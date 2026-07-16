"""v9.18.5 long-growth runtime: target stop, strict quality, refined corridor.

This revision removes the v9.18.4 sub-degree sliver workaround.  Exact edge-aware
insertion from v9.18.3 is retained, but mechanically separated bulk components
receive one minimal incremental x anchor so the full Mode-I body has no rigid
horizontal mode after complete crack-face separation.

The revision also:

* exits the accepted-step loop immediately after the committed physical target;
* rejects topology events below the production triangle-quality/child-area floors;
* rejects events whose tip-local resolution exceeds a configured fraction of da;
* pre-refines a Mode-I corridor through the requested committed growth length.

No barrier, hazard, cohesive-opening, source-refresh, MPZ, wake, or shielding law
is changed.
"""
from __future__ import annotations

import json
import math
import os
from pathlib import Path
import sys
from typing import Any

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import spsolve

from . import crack_backend as _cb
from . import fem as _fem
from . import mesh as _mesh
from . import mode_i_first_passage_v9_18_1 as _v9181
from . import mode_i_first_passage_v9_18_3 as _v9183
from . import sharp_front as _sharp_front


_RUNTIME: dict[str, Any] = {
    "mesh": None,
    "controller": None,
    "corridor_centers": [],
    "component_anchor_history": [],
    "quality_vetoes": [],
}


class _DynamicStepHorizon:
    """Numeric-like horizon whose loop comparison can stop on committed target.

    ``sharp_front.run_2d`` evaluates ``while step < args.steps`` every accepted
    step.  This object preserves the original numeric horizon for calculations
    and serialization while returning False from that comparison once the
    controller requests a committed-target stop.
    """

    def __init__(self, value: int, controller: Any):
        self.value = int(value)
        self.controller = controller

    def _running(self) -> bool:
        return not bool(getattr(self.controller, "v9185_stop_requested", False))

    def __int__(self):
        return self.value

    def __index__(self):
        return self.value

    def __float__(self):
        return float(self.value)

    def __repr__(self):
        return str(self.value)

    def __str__(self):
        return str(self.value)

    # Reflected comparison used by ``integer_step < horizon``.
    def __gt__(self, other):
        return self._running() and float(other) < float(self.value)

    def __ge__(self, other):
        return self._running() and float(other) <= float(self.value)

    def __lt__(self, other):
        return float(self.value) < float(other)

    def __le__(self, other):
        return float(self.value) <= float(other)

    def __eq__(self, other):
        try:
            return float(self.value) == float(other)
        except Exception:
            return False

    def __mul__(self, other):
        return self.value * other

    def __rmul__(self, other):
        return other * self.value

    def __truediv__(self, other):
        return self.value / other

    def __rtruediv__(self, other):
        return other / self.value

    def __floordiv__(self, other):
        return self.value // other

    def __rfloordiv__(self, other):
        return other // self.value

    def __sub__(self, other):
        return self.value - other

    def __rsub__(self, other):
        return other - self.value

    def __add__(self, other):
        return self.value + other

    def __radd__(self, other):
        return other + self.value


class TargetStopPersistentWakeController(
    _v9181.RenewalRollbackPersistentWakeController
):
    """v9.18.1 controller with a post-commit accepted-loop stop request."""

    def __init__(self) -> None:
        super().__init__()
        self.v9185_stop_requested = False
        self.v9185_target_stop_requests = 0

    def _commit_deferred_renewal(self) -> dict[str, float]:
        wake = super()._commit_deferred_renewal()
        if bool(getattr(self, "committed_target_reached", False)):
            self.v9185_stop_requested = True
            self.v9185_target_stop_requests += 1
            wake["v9185_committed_target_stop_requested"] = 1.0
        else:
            wake["v9185_committed_target_stop_requested"] = 0.0
        return wake

    def payload(self) -> dict[str, Any]:
        data = super().payload()
        data.update({
            "v9185_immediate_committed_target_stop_enabled": True,
            "v9185_stop_requested": bool(self.v9185_stop_requested),
            "v9185_target_stop_requests": int(self.v9185_target_stop_requests),
        })
        return data


def _corridor_centers(geom: Any, mesh_cfg: Any) -> np.ndarray:
    target_um = float(os.environ.get("ARRHENIUS_COMMITTED_TARGET_EXTENSION_UM", "0"))
    guard_um = max(float(os.environ.get("ARRHENIUS_CORRIDOR_GUARD_UM", "10")), 0.0)
    spacing_um = max(float(os.environ.get("ARRHENIUS_CORRIDOR_CENTER_SPACING_UM", "25")), 1.0)
    length_m = max(target_um + guard_um, 0.0) * 1.0e-6
    start = float(geom.a0)
    stop = min(float(geom.Lx), start + length_m)
    if stop <= start + 1.0e-15:
        return np.array([[start, 0.0]], dtype=float)
    n = max(int(math.ceil((stop - start) / (spacing_um * 1.0e-6))) + 1, 2)
    xs = np.linspace(start, stop, n)
    return np.column_stack([xs, np.zeros_like(xs)])


def _make_corridor_mesh(geom, mesh_cfg, seed=None, tip_center=None):
    original = _make_corridor_mesh._original
    graded = float(getattr(mesh_cfg, "tip_h_fine", 0.0) or 0.0) > 0.0
    enabled = os.environ.get("ARRHENIUS_PREFINED_MODE_I_CORRIDOR", "1") not in {
        "0", "false", "False", "no", "NO"
    }
    if tip_center is not None or not graded or not enabled:
        return original(geom, mesh_cfg, seed=seed, tip_center=tip_center)
    centers = _corridor_centers(geom, mesh_cfg)
    _RUNTIME["corridor_centers"] = centers.tolist()
    return original(geom, mesh_cfg, seed=seed, tip_center=centers)


def _assemble_with_mesh_capture(mesh, *args, **kwargs):
    _RUNTIME["mesh"] = mesh
    return _assemble_with_mesh_capture._original(mesh, *args, **kwargs)


def _node_components(mesh) -> tuple[int, np.ndarray]:
    nn = int(mesh.nn)
    conn = np.asarray(mesh.elems, dtype=int)
    edges = np.vstack([conn[:, [0, 1]], conn[:, [1, 2]], conn[:, [2, 0]]])
    rows = np.concatenate([edges[:, 0], edges[:, 1]])
    cols = np.concatenate([edges[:, 1], edges[:, 0]])
    graph = sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(nn, nn))
    return connected_components(graph, directed=False, return_labels=True)


def _component_anchored_solve(
    K: sparse.csr_matrix,
    Rint: np.ndarray,
    u: np.ndarray,
    bnd: Any,
    Uy_top: float,
    Uy_bot: float,
):
    """Dirichlet solve with one incremental x anchor per separated component."""
    mesh = _RUNTIME.get("mesh")
    if mesh is None or int(mesh.nn) != len(u) // 2:
        raise RuntimeError("v9.18.5 mechanics solve has no matching captured mesh")

    nn = len(u) // 2
    ndof = 2 * nn
    prescribed = np.zeros(ndof, dtype=bool)
    u_pres = np.zeros(ndof, dtype=float)

    top = np.asarray(bnd.top_nodes, dtype=int)
    bot = np.asarray(bnd.bot_nodes, dtype=int)
    prescribed[2 * top + 1] = True
    u_pres[2 * top + 1] = Uy_top
    prescribed[2 * bot + 1] = True
    u_pres[2 * bot + 1] = Uy_bot

    prescribed[2 * int(bnd.left_bot)] = True
    u_pres[2 * int(bnd.left_bot)] = 0.0
    prescribed[2 * int(bnd.left_bot) + 1] = True
    u_pres[2 * int(bnd.left_bot) + 1] = Uy_bot
    prescribed[2 * int(bnd.right_bot)] = True
    u_pres[2 * int(bnd.right_bot)] = 0.0

    ncomp, labels = _node_components(mesh)
    anchors = []
    boundary_nodes = np.unique(np.concatenate([top, bot]))
    for comp in range(int(ncomp)):
        nodes = np.where(labels == comp)[0]
        if nodes.size == 0:
            continue
        y_dofs = 2 * nodes + 1
        x_dofs = 2 * nodes
        if not bool(np.any(prescribed[y_dofs])):
            raise RuntimeError(
                f"v9.18.5 floating bulk component {comp} has no vertical grip anchor"
            )
        if bool(np.any(prescribed[x_dofs])):
            continue
        candidates = np.intersect1d(nodes, boundary_nodes, assume_unique=False)
        if candidates.size == 0:
            raise RuntimeError(
                f"v9.18.5 bulk component {comp} has no grip-boundary node for x anchor"
            )
        xy = np.asarray(mesh.nodes[candidates], float)
        order = np.lexsort((np.abs(xy[:, 1]), xy[:, 0]))
        nid = int(candidates[int(order[0])])
        prescribed[2 * nid] = True
        # Hold the current x position.  This removes only the incremental rigid
        # translation and introduces no displacement jump at component separation.
        u_pres[2 * nid] = float(u[2 * nid])
        anchors.append({
            "component": int(comp),
            "node": nid,
            "x_m": float(mesh.nodes[nid, 0]),
            "y_m": float(mesh.nodes[nid, 1]),
            "held_ux_m": float(u_pres[2 * nid]),
        })

    free = ~prescribed
    K_csr = K.tocsr()
    du_pres = u_pres[prescribed] - u[prescribed]
    rhs = -Rint[free] - K_csr[np.ix_(free, prescribed)] @ du_pres
    u_new = u.copy()
    u_new[free] = u[free] + spsolve(K_csr[np.ix_(free, free)], rhs)
    u_new[prescribed] = u_pres[prescribed]
    Rfull = Rint + K_csr @ (u_new - u)
    Ftop = float(np.sum(Rfull[2 * top + 1]))

    if not np.all(np.isfinite(u_new)) or not math.isfinite(Ftop):
        raise RuntimeError("v9.18.5 non-finite component-anchored FEM solution")
    if anchors:
        _RUNTIME["component_anchor_history"].append(anchors)
    return u_new, Ftop


def _affected_elements(old_mesh, new_mesh) -> np.ndarray:
    affected: set[int] = set(range(int(old_mesh.ne), int(new_mesh.ne)))
    common_e = min(int(old_mesh.ne), int(new_mesh.ne))
    if common_e:
        changed_conn = np.where(np.any(
            np.asarray(old_mesh.elems[:common_e]) != np.asarray(new_mesh.elems[:common_e]),
            axis=1,
        ))[0]
        affected.update(int(x) for x in changed_conn)
    common_n = min(int(old_mesh.nn), int(new_mesh.nn))
    if common_n:
        moved = np.where(np.linalg.norm(
            np.asarray(new_mesh.nodes[:common_n]) - np.asarray(old_mesh.nodes[:common_n]),
            axis=1,
        ) > 1.0e-14)[0]
        if moved.size:
            mask = np.any(np.isin(np.asarray(new_mesh.elems[:common_e]), moved), axis=1)
            affected.update(int(x) for x in np.where(mask)[0])
    if not affected:
        affected.update(range(int(new_mesh.ne)))
    return np.asarray(sorted(affected), dtype=int)


def _strict_quality_advance(self, *args, **kwargs):
    original = _strict_quality_advance._original
    old_mesh = kwargs["mesh"]
    snap = self._transaction_snapshot()
    result = original(self, *args, **kwargs)
    if not bool(getattr(result, "inserted", False)):
        return result

    new_mesh = result.mesh
    affected = _affected_elements(old_mesh, new_mesh)
    quality = self._triangle_quality(new_mesh.nodes, new_mesh.elems[affected])
    qmin = float(np.min(quality)) if quality.size else 1.0

    if not hasattr(self, "_v9185_production_q_floor"):
        self._v9185_production_q_floor = float(self.min_triangle_quality)
        self._v9185_production_area_floor = float(self.min_area_ratio)
    qfloor = float(os.environ.get(
        "ARRHENIUS_MIN_ACCEPTED_TRIANGLE_QUALITY",
        str(self._v9185_production_q_floor),
    ))
    afloor = float(os.environ.get(
        "ARRHENIUS_MIN_ACCEPTED_CHILD_AREA_RATIO",
        str(self._v9185_production_area_floor),
    ))

    parent_map = getattr(result, "elem_parent_map", None)
    ratios = []
    if parent_map is not None:
        pm = np.asarray(parent_map, dtype=int)
        for e in affected:
            if e < len(pm) and 0 <= pm[e] < int(old_mesh.ne):
                ratios.append(float(new_mesh.area_e[e]) / max(float(old_mesh.area_e[pm[e]]), 1.0e-300))
    else:
        for e in affected:
            if e < int(old_mesh.ne):
                ratios.append(float(new_mesh.area_e[e]) / max(float(old_mesh.area_e[e]), 1.0e-300))
    amin = float(min(ratios)) if ratios else 1.0

    p0 = np.asarray(kwargs.get("p0"), float)
    p1 = np.asarray(kwargs.get("p1"), float)
    da = max(float(np.linalg.norm(p1 - p0)), 1.0e-300)
    tip_ratio = float(new_mesh.hbar_tip) / da
    max_tip_ratio = float(os.environ.get("ARRHENIUS_MAX_TIP_H_OVER_DA", "0.75"))

    incidence = np.bincount(np.asarray(new_mesh.elems, int).ravel(), minlength=int(new_mesh.nn))
    orphan = np.where(incidence <= 0)[0]
    bad_endpoint = []
    for elem in self.cohesive_network.elements:
        for nid in elem.nodes4:
            if int(nid) < 0 or int(nid) >= int(new_mesh.nn) or incidence[int(nid)] <= 0:
                bad_endpoint.append(int(nid))

    issues = []
    if not np.all(np.isfinite(new_mesh.area_e)) or np.any(new_mesh.area_e <= 0.0):
        issues.append("nonpositive_or_nonfinite_area")
    if qmin < qfloor:
        issues.append(f"triangle_quality={qmin:.6e}<{qfloor:.6e}")
    if amin < afloor:
        issues.append(f"child_area_ratio={amin:.6e}<{afloor:.6e}")
    if tip_ratio > max_tip_ratio:
        issues.append(f"tip_h_over_da={tip_ratio:.6e}>{max_tip_ratio:.6e}")
    if orphan.size:
        issues.append(f"orphan_bulk_nodes={orphan[:10].tolist()}")
    if bad_endpoint:
        issues.append(f"unsupported_cohesive_endpoints={bad_endpoint[:10]}")

    row = {
        "front_id": int(kwargs.get("front_id", -1)),
        "min_triangle_quality": qmin,
        "triangle_quality_floor": qfloor,
        "min_child_area_ratio": amin,
        "child_area_ratio_floor": afloor,
        "tip_h_over_da": tip_ratio,
        "max_tip_h_over_da": max_tip_ratio,
        "affected_element_count": int(len(affected)),
        "accepted": not bool(issues),
        "issues": issues,
    }
    if not issues:
        if getattr(self, "advance_log", None):
            self.advance_log[-1].update({
                "v9185_quality_gate_passed": True,
                "v9185_min_triangle_quality": qmin,
                "v9185_min_child_area_ratio": amin,
                "v9185_tip_h_over_da": tip_ratio,
            })
        return result

    _RUNTIME["quality_vetoes"].append(row)
    self._transaction_rollback(snap)
    return _cb.CrackAdvanceResult(
        mesh=old_mesh,
        boundary=kwargs["boundary"],
        damage=kwargs["damage"],
        displacement=kwargs["displacement"],
        moved=0.0,
        inserted=False,
        angle_error_deg=float(getattr(result, "angle_error_deg", 0.0)),
        reason="v9185_quality_veto:" + ";".join(issues),
        elem_parent_map=None,
    )


def _run_2d_with_dynamic_target_stop(args):
    original = _run_2d_with_dynamic_target_stop._original
    controller = _RUNTIME.get("controller")
    original_steps = args.steps
    args.steps = _DynamicStepHorizon(int(original_steps), controller)
    try:
        return original(args)
    finally:
        args.steps = original_steps


def _option_value(argv: list[str], name: str) -> str | None:
    for i, token in enumerate(argv):
        if token == name and i + 1 < len(argv):
            return argv[i + 1]
        if token.startswith(name + "="):
            return token.split("=", 1)[1]
    return None


def main(argv=None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    original_controller_cls = _v9181.RenewalRollbackPersistentWakeController
    original_make_mesh = _mesh.make_tri_mesh
    original_assemble = _fem.assemble_mechanics
    original_solve = _fem.solve_dirichlet
    original_insert = _cb.AdaptiveCZMBackend._insert_target_in_incident_triangle
    original_advance = _cb.AdaptiveCZMBackend.advance
    original_run_2d = _sharp_front.run_2d

    controller_holder: dict[str, Any] = {}

    class _Controller(TargetStopPersistentWakeController):
        def __init__(self):
            super().__init__()
            controller_holder["controller"] = self
            _RUNTIME["controller"] = self

    _make_corridor_mesh._original = original_make_mesh
    _assemble_with_mesh_capture._original = original_assemble
    _strict_quality_advance._original = original_advance
    _run_2d_with_dynamic_target_stop._original = original_run_2d

    _v9181.RenewalRollbackPersistentWakeController = _Controller
    _mesh.make_tri_mesh = _make_corridor_mesh
    _fem.assemble_mechanics = _assemble_with_mesh_capture
    _fem.solve_dirichlet = _component_anchored_solve
    _cb.AdaptiveCZMBackend._insert_target_in_incident_triangle = (
        _v9183._edge_aware_insert_target_in_incident_triangle
    )
    _cb.AdaptiveCZMBackend.advance = _strict_quality_advance
    _sharp_front.run_2d = _run_2d_with_dynamic_target_stop

    try:
        results = _v9181.main(user_args)
    finally:
        _v9181.RenewalRollbackPersistentWakeController = original_controller_cls
        _mesh.make_tri_mesh = original_make_mesh
        _fem.assemble_mechanics = original_assemble
        _fem.solve_dirichlet = original_solve
        _cb.AdaptiveCZMBackend._insert_target_in_incident_triangle = original_insert
        _cb.AdaptiveCZMBackend.advance = original_advance
        _sharp_front.run_2d = original_run_2d

    out_value = _option_value(user_args, "--out")
    if out_value is not None:
        out = Path(out_value)
        out.mkdir(parents=True, exist_ok=True)
        controller = controller_holder.get("controller")
        payload = {
            "schema": "target_stop_quality_corridor_v9185_v1",
            "immediate_committed_target_stop_enabled": True,
            "target_stop_requested": bool(getattr(controller, "v9185_stop_requested", False)),
            "target_stop_requests": int(getattr(controller, "v9185_target_stop_requests", 0)),
            "exact_edge_aware_insertion_enabled": True,
            "same_length_angular_regularization_enabled": False,
            "component_wise_incremental_x_anchor_enabled": True,
            "prefined_mode_i_corridor_enabled": True,
            "corridor_centers_m": _RUNTIME.get("corridor_centers", []),
            "component_anchor_history": _RUNTIME.get("component_anchor_history", []),
            "strict_quality_gate_enabled": True,
            "quality_vetoes": _RUNTIME.get("quality_vetoes", []),
            "constitutive_physics_changed": False,
        }
        (out / "target_stop_quality_corridor_v9185.json").write_text(
            json.dumps(payload, indent=2, default=str)
        )
    return results


if __name__ == "__main__":
    main()
