"""v9.18.4 mechanically valid crack-front regularization.

v9.18.3 correctly avoided a zero-area child when a requested crack endpoint lay
on an existing mesh edge.  The accepted exact edge endpoint nevertheless left a
mechanically singular bulk topology after the 50 um event in all three 700 K
classes.  The next FEM solve returned NaN displacement, which subsequently
entered MPZ advection.

v9.18.4 preserves the requested physical advance length but, only when the
endpoint lies numerically on an edge, rotates the local topology segment by a
small bounded angle into an adjacent triangle.  Candidate endpoints are ranked
by child-element quality and tried transactionally.  An accepted event must also
pass bulk-incidence and boundary-anchor connectivity checks before the backend
transaction is committed.

No hazard, barrier, source-refresh, wake, shielding, or opening law is changed.
"""
from __future__ import annotations

import json
import math
import os
from pathlib import Path
import sys
from typing import Any

import numpy as np

from . import crack_backend as _cb
from . import fem as _fem
from . import mode_i_first_passage_v9_18_1 as _v9181
from . import mode_i_first_passage_v9_18_3 as _v9183


def _angle_schedule_deg() -> list[float]:
    raw = os.environ.get(
        "ARRHENIUS_EDGE_FRONT_REGULARIZATION_ANGLES_DEG",
        "0.01 0.025 0.05 0.1 0.2 0.35",
    )
    values: list[float] = []
    for token in str(raw).replace(",", " ").split():
        try:
            value = abs(float(token))
        except (TypeError, ValueError):
            continue
        if math.isfinite(value) and value > 0.0:
            values.append(value)
    return sorted(set(values)) or [0.01, 0.025, 0.05, 0.1, 0.2, 0.35]


def _mechanically_regularized_insert_target(
    self,
    mesh,
    displacement: np.ndarray,
    p0: np.ndarray,
    target: np.ndarray,
    front_id: int,
):
    """Use an interior same-length endpoint instead of an edge-front endpoint."""
    base = _v9183._edge_aware_insert_target_in_incident_triangle(
        self, mesh, displacement, p0, target, front_id
    )
    refined, u1, q, reason, meta, parent_map = base
    location = str((meta or {}).get("target_location_case", ""))
    if location != "existing_edge_exact_split":
        return base

    p0 = np.asarray(p0, float)
    target = np.asarray(target, float)
    vec = target - p0
    length = float(np.linalg.norm(vec))
    if not math.isfinite(length) or length <= 1.0e-14:
        return None, None, None, "edge_regularization_zero_length", {}, None
    tangent = vec / length
    normal = np.array([-tangent[1], tangent[0]], dtype=float)

    candidates: list[tuple[float, float, float, tuple]] = []
    for angle_deg in _angle_schedule_deg():
        theta = math.radians(float(angle_deg))
        c = math.cos(theta)
        s = math.sin(theta)
        for sign in (-1.0, 1.0):
            direction = c * tangent + sign * s * normal
            candidate = p0 + length * direction
            trial = _v9183._edge_aware_insert_target_in_incident_triangle(
                self, mesh, displacement, p0, candidate, front_id
            )
            tmesh, tu, tq, treason, tmeta, tmap = trial
            if tmesh is None or treason != "ok":
                continue
            if str((tmeta or {}).get("target_location_case", "")) != "strict_triangle_interior":
                continue
            quality = float((tmeta or {}).get("min_triangle_quality", float("nan")))
            area_ratio = float((tmeta or {}).get("min_area_ratio", float("nan")))
            if not math.isfinite(quality) or not math.isfinite(area_ratio):
                continue
            enriched = dict(tmeta)
            enriched.update({
                "target_location_case": "same_length_interior_regularization",
                "edge_front_regularization_used": True,
                "requested_target_x": float(target[0]),
                "requested_target_y": float(target[1]),
                "regularized_target_x": float(candidate[0]),
                "regularized_target_y": float(candidate[1]),
                "regularization_angle_deg": float(sign * angle_deg),
                "requested_advance_length_m": float(length),
                "regularized_advance_length_m": float(np.linalg.norm(candidate - p0)),
                "advance_length_error_m": float(np.linalg.norm(candidate - p0) - length),
                "replaced_v9183_edge_split": True,
            })
            trial = (tmesh, tu, tq, treason, enriched, tmap)
            # Prefer best quality; then the smallest absolute angular perturbation.
            candidates.append((-quality, abs(float(angle_deg)), -area_ratio, trial))

    if not candidates:
        return None, None, None, "no_quality_safe_same_length_edge_regularization", {
            "target_location_case": "edge_regularization_failed",
            "requested_target_x": float(target[0]),
            "requested_target_y": float(target[1]),
        }, None

    candidates.sort(key=lambda item: (item[0], item[1], item[2]))
    attempt = max(int(getattr(self, "_v9184_regularization_attempt", 0)), 0)
    chosen = candidates[min(attempt, len(candidates) - 1)][3]
    chosen[4]["regularization_candidate_index"] = int(min(attempt, len(candidates) - 1))
    chosen[4]["regularization_candidate_count"] = int(len(candidates))
    return chosen


def _bulk_components(mesh) -> tuple[np.ndarray, int, np.ndarray]:
    parent = np.arange(mesh.nn, dtype=int)
    rank = np.zeros(mesh.nn, dtype=np.int8)

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = int(parent[x])
        return x

    def union(a: int, b: int) -> None:
        ra = find(int(a)); rb = find(int(b))
        if ra == rb:
            return
        if rank[ra] < rank[rb]:
            parent[ra] = rb
        elif rank[ra] > rank[rb]:
            parent[rb] = ra
        else:
            parent[rb] = ra
            rank[ra] += 1

    for tri in np.asarray(mesh.elems, dtype=int):
        union(int(tri[0]), int(tri[1]))
        union(int(tri[1]), int(tri[2]))
    roots = np.array([find(i) for i in range(mesh.nn)], dtype=int)
    unique, labels = np.unique(roots, return_inverse=True)
    incidence = np.bincount(np.asarray(mesh.elems, dtype=int).ravel(), minlength=mesh.nn)
    return labels, int(len(unique)), incidence


def _mechanical_topology_issues(self, result) -> list[str]:
    mesh = result.mesh
    issues: list[str] = []
    if not np.all(np.isfinite(mesh.nodes)):
        issues.append("nonfinite_node_coordinates")
    if not np.all(np.isfinite(mesh.area_e)) or np.any(mesh.area_e <= 0.0):
        issues.append("nonpositive_or_nonfinite_bulk_area")

    labels, ncomp, incidence = _bulk_components(mesh)
    orphan = np.where(incidence <= 0)[0]
    if len(orphan):
        issues.append("orphan_bulk_nodes:" + ",".join(str(int(x)) for x in orphan[:12]))

    network = getattr(self, "cohesive_network", None)
    if network is not None:
        bad_surface_nodes: list[int] = []
        for elem in network.elements:
            for node in elem.nodes4:
                if int(node) < 0 or int(node) >= mesh.nn or incidence[int(node)] <= 0:
                    bad_surface_nodes.append(int(node))
        if bad_surface_nodes:
            unique_bad = sorted(set(bad_surface_nodes))[:12]
            issues.append("cohesive_nodes_without_bulk_support:" + ",".join(map(str, unique_bad)))

    # Failed cohesive links have zero tensile/shear stiffness.  Therefore every
    # disconnected bulk component must be independently anchored by the actual
    # Dirichlet boundary conditions; otherwise the next stiffness matrix has a
    # rigid-body null mode.
    bnd = result.boundary
    x_pins = {int(bnd.left_bot), int(bnd.right_bot)}
    y_pins = set(int(x) for x in np.asarray(bnd.top_nodes, dtype=int))
    y_pins.update(int(x) for x in np.asarray(bnd.bot_nodes, dtype=int))
    y_pins.add(int(bnd.left_bot))
    unanchored: list[str] = []
    for comp in range(ncomp):
        nodes = np.where(labels == comp)[0]
        if not len(nodes):
            continue
        node_set = set(int(x) for x in nodes)
        has_x = bool(node_set & x_pins)
        has_y = bool(node_set & y_pins)
        if not (has_x and has_y):
            unanchored.append(
                f"component={comp},nodes={len(nodes)},x_pin={int(has_x)},y_pin={int(has_y)}"
            )
    if unanchored:
        issues.append("unanchored_bulk_components:" + "|".join(unanchored[:8]))
    return issues


def _advance_with_mechanical_validation(self, *args, **kwargs):
    original = _advance_with_mechanical_validation._original
    snap = self._transaction_snapshot()
    result = original(self, *args, **kwargs)
    if bool(getattr(result, "inserted", False)):
        issues = _mechanical_topology_issues(self, result)
        if not issues:
            self._v9184_regularization_attempt = 0
            self._v9184_last_veto_signature = None
            self._v9184_identical_veto_count = 0
            return result

        self._transaction_rollback(snap)
        self._v9184_regularization_attempt = int(
            getattr(self, "_v9184_regularization_attempt", 0)
        ) + 1
        result = _cb.CrackAdvanceResult(
            mesh=kwargs["mesh"],
            boundary=kwargs["boundary"],
            damage=kwargs["damage"],
            displacement=kwargs["displacement"],
            moved=0.0,
            inserted=False,
            angle_error_deg=float(getattr(result, "angle_error_deg", 0.0)),
            reason="mechanical_topology_veto:" + ";".join(issues),
            elem_parent_map=None,
        )

    reason = str(getattr(result, "reason", "unknown"))
    front_id = int(kwargs.get("front_id", -1))
    p0 = np.asarray(kwargs.get("p0", [math.nan, math.nan]), float)
    p1 = np.asarray(kwargs.get("p1", [math.nan, math.nan]), float)
    scale = max(float(getattr(kwargs.get("mesh", None), "hbar_tip", 1.0e-12)), 1.0e-12)
    signature = (
        front_id,
        tuple(np.round(p0 / scale, 8)),
        tuple(np.round(p1 / scale, 8)),
        reason,
        int(getattr(self, "_v9184_regularization_attempt", 0)),
    )
    if signature == getattr(self, "_v9184_last_veto_signature", None):
        count = int(getattr(self, "_v9184_identical_veto_count", 0)) + 1
    else:
        count = 1
    self._v9184_last_veto_signature = signature
    self._v9184_identical_veto_count = count
    limit = max(int(os.environ.get("ARRHENIUS_MAX_IDENTICAL_GEOMETRY_VETOES", "12")), 1)
    if count >= limit:
        raise RuntimeError(
            "v9.18.4 repeated mechanically invalid geometry; physical event "
            f"remains unconsumed: front={front_id} count={count}/{limit} "
            f"attempt={getattr(self, '_v9184_regularization_attempt', 0)} "
            f"reason={reason} p0={p0.tolist()} p1={p1.tolist()}"
        )
    return result


def _finite_solve_dirichlet(*args, **kwargs):
    u_new, reaction = _finite_solve_dirichlet._original(*args, **kwargs)
    if not np.all(np.isfinite(u_new)) or not math.isfinite(float(reaction)):
        raise RuntimeError(
            "non-finite FEM solution detected immediately after mechanics solve; "
            "a singular or ill-conditioned topology was not allowed to propagate "
            "into the MPZ kinetics"
        )
    return u_new, reaction


def _option_value(args: list[str], name: str) -> str | None:
    for i, token in enumerate(args):
        if token == name and i + 1 < len(args):
            return args[i + 1]
        if token.startswith(name + "="):
            return token.split("=", 1)[1]
    return None


def main(argv=None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    original_insert = _cb.AdaptiveCZMBackend._insert_target_in_incident_triangle
    original_advance = _cb.AdaptiveCZMBackend.advance
    original_solve = _fem.solve_dirichlet
    _advance_with_mechanical_validation._original = original_advance
    _finite_solve_dirichlet._original = original_solve
    _cb.AdaptiveCZMBackend._insert_target_in_incident_triangle = (
        _mechanically_regularized_insert_target
    )
    _cb.AdaptiveCZMBackend.advance = _advance_with_mechanical_validation
    _fem.solve_dirichlet = _finite_solve_dirichlet
    try:
        results = _v9181.main(user_args)
    finally:
        _cb.AdaptiveCZMBackend._insert_target_in_incident_triangle = original_insert
        _cb.AdaptiveCZMBackend.advance = original_advance
        _fem.solve_dirichlet = original_solve

    out_value = _option_value(user_args, "--out")
    if out_value is not None:
        out = Path(out_value)
        out.mkdir(parents=True, exist_ok=True)
        payload = {
            "schema": "mechanically_valid_front_regularization_v9184_v1",
            "same_length_interior_regularization_enabled": True,
            "regularization_angles_deg": _angle_schedule_deg(),
            "bulk_incidence_validation_enabled": True,
            "cohesive_endpoint_bulk_support_required": True,
            "independent_component_boundary_anchoring_required": True,
            "nonfinite_mechanics_fail_fast_enabled": True,
            "constitutive_physics_changed": False,
        }
        (out / "mechanical_topology_validation_v9184.json").write_text(
            json.dumps(payload, indent=2, default=str)
        )
    return results


if __name__ == "__main__":
    main()
