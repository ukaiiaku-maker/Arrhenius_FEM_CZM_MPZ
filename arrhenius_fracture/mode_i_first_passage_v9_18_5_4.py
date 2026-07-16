"""v9.18.5.4 active-tip resolution gate and integrated veto guard.

v9.18.5 evaluated ``h_tip / da`` using ``TriMesh.hbar_tip``.  After topology
updates that stored quantity is an average over the nearest two percent of all
bulk elements; it is not a one-ring measure at the active cohesive endpoint.
The 700 K ceramic gate therefore rejected a valid third event with
``hbar_tip/da = 1.097`` and retried it thousands of times.

This revision evaluates resolution from the unique, nonzero bulk edges in the
active endpoint one-ring.  The accepted mesh stores that local value as its new
``hbar_tip``.  Repeated deterministic vetoes are counted in the same active
strict-quality wrapper and fail explicitly at the configured limit.

No barrier, hazard, cohesive, MPZ, wake, shielding, or material law changes.
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
from . import mode_i_first_passage_v9_18_5 as _v9185
from . import mode_i_first_passage_v9_18_5_3 as _v91853


_AUDIT: dict[str, Any] = {
    "accepted_events": [],
    "vetoes": [],
    "identical_veto_abort": None,
}


def _active_endpoint(self: Any, result: Any, kwargs: dict[str, Any]) -> np.ndarray:
    front_id = int(kwargs.get("front_id", -1))
    tips = getattr(self, "tip_nodes", {})
    if front_id in tips:
        return np.asarray(tips[front_id][2], float).copy()
    return np.asarray(kwargs.get("p1"), float).copy()


def _active_tip_one_ring_resolution(mesh: Any, point: np.ndarray) -> dict[str, Any]:
    """Return a geometric one-ring edge-size audit at ``point``.

    Coincident CZM node copies are grouped geometrically.  Zero-length edges are
    excluded, and each geometric bulk edge is counted once.
    """
    nodes = np.asarray(mesh.nodes, float)
    elems = np.asarray(mesh.elems, int)
    point = np.asarray(point, float)
    scale = max(float(getattr(mesh, "hbar", 0.0) or 0.0), 1.0e-12)
    tol = max(1.0e-12, 1.0e-7 * scale)
    distance = np.linalg.norm(nodes - point[None, :], axis=1)
    dmin = float(np.min(distance)) if len(distance) else float("inf")
    endpoint_nodes = np.where(distance <= dmin + tol)[0]
    if endpoint_nodes.size == 0:
        raise RuntimeError("active cohesive endpoint has no represented bulk node")

    incident_mask = np.any(np.isin(elems, endpoint_nodes), axis=1)
    incident = np.where(incident_mask)[0]
    if incident.size == 0:
        raise RuntimeError("active cohesive endpoint has no incident bulk elements")

    geom_edges: dict[tuple[tuple[int, int], tuple[int, int]], float] = {}
    quant = max(tol, 1.0e-14)
    for tri in elems[incident]:
        for ia, ib in ((0, 1), (1, 2), (2, 0)):
            a = nodes[int(tri[ia])]
            b = nodes[int(tri[ib])]
            length = float(np.linalg.norm(b - a))
            if length <= tol * 1.0e-3:
                continue
            ka = tuple(np.round(a / quant).astype(np.int64))
            kb = tuple(np.round(b / quant).astype(np.int64))
            key = tuple(sorted((ka, kb)))
            geom_edges[key] = length

    if not geom_edges:
        raise RuntimeError("active cohesive endpoint one-ring has no nonzero bulk edges")
    lengths = np.asarray(list(geom_edges.values()), float)
    return {
        "active_tip_point_m": point.tolist(),
        "active_tip_endpoint_node_count": int(endpoint_nodes.size),
        "active_tip_incident_element_count": int(incident.size),
        "active_tip_unique_edge_count": int(lengths.size),
        "active_tip_h_mean_m": float(np.mean(lengths)),
        "active_tip_h_median_m": float(np.median(lengths)),
        "active_tip_h_p90_m": float(np.quantile(lengths, 0.90)),
        "active_tip_h_max_m": float(np.max(lengths)),
        "legacy_stored_hbar_tip_m": float(getattr(mesh, "hbar_tip", float("nan"))),
    }


def _veto_signature(kwargs: dict[str, Any], reason: str) -> tuple[Any, ...]:
    p0 = np.asarray(kwargs.get("p0", [math.nan, math.nan]), float)
    p1 = np.asarray(kwargs.get("p1", [math.nan, math.nan]), float)
    return (
        int(kwargs.get("front_id", -1)),
        tuple(np.round(p0, 12)),
        tuple(np.round(p1, 12)),
        str(reason),
    )


def _record_veto_or_raise(self: Any, kwargs: dict[str, Any], result: Any) -> Any:
    if bool(getattr(result, "inserted", False)):
        self._v91854_last_veto_signature = None
        self._v91854_identical_veto_count = 0
        return result

    reason = str(getattr(result, "reason", "unknown"))
    signature = _veto_signature(kwargs, reason)
    if signature == getattr(self, "_v91854_last_veto_signature", None):
        count = int(getattr(self, "_v91854_identical_veto_count", 0)) + 1
    else:
        count = 1
    self._v91854_last_veto_signature = signature
    self._v91854_identical_veto_count = count
    limit = max(int(os.environ.get("ARRHENIUS_MAX_IDENTICAL_GEOMETRY_VETOES", "12")), 1)
    if count >= limit:
        payload = {
            "front_id": int(kwargs.get("front_id", -1)),
            "count": count,
            "limit": limit,
            "reason": reason,
            "p0_m": np.asarray(kwargs.get("p0"), float).tolist(),
            "p1_m": np.asarray(kwargs.get("p1"), float).tolist(),
        }
        _AUDIT["identical_veto_abort"] = payload
        raise RuntimeError(
            "v9.18.5.4 repeated identical geometry veto; physical renewal remains "
            f"unconsumed: front={payload['front_id']} count={count}/{limit} reason={reason}"
        )
    return result


def _strict_quality_advance_v91854(self: Any, *args, **kwargs):
    original = _strict_quality_advance_v91854._original
    old_mesh = kwargs["mesh"]
    snap = self._transaction_snapshot()
    result = original(self, *args, **kwargs)
    if not bool(getattr(result, "inserted", False)):
        return _record_veto_or_raise(self, kwargs, result)

    new_mesh = result.mesh
    affected = _v9185._affected_elements(old_mesh, new_mesh)
    quality = self._triangle_quality(new_mesh.nodes, new_mesh.elems[affected])
    qmin = float(np.min(quality)) if quality.size else 1.0

    if not hasattr(self, "_v91854_production_q_floor"):
        self._v91854_production_q_floor = float(self.min_triangle_quality)
        self._v91854_production_area_floor = float(self.min_area_ratio)
    qfloor = float(os.environ.get(
        "ARRHENIUS_MIN_ACCEPTED_TRIANGLE_QUALITY",
        str(self._v91854_production_q_floor),
    ))
    afloor = float(os.environ.get(
        "ARRHENIUS_MIN_ACCEPTED_CHILD_AREA_RATIO",
        str(self._v91854_production_area_floor),
    ))

    parent_map = getattr(result, "elem_parent_map", None)
    ratios: list[float] = []
    if parent_map is not None:
        pm = np.asarray(parent_map, dtype=int)
        for e in affected:
            if e < len(pm) and 0 <= pm[e] < int(old_mesh.ne):
                ratios.append(
                    float(new_mesh.area_e[e]) /
                    max(float(old_mesh.area_e[pm[e]]), 1.0e-300)
                )
    else:
        for e in affected:
            if e < int(old_mesh.ne):
                ratios.append(
                    float(new_mesh.area_e[e]) /
                    max(float(old_mesh.area_e[e]), 1.0e-300)
                )
    amin = float(min(ratios)) if ratios else 1.0

    p0 = np.asarray(kwargs.get("p0"), float)
    p1 = np.asarray(kwargs.get("p1"), float)
    da = max(float(np.linalg.norm(p1 - p0)), 1.0e-300)
    endpoint = _active_endpoint(self, result, kwargs)
    local = _active_tip_one_ring_resolution(new_mesh, endpoint)
    tip_ratio = float(local["active_tip_h_mean_m"]) / da
    max_tip_ratio = float(os.environ.get("ARRHENIUS_MAX_TIP_H_OVER_DA", "0.75"))

    incidence = np.bincount(
        np.asarray(new_mesh.elems, int).ravel(), minlength=int(new_mesh.nn)
    )
    orphan = np.where(incidence <= 0)[0]
    bad_endpoint: list[int] = []
    for elem in self.cohesive_network.elements:
        for nid in elem.nodes4:
            if int(nid) < 0 or int(nid) >= int(new_mesh.nn) or incidence[int(nid)] <= 0:
                bad_endpoint.append(int(nid))

    issues: list[str] = []
    if not np.all(np.isfinite(new_mesh.area_e)) or np.any(new_mesh.area_e <= 0.0):
        issues.append("nonpositive_or_nonfinite_area")
    if qmin < qfloor:
        issues.append(f"triangle_quality={qmin:.6e}<{qfloor:.6e}")
    if amin < afloor:
        issues.append(f"child_area_ratio={amin:.6e}<{afloor:.6e}")
    if tip_ratio > max_tip_ratio:
        issues.append(f"active_tip_h_over_da={tip_ratio:.6e}>{max_tip_ratio:.6e}")
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
        "active_tip_h_over_da": tip_ratio,
        "max_tip_h_over_da": max_tip_ratio,
        "requested_da_m": da,
        "affected_element_count": int(len(affected)),
        "accepted": not bool(issues),
        "issues": issues,
        "resolution_metric": "active_endpoint_one_ring_unique_nonzero_edge_mean",
        **local,
    }

    if not issues:
        # The active endpoint one-ring is the relevant process-zone resolution for
        # subsequent topology tolerances and diagnostics.
        new_mesh.hbar_tip = float(local["active_tip_h_mean_m"])
        if getattr(self, "advance_log", None):
            self.advance_log[-1].update({
                "v91854_quality_gate_passed": True,
                "v91854_min_triangle_quality": qmin,
                "v91854_min_child_area_ratio": amin,
                "v91854_active_tip_h_over_da": tip_ratio,
                "v91854_active_tip_h_mean_m": local["active_tip_h_mean_m"],
                "v91854_legacy_stored_hbar_tip_m": local["legacy_stored_hbar_tip_m"],
            })
        _AUDIT["accepted_events"].append(row)
        return _record_veto_or_raise(self, kwargs, result)

    _AUDIT["vetoes"].append(row)
    _v9185._RUNTIME["quality_vetoes"].append(row)
    self._transaction_rollback(snap)
    veto = _cb.CrackAdvanceResult(
        mesh=old_mesh,
        boundary=kwargs["boundary"],
        damage=kwargs["damage"],
        displacement=kwargs["displacement"],
        moved=0.0,
        inserted=False,
        angle_error_deg=float(getattr(result, "angle_error_deg", 0.0)),
        reason="v91854_quality_veto:" + ";".join(issues),
        elem_parent_map=None,
    )
    return _record_veto_or_raise(self, kwargs, veto)


def _option_value(argv: list[str], name: str) -> str | None:
    for i, token in enumerate(argv):
        if token == name and i + 1 < len(argv):
            return argv[i + 1]
        if token.startswith(name + "="):
            return token.split("=", 1)[1]
    return None


def _write_audit(argv: list[str], error: BaseException | None) -> None:
    out_value = _option_value(argv, "--out")
    if out_value is None:
        return
    out = Path(out_value)
    out.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema": "active_tip_resolution_veto_guard_v91854_v1",
        "active_tip_resolution_metric": "one_ring_unique_nonzero_bulk_edge_mean",
        "legacy_nearest_two_percent_hbar_tip_not_used_for_event_gate": True,
        "identical_veto_guard_integrated_with_active_quality_wrapper": True,
        "max_identical_geometry_vetoes": int(
            os.environ.get("ARRHENIUS_MAX_IDENTICAL_GEOMETRY_VETOES", "12")
        ),
        "run_completed_without_exception": error is None,
        "runtime_error_type": None if error is None else type(error).__name__,
        "runtime_error": None if error is None else str(error),
        "constitutive_physics_changed": False,
        **_AUDIT,
    }
    (out / "active_tip_resolution_veto_guard_v91854.json").write_text(
        json.dumps(payload, indent=2, default=str)
    )


def main(argv=None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    _AUDIT["accepted_events"] = []
    _AUDIT["vetoes"] = []
    _AUDIT["identical_veto_abort"] = None
    original = _v9185._strict_quality_advance
    _v9185._strict_quality_advance = _strict_quality_advance_v91854
    error: BaseException | None = None
    try:
        return _v91853.main(user_args)
    except BaseException as exc:
        error = exc
        raise
    finally:
        _write_audit(user_args, error)
        _v9185._strict_quality_advance = original


if __name__ == "__main__":
    main()
