"""Event-resolved exact-ray local refinement for v10.0.5.18.3.7.

One hazard-selected physical event is preserved exactly. Short events receive a
bounded local endpoint refinement. Long events are decomposed only at actual
mesh-ray crossings already implied by the exact segment; no equal-length
physical partition retry is permitted. The complete crossing sequence is one
atomic transaction.
"""
from __future__ import annotations

import copy
import math
import os
from typing import Any, Callable

import numpy as np

from . import crack_backend as _cb
from . import committed_tip_resolution_audit_v10051836 as _audit36
from . import mode_i_first_passage_v9_18_5 as _v9185
from . import mode_i_first_passage_v9_18_5_6 as _v91856
from . import quality_aware_czm_patch_v10051835 as _patch
from . import quality_aware_czm_ray_handoff_v10051835 as _ray
from . import quality_aware_czm_retry_v10051835 as _local


MODEL_ID = "event_resolved_exact_ray_local_cavity_retry_v10_0_5_18_3_7"
SCHEMA = "v10.0.5.18.3.7_atomic_exact_ray_crossing_local_cavity_retry"
_AUDIT: dict[str, Any] = {
    "schema": SCHEMA,
    "model_id": MODEL_ID,
    "events": [],
}


def configure_legacy_quality_wrapper(wrapper: Callable) -> None:
    _local.configure_legacy_quality_wrapper(wrapper)


def reset_audit() -> None:
    _local.reset_audit()
    _AUDIT["events"] = []


def _max_crossings(self: Any) -> int:
    default = int(getattr(self, "max_hrefine_subsegments", 512) or 512)
    try:
        value = int(os.environ.get("ARRHENIUS_EVENT_RESOLVED_MAX_RAY_SEGMENTS", default))
    except ValueError:
        value = default
    return max(value, 1)


def _point_tolerance(mesh: Any, requested: float) -> float:
    return max(
        1.0e-12,
        1.0e-8 * max(
            float(requested),
            float(getattr(mesh, "hbar_tip", 0.0) or 0.0),
            1.0e-12,
        ),
    )


def _local_endpoint_reachable(self: Any, mesh: Any, p0: np.ndarray, p1: np.ndarray, front_id: int) -> bool:
    if _patch._target_tip_triangle(self, mesh, p0, p1, front_id) is not None:
        return True
    adjacent, _ = _ray._adjacent_target_triangle(self, mesh, p0, p1, front_id)
    return adjacent is not None


def _clear_transaction_audit(a0: int, w0: int, r0: int) -> None:
    audit = _v91856._AUDIT
    del audit["accepted_events"][a0:]
    del audit["resolution_warnings"][w0:]
    del _v9185._RUNTIME["quality_vetoes"][r0:]


def _state_from_result(result: Any) -> _local._State:
    parent = getattr(result, "elem_parent_map", None)
    if parent is None:
        parent = np.arange(int(result.mesh.ne), dtype=int)
    return _local._State(
        result.mesh,
        result.boundary,
        np.asarray(result.damage),
        np.asarray(result.displacement),
        np.asarray(parent, dtype=int),
    )


def _tip_point(self: Any, front_id: int, fallback: np.ndarray) -> np.ndarray:
    tips = getattr(self, "tip_nodes", {})
    if int(front_id) in tips:
        return np.asarray(tips[int(front_id)][2], dtype=float).copy()
    return np.asarray(fallback, dtype=float).copy()


def _mark_accepted_rows(
    *,
    a0: int,
    physical_event_index: int,
    requested: float,
    segment_count: int,
) -> None:
    rows = _v91856._AUDIT["accepted_events"][a0:]
    for index, row in enumerate(rows):
        row.update(
            v10051837_event_resolved=True,
            physical_event_index=int(physical_event_index),
            exact_ray_segment_index=int(index),
            exact_ray_segment_count=int(segment_count),
            physical_event_requested_length_m=float(requested),
            equal_length_partition_retry=False,
        )


def _mark_advance_log(
    self: Any,
    *,
    log0: int,
    physical_event_index: int,
    requested: float,
    segment_count: int,
    patch_levels: int,
) -> None:
    for index, row in enumerate(self.advance_log[log0:]):
        row.update(
            v10051837_event_resolved=True,
            v10051837_retry_kind="exact_ray_crossing_local_cavity",
            v10051837_physical_event_index=int(physical_event_index),
            v10051837_exact_ray_segment_index=int(index),
            v10051837_exact_ray_segment_count=int(segment_count),
            v10051837_patch_refinement_levels=int(patch_levels),
            v10051837_physical_event_requested_length_m=float(requested),
            v10051837_equal_length_partition_retry=False,
        )


def _event_record_base(p0: np.ndarray, p1: np.ndarray, direction: np.ndarray, front_id: int) -> dict[str, Any]:
    return {
        "physical_event_index": int(len(_AUDIT["events"])),
        "front_id": int(front_id),
        "p0_m": np.asarray(p0, dtype=float).tolist(),
        "p1_m": np.asarray(p1, dtype=float).tolist(),
        "direction": np.asarray(direction, dtype=float).tolist(),
        "requested_length_m": float(np.linalg.norm(np.asarray(p1) - np.asarray(p0))),
        "equal_length_partition_retry": False,
        "exact_physical_endpoint_preserved": True,
        "exact_selected_direction_preserved": True,
        "atomic_transaction": True,
        "mesh_only_handoffs": [],
        "ray_segments": [],
    }


def _exact_ray_transaction(
    original: Callable,
    self: Any,
    kwargs: dict[str, Any],
    p0: np.ndarray,
    p1: np.ndarray,
    direction: np.ndarray,
):
    requested = float(np.linalg.norm(p1 - p0))
    front_id = int(kwargs.get("front_id", 0))
    root_snap = self._transaction_snapshot()
    audit = _v91856._AUDIT
    a0 = len(audit["accepted_events"])
    w0 = len(audit["resolution_warnings"])
    r0 = len(_v9185._RUNTIME["quality_vetoes"])
    log0 = len(self.advance_log)

    root = _local._State(
        kwargs["mesh"],
        kwargs["boundary"],
        np.asarray(kwargs["damage"]),
        np.asarray(kwargs["displacement"]),
        np.arange(int(kwargs["mesh"].ne), dtype=int),
    )
    current = root
    point = np.asarray(p0, dtype=float).copy()
    final = np.asarray(p1, dtype=float).copy()
    moved_total = 0.0
    max_angle = 0.0
    total_patch_levels = 0
    event = _event_record_base(p0, p1, direction, front_id)
    event_index = int(event["physical_event_index"])
    failure = None
    force_ray_crossing = False

    # Mesh-only handoffs do not consume physical ray segments.  Permit bounded
    # additional controller iterations while retaining the independent hard
    # limit on actual exact-ray crossing segments.
    max_iterations = 2 * _max_crossings(self) + _local._patch_levels()
    for transaction_iteration in range(max_iterations):
        remaining = float(np.linalg.norm(final - point))
        tol = _point_tolerance(current.mesh, requested)
        if remaining <= tol:
            break
        if len(event["ray_segments"]) >= _max_crossings(self):
            failure = f"exact_ray_segment_limit:{_max_crossings(self)}"
            break

        local_final = False if force_ray_crossing else _local_endpoint_reachable(
            self, current.mesh, point, final, front_id
        )
        if local_final:
            segment_target = final
            segment_kind = "exact_physical_endpoint"
        else:
            hit = self._ray_exit_edge(
                current.mesh,
                point,
                direction,
                front_id,
                remaining,
            )
            if hit is None:
                failure = (
                    "no_exact_ray_exit_after_mesh_only_handoff"
                    if force_ray_crossing
                    else "no_exact_ray_exit_to_distant_endpoint"
                )
                break
            t_hit, edge_i, edge_j, q_hit, elem_id, xi_hit = hit
            if float(t_hit) >= remaining - tol:
                segment_target = final
                segment_kind = "exact_physical_endpoint"
            else:
                segment_target = np.asarray(q_hit, dtype=float)
                segment_kind = "exact_mesh_ray_crossing"

        segment_index = int(len(event["ray_segments"]))
        result, patch_levels, segment_failure = _local._segment_retry(
            original,
            self,
            current,
            kwargs,
            point,
            segment_target,
            direction,
            segment_index=segment_index,
            segment_count=-1,
            segment_kind=segment_kind,
        )

        if isinstance(result, _local.MeshOnlyRayHandoff):
            total_patch_levels += int(patch_levels)
            current = result.state
            event["mesh_only_handoffs"].append(
                {
                    "handoff_index": int(len(event["mesh_only_handoffs"])),
                    "transaction_iteration": int(transaction_iteration),
                    "physical_tip_point_m": point.tolist(),
                    "physical_endpoint_m": final.tolist(),
                    "patch_levels": int(patch_levels),
                    "physical_motion_m": 0.0,
                    "force_exact_ray_crossing_restart": True,
                    "refinement_record": copy.deepcopy(result.record),
                }
            )
            force_ray_crossing = True
            continue

        event["ray_segments"].append(
            {
                "segment_index": int(segment_index),
                "segment_kind": segment_kind,
                "p0_m": point.tolist(),
                "p1_m": np.asarray(segment_target, dtype=float).tolist(),
                "requested_length_m": float(np.linalg.norm(segment_target - point)),
                "patch_levels": int(patch_levels),
                "success": result is not None,
                "failure": segment_failure,
            }
        )
        if result is None:
            failure = str(segment_failure or "local_segment_retry_failed")
            break
        if float(result.moved) <= tol:
            failure = "zero_or_subtolerance_exact_ray_segment"
            break

        force_ray_crossing = False
        total_patch_levels += int(patch_levels)
        moved_total += float(result.moved)
        max_angle = max(max_angle, abs(float(result.angle_error_deg)))
        current = _state_from_result(result)
        point = _tip_point(self, front_id, segment_target)

        if segment_kind == "exact_physical_endpoint":
            break
    else:
        failure = f"exact_ray_transaction_iteration_limit:{max_iterations}"

    length_error = abs(moved_total - requested)
    length_tol = max(1.0e-12, 1.0e-6 * requested)
    endpoint_error = float(np.linalg.norm(point - final))
    success = failure is None and length_error <= length_tol and endpoint_error <= length_tol

    event.update(
        {
            "success": bool(success),
            "failure": failure,
            "segment_count": int(len(event["ray_segments"])),
            "mesh_only_handoff_count": int(len(event["mesh_only_handoffs"])),
            "moved_length_m": float(moved_total),
            "length_error_m": float(length_error),
            "endpoint_error_m": float(endpoint_error),
            "patch_refinement_levels_total": int(total_patch_levels),
        }
    )

    if not success:
        self._transaction_rollback(root_snap)
        _clear_transaction_audit(a0, w0, r0)
        _AUDIT["events"].append(copy.deepcopy(event))
        return None, str(failure or "exact_ray_transaction_mismatch"), event

    _mark_accepted_rows(
        a0=a0,
        physical_event_index=event_index,
        requested=requested,
        segment_count=len(event["ray_segments"]),
    )
    _mark_advance_log(
        self,
        log0=log0,
        physical_event_index=event_index,
        requested=requested,
        segment_count=len(event["ray_segments"]),
        patch_levels=total_patch_levels,
    )
    _AUDIT["events"].append(copy.deepcopy(event))
    return _cb.CrackAdvanceResult(
        current.mesh,
        current.boundary,
        current.damage,
        current.displacement,
        moved_total,
        True,
        angle_error_deg=max_angle,
        selected_edge_length=moved_total,
        reason="v10051837_event_resolved_exact_ray_local_cavity",
        elem_parent_map=current.parent_to_root,
    ), None, event


def event_resolved_strict_advance_v10051837(self: Any, *args, **kwargs):
    if args:
        raise TypeError("v10.0.5.18.3.7 advance requires keyword arguments")
    original = event_resolved_strict_advance_v10051837._original

    # Recursive calls produced internally by the raw exact-ray marcher are
    # implementation segments, not new physical-event transactions.
    if int(kwargs.get("_subdepth", 0) or 0) > 0:
        return _local._checked_call(original, self, kwargs)

    p0 = np.asarray(kwargs["p0"], dtype=float)
    p1 = np.asarray(kwargs["p1"], dtype=float)
    direction = np.asarray(kwargs["direction"], dtype=float)
    norm = float(np.linalg.norm(direction))
    requested = float(np.linalg.norm(p1 - p0))
    if norm <= 0.0 or requested <= 0.0:
        return _local._checked_call(original, self, kwargs)
    direction = direction / norm

    result, failure, event = _exact_ray_transaction(
        original, self, kwargs, p0, p1, direction
    )
    if result is not None:
        _v91856._AUDIT["quality_aware_retry_successes"].append(
            {
                "front_id": int(kwargs.get("front_id", -1)),
                "retry_kind": "event_resolved_exact_ray_local_cavity",
                "patch_levels": int(event["patch_refinement_levels_total"]),
                "ray_segment_count": int(event["segment_count"]),
                "mesh_only_handoff_count": int(event.get("mesh_only_handoff_count", 0)),
                "partitions": 1,
                "equal_length_partition_retry": False,
                "requested_length_m": requested,
            }
        )
        return result

    failure_row = {
        "front_id": int(kwargs.get("front_id", -1)),
        "requested_length_m": requested,
        "ray_segment_count": int(event.get("segment_count", 0)),
        "mesh_only_handoff_count": int(event.get("mesh_only_handoff_count", 0)),
        "patch_levels_attempted": _local._patch_levels(),
        "equal_length_partition_retry": False,
        "last_reason": str(failure),
    }
    _v91856._AUDIT["quality_aware_retry_failures"].append(failure_row)
    veto = _cb.CrackAdvanceResult(
        kwargs["mesh"],
        kwargs["boundary"],
        kwargs["damage"],
        kwargs["displacement"],
        0.0,
        False,
        reason=f"v10051837_no_feasible_local_refinement_path:{failure}",
    )
    return _v91856._record_or_raise(self, kwargs, veto)


class EventResolvedAuditedAdaptiveCZMBackendV10051837(_cb.AdaptiveCZMBackend):
    """Telemetry wrapper without the historical equal-partition backend."""

    name = "adaptive_czm_v10051837_event_resolved_audited"

    def advance(self, **kwargs):
        nested = int(kwargs.get("_subdepth", 0) or 0) > 0 or bool(
            kwargs.get("_partition_retry_active", False)
        )
        if nested:
            return super().advance(**kwargs)

        mesh = kwargs["mesh"]
        p0 = np.asarray(kwargs["p0"], dtype=float)
        p1 = np.asarray(kwargs["p1"], dtype=float)
        front_id = int(kwargs.get("front_id", 0))
        requested = float(np.linalg.norm(p1 - p0))
        committed = _audit36._tip_one_ring(self, mesh, p0, front_id)
        record: dict[str, Any] = {
            "physical_event_index": int(len(_audit36._AUDIT["events"])),
            "front_id": front_id,
            "committed_tip_point_m": p0.tolist(),
            "requested_endpoint_m": p1.tolist(),
            "requested_event_length_m": requested,
            "requested_direction": np.asarray(kwargs.get("direction", []), float).tolist(),
            "committed_tip_one_ring": committed,
            "committed_tip_h_over_requested_da": (
                None if not committed.get("available") or requested <= 0.0
                else float(committed["h_mean_m"]) / requested
            ),
            "target_parent": _audit36._target_metrics(self, mesh, p1),
            "ray_corridor": _audit36._ray_corridor(self, mesh, p0, p1),
        }
        result = super().advance(**kwargs)
        endpoint = p1
        if bool(getattr(result, "inserted", False)) and front_id in getattr(self, "tip_nodes", {}):
            endpoint = np.asarray(self.tip_nodes[front_id][2], dtype=float)
        record["backend_result"] = _audit36._candidate_quality(
            self, mesh, result, endpoint, front_id
        )
        record["event_resolved_retry"] = (
            copy.deepcopy(_AUDIT["events"][-1]) if _AUDIT["events"] else None
        )
        _audit36._AUDIT["events"].append(record)
        return result


def audit_payload() -> dict[str, Any]:
    successes = sum(bool(row.get("success")) for row in _AUDIT["events"])
    return {
        "schema": SCHEMA,
        "model_id": MODEL_ID,
        "quality_gate_inside_atomic_transaction": True,
        "exact_ray_crossing_decomposition": True,
        "mesh_only_refinement_handoff_restart": True,
        "mesh_only_handoff_consumes_physical_length": False,
        "equal_length_partition_retry": False,
        "exact_physical_event_endpoint_preserved": True,
        "exact_selected_crack_direction_preserved": True,
        "local_tip_radial_refinement": True,
        "adjacent_target_cavity_refinement": True,
        "triangle_quality_floor_relaxed": False,
        "child_area_ratio_floor_relaxed": False,
        "tip_h_over_da_role": "refinement_trigger_and_audit_only_not_event_veto",
        "constitutive_physics_changed": False,
        "event_count": int(len(_AUDIT["events"])),
        "success_count": int(successes),
        "failure_count": int(len(_AUDIT["events"]) - successes),
        "events": copy.deepcopy(_AUDIT["events"]),
        "candidate_vetoes": copy.deepcopy(
            _v91856._AUDIT.get("quality_aware_candidate_vetoes", [])
        ),
        "patch_refinements": copy.deepcopy(
            _v91856._AUDIT.get("quality_aware_patch_refinements", [])
        ),
        "retry_successes": copy.deepcopy(
            _v91856._AUDIT.get("quality_aware_retry_successes", [])
        ),
        "retry_failures": copy.deepcopy(
            _v91856._AUDIT.get("quality_aware_retry_failures", [])
        ),
    }


__all__ = [
    "EventResolvedAuditedAdaptiveCZMBackendV10051837",
    "MODEL_ID",
    "SCHEMA",
    "_exact_ray_transaction",
    "audit_payload",
    "configure_legacy_quality_wrapper",
    "event_resolved_strict_advance_v10051837",
    "reset_audit",
]
