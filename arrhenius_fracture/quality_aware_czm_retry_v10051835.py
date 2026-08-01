"""Quality-aware adaptive-CZM transaction retry for v10.0.5.18.3.5."""
from __future__ import annotations

import copy
import os
from typing import Any, Callable

import numpy as np

from . import crack_backend as _cb
from . import mode_i_first_passage_v9_18_5 as _v9185
from . import mode_i_first_passage_v9_18_5_6 as _v91856
from .quality_aware_czm_patch_v10051835 import (
    State as _State,
    compose_parent_map as _compose_parent_map,
    floor_issues as _floor_issues,
    refine_once as _refine_once,
)


MODEL_ID = "quality_aware_adaptive_CZM_retry_v10_0_5_18_3_5"
SCHEMA = "v10.0.5.18.3.5_quality_aware_local_patch_transaction_retry"
_LEGACY_QUALITY_WRAPPER: Callable | None = None
_State = _State


def configure_legacy_quality_wrapper(wrapper: Callable) -> None:
    global _LEGACY_QUALITY_WRAPPER
    _LEGACY_QUALITY_WRAPPER = wrapper


def reset_audit() -> None:
    audit = _v91856._AUDIT
    audit["accepted_events"] = []
    audit["resolution_warnings"] = []
    audit["quality_vetoes"] = []
    audit["consecutive_veto_abort"] = None
    audit["quality_aware_candidate_vetoes"] = []
    audit["quality_aware_patch_refinements"] = []
    audit["quality_aware_retry_successes"] = []
    audit["quality_aware_retry_failures"] = []


def _patch_levels() -> int:
    try:
        value = int(os.environ.get("ARRHENIUS_QUALITY_AWARE_PATCH_LEVELS", "8"))
    except ValueError:
        value = 8
    return max(value, 1)


def _partition_counts() -> tuple[int, ...]:
    raw = os.environ.get("ARRHENIUS_QUALITY_AWARE_PARTITIONS", "2,4,8")
    values: list[int] = []
    for token in str(raw).replace(" ", ",").split(","):
        try:
            value = int(token)
        except ValueError:
            continue
        if value >= 2 and value not in values:
            values.append(value)
    return tuple(values or (2, 4, 8))


def _reset_veto_counter(self: Any) -> None:
    self._v91856_consecutive_geometry_vetoes = 0
    self._v91856_last_veto_reason = None


def _checked_call(original: Callable, self: Any, kwargs: dict[str, Any]):
    if _LEGACY_QUALITY_WRAPPER is None:
        raise RuntimeError("v10.0.5.18.3.5 quality wrapper not configured")
    _LEGACY_QUALITY_WRAPPER._original = original
    audit = _v91856._AUDIT
    q0 = len(audit["quality_vetoes"])
    a0 = len(audit["accepted_events"])
    w0 = len(audit["resolution_warnings"])
    r0 = len(_v9185._RUNTIME["quality_vetoes"])
    call = dict(kwargs)
    call["_partition_retry_active"] = True
    result = _LEGACY_QUALITY_WRAPPER(self, **call)
    if not bool(getattr(result, "inserted", False)):
        rows = audit["quality_vetoes"][q0:]
        if rows:
            audit["quality_aware_candidate_vetoes"].extend(copy.deepcopy(rows))
            del audit["quality_vetoes"][q0:]
        # A failed raw adaptive-CZM transaction can internally march across
        # several ray crossings before its final candidate is vetoed. Those
        # intermediate quality records are trial state, not committed events.
        del audit["accepted_events"][a0:]
        del audit["resolution_warnings"][w0:]
        del _v9185._RUNTIME["quality_vetoes"][r0:]
        _reset_veto_counter(self)
    return result


def _call_on_state(
    original: Callable,
    self: Any,
    state: _State,
    root_kwargs: dict[str, Any],
    p0: np.ndarray,
    p1: np.ndarray,
    direction: np.ndarray,
):
    call = dict(root_kwargs)
    call.update(
        mesh=state.mesh,
        boundary=state.boundary,
        damage=state.damage,
        displacement=state.displacement,
        p0=p0,
        p1=p1,
        direction=direction,
    )
    result = _checked_call(original, self, call)
    if result.inserted:
        result.elem_parent_map = _compose_parent_map(state, result)
    return result


def _append_patch_record(
    record: dict[str, Any],
    *,
    partitions: int,
    segment_index: int,
    segment_count: int,
    segment_requested_length_m: float,
) -> None:
    row = copy.deepcopy(record)
    row.update(
        physical_event_partition_count=int(partitions),
        physical_event_segment_index=int(segment_index),
        physical_event_segment_count=int(segment_count),
        segment_requested_length_m=float(segment_requested_length_m),
    )
    _v91856._AUDIT["quality_aware_patch_refinements"].append(row)


def _segment_retry(
    original: Callable,
    self: Any,
    state: _State,
    root_kwargs: dict[str, Any],
    p0: np.ndarray,
    p1: np.ndarray,
    direction: np.ndarray,
    *,
    partitions: int,
    segment_index: int,
    segment_count: int,
):
    """Solve one collinear segment with endpoint-local patch retry.

    This function never partitions again. The physical-event partition driver
    calls it for each segment, so a long event can traverse several tip
    triangles while still applying the same quality-aware patch logic at the
    precise subsegment that needs refinement.
    """
    segment_snap = self._transaction_snapshot()
    requested = float(np.linalg.norm(np.asarray(p1) - np.asarray(p0)))
    result = _call_on_state(
        original, self, state, root_kwargs, p0, p1, direction
    )
    if result.inserted:
        return result, 0, None

    last_reason = str(result.reason)
    self._transaction_rollback(segment_snap)
    candidate_state = state
    for level in range(1, _patch_levels() + 1):
        next_state, record = _refine_once(
            self,
            candidate_state,
            p0,
            p1,
            int(root_kwargs.get("front_id", 0)),
            level,
        )
        _append_patch_record(
            record,
            partitions=partitions,
            segment_index=segment_index,
            segment_count=segment_count,
            segment_requested_length_m=requested,
        )
        if next_state is None:
            self._transaction_rollback(segment_snap)
            return None, level - 1, str(
                record.get("reason", "patch_refinement_failed")
            )

        candidate_state = next_state
        result = _call_on_state(
            original,
            self,
            candidate_state,
            root_kwargs,
            p0,
            p1,
            direction,
        )
        if result.inserted:
            return result, level, None

        last_reason = str(result.reason)
        self._transaction_rollback(segment_snap)

    self._transaction_rollback(segment_snap)
    return None, _patch_levels(), last_reason


def _partition_retry(
    original: Callable,
    self: Any,
    state: _State,
    root_kwargs: dict[str, Any],
    p0: np.ndarray,
    direction: np.ndarray,
    requested: float,
    partitions: int,
):
    """Retry one physical event as atomic collinear geometry segments.

    Each segment receives its own endpoint-centered patch retry. Thus a distant
    event endpoint need not lie inside the initial tip one-ring.
    """
    snap = self._transaction_snapshot()
    audit = _v91856._AUDIT
    a0 = len(audit["accepted_events"])
    w0 = len(audit["resolution_warnings"])
    log0 = len(self.advance_log)
    current = state
    point = p0.copy()
    moved = 0.0
    max_angle = 0.0
    max_patch_levels = 0
    patched_segments = 0

    for index in range(partitions):
        target = p0 + requested * (index + 1) / partitions * direction
        result, patch_levels, failure = _segment_retry(
            original,
            self,
            current,
            root_kwargs,
            point,
            target,
            direction,
            partitions=partitions,
            segment_index=index,
            segment_count=partitions,
        )
        if result is None:
            self._transaction_rollback(snap)
            del audit["accepted_events"][a0:]
            del audit["resolution_warnings"][w0:]
            return None, (
                f"partition_{partitions}_segment_{index}:"
                f"{failure or 'segment_retry_failed'}"
            ), {
                "max_patch_levels": max_patch_levels,
                "patched_segments": patched_segments,
                "failed_segment_index": index,
            }

        if patch_levels > 0:
            patched_segments += 1
            max_patch_levels = max(max_patch_levels, patch_levels)

        current = _State(
            result.mesh,
            result.boundary,
            np.asarray(result.damage),
            np.asarray(result.displacement),
            np.asarray(result.elem_parent_map, dtype=int),
        )
        moved += float(result.moved)
        max_angle = max(max_angle, abs(float(result.angle_error_deg)))
        tips = getattr(self, "tip_nodes", {})
        front_id = int(root_kwargs.get("front_id", 0))
        point = (
            np.asarray(tips[front_id][2], dtype=float).copy()
            if front_id in tips
            else target.copy()
        )

    if abs(moved - requested) > max(1.0e-12, 1.0e-6 * requested):
        self._transaction_rollback(snap)
        del audit["accepted_events"][a0:]
        del audit["resolution_warnings"][w0:]
        return None, (
            f"partition_length_mismatch:{moved:.12e}!={requested:.12e}"
        ), {
            "max_patch_levels": max_patch_levels,
            "patched_segments": patched_segments,
            "failed_segment_index": None,
        }

    retry_kind = (
        "partition_with_segment_patch"
        if patched_segments
        else "partition"
    )
    for row in self.advance_log[log0:]:
        row.update(
            v10051835_quality_aware_retry=True,
            v10051835_retry_kind=retry_kind,
            v10051835_partition_count=partitions,
            v10051835_segment_patch_count=patched_segments,
            v10051835_max_segment_patch_levels=max_patch_levels,
            v10051835_physical_event_requested_length_m=requested,
        )
    return (
        _cb.CrackAdvanceResult(
            current.mesh,
            current.boundary,
            current.damage,
            current.displacement,
            moved,
            True,
            angle_error_deg=max_angle,
            selected_edge_length=moved,
            reason=f"quality_aware_{retry_kind}_{partitions}",
            elem_parent_map=current.parent_to_root,
        ),
        None,
        {
            "max_patch_levels": max_patch_levels,
            "patched_segments": patched_segments,
            "failed_segment_index": None,
        },
    )


def quality_aware_strict_advance_v10051835(self: Any, *args, **kwargs):
    if args:
        raise TypeError("v10.0.5.18.3.5 advance requires keyword arguments")
    original = quality_aware_strict_advance_v10051835._original
    p0 = np.asarray(kwargs["p0"], dtype=float)
    p1 = np.asarray(kwargs["p1"], dtype=float)
    direction = np.asarray(kwargs["direction"], dtype=float)
    norm = float(np.linalg.norm(direction))
    requested = float(np.linalg.norm(p1 - p0))
    if norm <= 0.0 or requested <= 0.0:
        return _checked_call(original, self, kwargs)
    direction /= norm

    root_snap = self._transaction_snapshot()
    root = _State(
        kwargs["mesh"],
        kwargs["boundary"],
        np.asarray(kwargs["damage"]),
        np.asarray(kwargs["displacement"]),
        np.arange(int(kwargs["mesh"].ne), dtype=int),
    )
    log0 = len(self.advance_log)
    result, patch_levels, last_reason = _segment_retry(
        original,
        self,
        root,
        kwargs,
        p0,
        p1,
        direction,
        partitions=1,
        segment_index=0,
        segment_count=1,
    )
    if result is not None:
        if patch_levels > 0:
            for row in self.advance_log[log0:]:
                row.update(
                    v10051835_quality_aware_retry=True,
                    v10051835_retry_kind="shape_regular_patch",
                    v10051835_patch_refinement_levels=patch_levels,
                    v10051835_physical_event_requested_length_m=requested,
                )
            _v91856._AUDIT["quality_aware_retry_successes"].append(
                {
                    "front_id": int(kwargs.get("front_id", -1)),
                    "retry_kind": "shape_regular_patch",
                    "patch_levels": patch_levels,
                    "partitions": 1,
                    "patched_segments": 1,
                    "requested_length_m": requested,
                }
            )
        return result

    self._transaction_rollback(root_snap)
    for partitions in _partition_counts():
        result, failure, summary = _partition_retry(
            original,
            self,
            root,
            kwargs,
            p0,
            direction,
            requested,
            partitions,
        )
        if result is not None:
            retry_kind = (
                "partition_with_segment_patch"
                if summary["patched_segments"]
                else "partition"
            )
            _v91856._AUDIT["quality_aware_retry_successes"].append(
                {
                    "front_id": int(kwargs.get("front_id", -1)),
                    "retry_kind": retry_kind,
                    "patch_levels": int(summary["max_patch_levels"]),
                    "partitions": partitions,
                    "patched_segments": int(summary["patched_segments"]),
                    "requested_length_m": requested,
                }
            )
            return result
        last_reason = str(failure)

    self._transaction_rollback(root_snap)
    failure = {
        "front_id": int(kwargs.get("front_id", -1)),
        "requested_length_m": requested,
        "patch_levels_attempted": _patch_levels(),
        "partition_counts_attempted": list(_partition_counts()),
        "last_reason": last_reason,
    }
    _v91856._AUDIT["quality_aware_retry_failures"].append(failure)
    veto = _cb.CrackAdvanceResult(
        kwargs["mesh"],
        kwargs["boundary"],
        kwargs["damage"],
        kwargs["displacement"],
        0.0,
        False,
        reason=f"v10051835_quality_aware_retry_exhausted:{last_reason}",
    )
    return _v91856._record_or_raise(self, kwargs, veto)


def audit_payload() -> dict[str, Any]:
    audit = _v91856._AUDIT
    return {
        "schema": SCHEMA,
        "model_id": MODEL_ID,
        "quality_gate_inside_retry_transaction": True,
        "shape_regular_tip_patch_midpoint_bisection": True,
        "exact_physical_event_endpoint_retried": True,
        "collinear_partition_retry_after_patch_refinement": True,
        "long_event_partition_segments_receive_patch_retry": True,
        "failed_trial_subsegments_removed_from_accepted_event_audit": True,
        "triangle_quality_floor_relaxed": False,
        "child_area_ratio_floor_relaxed": False,
        "immediate_parent_child_area_ratio_enforced": True,
        "cumulative_area_ratio_is_audit_only": True,
        "tip_h_over_da_role": "audit_warning_only",
        "constitutive_physics_changed": False,
        "candidate_vetoes": copy.deepcopy(
            audit.get("quality_aware_candidate_vetoes", [])
        ),
        "patch_refinements": copy.deepcopy(
            audit.get("quality_aware_patch_refinements", [])
        ),
        "retry_successes": copy.deepcopy(
            audit.get("quality_aware_retry_successes", [])
        ),
        "retry_failures": copy.deepcopy(
            audit.get("quality_aware_retry_failures", [])
        ),
    }


__all__ = [
    "MODEL_ID",
    "SCHEMA",
    "audit_payload",
    "configure_legacy_quality_wrapper",
    "quality_aware_strict_advance_v10051835",
    "reset_audit",
]
