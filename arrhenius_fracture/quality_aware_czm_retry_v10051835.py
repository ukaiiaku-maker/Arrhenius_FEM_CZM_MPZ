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
    r0 = len(_v9185._RUNTIME["quality_vetoes"])
    call = dict(kwargs)
    call["_partition_retry_active"] = True
    result = _LEGACY_QUALITY_WRAPPER(self, **call)
    if not bool(getattr(result, "inserted", False)):
        rows = audit["quality_vetoes"][q0:]
        if rows:
            audit["quality_aware_candidate_vetoes"].extend(copy.deepcopy(rows))
            del audit["quality_vetoes"][q0:]
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
    snap = self._transaction_snapshot()
    audit = _v91856._AUDIT
    a0 = len(audit["accepted_events"])
    w0 = len(audit["resolution_warnings"])
    log0 = len(self.advance_log)
    current = state
    point = p0.copy()
    moved = 0.0
    max_angle = 0.0
    for index in range(partitions):
        target = p0 + requested * (index + 1) / partitions * direction
        result = _call_on_state(
            original, self, current, root_kwargs, point, target, direction
        )
        if not result.inserted:
            self._transaction_rollback(snap)
            del audit["accepted_events"][a0:]
            del audit["resolution_warnings"][w0:]
            return None, str(result.reason)
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
        return None, f"partition_length_mismatch:{moved:.12e}!={requested:.12e}"
    for row in self.advance_log[log0:]:
        row.update(
            v10051835_quality_aware_retry=True,
            v10051835_retry_kind="partition",
            v10051835_partition_count=partitions,
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
            reason=f"quality_aware_partition_{partitions}",
            elem_parent_map=current.parent_to_root,
        ),
        None,
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
    result = _call_on_state(original, self, root, kwargs, p0, p1, direction)
    if result.inserted:
        return result
    last_reason = str(result.reason)
    self._transaction_rollback(root_snap)
    state = root
    for level in range(1, _patch_levels() + 1):
        state, record = _refine_once(
            self,
            state,
            p0,
            p1,
            int(kwargs.get("front_id", 0)),
            level,
        )
        _v91856._AUDIT["quality_aware_patch_refinements"].append(
            copy.deepcopy(record)
        )
        if state is None:
            last_reason = str(record.get("reason", "patch_refinement_failed"))
            break
        result = _call_on_state(original, self, state, kwargs, p0, p1, direction)
        if result.inserted:
            for row in self.advance_log:
                if row.get("v91856_quality_gate_passed") and not row.get(
                    "v10051835_quality_aware_retry"
                ):
                    row.update(
                        v10051835_quality_aware_retry=True,
                        v10051835_retry_kind="shape_regular_patch",
                        v10051835_patch_refinement_levels=level,
                        v10051835_physical_event_requested_length_m=requested,
                    )
            _v91856._AUDIT["quality_aware_retry_successes"].append(
                {
                    "front_id": int(kwargs.get("front_id", -1)),
                    "retry_kind": "shape_regular_patch",
                    "patch_levels": level,
                    "partitions": 1,
                    "requested_length_m": requested,
                }
            )
            return result
        last_reason = str(result.reason)
        for partitions in _partition_counts():
            result, failure = _partition_retry(
                original,
                self,
                state,
                kwargs,
                p0,
                direction,
                requested,
                partitions,
            )
            if result is not None:
                _v91856._AUDIT["quality_aware_retry_successes"].append(
                    {
                        "front_id": int(kwargs.get("front_id", -1)),
                        "retry_kind": "partition_after_patch",
                        "patch_levels": level,
                        "partitions": partitions,
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
