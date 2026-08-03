"""Compatibility layer for validated v10.0.5.18.3.5 local patch helpers.

The v10.0.5.18.3.7 controller imports these transaction primitives but does not
invoke the historical equal-length 2/4/8 partition policy.
"""
from __future__ import annotations

import copy
import os
from typing import Any, Callable

import numpy as np

from . import mode_i_first_passage_v9_18_5 as _v9185
from . import mode_i_first_passage_v9_18_5_6 as _v91856
from .quality_aware_czm_patch_v10051835 import (
    State as _State,
    compose_parent_map as _compose_parent_map,
    refine_once as _refine_once,
)

_LEGACY_QUALITY_WRAPPER: Callable | None = None
_State = _State


class MeshOnlyRayHandoff:
    """Accepted mesh refinement requiring exact-ray routing before tip motion.

    The backend transaction has been restored to the current physical tip, but
    ``state`` contains a conforming numerical refinement that must be retained.
    The outer event-resolved controller must recompute the next exact mesh-ray
    crossing on this refined mesh.  This object is deliberately not a crack
    advance result and carries no physical motion.
    """

    def __init__(self, state: _State, record: dict[str, Any]):
        self.state = state
        self.record = copy.deepcopy(record)


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
        value = int(os.environ.get("ARRHENIUS_EVENT_RESOLVED_PATCH_LEVELS", os.environ.get("ARRHENIUS_QUALITY_AWARE_PATCH_LEVELS", "12")))
    except ValueError:
        value = 12
    return max(value, 1)


def _reset_veto_counter(self: Any) -> None:
    self._v91856_consecutive_geometry_vetoes = 0
    self._v91856_last_veto_reason = None


def _checked_call(original: Callable, self: Any, kwargs: dict[str, Any]):
    if _LEGACY_QUALITY_WRAPPER is None:
        raise RuntimeError("quality wrapper not configured")
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
    segment_index: int,
    segment_count: int,
    segment_requested_length_m: float,
    segment_kind: str,
) -> None:
    row = copy.deepcopy(record)
    row.update(
        physical_event_partition_count=1,
        physical_event_segment_index=int(segment_index),
        physical_event_segment_count=int(segment_count),
        segment_requested_length_m=float(segment_requested_length_m),
        segment_kind=str(segment_kind),
        equal_length_partition_retry=False,
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
    segment_index: int = 0,
    segment_count: int = 1,
    segment_kind: str = "exact_physical_endpoint",
):
    snap = self._transaction_snapshot()
    requested = float(np.linalg.norm(np.asarray(p1) - np.asarray(p0)))
    result = _call_on_state(original, self, state, root_kwargs, p0, p1, direction)
    if result.inserted:
        return result, 0, None

    last_reason = str(result.reason)
    self._transaction_rollback(snap)
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
            segment_index=segment_index,
            segment_count=segment_count,
            segment_requested_length_m=requested,
            segment_kind=segment_kind,
        )
        if next_state is None:
            self._transaction_rollback(snap)
            return None, level - 1, str(record.get("reason", "patch_refinement_failed"))
        candidate_state = next_state

        # Some accepted refinements intentionally move the target out of the
        # directly reachable tip cavity.  Retrying the full endpoint from the
        # unchanged tip after such a refinement is topologically incorrect: the
        # refined mesh now contains an intervening exact ray crossing.  Preserve
        # the numerical refinement, restore the backend transaction, and hand
        # control to the outer exact-ray marcher before any physical motion.
        if bool(record.get("exact_ray_march_required_after_refinement", False)):
            self._transaction_rollback(snap)
            return MeshOnlyRayHandoff(candidate_state, record), level, None

        result = _call_on_state(original, self, candidate_state, root_kwargs, p0, p1, direction)
        if result.inserted:
            return result, level, None
        last_reason = str(result.reason)
        self._transaction_rollback(snap)

    self._transaction_rollback(snap)
    return None, _patch_levels(), last_reason


__all__ = [
    "MeshOnlyRayHandoff",
    "_State",
    "_call_on_state",
    "_checked_call",
    "_patch_levels",
    "_refine_once",
    "_segment_retry",
    "configure_legacy_quality_wrapper",
    "reset_audit",
]
