"""Adaptive quality-safe partitioning for arbitrary stochastic CZM events.

The v10.0.5.18 repair first tries a small set of uniform collinear
subdivisions.  A particular endpoint can still be incompatible with every
count in that finite list even though a nonuniform collinear partition is
admissible.  This module adds a bounded recursive bisection fallback that
splits only a rejected segment, preserves the exact physical endpoint and
length, and keeps the complete physical event atomic.

No triangle-quality or child-area floor is relaxed.  No constitutive, hazard,
loading, material, or stochastic event-length law is changed.
"""
from __future__ import annotations

from contextlib import contextmanager
import math
import os
from typing import Any, Iterator

import numpy as np

from . import crack_backend as _cb
from . import adaptive_czm_quality_subdivision_v100518 as _v100518

MODEL_ID = "adaptive_CZM_recursive_quality_partition_v10_0_5_19"


def _int_env(name: str, default: int, minimum: int, maximum: int) -> int:
    try:
        value = int(os.environ.get(name, str(default)))
    except (TypeError, ValueError):
        value = int(default)
    return min(max(value, minimum), maximum)


def _float_env(name: str, default: float, minimum: float) -> float:
    try:
        value = float(os.environ.get(name, str(default)))
    except (TypeError, ValueError):
        value = float(default)
    if not math.isfinite(value):
        value = float(default)
    return max(value, minimum)


def _max_depth() -> int:
    return _int_env("ARRHENIUS_QUALITY_PARTITION_MAX_DEPTH", 9, 1, 14)


def _min_segment_m() -> float:
    return _float_env("ARRHENIUS_QUALITY_PARTITION_MIN_SEGMENT_M", 2.0e-8, 1.0e-12)


def _parent_map_or_identity(result: Any, input_ne: int) -> np.ndarray:
    parent = getattr(result, "elem_parent_map", None)
    if parent is None:
        if int(result.mesh.ne) != int(input_ne):
            raise RuntimeError(
                "adaptive quality partition cannot compose a missing parent map "
                "across a topology-changing segment"
            )
        return np.arange(int(input_ne), dtype=int)
    parent = np.asarray(parent, dtype=int)
    if np.any(parent < 0) or np.any(parent >= int(input_ne)):
        raise RuntimeError("adaptive quality partition received an invalid parent map")
    return parent


def _rollback_attempt(
    self: Any,
    transaction_snapshot: Any,
    audit_snapshot: dict[str, int],
    counter_snapshot: tuple[int, Any],
) -> None:
    self._transaction_rollback(transaction_snapshot)
    _v100518._truncate_audits(audit_snapshot)
    _v100518._restore_counter(self, counter_snapshot)


def _attempt_interval(
    self: Any,
    original: Any,
    base_kwargs: dict[str, Any],
    event_p0: np.ndarray,
    direction: np.ndarray,
    event_length: float,
    start_fraction: float,
    end_fraction: float,
    state: dict[str, Any],
    depth: int,
    stats: dict[str, int],
):
    transaction_snapshot = self._transaction_snapshot()
    audit_snapshot = _v100518._audit_lengths()
    counter_snapshot = _v100518._counter_snapshot(self)

    segment_p0 = event_p0 + direction * (event_length * start_fraction)
    segment_p1 = event_p0 + direction * (event_length * end_fraction)
    segment_length = float(np.linalg.norm(segment_p1 - segment_p0))
    stats["attempts"] += 1
    stats["maximum_depth"] = max(stats["maximum_depth"], int(depth))

    local = dict(base_kwargs)
    local.update(
        {
            "mesh": state["mesh"],
            "boundary": state["boundary"],
            "damage": state["damage"],
            "displacement": state["displacement"],
            "p0": segment_p0,
            "p1": segment_p1,
            "direction": direction,
            "_subdepth": 1,
        }
    )
    log_start = len(getattr(self, "advance_log", []) or [])
    input_ne = int(state["mesh"].ne)
    result = original(self, **local)
    if bool(getattr(result, "inserted", False)):
        issues, row, warning = _v100518._assess(
            self, state["mesh"], result, local, log_start
        )
        row.update(
            {
                "adaptive_quality_partition_active": True,
                "adaptive_quality_partition_depth": int(depth),
                "adaptive_quality_partition_start_fraction": float(start_fraction),
                "adaptive_quality_partition_end_fraction": float(end_fraction),
                "physical_event_total_length_m": float(event_length),
            }
        )
        if not issues:
            _v100518._accept(self, result, row, warning)
            stats["accepted_leaves"] += 1
            return result

    _rollback_attempt(self, transaction_snapshot, audit_snapshot, counter_snapshot)

    if depth >= _max_depth() or segment_length <= _min_segment_m():
        return None

    midpoint = 0.5 * (start_fraction + end_fraction)
    first = _attempt_interval(
        self,
        original,
        base_kwargs,
        event_p0,
        direction,
        event_length,
        start_fraction,
        midpoint,
        state,
        depth + 1,
        stats,
    )
    if first is None:
        _rollback_attempt(self, transaction_snapshot, audit_snapshot, counter_snapshot)
        return None

    first_state = {
        "mesh": first.mesh,
        "boundary": first.boundary,
        "damage": first.damage,
        "displacement": first.displacement,
    }
    second = _attempt_interval(
        self,
        original,
        base_kwargs,
        event_p0,
        direction,
        event_length,
        midpoint,
        end_fraction,
        first_state,
        depth + 1,
        stats,
    )
    if second is None:
        _rollback_attempt(self, transaction_snapshot, audit_snapshot, counter_snapshot)
        return None

    first_map = _parent_map_or_identity(first, input_ne)
    combined_map = _v100518._compose_parent_map(
        first_map,
        second,
        int(first.mesh.ne),
    )
    return _cb.CrackAdvanceResult(
        mesh=second.mesh,
        boundary=second.boundary,
        damage=second.damage,
        displacement=second.displacement,
        moved=float(first.moved) + float(second.moved),
        inserted=True,
        angle_error_deg=0.0,
        selected_edge_length=float(first.moved) + float(second.moved),
        reason="ok_adaptive_quality_partition",
        elem_parent_map=combined_map,
    )


def _try_adaptive_partition(
    self: Any,
    original: Any,
    args: tuple[Any, ...],
    kwargs: dict[str, Any],
    outer_snapshot: Any,
    audit_snapshot: dict[str, int],
    counter_snapshot: tuple[int, Any],
):
    if args:
        return None
    event_p0 = np.asarray(kwargs.get("p0"), float).reshape(2)
    event_p1 = np.asarray(kwargs.get("p1"), float).reshape(2)
    direction = np.asarray(kwargs.get("direction"), float).reshape(2)
    event_length = float(np.linalg.norm(event_p1 - event_p0))
    norm = float(np.linalg.norm(direction))
    if not math.isfinite(event_length) or event_length <= 0.0 or norm <= 1.0e-30:
        return None
    direction = direction / norm

    _rollback_attempt(self, outer_snapshot, audit_snapshot, counter_snapshot)
    stats = {"attempts": 0, "accepted_leaves": 0, "maximum_depth": 0}
    state = {
        "mesh": kwargs["mesh"],
        "boundary": kwargs["boundary"],
        "damage": kwargs["damage"],
        "displacement": kwargs["displacement"],
    }
    result = _attempt_interval(
        self,
        original,
        kwargs,
        event_p0,
        direction,
        event_length,
        0.0,
        1.0,
        state,
        0,
        stats,
    )
    if result is None:
        _rollback_attempt(self, outer_snapshot, audit_snapshot, counter_snapshot)
        return None

    if not math.isclose(
        float(result.moved),
        event_length,
        rel_tol=1.0e-8,
        abs_tol=5.0e-10,
    ):
        _rollback_attempt(self, outer_snapshot, audit_snapshot, counter_snapshot)
        raise RuntimeError(
            "adaptive quality partition changed the physical event length: "
            f"requested={event_length:.9e} m moved={float(result.moved):.9e} m"
        )

    if getattr(self, "advance_log", None):
        self.advance_log[-1].update(
            {
                "v100519_adaptive_quality_partition_recovery": True,
                "v100519_partition_attempts": int(stats["attempts"]),
                "v100519_partition_leaf_count": int(stats["accepted_leaves"]),
                "v100519_partition_maximum_depth": int(stats["maximum_depth"]),
                "v100519_physical_event_total_length_m": float(event_length),
            }
        )
    result.reason = (
        "ok_adaptive_quality_partition_"
        f"{int(stats['accepted_leaves'])}_leaves_depth_{int(stats['maximum_depth'])}"
    )
    self._v91856_consecutive_geometry_vetoes = 0
    self._v91856_last_veto_reason = None
    return result


@contextmanager
def installed_adaptive_quality_partition_v100519() -> Iterator[None]:
    uniform_try = _v100518._try_subdivisions

    def combined_try(
        self: Any,
        original: Any,
        args: tuple[Any, ...],
        kwargs: dict[str, Any],
        outer_snapshot: Any,
        audit_snapshot: dict[str, int],
        counter_snapshot: tuple[int, Any],
    ):
        result = uniform_try(
            self,
            original,
            args,
            kwargs,
            outer_snapshot,
            audit_snapshot,
            counter_snapshot,
        )
        if result is not None:
            return result
        return _try_adaptive_partition(
            self,
            original,
            args,
            kwargs,
            outer_snapshot,
            audit_snapshot,
            counter_snapshot,
        )

    _v100518._try_subdivisions = combined_try
    try:
        yield
    finally:
        _v100518._try_subdivisions = uniform_try


__all__ = [
    "MODEL_ID",
    "installed_adaptive_quality_partition_v100519",
]
