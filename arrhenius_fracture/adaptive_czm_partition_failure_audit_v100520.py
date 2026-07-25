"""Persist terminal adaptive-CZM partition failures after atomic rollback.

v10.0.5.19 correctly rolls back every rejected geometry trial, including its
per-trial audit rows.  That is required for atomicity, but it also erased the
raw backend reasons needed to diagnose a sustained-growth failure.  This
additive diagnostic wrapper records rejection metadata outside the mutable
geometry transaction and writes one compact physical-event summary only after
all uniform and recursive recovery paths have failed.

No geometry, quality threshold, constitutive state, hazard clock, loading law,
or stochastic event length is modified.
"""
from __future__ import annotations

from collections import Counter
from contextlib import contextmanager
import math
from typing import Any, Iterator

import numpy as np

from . import adaptive_czm_quality_subdivision_v100518 as _v100518

MODEL_ID = "adaptive_CZM_partition_failure_audit_v10_0_5_20"


def _segment_payload(kwargs: dict[str, Any]) -> dict[str, Any]:
    p0 = np.asarray(kwargs.get("p0", []), float).reshape(-1)
    p1 = np.asarray(kwargs.get("p1", []), float).reshape(-1)
    length = float(np.linalg.norm(p1 - p0)) if p0.size == p1.size and p0.size else 0.0
    return {
        "p0_m": p0.tolist(),
        "p1_m": p1.tolist(),
        "segment_length_m": length,
        "subdepth": int(kwargs.get("_subdepth", 0) or 0),
        "front_id": int(kwargs.get("front_id", -1)),
    }


def _bounded_examples(rows: list[dict[str, Any]], limit: int = 24) -> list[dict[str, Any]]:
    if len(rows) <= limit:
        return rows
    head = rows[: limit // 2]
    tail = rows[-(limit - len(head)) :]
    return head + tail


@contextmanager
def installed_partition_failure_audit_v100520() -> Iterator[None]:
    """Wrap the active v10.0.5.19 recovery function with nonmutating diagnostics."""
    active_try = _v100518._try_subdivisions

    def audited_try(
        self: Any,
        original: Any,
        args: tuple[Any, ...],
        kwargs: dict[str, Any],
        outer_snapshot: Any,
        audit_snapshot: dict[str, int],
        counter_snapshot: tuple[int, Any],
    ):
        failures: list[dict[str, Any]] = []

        def recording_original(inner_self: Any, *inner_args, **inner_kwargs):
            result = original(inner_self, *inner_args, **inner_kwargs)
            if not bool(getattr(result, "inserted", False)):
                row = {
                    "kind": "raw_backend_rejection",
                    "reason": str(getattr(result, "reason", "unknown")),
                    "angle_error_deg": float(getattr(result, "angle_error_deg", 0.0)),
                    **_segment_payload(inner_kwargs),
                }
                failures.append(row)
            return result

        result = active_try(
            self,
            recording_original,
            args,
            kwargs,
            outer_snapshot,
            audit_snapshot,
            counter_snapshot,
        )
        if result is not None:
            return result

        # active_try has already restored the authoritative pre-event transaction.
        # Persist one diagnostic row outside that transaction so it survives the
        # v9.18.5.6 audit writer.
        reasons = Counter(row["reason"] for row in failures)
        examples = _bounded_examples(failures)
        requested = _segment_payload(kwargs)
        summary = {
            "schema": MODEL_ID,
            "accepted": False,
            "issues": ["adaptive_quality_partition_exhausted"],
            "front_id": requested["front_id"],
            "requested_da_m": requested["segment_length_m"],
            "requested_p0_m": requested["p0_m"],
            "requested_p1_m": requested["p1_m"],
            "raw_backend_rejection_count": int(len(failures)),
            "raw_backend_reason_counts": dict(sorted(reasons.items())),
            "dominant_raw_backend_reason": (
                reasons.most_common(1)[0][0] if reasons else "no_raw_backend_reason_captured"
            ),
            "maximum_observed_subdepth": int(
                max((row["subdepth"] for row in failures), default=0)
            ),
            "minimum_rejected_segment_m": float(
                min((row["segment_length_m"] for row in failures), default=math.nan)
            ),
            "maximum_rejected_segment_m": float(
                max((row["segment_length_m"] for row in failures), default=math.nan)
            ),
            "rejection_examples": examples,
            "geometry_transaction_rolled_back": True,
            "physical_event_consumed": False,
            "triangle_quality_floor_relaxed": False,
            "child_area_ratio_floor_relaxed": False,
            "constitutive_physics_changed": False,
        }
        _v100518._v91856._AUDIT.setdefault("quality_vetoes", []).append(summary)
        _v100518._v9185._RUNTIME.setdefault("quality_vetoes", []).append(summary.copy())
        self._v100520_last_partition_failure = summary
        return None

    _v100518._try_subdivisions = audited_try
    try:
        yield
    finally:
        _v100518._try_subdivisions = active_try


__all__ = [
    "MODEL_ID",
    "installed_partition_failure_audit_v100520",
]
