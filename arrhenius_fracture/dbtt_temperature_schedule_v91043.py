"""Dynamic temperature schedules for the v9.10.4.3 narrow-DBTT campaign.

The coarse search is always evaluated on a broad 100 K grid.  Each promoted
candidate then carries its own adjacent coarse transition bracket.  Four points
are placed across that bracket while low- and high-shelf anchors are retained
from the coarse grid.  This keeps the factor-of-two shelf test meaningful while
resolving the transition more finely without assuming its absolute temperature.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import json
from typing import Iterable, Mapping, Any

import numpy as np


@dataclass(frozen=True)
class DynamicTemperatureSchedule:
    coarse_temperatures_K: tuple[float, ...]
    transition_low_K: float
    transition_high_K: float
    transition_temperatures_K: tuple[float, ...]
    low_anchor_temperatures_K: tuple[float, ...]
    high_anchor_temperatures_K: tuple[float, ...]
    evaluation_temperatures_K: tuple[float, ...]
    refinement_points: int
    shelf_anchor_count: int

    def to_columns(self) -> dict[str, Any]:
        return {
            "coarse_transition_low_T_K": float(self.transition_low_K),
            "coarse_transition_high_T_K": float(self.transition_high_K),
            "refinement_transition_temperatures_K": json.dumps(list(self.transition_temperatures_K)),
            "refinement_low_anchor_temperatures_K": json.dumps(list(self.low_anchor_temperatures_K)),
            "refinement_high_anchor_temperatures_K": json.dumps(list(self.high_anchor_temperatures_K)),
            "refinement_evaluation_temperatures_K": json.dumps(list(self.evaluation_temperatures_K)),
            "refinement_points": int(self.refinement_points),
            "shelf_anchor_count": int(self.shelf_anchor_count),
        }

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _spanning_anchors(values: np.ndarray, count: int) -> tuple[float, ...]:
    values = np.unique(np.asarray(values, dtype=float))
    if values.size == 0 or count <= 0:
        return ()
    if values.size <= count:
        return tuple(float(x) for x in values)
    indices = np.linspace(0, values.size - 1, count)
    indices = np.unique(np.rint(indices).astype(int))
    return tuple(float(values[i]) for i in indices)


def schedule_from_bracket(
    coarse_temperatures_K: Iterable[float],
    transition_low_K: float,
    transition_high_K: float,
    *,
    refinement_points: int = 4,
    shelf_anchor_count: int = 2,
) -> DynamicTemperatureSchedule:
    coarse = np.unique(np.asarray(list(coarse_temperatures_K), dtype=float))
    if coarse.size < 2:
        raise ValueError("at least two coarse temperatures are required")
    low = float(transition_low_K)
    high = float(transition_high_K)
    if not np.isfinite(low) or not np.isfinite(high) or high <= low:
        raise ValueError(f"invalid transition bracket: {low}, {high}")
    if refinement_points < 2:
        raise ValueError("refinement_points must be at least two")

    transition = tuple(float(x) for x in np.linspace(low, high, int(refinement_points)))
    low_anchors = _spanning_anchors(coarse[coarse < low], int(shelf_anchor_count))
    high_anchors = _spanning_anchors(coarse[coarse > high], int(shelf_anchor_count))
    evaluation = tuple(
        float(x)
        for x in np.unique(
            np.asarray(low_anchors + transition + high_anchors, dtype=float)
        )
    )
    if len(low_anchors) < 1 or len(high_anchors) < 1:
        raise ValueError(
            "the selected transition bracket must leave at least one coarse "
            "temperature on each shelf"
        )
    return DynamicTemperatureSchedule(
        coarse_temperatures_K=tuple(float(x) for x in coarse),
        transition_low_K=low,
        transition_high_K=high,
        transition_temperatures_K=transition,
        low_anchor_temperatures_K=low_anchors,
        high_anchor_temperatures_K=high_anchors,
        evaluation_temperatures_K=evaluation,
        refinement_points=int(refinement_points),
        shelf_anchor_count=int(shelf_anchor_count),
    )


def _parse_json_temperatures(value: Any) -> tuple[float, ...]:
    if value is None:
        return ()
    if isinstance(value, (list, tuple, np.ndarray)):
        return tuple(float(x) for x in value)
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return ()
    parsed = json.loads(text)
    return tuple(float(x) for x in parsed)


def schedule_from_candidate_row(
    row: Mapping[str, Any],
    coarse_temperatures_K: Iterable[float] = tuple(range(300, 1101, 100)),
    *,
    refinement_points: int = 4,
    shelf_anchor_count: int = 2,
) -> DynamicTemperatureSchedule:
    existing = _parse_json_temperatures(row.get("refinement_evaluation_temperatures_K"))
    low = row.get("coarse_transition_low_T_K", row.get("transition_transition_low_K", np.nan))
    high = row.get("coarse_transition_high_T_K", row.get("transition_transition_high_K", np.nan))
    schedule = schedule_from_bracket(
        coarse_temperatures_K,
        float(low),
        float(high),
        refinement_points=int(row.get("refinement_points", refinement_points)),
        shelf_anchor_count=int(row.get("shelf_anchor_count", shelf_anchor_count)),
    )
    if existing and not np.allclose(existing, schedule.evaluation_temperatures_K):
        raise ValueError(
            "stored candidate temperature schedule is inconsistent with its "
            "coarse transition bracket"
        )
    return schedule


__all__ = [
    "DynamicTemperatureSchedule",
    "schedule_from_bracket",
    "schedule_from_candidate_row",
]
