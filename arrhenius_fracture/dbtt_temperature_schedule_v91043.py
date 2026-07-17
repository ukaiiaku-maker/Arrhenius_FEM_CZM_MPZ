"""Dynamic temperature schedules and fixed-bracket DBTT metrics for v9.10.4.3.

The coarse search is always evaluated on a broad 100 K grid. Each promoted
candidate then carries its own adjacent coarse transition bracket. Four points
are placed across that bracket while low- and high-shelf anchors are retained.
The refined metric scores the complete 100 K bracket rather than demanding that
75 percent of the rise occur in one of the new ~33 K subintervals.
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
        for x in np.unique(np.asarray(low_anchors + transition + high_anchors, dtype=float))
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
    stored_points = row.get("refinement_points", refinement_points)
    stored_anchors = row.get("shelf_anchor_count", shelf_anchor_count)
    if stored_points is None or not np.isfinite(float(stored_points)):
        stored_points = refinement_points
    if stored_anchors is None or not np.isfinite(float(stored_anchors)):
        stored_anchors = shelf_anchor_count
    schedule = schedule_from_bracket(
        coarse_temperatures_K,
        float(low),
        float(high),
        refinement_points=int(stored_points),
        shelf_anchor_count=int(stored_anchors),
    )
    if existing and not np.allclose(existing, schedule.evaluation_temperatures_K):
        raise ValueError(
            "stored candidate temperature schedule is inconsistent with its "
            "coarse transition bracket"
        )
    return schedule


def _indices_for_temperatures(T: np.ndarray, requested: Iterable[float]) -> np.ndarray:
    indices = []
    for value in requested:
        matches = np.flatnonzero(np.isclose(T, float(value), rtol=0.0, atol=1.0e-7))
        if matches.size != 1:
            raise ValueError(f"temperature {value} K is missing or duplicated")
        indices.append(int(matches[0]))
    return np.asarray(indices, dtype=int)


def _span_fraction(values: np.ndarray, center: float) -> float:
    if values.size <= 1:
        return 0.0
    return float((np.max(values) - np.min(values)) / max(abs(center), 1.0e-12))


def _crossing_temperature(T: np.ndarray, y: np.ndarray, level: float) -> float:
    if y[0] >= level:
        return float(T[0])
    for i in range(len(T) - 1):
        y0 = float(y[i])
        y1 = float(y[i + 1])
        if y0 < level <= y1 and y1 > y0:
            fraction = (level - y0) / (y1 - y0)
            return float(T[i] + fraction * (T[i + 1] - T[i]))
    return float("nan")


def fixed_bracket_transition_metrics(
    temperatures_K: Iterable[float],
    toughness: Iterable[float],
    schedule: DynamicTemperatureSchedule,
    *,
    plasticity_off_toughness: Iterable[float] | None = None,
    min_ratio: float = 2.0,
    robust_ratio_min: float = 1.8,
    max_low_span_fraction: float = 0.15,
    max_high_span_fraction: float = 0.20,
    min_bracket_concentration: float = 0.75,
    max_plasticity_off_ratio: float = 1.25,
    max_transition_width_K: float = 100.0,
    min_monotonic_fraction: float = 0.90,
) -> dict[str, Any]:
    """Score a candidate in its fixed coarse 100 K transition bracket."""
    T = np.asarray(list(temperatures_K), dtype=float)
    K = np.asarray(list(toughness), dtype=float)
    order = np.argsort(T)
    T = T[order]
    K = K[order]
    if not np.all(np.isfinite(K)):
        return {"valid": False, "loss": 1.0e12, "reason": "nonfinite_toughness"}

    low_i = _indices_for_temperatures(T, schedule.low_anchor_temperatures_K)
    high_i = _indices_for_temperatures(T, schedule.high_anchor_temperatures_K)
    trans_i = _indices_for_temperatures(T, schedule.transition_temperatures_K)
    low_values = K[low_i]
    high_values = K[high_i]
    trans_values = K[trans_i]
    trans_T = T[trans_i]

    KL = float(np.median(low_values))
    KH = float(np.median(high_values))
    delta = KH - KL
    if not np.isfinite(delta) or delta <= 0.0:
        return {"valid": False, "loss": 1.0e12, "reason": "nonpositive_shelf_rise"}

    ratio = KH / max(KL, 1.0e-12)
    robust_ratio = float(np.min(high_values) / max(np.max(low_values), 1.0e-12))
    bracket_rise = float(trans_values[-1] - trans_values[0])
    concentration = bracket_rise / max(delta, 1.0e-12)
    low_span = _span_fraction(low_values, KL)
    high_span = _span_fraction(high_values, KH)

    increments = np.diff(trans_values)
    total_variation = float(np.sum(np.abs(increments)))
    positive_variation = float(np.sum(np.maximum(increments, 0.0)))
    monotonic_fraction = 1.0 if total_variation <= 1.0e-12 else positive_variation / total_variation
    max_negative_fraction = float(np.max(np.maximum(-increments, 0.0)) / max(delta, 1.0e-12)) if increments.size else 0.0

    normalized = (trans_values - KL) / max(delta, 1.0e-12)
    T10 = _crossing_temperature(trans_T, normalized, 0.10)
    T90 = _crossing_temperature(trans_T, normalized, 0.90)
    width = float(T90 - T10) if np.isfinite(T10) and np.isfinite(T90) and T90 >= T10 else float("inf")

    off_ratio = 1.0
    if plasticity_off_toughness is not None:
        Koff = np.asarray(list(plasticity_off_toughness), dtype=float)[order]
        if not np.all(np.isfinite(Koff)):
            return {"valid": False, "loss": 1.0e12, "reason": "nonfinite_plasticity_off_toughness"}
        off_ratio = float(np.median(Koff[high_i]) / max(np.median(Koff[low_i]), 1.0e-12))

    penalties = {
        "ratio": max(min_ratio - ratio, 0.0) / 0.25,
        "robust_ratio": max(robust_ratio_min - robust_ratio, 0.0) / 0.25,
        "low_flatness": max(low_span - max_low_span_fraction, 0.0) / 0.05,
        "high_flatness": max(high_span - max_high_span_fraction, 0.0) / 0.05,
        "bracket_concentration": max(min_bracket_concentration - concentration, 0.0) / 0.10,
        "transition_width": max(width - max_transition_width_K, 0.0) / 25.0 if np.isfinite(width) else 20.0,
        "monotonicity": max(min_monotonic_fraction - monotonic_fraction, 0.0) / 0.10,
        "plasticity_off_ratio": max(off_ratio - max_plasticity_off_ratio, 0.0) / 0.10,
    }
    loss = float(sum(value * value for value in penalties.values()))
    return {
        "valid": True,
        "loss": loss,
        "split_index": -1,
        "transition_low_K": float(schedule.transition_low_K),
        "transition_high_K": float(schedule.transition_high_K),
        "low_shelf": KL,
        "high_shelf": KH,
        "shelf_ratio": ratio,
        "robust_shelf_ratio": robust_ratio,
        "main_jump": bracket_rise,
        "jump_concentration": concentration,
        "bracket_rise": bracket_rise,
        "low_span_fraction": low_span,
        "high_span_fraction": high_span,
        "secondary_jump_ratio": max_negative_fraction,
        "plasticity_off_ratio": off_ratio,
        "transition_T10_K": T10,
        "transition_T90_K": T90,
        "transition_width_K": width,
        "transition_monotonic_fraction": monotonic_fraction,
        "transition_max_negative_fraction": max_negative_fraction,
        "transition_temperatures_K": list(float(x) for x in trans_T),
        "transition_toughness": list(float(x) for x in trans_values),
        "penalties": penalties,
    }


__all__ = [
    "DynamicTemperatureSchedule",
    "fixed_bracket_transition_metrics",
    "schedule_from_bracket",
    "schedule_from_candidate_row",
]
