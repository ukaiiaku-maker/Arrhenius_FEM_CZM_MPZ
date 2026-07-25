"""Long-extension peak analysis for the v9.13 autonomous DBTT model.

The functions in this module are deliberately independent of the solver. They
operate on event dictionaries written by ``RCurveResult.as_dict()`` and enforce
strict checkpoint semantics: a requested checkpoint is unavailable unless the
run actually reached it.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Iterable, Mapping, Sequence

import numpy as np


@dataclass(frozen=True)
class PeakMetrics:
    peak_temperature_K: float
    peak_temperature_quadratic_K: float
    peak_value: float
    peak_index: int
    peak_at_boundary: bool
    peak_rise: float
    post_peak_drop: float
    post_peak_minimum: float
    final_rebound: float
    peak_prominence: float
    n_finite: int

    def as_dict(self) -> dict[str, Any]:
        return {
            "peak_temperature_K": self.peak_temperature_K,
            "peak_temperature_quadratic_K": self.peak_temperature_quadratic_K,
            "peak_value": self.peak_value,
            "peak_index": self.peak_index,
            "peak_at_boundary": self.peak_at_boundary,
            "peak_rise": self.peak_rise,
            "post_peak_drop": self.post_peak_drop,
            "post_peak_minimum": self.post_peak_minimum,
            "final_rebound": self.final_rebound,
            "peak_prominence": self.peak_prominence,
            "n_finite": self.n_finite,
        }


def truthy(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return False
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return bool(value) and math.isfinite(float(value))
    return str(value).strip().lower() in {"1", "true", "yes", "y", "on"}


def achieved_extension_um(events: Sequence[Mapping[str, Any]]) -> float:
    if not events:
        return 0.0
    return 1.0e6 * float(events[-1]["cumulative_projected_extension_m"])


def checkpoint_from_events(
    events: Sequence[Mapping[str, Any]],
    extension_um: float,
    *,
    strict: bool = True,
    tolerance_um: float = 1.0e-9,
) -> float:
    """Return the first event at or beyond a checkpoint.

    When ``strict`` is true, return NaN if the event history did not reach the
    requested extension. This prevents a 52 micrometre loading-map endpoint
    from being silently reported as K100.
    """
    target_um = float(extension_um)
    if not math.isfinite(target_um) or target_um < 0.0:
        raise ValueError("checkpoint extension must be finite and nonnegative")
    if not events:
        return float("nan")
    target_m = target_um * 1.0e-6
    for event in events:
        if (
            float(event["cumulative_projected_extension_m"])
            + tolerance_um * 1.0e-6
            >= target_m
        ):
            return float(event["K_MPa_sqrt_m"])
    if strict:
        return float("nan")
    return float(events[-1]["K_MPa_sqrt_m"])


def checkpoint_reached(
    events: Sequence[Mapping[str, Any]],
    extension_um: float,
    *,
    tolerance_um: float = 1.0e-9,
) -> bool:
    return achieved_extension_um(events) + tolerance_um >= float(extension_um)


def _quadratic_peak_temperature(
    temperatures: np.ndarray,
    values: np.ndarray,
    peak_index: int,
) -> float:
    if peak_index <= 0 or peak_index >= len(values) - 1:
        return float("nan")
    x = temperatures[peak_index - 1 : peak_index + 2]
    y = values[peak_index - 1 : peak_index + 2]
    if np.any(~np.isfinite(x)) or np.any(~np.isfinite(y)):
        return float("nan")
    try:
        a, b, _ = np.polyfit(x, y, 2)
    except (TypeError, ValueError, np.linalg.LinAlgError):
        return float("nan")
    if not np.isfinite(a) or not np.isfinite(b) or a >= 0.0:
        return float("nan")
    vertex = -b / (2.0 * a)
    if vertex < x[0] or vertex > x[-1]:
        return float("nan")
    return float(vertex)


def peak_metrics(
    temperatures_K: Iterable[float],
    values: Iterable[float],
) -> PeakMetrics:
    temperature_array = np.asarray(tuple(temperatures_K), dtype=float)
    value_array = np.asarray(tuple(values), dtype=float)
    if temperature_array.shape != value_array.shape:
        raise ValueError("temperature and response arrays must have equal shape")
    finite = np.isfinite(temperature_array) & np.isfinite(value_array)
    temperature_array = temperature_array[finite]
    value_array = value_array[finite]
    if len(value_array) == 0:
        nan = float("nan")
        return PeakMetrics(nan, nan, nan, -1, True, nan, nan, nan, nan, nan, 0)
    order = np.argsort(temperature_array, kind="stable")
    temperature_array = temperature_array[order]
    value_array = value_array[order]
    peak_index = int(np.argmax(value_array))
    peak_temperature = float(temperature_array[peak_index])
    peak_value = float(value_array[peak_index])
    at_boundary = peak_index == 0 or peak_index == len(value_array) - 1
    before = value_array[:peak_index]
    after = value_array[peak_index + 1 :]
    pre_minimum = float(np.min(before)) if len(before) else float("nan")
    post_minimum = float(np.min(after)) if len(after) else float("nan")
    rise = peak_value - pre_minimum if np.isfinite(pre_minimum) else float("nan")
    drop = peak_value - post_minimum if np.isfinite(post_minimum) else float("nan")
    rebound = (
        float(value_array[-1]) - post_minimum
        if np.isfinite(post_minimum)
        else float("nan")
    )
    prominence = (
        min(rise, drop)
        if np.isfinite(rise) and np.isfinite(drop)
        else float("nan")
    )
    quadratic = _quadratic_peak_temperature(
        temperature_array,
        value_array,
        peak_index,
    )
    return PeakMetrics(
        peak_temperature_K=peak_temperature,
        peak_temperature_quadratic_K=quadratic,
        peak_value=peak_value,
        peak_index=peak_index,
        peak_at_boundary=at_boundary,
        peak_rise=rise,
        post_peak_drop=drop,
        post_peak_minimum=post_minimum,
        final_rebound=rebound,
        peak_prominence=prominence,
        n_finite=len(value_array),
    )


def choose_peak_temperature(metrics: PeakMetrics, estimator: str) -> float:
    mode = str(estimator).lower()
    if mode == "discrete":
        return float(metrics.peak_temperature_K)
    if mode == "quadratic":
        if np.isfinite(metrics.peak_temperature_quadratic_K):
            return float(metrics.peak_temperature_quadratic_K)
        return float(metrics.peak_temperature_K)
    raise ValueError("peak estimator must be 'discrete' or 'quadratic'")


def peak_drift_classification(
    temperatures_K: Iterable[float],
    *,
    stable_limit_K: float = 50.0,
    maximum_alignable_drift_K: float = 100.0,
) -> tuple[float, str]:
    values = np.asarray(tuple(temperatures_K), dtype=float)
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return float("nan"), "unresolved"
    drift = float(np.max(values) - np.min(values))
    if drift <= float(stable_limit_K):
        return drift, "stable"
    if drift <= float(maximum_alignable_drift_K):
        return drift, "moderate_drift"
    return drift, "extension_dependent"


def loading_map_coverage_um(payload: Mapping[str, Any]) -> float:
    projected = np.asarray(payload.get("projected_advances_m", ()), dtype=float)
    if (
        projected.size == 0
        or np.any(~np.isfinite(projected))
        or np.any(projected <= 0.0)
    ):
        raise ValueError("loading map projected advances must be positive and finite")
    return float(np.sum(projected) * 1.0e6)


__all__ = [
    "PeakMetrics",
    "achieved_extension_um",
    "checkpoint_from_events",
    "checkpoint_reached",
    "choose_peak_temperature",
    "loading_map_coverage_um",
    "peak_drift_classification",
    "peak_metrics",
    "truthy",
]
