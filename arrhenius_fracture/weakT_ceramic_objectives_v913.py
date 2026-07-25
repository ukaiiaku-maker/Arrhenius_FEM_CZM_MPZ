"""Class-specific response metrics for v9.13 weak-T and ceramic searches.

The metrics are intentionally based on both initiation resistance and developed
resistance.  This prevents a short-distance, high-temperature-only screen from
promoting candidates that hide a low-temperature transition or a large long-range
R-curve.
"""
from __future__ import annotations

import math
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

from .emergent_gnd_contract_v913 import ACTIVE_CANDIDATE_PARAMETER_FIELDS


def _array(values: Sequence[float]) -> np.ndarray:
    return np.asarray([float(value) for value in values], dtype=float)


def _median(values: np.ndarray) -> float:
    finite = values[np.isfinite(values)]
    return float(np.median(finite)) if finite.size else float("nan")


def local_peak_prominence(values: Sequence[float]) -> float:
    curve = _array(values)
    if curve.size < 3 or not np.isfinite(curve).all():
        return float("nan")
    prominence = 0.0
    for index in range(1, curve.size - 1):
        if curve[index] > curve[index - 1] and curve[index] > curve[index + 1]:
            prominence = max(
                prominence,
                min(curve[index] - curve[index - 1], curve[index] - curve[index + 1]),
            )
    return float(prominence)


def response_metrics(
    temperatures_K: Sequence[float],
    K_initial_MPa_sqrt_m: Sequence[float],
    K_developed_MPa_sqrt_m: Sequence[float],
) -> dict[str, float | bool]:
    temperatures = _array(temperatures_K)
    initial = _array(K_initial_MPa_sqrt_m)
    developed = _array(K_developed_MPa_sqrt_m)
    complete = bool(
        temperatures.size >= 5
        and temperatures.size == initial.size == developed.size
        and np.isfinite(temperatures).all()
        and np.isfinite(initial).all()
        and np.isfinite(developed).all()
    )
    if not complete:
        return {
            "complete_temperature_grid": False,
            "initial_min_MPa_sqrt_m": float("nan"),
            "initial_span_MPa_sqrt_m": float("nan"),
            "initial_low_mid_span_MPa_sqrt_m": float("nan"),
            "developed_span_MPa_sqrt_m": float("nan"),
            "developed_low_mid_span_MPa_sqrt_m": float("nan"),
            "low_temperature_initial_rise_MPa_sqrt_m": float("nan"),
            "initial_high_temperature_loss_MPa_sqrt_m": float("nan"),
            "developed_high_temperature_loss_MPa_sqrt_m": float("nan"),
            "developed_high_temperature_rebound_MPa_sqrt_m": float("nan"),
            "median_R_rise_MPa_sqrt_m": float("nan"),
            "median_abs_R_rise_MPa_sqrt_m": float("nan"),
            "max_abs_R_rise_MPa_sqrt_m": float("nan"),
            "positive_R_rise_temperature_fraction": float("nan"),
            "initial_peak_prominence_MPa_sqrt_m": float("nan"),
            "developed_peak_prominence_MPa_sqrt_m": float("nan"),
        }

    order = np.argsort(temperatures, kind="stable")
    temperatures = temperatures[order]
    initial = initial[order]
    developed = developed[order]
    low_mid = temperatures <= 1100.0
    high = temperatures >= 1200.0
    low = temperatures <= 700.0
    if not low_mid.any() or not high.any() or not low.any():
        return response_metrics([], [], [])

    initial_low_mid = initial[low_mid]
    developed_low_mid = developed[low_mid]
    initial_high = initial[high]
    developed_high = developed[high]
    R_rise = developed - initial
    initial_lowest = float(initial[0])
    low_temperature_rise = float(np.max(initial[low]) - initial_lowest)
    developed_high_min = float(np.min(developed_high))
    developed_high_rebound = float(np.max(developed_high) - developed_high_min)

    return {
        "complete_temperature_grid": True,
        "initial_min_MPa_sqrt_m": float(np.min(initial)),
        "initial_span_MPa_sqrt_m": float(np.ptp(initial)),
        "initial_low_mid_span_MPa_sqrt_m": float(np.ptp(initial_low_mid)),
        "developed_span_MPa_sqrt_m": float(np.ptp(developed)),
        "developed_low_mid_span_MPa_sqrt_m": float(np.ptp(developed_low_mid)),
        "low_temperature_initial_rise_MPa_sqrt_m": low_temperature_rise,
        "initial_high_temperature_loss_MPa_sqrt_m": (
            _median(initial_low_mid) - _median(initial_high)
        ),
        "developed_high_temperature_loss_MPa_sqrt_m": (
            _median(developed_low_mid) - _median(developed_high)
        ),
        "developed_high_temperature_rebound_MPa_sqrt_m": developed_high_rebound,
        "median_R_rise_MPa_sqrt_m": _median(R_rise),
        "median_abs_R_rise_MPa_sqrt_m": _median(np.abs(R_rise)),
        "max_abs_R_rise_MPa_sqrt_m": float(np.max(np.abs(R_rise))),
        "positive_R_rise_temperature_fraction": float(np.mean(R_rise > 0.25)),
        "initial_peak_prominence_MPa_sqrt_m": local_peak_prominence(initial),
        "developed_peak_prominence_MPa_sqrt_m": local_peak_prominence(developed),
    }


def weakT_score(metrics: Mapping[str, Any]) -> tuple[bool, float]:
    finite_names = (
        "initial_min_MPa_sqrt_m",
        "initial_span_MPa_sqrt_m",
        "developed_span_MPa_sqrt_m",
        "low_temperature_initial_rise_MPa_sqrt_m",
        "initial_high_temperature_loss_MPa_sqrt_m",
        "developed_high_temperature_loss_MPa_sqrt_m",
        "median_R_rise_MPa_sqrt_m",
        "max_abs_R_rise_MPa_sqrt_m",
        "positive_R_rise_temperature_fraction",
        "developed_peak_prominence_MPa_sqrt_m",
    )
    finite = bool(metrics.get("complete_temperature_grid", False)) and all(
        math.isfinite(float(metrics.get(name, float("nan")))) for name in finite_names
    )
    if not finite:
        return False, 1.0e6
    initial_min = float(metrics["initial_min_MPa_sqrt_m"])
    initial_span = float(metrics["initial_span_MPa_sqrt_m"])
    developed_span = float(metrics["developed_span_MPa_sqrt_m"])
    low_rise = float(metrics["low_temperature_initial_rise_MPa_sqrt_m"])
    initial_loss = float(metrics["initial_high_temperature_loss_MPa_sqrt_m"])
    developed_loss = float(metrics["developed_high_temperature_loss_MPa_sqrt_m"])
    median_rise = float(metrics["median_R_rise_MPa_sqrt_m"])
    max_abs_rise = float(metrics["max_abs_R_rise_MPa_sqrt_m"])
    positive_fraction = float(metrics["positive_R_rise_temperature_fraction"])
    peak = float(metrics["developed_peak_prominence_MPa_sqrt_m"])
    gate = bool(
        initial_min >= 8.0
        and initial_span <= 6.0
        and developed_span <= 8.0
        and low_rise <= 3.0
        and abs(initial_loss) <= 4.0
        and abs(developed_loss) <= 5.0
        and 1.0 <= median_rise <= 12.0
        and max_abs_rise <= 16.0
        and positive_fraction >= 0.60
        and peak < 4.0
    )
    score = (
        max(8.0 - initial_min, 0.0) / 2.0
        + initial_span / 6.0
        + developed_span / 8.0
        + low_rise / 3.0
        + abs(initial_loss) / 4.0
        + abs(developed_loss) / 5.0
        + abs(median_rise - 5.0) / 5.0
        + max(max_abs_rise - 12.0, 0.0) / 4.0
        + max(0.75 - positive_fraction, 0.0) * 4.0
        + peak / 4.0
    )
    return gate, float(score)


def ceramic_score(metrics: Mapping[str, Any]) -> tuple[bool, float]:
    finite_names = (
        "initial_min_MPa_sqrt_m",
        "initial_low_mid_span_MPa_sqrt_m",
        "developed_low_mid_span_MPa_sqrt_m",
        "low_temperature_initial_rise_MPa_sqrt_m",
        "initial_high_temperature_loss_MPa_sqrt_m",
        "developed_high_temperature_rebound_MPa_sqrt_m",
        "median_abs_R_rise_MPa_sqrt_m",
        "max_abs_R_rise_MPa_sqrt_m",
        "initial_peak_prominence_MPa_sqrt_m",
        "developed_peak_prominence_MPa_sqrt_m",
    )
    finite = bool(metrics.get("complete_temperature_grid", False)) and all(
        math.isfinite(float(metrics.get(name, float("nan")))) for name in finite_names
    )
    if not finite:
        return False, 1.0e6
    initial_min = float(metrics["initial_min_MPa_sqrt_m"])
    initial_span = float(metrics["initial_low_mid_span_MPa_sqrt_m"])
    developed_span = float(metrics["developed_low_mid_span_MPa_sqrt_m"])
    low_rise = float(metrics["low_temperature_initial_rise_MPa_sqrt_m"])
    high_loss = float(metrics["initial_high_temperature_loss_MPa_sqrt_m"])
    rebound = float(metrics["developed_high_temperature_rebound_MPa_sqrt_m"])
    median_abs_rise = float(metrics["median_abs_R_rise_MPa_sqrt_m"])
    max_abs_rise = float(metrics["max_abs_R_rise_MPa_sqrt_m"])
    initial_peak = float(metrics["initial_peak_prominence_MPa_sqrt_m"])
    developed_peak = float(metrics["developed_peak_prominence_MPa_sqrt_m"])
    gate = bool(
        initial_min >= 8.0
        and initial_span <= 5.0
        and developed_span <= 6.0
        and low_rise <= 2.5
        and 0.5 <= high_loss <= 6.0
        and rebound <= 2.0
        and median_abs_rise <= 3.0
        and max_abs_rise <= 6.0
        and initial_peak < 3.0
        and developed_peak < 3.0
    )
    score = (
        max(8.0 - initial_min, 0.0) / 2.0
        + initial_span / 5.0
        + developed_span / 6.0
        + low_rise / 2.5
        + abs(high_loss - 2.5) / 3.0
        + rebound / 2.0
        + median_abs_rise / 3.0
        + max_abs_rise / 6.0
        + initial_peak / 3.0
        + developed_peak / 3.0
    )
    return gate, float(score)


def add_class_scores(
    frame: pd.DataFrame,
    temperatures_K: Sequence[float],
    *,
    initial_prefix: str,
    developed_prefix: str,
    metric_prefix: str,
) -> pd.DataFrame:
    out = frame.copy()
    records: list[dict[str, Any]] = []
    tags = [f"{float(T):g}".replace(".", "p") for T in temperatures_K]
    for row in out.to_dict(orient="records"):
        initial = [float(row[f"{initial_prefix}{tag}"]) for tag in tags]
        developed = [float(row[f"{developed_prefix}{tag}"]) for tag in tags]
        metrics = response_metrics(temperatures_K, initial, developed)
        weak_gate, weak_score_value = weakT_score(metrics)
        ceramic_gate, ceramic_score_value = ceramic_score(metrics)
        records.append(
            {
                **metrics,
                "weakT_gate": weak_gate,
                "weakT_score": weak_score_value,
                "ceramic_gate": ceramic_gate,
                "ceramic_score": ceramic_score_value,
            }
        )
    metrics_frame = pd.DataFrame(records, index=out.index)
    metrics_frame = metrics_frame.rename(
        columns={name: f"{metric_prefix}{name}" for name in metrics_frame.columns}
    )
    return pd.concat([out, metrics_frame], axis=1)


def normalized_active_features(frame: pd.DataFrame, policy: Mapping[str, Any]) -> np.ndarray:
    columns: list[np.ndarray] = []
    dimensions = policy["search_dimensions"]
    for name in ACTIVE_CANDIDATE_PARAMETER_FIELDS:
        if name not in dimensions:
            continue
        spec = dimensions[name]
        values = pd.to_numeric(frame[name], errors="coerce").to_numpy(float)
        low = float(spec["low"])
        high = float(spec["high"])
        if str(spec["mode"]) == "log10_delta":
            values = np.log10(np.maximum(values, 1.0e-300))
            low = math.log10(low)
            high = math.log10(high)
        columns.append(np.clip((values - low) / max(high - low, 1.0e-30), 0.0, 1.0))
    if not columns:
        raise RuntimeError("search policy contains no active feature dimensions")
    return np.column_stack(columns)


def diverse_select(
    ranked: pd.DataFrame,
    policy: Mapping[str, Any],
    *,
    count: int,
    score_column: str,
    gate_column: str,
) -> pd.DataFrame:
    if count < 1:
        raise ValueError("selection count must be positive")
    ordered = ranked.sort_values(
        [gate_column, score_column, "candidate_id"],
        ascending=[False, True, True],
        kind="stable",
    ).reset_index(drop=True)
    if len(ordered) <= count:
        result = ordered.copy()
        result["diversity_rank"] = np.arange(1, len(result) + 1)
        return result
    pool_count = min(len(ordered), max(count * 16, count))
    pool = ordered.head(pool_count).copy().reset_index(drop=True)
    features = normalized_active_features(pool, policy)
    scores = pd.to_numeric(pool[score_column], errors="coerce").to_numpy(float)
    finite_scores = scores[np.isfinite(scores)]
    scale = max(float(np.percentile(finite_scores, 80)) if finite_scores.size else 1.0, 1.0e-12)
    gates = pool[gate_column].fillna(False).astype(bool).to_numpy()
    selected = [0]
    minimum_distance = np.sum(np.square(features - features[0]), axis=1)
    minimum_distance[0] = -np.inf
    for _ in range(1, count):
        quality = np.exp(-np.clip(scores / scale, 0.0, 50.0))
        gate_weight = np.where(gates, 2.0, 1.0)
        acquisition = minimum_distance * (0.20 + 0.80 * quality) * gate_weight
        acquisition[selected] = -np.inf
        index = int(np.argmax(acquisition))
        selected.append(index)
        distance = np.sum(np.square(features - features[index]), axis=1)
        minimum_distance = np.minimum(minimum_distance, distance)
    result = pool.iloc[selected].copy()
    result["diversity_rank"] = np.arange(1, len(result) + 1)
    return result


__all__ = [
    "add_class_scores",
    "ceramic_score",
    "diverse_select",
    "local_peak_prominence",
    "normalized_active_features",
    "response_metrics",
    "weakT_score",
]
