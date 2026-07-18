#!/usr/bin/env python3
"""Corrected four-temperature moving-interface DBTT evaluation for v9.10.4.8.

The four detailed temperatures all lie inside one analytically predicted 100 K
transition bracket.  They are transition samples, not two low-shelf and two
high-shelf samples.  Therefore acceptance is based on the lower and upper
bracket endpoints, transition monotonicity and width, and matched
plasticity-off/mechanistic diagnostics.  The legacy midpoint-separation ratio
is retained as a diagnostic but is not an acceptance gate.
"""
from __future__ import annotations

from typing import Any

import numpy as np

import evaluate_dynamic_1d_mpz_v9_10_4_7 as _base


LOW_K_MIN = 8.0
LOW_K_MAX = 25.0
HIGH_K_MAX = 70.0
MIN_ENDPOINT_RATIO = 2.0
MIN_MONOTONIC_FRACTION = 0.90
MAX_PLASTICITY_OFF_RATIO = 1.25
MIN_MECHANISTIC_FRACTION = 0.60
MAX_TRANSITION_WIDTH_K = 100.0


def corrected_transition_metrics(
    temperatures: np.ndarray,
    full_K: np.ndarray,
    off_K: np.ndarray,
) -> dict[str, Any]:
    """Score a four-point transition using the complete 100 K bracket."""
    T = np.asarray(temperatures, dtype=float)
    K = np.asarray(full_K, dtype=float)
    Koff = np.asarray(off_K, dtype=float)
    if T.size != 4 or K.size != 4 or Koff.size != 4:
        return {
            "moving_1d_valid": False,
            "moving_1d_objective": 1.0e12,
            "moving_1d_accept": False,
            "moving_1d_reason": "four_transition_temperatures_required",
            "moving_1d_gate_version": "V9_10_4_8_ENDPOINT_100K",
        }
    if not np.all(np.isfinite(T)) or not np.all(np.isfinite(K)) or not np.all(np.isfinite(Koff)):
        return {
            "moving_1d_valid": False,
            "moving_1d_objective": 1.0e12,
            "moving_1d_accept": False,
            "moving_1d_reason": "incomplete_first_passage",
            "moving_1d_gate_version": "V9_10_4_8_ENDPOINT_100K",
        }
    if not np.all(np.diff(T) > 0.0):
        return {
            "moving_1d_valid": False,
            "moving_1d_objective": 1.0e12,
            "moving_1d_accept": False,
            "moving_1d_reason": "nonmonotone_temperature_schedule",
            "moving_1d_gate_version": "V9_10_4_8_ENDPOINT_100K",
        }

    low_K = float(K[0])
    high_K = float(K[-1])
    rise = high_K - low_K
    ratio = high_K / max(low_K, 1.0e-12)
    midpoint_separation = float(np.min(K[2:]) / max(np.max(K[:2]), 1.0e-12))

    increments = np.diff(K)
    total_variation = float(np.sum(np.abs(increments)))
    positive_variation = float(np.sum(np.maximum(increments, 0.0)))
    monotonic_fraction = (
        1.0 if total_variation <= 1.0e-12 else positive_variation / total_variation
    )
    max_negative_fraction = float(
        np.max(np.maximum(-increments, 0.0)) / max(abs(rise), 1.0e-12)
    )

    off_ratio = float(Koff[-1] / max(Koff[0], 1.0e-12))
    plastic = K - Koff
    mechanistic_fraction = float(
        (plastic[-1] - plastic[0]) / max(rise, 1.0e-12)
    )

    T10 = float("nan")
    T90 = float("nan")
    width = float("inf")
    if rise > 0.0:
        normalized = (K - low_K) / rise

        def crossing(level: float) -> float:
            if normalized[0] >= level:
                return float(T[0])
            for i in range(len(T) - 1):
                y0 = float(normalized[i])
                y1 = float(normalized[i + 1])
                if y0 < level <= y1 and y1 > y0:
                    fraction = (level - y0) / (y1 - y0)
                    return float(T[i] + fraction * (T[i + 1] - T[i]))
            return float("nan")

        T10 = crossing(0.10)
        T90 = crossing(0.90)
        if np.isfinite(T10) and np.isfinite(T90) and T90 >= T10:
            width = float(T90 - T10)

    penalties = {
        "low_floor": max(LOW_K_MIN - low_K, 0.0) / 2.0,
        "low_ceiling": max(low_K - LOW_K_MAX, 0.0) / 3.0,
        "high_ceiling": max(high_K - HIGH_K_MAX, 0.0) / 5.0,
        "endpoint_ratio": max(MIN_ENDPOINT_RATIO - ratio, 0.0) / 0.25,
        "monotonicity": max(MIN_MONOTONIC_FRACTION - monotonic_fraction, 0.0) / 0.10,
        "plasticity_off_ratio": max(off_ratio - MAX_PLASTICITY_OFF_RATIO, 0.0) / 0.10,
        "mechanistic_fraction": max(MIN_MECHANISTIC_FRACTION - mechanistic_fraction, 0.0) / 0.15,
        "transition_width": (
            max(width - MAX_TRANSITION_WIDTH_K, 0.0) / 25.0
            if np.isfinite(width)
            else 20.0
        ),
    }
    objective = float(sum(value * value for value in penalties.values()))

    checks = [
        (low_K >= LOW_K_MIN, "moving_1d_low_toughness_below_floor"),
        (low_K <= LOW_K_MAX, "moving_1d_low_toughness_above_ceiling"),
        (high_K <= HIGH_K_MAX, "moving_1d_high_toughness_above_ceiling"),
        (rise > 0.0, "moving_1d_nonpositive_transition_rise"),
        (ratio >= MIN_ENDPOINT_RATIO, "moving_1d_endpoint_ratio_below_two"),
        (monotonic_fraction >= MIN_MONOTONIC_FRACTION, "moving_1d_transition_nonmonotone"),
        (off_ratio <= MAX_PLASTICITY_OFF_RATIO, "moving_1d_cleavage_only_T_dependence"),
        (mechanistic_fraction >= MIN_MECHANISTIC_FRACTION, "moving_1d_plastic_increment_too_small"),
        (width <= MAX_TRANSITION_WIDTH_K, "moving_1d_transition_wider_than_bracket"),
    ]
    accepted = True
    reason = "moving_1d_endpoint_transition_passed"
    for passed, failure in checks:
        if not passed:
            accepted = False
            reason = failure
            break

    return {
        "moving_1d_valid": True,
        "moving_1d_objective": objective,
        "moving_1d_accept": accepted,
        "moving_1d_reason": reason,
        "moving_1d_gate_version": "V9_10_4_8_ENDPOINT_100K",
        "moving_1d_low_endpoint_K": low_K,
        "moving_1d_high_endpoint_K": high_K,
        "moving_1d_edge_ratio": ratio,
        "moving_1d_robust_ratio": midpoint_separation,
        "moving_1d_midpoint_separation_ratio": midpoint_separation,
        "moving_1d_rise_MPa_sqrt_m": rise,
        "moving_1d_monotonic_fraction": monotonic_fraction,
        "moving_1d_max_negative_fraction": max_negative_fraction,
        "moving_1d_plasticity_off_ratio": off_ratio,
        "moving_1d_mechanistic_fraction": mechanistic_fraction,
        "moving_1d_T10_K": T10,
        "moving_1d_T90_K": T90,
        "moving_1d_transition_width_K": width,
        **{f"moving_1d_penalty_{key}": value for key, value in penalties.items()},
    }


_base._transition_metrics = corrected_transition_metrics


if __name__ == "__main__":
    _base.main()
