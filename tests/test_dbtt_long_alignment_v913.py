from __future__ import annotations

import math

import pytest

from arrhenius_fracture.dbtt_long_alignment_v913 import (
    checkpoint_from_events,
    checkpoint_reached,
    peak_drift_classification,
    peak_metrics,
)


def events():
    return [
        {"cumulative_projected_extension_m": 25.0e-6, "K_MPa_sqrt_m": 20.0},
        {"cumulative_projected_extension_m": 52.0e-6, "K_MPa_sqrt_m": 30.0},
    ]


def test_checkpoint_is_strict_when_loading_map_is_exhausted():
    assert checkpoint_from_events(events(), 50.0) == pytest.approx(30.0)
    assert checkpoint_reached(events(), 50.0)
    assert math.isnan(checkpoint_from_events(events(), 100.0))
    assert not checkpoint_reached(events(), 100.0)
    assert checkpoint_from_events(events(), 100.0, strict=False) == pytest.approx(30.0)


def test_peak_metrics_detect_interior_peak_drop_and_rebound():
    metrics = peak_metrics(
        [700.0, 800.0, 900.0, 1000.0, 1100.0, 1200.0],
        [20.0, 30.0, 45.0, 60.0, 50.0, 55.0],
    )
    assert metrics.peak_temperature_K == pytest.approx(1000.0)
    assert not metrics.peak_at_boundary
    assert metrics.peak_rise == pytest.approx(40.0)
    assert metrics.post_peak_drop == pytest.approx(10.0)
    assert metrics.final_rebound == pytest.approx(5.0)
    assert metrics.peak_prominence == pytest.approx(10.0)


def test_peak_drift_classification():
    assert peak_drift_classification([900.0, 925.0, 950.0]) == (50.0, "stable")
    assert peak_drift_classification([900.0, 950.0, 1000.0]) == (
        100.0,
        "moderate_drift",
    )
    assert peak_drift_classification([800.0, 1000.0])[1] == "extension_dependent"
