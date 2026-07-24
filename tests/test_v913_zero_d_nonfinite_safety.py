from __future__ import annotations

import json
import warnings

import numpy as np

import arrhenius_fracture.zero_d_search_v913 as search_helpers
from scripts.run_v913_zero_d_large_search_safe import (
    curve_metrics_matrix_safe,
    install_safety_patch,
    json_safe,
    local_peak_metrics_safe,
)


def test_json_safe_converts_nested_nonfinite_values() -> None:
    payload = {
        "positive": np.float64(2.5),
        "nan": float("nan"),
        "positive_inf": np.float64(float("inf")),
        "negative_inf": float("-inf"),
        "nested": [np.float32(1.0), float("nan")],
    }
    safe = json_safe(payload)
    assert safe == {
        "positive": 2.5,
        "nan": None,
        "positive_inf": None,
        "negative_inf": None,
        "nested": [1.0, None],
    }
    json.dumps(safe, allow_nan=False)


def test_curve_metrics_all_nan_row_is_invalid_without_warning() -> None:
    temperatures = np.asarray([700.0, 900.0, 1000.0, 1100.0, 1300.0])
    curves = np.asarray(
        [
            [np.nan, np.nan, np.nan, np.nan, np.nan],
            [20.0, 30.0, 45.0, 35.0, 38.0],
        ]
    )
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        result = curve_metrics_matrix_safe(
            temperatures,
            curves,
            peak_min=850.0,
            peak_max=1100.0,
        )
    assert np.isnan(result["peak_temperature_K"][0])
    assert np.isnan(result["two_sided_prominence_MPa_sqrt_m"][0])
    assert result["peak_internal"][0] == 0
    assert result["peak_temperature_K"][1] == 1000.0
    assert result["two_sided_prominence_MPa_sqrt_m"][1] == 10.0
    assert result["post_peak_drop_MPa_sqrt_m"][1] == 10.0
    assert result["high_temperature_rebound_MPa_sqrt_m"][1] == -7.0


def test_install_patch_replaces_proxy_module_metric() -> None:
    install_safety_patch()
    assert search_helpers._curve_metrics_matrix is curve_metrics_matrix_safe

    temperatures = np.asarray([700.0, 900.0, 1000.0, 1100.0, 1300.0])
    curves = np.full((1, temperatures.size), np.nan)
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        result = search_helpers._curve_metrics_matrix(
            temperatures,
            curves,
            peak_min=850.0,
            peak_max=1100.0,
        )
    assert np.isnan(result["peak_temperature_K"][0])
    assert result["peak_internal"][0] == 0


def test_scalar_boundary_peak_uses_nan_not_infinity() -> None:
    result = local_peak_metrics_safe(
        [700.0, 900.0, 1100.0, 1300.0],
        [20.0, 30.0, 40.0, 50.0],
    )
    assert result["peak_temperature_K"] == 1300.0
    assert result["peak_internal"] is False
    assert np.isnan(result["two_sided_prominence"])
    assert np.isnan(result["post_peak_drop"])
    assert np.isnan(result["high_temperature_rebound"])
