from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from arrhenius_fracture.emergent_gnd_contract_v913 import (
    ACTIVE_CANDIDATE_PARAMETER_FIELDS,
)
from scripts.run_v913_zero_d_large_search import (
    FIXED_ACTIVE_FIELDS,
    VARIABLE_FIELDS,
    _curve_metrics_matrix,
    _load_policy,
    _sample_rows,
)


def _policy_path() -> Path:
    return Path("mpz_v9_13_zero_d_large_search_policy.json")


def test_policy_covers_current_active_contract() -> None:
    policy = _load_policy(_policy_path())
    assert set(policy["search_dimensions"]) == set(VARIABLE_FIELDS)
    assert set(VARIABLE_FIELDS) | set(FIXED_ACTIVE_FIELDS) == set(
        ACTIVE_CANDIDATE_PARAMETER_FIELDS
    )
    assert "source_refresh_length_um" not in policy["search_dimensions"]
    assert "recovery_H0_eV" not in policy["search_dimensions"]


def test_sobol_sampling_is_deterministic_and_bounded() -> None:
    policy = _load_policy(_policy_path())
    anchor = {"candidate_id": policy["anchor_candidate_ids"][0]}
    for name, value in FIXED_ACTIVE_FIELDS.items():
        anchor[name] = value
    for name, spec in policy["search_dimensions"].items():
        if spec["mode"] == "log10_delta":
            anchor[name] = np.sqrt(float(spec["low"]) * float(spec["high"]))
        else:
            anchor[name] = 0.5 * (float(spec["low"]) + float(spec["high"]))
    anchors = pd.DataFrame([anchor])
    a = _sample_rows(
        start=0,
        count=16,
        total_samples=16,
        seed=913100,
        anchors=anchors,
        policy=policy,
    )
    b = _sample_rows(
        start=0,
        count=16,
        total_samples=16,
        seed=913100,
        anchors=anchors,
        policy=policy,
    )
    pd.testing.assert_frame_equal(a, b)
    for name, spec in policy["search_dimensions"].items():
        assert a[name].between(float(spec["low"]), float(spec["high"])).all()


def test_curve_metrics_detect_local_peak_and_high_temperature_rebound() -> None:
    temperatures = np.asarray([800, 900, 1000, 1100, 1200, 1300], dtype=float)
    curves = np.asarray(
        [
            [30, 40, 55, 44, 48, 58],
            [30, 35, 40, 45, 50, 55],
        ],
        dtype=float,
    )
    metrics = _curve_metrics_matrix(
        temperatures,
        curves,
        peak_min=850.0,
        peak_max=1100.0,
    )
    assert metrics["peak_temperature_K"][0] == 1000.0
    assert metrics["two_sided_prominence_MPa_sqrt_m"][0] == 11.0
    assert metrics["high_temperature_rebound_MPa_sqrt_m"][0] == 3.0
    assert metrics["peak_internal"][0] == 1
    assert metrics["peak_internal"][1] == 0
