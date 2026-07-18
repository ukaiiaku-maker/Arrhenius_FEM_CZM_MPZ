from __future__ import annotations

import numpy as np
import pandas as pd

from arrhenius_fracture.reduced_campaign_front_v9104 import ReducedFrontSettings
from evaluate_mechanism_ablation_mpz_v9_10_4_9 import (
    OPTIONAL_MODE,
    _solver_configuration,
    candidate_sensitivity,
    curve_metrics,
)
from select_mechanism_ablation_candidates_mpz_v9_10_4_9 import select_candidates


def test_selects_strict_pass_then_best_near_miss_in_each_bracket():
    rows = [
        {
            "candidate_id": "a0",
            "transition_bracket": "T0400_0500K",
            "coarse_transition_low_T_K": 400.0,
            "moving_1d_valid": True,
            "moving_1d_accept": False,
            "moving_1d_objective": 0.1,
        },
        {
            "candidate_id": "a1",
            "transition_bracket": "T0400_0500K",
            "coarse_transition_low_T_K": 400.0,
            "moving_1d_valid": True,
            "moving_1d_accept": True,
            "moving_1d_objective": 10.0,
        },
        {
            "candidate_id": "b0",
            "transition_bracket": "T0500_0600K",
            "coarse_transition_low_T_K": 500.0,
            "moving_1d_valid": True,
            "moving_1d_accept": False,
            "moving_1d_objective": 0.2,
        },
    ]
    selected = select_candidates(
        pd.DataFrame(rows), per_bracket=1, expected_brackets=2
    )
    assert selected.candidate_id.tolist() == ["a1", "b0"]
    assert selected.ablation_selection_basis.tolist() == [
        "strict_1d_pass",
        "best_valid_near_miss",
    ]


def test_curve_metrics_uses_complete_four_point_bracket():
    T = np.asarray([700.0, 733.3333333333, 766.6666666667, 800.0])
    K = np.asarray([8.429344, 8.311141, 9.396201, 19.270259])
    result = curve_metrics(T, K)
    assert result["curve_valid"] is True
    assert result["endpoint_ratio"] > 2.0
    assert result["transition_width_K"] < 100.0


def test_background_field_off_preserves_full_mode_but_zeroes_rho0():
    settings = ReducedFrontSettings(rho0_m2=5.0e12)
    solver_mode, changed, description = _solver_configuration(OPTIONAL_MODE, settings)
    assert solver_mode == "full"
    assert changed.rho0_m2 == 0.0
    assert settings.rho0_m2 == 5.0e12
    assert description == "full_with_rho0_zero"


def test_emission_opening_priority_requires_blunting_not_shielding_or_background():
    row = pd.Series(
        {
            "candidate_id": "candidate",
            "transition_bracket": "T0700_0800K",
            "coarse_transition_low_T_K": 700.0,
        }
    )
    mode_rows = {
        "full": {
            "rise_MPa_sqrt_m": 12.0,
            "endpoint_ratio": 2.2,
            "low_endpoint_K": 10.0,
            "high_endpoint_K": 22.0,
            "monotonic_fraction": 1.0,
            "transition_width_K": 40.0,
            "max_K_shield_MPa_sqrt_m": 0.1,
            "low_cumulative_emitted": 2.0,
            "high_cumulative_emitted": 100.0,
        },
        "plasticity_off": {
            "rise_MPa_sqrt_m": -0.5,
            "endpoint_ratio": 0.95,
        },
        "blunting_off": {"rise_MPa_sqrt_m": 2.0},
        "backstress_off": {"rise_MPa_sqrt_m": 14.0},
        "shielding_off": {"rise_MPa_sqrt_m": 11.5},
        OPTIONAL_MODE: {"rise_MPa_sqrt_m": 10.0},
    }
    result = candidate_sensitivity(row, mode_rows)
    assert result["two_d_emission_opening_priority"] is True
    assert result["ablation_blunting_sensitivity_fraction"] > 0.5
    assert abs(result["ablation_shielding_sensitivity_fraction"]) < 0.2
    assert result["ablation_background_off_retained_rise_fraction"] > 0.75
