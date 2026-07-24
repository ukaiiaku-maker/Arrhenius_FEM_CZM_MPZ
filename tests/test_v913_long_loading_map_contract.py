from __future__ import annotations

import pytest

from scripts.extract_v10222_long_rcurve_loading_map import (
    validate_calibrated_stochastic_contract,
    validate_expected_prefix,
)


def calibrated_stack() -> dict:
    return {
        "event_length_uses_same_integrated_hazard_threshold": True,
        "stochastic_hazard": {
            "mode": "exponential",
            "distribution": "exponential_unit_mean",
        },
        "stochastic_avalanche": {
            "mode": "threshold_scaled",
            "minimum_factor": 0.5,
            "maximum_factor": 4.0,
            "geometry_subsegment_fraction": 0.1,
        },
    }


def loading_map() -> dict:
    return {
        "K_per_U_MPa_sqrt_m_per_m": [1.0, 2.0],
        "threshold_actions": [0.5, 1.25],
        "path_advances_m": [2.5e-6, 6.25e-6],
        "projected_advances_m": [2.5e-6, 6.0e-6],
        "seed": 3621,
        "nominal_dU_m": 2.0e-7,
        "nominal_dt_s": 8.4,
    }


def test_calibrated_stochastic_contract_rejects_deterministic_reference() -> None:
    accepted = validate_calibrated_stochastic_contract(calibrated_stack())
    assert accepted["hazard_mode"] == "exponential"
    assert accepted["event_length_mode"] == "threshold_scaled"

    invalid = calibrated_stack()
    invalid["stochastic_hazard"] = {
        "mode": "deterministic",
        "distribution": "delta_at_one",
    }
    invalid["stochastic_avalanche"] = {
        "mode": "fixed",
        "minimum_factor": 0.5,
        "maximum_factor": 4.0,
        "geometry_subsegment_fraction": 0.1,
    }
    with pytest.raises(RuntimeError, match="stochastic_hazard.mode"):
        validate_calibrated_stochastic_contract(invalid)


def test_expected_prefix_requires_exact_calibrated_history() -> None:
    expected = loading_map()
    current = {
        **loading_map(),
        "K_per_U_MPa_sqrt_m_per_m": [1.0, 2.0, 3.0],
        "threshold_actions": [0.5, 1.25, 0.75],
        "path_advances_m": [2.5e-6, 6.25e-6, 3.75e-6],
        "projected_advances_m": [2.5e-6, 6.0e-6, 3.75e-6],
    }
    audit = validate_expected_prefix(current, expected)
    assert audit["prefix_events"] == 2
    assert audit["maximum_absolute_errors"]["threshold_actions"] == 0.0

    current["threshold_actions"][1] = 1.0
    with pytest.raises(RuntimeError, match="does not reproduce calibrated prefix"):
        validate_expected_prefix(current, expected)
