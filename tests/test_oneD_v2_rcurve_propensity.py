from __future__ import annotations

import json

import pytest

from reduced_fracture_v2.rcurve_propensity import (
    reload_separated_onsets,
    summarize_rcurve_propensity,
)


def _event(index, avalanche, opening_um, kj, *, reload=False, length_um=5.0):
    return {
        "event_index": index,
        "physical_avalanche_index": avalanche,
        "event_time_s": float(index),
        "event_opening_m": opening_um * 1.0e-6,
        "reload_opening_m": None if index == 0 else (opening_um - 10.0) * 1.0e-6,
        "reload_separated": reload,
        "extension_before_m": index * 5.0e-6,
        "event_length_m": length_um * 1.0e-6,
        "native_KJ_MPa_sqrt_m": kj,
        "qualified_G_J_per_m2": kj * kj,
        "tip_radius_m": 1.0e-6,
        "front_width_m": 2.0e-6,
        "mobile_density_m2": 3.0,
        "retained_density_m2": 4.0,
        "source_multiplicity": 5.0,
        "backstress_Pa": 6.0,
    }


def test_only_reload_separated_pre_event_states_are_resistance_candidates():
    result = {
        "events": [
            _event(0, 0, 10.0, 20.0),
            _event(1, 0, 10.1, 45.0),
            _event(2, 1, 25.0, 30.0, reload=True),
            _event(3, 1, 25.1, 80.0),
        ],
        "physical_avalanche_count": 2,
        "largest_avalanche_fraction": 0.5,
    }
    onsets = reload_separated_onsets(result)
    assert [item["event_index"] for item in onsets] == [0, 2]
    summary = summarize_rcurve_propensity(result)
    assert summary["deltaK_reinit_MPa_sqrt_m"] == pytest.approx(10.0)
    assert summary["N_reinit"] == 1
    assert summary["onset_sequence_monotonicity"] == "NONDECREASING"
    assert json.loads(summary["in_avalanche_native_drive_json"]) == [20.0, 45.0, 30.0, 80.0]
    assert summary["in_avalanche_drive_interpretation"].endswith("NOT_RESISTANCE")


def test_single_avalanche_has_zero_reinitiation_increment():
    result = {
        "events": [_event(0, 0, 10.0, 20.0), _event(1, 0, 10.1, 50.0)],
        "physical_avalanche_count": 1,
        "largest_avalanche_fraction": 1.0,
    }
    summary = summarize_rcurve_propensity(result)
    assert summary["N_reinit"] == 0
    assert summary["deltaK_reinit_MPa_sqrt_m"] == 0.0
    assert summary["onset_sequence_monotonicity"] == "NOT_APPLICABLE"


def test_target_completion_is_never_relabelled_physical_arrest():
    result = {
        "events": [_event(0, 0, 10.0, 20.0)],
        "physical_avalanche_count": 1,
        "largest_avalanche_fraction": 1.0,
        "status": "TARGET_RIGHT_CENSORED",
    }
    summary = summarize_rcurve_propensity(result)
    assert summary["target_termination_semantics"] == "RIGHT_CENSORED_NOT_PHYSICAL_ARREST"
