import numpy as np

from arrhenius_fracture.dbtt_temperature_schedule_v91043 import (
    fixed_bracket_transition_metrics,
    schedule_from_bracket,
)


def test_refined_metric_scores_complete_100K_bracket_not_one_33K_interval():
    coarse = np.arange(300.0, 1101.0, 100.0)
    schedule = schedule_from_bracket(coarse, 600.0, 700.0)
    T = np.asarray(schedule.evaluation_temperatures_K)
    K = np.array([15.0, 15.0, 15.0, 20.0, 25.0, 30.0, 30.0, 30.0])
    Koff = np.full_like(K, 15.0)
    result = fixed_bracket_transition_metrics(T, K, schedule, plasticity_off_toughness=Koff)
    assert result["valid"]
    assert result["shelf_ratio"] == 2.0
    assert result["jump_concentration"] == 1.0
    assert result["transition_width_K"] <= 100.0
    assert result["transition_monotonic_fraction"] == 1.0
    # No individual 33 K jump contains 75% of the 15-unit shelf change, but the
    # full 100 K bracket does. This is the intended refined-stage interpretation.
    assert np.max(np.diff(K[2:6])) / 15.0 < 0.75


def test_broad_ramp_outside_selected_bracket_has_low_concentration():
    coarse = np.arange(300.0, 1101.0, 100.0)
    schedule = schedule_from_bracket(coarse, 600.0, 700.0)
    T = np.asarray(schedule.evaluation_temperatures_K)
    K = np.array([15.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0, 35.0])
    result = fixed_bracket_transition_metrics(T, K, schedule)
    assert result["valid"]
    assert result["jump_concentration"] < 0.75
    assert result["penalties"]["bracket_concentration"] > 0.0
