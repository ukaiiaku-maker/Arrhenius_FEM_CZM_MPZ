import json

import numpy as np
import pytest

from arrhenius_fracture.dbtt_temperature_schedule_v91043 import (
    schedule_from_bracket,
    schedule_from_candidate_row,
)


def test_four_points_span_candidate_specific_100K_bracket_with_shelf_anchors():
    coarse = np.arange(300.0, 1101.0, 100.0)
    schedule = schedule_from_bracket(
        coarse,
        600.0,
        700.0,
        refinement_points=4,
        shelf_anchor_count=2,
    )
    assert np.allclose(schedule.transition_temperatures_K, [600.0, 633.3333333333, 666.6666666667, 700.0])
    assert schedule.low_anchor_temperatures_K == (300.0, 500.0)
    assert schedule.high_anchor_temperatures_K == (800.0, 1100.0)
    assert np.allclose(
        schedule.evaluation_temperatures_K,
        [300.0, 500.0, 600.0, 633.3333333333, 666.6666666667, 700.0, 800.0, 1100.0],
    )


def test_schedule_round_trips_through_candidate_columns():
    coarse = np.arange(300.0, 1101.0, 100.0)
    original = schedule_from_bracket(coarse, 800.0, 900.0)
    row = original.to_columns()
    recovered = schedule_from_candidate_row(row, coarse)
    assert recovered == original
    assert json.loads(row["refinement_transition_temperatures_K"])[0] == 800.0


def test_edge_bracket_without_a_shelf_anchor_is_rejected():
    coarse = np.arange(300.0, 1101.0, 100.0)
    with pytest.raises(ValueError, match="at least one coarse temperature on each shelf"):
        schedule_from_bracket(coarse, 300.0, 400.0)
