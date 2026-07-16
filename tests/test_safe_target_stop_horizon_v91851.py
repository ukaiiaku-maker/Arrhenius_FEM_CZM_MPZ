from __future__ import annotations

import json
import operator
from types import SimpleNamespace

from arrhenius_fracture.mode_i_first_passage_v9_18_5_1 import (
    SafeDynamicStepHorizon,
)


def test_horizon_is_strict_integer_compatible_before_stop():
    controller = SimpleNamespace(v9185_stop_requested=False)
    horizon = SafeDynamicStepHorizon(15, controller)

    assert isinstance(horizon, int)
    assert int(horizon) == 15
    assert operator.index(horizon) == 15
    assert list(range(horizon)) == list(range(15))
    assert horizon * 8.4 == 126.0
    assert 2 * horizon == 30
    assert horizon + 3 == 18
    assert 3 + horizon == 18
    assert horizon - 3 == 12
    assert json.dumps({"steps": horizon}) == '{"steps": 15}'
    assert 0 < horizon
    assert 14 < horizon
    assert not (15 < horizon)


def test_horizon_stops_only_reflected_loop_comparison_after_commit():
    controller = SimpleNamespace(v9185_stop_requested=False)
    horizon = SafeDynamicStepHorizon(15000, controller)

    accepted_step = 1393
    assert accepted_step < horizon

    controller.v9185_stop_requested = True
    assert not (accepted_step < horizon)

    # The nominal horizon remains an ordinary integer everywhere else.
    assert int(horizon) == 15000
    assert operator.index(horizon) == 15000
    assert horizon * 2 == 30000
    assert str(horizon) == "15000"


def test_zero_step_and_boundary_comparisons_remain_well_defined():
    controller = SimpleNamespace(v9185_stop_requested=False)
    horizon = SafeDynamicStepHorizon(0, controller)
    assert not (-1 >= horizon)
    assert not (0 < horizon)
    controller.v9185_stop_requested = True
    assert not (-1 < horizon)
