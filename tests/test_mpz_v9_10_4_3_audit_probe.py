import numpy as np

import audit_mpz_v9_10_4_2_current_dbtt as guard


def test_three_temperature_probe_is_invalid_but_does_not_raise():
    guard._RECORDED_TRANSITIONS.clear()
    result = guard.guarded_best_adjacent_transition(
        np.array([300.0, 700.0, 1100.0]),
        np.array([15.0, 28.0, 39.0]),
    )
    assert result["valid"] is False
    assert result["reason"] == "insufficient_temperatures_for_two_shelves"
    assert np.isnan(result["shelf_ratio"])
