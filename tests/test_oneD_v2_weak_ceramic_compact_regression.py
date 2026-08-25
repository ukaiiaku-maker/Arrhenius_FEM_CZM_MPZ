from scripts.run_oneD_v2_weak_ceramic_compact_regression import ROWS, TEMPERATURES


def test_compact_regression_uses_frozen_selected_rows():
    assert ROWS == {
        "weak-T": "oneD_v2_focused_weak_T_0016",
        "ceramic-like": "oneD_v2_focused_ceramic_like_0018",
    }


def test_compact_regression_uses_required_temperatures():
    assert TEMPERATURES == (300.0, 1000.0, 1200.0)
