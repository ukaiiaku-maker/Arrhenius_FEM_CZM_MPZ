import pytest

from scripts.check_oneD_v2_peak_dbtt_R_exact_finalists import FINALISTS, _segments


def test_exact_finalists_use_canonical_seeds():
    assert FINALISTS["Peak"][2] == 8666
    assert FINALISTS["DBTT"][2] == 1008666


def test_pf_reconstruction_uses_exact_five_micron_straight_segments():
    events = _segments(12.0e-6)
    lengths = [event["p1_m"][0] - event["p0_m"][0] for event in events]
    assert lengths == pytest.approx([5.0e-6, 5.0e-6, 2.0e-6], rel=0.0, abs=1.0e-18)
    assert all(event["p0_m"][1] == event["p1_m"][1] == 0.0 for event in events)
