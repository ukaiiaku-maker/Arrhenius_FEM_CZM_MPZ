from scripts.prepare_oneD_v2_peak_dbtt_R_pf_transfer import ROWS


def test_pf_transfer_input_has_two_controls_and_two_finalists():
    assert len(ROWS) == 4
    assert {row[1] for row in ROWS} == {
        "Peak_control", "Peak_R", "DBTT_control", "DBTT_R",
    }
    assert len({row[2] for row in ROWS}) == 4
