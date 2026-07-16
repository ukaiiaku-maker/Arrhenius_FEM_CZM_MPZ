from pathlib import Path


def test_mode_i_runner_routes_to_v91854():
    text = Path("run_mpz_v9_18_5_4_mode_i_rcurve.py").read_text()
    assert "arrhenius_fracture.mode_i_first_passage_v9_18_5_4" in text


def test_campaign_routes_to_v91854_runner():
    text = Path("run_mpz_v9_18_5_4_persistent_plastic_wake.py").read_text()
    assert "run_mpz_v9_18_5_4_mode_i_rcurve.py" in text


def test_v91854_wraps_quality_selected_corridor():
    text = Path("arrhenius_fracture/mode_i_first_passage_v9_18_5_4.py").read_text()
    assert "_v91853.main(user_args)" in text
    assert "_v9185._strict_quality_advance = _strict_quality_advance_v91854" in text
