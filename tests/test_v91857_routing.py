from pathlib import Path


def test_v91857_mode_i_routes_to_v91856_quality_module():
    text = Path("run_mpz_v9_18_5_7_mode_i_rcurve.py").read_text()
    assert "arrhenius_fracture.mode_i_first_passage_v9_18_5_6" in text


def test_v91857_wake_routes_to_v91857_mode_i_driver():
    text = Path("run_mpz_v9_18_5_7_persistent_plastic_wake.py").read_text()
    assert "run_mpz_v9_18_5_7_mode_i_rcurve.py" in text


def test_v91857_sweep_uses_subsegment_aware_auditor():
    text = Path("run_mpz_v9_18_5_7_persistent_plastic_wake_sweep.sh").read_text()
    assert "audit_v91857_subsegment_quality.py" in text
    assert "len(accepted) != expected" not in text
