from pathlib import Path


def test_v91856_mode_i_routing():
    text = Path("run_mpz_v9_18_5_6_mode_i_rcurve.py").read_text()
    assert "arrhenius_fracture.mode_i_first_passage_v9_18_5_6" in text


def test_v91856_persistent_wake_routing():
    text = Path("run_mpz_v9_18_5_6_persistent_plastic_wake.py").read_text()
    assert "run_mpz_v9_18_5_6_mode_i_rcurve.py" in text


def test_v91856_sweep_routing():
    text = Path("run_mpz_v9_18_5_6_persistent_plastic_wake_sweep.sh").read_text()
    assert "run_mpz_v9_18_5_6_persistent_plastic_wake.py" in text
    assert "explicit_quality_wrapper_chain_v91856.json" in text
