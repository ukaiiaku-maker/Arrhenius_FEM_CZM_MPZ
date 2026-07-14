import search_mpz_peierls_taylor_parameters_v94 as v94


def test_wrapper_preserves_original_evaluator_before_monkeypatch():
    assert v94._LEGACY_EVALUATE_ONE is not v94.evaluate_one


def test_v94_entropy_prior_is_bounded():
    rows = v94.sample_transport_parameters(32, 94017)
    assert rows["pt_entropy_multiplier"].min() >= 0.25
    assert rows["pt_entropy_multiplier"].max() <= 8.0
