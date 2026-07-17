from analyze_v10_0_4_parameter_matrix import MatrixCase, certify


def row(material, T, Kc, fingerprint):
    return MatrixCase(
        case_id=f"{material}_{T:g}",
        material=material,
        temperature_K=T,
        Kc_MPa_sqrt_m=Kc,
        event_Uapp_m=1e-4,
        event_Ftop_N=1.0,
        max_N_em=3.0,
        source_population_bound=6.0,
        parameter_fingerprint_sha256=fingerprint,
        mode_classification="brittle",
        certified=True,
    )


def valid_cases():
    return [
        row("ceramic", 300.0, 14.0, "ceramic-fp"),
        row("ceramic", 700.0, 14.4, "ceramic-fp"),
        row("ceramic", 1100.0, 14.8, "ceramic-fp"),
        row("weakT", 300.0, 15.0, "weak-fp"),
        row("weakT", 700.0, 15.8, "weak-fp"),
        row("weakT", 1100.0, 16.5, "weak-fp"),
        row("DBTT", 300.0, 14.0, "dbtt-fp"),
        row("DBTT", 700.0, 20.0, "dbtt-fp"),
        row("DBTT", 1100.0, 28.0, "dbtt-fp"),
    ]


def test_three_parameterizations_are_certified_when_class_responses_separate():
    report = certify(
        valid_cases(),
        {"ceramic", "weakT", "DBTT"},
        [300.0, 700.0, 1100.0],
        0.35,
        0.35,
        1.5,
        2.0,
        0.15,
    )
    assert report["certified"] is True
    assert report["short_growth_authorized"] is True
    assert report["parameter_fingerprints"]["DBTT"] == "dbtt-fp"


def test_matrix_rejects_a_dbtt_parameterization_without_transition_contrast():
    cases = valid_cases()
    cases[-3:] = [
        row("DBTT", 300.0, 16.0, "dbtt-fp"),
        row("DBTT", 700.0, 16.2, "dbtt-fp"),
        row("DBTT", 1100.0, 16.4, "dbtt-fp"),
    ]
    report = certify(
        cases,
        {"ceramic", "weakT", "DBTT"},
        [300.0, 700.0, 1100.0],
        0.35,
        0.35,
        1.5,
        2.0,
        0.15,
    )
    assert report["certified"] is False
    assert "DBTT_high_low_ratio_too_small" in report["failures"]
