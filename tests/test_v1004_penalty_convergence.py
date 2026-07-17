from analyze_v10_0_4_penalty_convergence import PenaltyCase, certify


def case(case_id, normal, tangent, Kc, U):
    return PenaltyCase(
        case_id=case_id,
        normal_penalty=normal,
        tangent_penalty=tangent,
        material="weakT",
        temperature_K=700.0,
        Kc_MPa_sqrt_m=Kc,
        event_Uapp_m=U,
        event_Ftop_N=100.0 * U,
        event_sigma_tip_Pa=6.0e9,
        max_N_em=3.0,
        source_population_bound=5.7,
        certified=True,
        quality_veto_count=0,
        full_rollbacks=0,
        committed_events=1,
    )


def converged_cases():
    return [
        case("n5", 5e17, 1e18, 15.90, 7.48e-5),
        case("base", 1e18, 1e18, 15.94, 7.50e-5),
        case("n2", 2e18, 1e18, 15.96, 7.51e-5),
        case("t5", 1e18, 5e17, 15.93, 7.50e-5),
        case("t2", 1e18, 2e18, 15.95, 7.50e-5),
    ]


def test_penalty_matrix_certifies_small_normal_and_tangent_sensitivity():
    report = certify(converged_cases(), 1e18, 0.01, 0.01, 0.0025, 0.005)
    assert report["certified"] is True
    assert report["parameterization_matrix_authorized"] is True
    assert report["recommended_normal_penalty_Pa_per_m"] == 1e18


def test_penalty_matrix_rejects_large_normal_stiffness_dependence():
    cases = converged_cases()
    cases[0] = case("n5", 5e17, 1e18, 14.0, 6.5e-5)
    report = certify(cases, 1e18, 0.01, 0.01, 0.0025, 0.005)
    assert report["certified"] is False
    assert "normal_Kc_relative_span" in report["failures"]
