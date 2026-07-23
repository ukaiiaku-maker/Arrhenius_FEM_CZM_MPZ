from dataclasses import asdict

import pytest

from arrhenius_fracture.emergent_gnd_rcurve_v913 import (
    RCurveEvent,
    RCurveLoadingMap,
    RCurveResult,
    run_autonomous_rcurve,
)
from arrhenius_fracture.emergent_gnd_types_v912 import (
    ExpFloorSurface,
    PTMechanism,
)
from arrhenius_fracture.emergent_gnd_types_v913 import (
    CandidateParameters,
    CommonPhysics,
)
from scripts.calibrate_v913_rcurve_to_v10222_top5 import _select_targets


def candidate() -> CandidateParameters:
    surface = ExpFloorSurface(
        G00_eV=0.4,
        gT_eV_per_K=0.0,
        sigc0_Pa=2.0e9,
        sT_Pa_per_K=0.0,
        exp_a=1.0,
        exp_n=1.0,
        floor_fraction=0.01,
    )
    return CandidateParameters(
        candidate_id="autonomous-test",
        cleavage=surface,
        emission=surface,
        peierls=PTMechanism(0.05, 0.0, 1.0, 1.0, 1.0e8),
        taylor=PTMechanism(0.10, 0.0, 1.0, 1.0, 1.0e8),
        rho_source0_m2=1.4e16,
        source_refresh_length_m=1.0e-6,
        taylor_corr_rho_c_m2=1.0e14,
        taylor_corr_scale=1.0,
        recovery_nu0_s=0.0,
        c_blunt=1.4,
    )


def physics() -> CommonPhysics:
    result = CommonPhysics(
        n_bins=8,
        mpz_length_m=8.0e-6,
        source_zone_length_m=2.0e-6,
        blunting_length_m=1.0e-6,
        activation_to_line_content_per_system=(1.0, 1.0),
    )
    result.validate()
    return result


def loading_map() -> RCurveLoadingMap:
    return RCurveLoadingMap(
        K_per_U_MPa_sqrt_m_per_m=(1.0e8, 1.0e8),
        threshold_actions=(0.2, 0.3),
        path_advances_m=(1.0e-6, 1.5e-6),
        projected_advances_m=(1.0e-6, 1.0e-6),
        nominal_dU_m=1.0e-8,
        nominal_dt_s=1.0e-7,
        seed=3621,
        reference_candidate_id="plastic-free-reference",
        reference_temperature_K=300.0,
    )


def test_loading_map_round_trip_and_validation():
    original = loading_map()
    restored = RCurveLoadingMap.from_dict(original.as_dict())
    assert restored == original

    invalid = original.as_dict()
    invalid["threshold_actions"] = [0.2]
    with pytest.raises(ValueError, match="equal nonzero length"):
        RCurveLoadingMap.from_dict(invalid)


def test_checkpoint_uses_first_accepted_event_at_or_beyond_extension():
    result = RCurveResult(
        candidate_id="test",
        temperature_K=900.0,
        status="complete",
        seed=3621,
        target_projected_extension_m=10.0e-6,
        events=[
            RCurveEvent(
                event_index=index,
                threshold_action=1.0,
                applied_displacement_m=1.0e-6,
                elapsed_time_s=1.0,
                K_MPa_sqrt_m=K,
                path_advance_m=extension,
                projected_advance_m=extension,
                cumulative_path_extension_m=extension,
                cumulative_projected_extension_m=extension,
                tip_radius_pre_advance_m=1.0e-6,
                tip_radius_post_advance_m=1.0e-6,
                front_width_pre_advance_m=1.0e-6,
                backstress_pre_advance_Pa=0.0,
                source_multiplicity_pre_advance=1.0,
                cumulative_source_activations=0.0,
                cumulative_line_content=0.0,
                integration_substeps=1,
            )
            for index, (extension, K) in enumerate(
                ((2.0e-6, 20.0), (7.0e-6, 30.0), (12.0e-6, 40.0))
            )
        ],
    )
    assert result.checkpoint_K(0.0) == 20.0
    assert result.checkpoint_K(7.0e-6) == 30.0
    assert result.checkpoint_K(10.0e-6) == 40.0


def test_autonomous_driver_completes_without_changing_candidate_parameters():
    fixed_candidate = candidate()
    original = asdict(fixed_candidate)
    result = run_autonomous_rcurve(
        fixed_candidate,
        physics(),
        loading_map(),
        900.0,
        target_projected_extension_m=2.0e-6,
        maximum_integration_substeps=10_000,
        translation_mode="hazard_coupled",
        translation_action_exponent=0.95,
    )
    assert result.status == "complete"
    assert len(result.events) == 2
    assert result.achieved_projected_extension_m == pytest.approx(2.0e-6)
    assert result.events[-1].cumulative_path_extension_m == pytest.approx(2.5e-6)
    assert (
        result.numerical_integration["candidate_parameters_modified_by_driver"] is False
    )
    assert result.numerical_integration["translation_action_exponent"] == 0.95
    assert asdict(fixed_candidate) == original


def test_exact_case_selection_does_not_create_a_cross_product():
    targets = [
        {"candidate_id": candidate_id, "temperature_K": temperature}
        for candidate_id in ("a", "b")
        for temperature in (800.0, 900.0)
    ]
    selected = _select_targets(
        targets,
        candidate_ids=(),
        temperatures=(),
        case_specs=("a:800", "b:900"),
    )
    assert [(row["candidate_id"], row["temperature_K"]) for row in selected] == [
        ("a", 800.0),
        ("b", 900.0),
    ]

    with pytest.raises(ValueError, match="CANDIDATE_ID:T_K"):
        _select_targets(
            targets,
            candidate_ids=(),
            temperatures=(),
            case_specs=("not-a-case",),
        )
