from __future__ import annotations

import numpy as np

from arrhenius_fracture.emergent_gnd_rcurve_v913 import RCurveLoadingMap
from arrhenius_fracture.emergent_gnd_types_v912 import ExpFloorSurface, PTMechanism
from arrhenius_fracture.emergent_gnd_types_v913 import CandidateParameters, CommonPhysics
from arrhenius_fracture.zero_d_persistent_v913 import (
    ZeroDRunSettings,
    local_peak_metrics,
    reduction_geometry,
    run_zero_d_rcurve,
)


def candidate() -> CandidateParameters:
    cleavage = ExpFloorSurface(
        G00_eV=0.20,
        gT_eV_per_K=0.0,
        sigc0_Pa=2.0e9,
        sT_Pa_per_K=0.0,
        exp_a=1.0,
        exp_n=1.0,
        floor_fraction=0.05,
    )
    emission = ExpFloorSurface(
        G00_eV=4.0,
        gT_eV_per_K=0.0,
        sigc0_Pa=5.0e9,
        sT_Pa_per_K=0.0,
        exp_a=1.0,
        exp_n=1.0,
        floor_fraction=0.05,
    )
    return CandidateParameters(
        candidate_id="unit",
        cleavage=cleavage,
        emission=emission,
        peierls=PTMechanism(2.0, 0.0, 1.0, 1.0, 1.0e12),
        taylor=PTMechanism(1.0, 0.0, 1.0, 1.0, 1.0e11),
        rho_source0_m2=1.0e15,
        source_refresh_length_m=0.0,
        taylor_corr_rho_c_m2=1.0e15,
        taylor_corr_scale=1.0,
        c_blunt=1.0,
    )


def physics() -> CommonPhysics:
    return CommonPhysics(
        n_bins=16,
        minimum_front_width_m=2.74e-10,
        maximum_front_width_m=50.0e-6,
        shielding_orientation_factors=(0.0, 0.0),
        emission_schmid_factors=(0.1, 0.02),
        activation_to_line_content_per_system=(1.0, 1.0),
        encounter_efficiency=1.0,
        taylor_phi_max=20.0,
        mobile_transport_velocity_scale=0.0,
    )


def loading_map() -> RCurveLoadingMap:
    return RCurveLoadingMap(
        K_per_U_MPa_sqrt_m_per_m=(1.0e8, 1.0e8),
        threshold_actions=(0.01, 0.01),
        path_advances_m=(5.0e-6, 5.0e-6),
        projected_advances_m=(5.0e-6, 5.0e-6),
        nominal_dU_m=1.0e-7,
        nominal_dt_s=1.0,
        seed=1,
        reference_candidate_id="unit",
        reference_temperature_K=1000.0,
    )


def test_reduction_geometry_is_physical() -> None:
    reduced = reduction_geometry(physics())
    assert reduced.source_bins >= 1
    assert reduced.cell_area_m2 > 0.0
    assert np.all(np.asarray(reduced.density_increment_per_activation_m2) > 0.0)
    assert np.all(np.asarray(reduced.slip_count_increment_per_activation) > 0.0)


def test_zero_d_run_uses_persistent_contract() -> None:
    result = run_zero_d_rcurve(
        candidate(),
        physics(),
        loading_map(),
        1000.0,
        settings=ZeroDRunSettings(
            target_projected_extension_m=10.0e-6,
            load_increment_factor=1.0,
            maximum_applied_displacement_m=1.0e-3,
        ),
    )
    assert result.status == "complete"
    assert len(result.events) == 2
    contract = result.numerical_contract
    assert contract["finite_source_inventory"] is False
    assert contract["source_refresh_on_crack_advance"] is False
    assert contract["explicit_recovery"] is False
    assert contract["persistent_multiplicity"] is True
    assert contract["implicit_backstress_limited_emission"] is True


def test_local_peak_metrics_separate_sharpness_and_rebound() -> None:
    metrics = local_peak_metrics(
        [800, 900, 1000, 1100, 1200, 1300],
        [30, 40, 55, 44, 48, 58],
    )
    assert metrics["peak_temperature_K"] == 1000.0
    assert metrics["two_sided_prominence"] == 11.0
    assert metrics["post_peak_drop"] == 11.0
    assert metrics["high_temperature_rebound"] == 3.0
    assert metrics["peak_internal"] is True
