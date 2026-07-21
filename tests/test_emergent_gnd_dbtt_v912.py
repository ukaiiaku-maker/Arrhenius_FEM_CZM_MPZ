import math

import numpy as np

from arrhenius_fracture.emergent_gnd_dbtt_v912 import (
    CandidateParameters,
    CommonPhysics,
    EmergentGNDState,
    ExpFloorSurface,
    PTMechanism,
    ProtocolSegment,
    developed_delta_K,
    run_temperature_protocol,
    score_microstructural_transition,
)


def candidate() -> CandidateParameters:
    cleavage = ExpFloorSurface(
        G00_eV=2.0,
        gT_eV_per_K=0.0,
        sigc0_Pa=4.0e9,
        sT_Pa_per_K=0.0,
        exp_a=1.0,
        exp_n=1.0,
        floor_fraction=0.05,
    )
    emission = ExpFloorSurface(
        G00_eV=1.2,
        gT_eV_per_K=0.0,
        sigc0_Pa=3.0e9,
        sT_Pa_per_K=0.0,
        exp_a=1.0,
        exp_n=1.0,
        floor_fraction=0.05,
    )
    return CandidateParameters(
        candidate_id="test",
        cleavage=cleavage,
        emission=emission,
        peierls=PTMechanism(0.5, 0.0, 1.0, 1.0, 1.0e8),
        taylor=PTMechanism(0.8, 0.0, 1.0, 1.0, 1.0e7),
        rho_source0_m2=1.0e13,
        source_refresh_length_m=5.0e-6,
        taylor_corr_rho_c_m2=1.0e14,
        taylor_corr_scale=1.0,
    )


def test_zero_state_has_zero_gnd_feedback():
    state = EmergentGNDState(candidate(), CommonPhysics(n_bins=16))
    assert np.allclose(state.tau_gnd_Pa(), 0.0)
    assert state.K_shield_MPa_sqrt_m() == 0.0


def test_equal_opposite_sign_retained_content_cancels():
    state = EmergentGNDState(candidate(), CommonPhysics(n_bins=16))
    profile = np.linspace(1.0e12, 2.0e12, state.c.n_bins)
    state.retained_m2[:, 0, :] = profile
    state.retained_m2[:, 1, :] = profile
    assert np.allclose(state.signed_gnd_m2(), 0.0)
    assert np.allclose(state.tau_gnd_Pa(), 0.0)
    assert abs(state.K_shield_MPa_sqrt_m()) < 1.0e-14


def test_physical_source_inventory_is_grid_invariant():
    c = candidate()
    p40 = CommonPhysics(n_bins=40)
    p80 = CommonPhysics(n_bins=80)
    s40 = EmergentGNDState(c, p40)
    s80 = EmergentGNDState(c, p80)
    total40 = np.sum(s40.source_capacity_m2) * s40.cell_area_m2
    total80 = np.sum(s80.source_capacity_m2) * s80.cell_area_m2
    assert math.isclose(total40, total80, rel_tol=1.0e-12)


def test_homogeneous_0d_translation_preserves_state_and_refreshes_sources():
    state = EmergentGNDState(candidate(), CommonPhysics(n_bins=1))
    state.retained_m2[0, 1, 0] = 2.5e12
    state.source_available_m2[:] = 0.25 * state.source_capacity_m2
    retained_before = state.retained_m2.copy()
    source_before = state.source_available_m2.copy()

    state.translate_tip(1.0e-6)

    assert np.array_equal(state.retained_m2, retained_before)
    assert np.all(state.source_available_m2 > source_before)
    assert np.all(state.source_available_m2 <= state.source_capacity_m2)
    assert math.isclose(state.extension_m, 1.0e-6)


def test_objective_uses_only_developed_microstructural_increment():
    T = [300, 400, 500, 600, 700]
    delta = [0.0, 0.1, 0.2, 9.0, 9.2]
    score = score_microstructural_transition(T, delta)
    assert score["pass"]
    assert score["largest_jump_localization"] > 0.9


def test_short_protocol_smoke():
    physics = CommonPhysics(
        n_bins=8,
        mpz_length_m=8.0e-6,
        active_strip_width_m=2.0e-6,
        source_zone_length_m=1.0e-6,
        max_fractional_state_change=0.2,
    )
    protocol = [
        ProtocolSegment(0.0, 1.0e-6, 2.0, 3.0, 1.0e-5),
        ProtocolSegment(1.0e-6, 2.0e-6, 3.0, 4.0, 1.0e-5),
    ]
    result = run_temperature_protocol(
        candidate(), physics, protocol, 700.0, target_cleavage_rate_s=1.0e-6
    )
    assert len(result.extensions_um) == 2
    assert np.all(np.isfinite(result.delta_K_micro_MPa_sqrt_m))
    assert developed_delta_K(result, (1.0, 2.0)) >= -1.0e-8
