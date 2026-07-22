"""Regression tests for stiff-safe v9.12 Arrhenius state integration."""
import numpy as np

from arrhenius_fracture.emergent_gnd_dbtt_v912 import (
    CandidateParameters,
    CommonPhysics,
    EmergentGNDState,
    ExpFloorSurface,
    PTMechanism,
)


def fast_candidate() -> CandidateParameters:
    cleavage = ExpFloorSurface(2.0, 0.0, 4.0e9, 0.0, 1.0, 1.0, 0.05)
    emission = ExpFloorSurface(0.05, 0.0, 1.0e9, 0.0, 2.0, 1.0, 0.01)
    return CandidateParameters(
        candidate_id="near_prefactor",
        cleavage=cleavage,
        emission=emission,
        peierls=PTMechanism(0.05, 0.0, 2.0, 1.0, 1.0e12),
        taylor=PTMechanism(0.10, 0.0, 2.0, 1.0, 1.0e11),
        rho_source0_m2=1.0e14,
        source_refresh_length_m=5.0e-6,
        taylor_corr_rho_c_m2=1.0e14,
        taylor_corr_scale=1.0,
        recovery_nu0_s=1.0e12,
        recovery_H0_eV=0.05,
    )


def test_near_prefactor_homogeneous_segment_completes():
    state = EmergentGNDState(fast_candidate(), CommonPhysics(n_bins=1))
    totals = state.advance_time(8.4, 20.0, 1200.0)
    assert np.isclose(state.time_s, 8.4)
    assert np.all(np.isfinite(state.mobile_m2))
    assert np.all(np.isfinite(state.retained_m2))
    assert np.all(state.mobile_m2 >= 0.0)
    assert np.all(state.retained_m2 >= 0.0)
    assert totals["emitted_per_m"] >= 0.0


def test_coupled_transport_storage_handles_large_rates():
    physics = CommonPhysics(
        n_bins=16,
        n_systems=2,
        mpz_length_m=16.0e-6,
        source_zone_length_m=2.0e-6,
    )
    state = EmergentGNDState(fast_candidate(), physics)
    state.mobile_m2[0, 1, :4] = 1.0e14
    initial = float(np.sum(state.mobile_m2 + state.retained_m2))

    rates = {
        "velocity_m_s": np.vstack(
            [np.full(16, 1.0e3), np.full(16, 1.0e3)]
        ),
        "encounter_s": np.vstack(
            [np.full(16, 2.0e9), np.full(16, 2.0e9)]
        ),
        "taylor_completion_s": np.vstack(
            [np.full(16, 1.0e8), np.full(16, 1.0e8)]
        ),
        "recovery_rate_s": np.asarray(0.0),
    }
    state._coupled_mobile_retained(rates, 0.1)

    total = float(np.sum(state.mobile_m2 + state.retained_m2))
    assert np.all(np.isfinite(state.mobile_m2))
    assert np.all(np.isfinite(state.retained_m2))
    assert np.all(state.mobile_m2 >= 0.0)
    assert np.all(state.retained_m2 >= 0.0)
    assert total <= initial * (1.0 + 1.0e-12)
    assert float(np.sum(state.retained_m2)) > 0.0


def test_exact_annihilation_preserves_signed_gnd():
    physics = CommonPhysics(
        n_bins=4,
        n_systems=2,
        mpz_length_m=4.0e-6,
        source_zone_length_m=1.0e-6,
    )
    state = EmergentGNDState(fast_candidate(), physics)
    state.retained_m2[0, 0] = np.asarray([4.0, 8.0, 2.0, 10.0]) * 1.0e13
    state.retained_m2[0, 1] = np.asarray([7.0, 3.0, 2.0, 15.0]) * 1.0e13
    signed_before = state.signed_gnd_m2().copy()

    rates = {
        "velocity_m_s": np.vstack(
            [np.full(4, 1.0e3), np.zeros(4)]
        )
    }
    removed = state._annihilate_exact(rates, 0.1)

    assert removed > 0.0
    assert np.allclose(state.signed_gnd_m2(), signed_before)
    assert np.all(state.retained_m2 >= 0.0)


def test_near_prefactor_spatial_segment_completes():
    physics = CommonPhysics(
        n_bins=8,
        mpz_length_m=5.0e-5,
        source_zone_length_m=2.0e-6,
    )
    state = EmergentGNDState(fast_candidate(), physics)
    totals = state.advance_time(8.4, 20.0, 1200.0)

    assert np.isclose(state.time_s, 8.4)
    assert np.all(np.isfinite(state.mobile_m2))
    assert np.all(np.isfinite(state.retained_m2))
    assert np.all(state.mobile_m2 >= 0.0)
    assert np.all(state.retained_m2 >= 0.0)
    assert totals["emitted_per_m"] >= 0.0
    metadata = state.integration_metadata()
    assert metadata["spatial_integrator"] == (
        "coupled_mobile_retained_backward_euler_v2_post_emit_refresh"
    )
    assert metadata["constitutive_feedback_update"] == (
        "refresh_after_midpoint_emission"
    )


def test_second_coupled_half_uses_post_emission_rates(monkeypatch):
    physics = CommonPhysics(
        n_bins=4,
        n_systems=2,
        mpz_length_m=4.0e-6,
        source_zone_length_m=1.0e-6,
    )
    state = EmergentGNDState(fast_candidate(), physics)
    rate_snapshots = []
    coupled_mobile_markers = []

    def fake_local_rates(K_MPa_sqrt_m, T_K):
        del K_MPa_sqrt_m, T_K
        mobile_total = float(np.sum(state.mobile_m2))
        rate_snapshots.append(mobile_total)
        emission = np.zeros(
            (state.c.n_systems, 2, state.c.n_bins),
            dtype=float,
        )
        # The second local-rate evaluation is the midpoint, before emission.
        if len(rate_snapshots) == 2:
            for system, sign in enumerate(state.c.emission_signs):
                q = 1 if sign > 0 else 0
                emission[system, q, 0] = 1.0e6
        marker = np.full(
            (state.c.n_systems, state.c.n_bins),
            mobile_total,
            dtype=float,
        )
        zeros = np.zeros_like(marker)
        return {
            "velocity_m_s": zeros,
            "encounter_s": marker,
            "taylor_completion_s": zeros,
            "emission_rate_s": emission,
            "recovery_rate_s": np.asarray(0.0),
            "tau_external_Pa": zeros,
            "tau_nonlocal_shielding_Pa": zeros,
            "tau_gnd_Pa": zeros,
            "tau_eff_Pa": zeros,
        }

    def fake_coupled(rates, dt):
        del dt
        coupled_mobile_markers.append(float(rates["encounter_s"][0, 0]))

    monkeypatch.setattr(state, "local_rates", fake_local_rates)
    monkeypatch.setattr(state, "_coupled_mobile_retained", fake_coupled)
    monkeypatch.setattr(state, "_annihilate_exact", lambda rates, dt: 0.0)

    state._advance_spatial_step(0.1, 20.0, 900.0)

    assert len(rate_snapshots) == 4
    assert len(coupled_mobile_markers) == 2
    assert coupled_mobile_markers[0] == 0.0
    assert coupled_mobile_markers[1] > 0.0
    assert rate_snapshots[2] > rate_snapshots[1]
