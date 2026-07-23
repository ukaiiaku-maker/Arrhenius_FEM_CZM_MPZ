from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from arrhenius_fracture.anisotropic_two_channel_drive_v100514 import (
    resolve_two_channel_drive,
    resolve_two_channel_drive_from_tensors,
)
from arrhenius_fracture.persistent_site_registry_v100514 import (
    ROWS,
    select_persistent_site_row,
)
from arrhenius_fracture.persistent_site_signed_mpz_v100514 import (
    PersistentSiteSignedMPZStateV100514,
    SignedShieldingKernelV100514,
    effective_front_width_m,
    solve_backstress_limited_activations,
)


def kernel(n_bins=80):
    x = np.arange(n_bins, dtype=float)
    active = np.vstack(
        (1.0e3 / np.sqrt(x + 1.0), 0.8e3 / np.sqrt(x + 1.0))
    )
    return SignedShieldingKernelV100514(
        active_kernel_Pa_sqrt_m_per_signed_line=active,
        wake_kernel_Pa_sqrt_m_per_signed_line=np.zeros_like(active),
        activation_to_line_content_by_system=np.array([1.0, 1.0]),
        metadata={
            "candidate_independent": True,
            "counts_are_signed_burgers_lines": True,
            "normalization_is_mechanically_derived": True,
        },
        source_path="synthetic_test_kernel",
    )


def state(candidate=None):
    candidate = candidate or select_persistent_site_row(
        "v912_peak_0118_persistent_sites"
    )
    return PersistentSiteSignedMPZStateV100514(
        candidate,
        kernel(candidate.n_bins_recommended),
        G_Pa=160.0e9,
        nu=0.28,
        b_m=2.74e-10,
        r0_m=1.0e-6,
        blunting_length_m=0.5e-6,
    )


def test_exact_registry_contract_and_candidate_0118():
    assert len(ROWS) == 5
    candidate = select_persistent_site_row("v912_peak_0118_persistent_sites")
    assert candidate.candidate_id == "v912_targeted_local_peak_005518_0118"
    assert candidate.emit_G00_eV == pytest.approx(4.654189503333349)
    assert candidate.rho_source0_m2 == pytest.approx(2.3012695321899512e16)
    assert candidate.source_recovery_rate_s == 0.0
    assert candidate.retained_recovery_rate_s == 0.0
    assert candidate.recovery_nu0_s == 0.0


def test_persistent_source_invariance_after_emission():
    model = state()
    before = model.available_sites.copy()
    result = model.emit_persistent(
        dt_s=1.0e-6,
        T_K=900.0,
        opening_stress_Pa=4.0e9,
        drive_factors=np.array([0.5, 0.5]),
        tau_signed_Pa=np.array([1.0e9, -1.0e9]),
        rate_function=lambda stress, T: 1.0e3,
    )
    assert result["dN_emit"] > 0.0
    assert np.sum(model.mobile_positive) > 0.0
    assert np.sum(model.mobile_negative) > 0.0
    assert np.allclose(model.available_sites, before)
    assert model.available_site_fraction == 1.0
    assert result["source_depletion_active"] is False


def test_grid_independent_width_and_multiplicity():
    candidate = select_persistent_site_row("v912_peak_0118_persistent_sites")
    state80 = state(candidate)
    candidate160 = replace(candidate, n_bins_recommended=160)
    state160 = PersistentSiteSignedMPZStateV100514(
        candidate160,
        kernel(160),
        G_Pa=160.0e9,
        nu=0.28,
        b_m=2.74e-10,
        r0_m=1.0e-6,
        blunting_length_m=0.5e-6,
    )
    geometry80 = state80.source_geometry()
    geometry160 = state160.source_geometry()
    assert geometry80["front_width_m"] == pytest.approx(
        geometry160["front_width_m"], rel=1e-14
    )
    assert geometry80["multiplicity_per_system"] == pytest.approx(
        geometry160["multiplicity_per_system"], rel=1e-14
    )
    width = effective_front_width_m(
        1.0e30,
        reference_width_m=10e-6,
        reference_density_m2=5e12,
        minimum_physical_width_m=0.0,
        burgers_m=2.74e-10,
        maximum_width_m=50e-6,
    )
    assert width >= 2.74e-10


def test_interior_implicit_root():
    value, block, blocked = solve_backstress_limited_activations(
        multiplicity=10.0,
        dt_s=1.0,
        drive_stress_Pa=100.0,
        rho_initial_m2=0.0,
        rho_increment_per_activation_m2=1.0,
        backstress_prefactor_Pa_sqrt_m2=1.0,
        rate_function=lambda stress: 0.2,
    )
    assert block > 2.0
    assert blocked is False
    assert value == pytest.approx(2.0, rel=1e-9)


def test_mechanical_blocking_complementarity():
    value, block, blocked = solve_backstress_limited_activations(
        multiplicity=1.0e12,
        dt_s=1.0,
        drive_stress_Pa=10.0,
        rho_initial_m2=0.0,
        rho_increment_per_activation_m2=1.0,
        backstress_prefactor_Pa_sqrt_m2=1.0,
        rate_function=lambda stress: 1.0e12,
    )
    assert blocked is True
    assert value == pytest.approx(block)
    assert block == pytest.approx(100.0)


def test_zero_drive_gate():
    value, block, blocked = solve_backstress_limited_activations(
        multiplicity=1.0e6,
        dt_s=1.0,
        drive_stress_Pa=1.0,
        rho_initial_m2=4.0,
        rho_increment_per_activation_m2=1.0,
        backstress_prefactor_Pa_sqrt_m2=1.0,
        rate_function=lambda stress: 1.0e6,
    )
    assert value == 0.0
    assert block == 0.0
    assert blocked is False


def test_burgers_sign_reversal_equal_magnitude():
    positive = state()
    negative = state()
    args = dict(
        dt_s=1.0e-7,
        T_K=900.0,
        opening_stress_Pa=3.0e9,
        drive_factors=np.array([0.5, 0.0]),
        rate_function=lambda stress, T: 100.0,
    )
    positive.emit_persistent(tau_signed_Pa=np.array([1.0, 0.0]), **args)
    negative.emit_persistent(tau_signed_Pa=np.array([-1.0, 0.0]), **args)
    assert np.sum(positive.mobile_positive) == pytest.approx(
        np.sum(negative.mobile_negative)
    )
    assert np.sum(positive.mobile_negative) == 0.0
    assert np.sum(negative.mobile_positive) == 0.0


def test_crack_advance_resharpens_without_refresh():
    model = state()
    model.accumulated_slip_positive[0, :4] = 10.0
    radius_before = model.blunted_radius()
    available_before = model.available_sites.copy()
    result = model.advance(2.0e-6)
    assert result["wake_slip"] > 0.0
    assert result["tip_radius_after_advance_m"] < radius_before
    assert result["source_sites_refreshed"] == 0.0
    assert np.allclose(model.available_sites, available_before)
    assert model.available_site_fraction == 1.0


def test_unsigned_backstress_survives_signed_shielding_cancellation():
    model = state()
    model.retained_positive[0, :3] = 4.0
    model.retained_negative[0, :3] = 4.0
    _, backstress = model.backstress()
    assert backstress[0] > 0.0
    assert model.shielding_K() == pytest.approx(0.0, abs=1e-14)


def test_trial_copy_does_not_mutate_committed_state():
    committed = state()
    trial = committed.copy()
    trial.emit_persistent(
        dt_s=1.0e-6,
        T_K=900.0,
        opening_stress_Pa=4.0e9,
        drive_factors=np.array([0.5, 0.5]),
        tau_signed_Pa=np.array([1.0, -1.0]),
        rate_function=lambda stress, T: 1.0e3,
    )
    assert committed.mobile_count == 0.0
    assert trial.mobile_count > 0.0


def test_two_channel_tensor_drive_has_signed_channels():
    stress = np.array([[2.0e9, 0.4e9], [0.4e9, 3.0e9]])
    result = resolve_two_channel_drive(stress, 3.0e9, 45.0)
    assert result["two_channel_drive_reliable"] is True
    assert len(result["two_channel_tau_signed_Pa"]) == 2
    assert len(result["two_channel_drive_factors"]) == 2
    assert all(value >= 0.0 for value in result["two_channel_drive_factors"])


def test_restart_round_trip_preserves_signed_state_and_full_sources():
    original = state()
    original.mobile_positive[0, 0] = 3.0
    original.retained_negative[1, 2] = 4.0
    original.accumulated_slip_positive[0, 1] = 5.0
    original.advance_total_m = 2.0e-6
    original.time_s = 7.0
    restored = PersistentSiteSignedMPZStateV100514.from_state_dict(
        original.state_dict(), kernel(original.n_bins)
    )
    assert np.array_equal(restored.mobile_positive, original.mobile_positive)
    assert np.array_equal(restored.retained_negative, original.retained_negative)
    assert np.array_equal(
        restored.accumulated_slip_positive,
        original.accumulated_slip_positive,
    )
    assert restored.advance_total_m == pytest.approx(original.advance_total_m)
    assert restored.time_s == pytest.approx(original.time_s)
    assert restored.available_site_fraction == 1.0
    assert np.allclose(restored.available_sites, restored.site_capacity)


def test_two_channel_resolution_uses_distinct_channel_tensors():
    opening = np.array([[1.0e9, 0.0], [0.0, 3.0e9]])
    channel_a = np.array([[2.0e9, 0.8e9], [0.8e9, 1.0e9]])
    channel_b = np.array([[1.0e9, -0.3e9], [-0.3e9, 2.0e9]])
    result = resolve_two_channel_drive_from_tensors(
        opening, [channel_a, channel_b], [1.0, 0.0], 45.0
    )
    tau = result["two_channel_tau_signed_Pa"]
    assert len(tau) == 2
    assert tau[0] != pytest.approx(tau[1])
    assert len(result["two_channel_drive_factors"]) == 2
