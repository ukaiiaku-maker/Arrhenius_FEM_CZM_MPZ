from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from arrhenius_fracture.persistent_site_registry_v100514 import (
    select_persistent_site_row,
)
from arrhenius_fracture.persistent_site_signed_mpz_v100514 import (
    PersistentSiteSignedMPZStateV100514,
    SignedShieldingKernelV100514,
)
from arrhenius_fracture.persistent_site_transport_v1005144 import (
    TRANSPORT_INTEGRATOR,
    _advect_mobile_exact,
    _exact_exchange_pair,
    installed_split_transport_v1005144,
)


class FixedPTRates:
    def __init__(
        self,
        n_bins: int,
        *,
        peierls_s: float,
        release_s: float,
        jump_m: float,
    ):
        self.n_bins = int(n_bins)
        self.peierls_s = float(peierls_s)
        self.release_s = float(release_s)
        self.jump_m = float(jump_m)

    def rates(self, sigma, forest, T_K, b_m):
        n = self.n_bins
        return {
            "peierls_rate_s": np.full(n, self.peierls_s),
            "taylor_completion_rate_s": np.full(n, self.release_s),
            "jump_length_m": np.full(n, self.jump_m),
        }


def make_state(*, encounter_efficiency: float = 9.160246308716648):
    candidate = select_persistent_site_row("v912_peak_0118_persistent_sites")
    candidate = replace(candidate, encounter_efficiency=encounter_efficiency)
    n = candidate.n_bins_recommended
    active = np.zeros((2, n))
    kernel = SignedShieldingKernelV100514(
        active_kernel_Pa_sqrt_m_per_signed_line=active,
        wake_kernel_Pa_sqrt_m_per_signed_line=np.zeros_like(active),
        activation_to_line_content_by_system=np.ones(2),
        metadata={
            "candidate_independent": True,
            "counts_are_signed_burgers_lines": True,
            "normalization_is_mechanically_derived": True,
        },
        source_path="synthetic_transport_kernel_v1005144",
    )
    return PersistentSiteSignedMPZStateV100514(
        candidate,
        kernel,
        G_Pa=160.0e9,
        nu=0.28,
        b_m=2.74e-10,
        r0_m=1.0e-6,
        blunting_length_m=0.5e-6,
    )


def total_active_content(state) -> float:
    return float(
        np.sum(state.mobile_positive)
        + np.sum(state.mobile_negative)
        + np.sum(state.retained_positive)
        + np.sum(state.retained_negative)
    )


def test_exact_exchange_is_nonnegative_and_conservative_at_extreme_rates():
    mobile = np.array([[5.0e-6, 0.0], [0.0, 2.0e-7]])
    retained = np.array([[0.0, 3.0e-7], [1.0e-7, 0.0]])
    before = float(np.sum(mobile) + np.sum(retained))
    m1, r1, trapped, released = _exact_exchange_pair(
        mobile,
        retained,
        encounter_rate_s=np.array([1.0e15, 1.0e12]),
        release_rate_s=np.array([1.0e14, 1.0e13]),
        dt_s=420.0,
    )
    assert np.min(m1) >= 0.0
    assert np.min(r1) >= 0.0
    assert float(np.sum(m1) + np.sum(r1)) == pytest.approx(
        before, rel=2.0e-14, abs=1.0e-20
    )
    assert trapped >= 0.0
    assert released >= 0.0


def test_absorbing_advection_has_exact_semigroup_and_mass_balance():
    mobile = np.zeros((4, 80))
    mobile[0, 0] = 5.0e-6
    mobile[1, 5] = 2.0e-7
    velocity = np.linspace(1.0e-7, 5.0e-4, 80)
    full, escape_full, _ = _advect_mobile_exact(
        mobile, velocity, 0.625e-6, 840.0
    )
    half, escape_a, _ = _advect_mobile_exact(
        mobile, velocity, 0.625e-6, 420.0
    )
    two_half, escape_b, _ = _advect_mobile_exact(
        half, velocity, 0.625e-6, 420.0
    )
    assert np.allclose(full, two_half, rtol=2.0e-11, atol=1.0e-20)
    assert escape_full == pytest.approx(
        escape_a + escape_b, rel=2.0e-11, abs=1.0e-20
    )
    initial = float(np.sum(mobile))
    assert float(np.sum(full)) + escape_full == pytest.approx(
        initial, rel=2.0e-11, abs=1.0e-20
    )


def test_tiny_stiff_population_finishes_without_augmented_matrix_error():
    state = make_state(encounter_efficiency=20.0)
    state.mobile_positive[0, 0] = 5.092917e-6
    state.retained_negative[1, 2] = 2.984848e-7
    state._pt_model = lambda: FixedPTRates(
        state.n_bins,
        peierls_s=1.0e13,
        release_s=2.0e12,
        jump_m=1.0e-6,
    )
    state.max_transport_substeps = 300
    initial = total_active_content(state)
    with installed_split_transport_v1005144():
        result = state.transport(
            dt_s=840.0, T_K=1000.0, opening_stress_Pa=4.0e9
        )
    assert result["transport_integrator"] == TRANSPORT_INTEGRATOR
    assert result["transport_attempted_physical_solves"] <= 300
    assert np.min(state.mobile_positive) >= 0.0
    assert np.min(state.mobile_negative) >= 0.0
    assert np.min(state.retained_positive) >= 0.0
    assert np.min(state.retained_negative) >= 0.0
    assert total_active_content(state) + result["dN_escaped"] == pytest.approx(
        initial, rel=2.0e-9, abs=1.0e-18
    )


@pytest.mark.parametrize("temperature_K", [900.0, 1000.0, 1100.0, 1200.0])
def test_candidate_0118_real_high_temperature_rates_conserve_tiny_content(
    temperature_K,
):
    state = make_state()
    state.mobile_positive[0, :4] = 5.0e-7
    state.mobile_negative[1, :4] = 7.5e-8
    state.retained_positive[0, 4:8] = 2.5e-8
    state.max_transport_substeps = 2000
    initial = total_active_content(state)
    with installed_split_transport_v1005144():
        result = state.transport(
            dt_s=840.0,
            T_K=temperature_K,
            opening_stress_Pa=4.0e9,
        )
    assert result["transport_integrator"] == TRANSPORT_INTEGRATOR
    assert result["transport_attempted_physical_solves"] <= 2000
    assert np.isfinite(result["transport_nonlinear_error_max"])
    assert result.get("physical_generator_mass_gain", 0.0) <= (
        2.0e-9 * initial + 1.0e-20
    )
    final = total_active_content(state)
    assert final + result["dN_escaped"] == pytest.approx(
        initial, rel=2.0e-8, abs=1.0e-17
    )


def test_zero_content_advances_clock_without_transport_work():
    state = make_state()
    with installed_split_transport_v1005144():
        result = state.transport(
            dt_s=840.0, T_K=1200.0, opening_stress_Pa=4.0e9
        )
    assert result["transport_substeps"] == 0
    assert result["transport_attempted_physical_solves"] == 0
    assert state.time_s == pytest.approx(840.0)
