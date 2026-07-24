from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from arrhenius_fracture.persistent_site_pf_update_v100515 import (
    PF_UPDATE_MAP,
    exact_exchange,
    evolve_pf_v10222,
    fractional_advect_forward,
)
from arrhenius_fracture.persistent_site_registry_v100514 import (
    select_persistent_site_row,
)
from arrhenius_fracture.persistent_site_signed_mpz_v100514 import (
    PersistentSiteSignedMPZStateV100514,
    SignedShieldingKernelV100514,
)


class FixedRates:
    def __init__(self, n_bins, peierls, release, jump):
        self.n_bins = int(n_bins)
        self.peierls = np.asarray(peierls, dtype=float)
        self.release = np.asarray(release, dtype=float)
        self.jump = np.asarray(jump, dtype=float)

    def rates(self, stress, forest, T_K, b_m):
        return {
            "peierls_rate_s": np.broadcast_to(self.peierls, (self.n_bins,)).copy(),
            "taylor_completion_rate_s": np.broadcast_to(
                self.release, (self.n_bins,)
            ).copy(),
            "jump_length_m": np.broadcast_to(self.jump, (self.n_bins,)).copy(),
        }


def make_state(*, encounter_efficiency=1.0):
    candidate = select_persistent_site_row("v912_peak_0118_persistent_sites")
    candidate = replace(candidate, encounter_efficiency=float(encounter_efficiency))
    n = candidate.n_bins_recommended
    kernel = SignedShieldingKernelV100514(
        active_kernel_Pa_sqrt_m_per_signed_line=np.zeros((2, n)),
        wake_kernel_Pa_sqrt_m_per_signed_line=np.zeros((2, n)),
        activation_to_line_content_by_system=np.ones(2),
        metadata={
            "candidate_independent": True,
            "counts_are_signed_burgers_lines": True,
            "normalization_is_mechanically_derived": True,
        },
        source_path="synthetic_pf_update_map",
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


def inventory(state):
    return float(
        np.sum(state.mobile_positive)
        + np.sum(state.mobile_negative)
        + np.sum(state.retained_positive)
        + np.sum(state.retained_negative)
        + np.sum(state.wake_mobile_positive)
        + np.sum(state.wake_mobile_negative)
        + np.sum(state.wake_retained_positive)
        + np.sum(state.wake_retained_negative)
        + state.escaped_total
        + state.wake_discarded_total
    )


def test_exact_exchange_matches_pf_closed_form():
    mobile = np.array([[3.0, 1.0], [0.5, 0.0]])
    retained = np.array([[1.0, 2.0], [0.5, 4.0]])
    ke = np.array([2.0, 0.5])
    kr = np.array([1.0, 1.5])
    dt = 0.7
    total = mobile + retained
    rate = ke[None, :] + kr[None, :]
    req = ke[None, :] / rate * total
    expected_r = req + (retained - req) * np.exp(-rate * dt)
    expected_m = total - expected_r
    new_m, new_r, trapped, released = exact_exchange(
        mobile, retained, ke, kr, dt
    )
    assert np.allclose(new_m, expected_m, rtol=2e-15, atol=2e-15)
    assert np.allclose(new_r, expected_r, rtol=2e-15, atol=2e-15)
    assert np.sum(new_m + new_r) == pytest.approx(np.sum(total), rel=2e-15)
    assert trapped >= 0.0
    assert released >= 0.0


def test_fractional_advection_matches_pf_cell_remap():
    field = np.zeros((2, 5))
    field[0, 0] = 4.0
    field[1, 3] = 2.0
    moved, escaped = fractional_advect_forward(field, 1.25, 1.0)
    expected = np.zeros_like(field)
    expected[0, 1] = 3.0
    expected[0, 2] = 1.0
    expected[1, 4] = 1.5
    assert np.allclose(moved, expected)
    assert escaped == pytest.approx(0.5)
    assert np.sum(moved) + escaped == pytest.approx(np.sum(field))


def test_complete_map_uses_retained_forest_exact_exchange_and_scalar_velocity():
    state = make_state(encounter_efficiency=0.5)
    state.mobile_positive[0, 0] = 2.0
    state.mobile_negative[1, 2] = 1.0
    state.retained_positive[0, 1] = 3.0
    state.emit_persistent = lambda **kwargs: {
        "dN_emit": 0.0,
        "aggregate_hazard_initial_by_system_s": np.zeros(2),
    }
    peierls = np.linspace(1.0, 2.0, state.n_bins)
    state._pt_model = lambda: FixedRates(
        state.n_bins,
        peierls=peierls,
        release=np.full(state.n_bins, 0.25),
        jump=np.full(state.n_bins, state.dx),
    )
    before = inventory(state)
    result = evolve_pf_v10222(
        state,
        dt_s=0.2,
        T_K=700.0,
        opening_stress_Pa=2.0e9,
        drive_factors=np.array([0.2, 0.1]),
        tau_signed_Pa=np.array([1.0, -1.0]),
        emission_rate_function=lambda stress, T: 0.0,
    )
    assert result["transport_integrator"] == PF_UPDATE_MAP
    assert result["explicit_recovery_active"] is False
    assert result["glide_velocity_m_s"] >= 0.0
    assert inventory(state) == pytest.approx(before, rel=5e-12, abs=1e-14)
    assert result["line_content_balance_relative_error"] < 5e-11


@pytest.mark.parametrize(
    "temperature_K", [300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0, 1100.0, 1200.0]
)
def test_candidate_0118_pf_map_is_finite_and_conservative_all_temperatures(
    temperature_K,
):
    state = make_state()
    state.mobile_positive[0, :4] = 5.0e-7
    state.mobile_negative[1, :4] = 7.5e-8
    state.retained_positive[0, 4:8] = 2.5e-8
    state.emit_persistent = lambda **kwargs: {
        "dN_emit": 0.0,
        "aggregate_hazard_initial_by_system_s": np.zeros(2),
    }
    before = inventory(state)
    result = evolve_pf_v10222(
        state,
        dt_s=840.0,
        T_K=temperature_K,
        opening_stress_Pa=4.0e9,
        drive_factors=np.array([0.2, 0.1]),
        tau_signed_Pa=np.array([1.0, -1.0]),
        emission_rate_function=lambda stress, T: 0.0,
    )
    assert result["transport_integrator"] == PF_UPDATE_MAP
    assert np.isfinite(result["peierls_rate_s"])
    assert np.isfinite(result["taylor_completion_rate_s"])
    assert np.min(state.mobile_positive) >= 0.0
    assert np.min(state.mobile_negative) >= 0.0
    assert np.min(state.retained_positive) >= 0.0
    assert np.min(state.retained_negative) >= 0.0
    assert inventory(state) == pytest.approx(before, rel=5e-11, abs=1e-18)
