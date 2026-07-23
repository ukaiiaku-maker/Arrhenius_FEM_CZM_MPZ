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


class FixedPTRates:
    def __init__(self, n_bins: int, *, peierls_s: float, release_s: float, jump_m: float):
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


def make_state(*, encounter_efficiency: float = 0.0):
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
        source_path="synthetic_transport_kernel",
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


def test_stiff_macro_transport_completes_without_cfl_microsteps():
    state = make_state(encounter_efficiency=0.0)
    state.mobile_positive[0, 0] = 1.0
    state._pt_model = lambda: FixedPTRates(
        state.n_bins,
        peierls_s=1.0e12,
        release_s=0.0,
        jump_m=1.0e-6,
    )
    state.max_transport_substeps = 120

    initial = total_active_content(state)
    result = state.transport(dt_s=840.0, T_K=700.0, opening_stress_Pa=1.0e9)

    assert result["transport_integrator"] == (
        "adaptive_backward_euler_upwind_v10_0_5_14_2"
    )
    assert result["transport_cfl_limited"] is False
    assert result["max_frozen_courant"] > 1.0e12
    assert result["transport_attempted_linear_solves"] <= 120
    assert result["dN_escaped"] > 0.999999 * initial
    assert total_active_content(state) < 1.0e-6 * initial
    assert total_active_content(state) + result["dN_escaped"] == pytest.approx(
        initial, rel=1.0e-8, abs=1.0e-12
    )


def test_stiff_release_and_escape_remain_nonnegative_and_conservative():
    state = make_state(encounter_efficiency=0.0)
    state.retained_negative[1, 0] = 3.0
    state._pt_model = lambda: FixedPTRates(
        state.n_bins,
        peierls_s=1.0e10,
        release_s=1.0e8,
        jump_m=1.0e-6,
    )
    state.max_transport_substeps = 240

    initial = total_active_content(state)
    result = state.transport(dt_s=840.0, T_K=700.0, opening_stress_Pa=1.0e9)

    assert result["dN_detrapped"] > 0.0
    assert result["dN_escaped"] > 0.999 * initial
    assert np.min(state.mobile_negative) >= 0.0
    assert np.min(state.retained_negative) >= 0.0
    assert total_active_content(state) + result["dN_escaped"] == pytest.approx(
        initial, rel=1.0e-8, abs=1.0e-12
    )


def test_small_courant_limit_matches_first_order_upwind():
    state = make_state(encounter_efficiency=0.0)
    state.mobile_positive[0, 0] = 1.0
    velocity = 1.0e-9
    state._pt_model = lambda: FixedPTRates(
        state.n_bins,
        peierls_s=1.0,
        release_s=0.0,
        jump_m=velocity,
    )
    dt = 1.0e-3
    courant = velocity * dt / state.dx
    assert courant < 1.0e-5

    result = state.transport(dt_s=dt, T_K=700.0, opening_stress_Pa=1.0e9)

    expected_cell0 = 1.0 - courant
    expected_cell1 = courant
    assert state.mobile_positive[0, 0] == pytest.approx(
        expected_cell0, rel=1.0e-8, abs=1.0e-12
    )
    assert state.mobile_positive[0, 1] == pytest.approx(
        expected_cell1, rel=1.0e-8, abs=1.0e-12
    )
    assert result["dN_escaped"] == pytest.approx(0.0, abs=1.0e-12)


def test_zero_content_transport_is_a_noop():
    state = make_state(encounter_efficiency=0.0)
    result = state.transport(dt_s=840.0, T_K=700.0, opening_stress_Pa=1.0e9)
    assert result["transport_substeps"] == 0
    assert result["dN_escaped"] == 0.0
    assert state.time_s == 0.0
