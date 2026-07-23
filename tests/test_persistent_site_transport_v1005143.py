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
from arrhenius_fracture.persistent_site_transport_v1005143 import (
    TRANSPORT_INTEGRATOR,
    installed_exponential_transport_v1005143,
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


def test_frozen_exponential_has_stiff_semigroup_property():
    state = make_state(encounter_efficiency=0.0)
    state.mobile_positive[0, 0] = 1.0
    state._pt_model = lambda: FixedPTRates(
        state.n_bins,
        peierls_s=1.0e12,
        release_s=0.0,
        jump_m=1.0e-6,
    )
    with installed_exponential_transport_v1005143():
        initial = state._transport_snapshot()
        full, full_diag = state._frozen_transport_step(
            initial, dt_s=840.0, T_K=700.0, opening_stress_Pa=1.0e9
        )
        half, first_diag = state._frozen_transport_step(
            initial, dt_s=420.0, T_K=700.0, opening_stress_Pa=1.0e9
        )
        two_half, second_diag = state._frozen_transport_step(
            half, dt_s=420.0, T_K=700.0, opening_stress_Pa=1.0e9
        )
    assert state._snapshot_difference(full, two_half) < 1.0e-11
    assert full_diag["dN_escaped"] == pytest.approx(
        first_diag["dN_escaped"] + second_diag["dN_escaped"],
        rel=1.0e-10,
        abs=1.0e-12,
    )


def test_840_second_stiff_transport_finishes_without_tail_refinement():
    state = make_state(encounter_efficiency=0.0)
    state.mobile_positive[0, 0] = 1.0
    state._pt_model = lambda: FixedPTRates(
        state.n_bins,
        peierls_s=1.0e12,
        release_s=0.0,
        jump_m=1.0e-6,
    )
    state.max_transport_substeps = 24
    initial = total_active_content(state)
    with installed_exponential_transport_v1005143():
        result = state.transport(
            dt_s=840.0, T_K=700.0, opening_stress_Pa=1.0e9
        )
    assert result["transport_integrator"] == TRANSPORT_INTEGRATOR
    assert result["transport_attempted_exponentials"] <= 3
    assert result["transport_rejected_intervals"] == 0
    assert result["max_frozen_courant"] > 1.0e12
    assert result["dN_escaped"] > 0.999999 * initial
    assert total_active_content(state) < 1.0e-8 * initial
    assert total_active_content(state) + result["dN_escaped"] == pytest.approx(
        initial, rel=2.0e-8, abs=1.0e-12
    )


def test_real_candidate_0118_rates_complete_840_second_transport():
    state = make_state(encounter_efficiency=9.160246308716648)
    state.mobile_positive[0, :4] = 0.25
    state.mobile_negative[1, :4] = 0.25
    state.max_transport_substeps = 240
    initial = total_active_content(state)
    with installed_exponential_transport_v1005143():
        result = state.transport(
            dt_s=840.0, T_K=700.0, opening_stress_Pa=2.0e9
        )
    assert result["transport_attempted_exponentials"] <= 240
    assert result["transport_integrator"] == TRANSPORT_INTEGRATOR
    assert np.min(state.mobile_positive) >= 0.0
    assert np.min(state.mobile_negative) >= 0.0
    assert np.min(state.retained_positive) >= 0.0
    assert np.min(state.retained_negative) >= 0.0
    assert total_active_content(state) + result["dN_escaped"] == pytest.approx(
        initial, rel=2.0e-8, abs=1.0e-11
    )


def test_zero_content_advances_clock_without_transport_work():
    state = make_state(encounter_efficiency=0.0)
    with installed_exponential_transport_v1005143():
        result = state.transport(
            dt_s=840.0, T_K=700.0, opening_stress_Pa=1.0e9
        )
    assert result["transport_substeps"] == 0
    assert result["transport_attempted_exponentials"] == 0
    assert state.time_s == pytest.approx(840.0)
