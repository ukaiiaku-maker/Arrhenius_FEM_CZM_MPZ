from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest
from scipy.linalg import expm

from arrhenius_fracture.persistent_site_registry_v100514 import (
    select_persistent_site_row,
)
from arrhenius_fracture.persistent_site_signed_mpz_v100514 import (
    PersistentSiteSignedMPZStateV100514,
    SignedShieldingKernelV100514,
)
from arrhenius_fracture.persistent_site_transport_v1005145 import (
    ASYMPTOTIC_MODEL,
    TRANSPORT_INTEGRATOR,
    _mobile_absorption_probabilities,
    _retained_generator,
    installed_asymptotic_transport_v1005145,
)


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
        source_path="synthetic_transport_kernel_v1005145",
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


def test_mobile_absorption_probabilities_are_substochastic_and_directional():
    outflow = np.array([4.0, 3.0, 2.0, 1.0])
    encounter = np.array([1.0, 2.0, 3.0, 4.0])
    P, escape = _mobile_absorption_probabilities(outflow, encounter)
    assert np.min(P) >= 0.0
    assert np.min(escape) >= 0.0
    assert np.allclose(np.sum(P, axis=1) + escape, 1.0, rtol=0.0, atol=2.0e-14)
    assert np.allclose(np.tril(P, -1), 0.0)


def test_asymptotic_reduction_matches_exact_fast_slow_generator():
    n = 5
    fast_scale = 1.0e7
    outflow = fast_scale * np.array([3.0, 2.5, 2.0, 1.5, 1.0])
    encounter = fast_scale * np.array([1.0, 1.5, 2.0, 2.5, 3.0])
    release = np.array([0.2, 0.4, 0.7, 1.0, 1.3])
    P, escape = _mobile_absorption_probabilities(outflow, encounter)
    Q = _retained_generator(P, escape, release)

    A = np.zeros((2 * n, 2 * n))
    for i in range(n):
        A[i, i] -= outflow[i] + encounter[i]
        if i + 1 < n:
            A[i + 1, i] += outflow[i]
        A[n + i, i] += encounter[i]
        A[i, n + i] += release[i]
        A[n + i, n + i] -= release[i]

    mobile0 = np.zeros(n)
    mobile0[0] = 0.8
    mobile0[2] = 0.2
    retained0 = np.array([0.1, 0.2, 0.3, 0.1, 0.0])
    y0 = np.concatenate([mobile0, retained0])
    dt = 2.0
    exact = expm(dt * A) @ y0

    retained_start = retained0 + mobile0 @ P
    reduced_retained = expm(dt * Q) @ retained_start
    reduced = np.concatenate([np.zeros(n), reduced_retained])
    assert np.sum(np.abs(exact - reduced)) / np.sum(y0) < 2.0e-6
    exact_escape = np.sum(y0) - np.sum(exact)
    reduced_escape = np.sum(y0) - np.sum(reduced)
    assert reduced_escape == pytest.approx(exact_escape, rel=2.0e-6, abs=1.0e-10)


@pytest.mark.parametrize("temperature_K", [300.0, 400.0, 500.0, 600.0, 700.0, 800.0])
def test_candidate_0118_low_temperature_production_scale_finishes(
    temperature_K,
):
    state = make_state()
    # Production logs show active retained populations ranging from fractions of
    # a line to hundreds of lines.  Exercise both signs and both mobile/retained
    # populations at the corresponding forest-density scale.
    state.mobile_positive[0, 0] = 180.0
    state.mobile_negative[1, 0] = 35.0
    state.retained_positive[0, 1] = 80.0
    state.retained_negative[1, 2] = 20.0
    state.max_transport_substeps = 2000
    initial = total_active_content(state)
    with installed_asymptotic_transport_v1005145():
        result = state.transport(
            dt_s=840.0,
            T_K=temperature_K,
            opening_stress_Pa=8.0e9,
        )
    assert result["transport_integrator"] == TRANSPORT_INTEGRATOR
    assert result["transport_asymptotic_active"] is True
    assert result["physical_generator_action"] == ASYMPTOTIC_MODEL
    assert result["transport_asymptotic_iterations"] <= 80
    assert np.isfinite(result["transport_nonlinear_error_final"])
    assert result["transport_nonlinear_error_final"] <= result["transport_nonlinear_rtol"]
    assert np.min(state.mobile_positive) >= 0.0
    assert np.min(state.mobile_negative) >= 0.0
    assert np.min(state.retained_positive) >= 0.0
    assert np.min(state.retained_negative) >= 0.0
    final = total_active_content(state)
    assert final + result["dN_escaped"] == pytest.approx(
        initial, rel=2.0e-9, abs=1.0e-10
    )


def test_high_temperature_path_remains_available_and_conservative():
    state = make_state()
    state.mobile_positive[0, :4] = 5.0e-7
    state.mobile_negative[1, :4] = 7.5e-8
    state.retained_positive[0, 4:8] = 2.5e-8
    state.max_transport_substeps = 2000
    initial = total_active_content(state)
    with installed_asymptotic_transport_v1005145():
        result = state.transport(
            dt_s=840.0,
            T_K=1100.0,
            opening_stress_Pa=4.0e9,
        )
    assert result["transport_integrator"] == TRANSPORT_INTEGRATOR
    final = total_active_content(state)
    assert final + result["dN_escaped"] == pytest.approx(
        initial, rel=2.0e-8, abs=1.0e-17
    )
