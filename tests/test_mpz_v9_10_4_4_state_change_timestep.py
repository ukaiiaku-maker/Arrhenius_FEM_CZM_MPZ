import math

import numpy as np

from arrhenius_fracture.reduced_campaign_front_v9104 import (
    ReducedCampaignFront,
    ReducedFrontSettings,
    simulate_reduced_response,
)
from arrhenius_fracture.reduced_campaign_front_v91044_timestep_guard import (
    _limiter_diagnostics,
    _time_for_exact_change,
)


def parameters():
    return {
        "cleave_G00_eV": 2.0, "cleave_gT_eV_per_K": 0.0,
        "cleave_sigc0_GPa": 4.0, "cleave_sT_GPa_per_K": 0.0,
        "cleave_exp_a": 1.0, "cleave_exp_n": 1.0, "cleave_floor_frac": 0.05,
        "emit_G00_eV": 2.0, "emit_gT_eV_per_K": 0.0,
        "emit_sigc0_GPa": 4.0, "emit_sT_GPa_per_K": 0.0,
        "emit_exp_a": 1.0, "emit_exp_n": 1.0, "emit_floor_frac": 0.05,
        "peierls_H0_eV": 1.0, "peierls_exp_a": 1.0, "peierls_exp_n": 1.0,
        "delta_H_PT_eV": 1.0, "taylor_H0_eV": 2.0,
        "taylor_exp_a": 1.0, "taylor_exp_n": 1.0,
        "peierls_activation_entropy_kB": 0.0,
        "taylor_activation_entropy_kB": 0.0,
        "taylor_corr_rho_c_m2": 1.0e14, "taylor_corr_scale": 1.0,
        "source_sites_per_system": 10.0, "encounter_efficiency": 1.0,
        "retained_recovery_rate_s": 0.0, "source_refresh_length_um": 50.0,
        "peierls_nu0_s": 1.0e12, "taylor_nu0_s": 1.0e11, "c_blunt": 1.0,
    }


def rates(**overrides):
    base = {
        "lambda_c_s": 0.0,
        "lambda_e_s": np.zeros(2),
        "encounter_rate_s": np.zeros(2),
        "taylor_rate_s": np.zeros(2),
        "velocity_m_s": np.zeros(2),
    }
    base.update(overrides)
    return base


def test_exact_change_has_no_limit_when_asymptotic_change_is_small():
    assert math.isinf(_time_for_exact_change([0.5], [1.0e30], [1.0]))


def test_exact_change_limit_reproduces_requested_change():
    h = _time_for_exact_change([10.0], [4.0], [2.0])
    changed = 10.0 * (1.0 - math.exp(-4.0 * h))
    assert np.isclose(changed, 2.0)


def test_huge_exchange_rate_at_equilibrium_does_not_limit_loading_step():
    front = ReducedCampaignFront(parameters(), ReducedFrontSettings(), mode="full")
    front.available[:] = 0.0
    front.mobile[:] = 5.0
    front.retained[:] = 5.0
    diag = _limiter_diagnostics(
        front,
        rates(
            encounter_rate_s=np.full(2, 1.0e30),
            taylor_rate_s=np.full(2, 1.0e30),
        ),
    )
    expected = front.s.max_dK_substep_MPa_sqrt_m / front.s.Kdot_MPa_sqrt_m_s
    assert diag["dt_limiter"] == "loading"
    assert diag["dt_selected_s"] == expected
    assert math.isinf(diag["dt_limit_exchange_s"])


def test_escape_does_not_limit_when_all_content_is_retained():
    front = ReducedCampaignFront(parameters(), ReducedFrontSettings(), mode="full")
    front.available[:] = 0.0
    front.mobile[:] = 0.0
    front.retained[:] = 10.0
    diag = _limiter_diagnostics(
        front,
        rates(velocity_m_s=np.full(2, 1.0e20)),
    )
    assert math.isinf(diag["dt_limit_escape_s"])


def test_incomplete_run_summary_contains_terminal_limiter_diagnostics():
    settings = ReducedFrontSettings(Kmax_MPa_sqrt_m=0.01, max_internal_steps=2)
    result = simulate_reduced_response(parameters(), 300.0, settings, mode="plasticity_off")
    assert "termination_reason" in result
    assert "terminal_K_MPa_sqrt_m" in result
    assert "terminal_B" in result
    assert "dominant_dt_limiter" in result
    assert "dt_limiter_counts" in result
