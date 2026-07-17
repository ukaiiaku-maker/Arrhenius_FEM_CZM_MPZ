import math

import numpy as np

from arrhenius_fracture.reduced_campaign_front_v9104 import (
    ReducedCampaignFront,
    ReducedFrontSettings,
    exact_depletion,
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


def rates(lambda_e=0.0, kinetic=0.0):
    return {
        "lambda_c_s": 0.0,
        "lambda_e_s": np.full(2, lambda_e),
        "encounter_rate_s": np.full(2, kinetic),
        "taylor_rate_s": np.zeros(2),
        "velocity_m_s": np.zeros(2),
    }


def test_trace_available_source_does_not_limit_global_step():
    front = ReducedCampaignFront(parameters(), ReducedFrontSettings(), mode="full")
    front.available[:] = 1.0e-6
    front.mobile[:] = 0.0
    front.retained[:] = 0.0
    expected = front.s.max_dK_substep_MPa_sqrt_m / front.s.Kdot_MPa_sqrt_m_s
    assert front._choose_dt(rates(lambda_e=1.0e30)) == expected


def test_emission_limit_is_based_on_absolute_fraction_of_capacity():
    front = ReducedCampaignFront(parameters(), ReducedFrontSettings(), mode="full")
    front.available[:] = front.capacity
    rate = 100.0
    h = front._choose_dt(rates(lambda_e=rate))
    allowed = front.s.max_emit_fraction_substep * front.capacity
    expected = -math.log1p(-allowed[0] / front.available[0]) / rate
    assert math.isclose(h, expected, rel_tol=1.0e-12)
    emitted = exact_depletion(front.available, np.full(2, rate), h)
    assert np.all(emitted <= allowed * (1.0 + 1.0e-12))


def test_trace_active_population_does_not_limit_global_step():
    front = ReducedCampaignFront(parameters(), ReducedFrontSettings(), mode="full")
    front.available[:] = 0.0
    front.mobile[:] = 1.0e-6
    front.retained[:] = 0.0
    expected = front.s.max_dK_substep_MPa_sqrt_m / front.s.Kdot_MPa_sqrt_m_s
    assert front._choose_dt(rates(kinetic=1.0e30)) == expected
