import math

import numpy as np

from arrhenius_fracture.reduced_campaign_front_v9104 import (
    ReducedCampaignFront,
    ReducedFrontSettings,
    best_adjacent_transition,
    exact_depletion,
    exact_refresh,
)


def parameters():
    return {
        "cleave_G00_eV": 2.0,
        "cleave_gT_eV_per_K": 0.0,
        "cleave_sigc0_GPa": 4.0,
        "cleave_sT_GPa_per_K": 0.0,
        "cleave_exp_a": 1.0,
        "cleave_exp_n": 1.0,
        "cleave_floor_frac": 0.05,
        "emit_G00_eV": 2.0,
        "emit_gT_eV_per_K": 0.0,
        "emit_sigc0_GPa": 4.0,
        "emit_sT_GPa_per_K": 0.0,
        "emit_exp_a": 1.0,
        "emit_exp_n": 1.0,
        "emit_floor_frac": 0.05,
        "peierls_H0_eV": 1.0,
        "peierls_exp_a": 1.0,
        "peierls_exp_n": 1.0,
        "delta_H_PT_eV": 1.0,
        "taylor_H0_eV": 2.0,
        "taylor_exp_a": 1.0,
        "taylor_exp_n": 1.0,
        "peierls_activation_entropy_kB": 0.0,
        "taylor_activation_entropy_kB": 0.0,
        "taylor_corr_rho_c_m2": 1.0e14,
        "taylor_corr_scale": 1.0,
        "source_sites_per_system": 10.0,
        "encounter_efficiency": 1.0,
        "retained_recovery_rate_s": 0.0,
        "source_refresh_length_um": 50.0,
        "peierls_nu0_s": 1.0e12,
        "taylor_nu0_s": 1.0e11,
        "c_blunt": 1.0,
    }


def test_exact_source_depletion_is_partition_invariant():
    available = np.asarray([10.0, 20.0])
    rate = np.asarray([0.3, 1.7])
    one = available - exact_depletion(available, rate, 2.0)
    half = available - exact_depletion(available, rate, 1.0)
    two_halves = half - exact_depletion(half, rate, 1.0)
    np.testing.assert_allclose(one, two_halves, rtol=1e-13, atol=1e-13)


def test_source_refresh_requires_crack_advance_and_is_partition_invariant():
    available = np.asarray([1.0, 4.0])
    capacity = np.asarray([10.0, 10.0])
    np.testing.assert_allclose(exact_refresh(available, capacity, 0.0, 50e-6), available)
    one = exact_refresh(available, capacity, 10e-6, 50e-6)
    half = exact_refresh(available, capacity, 5e-6, 50e-6)
    two_halves = exact_refresh(half, capacity, 5e-6, 50e-6)
    np.testing.assert_allclose(one, two_halves, rtol=1e-13, atol=1e-13)


def test_shielding_changes_cleavage_but_not_opening_stress():
    settings = ReducedFrontSettings(max_K_shield_MPa_sqrt_m=20.0)
    full = ReducedCampaignFront(parameters(), settings, mode="full")
    no_shield = ReducedCampaignFront(parameters(), settings, mode="shielding_off")
    full.retained[:] = 50.0
    no_shield.retained[:] = 50.0
    a = full.stress_channels(20.0)
    b = no_shield.stress_channels(20.0)
    assert a["K_shield_MPa_sqrt_m"] > 0.0
    assert math.isclose(a["sigma_open_Pa"], b["sigma_open_Pa"], rel_tol=1e-14)
    assert a["sigma_cleave_Pa"] < b["sigma_cleave_Pa"]


def test_backstress_changes_emission_but_not_cleavage_stress():
    settings = ReducedFrontSettings(max_K_shield_MPa_sqrt_m=20.0)
    full = ReducedCampaignFront(parameters(), settings, mode="full")
    no_back = ReducedCampaignFront(parameters(), settings, mode="backstress_off")
    for front in (full, no_back):
        front.mobile[:] = 100.0
        front.retained[:] = 20.0
    a = full.stress_channels(20.0)
    b = no_back.stress_channels(20.0)
    assert np.max(a["sigma_back_Pa"]) > 0.0
    assert math.isclose(a["sigma_cleave_Pa"], b["sigma_cleave_Pa"], rel_tol=1e-14)
    assert np.max(a["sigma_emit_Pa"]) < np.max(b["sigma_emit_Pa"])


def test_free_transition_finder_selects_one_dominant_adjacent_jump():
    temperatures = np.arange(300.0, 1200.0, 100.0)
    full = np.asarray([15.0, 15.2, 14.9, 15.1, 15.3, 31.0, 31.5, 31.2, 31.7])
    off = np.asarray([15.0, 14.8, 14.5, 14.1, 13.8, 13.5, 13.2, 13.0, 12.8])
    result = best_adjacent_transition(
        temperatures,
        full,
        plasticity_off_toughness=off,
    )
    assert result["transition_low_K"] == 700.0
    assert result["transition_high_K"] == 800.0
    assert result["shelf_ratio"] >= 2.0
    assert result["jump_concentration"] >= 0.9
    assert result["plasticity_off_ratio"] < 1.0
