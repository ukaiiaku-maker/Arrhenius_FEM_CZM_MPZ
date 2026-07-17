from __future__ import annotations

import copy
from types import SimpleNamespace

import numpy as np
import pytest

from arrhenius_fracture.kinetic_campaign_czm import (
    CampaignCalibratedCZMFrontEngine,
    CampaignKineticMPZState,
    KineticCampaignCZMConfig,
    apply_pf_manifest_to_mpz_config,
)
from arrhenius_fracture.moving_process_zone import MovingProcessZoneConfig
from arrhenius_fracture.pf_equivalent_material_manifest import load_material_manifest


def make_state(material="weakT", n_bins=100):
    manifest = load_material_manifest(material)
    cfg = MovingProcessZoneConfig(length_m=100e-6, n_bins=n_bins, n_systems=2)
    kinetic = KineticCampaignCZMConfig()
    apply_pf_manifest_to_mpz_config(cfg, manifest, kinetic)
    state = CampaignKineticMPZState(
        cfg,
        manifest,
        b=2.48e-10,
        G_Pa=80e9,
        kinetic_cfg=kinetic,
    )
    return state, manifest, kinetic


def make_channel_engine(state, manifest, kinetic):
    eng = object.__new__(CampaignCalibratedCZMFrontEngine)
    eng.mpz_state = state
    eng.manifest = manifest
    eng.kinetic_config = kinetic
    eng.G = 80e9
    eng.nu = 0.28
    eng.b = 2.48e-10
    eng.f = SimpleNamespace(r0=1e-6, c_blunt=manifest.c_blunt)
    eng._last_channels = {}
    return eng


def test_finite_source_depletion_is_timestep_partition_invariant():
    state1, _, _ = make_state()
    state2 = copy.deepcopy(state1)
    sigma = 6e9
    T = 700.0
    total = 1e-3
    state1.emit_exact(total, sigma, T, np.ones(2))
    state2.emit_exact(0.5 * total, sigma, T, np.ones(2))
    state2.emit_exact(0.5 * total, sigma, T, np.ones(2))
    assert np.allclose(state1.available_sites, state2.available_sites, rtol=1e-12, atol=1e-12)
    assert np.all(state1.available_sites <= state1.site_capacity)
    assert np.all(state1.available_sites >= 0.0)


def test_source_refresh_is_advance_only_and_partition_invariant():
    state1, _, _ = make_state()
    state2 = copy.deepcopy(state1)
    state1.available_sites *= 0.2
    state2.available_sites *= 0.2
    before = state1.available_sites.copy()
    state1.advance_campaign(0.0)
    assert np.array_equal(state1.available_sites, before)

    state1.advance_campaign(5e-6)
    state2.advance_campaign(2e-6)
    state2.advance_campaign(3e-6)
    assert np.allclose(state1.available_sites, state2.available_sites, rtol=1e-12, atol=1e-12)
    assert np.all(state1.available_sites <= state1.site_capacity)


def test_shielding_changes_cleavage_but_not_opening():
    state, manifest, kinetic = make_state()
    eng = make_channel_engine(state, manifest, kinetic)
    K = 20e6
    opening0 = eng.sigma_opening_tip(K)
    cleavage0 = eng.sigma_cleavage_tip(K)
    state.retained[:, :4] = 50.0
    opening1 = eng.sigma_opening_tip(K)
    cleavage1 = eng.sigma_cleavage_tip(K)
    assert opening1 == pytest.approx(opening0)
    assert cleavage1 < cleavage0


def test_taylor_backstress_changes_emission_but_not_cleavage():
    state, manifest, kinetic = make_state()
    eng = make_channel_engine(state, manifest, kinetic)
    K = 20e6
    base = eng.stress_channels(K, K, np.ones(2))
    state.mobile[:, :4] = 100.0
    changed = eng.stress_channels(K, K, np.ones(2))
    assert changed["sigma_emission_effective_Pa"] < base["sigma_emission_effective_Pa"]
    back = state.taylor_backstress_Pa().copy()
    state._kinetic_cfg.backstress_scale = 0.0
    no_back = eng.stress_channels(K, K, np.ones(2))
    assert no_back["sigma_cleave_eff_Pa"] == pytest.approx(changed["sigma_cleave_eff_Pa"])
    assert np.mean(back) > 0.0


def test_shielding_is_applied_once():
    state, manifest, kinetic = make_state()
    eng = make_channel_engine(state, manifest, kinetic)
    state.retained[:, :5] = 10.0
    K = 25e6
    Ksh = eng.active_K_shielding()
    expected = max(K - Ksh, 0.0) / np.sqrt(2.0 * np.pi * eng.r_eff())
    assert eng.sigma_cleavage_tip(K) == pytest.approx(expected)


def test_local_backstress_averaging_length_does_not_expand_with_blunting():
    state, _, _ = make_state()
    state.mobile[:, :5] = 10.0
    rho0 = state.local_backstress_density_m2()
    state.accumulated_slip[:, :5] = 1e8
    rho1 = state.local_backstress_density_m2()
    assert np.array_equal(rho0, rho1)
