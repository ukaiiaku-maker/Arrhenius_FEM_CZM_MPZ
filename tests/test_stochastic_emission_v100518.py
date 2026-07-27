from __future__ import annotations

import copy
import json
from types import MethodType, SimpleNamespace

import numpy as np
import pytest

from arrhenius_fracture.persistent_site_registry_v100514 import (
    select_persistent_site_row,
)
from arrhenius_fracture.persistent_site_signed_mpz_v100514 import (
    SignedShieldingKernelV100514,
)
from arrhenius_fracture.persistent_site_stochastic_emission_v100518 import (
    PersistentSiteStochasticEmissionMovingTipFrontEngineV100518,
    derive_stream_seed,
)
from arrhenius_fracture.sharp_front import (
    FrontConfig,
    default_cleavage_barrier,
    default_emission_barrier,
)


def make_engine(monkeypatch, seed=421):
    monkeypatch.delenv("EMISSION_HAZARD_SEED", raising=False)
    monkeypatch.delenv("EMISSION_HAZARD_MIN_THRESHOLD", raising=False)
    monkeypatch.delenv("EMISSION_MAX_ACTION_SUBSTEP", raising=False)
    monkeypatch.delenv("EMISSION_MAX_EVENTS_PER_HALF_STEP", raising=False)
    candidate = select_persistent_site_row("v912_peak_0118_persistent_sites")
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
        source_path="synthetic_stochastic_emission_kernel",
    )
    PersistentSiteStochasticEmissionMovingTipFrontEngineV100518.configure(
        candidate, kernel
    )
    PersistentSiteStochasticEmissionMovingTipFrontEngineV100518.configure_stochastic(
        hazard_mode="exponential",
        hazard_seed=seed,
        hazard_minimum_threshold=1.0e-12,
        event_length_mode="threshold_scaled",
        event_minimum_factor=0.5,
        event_maximum_factor=4.0,
    )
    front = FrontConfig()
    front.r0 = 1.0e-6
    front.L_pz = 50.0e-6
    front.da = 5.0e-6
    front.sigma_cap = 0.0
    mpz_cfg = SimpleNamespace(
        blunting_length_m=0.5e-6,
        max_transport_cfl=0.35,
        max_transport_substeps=2000,
    )
    engine = PersistentSiteStochasticEmissionMovingTipFrontEngineV100518(
        front,
        default_cleavage_barrier(),
        default_emission_barrier(2.74e-10),
        160.0e9,
        0.28,
        2.74e-10,
        mpz_cfg,
    )
    engine._mm = SimpleNamespace(
        latest={
            "two_channel_drive_reliable": True,
            "two_channel_drive_factors": [0.2, 0.1],
            "two_channel_tau_signed_Pa": [1.0e8, -1.0e8],
            "two_channel_names": ["positive", "negative"],
        }
    )
    engine.max_action_substep = 1.0
    engine.max_translation_substep_m = 1.0
    return engine


def install_constant_emission_rate(engine, aggregate_hazard_per_system_s=1.0):
    multiplicity = engine.mpz_state.source_geometry()["multiplicity_per_system"]
    rate = float(aggregate_hazard_per_system_s) / float(multiplicity)
    engine._emission_rate_per_site = MethodType(
        lambda self, stress, temperature: rate if stress > 0.0 else 0.0,
        engine,
    )
    return rate


def install_no_emission(engine):
    engine._emission_rate_per_site = MethodType(
        lambda self, stress, temperature: 0.0,
        engine,
    )


def test_domain_separated_emission_seeds_are_stable_and_distinct():
    first = derive_stream_seed(1234, "signed_emission", 0)
    second = derive_stream_seed(1234, "signed_emission", 1)
    assert first == derive_stream_seed(1234, "signed_emission", 0)
    assert first != second
    assert first != derive_stream_seed(1235, "signed_emission", 0)


def test_zero_emission_hazard_produces_no_event(monkeypatch):
    engine = make_engine(monkeypatch)
    install_no_emission(engine)
    before = engine.mpz_state.emitted_total
    result, _ = engine._plastic_half_step(
        dt_s=1.0,
        T_K=700.0,
        K_emit=20.0e6,
        drive_factors=np.array([0.2, 0.1]),
        tau_signed_Pa=np.array([1.0e8, -1.0e8]),
    )
    assert result["stochastic_emission_event_count"] == 0
    assert result["dN_emit"] == 0.0
    assert engine.mpz_state.emitted_total == before


def test_one_stochastic_event_adds_one_packet_and_blunts_tip(monkeypatch):
    engine = make_engine(monkeypatch)
    install_constant_emission_rate(engine, aggregate_hazard_per_system_s=1.0)
    engine.emission_threshold_action[:] = [0.05, 1.0e9]
    engine.emission_action_current[:] = 0.0
    engine._emission_rngs[0] = np.random.default_rng(12345)
    radius_before = engine.mpz_state.blunted_radius()
    result, _ = engine._plastic_half_step(
        dt_s=0.051,
        T_K=700.0,
        K_emit=20.0e6,
        drive_factors=np.array([0.2, 0.0]),
        tau_signed_Pa=np.array([1.0e8, -1.0e8]),
    )
    assert result["stochastic_emission_event_count"] == 1
    assert result["source_activations"] == 1.0
    assert result["dN_emit"] == pytest.approx(1.0)
    assert engine.mpz_state.emitted_total == pytest.approx(1.0)
    assert engine.mpz_state.blunted_radius() > radius_before
    event = result["stochastic_emission_events"][0]
    assert event["system"] == 0
    assert event["sign"] == 1.0
    assert event["line_content_added"] == pytest.approx(1.0)


def test_trial_restore_recovers_emission_rng_clock_and_microstructure(monkeypatch):
    engine = make_engine(monkeypatch)
    install_constant_emission_rate(engine, aggregate_hazard_per_system_s=1.0)
    engine.emission_threshold_action[:] = [0.05, 1.0e9]
    engine.emission_action_current[:] = 0.0
    before = engine._capture_state()
    before_state = engine.mpz_state.state_dict()
    engine._plastic_half_step(
        dt_s=0.051,
        T_K=700.0,
        K_emit=20.0e6,
        drive_factors=np.array([0.2, 0.0]),
        tau_signed_Pa=np.array([1.0e8, -1.0e8]),
    )
    assert engine.emission_event_count_total == 1
    engine._restore_state(before)
    assert engine.emission_event_count_total == 0
    assert np.array_equal(
        engine.emission_threshold_action,
        before["emission_threshold_action"],
    )
    assert np.array_equal(
        engine.emission_action_current,
        before["emission_action_current"],
    )
    assert engine.mpz_state.state_dict() == before_state
    restored_rng = json.dumps(
        engine._emission_rngs[0].bit_generator.state,
        sort_keys=True,
    )
    expected_rng = json.dumps(before["emission_rng_states"][0], sort_keys=True)
    assert restored_rng == expected_rng


def test_stochastic_crack_progress_advects_mpz_before_checkpoint(monkeypatch):
    engine = make_engine(monkeypatch)
    install_no_emission(engine)
    engine.lambda_cleave = MethodType(
        lambda self, sigma, temperature: (0.25, 0.25, 1.0e-19),
        engine,
    )
    threshold = engine.hazard_threshold_action
    expected_progress = 0.25 / threshold
    result = engine._integrate_coupled(
        K_cleave=20.0e6,
        K_emit=20.0e6,
        T_K=700.0,
        dt_s=1.0,
        drive_factors=np.array([0.2, 0.1]),
        tau_signed_Pa=np.array([1.0e8, -1.0e8]),
    )
    expected_da = min(expected_progress, 1.0) * engine.stochastic_event_advance_m
    assert result["da"] == pytest.approx(expected_da)
    assert engine.mpz_state.advance_total_m == pytest.approx(expected_da)
    if expected_progress < 1.0:
        assert result["fired"] is False
        assert engine.checkpoint_advance_total_m == 0.0


def test_nonmutating_prediction_preserves_both_rng_streams(monkeypatch):
    engine = make_engine(monkeypatch)
    install_no_emission(engine)
    before_hazard = copy.deepcopy(engine._hazard_rng.bit_generator.state)
    before_emission = [
        copy.deepcopy(rng.bit_generator.state) for rng in engine._emission_rngs
    ]
    before_mpz = engine.mpz_state.state_dict()
    engine.predict_clock_increment_drives(20.0e6, 20.0e6, 700.0, 1.0)
    assert engine._hazard_rng.bit_generator.state == before_hazard
    assert [rng.bit_generator.state for rng in engine._emission_rngs] == before_emission
    assert engine.mpz_state.state_dict() == before_mpz
