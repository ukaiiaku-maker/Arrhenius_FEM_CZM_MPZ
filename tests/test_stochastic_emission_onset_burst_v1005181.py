from __future__ import annotations

import copy
import json
from types import MethodType, SimpleNamespace

import numpy as np

from arrhenius_fracture.persistent_site_registry_v100514 import (
    select_persistent_site_row,
)
from arrhenius_fracture.persistent_site_signed_mpz_v100514 import (
    SignedShieldingKernelV100514,
)
from arrhenius_fracture.persistent_site_stochastic_onset_burst_v1005181 import (
    PersistentSiteStochasticOnsetBurstMovingTipFrontEngineV1005181,
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
        source_path="synthetic_stochastic_onset_burst_kernel",
    )
    cls = PersistentSiteStochasticOnsetBurstMovingTipFrontEngineV1005181
    cls.configure(candidate, kernel)
    cls.configure_stochastic(
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
    engine = cls(
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


def install_constant_aggregate_hazard(engine, aggregate_hazard_s):
    multiplicity = engine.mpz_state.source_geometry()["multiplicity_per_system"]
    rate = float(aggregate_hazard_s) / float(multiplicity)
    engine._emission_rate_per_site = MethodType(
        lambda self, stress, temperature: rate if stress > 0.0 else 0.0,
        engine,
    )


def test_high_rate_episode_localizes_one_onset_then_resolves_burst(monkeypatch):
    engine = make_engine(monkeypatch)
    install_constant_aggregate_hazard(engine, 1.0e8)
    engine.emission_threshold_action[:] = [0.01, 1.0e30]
    engine.emission_action_current[:] = 0.0
    radius_before = engine.mpz_state.blunted_radius()

    result, _ = engine._plastic_half_step(
        dt_s=1.0,
        T_K=700.0,
        K_emit=20.0e6,
        drive_factors=np.array([0.2, 0.0]),
        tau_signed_Pa=np.array([1.0e8, -1.0e8]),
    )

    assert result["stochastic_emission_event_count"] == 1
    assert result["trigger_activations_by_system"][0] == 1.0
    assert result["post_onset_burst_activations_by_system"][0] >= 0.0
    assert result["source_activations"] >= 1.0
    assert result["stochastic_onset_exact"] is True
    assert result["post_onset_burst_mean_field"] is True
    assert result["individual_post_onset_activations_explicitly_sampled"] is False
    assert engine.mpz_state.blunted_radius() > radius_before
    assert engine.emission_event_count_total == 1


def test_prediction_is_nonmutating_under_stiff_emission(monkeypatch):
    engine = make_engine(monkeypatch)
    install_constant_aggregate_hazard(engine, 1.0e8)
    engine.emission_threshold_action[:] = [0.01, 1.0e30]
    before = engine._capture_state()
    before_mpz = engine.mpz_state.state_dict()

    predicted = engine.predict_clock_increment_drives(
        20.0e6,
        20.0e6,
        700.0,
        840.0,
    )

    assert np.isfinite(predicted)
    assert engine.emission_event_count_total == before["emission_event_count_total"]
    assert np.array_equal(
        engine.emission_threshold_action,
        before["emission_threshold_action"],
    )
    assert np.array_equal(
        engine.emission_action_current,
        before["emission_action_current"],
    )
    assert engine.mpz_state.state_dict() == before_mpz
    assert [
        json.dumps(rng.bit_generator.state, sort_keys=True)
        for rng in engine._emission_rngs
    ] == [
        json.dumps(state, sort_keys=True)
        for state in before["emission_rng_states"]
    ]


def test_trial_restore_recovers_onset_rng_burst_and_microstructure(monkeypatch):
    engine = make_engine(monkeypatch)
    install_constant_aggregate_hazard(engine, 1.0e6)
    engine.emission_threshold_action[:] = [0.01, 1.0e30]
    before = engine._capture_state()
    before_mpz = engine.mpz_state.state_dict()

    engine._plastic_half_step(
        dt_s=1.0,
        T_K=700.0,
        K_emit=20.0e6,
        drive_factors=np.array([0.2, 0.0]),
        tau_signed_Pa=np.array([1.0e8, -1.0e8]),
    )
    assert engine.emission_event_count_total == 1

    engine._restore_state(before)
    assert engine.emission_event_count_total == 0
    assert engine.mpz_state.state_dict() == before_mpz
    assert np.array_equal(
        engine.emission_threshold_action,
        before["emission_threshold_action"],
    )
    assert np.array_equal(
        engine.emission_action_current,
        before["emission_action_current"],
    )


def test_audit_states_exact_onset_and_mean_field_burst(monkeypatch):
    engine = make_engine(monkeypatch)
    audit = type(engine).audit_payload()["stochastic_emission"]
    assert audit["stochastic_onset_exact"] is True
    assert audit["post_onset_burst_mean_field"] is True
    assert audit["individual_post_onset_activations_explicitly_sampled"] is False
    assert audit["finite_source_inventory"] is False
    assert audit["source_refresh"] is False
