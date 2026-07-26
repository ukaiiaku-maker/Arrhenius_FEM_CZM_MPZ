from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from arrhenius_fracture.four_class_parameter_bridge_v100518 import (
    load_four_class_parameter_option,
)
from arrhenius_fracture.persistent_site_load_ramp_stochastic_emission_v1005183 import (
    LOAD_RAMP_SCHEMA,
    PersistentSiteLoadRampStochasticEmissionFrontEngineV1005183,
    two_channel_absolute_opening_drives_v1005183,
)
from arrhenius_fracture.sharp_front import (
    FrontConfig,
    default_cleavage_barrier,
    default_emission_barrier,
)
from arrhenius_fracture.signed_kernel_family_v1005141 import (
    load_signed_shielding_artifact_v1005141,
)

ROOT = Path(__file__).resolve().parents[1]
PARAMETER_ROOT = ROOT / "runtime_inputs" / "v10_0_5_18_four_class"
FAMILY = (
    ROOT
    / "runtime_inputs"
    / "v10_0_5_17_frozen_pf_inputs"
    / "signed_kernel"
    / "v10_2_14_active_only_campaign_family.json"
)
OPTION = "v913_paper_weakT01_0129902_persistent_sites"


def _engine(monkeypatch):
    monkeypatch.setenv("EMISSION_HAZARD_SEED", "1278783511")
    monkeypatch.setenv("EMISSION_RAMP_MAX_LOG_RATE_CHANGE", "0.5")
    monkeypatch.setenv("EMISSION_RAMP_ACTION_FLOOR", "1e-6")
    candidate, _ = load_four_class_parameter_option(PARAMETER_ROOT, OPTION)
    family = load_signed_shielding_artifact_v1005141(FAMILY)
    cls = PersistentSiteLoadRampStochasticEmissionFrontEngineV1005183
    cls.configure(candidate, family)
    cls.configure_stochastic(
        hazard_mode="exponential",
        hazard_seed=1278783511,
        hazard_minimum_threshold=1.0e-12,
        event_length_mode="threshold_scaled",
        event_minimum_factor=0.5,
        event_maximum_factor=4.0,
    )
    front = FrontConfig()
    front.r0 = 1.0e-6
    front.L_pz = 50.0e-6
    front.da = 5.0e-6
    front.sigma_cap = 30.0e9
    engine = cls(
        front,
        default_cleavage_barrier(),
        default_emission_barrier(2.74e-10),
        160.0e9,
        0.28,
        2.74e-10,
        SimpleNamespace(
            blunting_length_m=0.5e-6,
            max_transport_cfl=0.35,
            max_transport_substeps=2000,
        ),
    )
    engine._mm = SimpleNamespace(
        latest={
            "two_channel_drive_reliable": True,
            "two_channel_drive_factors": [1.0, 0.5],
            "two_channel_tau_signed_Pa": [1.0e8, -1.0e8],
            "two_channel_names": ["positive", "negative"],
            "cleavage_factor": 1.0,
            "emission_factor": 4.0,
        }
    )
    return engine


def test_two_channel_absolute_opening_drive_applies_anisotropy_once():
    engine = SimpleNamespace(
        persistent_site_source_active=True,
        _mm=SimpleNamespace(
            latest={
                "two_channel_drive_reliable": True,
                "cleavage_factor": 1.25,
                "emission_factor": 4.0,
            }
        ),
    )
    Kc, Ke, metadata = two_channel_absolute_opening_drives_v1005183(
        engine, 20.0e6
    )
    assert Kc == pytest.approx(25.0e6)
    assert Ke == pytest.approx(20.0e6)
    assert metadata["scalar_emission_factor_applied_to_Kemit"] is False


def test_actual_weakT_startup_ramp_predictor_is_finite_and_nonmutating(monkeypatch):
    engine = _engine(monkeypatch)
    before = engine._capture_state()
    predicted = engine.predict_clock_increment_drives(
        80.0e6,
        80.0e6,
        300.0,
        840.0,
    )
    assert np.isfinite(predicted)
    assert 0.0 <= predicted <= 1.0 + 1.0e-12
    after = engine._capture_state()
    assert after["mpz_state"].state_dict() == before["mpz_state"].state_dict()
    assert after["emission_event_count_total"] == before["emission_event_count_total"]
    assert after["hazard_rng_state"] == before["hazard_rng_state"]


def test_linear_ramp_is_less_aggressive_than_endpoint_hold(monkeypatch):
    ramp = _engine(monkeypatch)
    held = _engine(monkeypatch)
    endpoint = 30.0e6
    held.load_ramp_K_cleave_prev_Pa_sqrt_m = endpoint
    held.load_ramp_K_emit_prev_Pa_sqrt_m = endpoint
    predicted_ramp = ramp.predict_clock_increment_drives(
        endpoint, endpoint, 300.0, 840.0
    )
    predicted_held = held.predict_clock_increment_drives(
        endpoint, endpoint, 300.0, 840.0
    )
    assert predicted_ramp <= predicted_held + 1.0e-12


def test_predictor_never_calls_deterministic_expected_emission(monkeypatch):
    engine = _engine(monkeypatch)

    def forbidden(*args, **kwargs):
        raise AssertionError("deterministic emit_persistent predictor was invoked")

    monkeypatch.setattr(engine.mpz_state, "emit_persistent", forbidden)
    predicted = engine.predict_clock_increment_drives(
        5.0e6, 5.0e6, 300.0, 840.0
    )
    assert np.isfinite(predicted)


def test_probe_fallback_is_bounded_and_audited(monkeypatch):
    engine = _engine(monkeypatch)
    factors, tau, _ = engine._two_channel_drive()
    engine._mm.latest = {
        "two_channel_drive_reliable": False,
        "traction_probe_reliable": False,
    }
    factors2, tau2, metadata = engine._two_channel_drive()
    assert np.allclose(factors2, factors)
    assert np.allclose(tau2, tau)
    assert metadata["two_channel_probe_fallback_active"] is True
    assert engine.two_channel_probe_fallback_total == 1


def test_audit_reports_local_ramp_without_global_state_predictor(monkeypatch):
    engine = _engine(monkeypatch)
    audit = type(engine).audit_payload()["stochastic_emission"]
    assert audit["load_ramp_schema"] == LOAD_RAMP_SCHEMA
    assert audit["global_FEM_state_change_predictor"] is False
    assert audit["local_linear_K_ramp"] is True
    assert audit["post_onset_mean_field_burst"] is False
