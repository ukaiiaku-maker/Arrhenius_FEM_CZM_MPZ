from __future__ import annotations

import faulthandler
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

import arrhenius_fracture.persistent_site_joint_K_ramp_emission_v10051832 as joint32
from arrhenius_fracture.four_class_parameter_bridge_v100518 import (
    load_four_class_parameter_option,
)
from arrhenius_fracture.persistent_site_joint_K_ramp_emission_v10051832 import (
    EVENT_DRIVEN_SCHEMA,
    LOAD_RAMP_SCHEMA,
    OUTER_EMISSION_LIMITER_SCHEMA,
    PersistentSiteJointKRampEventDrivenFrontEngineV10051832,
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
PEAK = "v913_paper_peak01_0242980_persistent_sites"
WEAKT = "v913_paper_weakT01_0129902_persistent_sites"


def _engine(monkeypatch, *, option=WEAKT, seed=331085649, factors=(1.0, 0.5)):
    monkeypatch.setenv("EMISSION_HAZARD_SEED", str(seed))
    monkeypatch.setenv("EMISSION_RAMP_MAX_LOG_RATE_CHANGE", "0.25")
    monkeypatch.setenv("EMISSION_RAMP_ACTION_FLOOR", "1e-6")
    monkeypatch.setenv("EMISSION_INNER_MAX_LOG_HAZARD_CHANGE", "0.25")
    monkeypatch.setenv("EMISSION_INNER_ACTION_FLOOR", "1e-6")
    monkeypatch.setenv("EMISSION_EVENT_HORIZON_FACTOR", "1.25")
    candidate, _ = load_four_class_parameter_option(PARAMETER_ROOT, option)
    family = load_signed_shielding_artifact_v1005141(FAMILY)
    cls = PersistentSiteJointKRampEventDrivenFrontEngineV10051832
    cls.configure(candidate, family)
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
            "two_channel_drive_factors": list(factors),
            "two_channel_tau_signed_Pa": [1.0e8, -1.0e8],
            "two_channel_names": ["positive", "negative"],
            "cleavage_factor": 1.0,
            "emission_factor": 4.0,
        }
    )
    return engine


def test_outer_cell_has_no_emission_action_or_rate_limiter(monkeypatch):
    engine = _engine(monkeypatch)
    engine.load_ramp_K_cleave_prev_Pa_sqrt_m = 5.0e6
    engine.load_ramp_K_emit_prev_Pa_sqrt_m = 0.0
    monkeypatch.setattr(engine, "lambda_cleave", lambda sigma, T: (0.0, 0.0, 0.0))

    h, diagnostics = engine._local_ramp_substep(
        remaining_s=10.0,
        consumed_s=0.0,
        total_s=10.0,
        K_cleave_end=5.0e6,
        K_emit_end=100.0e6,
        T_K=1000.0,
        drive_factors=np.asarray([1.0, 0.5]),
        tau_signed_Pa=np.asarray([1.0e8, -1.0e8]),
    )

    assert h == pytest.approx(10.0)
    assert diagnostics["emission_action"] == 0.0
    assert diagnostics["emission_log_change"] == 0.0
    assert diagnostics["outer_emission_action_limiter_active"] == 0.0
    assert diagnostics["outer_emission_rate_variation_limiter_active"] == 0.0


def test_zero_initial_hazard_is_localized_inside_linear_K_ramp(monkeypatch):
    engine = _engine(monkeypatch)
    engine.emission_threshold_action = np.asarray([0.25, 100.0], dtype=float)
    engine.emission_action_current = np.zeros(2, dtype=float)
    monkeypatch.setattr(engine, "_opening_stress", lambda K: (float(K), 0.0, 0.0, 0.0))

    def fake_hazards(state, *, opening_stress_Pa, **kwargs):
        geometry = state.source_geometry()
        multiplicity = float(geometry["multiplicity_per_system"])
        rate = max(float(opening_stress_Pa), 0.0)
        return {
            "geometry": geometry,
            "multiplicity": multiplicity,
            "rho_back_m2": np.zeros(2, dtype=float),
            "sigma_back_Pa": np.zeros(2, dtype=float),
            "rate_per_site_s": np.asarray([rate / multiplicity, 0.0]),
            "aggregate_hazard_s": np.asarray([rate, 0.0]),
            "signs": np.asarray([1.0, -1.0]),
        }

    def fake_transport(state, *, dt_s, **kwargs):
        state.time_s += float(dt_s)
        return {
            "dN_emit": 0.0,
            "dN_trapped": 0.0,
            "dN_released": 0.0,
            "dN_recovered": 0.0,
            "dN_escaped": 0.0,
            "transport_substeps": 1,
        }

    monkeypatch.setattr(joint32, "_emission_hazards", fake_hazards)
    monkeypatch.setattr(joint32, "advance_pf_transport_only_v10222", fake_transport)
    monkeypatch.setattr(joint32, "_draw_exponential", lambda rng, floor: 100.0)

    result = engine._stochastic_emission_transport_ramp(
        dt_s=1.0,
        T_K=1000.0,
        K_emit_start=0.0,
        K_emit_end=2.0,
        drive_factors=np.asarray([1.0, 0.5]),
        tau_signed_Pa=np.asarray([1.0e8, -1.0e8]),
    )

    assert result["stochastic_emission_event_count"] == 1
    event = result["stochastic_emission_events"][0]
    assert event["time_within_half_step_s"] == pytest.approx(0.5)
    assert event["K_emit_at_event_Pa_sqrt_m"] == pytest.approx(1.0)
    assert event["opening_stress_at_event_Pa"] == pytest.approx(1.0)
    assert result["joint_K_ramp_inside_event_horizon"] is True


def test_peak_1000K_physical_startup_predictor_is_finite_and_nonmutating(monkeypatch):
    engine = _engine(
        monkeypatch,
        option=PEAK,
        seed=3100913530,
        factors=(1.0, 0.5),
    )
    before = engine._capture_state()
    faulthandler.dump_traceback_later(30.0, repeat=True)
    try:
        predicted = engine.predict_clock_increment_drives(
            16.5e6,
            16.5e6,
            1000.0,
            840.0,
        )
    finally:
        faulthandler.cancel_dump_traceback_later()
    assert np.isfinite(predicted)
    assert 0.0 <= predicted <= 1.0 + 1.0e-12
    after = engine._capture_state()
    assert after["mpz_state"].state_dict() == before["mpz_state"].state_dict()
    assert after["emission_event_count_total"] == before["emission_event_count_total"]
    assert after["hazard_rng_state"] == before["hazard_rng_state"]


def test_audit_records_joint_K_event_horizon(monkeypatch):
    engine = _engine(monkeypatch)
    audit = type(engine).audit_payload()["stochastic_emission"]
    assert audit["load_ramp_schema"] == LOAD_RAMP_SCHEMA
    assert audit["event_driven_transport_schema"] == EVENT_DRIVEN_SCHEMA
    assert audit["outer_emission_action_limiter_active"] is False
    assert audit["outer_emission_rate_variation_limiter_active"] is False
    assert audit["K_ramp_evaluated_inside_event_horizon"] is True
    assert audit["hazard_drive_quadrature"] == "linear_K_endpoint_trapezoid"
    assert audit["events_batched"] is False
    assert audit["accepted_update"] == "exact_event_localized_one_activation_packets"
    assert OUTER_EMISSION_LIMITER_SCHEMA.startswith("v10.0.5.18.3.2")
