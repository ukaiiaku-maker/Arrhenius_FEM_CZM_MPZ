from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

import arrhenius_fracture.persistent_site_load_ramp_stochastic_emission_v10051831 as ramp31
from arrhenius_fracture.four_class_parameter_bridge_v100518 import (
    load_four_class_parameter_option,
)
from arrhenius_fracture.persistent_site_event_driven_emission_v10051831 import (
    EVENT_DRIVEN_SCHEMA,
    PersistentSiteEventDrivenSingleLayerLoadRampFrontEngineV10051831,
)
from arrhenius_fracture.persistent_site_load_ramp_stochastic_emission_v10051831 import (
    LOAD_RAMP_SCHEMA,
    OUTER_EMISSION_LIMITER_SCHEMA,
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


def _engine(monkeypatch, factors=(1.0, 0.5)):
    monkeypatch.setenv("EMISSION_HAZARD_SEED", "331085649")
    monkeypatch.setenv("EMISSION_RAMP_MAX_LOG_RATE_CHANGE", "0.5")
    monkeypatch.setenv("EMISSION_RAMP_ACTION_FLOOR", "1e-6")
    monkeypatch.setenv("EMISSION_INNER_MAX_LOG_HAZARD_CHANGE", "0.5")
    monkeypatch.setenv("EMISSION_INNER_ACTION_FLOOR", "1e-6")
    monkeypatch.setenv("EMISSION_EVENT_HORIZON_FACTOR", "1.25")
    candidate, _ = load_four_class_parameter_option(PARAMETER_ROOT, OPTION)
    family = load_signed_shielding_artifact_v1005141(FAMILY)
    cls = PersistentSiteEventDrivenSingleLayerLoadRampFrontEngineV10051831
    cls.configure(candidate, family)
    cls.configure_stochastic(
        hazard_mode="exponential",
        hazard_seed=331085649,
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


def test_outer_ramp_does_not_duplicate_inner_emission_action_limiter(monkeypatch):
    engine = _engine(monkeypatch)
    engine.load_ramp_K_cleave_prev_Pa_sqrt_m = 5.0e6
    engine.load_ramp_K_emit_prev_Pa_sqrt_m = 5.0e6
    monkeypatch.setattr(engine, "lambda_cleave", lambda sigma, T: (0.0, 0.0, 0.0))
    monkeypatch.setattr(
        ramp31,
        "_emission_hazards",
        lambda *args, **kwargs: {
            "aggregate_hazard_s": np.asarray([1.0e6, 1.0e6], dtype=float)
        },
    )

    h, diagnostics = engine._local_ramp_substep(
        remaining_s=10.0,
        consumed_s=0.0,
        total_s=10.0,
        K_cleave_end=5.0e6,
        K_emit_end=5.0e6,
        T_K=300.0,
        drive_factors=np.asarray([1.0, 0.5]),
        tau_signed_Pa=np.asarray([1.0e8, -1.0e8]),
    )

    assert h == pytest.approx(10.0)
    assert diagnostics["emission_action"] == pytest.approx(1.0e7)
    assert diagnostics["outer_emission_action_limiter_active"] == 0.0


def test_actual_weakT_startup_predictor_remains_finite_and_nonmutating(monkeypatch):
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


def test_bounded_maximum_two_channel_factor_predictor_is_finite(monkeypatch):
    engine = _engine(monkeypatch, factors=(5.0, 5.0))
    predicted = engine.predict_clock_increment_drives(
        80.0e6,
        80.0e6,
        300.0,
        840.0,
    )
    assert np.isfinite(predicted)
    assert 0.0 <= predicted <= 1.0 + 1.0e-12


def test_audit_records_single_layer_action_control(monkeypatch):
    engine = _engine(monkeypatch)
    audit = type(engine).audit_payload()["stochastic_emission"]
    assert audit["load_ramp_schema"] == LOAD_RAMP_SCHEMA
    assert audit["outer_emission_action_limiter_active"] is False
    assert audit["outer_emission_limiter_schema"] == OUTER_EMISSION_LIMITER_SCHEMA
    assert audit["inner_exact_emission_action_control"] is True
    assert audit["event_driven_transport_active"] is True
    assert audit["event_driven_transport_schema"] == EVENT_DRIVEN_SCHEMA
    assert audit["fixed_inner_action_substep_per_event"] is False
    assert audit["threshold_crossings_localized_individually"] is True
    assert audit["events_batched"] is False
    assert audit["global_FEM_state_change_predictor"] is False
