from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from arrhenius_fracture.four_class_parameter_bridge_v100518 import (
    load_four_class_parameter_option,
)
from arrhenius_fracture.persistent_site_state_coupled_stochastic_emission_v1005182 import (
    PersistentSiteStateCoupledStochasticEmissionFrontEngineV1005182,
    two_channel_absolute_opening_drives,
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
    monkeypatch.setenv("EMISSION_HAZARD_SEED", "444111166")
    monkeypatch.setenv("EMISSION_ADAPT_REFERENCE_TARGET", "0.05")
    monkeypatch.setenv("EMISSION_ADAPT_MAX_LOG_HAZARD_CHANGE", "0.25")
    monkeypatch.setenv("EMISSION_ADAPT_MAX_BACKSTRESS_FRACTION", "0.05")
    monkeypatch.setenv("EMISSION_ADAPT_MAX_GEOMETRY_FRACTION", "0.05")
    candidate, _ = load_four_class_parameter_option(PARAMETER_ROOT, OPTION)
    family = load_signed_shielding_artifact_v1005141(FAMILY)
    cls = PersistentSiteStateCoupledStochasticEmissionFrontEngineV1005182
    cls.configure(candidate, family)
    cls.configure_stochastic(
        hazard_mode="exponential",
        hazard_seed=444111166,
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
    return engine, candidate


def test_persistent_two_channel_drive_does_not_apply_scalar_emission_factor_twice():
    engine = SimpleNamespace(
        persistent_site_source_active=True,
        _mm=SimpleNamespace(
            latest={
                "cleavage_factor": 1.25,
                "emission_factor": 4.0,
                "two_channel_drive_reliable": True,
            }
        ),
    )
    Kc, Ke, metadata = two_channel_absolute_opening_drives(engine, 20.0e6)
    assert Kc == pytest.approx(25.0e6)
    assert Ke == pytest.approx(20.0e6)
    assert metadata["persistent_two_channel_absolute_opening_drive"] is True
    assert metadata["scalar_emission_factor_applied_to_Kemit"] is False
    assert metadata["two_channel_factor_applied_inside_emission_hazard"] is True


def test_nonpersistent_path_retains_inherited_scalar_emission_factor():
    engine = SimpleNamespace(
        persistent_site_source_active=False,
        _mm=SimpleNamespace(
            latest={
                "cleavage_factor": 1.25,
                "emission_factor": 4.0,
                "two_channel_drive_reliable": True,
            }
        ),
    )
    Kc, Ke, metadata = two_channel_absolute_opening_drives(engine, 20.0e6)
    assert Kc == pytest.approx(25.0e6)
    assert Ke == pytest.approx(80.0e6)
    assert metadata["scalar_emission_factor_applied_to_Kemit"] is True


def test_weakT_reference_site_multiplicity_uses_areal_density(monkeypatch):
    engine, candidate = _engine(monkeypatch)
    geometry = engine.mpz_state.source_geometry()
    expected = candidate.rho_source0_m2 * candidate.reference_source_area_um2 * 1.0e-12
    assert geometry["multiplicity_per_system"] == pytest.approx(expected)
    assert geometry["multiplicity_per_system"] == pytest.approx(472.38905390836175)
    assert candidate.source_sites_per_system_provenance == pytest.approx(
        141.0590567476921
    )
    assert geometry["multiplicity_per_system"] != pytest.approx(
        candidate.source_sites_per_system_provenance
    )


def test_large_weakT_startup_increment_requests_load_subdivision_without_mutation(
    monkeypatch,
):
    engine, _ = _engine(monkeypatch)
    before = engine.mpz_state.state_dict()
    event_count_before = engine.emission_event_count_total

    predicted = engine.predict_clock_increment_drives(
        80.0e6,
        80.0e6,
        300.0,
        840.0,
    )

    assert np.isfinite(predicted)
    assert predicted > engine.emission_adaptive_reference_target
    assert engine.mpz_state.state_dict() == before
    assert engine.emission_event_count_total == event_count_before
    audit = engine.emission_last_adaptive_prediction
    assert audit["predicted_source_activations"] > 1.0e4
    assert audit["backstress_fraction"] > 0.05
    assert audit["geometry_fraction"] > 0.05
    assert audit["adaptive_increment"] == pytest.approx(predicted)
    assert engine.emission_adaptive_rejection_requests == 1


def test_small_increment_can_reach_exact_stochastic_predictor(monkeypatch):
    engine, _ = _engine(monkeypatch)
    before = engine.mpz_state.state_dict()
    predicted = engine.predict_clock_increment_drives(
        1.0e3,
        1.0e3,
        300.0,
        1.0e-9,
    )
    assert np.isfinite(predicted)
    assert engine.mpz_state.state_dict() == before
    assert engine.emission_event_count_total == 0


def test_audit_rejects_post_onset_mean_field_burst(monkeypatch):
    engine, _ = _engine(monkeypatch)
    audit = type(engine).audit_payload()["stochastic_emission"]
    assert audit["accepted_update"] == "exact_event_localized_one_activation_packets"
    assert audit["post_onset_mean_field_burst"] is False
    assert audit["scalar_emission_factor_applied_to_opening_K"] is False
    assert audit["signed_channel_factors_applied_once"] is True
    assert audit["backstress_recomputed_after_each_event"] is True
    assert audit["front_width_recomputed_after_each_event"] is True
    assert audit["site_multiplicity_recomputed_after_each_event"] is True
