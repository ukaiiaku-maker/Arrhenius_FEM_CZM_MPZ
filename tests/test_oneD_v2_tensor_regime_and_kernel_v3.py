import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from reduced_fracture_v2.local_drive import (
    DriveQualificationError,
    PFProbeRegimeKey,
    PFTensorDriveProvider,
)


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_outputs/oneD_v2_local_drive"


def provider(**kwargs):
    return PFTensorDriveProvider.from_artifacts(
        OUT / "oneD_v2_pf_tensor_drive_node_summary.csv",
        OUT / "oneD_v2_pf_tensor_drive_map_manifest.json",
        allow_source_nodes_only=True,
        **kwargs,
    )


def test_interpolation_never_crosses_incompatible_pf_regime_keys():
    with pytest.raises(DriveQualificationError, match="exact deterministic query"):
        provider().evaluate(7.5e-6, 1e-5)


def test_regime_transition_uses_explicit_exact_value_and_cache():
    source = provider().evaluate(10e-6, 1e-5, exact_node_only=True)
    calls = []

    def exact(extension, opening, key, geometry):
        calls.append((extension, opening, key, geometry))
        return source

    p = provider(exact_query=exact)
    assert p.evaluate(7.5e-6, 1e-5, geometry_state_fingerprint="event-path-a") is source
    assert p.evaluate(7.5e-6, 1e-5, geometry_state_fingerprint="event-path-a") is source
    assert len(calls) == 1


def test_exact_cache_cannot_alias_two_continuous_event_geometries():
    source = provider().evaluate(10e-6, 1e-5, exact_node_only=True)
    calls = []

    def exact(extension, opening, key, geometry):
        calls.append(geometry)
        return source

    p = provider(exact_query=exact)
    p.evaluate(7.5e-6, 1e-5, geometry_state_fingerprint="event-path-a")
    p.evaluate(7.5e-6, 1e-5, geometry_state_fingerprint="event-path-b")
    assert calls == ["event-path-a", "event-path-b"]


def test_reliability_change_is_a_hard_boundary():
    p = provider()
    row = p.nodes.iloc[0].copy()
    key = PFProbeRegimeKey.from_row(row)
    incompatible = PFProbeRegimeKey(
        key.wake_topology_fingerprint, key.front_basis_fingerprint,
        key.probe_support_fingerprint, key.probe_weight_fingerprint,
        key.selected_system_signature, key.resolved_sign_signature, not key.reliable,
    )
    with pytest.raises(DriveQualificationError, match="exact deterministic query"):
        p.evaluate(float(row.actual_extension_um) * 1e-6, 1e-5, regime_key=incompatible)


def test_signed_shear_and_selected_system_match_exact_pf_source_fixture():
    shadow = pd.read_csv(OUT / "oneD_v2_pf_exact_wrapper_shadow_v2.csv")
    assert len(shadow) == 12
    assert shadow.tau_sign_exact.all() and shadow.reliability_exact.all()
    assert shadow.selected_system_exact.all()


def test_barrier_rate_probability_and_hazard_errors_are_validated():
    shadow = pd.read_csv(OUT / "oneD_v2_pf_exact_wrapper_shadow_v2.csv")
    assert (shadow.maximum_barrier_absolute_error_J == 0).all()
    assert (shadow.maximum_rate_relative_error == 0).all()
    assert (shadow.maximum_probability_absolute_error == 0).all()
    assert (shadow.aggregate_hazard_relative_error == 0).all()


def test_near_zero_sign_failure_has_absolute_classification():
    failures = pd.read_csv(OUT / "oneD_v2_pf_tensor_failure_classification.csv")
    assert len(failures) == 23
    assert failures.primary_classification.nunique() > 1
    near = failures[failures.primary_classification == "NEAR_ZERO_COMPONENT_SIGN_ARTIFACT"]
    assert len(near) > 0 and np.isfinite(near.absolute_error).all()


def test_out_of_map_and_out_of_profile_fail_closed():
    with pytest.raises(DriveQualificationError, match="outside the PF tensor map"):
        provider().evaluate(1001e-6, 1e-5)
    from tests.test_oneD_v2_local_drive import _bundle
    with pytest.raises(DriveQualificationError, match="outside the qualified profile"):
        _bundle(radius=3e-6).require_radius_in_domain()


def test_raw_profile_and_regime_audit_fields_remain_available():
    regime = pd.read_parquet(OUT / "oneD_v2_pf_regime_table.parquet")
    required = {
        "killed_node_ids_json", "probe_element_ids_json", "probe_element_weights_json",
        "opening_sigma_xx_Pa", "channel0_sigma_xy_Pa", "channel1_sigma_xy_Pa",
        "profile_s_min_m", "profile_s_max_m", "profile_n_min_m", "profile_n_max_m",
    }
    assert required <= set(regime.columns)
    assert len(regime) == 17 and regime.reliable.all()


def test_event_geometry_is_continuous_and_exact_lookup_policy_is_explicit():
    transitions = pd.read_csv(OUT / "oneD_v2_pf_tensor_transition_table.csv")
    assert not transitions.interpolation_allowed.any()
    assert (transitions.policy == "EXACT_DETERMINISTIC_QUERY_AND_CACHE").all()
    report = (OUT / "ONE_D_V2_PF_REGIME_AWARE_TENSOR_MAP.md").read_text()
    assert "clipped continuous exponential draw" in report


def test_fem_kernel_is_not_silently_substituted_and_run_provenance_wins():
    inventory = pd.read_csv(OUT / "oneD_v2_fem_kernel_inventory.csv")
    historical = inventory[inventory.artifact_role == "historical_PF_byte_identity_contract"].iloc[0]
    corrected = inventory[inventory.artifact_role == "authoritative_corrected_FEM_runtime_source"].iloc[0]
    assert not historical.present and historical.sha256.startswith("a85b57")
    assert corrected.present and corrected.sha256.startswith("d41b08")
    composition = json.loads((OUT / "oneD_v2_fem_factory_composition_v3.json").read_text())
    assert composition["historical_pin_modified"] is False
    assert composition["silent_substitution"] is False
    assert composition["configured_exact_engine_instance_constructed"] is True


def test_level3_stays_closed_without_fem_tensor_provider_or_common_parity():
    decision = json.loads((OUT / "oneD_v2_level3_decision_v3.json").read_text())
    assert decision["Level1"] == decision["Level2"] == "PASS_EXACT"
    assert decision["PF_exact_wrapper_full_domain"].startswith("UNQUALIFIED")
    assert decision["FEM_local_drive_provider"].startswith("UNQUALIFIED")
    assert decision["common_tensor_parity"].startswith("NOT_RUN")
    assert decision["Level3"] == "NATIVE_CLOSURE_UNQUALIFIED"
    assert not decision["forward_baselines_authorized"] and not decision["campaign_authorized"]
    assert decision["canonical_parameters_changed"] is False
    assert decision["new_stochastic_FEM_runs"] == 0
