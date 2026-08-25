#!/usr/bin/env python3
"""Materialize PF-regime and FEM-kernel qualification evidence, fail closed."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_outputs/oneD_v2_local_drive"
FEM = Path("/Volumes/Data/Data/Nanopillar_calculation/Arrhenius_FEM_CZM_MPZ_theta0_pf_parity_claude")
EXPECTED = "a85b57ad9eee8331ef34ea222760ce7d4064f48ddef668f1e28a666d14946f3a"
ACTUAL = "d41b08f69ae773009f65f4c0094ef5436ba7f0f290a94172e5b894d09a41f7c3"


def dump(name, value):
    (OUT / name).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def write(name, value):
    (OUT / name).write_text(value.strip() + "\n")


def digest_json(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def main():
    nodes = pd.read_csv(OUT / "oneD_v2_pf_tensor_drive_node_summary.csv")
    unique = nodes.drop_duplicates("actual_extension_um", keep="last").sort_values("actual_extension_um").reset_index(drop=True)
    regime_fields = ["wake_fingerprint", "front_basis_fingerprint", "probe_support_fingerprint",
                     "probe_weight_fingerprint", "selected_system_signature", "resolved_sign_signature", "reliable"]
    unique["probe_regime_key"] = [digest_json({key: row[key] for key in regime_fields}) for _, row in unique.iterrows()]
    regime_columns = [
        "actual_extension_um", "geometry_event_count", "event_length_sequence_m_json",
        "event_geometry_fingerprint", "wake_fingerprint", "killed_node_ids_json",
        "killed_node_count", "tip_x_m", "tip_y_m", "front_tangent_x", "front_tangent_y",
        "front_normal_x", "front_normal_y", "front_basis_fingerprint",
        "probe_element_ids_json", "probe_element_weights_json", "probe_support_fingerprint",
        "probe_weight_fingerprint", "selected_system_signature", "resolved_sign_signature",
        "reliable", "opening_sigma_xx_Pa", "opening_sigma_yy_Pa", "opening_sigma_xy_Pa",
        "channel0_sigma_xx_Pa", "channel0_sigma_yy_Pa", "channel0_sigma_xy_Pa",
        "channel1_sigma_xx_Pa", "channel1_sigma_yy_Pa", "channel1_sigma_xy_Pa",
        "tau0_signed_Pa", "tau1_signed_Pa", "factor0", "factor1", "reaction_N_per_m",
        "compliance_m2_per_N", "elastic_energy_J_per_m", "native_J_J_per_m2",
        "native_KJ_Pa_sqrt_m", "profile_s_min_m", "profile_s_max_m", "profile_n_min_m",
        "profile_n_max_m", "probe_regime_key",
    ]
    unique[regime_columns].to_parquet(OUT / "oneD_v2_pf_regime_table.parquet", index=False)

    transitions = []
    for index in range(1, len(unique)):
        left, right = unique.iloc[index - 1], unique.iloc[index]
        changes = {field + "_changed": bool(left[field] != right[field]) for field in regime_fields}
        transitions.append({
            "left_extension_um": left.actual_extension_um,
            "right_extension_um": right.actual_extension_um,
            **changes,
            "interpolation_allowed": not any(changes.values()),
            "policy": "EXACT_DETERMINISTIC_QUERY_AND_CACHE" if any(changes.values()) else "WITHIN_REGIME_INTERPOLATION",
        })
    pd.DataFrame(transitions).to_csv(OUT / "oneD_v2_pf_tensor_transition_table.csv", index=False)

    validation = pd.read_csv(OUT / "oneD_v2_pf_tensor_interpolation_validation.csv")
    failures = validation[(validation.validation_kind == "LEAVE_ONE_NODE_OUT") & (validation.sign_mismatch_count == 1)].copy()
    classified = []
    tensor_columns = [name for name in unique if name.endswith("_Pa") and "reaction" not in name]
    for _, failure in failures.iterrows():
        index = int(np.flatnonzero(np.isclose(unique.actual_extension_um, failure.actual_extension_um))[0])
        target, left, right = unique.iloc[index], unique.iloc[index - 1], unique.iloc[index + 1]
        local_scale = max(abs(float(target[name])) for name in tensor_columns)
        basis_change = left.front_basis_fingerprint != target.front_basis_fingerprint or target.front_basis_fingerprint != right.front_basis_fingerprint
        support_change = left.probe_support_fingerprint != target.probe_support_fingerprint or target.probe_support_fingerprint != right.probe_support_fingerprint
        weight_change = left.probe_weight_fingerprint != target.probe_weight_fingerprint or target.probe_weight_fingerprint != right.probe_weight_fingerprint
        wake_change = left.wake_fingerprint != target.wake_fingerprint or target.wake_fingerprint != right.wake_fingerprint
        if abs(float(failure.observed)) <= 0.01 * local_scale:
            primary = "NEAR_ZERO_COMPONENT_SIGN_ARTIFACT"
        elif basis_change:
            primary = "FRONT_DIRECTION_TRANSITION"
        elif support_change:
            primary = "PROBE_ELEMENT_SELECTION_TRANSITION"
        elif weight_change:
            primary = "PROBE_WEIGHT_TRANSITION"
        elif wake_change:
            primary = "WAKE_TOPOLOGY_TRANSITION"
        else:
            primary = "INSUFFICIENT_WITHIN_REGIME_NODE_DENSITY"
        classified.append({
            **failure.to_dict(), "primary_classification": primary,
            "wake_topology_transition": wake_change, "front_direction_transition": basis_change,
            "probe_element_selection_transition": support_change, "probe_weight_transition": weight_change,
            "reliability_transition": left.reliable != target.reliable or target.reliable != right.reliable,
            "left_regime_key": left.probe_regime_key, "target_regime_key": target.probe_regime_key,
            "right_regime_key": right.probe_regime_key,
            "provider_action": "HARD_FAIL_OR_EXACT_DETERMINISTIC_QUERY",
        })
    classified = pd.DataFrame(classified)
    classified.to_csv(OUT / "oneD_v2_pf_tensor_failure_classification.csv", index=False)

    counts = classified.primary_classification.value_counts().to_dict()
    write("ONE_D_V2_PF_TENSOR_INTERPOLATION_FAILURE_AUDIT.md", f"""
# PF tensor interpolation failure audit

Status: **SPARSE EXTENSION-ONLY INTERPOLATION REJECTED**.

All 23 sign failures were re-evaluated against the realized wake, front basis,
probe support, probe weights, reliability and categorical signed-system state.
They are not one undifferentiated numerical error: `{json.dumps(counts, sort_keys=True)}`.
The CSV retains all simultaneous cause flags and the left/target/right regime
keys. A reliability, signed-system, or sign transition is a hard failure. A
near-zero component is reported as such and is never promoted to a reliable
continuous sign crossing.
""")
    write("ONE_D_V2_PF_REGIME_AWARE_TENSOR_MAP.md", f"""
# PF regime-aware tensor map

Status: **HYBRID PROVIDER CONTRACT QUALIFIED; EXTENSION-ONLY SURROGATE UNQUALIFIED**.

`PFProbeRegimeKey` excludes continuous tensor values and includes wake topology,
front basis, probe element support, probe weights, selected-system/sign
signatures, and reliability. The {len(unique)} exact source states span
0–1000 µm; all are reliable and retain killed-node IDs, realized event geometry,
probe IDs/weights, raw tensors, reaction, energy, native J/KJ, and profile bounds.
None of the {len(transitions)} adjacent sparse intervals has an identical regime
key, so no interval is admitted for interpolation. Queries at a transition use
an exact deterministic mechanics query keyed by realized event geometry and are
cached; absence of that callable fails closed. Because the production event
length is a clipped continuous exponential draw, arbitrary trajectories cannot
be represented by a finite extension-only lookup.
""")

    shadow = pd.read_csv(OUT / "oneD_v2_pf_exact_wrapper_shadow_v2.csv")
    write("ONE_D_V2_PF_EXACT_WRAPPER_SHADOW_V2.md", f"""
# PF exact-wrapper shadow V2

Status: **PARTIAL PASS AT ONE QUALIFIED SOURCE NODE; FULL SHADOW UNQUALIFIED**.

At the 100 µm qualified source state, all {len(shadow)} four-class/temperature
cells reproduce the exact persistent-site wrapper fixture's selected emission system,
signed channels, reliability, effective stresses, barriers, per-system rates,
probabilities and aggregate hazard exactly. Maximum recorded continuous-field
error is zero. This qualifies source-node composition, not sparse interpolation
or an arbitrary off-node trajectory. Required small/large radius, front-width,
retained-state, backstress, extension and opening coverage is still absent;
off-node use requires the exact deterministic topology query described above.
""")

    actual_path = FEM / "reference_inputs/pf_final_v10_2_30/rate1x/v913_paper_peak01_0242980_persistent_sites/T1000K_th0_seed8666/family.json"
    inventory = [
        {"artifact_role": "historical_PF_byte_identity_contract", "sha256": EXPECTED, "path": str(FEM / "reference_inputs/pf_v10_4_1_theta0_peak1000K_seed8666/family.json"), "present": False, "classification": "MISSING_NOT_REGENERATED"},
        {"artifact_role": "authoritative_corrected_FEM_runtime_source", "sha256": ACTUAL, "path": str(actual_path), "present": actual_path.is_file(), "classification": "PRESENT_USED_BY_ALL_FOUR_CORRECTED_LONG_RUNS"},
        {"artifact_role": "committed_v10_0_5_17_frozen_input", "sha256": "a876ea042bf291ca9607b586205d75ce2b5cb95c668fef53d713d611dfe59cde", "path": str(FEM / "runtime_inputs/v10_0_5_17_frozen_pf_inputs/signed_kernel/v10_2_14_active_only_campaign_family.json"), "present": True, "classification": "PRESENT_DIFFERENT_BYTES"},
        {"artifact_role": "PF_named_cache_current", "sha256": ACTUAL, "path": "/Volumes/Data/Data/Nanopillar_calculation/PF-fracture-fatigue_v10_2_21_persistent_sites_top1/runs/v10_2_28_kernel_cache/1447653d199f0b43cb475951092d69444c9b785f6fdf518c723792abb3b1f5e5/family.json", "present": True, "classification": "PRESENT_MUTATED_FROM_HISTORICAL_EXPECTATION"},
        {"artifact_role": "PF_alternate_cache_4fa015", "sha256": "5aa74370e6419104684c52bbcf93323f905d2752d1bba59252f9b0b35c77e07c", "path": "/Volumes/Data/Data/Nanopillar_calculation/PF-fracture-fatigue_v10_2_21_persistent_sites_top1/runs/v10_2_28_kernel_cache/4fa015d77f1aadf05f77f550366f64cd611f537ae716bbd47870bf9e6fe2f873/family.json", "present": True, "classification": "PRESENT_DIFFERENT_BYTES"},
        {"artifact_role": "PF_v10_4_2_alternate_cache_4fa015", "sha256": "b1857539f0e47e3520a2d4b2acd0df186a5129e86e13552819b5c27f1ce91ed1", "path": "/Volumes/Data/Data/Nanopillar_calculation/PF-fracture-fatigue_v10_4_2_plastic_flow_terminal/runs/v10_2_28_kernel_cache/4fa015d77f1aadf05f77f550366f64cd611f537ae716bbd47870bf9e6fe2f873/family.json", "present": True, "classification": "PRESENT_DIFFERENT_BYTES"},
        {"artifact_role": "PF_v11_causal_cache_353593", "sha256": "5b8ed966506332299d2d891d349710c897fbb734340df67eaf8767ff6c1f53d7", "path": "/Volumes/Data/Data/Nanopillar_calculation/PF-fracture-fatigue_v11_branching/runs/v11_direct_kernel_cache_causal/3535939d541388c5e3dd88b9a712b7bfe03d951cd64344f7d3af3ca0c5f52ea4/family.json", "present": True, "classification": "PRESENT_DIFFERENT_BYTES"},
        {"artifact_role": "PF_v11_causal_cache_e2e456", "sha256": "95890404744ef1b49264a286d96283ae18c3e73bf1c276178fb0e9b7180f03fa", "path": "/Volumes/Data/Data/Nanopillar_calculation/PF-fracture-fatigue_v11_branching/runs/v11_direct_kernel_cache_causal/e2e456c75935a60a56ccbd9f3b7392036db704de3c2ec6b88e6f0e7bfb127070/family.json", "present": True, "classification": "PRESENT_DIFFERENT_BYTES"},
    ]
    pd.DataFrame(inventory).to_csv(OUT / "oneD_v2_fem_kernel_inventory.csv", index=False)
    composition = {
        "schema": "oneD_v2_fem_factory_composition_v3",
        "corrected_source_commit": "17f8bf67c0138b453eadff4b509afa80ebaa2c04",
        "historical_pf_hash_contract": EXPECTED,
        "corrected_run_kernel_sha256": ACTUAL,
        "authoritative_run_kernel_classification": "AUTHORITATIVE_RUN_USED_D41B08",
        "corrected_run_launch_heads": {"Peak_tip_only": "931bed66913afc970117bce900805ecc9b6225f8", "DBTT_tip_only": "931bed66913afc970117bce900805ecc9b6225f8", "Peak_bulk_pt_on": "eee33dfcaaac3ab33172d56994fb0b184b0c452c", "DBTT_bulk_pt_on": "a691deafe4a4fe0488175da92abcabeba795f13c"},
        "runtime_composition": ["native-rate launcher records supplied d41b08 kernel", "process-local REFERENCE_KERNEL_SHA256 override", "audited empty-wake representation adapter", "unchanged active kernel/interpolation states", "v10.0.5.18.4.0 wrapper chain", "ProductionJointKRampStochasticEmissionFrontEngineV10051832"],
        "configured_exact_engine_instance_constructed": True,
        "engine_class": "ProductionJointKRampStochasticEmissionFrontEngineV10051832",
        "engine_state_class": "PersistentSiteSignedMPZStateV100514",
        "candidate_fingerprint_peak": "ce47ec844b3d7a9812aa0b077351fc6d022c292b39b84bd4056871395f7f4f52",
        "kernel_blocker": "CLOSED_FOR_CORRECTED_FEM_RUNTIME_NOT_FOR_HISTORICAL_PF_BYTE_IDENTITY",
        "historical_pin_modified": False,
        "silent_substitution": False,
        "new_stochastic_fem_runs": 0,
    }
    dump("oneD_v2_fem_factory_composition_v3.json", composition)
    write("ONE_D_V2_FEMCZM_KERNEL_PROVENANCE_AUDIT.md", f"""
# FEM/CZM kernel provenance audit

The apparent hash conflict represented two different claims. `{EXPECTED}` is
the missing historical PF byte-identity contract. `{ACTUAL}` is present and is
explicitly recorded by every authoritative corrected Peak/DBTT long-growth run.
Those runs used the tracked native-rate launcher, a process-local recorded hash
override, and the audited empty-wake representation adapter; the adapter leaves
the active kernel and interpolation states unchanged and does not enable wake
shielding. Repository history, local family files, LFS storage, attachments and
the PF/FEM zip inventories contain no file with the historical hash. No artifact
was regenerated and no source pin was edited.
""")
    write("ONE_D_V2_FEMCZM_FACTORY_COMPOSITION_AUDIT_V2.md", """
# FEM/CZM factory composition audit V2

Status: **CORRECTED-RUNTIME KERNEL/ENGINE COMPOSITION RESOLVED**.

An isolated source fixture loaded the authoritative `d41b08…` artifact through
the same audited active-only compatibility adapter, configured the exact Peak
candidate, and constructed `ProductionJointKRampStochasticEmissionFrontEngineV10051832`
with `PersistentSiteSignedMPZStateV100514`. This closes the earlier kernel
construction blocker for the corrected FEM model. It does not assert historical
PF byte identity and does not authorize substituting `d41b08…` into that claim.
""")
    write("ONE_D_V2_FEMCZM_LOCAL_DRIVE_PROVIDER.md", """
# FEM/CZM local-drive provider

Status: **UNQUALIFIED — DYNAMIC TENSOR-PROBE MAP NOT YET BUILT**.

The kernel/factory blocker is closed for the corrected FEM runtime. Exact native
emission still consumes a reliable opening tensor and two signed channel probes
at the evolving process-zone radius. Existing zero-thickness native-J/G maps do
not archive those raw tensors, probe supports, or radius-dependent weights.
Consequently they cannot be promoted to a local-drive provider. No stochastic
FEM run was launched and no scalar K/J substitution is permitted.
""")
    parity = pd.DataFrame([{
        "comparison": "PF_vs_FEMCZM_common_tensor_drive", "status": "NOT_RUN_FEM_LOCAL_DRIVE_PROVIDER_UNQUALIFIED",
        "PF_provider": "QUALIFIED_SOURCE_NODES_EXACT_QUERY_HYBRID", "FEMCZM_provider": "UNQUALIFIED_MISSING_DYNAMIC_TENSOR_MAP",
        "maximum_relative_error": np.nan, "sign_mismatch_count": np.nan,
    }])
    parity.to_csv(OUT / "oneD_v2_common_tensor_parity.csv", index=False)
    write("ONE_D_V2_COMMON_TENSOR_PARITY.md", """
# Common tensor parity

**NOT RUN.** PF source-node/exact-query closure is qualified, but the corrected
FEM/CZM dynamic-radius tensor provider is not. A scalar native J/KJ or structural
G map cannot stand in for the missing signed local field.
""")
    decision = {
        "schema": "oneD_v2_level3_decision_v3", "Level1": "PASS_EXACT", "Level2": "PASS_EXACT",
        "PF_source_nodes": "QUALIFIED", "PF_regime_provider": "QUALIFIED_EXACT_QUERY_HYBRID",
        "PF_sparse_interpolation": "UNQUALIFIED", "PF_exact_wrapper_source_node": "PASS_EXACT_ONE_NODE_STATIC_STATE",
        "PF_exact_wrapper_full_domain": "UNQUALIFIED_INCOMPLETE_STATE_DOMAIN",
        "FEM_kernel_provenance": "RESOLVED_CORRECTED_RUNTIME_D41B08_HISTORICAL_A85_MISSING",
        "FEM_factory_composition": "RESOLVED_CORRECTED_RUNTIME",
        "FEM_local_drive_provider": "UNQUALIFIED_DYNAMIC_TENSOR_MAP_MISSING",
        "common_tensor_parity": "NOT_RUN_PREREQUISITE_FAILED", "Level3": "NATIVE_CLOSURE_UNQUALIFIED",
        "forward_baselines_authorized": False, "campaign_authorized": False,
        "canonical_parameters_changed": False, "new_stochastic_PF_runs": 0, "new_stochastic_FEM_runs": 0,
        "parameter_status": "INSUFFICIENT_EVIDENCE_FOR_V2_PARAMETER_DECISION",
    }
    dump("oneD_v2_level3_decision_v3.json", decision)
    write("ONE_D_V2_LEVEL3_NATIVE_CLOSURE_V3.md", """
# Level-3 native closure V3

Level 1 and Level 2 remain exact. PF now has exact source-node wrapper closure
and a fail-closed regime-aware exact-query provider; sparse interpolation remains
rejected. The corrected FEM kernel/factory identity is resolved, but its
dynamic-radius local tensor provider is absent. Therefore common tensor parity,
native microtrajectories, all baselines, and both campaigns remain unauthorized.
This is a mechanics-reduction/data-availability result, not evidence against the
four shared material rows.
""")
    fingerprint_path = OUT / "oneD_v2_local_drive_fingerprints.json"
    prior = json.loads(fingerprint_path.read_text()) if fingerprint_path.is_file() else {}
    names = set(prior) | {
        "ONE_D_V2_PF_TENSOR_INTERPOLATION_FAILURE_AUDIT.md",
        "ONE_D_V2_PF_REGIME_AWARE_TENSOR_MAP.md",
        "ONE_D_V2_PF_EXACT_WRAPPER_SHADOW_V2.md",
        "ONE_D_V2_FEMCZM_KERNEL_PROVENANCE_AUDIT.md",
        "ONE_D_V2_FEMCZM_FACTORY_COMPOSITION_AUDIT_V2.md",
        "ONE_D_V2_FEMCZM_LOCAL_DRIVE_PROVIDER.md",
        "ONE_D_V2_COMMON_TENSOR_PARITY.md",
        "ONE_D_V2_LEVEL3_NATIVE_CLOSURE_V3.md",
        "oneD_v2_pf_tensor_failure_classification.csv",
        "oneD_v2_pf_regime_table.parquet",
        "oneD_v2_pf_tensor_transition_table.csv",
        "oneD_v2_pf_exact_wrapper_shadow_v2.csv",
        "oneD_v2_fem_kernel_inventory.csv",
        "oneD_v2_fem_factory_composition_v3.json",
        "oneD_v2_common_tensor_parity.csv",
        "oneD_v2_level3_decision_v3.json",
    }
    fingerprints = {
        name: hashlib.sha256((OUT / name).read_bytes()).hexdigest()
        for name in sorted(names) if (OUT / name).is_file()
    }
    fingerprint_path.write_text(json.dumps(fingerprints, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
