from __future__ import annotations

import json
from pathlib import Path

import pytest

from reduced_fracture_v3 import (
    COMMON_VOID_KINETICS_ROW_ID,
    FRACTURE_ROWS,
    MaterialBundle,
    load_exact_fracture_rows,
    material_field_mapping_audit,
    paired_case_ledger,
    pilot_material_bundles,
    require_complete_material_mapping,
)


ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "analysis_outputs/oneD_v2_terminal_predictive_program/oneD_v2_new_four_class_registry.csv"
OUT = ROOT / "analysis_outputs/oneD_v3_aligned_temperature_transfer"


def test_01_exact_four_family_rows_are_frozen():
    rows = load_exact_fracture_rows(REGISTRY)
    assert {family: row["candidate_id"] for family, row in rows.items()} == FRACTURE_ROWS


def test_02_material_bundle_has_five_separate_versioned_owners():
    bundles = pilot_material_bundles()
    assert set(bundles) == {"Peak", "DBTT", "weak-T", "ceramic-like"}
    assert {bundle.void_kinetics_row_id for bundle in bundles.values()} == {
        COMMON_VOID_KINETICS_ROW_ID
    }
    for bundle in bundles.values():
        assert bundle.fracture_material_row_id != bundle.void_kinetics_row_id
        assert bundle.schema == "oneD.v3.material-bundle/1"


def test_03_bundle_rejects_conflated_fracture_and_void_ownership():
    with pytest.raises(ValueError, match="ownership"):
        MaterialBundle("same", "same", "elastic/1", "sites/1", "loading/1")


def test_04_mapping_audits_every_field_of_every_row_at_full_precision():
    rows = load_exact_fracture_rows(REGISTRY)
    audit = material_field_mapping_audit(
        rows, source_registry="registry.csv", v5_driver_identity="sha:driver"
    )
    expected = sum(len(row) for row in rows.values())
    assert len(audit["records"]) == expected
    peak = next(
        row for row in audit["records"]
        if row["material_class"] == "Peak" and row["source_field"] == "cleave_G00_eV"
    )
    assert peak["full_precision_value"] == "4.011803912930191"
    assert peak["units"] == "eV"
    assert peak["conversion"] == "identity"


def test_05_unmapped_active_v5_fields_fail_closed_before_execution():
    rows = load_exact_fracture_rows(REGISTRY)
    audit = material_field_mapping_audit(
        rows, source_registry="registry.csv", v5_driver_identity="sha:driver"
    )
    assert audit["summary"]["gate"] == "BLOCKED_UNMAPPED_ACTIVE_FIELD"
    assert not audit["summary"]["all_active_fields_mapped"]
    with pytest.raises(RuntimeError, match="BLOCKED_UNMAPPED_ACTIVE_FIELD"):
        require_complete_material_mapping(audit)


def test_06_unified_runtime_binding_is_exact_and_active():
    audit = json.loads((OUT / "material_field_mapping_audit.json").read_text())
    record = next(
        row for row in audit["records"]
        if row["material_class"] == "DBTT" and row["source_field"] == "c_blunt"
    )
    assert record["two_d_target"] == "MaterialManifest.c_blunt and FrontConfig.c_blunt"
    assert record["classification"] == "active"
    assert record["reason"] == "exact unified 2-D material-bundle runtime binding"
    assert audit["summary"]["gate"] == "PASS_EXACT_RUNTIME_BINDING"


def test_07_common_void_scout_is_exact_v5_code_and_material_independent():
    scout = json.loads((OUT / "common_void_rate_scout.json").read_text())
    assert scout["source_function"] == "arrhenius_fracture.voiding_v5.arrhenius_rates"
    assert scout["temperature_grid_K"] == list(range(300, 1201, 25))
    assert scout["material_dependence"] == "NONE_COMMON_REFERENCE_VOID_KINETICS_ROW"
    assert len(scout["rows"]) == 37


def test_08_twelve_case_ledger_is_present_and_none_were_run():
    ledger = paired_case_ledger(
        {family: (300.0, None, 1200.0) for family in FRACTURE_ROWS},
        gate="BLOCKED_UNMAPPED_ACTIVE_FIELD",
    )
    assert len(ledger) == 12
    assert {row["material_class"] for row in ledger} == set(FRACTURE_ROWS)
    assert all(row["void_2d_status"] == "NOT_RUN_M2_GATE_BLOCKED" for row in ledger)


def test_09_ceramic_like_is_complete_holdout_and_fatigue_remains_unstarted():
    decision = json.loads((OUT / "decision.json").read_text())
    assert decision["held_out_policy"] == {
        "ceramic_like_used_for_tuning": False,
        "complete_held_out_material_family": "ceramic-like",
        "development_sentinels": ["Peak", "DBTT", "weak-T"],
    }
    assert decision["pilot_terminal"]["FATIGUE_IMPLEMENTATION"] == "NOT_STARTED_BY_CONTRACT"
    assert decision["execution_counts"]["new_2d_mechanics_solves"] == 1


def test_10_preserved_geometric_and_no_inference_contracts():
    decision = json.loads((OUT / "decision.json").read_text())
    assert decision["preserved"]["r_tip_equals_R_void"] is False
    assert decision["preserved"]["fractured_length_equals_free_span_equals_front_coordinate"] is False
    assert decision["preserved"]["missing_mechanics_fields_inferred"] is False
    assert decision["preserved"]["fracture_rows_changed"] is False


def test_11_later_phases_remain_bounded_after_oracle_readiness_block():
    decision = json.loads((OUT / "decision.json").read_text())
    contract = decision["prospective_execution_contract"]
    assert contract["oracle_matrix"]["planned_unique_states"] == 18
    assert contract["oracle_matrix"]["status"] == "BLOCKED_CENTRAL_DBTT_V4_RESOLUTION_AND_QUALITY"
    assert contract["load_mapping"]["raw_2d_opening_used_as_1d_K"] is False
    assert contract["baseline_comparison"]["void_increment"] == (
        "Delta O_void = O_void - O_no_void"
    )
    assert contract["paired_state_exact_equality"] == [
        "thresholds",
        "RNG state",
        "site state",
        "initial void state",
        "fracture state",
        "r_tip/emission/shielding state",
        "load history",
        "temperature history",
        "stop condition",
    ]
    assert contract["complete_pass_rule"].startswith("all four material families")


def test_12_feature_temperatures_were_frozen_before_paired_results():
    scout = json.loads((OUT / "fracture_rate_feature_scout.json").read_text())
    decision = json.loads((OUT / "decision.json").read_text())
    assert scout["selected_feature_temperatures_K"] == {
        "Peak": 1100,
        "DBTT": 875,
        "weak-T": 950,
        "ceramic-like": 1100,
    }
    assert scout["selection_executed_before_paired_results"] is True
    assert scout["paired_material_temperature_cells_present_at_selection"] == 0
    assert scout["held_out_results_used_for_temperature_selection"] is False
    assert decision["feature_temperature_gate"] == "PASS_FROZEN_PROSPECTIVE_ANCHORS"
    assert all(
        value["T_feature_status"] == "FROZEN_PROSPECTIVELY_BEFORE_PAIRED_RESULTS"
        for value in decision["temperature_anchors"].values()
    )


def test_13_passed_m2_ledger_waits_only_for_aligned_oracle():
    anchors = {
        family: (
            300.0,
            json.loads((OUT / "decision.json").read_text())["temperature_anchors"][family]["T_feature_K"],
            1200.0,
        )
        for family in FRACTURE_ROWS
    }
    ledger = paired_case_ledger(anchors, gate="PASS_EXACT_RUNTIME_BINDING")
    assert all(row["temperature_K"] is not None for row in ledger)
    assert all(row["void_2d_status"] == "NOT_RUN_PENDING_ALIGNED_ORACLE" for row in ledger)


def test_14_oracle_readiness_fails_closed_on_expected_source_noncertification():
    readiness = json.loads((OUT / "aligned_oracle_readiness.json").read_text())
    decision = json.loads((OUT / "decision.json").read_text())
    assert readiness["source_head_matches_pinned_unified_commit"] is True
    assert readiness["schema"] == "oneD.v3.aligned-oracle-readiness/2"
    assert readiness["source_recovery_operator"] == "CAVITY_FIXED_ARC_PATCH_RECOVERY_V2"
    assert readiness["requested_refinement_levels"] == 3
    assert readiness["failure_classification"] == [
        "TANGENTIAL_STRESS_CONVERGENCE",
        "NORMAL_DIRECTION_RESOLUTION",
    ]
    assert readiness["readiness_gate_taxonomy"] == {
        "decision": "B_STILL_MANDATORY_V2_GATES",
        "eta_n": "MANDATORY_V2_GATE",
        "eta_t": "MANDATORY_V2_GATE",
    }
    assert readiness["accepted_oracle_states"] == 0
    assert readiness["oracle_matrix_status"] == "BLOCKED_DOWNSTREAM_SOURCE_TENSOR_UNQUALIFIED"
    assert readiness["expected_scientific_noncertification_converted_to_unavailable"] is True
    assert readiness["unexpected_programming_exceptions_caught"] is False
    assert readiness["missing_fields_inferred_or_synthesized"] is False
    assert readiness["final_governing_predicates"] == {
        "fixed_arc_tensor_convergence": False,
        "cavity_traction": True,
        "normal_direction_resolution": False,
        "tangential_direction_resolution": True,
    }
    assert decision["historical_v2_readiness"] == {
        "artifact": "aligned_oracle_readiness.json",
        "status": "BLOCKED_DOWNSTREAM_SOURCE_TENSOR_UNQUALIFIED",
        "taxonomy": "B_STILL_MANDATORY_V2_GATES",
        "failure_classification": [
            "TANGENTIAL_STRESS_CONVERGENCE",
            "NORMAL_DIRECTION_RESOLUTION",
        ],
    }
    assert decision["execution_counts"]["paired_material_temperature_cells_run"] == 0


def test_15_v3_source_readiness_stops_at_exact_geometry_failure_without_oracle():
    readiness = json.loads((OUT / "v3_source_readiness.json").read_text())
    decision = json.loads((OUT / "decision.json").read_text())
    assert readiness["source_head"] == "67528be7e38deb7047ad630ac3497bf9926543e1"
    assert readiness["source_recovery_operator"] == (
        "CAVITY_FIXED_PHYSICAL_ARC_PATCH_RECOVERY_V3"
    )
    central = readiness["central_dbtt_v3"]
    assert central["DBTT_SOURCE_READINESS"] == "BLOCKED_WITH_EXACT_V3_FAILURE_CLASS"
    assert central["exact_v3_failure_class"] == ["SOURCE_GEOMETRY_IDENTITY"]
    assert central["refinement_levels_run"] == 0
    assert central["accepted_pre_source_state_unchanged_on_noncertification"] is True
    assert central["thresholds_and_rng_unchanged_on_noncertification"] is True
    assert central["unexpected_programming_exceptions_caught"] is False
    assert readiness["oracle_states_accepted"] == 0
    assert readiness["paired_trajectories_run"] == 0
    assert readiness["fatigue_started"] is False
    assert readiness["missing_fields_inferred_or_synthesized"] is False
    assert readiness["next_bounded_step"] == "DERIVE_FINITE_ACTIVATION_ZONE_WORK_OBSERVABLE"
    assert readiness["bounded_worker"]["run_id"] == 34719718024
    assert readiness["bounded_worker"]["tests_passed"] == 22
    assert decision["v3_source_readiness"]["DBTT_SOURCE_READINESS"] == (
        "BLOCKED_WITH_EXACT_V3_FAILURE_CLASS"
    )
    assert decision["v3_source_readiness"]["exact_v3_failure_class"] == [
        "SOURCE_GEOMETRY_IDENTITY"
    ]
    assert decision["execution_counts"]["v3_source_refinement_levels_run"] == 0


def test_16_v4_source_geometry_passes_but_resolution_and_quality_block_oracle():
    readiness = json.loads((OUT / "v4_source_readiness.json").read_text())
    decision = json.loads((OUT / "decision.json").read_text())
    assert readiness["source_head"] == "9d2cb7cb63fc74f02e5f0b54b5a2d1a394fe7d05"
    assert readiness["source_geometry_contract"] == "CAVITY_SOURCE_CONFORMING_GEOMETRY_V4"
    central = readiness["central_dbtt_v4"]
    assert central["DBTT_SOURCE_READINESS"] == "BLOCKED_WITH_EXACT_V4_FAILURE_CLASS"
    assert central["exact_v4_failure_class"] == [
        "NORMAL_DIRECTION_RESOLUTION", "TANGENTIAL_DIRECTION_RESOLUTION", "MESH_QUALITY"
    ]
    assert central["levels_run"] == 3
    final = central["level_records"][-1]
    assert final["predicates"]["source_geometry"] is True
    assert final["predicates"]["tensor_convergence"] is True
    assert final["predicates"]["cavity_traction"] is True
    assert final["predicates"]["patch_conditioning"] is True
    assert final["predicates"]["normal_direction_resolution"] is False
    assert final["predicates"]["tangential_direction_resolution"] is False
    assert final["predicates"]["minimum_mesh_quality"] is False
    assert central["accepted_material_identity_exact_across_levels"] is True
    assert central["physical_geometry_exact_across_levels"] is True
    assert central["thresholds_and_rng_exact_across_levels"] is True
    assert central["accepted_input_state_unchanged_on_noncertification"] is True
    assert readiness["oracle_states_accepted"] == 0
    assert readiness["paired_trajectories_run"] == 0
    assert readiness["fatigue_started"] is False
    assert readiness["mechanics_map_fitting_complete"] is False
    assert readiness["monotonic_2d_transfer_complete"] is False
    assert readiness["finite_activation_zone_observable_derived"] is False
    assert readiness["bounded_worker"]["run_id"] == 34723836462
    assert readiness["bounded_worker"]["artifact_id"] == 10307550298
    assert readiness["bounded_worker"]["tests_passed"] == 36
    assert decision["blocking_reason"]["code"] == "CENTRAL_DBTT_V4_RESOLUTION_AND_QUALITY"
    assert decision["execution_counts"]["v4_angular_levels_run"] == 3
    assert decision["prospective_execution_contract"]["oracle_matrix"]["status"] == (
        "BLOCKED_CENTRAL_DBTT_V4_RESOLUTION_AND_QUALITY"
    )
