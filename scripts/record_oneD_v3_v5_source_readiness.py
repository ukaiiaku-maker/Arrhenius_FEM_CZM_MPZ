#!/usr/bin/env python3
"""Import the bounded V5 source decision without launching a 2-D solve."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess


ROOT = Path(__file__).resolve().parents[1]
V5_HEAD = "b6224cbdc41b76048d7ad4b459046fd8461f53e5"
V5_BRANCH = "codex/v5-unified-fracture-multitip-voiding"
SOURCE = Path("artifacts/v5_cavity_source_recovery_v5/central_dbtt_v5_readiness.json")
ATTESTATION = Path("artifacts/v5_cavity_source_recovery_v5/attestation.json")
OUTPUT = Path("analysis_outputs/oneD_v3_aligned_temperature_transfer/v5_source_readiness.json")
DECISION = Path("analysis_outputs/oneD_v3_aligned_temperature_transfer/decision.json")
FAILURE = ["RAW_ADJACENT_ELEMENT_TRACTION"]


def _git(repo: Path, *args: str) -> str:
    return subprocess.check_output(("git", *args), cwd=repo, text=True).strip()


def _level_summary(row):
    keys = (
        "N_theta", "local_level", "scientific_unavailability",
        "accepted_input_state_fingerprint", "accepted_input_state_unchanged",
        "material_identity_matches_reference", "thresholds_match_reference",
        "rng_state_matches_reference", "load_history_matches_reference",
        "source_coordinate_matches_reference", "nominal_cavity_geometry_matches_reference",
        "source_coordinate_m", "discrete_cavity_boundary_fingerprints",
        "constrained_boundary_limit_tensor_global_xy_Pa",
        "constrained_boundary_limit_tensor_local_nt_Pa",
        "raw_adjacent_element_traction_normalized",
        "assembled_weak_cavity_boundary_residual_normalized",
        "constrained_boundary_limit_traction_residual", "constrained_fit_residual",
        "patch_rank", "patch_condition", "eta_n", "eta_t",
        "local_minimum_quality", "global_minimum_quality", "physical_window_identity",
        "source_node_certificate", "sigma_zz_and_mean_stress_audit",
        "observables", "ligament_energy_gate", "candidate_kinetic_drive",
        "candidate_kinetic_drive_status",
    )
    return {key: row.get(key) for key in keys}


def build(v5_repo: Path):
    head = _git(v5_repo, "rev-parse", "HEAD")
    branch = _git(v5_repo, "branch", "--show-current")
    if head != V5_HEAD or branch != V5_BRANCH:
        raise RuntimeError("2-D source is not at the exact final V5 branch head")
    if _git(v5_repo, "status", "--porcelain=v1"):
        raise RuntimeError("2-D source worktree is not clean")
    source = json.loads((v5_repo / SOURCE).read_text())
    attestation = json.loads((v5_repo / ATTESTATION).read_text())
    if source["DBTT_SOURCE_READINESS"] != "BLOCKED_WITH_EXACT_V5_FAILURE_CLASS":
        raise RuntimeError("unexpected central DBTT V5 decision")
    if source["exact_v5_failure_class"] != FAILURE or source["oracle_states_accepted"] != 0:
        raise RuntimeError("unexpected central DBTT V5 failure or oracle count")
    worker = attestation["bounded_worker"]
    if (worker["conclusion"], worker["job_count"], worker["tests_passed"]) != (
        "success", 1, 46,
    ):
        raise RuntimeError("bounded V5 worker is not successful")
    angular = source["angular_family"]
    local = source["fixed_geometry_local_family"]
    readiness = {
        "schema": "oneD.v3.cavity-source-v5-readiness/1",
        "source_repository": "ukaiiaku-maker/PF-fracture-fatigue",
        "source_branch": V5_BRANCH,
        "source_head": head,
        "source_implementation_head": source["executed_code_sha"],
        "source_mesh_contract": source["mesh_contract"],
        "source_recovery_operator": "CAVITY_FIXED_PHYSICAL_ARC_PATCH_RECOVERY_V3",
        "preserved_v1_v2_v3": source["preserved_v1_v2_v3"],
        "preserved_v4": source["preserved_v4"],
        "v4_mesh_failure_cause": attestation["scientific_ledger"]["V4_MESH_FAILURE_CAUSE"],
        "central_dbtt_v5": {
            "DBTT_SOURCE_READINESS": source["DBTT_SOURCE_READINESS"],
            "exact_v5_failure_class": source["exact_v5_failure_class"],
            "required_angular_rows": [_level_summary(row) for row in angular["required_rows"]],
            "required_64_to_128": angular["required_64_to_128"],
            "conditional_256_row": _level_summary(angular["conditional_256_row"]),
            "conditional_128_to_256": angular["conditional_128_to_256"],
            "fixed_polygon_sectors": local["polygon_sectors"],
            "fixed_geometry_local_rows": [_level_summary(row) for row in local["rows"]],
            "fixed_geometry_B_to_C": local["B_to_C"],
        },
        "oracle_states_accepted": 0,
        "oracle_matrix_status": "BLOCKED_CENTRAL_DBTT_V5_RAW_ADJACENT_ELEMENT_TRACTION",
        "paired_trajectories_run": 0,
        "fatigue_started": False,
        "mechanics_map_fitting_complete": False,
        "monotonic_2d_transfer_complete": False,
        "missing_fields_inferred_or_synthesized": False,
        "finite_activation_zone_observable_derived": False,
        "point_source_mesh_development": "STOPPED_AFTER_V5_FAILURE",
        "next_bounded_step": "STOP_POINT_SOURCE_MESH_DEVELOPMENT",
        "bounded_worker": worker,
        "scope": attestation["scope"],
    }
    decision = json.loads((ROOT / DECISION).read_text())
    decision["read_only_v5_head"] = head
    decision["blocking_reason"] = {
        "artifact": "v5_source_readiness.json",
        "code": "CENTRAL_DBTT_V5_RAW_ADJACENT_ELEMENT_TRACTION",
        "exact_v5_failure_class": FAILURE,
        "required_action": "stop further point-source mesh development",
    }
    decision["next_bounded_step"] = readiness["next_bounded_step"]
    decision["execution_counts"].update({
        "bounded_2d_oracle_readiness_cases": 4,
        "v5_bounded_2d_oracle_readiness_cases": 1,
        "v5_required_angular_levels_run": 2,
        "v5_conditional_angular_levels_run": 1,
        "v5_fixed_geometry_local_levels_run": 3,
        "new_2d_oracle_states_accepted": 0,
        "new_2d_trajectories": 0,
        "paired_material_temperature_cells_run": 0,
    })
    decision["pilot_terminal"].update({
        "ONE_D_V3_COMPLETE_ALIGNED_MONOTONIC_TRANSFER": "BLOCKED_ZERO_OF_18_ORACLE_STATES",
        "ONE_D_V3_FIXED_STATE_2D_TRANSFER": "BLOCKED_ZERO_OF_18_ORACLE_STATES",
        "ONE_D_V3_MECHANICS_MAP_FIT": "BLOCKED_CENTRAL_DBTT_V5_RAW_ADJACENT_ELEMENT_TRACTION",
    })
    oracle = decision["prospective_execution_contract"]["oracle_matrix"]
    oracle.update({
        "fixed_arc_source": "CAVITY_FIXED_PHYSICAL_ARC_PATCH_RECOVERY_V3",
        "source_geometry_contract": "CAVITY_SOURCE_SHAPE_REGULAR_LOCAL_PATCH_V5",
        "status": "BLOCKED_CENTRAL_DBTT_V5_RAW_ADJACENT_ELEMENT_TRACTION",
    })
    decision["v5_source_readiness"] = {
        "artifact": "v5_source_readiness.json",
        "DBTT_SOURCE_READINESS": source["DBTT_SOURCE_READINESS"],
        "exact_v5_failure_class": FAILURE,
        "oracle_states_accepted": 0,
        "bounded_worker_run_id": worker["run_id"],
        "bounded_worker_artifact_id": worker["artifact_id"],
        "bounded_worker_artifact_digest": worker["artifact_digest"],
    }
    decision["preserved"]["finite_activation_zone_observable_derived"] = False
    return readiness, decision


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--v5-repo", type=Path, required=True)
    args = parser.parse_args()
    readiness, decision = build(args.v5_repo.resolve())
    (ROOT / OUTPUT).write_text(json.dumps(readiness, indent=2, sort_keys=True) + "\n")
    (ROOT / DECISION).write_text(json.dumps(decision, indent=2, sort_keys=True) + "\n")
    print(readiness["oracle_matrix_status"])
    print("ORACLE_STATES_ACCEPTED=0/18")
    print("PAIRED_TRAJECTORIES_RUN=0/12")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
