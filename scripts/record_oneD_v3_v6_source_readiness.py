#!/usr/bin/env python3
"""Import the bounded V6 source decision without launching a 2-D solve."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess


ROOT = Path(__file__).resolve().parents[1]
V6_HEAD = "c6c8fd500cf506258d70e56a8f391d52908f03d4"
V6_BRANCH = "codex/v5-unified-fracture-multitip-voiding"
SOURCE = Path("artifacts/v6_cavity_source_raw_traction_closure/central_dbtt_v6_readiness.json")
ATTESTATION = Path("artifacts/v6_cavity_source_raw_traction_closure/attestation.json")
OUTPUT = Path("analysis_outputs/oneD_v3_aligned_temperature_transfer/v6_source_readiness.json")
DECISION = Path("analysis_outputs/oneD_v3_aligned_temperature_transfer/decision.json")
FAILURE = [
    "RAW_ADJACENT_ELEMENT_TRACTION_E", "MESH_QUALITY", "FIXED_GEOMETRY_IDENTITY",
]


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
        "source_node_certificate", "sigma_zz_and_mean_stress_audit", "observables",
        "physical_equilibrium_observable_predicates", "ligament_energy_gate",
        "candidate_kinetic_drive", "candidate_kinetic_drive_status",
    )
    return {key: row.get(key) for key in keys}


def build(v6_repo: Path):
    head = _git(v6_repo, "rev-parse", "HEAD")
    branch = _git(v6_repo, "branch", "--show-current")
    if head != V6_HEAD or branch != V6_BRANCH:
        raise RuntimeError("2-D source is not at the exact final V6 branch head")
    if _git(v6_repo, "status", "--porcelain=v1"):
        raise RuntimeError("2-D source worktree is not clean")
    source = json.loads((v6_repo / SOURCE).read_text())
    attestation = json.loads((v6_repo / ATTESTATION).read_text())
    if source["DBTT_SOURCE_READINESS"] != "BLOCKED_WITH_EXACT_V6_FAILURE_CLASS":
        raise RuntimeError("unexpected central DBTT V6 decision")
    if source["exact_v6_failure_class"] != FAILURE or source["oracle_states_accepted"] != 0:
        raise RuntimeError("unexpected central DBTT V6 failure or oracle count")
    worker = attestation["bounded_worker"]
    if (worker["conclusion"], worker["job_count"], worker["tests_passed"]) != (
        "success", 1, 34,
    ):
        raise RuntimeError("bounded V6 worker is not successful")
    family = source["fixed_geometry_local_family"]
    readiness = {
        "schema": "oneD.v3.cavity-source-v6-readiness/1",
        "source_repository": "ukaiiaku-maker/PF-fracture-fatigue",
        "source_branch": V6_BRANCH,
        "source_head": head,
        "source_implementation_head": source["executed_code_sha"],
        "source_workflow_head": worker["workflow_head"],
        "source_contract": source["contract"],
        "source_recovery_operator": "CAVITY_FIXED_PHYSICAL_ARC_PATCH_RECOVERY_V3",
        "preserved_v5": source["preserved_v5"],
        "central_dbtt_v6": {
            "DBTT_SOURCE_READINESS": source["DBTT_SOURCE_READINESS"],
            "exact_v6_failure_class": source["exact_v6_failure_class"],
            "retained_A_B_C_raw_traction": family["retained_A_B_C_raw_traction"],
            "new_D_E_rows": [_level_summary(row) for row in family["new_rows"]],
            "D_to_E_comparisons": family["D_to_E_comparisons"],
            "predicates": family["predicates"],
        },
        "oracle_states_accepted": 0,
        "oracle_matrix_status": "BLOCKED_CENTRAL_DBTT_V6_RAW_TRACTION_QUALITY_AND_GEOMETRY_IDENTITY",
        "paired_trajectories_run": 0,
        "fatigue_started": False,
        "mechanics_map_fitting_complete": False,
        "monotonic_2d_transfer_complete": False,
        "missing_fields_inferred_or_synthesized": False,
        "finite_activation_zone_observable_derived": False,
        "point_source_mesh_development": "CLOSED_AFTER_FINAL_V6_FAILURE",
        "next_bounded_step": "EQUILIBRATED_BOUNDARY_STRESS_RECONSTRUCTION",
        "bounded_worker": worker,
        "cancelled_scope_control_runs": attestation["cancelled_scope_control_runs"],
    }
    decision = json.loads((ROOT / DECISION).read_text())
    decision["read_only_v6_head"] = head
    decision["blocking_reason"] = {
        "artifact": "v6_source_readiness.json",
        "code": "CENTRAL_DBTT_V6_RAW_TRACTION_QUALITY_AND_GEOMETRY_IDENTITY",
        "exact_v6_failure_class": FAILURE,
        "required_action": "open equilibrated boundary-stress reconstruction as a separate numerical-method task",
    }
    decision["next_bounded_step"] = readiness["next_bounded_step"]
    decision["execution_counts"].update({
        "bounded_2d_oracle_readiness_cases": 5,
        "v6_bounded_2d_oracle_readiness_cases": 1,
        "v6_fixed_geometry_local_levels_run": 2,
        "new_2d_oracle_states_accepted": 0,
        "new_2d_trajectories": 0,
        "paired_material_temperature_cells_run": 0,
    })
    decision["pilot_terminal"].update({
        "ONE_D_V3_COMPLETE_ALIGNED_MONOTONIC_TRANSFER": "BLOCKED_ZERO_OF_18_ORACLE_STATES",
        "ONE_D_V3_FIXED_STATE_2D_TRANSFER": "BLOCKED_ZERO_OF_18_ORACLE_STATES",
        "ONE_D_V3_MECHANICS_MAP_FIT": (
            "BLOCKED_CENTRAL_DBTT_V6_RAW_TRACTION_QUALITY_AND_GEOMETRY_IDENTITY"
        ),
    })
    decision["prospective_execution_contract"]["oracle_matrix"].update({
        "fixed_arc_source": "CAVITY_FIXED_PHYSICAL_ARC_PATCH_RECOVERY_V3",
        "source_geometry_contract": "CAVITY_SOURCE_RAW_TRACTION_CLOSURE_V6",
        "status": "BLOCKED_CENTRAL_DBTT_V6_RAW_TRACTION_QUALITY_AND_GEOMETRY_IDENTITY",
    })
    decision["v6_source_readiness"] = {
        "artifact": "v6_source_readiness.json",
        "DBTT_SOURCE_READINESS": source["DBTT_SOURCE_READINESS"],
        "exact_v6_failure_class": FAILURE,
        "oracle_states_accepted": 0,
        "bounded_worker_run_id": worker["run_id"],
        "bounded_worker_artifact_id": worker["artifact_id"],
        "bounded_worker_artifact_digest": worker["artifact_digest"],
    }
    decision["preserved"]["finite_activation_zone_observable_derived"] = False
    return readiness, decision


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--v6-repo", type=Path, required=True)
    args = parser.parse_args()
    readiness, decision = build(args.v6_repo.resolve())
    (ROOT / OUTPUT).write_text(json.dumps(readiness, indent=2, sort_keys=True) + "\n")
    (ROOT / DECISION).write_text(json.dumps(decision, indent=2, sort_keys=True) + "\n")
    print(readiness["oracle_matrix_status"])
    print("ORACLE_STATES_ACCEPTED=0/18")
    print("PAIRED_TRAJECTORIES_RUN=0/12")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
