#!/usr/bin/env python3
"""Import the bounded V4 source decision without launching a 2-D solve."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess


ROOT = Path(__file__).resolve().parents[1]
V5_HEAD = "9d2cb7cb63fc74f02e5f0b54b5a2d1a394fe7d05"
V5_BRANCH = "codex/v5-unified-fracture-multitip-voiding"
SOURCE = Path("artifacts/v5_cavity_source_recovery_v4/central_dbtt_v4_readiness.json")
ATTESTATION = Path("artifacts/v5_cavity_source_recovery_v4/attestation.json")
OUTPUT = Path("analysis_outputs/oneD_v3_aligned_temperature_transfer/v4_source_readiness.json")
DECISION = Path("analysis_outputs/oneD_v3_aligned_temperature_transfer/decision.json")
FAILURE = ["NORMAL_DIRECTION_RESOLUTION", "TANGENTIAL_DIRECTION_RESOLUTION", "MESH_QUALITY"]


def _git(repo: Path, *args: str) -> str:
    return subprocess.check_output(("git", *args), cwd=repo, text=True).strip()


def build(v5_repo: Path):
    head = _git(v5_repo, "rev-parse", "HEAD")
    branch = _git(v5_repo, "branch", "--show-current")
    if head != V5_HEAD or branch != V5_BRANCH:
        raise RuntimeError("2-D source is not at the exact final V4 branch head")
    if _git(v5_repo, "status", "--porcelain=v1"):
        raise RuntimeError("2-D source worktree is not clean")
    source = json.loads((v5_repo / SOURCE).read_text())
    attestation = json.loads((v5_repo / ATTESTATION).read_text())
    if source["DBTT_SOURCE_READINESS"] != "BLOCKED_WITH_EXACT_V4_FAILURE_CLASS":
        raise RuntimeError("unexpected central DBTT V4 decision")
    if source["exact_v4_failure_class"] != FAILURE or source["oracle_states_accepted"] != 0:
        raise RuntimeError("unexpected central DBTT V4 failure or oracle count")
    worker = attestation["bounded_worker"]
    if worker["conclusion"] != "success" or worker["tests_passed"] != 36:
        raise RuntimeError("bounded V4 worker is not successful")
    readiness = {
        "schema": "oneD.v3.cavity-source-v4-readiness/1",
        "source_repository": "ukaiiaku-maker/PF-fracture-fatigue",
        "source_branch": V5_BRANCH,
        "source_head": head,
        "source_geometry_contract": source["geometry_contract"],
        "source_recovery_operator": "CAVITY_FIXED_PHYSICAL_ARC_PATCH_RECOVERY_V3",
        "preserved_v1_v2_v3": source["preserved_v1_v2_v3"],
        "central_dbtt_v4": {
            "DBTT_SOURCE_READINESS": source["DBTT_SOURCE_READINESS"],
            "exact_v4_failure_class": source["exact_v4_failure_class"],
            "geometry": source["geometry"],
            "levels_requested": source["levels_requested"],
            "levels_run": source["levels_run"],
            "level_records": source["level_records"],
            "accepted_material_identity_exact_across_levels": source[
                "accepted_material_identity_exact_across_levels"
            ],
            "physical_geometry_exact_across_levels": source[
                "physical_geometry_exact_across_levels"
            ],
            "thresholds_and_rng_exact_across_levels": source[
                "thresholds_and_rng_exact_across_levels"
            ],
            "accepted_input_state_unchanged_on_noncertification": source[
                "accepted_input_state_unchanged_on_noncertification"
            ],
            "unexpected_programming_exceptions_caught": source[
                "unexpected_programming_exceptions_caught"
            ],
        },
        "oracle_states_accepted": 0,
        "oracle_matrix_status": "BLOCKED_CENTRAL_DBTT_V4_RESOLUTION_AND_QUALITY",
        "paired_trajectories_run": 0,
        "fatigue_started": False,
        "mechanics_map_fitting_complete": False,
        "monotonic_2d_transfer_complete": False,
        "missing_fields_inferred_or_synthesized": False,
        "finite_activation_zone_observable_derived": False,
        "next_bounded_step": "SEPARATELY_FORMULATE_FINITE_ACTIVATION_ZONE_WORK_OBSERVABLE",
        "bounded_worker": worker,
        "scope": attestation["scope"],
    }
    decision = json.loads((ROOT / DECISION).read_text())
    decision["read_only_v5_head"] = head
    decision["blocking_reason"] = {
        "artifact": "v4_source_readiness.json",
        "code": "CENTRAL_DBTT_V4_RESOLUTION_AND_QUALITY",
        "exact_v4_failure_class": FAILURE,
        "required_action": "separately formulate a finite activation-zone work observable",
    }
    decision["next_bounded_step"] = readiness["next_bounded_step"]
    decision["execution_counts"].update({
        "bounded_2d_oracle_readiness_cases": 3,
        "v4_bounded_2d_oracle_readiness_cases": 1,
        "v4_angular_levels_run": 3,
        "new_2d_oracle_states_accepted": 0,
        "new_2d_trajectories": 0,
        "paired_material_temperature_cells_run": 0,
    })
    decision["pilot_terminal"].update({
        "ONE_D_V3_COMPLETE_ALIGNED_MONOTONIC_TRANSFER": "BLOCKED_ZERO_OF_18_ORACLE_STATES",
        "ONE_D_V3_FIXED_STATE_2D_TRANSFER": "BLOCKED_ZERO_OF_18_ORACLE_STATES",
        "ONE_D_V3_MECHANICS_MAP_FIT": "BLOCKED_CENTRAL_DBTT_V4_RESOLUTION_AND_QUALITY",
    })
    oracle = decision["prospective_execution_contract"]["oracle_matrix"]
    oracle.update({
        "fixed_arc_source": "CAVITY_FIXED_PHYSICAL_ARC_PATCH_RECOVERY_V3",
        "source_geometry_contract": "CAVITY_SOURCE_CONFORMING_GEOMETRY_V4",
        "status": "BLOCKED_CENTRAL_DBTT_V4_RESOLUTION_AND_QUALITY",
    })
    decision["v4_source_readiness"] = {
        "artifact": "v4_source_readiness.json",
        "DBTT_SOURCE_READINESS": source["DBTT_SOURCE_READINESS"],
        "exact_v4_failure_class": FAILURE,
        "oracle_states_accepted": 0,
        "bounded_worker_run_id": worker["run_id"],
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
