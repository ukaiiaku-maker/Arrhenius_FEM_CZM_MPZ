#!/usr/bin/env python3
"""Import the bounded V3 source decision without launching a 2-D solve."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess


ROOT = Path(__file__).resolve().parents[1]
V5_HEAD = "67528be7e38deb7047ad630ac3497bf9926543e1"
V5_BRANCH = "codex/v5-unified-fracture-multitip-voiding"
SOURCE = Path("artifacts/v5_cavity_source_recovery_v3/central_dbtt_v3_readiness.json")
ATTESTATION = Path("artifacts/v5_cavity_source_recovery_v3/attestation.json")
OUTPUT = Path("analysis_outputs/oneD_v3_aligned_temperature_transfer/v3_source_readiness.json")


def _git(repo: Path, *args: str) -> str:
    return subprocess.check_output(("git", *args), cwd=repo, text=True).strip()


def build(v5_repo: Path) -> dict[str, object]:
    head = _git(v5_repo, "rev-parse", "HEAD")
    branch = _git(v5_repo, "branch", "--show-current")
    if head != V5_HEAD or branch != V5_BRANCH:
        raise RuntimeError("2-D source is not at the exact final V3 branch head")
    if _git(v5_repo, "status", "--porcelain=v1"):
        raise RuntimeError("2-D source worktree is not clean")
    source = json.loads((v5_repo / SOURCE).read_text())
    attestation = json.loads((v5_repo / ATTESTATION).read_text())
    if source["DBTT_SOURCE_READINESS"] != "BLOCKED_WITH_EXACT_V3_FAILURE_CLASS":
        raise RuntimeError("unexpected central DBTT V3 decision")
    if source["exact_v3_failure_class"] != ["SOURCE_GEOMETRY_IDENTITY"]:
        raise RuntimeError("unexpected central DBTT V3 failure class")
    if source["oracle_states_accepted"] != 0:
        raise RuntimeError("blocked V3 source cannot provide an oracle row")
    if attestation["workflow"]["conclusion"] != "success":
        raise RuntimeError("bounded V3 worker is not successful")
    return {
        "schema": "oneD.v3.cavity-source-v3-readiness/1",
        "source_repository": "ukaiiaku-maker/PF-fracture-fatigue",
        "source_branch": V5_BRANCH,
        "source_head": head,
        "source_recovery_operator": source["operator"],
        "historical_v2_readiness_artifact": "aligned_oracle_readiness.json",
        "preserved_v1_v2": source["preserved_v1_v2"],
        "reconciled_v2_taxonomy": {
            "decision": "B_STILL_MANDATORY_V2_GATES",
            "failure_classification": [
                "TANGENTIAL_STRESS_CONVERGENCE",
                "NORMAL_DIRECTION_RESOLUTION",
            ],
        },
        "central_dbtt_v3": {
            "DBTT_SOURCE_READINESS": source["DBTT_SOURCE_READINESS"],
            "exact_v3_failure_class": source["exact_v3_failure_class"],
            "scientific_unavailability": source["scientific_unavailability"],
            "geometry": source["geometry"],
            "source_geometry_diagnostic": source["source_geometry_diagnostic"],
            "material_identity": source["material_identity"],
            "accepted_material_identity_complete": source[
                "accepted_material_identity_complete"
            ],
            "ligament_energy_gate": source["ligament_energy_gate"],
            "requested_refinement_levels": source["requested_refinement_levels"],
            "refinement_levels_run": source["refinement_levels_run"],
            "accepted_pre_source_state_unchanged_on_noncertification": source[
                "accepted_pre_source_state_unchanged_on_noncertification"
            ],
            "thresholds_and_rng_unchanged_on_noncertification": source[
                "thresholds_and_rng_unchanged_on_noncertification"
            ],
            "unexpected_programming_exceptions_caught": source[
                "unexpected_programming_exceptions_caught"
            ],
        },
        "oracle_states_accepted": 0,
        "oracle_matrix_status": "BLOCKED_CENTRAL_DBTT_V3_SOURCE_GEOMETRY_IDENTITY",
        "paired_trajectories_run": 0,
        "fatigue_started": False,
        "missing_fields_inferred_or_synthesized": False,
        "next_bounded_step": "DERIVE_FINITE_ACTIVATION_ZONE_WORK_OBSERVABLE",
        "bounded_worker": attestation["workflow"],
        "scope": source["scope"],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--v5-repo", type=Path, required=True)
    args = parser.parse_args()
    payload = build(args.v5_repo.resolve())
    destination = ROOT / OUTPUT
    destination.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(payload["oracle_matrix_status"])
    print("ORACLE_STATES_ACCEPTED=0/18")
    print("PAIRED_TRAJECTORIES_RUN=0/12")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
