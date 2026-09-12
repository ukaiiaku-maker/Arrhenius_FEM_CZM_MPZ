#!/usr/bin/env python3
"""Run the bounded material-transfer preflight without launching 2-D solves."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from reduced_fracture_v3.materials import (  # noqa: E402
    FRACTURE_ROWS,
    load_exact_fracture_rows,
    material_field_mapping_audit,
    paired_case_ledger,
    pilot_material_bundles,
)


PINNED_V5_RUNTIME_SHA = "b58997bdb18cf4e9a32c251c073115d8b405bb27"
PINNED_V5_ATTESTATION_SHA = "c7583ecd0a259f28ce92780833d21a358a920f45"
FROZEN_TENSOR_PA = ((3.0e9, 0.0), (0.0, 3.0e9))
TEMPERATURES_K = tuple(range(300, 1201, 25))


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _git(repo: Path, *args: str) -> str:
    return subprocess.check_output(("git", *args), cwd=repo, text=True).strip()


def _verify_v5(repo: Path) -> str:
    head = _git(repo, "rev-parse", "HEAD")
    for required in (PINNED_V5_RUNTIME_SHA, PINNED_V5_ATTESTATION_SHA):
        result = subprocess.run(
            ("git", "merge-base", "--is-ancestor", required, head), cwd=repo
        )
        if result.returncode:
            raise RuntimeError(f"pinned V5 identity {required} is not an ancestor of {head}")
    if _git(repo, "status", "--porcelain=v1"):
        raise RuntimeError("read-only V5 oracle worktree is not clean")
    if subprocess.run(
        ("git", "diff", "--quiet", PINNED_V5_ATTESTATION_SHA, head, "--", "arrhenius_fracture"),
        cwd=repo,
    ).returncode:
        raise RuntimeError("V5 source tree differs from the pinned attestation")
    return head


def _void_rate_scout(v5_repo: Path, v5_head: str) -> dict[str, object]:
    request = {
        "temperatures_K": TEMPERATURES_K,
        "stress_tensor_Pa": FROZEN_TENSOR_PA,
    }
    code = r'''
import dataclasses
import json
import sys
import numpy as np
from arrhenius_fracture.voiding_v5 import VoidingConfig, arrhenius_rates

request = json.loads(sys.stdin.read())
config = VoidingConfig(enabled=True)
tensor = np.asarray(request["stress_tensor_Pa"], dtype=float)
rows = []
for temperature in request["temperatures_K"]:
    rates = arrhenius_rates(config, temperature_K=float(temperature), stress_tensor_Pa=tensor)
    rows.append({"temperature_K": temperature, **rates})
print(json.dumps({"void_kinetics_row": dataclasses.asdict(config), "rows": rows}, sort_keys=True))
'''
    env = {"PYTHONPATH": str(v5_repo)}
    completed = subprocess.run(
        (sys.executable, "-c", code),
        cwd=v5_repo,
        env=env,
        input=json.dumps(request),
        text=True,
        capture_output=True,
        check=True,
    )
    result = json.loads(completed.stdout)
    return {
        "schema": "oneD.v3.common-void-rate-scout/1",
        "source_repository": "ukaiiaku-maker/PF-fracture-fatigue",
        "source_head": v5_head,
        "runtime_hardening_sha": PINNED_V5_RUNTIME_SHA,
        "attestation_sha": PINNED_V5_ATTESTATION_SHA,
        "source_function": "arrhenius_fracture.voiding_v5.arrhenius_rates",
        "frozen_representative_tensor_Pa": FROZEN_TENSOR_PA,
        "temperature_grid_K": TEMPERATURES_K,
        "void_kinetics_row": result["void_kinetics_row"],
        "rows": result["rows"],
        "material_dependence": "NONE_COMMON_REFERENCE_VOID_KINETICS_ROW",
        "feature_anchor_use": (
            "INSUFFICIENT_FOR_MATERIAL_SPECIFIC_FEATURE_ANCHORS; "
            "requires an exactly bound fracture row"
        ),
    }


def _write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError("cannot write empty transfer ledger")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--registry",
        type=Path,
        default=ROOT / "analysis_outputs/oneD_v2_terminal_predictive_program/oneD_v2_new_four_class_registry.csv",
    )
    parser.add_argument("--v5-repo", type=Path, required=True)
    parser.add_argument(
        "--out",
        type=Path,
        default=ROOT / "analysis_outputs/oneD_v3_aligned_temperature_transfer",
    )
    args = parser.parse_args()

    v5_head = _verify_v5(args.v5_repo.resolve())
    rows = load_exact_fracture_rows(args.registry)
    bundles = pilot_material_bundles()
    audit = material_field_mapping_audit(
        rows,
        source_registry="analysis_outputs/oneD_v2_terminal_predictive_program/oneD_v2_new_four_class_registry.csv",
        v5_driver_identity=f"{v5_head}:arrhenius_fracture.voiding_production_v5",
    )
    scout = _void_rate_scout(args.v5_repo.resolve(), v5_head)
    anchors = {family: (300.0, None, 1200.0) for family in FRACTURE_ROWS}
    ledger = paired_case_ledger(anchors, gate=str(audit["summary"]["gate"]))

    bundle_payload = {
        "schema": "oneD.v3.material-bundle-registry/1",
        "fracture_rows_source_sha256": _sha256(args.registry),
        "common_reference_void_row_for_all_families": True,
        "separate_void_barrier_fit_per_family": False,
        "bundles": {family: bundle.__dict__ for family, bundle in bundles.items()},
    }
    decision = {
        "schema": "oneD.v3.aligned-temperature-transfer-decision/1",
        "mission": "PAIRED 2-D / 1-D ALIGNED-VOID TEMPERATURE TRANSFER PILOT",
        "amendment": "FOUR_RESPONSE_FAMILIES_WITH_COMMON_REFERENCE_VOID_KINETICS",
        "one_d_checkpoint_sha": _git(ROOT, "rev-parse", "HEAD"),
        "read_only_v5_head": v5_head,
        "read_only_v5_source_tree_matches_attestation": True,
        "material_mapping_gate": audit["summary"]["gate"],
        "temperature_scout": "PASS_COMMON_VOID_RATE_SCOUT_300_TO_1200K",
        "temperature_anchors": {
            family: {
                "T_low_K": 300.0,
                "T_feature_K": None,
                "T_feature_status": "NOT_FROZEN_M2_MATERIAL_BINDING_GATE_BLOCKED",
                "T_high_K": 1200.0,
            }
            for family in FRACTURE_ROWS
        },
        "execution_counts": {
            "new_2d_mechanics_solves": 0,
            "new_2d_trajectories": 0,
            "new_1d_trajectories": 0,
            "paired_material_temperature_cells_run": 0,
            "planned_material_temperature_cells": 12,
        },
        "material_transfer": {
            "MATERIAL_ROW_TRANSFER_PEAK": "BLOCKED_UNMAPPED_ACTIVE_FIELD",
            "MATERIAL_ROW_TRANSFER_DBTT": "BLOCKED_UNMAPPED_ACTIVE_FIELD",
            "MATERIAL_ROW_TRANSFER_WEAK_T": "BLOCKED_UNMAPPED_ACTIVE_FIELD",
            "MATERIAL_ROW_TRANSFER_CERAMIC_LIKE": "BLOCKED_UNMAPPED_ACTIVE_FIELD",
        },
        "temperature_transfer": {
            "ONE_D_V3_PEAK_TEMPERATURE_TRANSFER": "NOT_RUN_M2_GATE_BLOCKED",
            "ONE_D_V3_DBTT_TEMPERATURE_TRANSFER": "NOT_RUN_M2_GATE_BLOCKED",
            "ONE_D_V3_WEAK_T_TEMPERATURE_TRANSFER": "NOT_RUN_M2_GATE_BLOCKED",
            "ONE_D_V3_CERAMIC_LIKE_TEMPERATURE_TRANSFER": "NOT_RUN_M2_GATE_BLOCKED",
        },
        "pilot_terminal": {
            "ONE_D_V3_MECHANICS_MAP_FIT": "NOT_RUN_M2_GATE_BLOCKED",
            "ONE_D_V3_FIXED_STATE_2D_TRANSFER": "NOT_RUN_M2_GATE_BLOCKED",
            "ONE_D_V3_LOW_T_TRAJECTORY_TRANSFER": "NOT_RUN_M2_GATE_BLOCKED",
            "ONE_D_V3_TRANSITION_T_TRAJECTORY_TRANSFER": "NOT_RUN_M2_GATE_BLOCKED",
            "ONE_D_V3_HIGH_T_TRAJECTORY_TRANSFER": "NOT_RUN_M2_GATE_BLOCKED",
            "ONE_D_V3_COMPLETE_ALIGNED_MONOTONIC_TRANSFER": "BLOCKED",
            "FATIGUE_IMPLEMENTATION": "NOT_STARTED_BY_CONTRACT",
        },
        "held_out_policy": {
            "development_sentinels": ["Peak", "DBTT", "weak-T"],
            "complete_held_out_material_family": "ceramic-like",
            "ceramic_like_used_for_tuning": False,
        },
        "prospective_execution_contract": {
            "load_mapping": {
                "G0_2D": "candidate kinetic drive of matched no-void 2-D companion",
                "K0_1D": "sqrt(Eprime * max(G0_2D, 0))",
                "raw_2d_opening_used_as_1d_K": False,
            },
            "oracle_matrix": {
                "planned_unique_states": 18,
                "topology_regimes": [
                    "PRECONNECTION",
                    "CONNECTED_CAVITY",
                    "DOWNSTREAM_CHILD",
                ],
                "geometry_levels": ["low", "central", "high"],
                "load_levels_per_geometry": 2,
                "fixed_arc_source": "actual V5 fixed-arc recovery operator",
                "status": "NOT_RUN_M2_GATE_BLOCKED",
            },
            "fit_criteria": {
                "held_out_G_relative_error_max": 0.05,
                "held_out_tensor_relative_error_max": 0.05,
                "compliance_relative_error_max": 0.02,
                "energy_relative_error_max": 0.02,
                "correct_sign_and_event_eligibility": True,
                "exact_no_void_far_void_small_void_limits": True,
                "fail_closed_out_of_domain": True,
                "load_scaling_test_before_factorization": True,
            },
            "baseline_comparison": {
                "raw_observables": True,
                "void_increment": "Delta O_void = O_void - O_no_void",
                "attribute_no_void_model_difference_to_void_reduction": False,
            },
            "paired_state_exact_equality": [
                "thresholds",
                "RNG state",
                "site state",
                "initial void state",
                "fracture state",
                "r_tip/emission/shielding state",
                "load history",
                "temperature history",
                "stop condition",
            ],
            "complete_pass_rule": (
                "all four material families pass fixed-state and trajectory comparisons"
            ),
        },
        "preserved": {
            "r_tip_equals_R_void": False,
            "fractured_length_equals_free_span_equals_front_coordinate": False,
            "missing_mechanics_fields_inferred": False,
            "fracture_rows_changed": False,
            "separate_void_barriers_fit": False,
            "fatigue_started": False,
        },
        "blocking_reason": (
            "The pinned read-only V5 one-void runtime creates default FrontEngine and "
            "material state internally and exposes no exact fracture-material-row binding. "
            "Running the matrix would substitute defaults and violate M2."
        ),
        "next_bounded_step": (
            "Add and separately qualify an explicit material-bundle injection contract in "
            "the V5 runtime, then rerun this preflight before any oracle or trajectory solve."
        ),
    }

    _write_json(args.out / "material_bundles.json", bundle_payload)
    _write_json(args.out / "material_field_mapping_audit.json", audit)
    _write_json(args.out / "common_void_rate_scout.json", scout)
    _write_csv(args.out / "paired_case_ledger.csv", ledger)
    _write_json(args.out / "decision.json", decision)
    print(audit["summary"]["gate"])
    print("PAIRED_CASES_RUN=0/12")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
