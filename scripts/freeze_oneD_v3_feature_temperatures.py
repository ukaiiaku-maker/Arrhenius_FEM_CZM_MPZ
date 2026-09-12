#!/usr/bin/env python3
"""Freeze material feature temperatures before any paired trajectory exists."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from reduced_fracture_v3.materials import FRACTURE_ROWS, paired_case_ledger  # noqa: E402
from scripts.audit_oneD_v3_material_transfer import (  # noqa: E402
    FROZEN_TENSOR_PA,
    UNIFIED_V5_BRANCH,
    UNIFIED_V5_SHA,
    _git,
    _verify_v5,
    _write_csv,
    _write_json,
)

TEMPERATURES_K = tuple(range(300, 1201, 25))
FEATURE_GRID_K = tuple(range(450, 1101, 25))
WEAK_T_GRID_K = tuple(range(450, 1051, 25))
CERAMIC_LIKE_FEATURE_K = 1100


def _fracture_rate_rows(v5_repo: Path) -> list[dict[str, object]]:
    request = {
        "temperatures_K": TEMPERATURES_K,
        "stress_Pa": float(FROZEN_TENSOR_PA[0][0]),
    }
    code = r'''
import json
import math
import sys
from arrhenius_fracture.unified_fracture_material_v5 import all_material_bundles, material_manifest

request = json.loads(sys.stdin.read())
rows = []
for family, bundle in all_material_bundles().items():
    manifest = material_manifest(bundle)
    for temperature in request["temperatures_K"]:
        cleavage = float(manifest.cleavage.rate(request["stress_Pa"], temperature))
        emission = float(manifest.emission.rate(request["stress_Pa"], temperature))
        rows.append({
            "material_class": family,
            "fracture_material_row_id": bundle.fracture_material_row_id,
            "temperature_K": temperature,
            "cleavage_rate_s": cleavage,
            "emission_rate_s": emission,
            "log10_cleavage_rate": math.log10(max(cleavage, 1.0e-300)),
            "log10_emission_rate": math.log10(max(emission, 1.0e-300)),
            "log10_cleavage_to_emission": math.log10(max(cleavage, 1.0e-300)) - math.log10(max(emission, 1.0e-300)),
        })
print(json.dumps(rows, sort_keys=True))
'''
    completed = subprocess.run(
        (sys.executable, "-c", code), cwd=v5_repo,
        env={"PYTHONPATH": str(v5_repo)}, input=json.dumps(request),
        text=True, capture_output=True, check=True,
    )
    return json.loads(completed.stdout)


def _select(rows: list[dict[str, object]]) -> tuple[dict[str, int], dict[str, object]]:
    by_family = {
        family: [row for row in rows if row["material_class"] == family]
        for family in FRACTURE_ROWS
    }
    by_family_temperature = {
        family: {int(row["temperature_K"]): row for row in family_rows}
        for family, family_rows in by_family.items()
    }

    peak_candidates = [by_family_temperature["Peak"][T] for T in FEATURE_GRID_K]
    peak = max(
        peak_candidates,
        key=lambda row: 0.5 * (
            float(row["log10_cleavage_rate"]) + float(row["log10_emission_rate"])
        ),
    )
    dbtt_candidates = [by_family_temperature["DBTT"][T] for T in FEATURE_GRID_K]
    dbtt = min(
        dbtt_candidates,
        key=lambda row: abs(float(row["log10_cleavage_to_emission"])),
    )

    development = ("Peak", "DBTT", "weak-T")
    vectors = np.asarray([
        [
            float(by_family_temperature[family][T]["log10_cleavage_rate"]),
            float(by_family_temperature[family][T]["log10_emission_rate"]),
        ]
        for family in development for T in WEAK_T_GRID_K
    ])
    mean = vectors.mean(axis=0)
    scale = vectors.std(axis=0)
    weak_scores: list[tuple[float, int, list[float]]] = []
    for temperature in WEAK_T_GRID_K:
        weak = np.asarray([
            float(by_family_temperature["weak-T"][temperature]["log10_cleavage_rate"]),
            float(by_family_temperature["weak-T"][temperature]["log10_emission_rate"]),
        ])
        distances = []
        for comparator in ("Peak", "DBTT"):
            other = np.asarray([
                float(by_family_temperature[comparator][temperature]["log10_cleavage_rate"]),
                float(by_family_temperature[comparator][temperature]["log10_emission_rate"]),
            ])
            distances.append(float(np.linalg.norm((weak - other) / scale)))
        weak_scores.append((min(distances), temperature, distances))
    weak_score, weak_temperature, weak_distances = max(weak_scores)

    anchors = {
        "Peak": int(peak["temperature_K"]),
        "DBTT": int(dbtt["temperature_K"]),
        "weak-T": int(weak_temperature),
        # This value is declared independently of the held-out row results.
        "ceramic-like": CERAMIC_LIKE_FEATURE_K,
    }
    evidence = {
        "Peak": {
            "rule": "maximize geometric mean of cleavage and emission rates on 450..1100 K interior grid",
            "score_log10_rate_s": 0.5 * (
                float(peak["log10_cleavage_rate"]) + float(peak["log10_emission_rate"])
            ),
        },
        "DBTT": {
            "rule": "minimize absolute log10 cleavage/emission rate gap on 450..1100 K interior grid",
            "absolute_log10_rate_gap": abs(float(dbtt["log10_cleavage_to_emission"])),
        },
        "weak-T": {
            "rule": "maximize minimum standardized log-rate-vector distance from Peak and DBTT on 450..1050 K midrange grid",
            "minimum_standardized_distance": weak_score,
            "distances_to_peak_and_dbtt": weak_distances,
        },
        "ceramic-like": {
            "rule": "prospectively declared high-temperature interior sentinel",
            "held_out_rate_results_used_for_selection": False,
        },
    }
    return anchors, evidence


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--v5-repo", type=Path, required=True)
    parser.add_argument(
        "--out", type=Path,
        default=ROOT / "analysis_outputs/oneD_v3_aligned_temperature_transfer",
    )
    args = parser.parse_args()
    v5_repo = args.v5_repo.resolve()
    head = _verify_v5(v5_repo)
    if head != UNIFIED_V5_SHA or _git(v5_repo, "branch", "--show-current") != UNIFIED_V5_BRANCH:
        raise RuntimeError("unified 2-D source identity changed")

    decision_path = args.out / "decision.json"
    decision = json.loads(decision_path.read_text())
    if decision["material_mapping_gate"] != "PASS_EXACT_RUNTIME_BINDING":
        raise RuntimeError("feature temperatures cannot be frozen before all M2 gates pass")
    if any(
        int(value)
        for key, value in decision["execution_counts"].items()
        if key != "planned_material_temperature_cells"
    ):
        raise RuntimeError("paired or oracle execution already exists; prospective freeze refused")

    rows = _fracture_rate_rows(v5_repo)
    anchors, evidence = _select(rows)
    scout = {
        "schema": "oneD.v3.fracture-rate-feature-scout/1",
        "source_repository": "ukaiiaku-maker/PF-fracture-fatigue",
        "source_head": head,
        "source_function": "MaterialManifest.cleavage/emission.rate",
        "frozen_representative_stress_Pa": float(FROZEN_TENSOR_PA[0][0]),
        "temperature_grid_K": list(TEMPERATURES_K),
        "feature_grid_K": list(FEATURE_GRID_K),
        "weak_T_discrimination_grid_K": list(WEAK_T_GRID_K),
        "selection_executed_before_paired_results": True,
        "paired_material_temperature_cells_present_at_selection": 0,
        "held_out_material_family": "ceramic-like",
        "held_out_results_used_for_temperature_selection": False,
        "selected_feature_temperatures_K": anchors,
        "selection_evidence": evidence,
        "rows": rows,
    }
    _write_json(args.out / "fracture_rate_feature_scout.json", scout)

    for family, temperature in anchors.items():
        decision["temperature_anchors"][family]["T_feature_K"] = float(temperature)
        decision["temperature_anchors"][family]["T_feature_status"] = "FROZEN_PROSPECTIVELY_BEFORE_PAIRED_RESULTS"
    decision["feature_temperature_gate"] = "PASS_FROZEN_PROSPECTIVE_ANCHORS"
    decision["prospective_execution_contract"]["oracle_matrix"]["status"] = "NOT_RUN_PENDING_BOUNDED_ALIGNED_ORACLE"
    decision["next_bounded_step"] = "Generate the bounded 18-state aligned 2-D mechanics oracle without inferring absent fields."
    _write_json(decision_path, decision)
    ledger = paired_case_ledger(
        {family: (300.0, float(temperature), 1200.0) for family, temperature in anchors.items()},
        gate="PASS_EXACT_RUNTIME_BINDING",
    )
    _write_csv(args.out / "paired_case_ledger.csv", ledger)
    print("FEATURE_TEMPERATURE_GATE=PASS_FROZEN_PROSPECTIVE_ANCHORS")
    print(json.dumps(anchors, sort_keys=True))
    print("PAIRED_CASES_RUN=0/12")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
