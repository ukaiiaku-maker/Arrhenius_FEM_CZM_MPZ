#!/usr/bin/env python3
"""Long-extension and multiseed validation for Taylor/Peierls survivors."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
from typing import Any

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_outputs" / "oneD_v2_taylor_peierls_rcurve_search"
sys.path.insert(0, str(ROOT))

from reduced_fracture_v2.taylor_peierls_contract import SEARCH_WHITELIST
from scripts.run_oneD_v2_taylor_peierls_search import (
    CLASS_SEEDS,
    CONTROL_IDS,
    TEMPERATURES,
    _case_task,
    _worker_initialize,
)


SEEDS = {
    "Peak": (8666, 8667, 8668),
    "DBTT": (1008666, 1008667, 1008668),
}


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _head() -> str:
    return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()


def _selection_features(frame: pd.DataFrame, diagnostics: pd.DataFrame) -> np.ndarray:
    merged = frame.merge(diagnostics, on=["candidate_id", "target_class"], how="left")
    columns = []
    for field in SEARCH_WHITELIST:
        values = merged[field].to_numpy(float)
        if field not in {"peierls_activation_entropy_kB", "taylor_activation_entropy_kB"}:
            values = np.log10(values)
        low, high = np.nanpercentile(values, [1.0, 99.0])
        columns.append(np.clip((values - low) / max(high - low, 1e-300), 0.0, 1.0))
    for field in (
        "log10_peierls_rate_s_median", "log10_taylor_completion_rate_s_median",
        "log10_chi_taylor_completion_median", "retained_equilibrium_fraction_median",
    ):
        values = merged[field].to_numpy(float)
        low, high = np.nanpercentile(values, [1.0, 99.0])
        columns.append(np.clip((values - low) / max(high - low, 1e-300), 0.0, 1.0))
    return np.column_stack(columns)


def _maximin_indices(features: np.ndarray, first: int, count: int) -> list[int]:
    selected = [int(first)]
    distance = np.sum((features - features[first]) ** 2, axis=1)
    distance[first] = -np.inf
    while len(selected) < count:
        index = int(np.argmax(distance))
        selected.append(index)
        distance = np.minimum(distance, np.sum((features - features[index]) ** 2, axis=1))
        distance[np.asarray(selected)] = -np.inf
    return selected


def build_survivors() -> dict[str, Any]:
    selection = pd.read_csv(OUT / "oneD_v2_taylor_peierls_stage1_selection.csv")
    diagnostics = pd.read_csv(OUT / "oneD_v2_taylor_peierls_timescale_diagnostics.csv")
    cases = pd.read_parquet(OUT / "oneD_v2_taylor_peierls_stage1_100um_cases.parquet")
    payload: dict[str, Any] = {
        "schema": "oneD_v2_taylor_peierls_survivor_selection_v1",
        "selection_policy": (
            "MAXIMIN_TAYLOR_PEIERLS_PLUS_SOURCE_KINETIC_FEATURES; CONTROL_FORCED; "
            "NO_100UM_REINIT_REQUIRED_FOR_SLOW_EVOLUTION_GATE"
        ),
        "classes": {},
    }
    for material in ("Peak", "DBTT"):
        local = selection[selection.target_class == material].reset_index(drop=True)
        features = _selection_features(local, diagnostics[diagnostics.target_class == material])
        control_index = int(np.flatnonzero(local.candidate_id.eq(CONTROL_IDS[material]))[0])
        stage2_indices = _maximin_indices(features, control_index, 16)
        stage2 = local.iloc[stage2_indices].candidate_id.astype(str).tolist()
        stage3_indices_local = _maximin_indices(features[stage2_indices], 0, 6)
        stage3 = [stage2[index] for index in stage3_indices_local]
        # The 500-um lane samples the control plus the most distinct kinetic
        # regimes. It is deliberately bounded to five rows per class.
        stage500_indices_local = _maximin_indices(features[stage2_indices], 0, 5)
        stage500 = [stage2[index] for index in stage500_indices_local]
        class_cases = cases[cases.target_class == material]
        payload["classes"][material] = {
            "stage2_300um_candidates": stage2,
            "stage2_500um_candidates": stage500,
            "stage3_multiseed_candidates": stage3,
            "stage1_positive_reinit_case_count": int((class_cases.N_reinit > 0).sum()),
            "stage1_case_count": int(len(class_cases)),
        }
    return payload


def _tasks(selection: pd.DataFrame, survivors: dict[str, Any]) -> list[tuple]:
    lookup = selection.set_index("candidate_id", drop=False)
    tasks: list[tuple] = []
    for material in ("Peak", "DBTT"):
        group = survivors["classes"][material]
        for candidate_id in group["stage2_300um_candidates"]:
            row = lookup.loc[candidate_id].to_dict()
            for provider in ("PF", "FEMCZM"):
                for temperature in TEMPERATURES[material]:
                    tasks.append((row, material, provider, temperature, 300.0,
                                  CLASS_SEEDS[material], "STAGE2_300UM_SLOW_EVOLUTION"))
        for candidate_id in group["stage2_500um_candidates"]:
            row = lookup.loc[candidate_id].to_dict()
            for provider in ("PF", "FEMCZM"):
                for temperature in TEMPERATURES[material]:
                    tasks.append((row, material, provider, temperature, 500.0,
                                  CLASS_SEEDS[material], "STAGE2_500UM_BOUNDED_CONTINUATION"))
        for candidate_id in group["stage3_multiseed_candidates"]:
            row = lookup.loc[candidate_id].to_dict()
            for provider in ("PF", "FEMCZM"):
                for temperature in TEMPERATURES[material]:
                    for seed in SEEDS[material][1:]:
                        tasks.append((row, material, provider, temperature, 300.0, seed,
                                      "STAGE3_300UM_MULTISEED"))
    return tasks


def run_validation(workers: int, checkpoint_every: int) -> None:
    selection = pd.read_csv(OUT / "oneD_v2_taylor_peierls_stage1_selection.csv")
    survivors = build_survivors()
    selection_path = OUT / "oneD_v2_taylor_peierls_survivor_selection.json"
    selection_path.write_text(json.dumps(survivors, indent=2, sort_keys=True) + "\n")
    checkpoint = OUT / "oneD_v2_taylor_peierls_validation_checkpoint.parquet"
    records = pd.read_parquet(checkpoint).to_dict("records") if checkpoint.exists() else []
    done = {
        (str(row["candidate_id"]), str(row["provider"]), float(row["temperature_K"]),
         float(row["target_um"]), int(row["hazard_seed"]))
        for row in records
    }
    tasks = [
        task for task in _tasks(selection, survivors)
        if (str(task[0]["candidate_id"]), task[2], float(task[3]), float(task[4]), int(task[5]))
        not in done
    ]
    print(f"VALIDATION_START pending={len(tasks)} existing={len(records)} workers={workers}", flush=True)
    with ProcessPoolExecutor(max_workers=workers, initializer=_worker_initialize) as executor:
        futures = [executor.submit(_case_task, task) for task in tasks]
        for completed, future in enumerate(as_completed(futures), start=1):
            records.append(future.result())
            if completed % checkpoint_every == 0 or completed == len(futures):
                pd.DataFrame(records).to_parquet(checkpoint, index=False)
                print(f"VALIDATION_PROGRESS completed={completed}/{len(futures)}", flush=True)
    frame = pd.DataFrame(records).sort_values(
        ["target_class", "candidate_id", "target_um", "provider", "temperature_K", "hazard_seed"],
        kind="stable",
    ).reset_index(drop=True)
    validation_path = OUT / "oneD_v2_taylor_peierls_long_validation.parquet"
    frame.to_parquet(validation_path, index=False)
    # Canonical multiseed artifact includes the canonical Stage-2 seed and two
    # independent Stage-3 seeds at 300 um for the six-row diverse subset.
    multiseed_ids = {
        material: set(survivors["classes"][material]["stage3_multiseed_candidates"])
        for material in ("Peak", "DBTT")
    }
    multiseed = frame[
        frame.target_um.eq(300.0)
        & frame.apply(lambda row: row.candidate_id in multiseed_ids[row.target_class], axis=1)
    ].copy()
    multiseed_path = OUT / "oneD_v2_taylor_peierls_multiseed_results.parquet"
    multiseed.to_parquet(multiseed_path, index=False)
    manifest = {
        "schema": "oneD_v2_taylor_peierls_long_multiseed_validation_v1",
        "producer_commit": _head(), "survivor_selection_sha256": _sha(selection_path),
        "case_count": int(len(frame)), "multiseed_case_count": int(len(multiseed)),
        "seeds": SEEDS, "providers": ["PF", "FEMCZM"],
        "long_results_sha256": _sha(validation_path),
        "multiseed_results_sha256": _sha(multiseed_path),
        "heavy_PF_worker_count": 0, "new_FEMCZM_runs": 0,
        "status_counts": frame.status.value_counts(dropna=False).to_dict(),
    }
    (OUT / "oneD_v2_taylor_peierls_validation_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    checkpoint.unlink(missing_ok=True)
    print(f"VALIDATION_COMPLETE cases={len(frame)} multiseed={len(multiseed)}", flush=True)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=min(max(os.cpu_count() or 2, 2), 8))
    parser.add_argument("--checkpoint-every", type=int, default=64)
    args = parser.parse_args()
    run_validation(args.workers, args.checkpoint_every)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
