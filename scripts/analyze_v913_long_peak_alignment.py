#!/usr/bin/env python3
"""Analyze long-extension R-curves and create a common-temperature registry."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import numpy as np
import pandas as pd

from arrhenius_fracture.dbtt_long_alignment_v913 import (
    achieved_extension_um,
    checkpoint_from_events,
    checkpoint_reached,
    choose_peak_temperature,
    peak_drift_classification,
    peak_metrics,
)
from arrhenius_fracture.dbtt_transform_v913 import (
    active_parameter_row,
    temperature_scale_candidate_row,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case-root", type=Path, required=True)
    parser.add_argument("--candidate-registry", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument(
        "--checkpoints-um",
        type=float,
        nargs="+",
        default=(25.0, 50.0, 75.0, 100.0),
    )
    parser.add_argument("--target-peak-temperature-K", type=float, default=900.0)
    parser.add_argument(
        "--peak-estimator",
        choices=("discrete", "quadratic"),
        default="discrete",
    )
    parser.add_argument("--stable-drift-limit-K", type=float, default=50.0)
    parser.add_argument("--maximum-alignable-drift-K", type=float, default=100.0)
    parser.add_argument("--minimum-post-peak-drop-MPa-sqrt-m", type=float, default=1.0)
    parser.add_argument("--refinement-step-K", type=float, default=25.0)
    parser.add_argument("--refinement-half-width-K", type=float, default=100.0)
    return parser.parse_args()


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise RuntimeError(f"empty candidate registry: {path}")
    return rows


def load_payloads(case_root: Path) -> list[dict[str, Any]]:
    paths = sorted(case_root.glob("*/T*K.json"))
    if not paths:
        raise RuntimeError(f"no case JSON files under {case_root}")
    payloads: list[dict[str, Any]] = []
    for path in paths:
        payload = json.loads(path.read_text())
        payload["_source_path"] = str(path)
        payloads.append(payload)
    return payloads


def tag(value: float) -> str:
    rounded = round(float(value))
    if np.isclose(value, rounded):
        return str(int(rounded))
    return f"{float(value):g}".replace(".", "p")


def case_rows(
    payloads: list[dict[str, Any]],
    checkpoints_um: tuple[float, ...],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for payload in payloads:
        events = list(payload.get("events", []))
        row: dict[str, Any] = {
            "candidate_id": str(payload["candidate_id"]),
            "temperature_K": float(payload["temperature_K"]),
            "status": str(payload["status"]),
            "seed": int(payload["seed"]),
            "achieved_projected_extension_um": achieved_extension_um(events),
            "n_events": len(events),
            "max_backstress_GPa": float(payload.get("max_backstress_GPa", float("nan"))),
            "min_front_width_um": float(payload.get("min_front_width_um", float("nan"))),
            "max_tip_radius_um": float(payload.get("max_tip_radius_um", float("nan"))),
            "max_source_multiplicity": float(payload.get("max_source_multiplicity", float("nan"))),
            "source_case_json": payload.get("_source_path", ""),
        }
        for checkpoint in checkpoints_um:
            suffix = tag(checkpoint)
            reached = checkpoint_reached(events, checkpoint)
            row[f"reached_{suffix}um"] = reached
            row[f"K_{suffix}um_MPa_sqrt_m"] = checkpoint_from_events(
                events,
                checkpoint,
                strict=True,
            )
        rows.append(row)
    return rows


def refinement_plan(
    candidate_id: str,
    peak_temperature_K: float,
    existing_temperatures: list[float],
    *,
    step_K: float,
    half_width_K: float,
) -> list[dict[str, Any]]:
    if not np.isfinite(peak_temperature_K):
        return []
    low = peak_temperature_K - half_width_K
    high = peak_temperature_K + half_width_K
    first = math.ceil(low / step_K) * step_K
    requested: list[dict[str, Any]] = []
    existing = np.asarray(existing_temperatures, dtype=float)
    temperature = first
    while temperature <= high + 1.0e-9:
        if not np.any(np.isclose(existing, temperature, rtol=0.0, atol=1.0e-9)):
            requested.append(
                {
                    "candidate_id": candidate_id,
                    "temperature_K": float(temperature),
                    "reason": "long_extension_peak_refinement",
                    "coarse_peak_temperature_K": peak_temperature_K,
                }
            )
        temperature += step_K
    return requested


def main() -> int:
    args = parse_args()
    checkpoints = tuple(sorted({float(value) for value in args.checkpoints_um}))
    if not checkpoints or checkpoints[0] <= 0.0:
        raise ValueError("checkpoints must be positive")
    if args.target_peak_temperature_K <= 0.0:
        raise ValueError("target peak temperature must be positive")
    if args.stable_drift_limit_K < 0.0:
        raise ValueError("stable drift limit must be nonnegative")
    if args.maximum_alignable_drift_K < args.stable_drift_limit_K:
        raise ValueError("maximum alignable drift must exceed stable drift limit")
    if args.refinement_step_K <= 0.0 or args.refinement_half_width_K < 0.0:
        raise ValueError("refinement controls are invalid")

    args.out.mkdir(parents=True, exist_ok=True)
    registry_rows = read_csv_rows(args.candidate_registry)
    registry_by_id = {str(row["candidate_id"]): row for row in registry_rows}
    if len(registry_by_id) != len(registry_rows):
        raise RuntimeError("candidate IDs are not unique in the registry")

    payloads = load_payloads(args.case_root)
    unknown = sorted(
        {str(payload["candidate_id"]) for payload in payloads} - set(registry_by_id)
    )
    if unknown:
        raise RuntimeError(f"case payloads contain unknown candidates: {unknown}")
    cases = case_rows(payloads, checkpoints)
    case_table = pd.DataFrame(cases).sort_values(["candidate_id", "temperature_K"])
    case_table.to_csv(args.out / "long_case_checkpoint_table.csv", index=False)

    checkpoint_metrics: list[dict[str, Any]] = []
    candidate_metrics: list[dict[str, Any]] = []
    refinement_rows: list[dict[str, Any]] = []
    aligned_registry_rows: list[dict[str, Any]] = []
    rejected_rows: list[dict[str, Any]] = []
    aligned_case_rows: list[dict[str, Any]] = []

    for candidate_id in registry_by_id:
        local = case_table[case_table["candidate_id"] == candidate_id].copy()
        if local.empty:
            rejected_rows.append(
                {"candidate_id": candidate_id, "alignment_rejection_reason": "no_cases"}
            )
            continue
        local = local.sort_values("temperature_K")
        temperatures = local["temperature_K"].to_numpy(dtype=float)
        peak_temperatures: list[float] = []
        local_checkpoint_metrics: dict[float, dict[str, Any]] = {}

        for checkpoint in checkpoints:
            suffix = tag(checkpoint)
            values = local[f"K_{suffix}um_MPa_sqrt_m"].to_numpy(dtype=float)
            metrics = peak_metrics(temperatures, values)
            chosen_peak = choose_peak_temperature(metrics, args.peak_estimator)
            peak_temperatures.append(chosen_peak)
            record = {
                "candidate_id": candidate_id,
                "checkpoint_um": checkpoint,
                "peak_temperature_used_K": chosen_peak,
                **metrics.as_dict(),
            }
            checkpoint_metrics.append(record)
            local_checkpoint_metrics[checkpoint] = record

        drift, drift_class = peak_drift_classification(
            peak_temperatures,
            stable_limit_K=args.stable_drift_limit_K,
            maximum_alignable_drift_K=args.maximum_alignable_drift_K,
        )
        final_checkpoint = checkpoints[-1]
        final_metrics = local_checkpoint_metrics[final_checkpoint]
        final_peak = float(final_metrics["peak_temperature_used_K"])
        boundary = bool(final_metrics["peak_at_boundary"])
        drop = float(final_metrics["post_peak_drop"])
        all_reached = bool(
            local[[f"reached_{tag(value)}um" for value in checkpoints]]
            .to_numpy(dtype=bool)
            .all()
        )
        all_complete = bool((local["status"].astype(str) == "complete").all())

        reasons: list[str] = []
        if not all_complete:
            reasons.append("incomplete_case_grid")
        if not all_reached:
            reasons.append("checkpoint_not_reached")
        if boundary:
            reasons.append("long_peak_at_temperature_boundary")
        if not np.isfinite(final_peak):
            reasons.append("long_peak_unresolved")
        if not np.isfinite(drop) or drop < args.minimum_post_peak_drop_MPa_sqrt_m:
            reasons.append("insufficient_post_peak_drop")
        if drift_class == "extension_dependent":
            reasons.append("peak_temperature_drifts_with_extension")

        alignable = not reasons
        scale = (
            float(args.target_peak_temperature_K) / final_peak
            if alignable
            else float("nan")
        )
        summary: dict[str, Any] = {
            "candidate_id": candidate_id,
            "all_cases_complete": all_complete,
            "all_checkpoints_reached": all_reached,
            "peak_estimator": args.peak_estimator,
            "long_checkpoint_um": final_checkpoint,
            "long_peak_temperature_K": final_peak,
            "long_peak_value_MPa_sqrt_m": final_metrics["peak_value"],
            "long_peak_rise_MPa_sqrt_m": final_metrics["peak_rise"],
            "long_post_peak_drop_MPa_sqrt_m": drop,
            "long_final_rebound_MPa_sqrt_m": final_metrics["final_rebound"],
            "peak_temperature_drift_K": drift,
            "peak_drift_classification": drift_class,
            "target_peak_temperature_K": args.target_peak_temperature_K,
            "temperature_scale": scale,
            "alignable": alignable,
            "alignment_rejection_reason": ";".join(reasons),
        }
        for checkpoint in checkpoints:
            suffix = tag(checkpoint)
            metrics = local_checkpoint_metrics[checkpoint]
            summary[f"peak_temperature_{suffix}um_K"] = metrics[
                "peak_temperature_used_K"
            ]
            summary[f"peak_value_{suffix}um_MPa_sqrt_m"] = metrics["peak_value"]
            summary[f"post_peak_drop_{suffix}um_MPa_sqrt_m"] = metrics[
                "post_peak_drop"
            ]
            summary[f"peak_boundary_{suffix}um"] = metrics["peak_at_boundary"]
        candidate_metrics.append(summary)

        refinement_rows.extend(
            refinement_plan(
                candidate_id,
                final_peak,
                temperatures.tolist(),
                step_K=args.refinement_step_K,
                half_width_K=args.refinement_half_width_K,
            )
        )

        if not alignable:
            rejected_rows.append(summary)
            continue

        transformed = temperature_scale_candidate_row(
            registry_by_id[candidate_id],
            scale,
            candidate_id=(
                f"{candidate_id}__aligned_Tp{tag(args.target_peak_temperature_K)}K"
            ),
        )
        registry_record = active_parameter_row(transformed)
        registry_record.update(
            {
                "source_candidate_id": candidate_id,
                "source_long_peak_temperature_K": final_peak,
                "target_peak_temperature_K": args.target_peak_temperature_K,
                "peak_temperature_drift_K": drift,
                "peak_drift_classification": drift_class,
                "alignment_source_checkpoint_um": final_checkpoint,
            }
        )
        aligned_registry_rows.append(registry_record)

        for _, case in local.iterrows():
            record = case.to_dict()
            record["source_candidate_id"] = candidate_id
            record["candidate_id"] = registry_record["candidate_id"]
            record["source_temperature_K"] = float(case["temperature_K"])
            record["temperature_K"] = scale * float(case["temperature_K"])
            record["temperature_scale"] = scale
            record["derived_by_exact_temperature_axis_scale"] = True
            aligned_case_rows.append(record)

    metrics_table = pd.DataFrame(candidate_metrics)
    if not metrics_table.empty:
        metrics_table = metrics_table.sort_values(
            ["alignable", "long_post_peak_drop_MPa_sqrt_m", "candidate_id"],
            ascending=[False, False, True],
        )
    metrics_table.to_csv(args.out / "candidate_long_peak_metrics.csv", index=False)
    pd.DataFrame(checkpoint_metrics).to_csv(
        args.out / "checkpoint_peak_metrics.csv",
        index=False,
    )
    pd.DataFrame(aligned_registry_rows).to_csv(
        args.out / f"aligned_registry_Tp{tag(args.target_peak_temperature_K)}K.csv",
        index=False,
    )
    pd.DataFrame(rejected_rows).to_csv(
        args.out / "alignment_rejections.csv",
        index=False,
    )
    pd.DataFrame(aligned_case_rows).to_csv(
        args.out / "aligned_case_results_derived.csv",
        index=False,
    )
    if refinement_rows:
        refinement_table = pd.DataFrame(refinement_rows).drop_duplicates(
            subset=["candidate_id", "temperature_K"]
        ).sort_values(["candidate_id", "temperature_K"])
    else:
        refinement_table = pd.DataFrame(
            columns=[
                "candidate_id",
                "temperature_K",
                "reason",
                "coarse_peak_temperature_K",
            ]
        )
    refinement_table.to_csv(args.out / "refinement_case_plan.csv", index=False)

    manifest = {
        "schema": "v9.13_long_extension_peak_alignment_v1",
        "candidate_registry": str(args.candidate_registry.resolve()),
        "candidate_registry_sha256": sha256_path(args.candidate_registry),
        "case_root": str(args.case_root.resolve()),
        "candidate_count": len(registry_rows),
        "case_count": len(payloads),
        "checkpoints_um": list(checkpoints),
        "target_peak_temperature_K": args.target_peak_temperature_K,
        "peak_estimator": args.peak_estimator,
        "stable_drift_limit_K": args.stable_drift_limit_K,
        "maximum_alignable_drift_K": args.maximum_alignable_drift_K,
        "minimum_post_peak_drop_MPa_sqrt_m": (
            args.minimum_post_peak_drop_MPa_sqrt_m
        ),
        "alignable_count": len(aligned_registry_rows),
        "rejected_count": len(rejected_rows),
        "aligned_curves_are_derived_not_rerun": True,
    }
    (args.out / "alignment_manifest.json").write_text(
        json.dumps(manifest, indent=2, allow_nan=False) + "\n"
    )
    print(
        "V913_LONG_PEAK_ALIGNMENT_COMPLETE "
        f"candidates={len(registry_rows)} alignable={len(aligned_registry_rows)} "
        f"rejected={len(rejected_rows)} target_T={args.target_peak_temperature_K:g} "
        f"out={args.out}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
