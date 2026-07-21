#!/usr/bin/env python3
"""Run the v9.12 0-D/1-D emergent-GND DBTT parameterization screen."""
from __future__ import annotations

import argparse
import csv
from dataclasses import replace
import json
from pathlib import Path
from statistics import median
import sys
from typing import Any

from arrhenius_fracture.emergent_gnd_dbtt_v912 import (
    CommonPhysics,
    candidate_from_registry_row,
    developed_delta_K,
    dump_result_json,
    load_protocol_csv,
    run_temperature_protocol,
    score_microstructural_transition,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-registry", required=True)
    parser.add_argument("--protocol-csv", required=True)
    parser.add_argument("--physics-json", required=True)
    parser.add_argument("--temperatures", nargs="+", type=float, required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--stage", choices=("0d", "1d"), default="1d")
    parser.add_argument("--window-um", nargs=2, type=float, default=(10.0, 30.0))
    parser.add_argument("--target-cleavage-rate-s", type=float, default=1.0e-3)
    parser.add_argument("--min-amplitude", type=float, default=8.0)
    parser.add_argument("--target-localization", type=float, default=0.50)
    parser.add_argument("--max-width-K", type=float, default=200.0)
    parser.add_argument(
        "--compact-output",
        action="store_true",
        help=(
            "Write one candidate_summary.json per candidate but omit the "
            "per-temperature JSON files. Recommended for large training runs."
        ),
    )
    parser.add_argument(
        "--quiet-cases",
        action="store_true",
        help="Suppress per-temperature CASE_RESULT terminal lines.",
    )
    return parser.parse_args()


def load_physics(path: str | Path, stage: str) -> CommonPhysics:
    payload = json.loads(Path(path).read_text())
    if "common_physics" in payload:
        payload = payload["common_physics"]
    for key in (
        "emission_signs",
        "emission_schmid_factors",
        "shielding_orientation_factors",
    ):
        if key in payload:
            payload[key] = tuple(payload[key])
    for key in ("forest_interaction_matrix", "gnd_stress_projection_matrix"):
        if key in payload:
            payload[key] = tuple(tuple(row) for row in payload[key])
    physics = CommonPhysics(**payload)
    if stage == "0d":
        physics = replace(
            physics,
            n_bins=1,
            source_zone_length_m=physics.mpz_length_m,
        )
    physics.validate()
    return physics


def load_rows(path: str | Path) -> list[dict[str, Any]]:
    with Path(path).open(newline="") as fp:
        rows = list(csv.DictReader(fp))
    if not rows:
        raise RuntimeError(f"candidate registry is empty: {path}")
    return rows


def values_in_window(
    result: Any,
    field_name: str,
    window_um: tuple[float, float],
) -> list[float]:
    values = list(getattr(result, field_name))
    if len(values) != len(result.extensions_um):
        raise RuntimeError(
            f"{field_name} length differs from extension checkpoint length"
        )
    return [
        float(value)
        for extension, value in zip(result.extensions_um, values)
        if window_um[0] <= float(extension) <= window_um[1]
    ]


def write_ranking(path: Path, records: list[dict[str, Any]]) -> None:
    if not records:
        return
    fields = [
        "rank",
        "candidate_id",
        "stage",
        "score",
        "pass",
        "amplitude_MPa_sqrt_m",
        "largest_jump_localization",
        "transition_width_10_90_K",
        "linear_r2",
        "max_abs_K_shield_MPa_sqrt_m",
        "max_tau_gnd_tip_MPa",
        "max_gnd_abs_line_count_per_unit_thickness",
        "min_source_available_fraction",
        "min_source_available_fraction_pre_advance",
        "median_abs_K_shield_developed_window_MPa_sqrt_m",
        "median_abs_tau_gnd_tip_developed_window_MPa",
        "median_gnd_abs_line_count_developed_window_per_unit_thickness",
        "min_source_available_fraction_pre_advance_developed_window",
    ]
    ordered = sorted(records, key=lambda row: float(row["score"]), reverse=True)
    with path.open("w", newline="") as fp:
        writer = csv.DictWriter(fp, fieldnames=fields)
        writer.writeheader()
        for rank, row in enumerate(ordered, start=1):
            writer.writerow({"rank": rank, **{key: row[key] for key in fields[1:]}})


def main() -> int:
    args = parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    physics = load_physics(args.physics_json, args.stage)
    protocol = load_protocol_csv(args.protocol_csv)
    rows = load_rows(args.candidate_registry)
    window_um = tuple(float(value) for value in args.window_um)

    print(
        "CAMPAIGN_START "
        f"stage={args.stage} candidates={len(rows)} "
        f"temperatures={','.join(str(int(T)) for T in args.temperatures)}",
        flush=True,
    )
    ranking: list[dict[str, Any]] = []
    for index, row in enumerate(rows, start=1):
        candidate = candidate_from_registry_row(row)
        candidate_root = out / candidate.candidate_id
        candidate_root.mkdir(parents=True, exist_ok=True)
        print(
            f"CANDIDATE_START index={index}/{len(rows)} "
            f"candidate={candidate.candidate_id}",
            flush=True,
        )
        results = []
        developed = []
        for T in args.temperatures:
            result = run_temperature_protocol(
                candidate,
                physics,
                protocol,
                T,
                target_cleavage_rate_s=args.target_cleavage_rate_s,
            )
            results.append(result)
            developed_value = developed_delta_K(result, window_um)
            developed.append(developed_value)
            if not args.compact_output:
                dump_result_json(
                    candidate_root / f"T{int(round(T))}K.json",
                    result.as_dict(),
                )
            if not args.quiet_cases:
                print(
                    "CASE_RESULT "
                    f"candidate={candidate.candidate_id} T={T:g} "
                    f"deltaKmicro={developed_value:.8g} "
                    f"Kshield={result.K_shield_MPa_sqrt_m[-1]:.8g} "
                    f"tauGND_MPa={result.tau_gnd_tip_MPa[-1]:.8g}",
                    flush=True,
                )

        score = score_microstructural_transition(
            args.temperatures,
            developed,
            min_amplitude=args.min_amplitude,
            target_localization=args.target_localization,
            max_width_K=args.max_width_K,
        )
        pre_advance_source_values = [
            value
            for result in results
            for value in result.source_available_fraction_pre_advance
        ]
        developed_shield = [
            abs(value)
            for result in results
            for value in values_in_window(
                result, "K_shield_MPa_sqrt_m", window_um
            )
        ]
        developed_tau = [
            abs(value)
            for result in results
            for value in values_in_window(result, "tau_gnd_tip_MPa", window_um)
        ]
        developed_gnd = [
            value
            for result in results
            for value in values_in_window(
                result,
                "gnd_abs_line_count_per_unit_thickness",
                window_um,
            )
        ]
        developed_source_pre = [
            value
            for result in results
            for value in values_in_window(
                result,
                "source_available_fraction_pre_advance",
                window_um,
            )
        ]
        if not (
            developed_shield
            and developed_tau
            and developed_gnd
            and developed_source_pre
        ):
            raise RuntimeError(
                f"candidate {candidate.candidate_id} has no diagnostics in "
                f"developed window {window_um}"
            )

        record = {
            "candidate_id": candidate.candidate_id,
            "stage": args.stage,
            **score,
            "max_abs_K_shield_MPa_sqrt_m": float(
                max(abs(v) for result in results for v in result.K_shield_MPa_sqrt_m)
            ),
            "max_tau_gnd_tip_MPa": float(
                max(v for result in results for v in result.tau_gnd_tip_MPa)
            ),
            "max_gnd_abs_line_count_per_unit_thickness": float(
                max(
                    v
                    for result in results
                    for v in result.gnd_abs_line_count_per_unit_thickness
                )
            ),
            "min_source_available_fraction": float(
                min(v for result in results for v in result.source_available_fraction)
            ),
            "min_source_available_fraction_pre_advance": float(
                min(pre_advance_source_values)
                if pre_advance_source_values
                else min(
                    v for result in results for v in result.source_available_fraction
                )
            ),
            "median_abs_K_shield_developed_window_MPa_sqrt_m": float(
                median(developed_shield)
            ),
            "median_abs_tau_gnd_tip_developed_window_MPa": float(
                median(developed_tau)
            ),
            "median_gnd_abs_line_count_developed_window_per_unit_thickness": float(
                median(developed_gnd)
            ),
            "min_source_available_fraction_pre_advance_developed_window": float(
                min(developed_source_pre)
            ),
        }
        signed_state_active = (
            record["max_abs_K_shield_MPa_sqrt_m"] > 1.0e-6
            and record["max_gnd_abs_line_count_per_unit_thickness"] > 0.0
        )
        spatial_gnd_active = record["max_tau_gnd_tip_MPa"] > 1.0e-6
        mechanism_active = signed_state_active and (
            args.stage == "0d" or spatial_gnd_active
        )
        record["signed_state_gate_pass"] = signed_state_active
        record["spatial_gnd_gate_required"] = args.stage == "1d"
        record["spatial_gnd_gate_pass"] = spatial_gnd_active
        record["mechanism_active"] = mechanism_active
        record["pass"] = bool(record["pass"] and mechanism_active)
        ranking.append(record)
        dump_result_json(
            candidate_root / "candidate_summary.json",
            {
                "candidate_id": candidate.candidate_id,
                "stage": args.stage,
                "temperatures_K": list(args.temperatures),
                "developed_window_um": list(window_um),
                "developed_delta_K_micro_MPa_sqrt_m": developed,
                "objective": record,
                "K0_target_or_penalty_active": False,
                "explicit_N_sat_active": False,
                "independent_backstress_law_active": False,
                "compact_output": bool(args.compact_output),
            },
        )
        print(
            "CANDIDATE_COMPLETE "
            f"candidate={candidate.candidate_id} score={record['score']:.8g} "
            f"pass={int(record['pass'])}",
            flush=True,
        )

    write_ranking(out / "ranking.csv", ranking)
    dump_result_json(
        out / "campaign_summary.json",
        {
            "stage": args.stage,
            "candidate_registry": str(Path(args.candidate_registry).resolve()),
            "protocol_csv": str(Path(args.protocol_csv).resolve()),
            "physics_json": str(Path(args.physics_json).resolve()),
            "temperatures_K": list(args.temperatures),
            "developed_window_um": list(window_um),
            "K0_target_or_penalty_active": False,
            "compact_output": bool(args.compact_output),
            "ranked_candidates": sorted(
                ranking, key=lambda row: float(row["score"]), reverse=True
            ),
        },
    )
    print(
        f"CAMPAIGN_COMPLETE candidates={len(ranking)} out={out}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
