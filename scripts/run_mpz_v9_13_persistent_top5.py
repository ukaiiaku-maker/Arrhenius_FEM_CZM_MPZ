#!/usr/bin/env python3
"""Rerun the selected v9.12 candidates with the v9.13 persistent-site law."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys
import time
from typing import Any

import numpy as np

from arrhenius_fracture.emergent_gnd_dbtt_v913 import (
    CommonPhysics,
    candidate_from_registry_row,
    developed_delta_K,
    dump_result_json,
    load_protocol_csv,
    run_temperature_protocol,
    score_microstructural_transition,
)


OLD_REFERENCE = {
    "v912_targeted_local_peak_013476_0368": [
        0, 0, 0, 0, 0, 0, 15.89420517793896, 9.835585678115162,
        1.590275820632705, 0.8850894135232181,
    ],
    "v912_targeted_local_peak_013476_0314": [
        0, 0, 0, 0, 0, 5.08364767682389, 15.146514773430699,
        9.858534253225447, 2.226740819019838, 0.7105000696549979,
    ],
    "v912_targeted_local_peak_013476_0162": [
        0, 0, 0, 0, 0, 3.8031528226309277, 10.403619730387117,
        7.769957241582379, 5.2306921637390715, 1.7864366875737119,
    ],
    "v912_targeted_local_peak_005518_0118": [
        0, 0, 0, 0, 0, 0, 0.29613604492855217, 4.7521916442766745,
        5.352200533471027, 8.132296921811808,
    ],
    "v912_targeted_local_plateau_010759_0403": [
        0, 0, 0, 0, 0, 0, 13.063468568822184, 14.645020510375357,
        16.13723968585314, 17.809147095201837,
    ],
}
OLD_TEMPERATURES = np.arange(300.0, 1201.0, 100.0)

PERSISTENT_FIELDS = (
    "persistent_site_multiplicity_per_system",
    "persistent_site_source_area_m2",
    "persistent_site_front_width_m",
    "persistent_site_width_density_m2",
    "persistent_tip_radius_m",
    "persistent_rho_back_mean_m2",
    "persistent_sigma_back_mean_Pa",
    "persistent_backstress_drive_ratio_max",
    "persistent_last_source_activations",
    "persistent_last_line_content",
    "persistent_local_accumulated_slip_count",
    "tip_resharpening_by_advance_m",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-registry", required=True)
    parser.add_argument("--protocol-csv", required=True)
    parser.add_argument("--physics-json", required=True)
    parser.add_argument("--temperatures", nargs="+", type=float, required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--window-um", nargs=2, type=float, default=(10.0, 30.0))
    parser.add_argument("--target-cleavage-rate-s", type=float, default=1.0e-3)
    parser.add_argument("--min-amplitude", type=float, default=8.0)
    parser.add_argument("--target-localization", type=float, default=0.50)
    parser.add_argument("--max-width-K", type=float, default=200.0)
    parser.add_argument("--compact-output", action="store_true")
    return parser.parse_args()


def load_physics(path: Path) -> tuple[CommonPhysics, dict[str, Any]]:
    payload = json.loads(path.read_text())
    metadata = dict(payload)
    common = dict(payload.get("common_physics", payload))
    for key in (
        "emission_signs",
        "emission_schmid_factors",
        "shielding_orientation_factors",
        "activation_to_line_content_per_system",
    ):
        if key in common:
            common[key] = tuple(common[key])
    for key in ("forest_interaction_matrix", "gnd_stress_projection_matrix"):
        if key in common:
            common[key] = tuple(tuple(row) for row in common[key])
    physics = CommonPhysics(**common)
    physics.validate()
    return physics, metadata


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise RuntimeError(f"empty candidate registry: {path}")
    return rows


def window_values(result: Any, field: str, window: tuple[float, float]) -> list[float]:
    values = list(getattr(result, field))
    return [
        float(value)
        for extension, value in zip(result.extensions_um, values)
        if window[0] <= float(extension) <= window[1]
    ]


def response_class(temperatures: np.ndarray, values: np.ndarray) -> str:
    peak_index = int(np.argmax(values))
    peak = float(values[peak_index])
    peak_temperature = float(temperatures[peak_index])
    final_ratio = float(values[-1] / peak) if peak > 1.0e-12 else float("nan")
    low = (
        float(np.median(values[temperatures <= 700.0]))
        if np.any(temperatures <= 700.0)
        else float(values[0])
    )
    high = (
        float(np.median(values[temperatures >= 1000.0]))
        if np.any(temperatures >= 1000.0)
        else float(values[-1])
    )
    if peak < 0.01:
        return "negligible_shielding"
    if 800.0 <= peak_temperature <= 1000.0 and final_ratio < 0.8:
        return "transition_peak"
    if peak_temperature >= 1000.0 and high - low > 0.1 and final_ratio >= 0.8:
        return "late_high_temperature_response"
    if peak_temperature <= 500.0:
        return "low_temperature_inverse"
    return "other"


def old_values_for(candidate_id: str, temperatures: np.ndarray) -> np.ndarray:
    old = np.asarray(OLD_REFERENCE.get(candidate_id, []), dtype=float)
    if old.size != OLD_TEMPERATURES.size:
        return np.full(temperatures.shape, np.nan)
    return np.interp(temperatures, OLD_TEMPERATURES, old)


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    out = Path(args.out)
    if (out / "persistent_top5_summary.json").exists():
        raise RuntimeError(f"completed output already exists: {out}")
    out.mkdir(parents=True, exist_ok=True)
    physics, physics_payload = load_physics(Path(args.physics_json))
    rows = load_rows(Path(args.candidate_registry))
    protocol = load_protocol_csv(args.protocol_csv)
    temperatures = np.asarray(args.temperatures, dtype=float)
    window = (float(args.window_um[0]), float(args.window_um[1]))

    conversion = list(physics.activation_to_line_content_per_system)
    print(
        "V913_CAMPAIGN_START "
        f"candidates={len(rows)} temperatures={len(temperatures)} "
        f"line_conversion={conversion} "
        f"provenance={physics_payload.get('conversion_provenance', 'unspecified')}",
        flush=True,
    )

    temperature_rows: list[dict[str, Any]] = []
    candidate_rows: list[dict[str, Any]] = []
    failures: list[dict[str, str]] = []
    start = time.perf_counter()

    for index, row in enumerate(rows, start=1):
        candidate = candidate_from_registry_row(row)
        candidate_root = out / candidate.candidate_id
        candidate_root.mkdir(parents=True, exist_ok=True)
        print(
            f"V913_CANDIDATE_START index={index}/{len(rows)} candidate={candidate.candidate_id}",
            flush=True,
        )
        try:
            results = []
            developed = []
            for temperature in temperatures:
                result = run_temperature_protocol(
                    candidate,
                    physics,
                    protocol,
                    float(temperature),
                    target_cleavage_rate_s=args.target_cleavage_rate_s,
                )
                results.append(result)
                value = developed_delta_K(result, window)
                developed.append(value)
                if not args.compact_output:
                    dump_result_json(
                        candidate_root / f"T{int(round(temperature))}K.json",
                        result.as_dict(),
                    )

                temp_record: dict[str, Any] = {
                    "candidate_id": candidate.candidate_id,
                    "temperature_K": float(temperature),
                    "developed_delta_K_micro_MPa_sqrt_m": float(value),
                    "old_reference_delta_K_micro_MPa_sqrt_m": float(
                        old_values_for(candidate.candidate_id, np.asarray([temperature]))[0]
                    ),
                }
                for field in PERSISTENT_FIELDS:
                    vals = window_values(result, field, window)
                    temp_record[f"developed_median_{field}"] = (
                        float(np.median(vals)) if vals else float("nan")
                    )
                temperature_rows.append(temp_record)
                print(
                    "V913_CASE_RESULT "
                    f"candidate={candidate.candidate_id} T={temperature:g} "
                    f"deltaKmicro={value:.8g} "
                    f"r_tip_um={result.persistent_tip_radius_m[-1] * 1e6:.8g} "
                    f"multiplicity={result.persistent_site_multiplicity_per_system[-1]:.8g}",
                    flush=True,
                )

            values = np.asarray(developed, dtype=float)
            old_values = old_values_for(candidate.candidate_id, temperatures)
            peak_index = int(np.argmax(values))
            old_peak_index = (
                int(np.nanargmax(old_values))
                if np.any(np.isfinite(old_values))
                else 0
            )
            score = score_microstructural_transition(
                temperatures,
                values,
                min_amplitude=args.min_amplitude,
                target_localization=args.target_localization,
                max_width_K=args.max_width_K,
            )
            candidate_record: dict[str, Any] = {
                "candidate_id": candidate.candidate_id,
                "status": "complete",
                "response_class": response_class(temperatures, values),
                "peak_delta_K_micro_MPa_sqrt_m": float(values[peak_index]),
                "peak_temperature_K": float(temperatures[peak_index]),
                "final_delta_K_micro_MPa_sqrt_m": float(values[-1]),
                "old_reference_peak_delta_K_micro_MPa_sqrt_m": float(
                    old_values[old_peak_index]
                ),
                "old_reference_peak_temperature_K": float(
                    temperatures[old_peak_index]
                ),
                "peak_change_MPa_sqrt_m": float(
                    values[peak_index] - old_values[old_peak_index]
                ),
                "maximum_tip_radius_um": float(
                    max(
                        v
                        for result in results
                        for v in result.persistent_tip_radius_m
                    )
                    * 1e6
                ),
                "minimum_front_width_um": float(
                    min(
                        v
                        for result in results
                        for v in result.persistent_site_front_width_m
                    )
                    * 1e6
                ),
                "minimum_multiplicity_per_system": float(
                    min(
                        v
                        for result in results
                        for v in result.persistent_site_multiplicity_per_system
                    )
                ),
                "maximum_multiplicity_per_system": float(
                    max(
                        v
                        for result in results
                        for v in result.persistent_site_multiplicity_per_system
                    )
                ),
                "maximum_backstress_MPa": float(
                    max(
                        v
                        for result in results
                        for v in result.persistent_sigma_back_mean_Pa
                    )
                    / 1e6
                ),
                "maximum_resharpening_um": float(
                    max(
                        v
                        for result in results
                        for v in result.tip_resharpening_by_advance_m
                    )
                    * 1e6
                ),
                "finite_source_inventory": False,
                "source_refresh_active": False,
                "explicit_recovery_active": False,
                **score,
            }
            candidate_rows.append(candidate_record)
            dump_result_json(
                candidate_root / "candidate_summary.json",
                {
                    "candidate_id": candidate.candidate_id,
                    "model": "v9.13_persistent_site_backstress_blunting",
                    "temperatures_K": temperatures.tolist(),
                    "developed_window_um": list(window),
                    "developed_delta_K_micro_MPa_sqrt_m": values.tolist(),
                    "old_reference_delta_K_micro_MPa_sqrt_m": old_values.tolist(),
                    "objective": candidate_record,
                    "physics": physics_payload,
                },
            )
            print(
                "V913_CANDIDATE_COMPLETE "
                f"candidate={candidate.candidate_id} peak={values[peak_index]:.8g} "
                f"Tpeak={temperatures[peak_index]:g} "
                f"class={candidate_record['response_class']}",
                flush=True,
            )
        except Exception as exc:
            failures.append(
                {"candidate_id": candidate.candidate_id, "error": repr(exc)}
            )
            candidate_rows.append(
                {
                    "candidate_id": candidate.candidate_id,
                    "status": "unresolved",
                    "error": repr(exc),
                }
            )
            (candidate_root / "failure.json").write_text(
                json.dumps(failures[-1], indent=2) + "\n"
            )
            print(
                f"V913_CANDIDATE_FAILED candidate={candidate.candidate_id} error={exc!r}",
                flush=True,
            )

    candidate_rows.sort(
        key=lambda row: float(row.get("peak_delta_K_micro_MPa_sqrt_m", -1.0)),
        reverse=True,
    )
    for rank, row in enumerate(candidate_rows, start=1):
        row["rank"] = rank
    write_csv(out / "persistent_top5_temperature.csv", temperature_rows)
    write_csv(out / "persistent_top5_candidate.csv", candidate_rows)
    write_csv(out / "ranking.csv", candidate_rows)
    summary = {
        "model": "v9.13_persistent_site_backstress_blunting",
        "candidate_count": len(rows),
        "complete_count": len(rows) - len(failures),
        "unresolved_count": len(failures),
        "elapsed_s": time.perf_counter() - start,
        "conversion_provenance": physics_payload.get("conversion_provenance"),
        "activation_to_line_content_per_system": conversion,
        "finite_source_inventory": False,
        "source_refresh_active": False,
        "explicit_recovery_active": False,
        "records": candidate_rows,
        "failures": failures,
    }
    (out / "persistent_top5_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    print(
        "V913_CAMPAIGN_COMPLETE "
        f"complete={summary['complete_count']} "
        f"unresolved={summary['unresolved_count']} out={out}",
        flush=True,
    )
    return 0 if not failures else 1


if __name__ == "__main__":
    sys.exit(main())
