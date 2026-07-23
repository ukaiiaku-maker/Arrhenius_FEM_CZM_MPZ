#!/usr/bin/env python3
"""Run exact temperature translations and constrained low-shelf diagnostics."""
from __future__ import annotations

import argparse
import csv
import itertools
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import pandas as pd

from arrhenius_fracture.dbtt_transform_v913 import (
    active_parameter_row,
    anchored_cleavage_pivot_row,
    scale_cleavage_stress_axis_row,
    temperature_scale_candidate_row,
    validate_positive_barrier_domain,
)
from arrhenius_fracture.emergent_gnd_campaign_v913 import candidate_from_registry_row
from arrhenius_fracture.emergent_gnd_contract_v913 import (
    ACTIVE_CANDIDATE_PARAMETER_FIELDS,
)
from arrhenius_fracture.emergent_gnd_rcurve_v913 import (
    RCurveLoadingMap,
    RCurveResult,
    run_autonomous_rcurve,
)
from scripts.run_mpz_v9_13_persistent_top5 import load_physics


IDENTITY_FIELDS = (
    "threshold_action",
    "applied_displacement_m",
    "elapsed_time_s",
    "K_MPa_sqrt_m",
    "path_advance_m",
    "projected_advance_m",
    "cumulative_path_extension_m",
    "cumulative_projected_extension_m",
    "tip_radius_pre_advance_m",
    "tip_radius_post_advance_m",
    "front_width_pre_advance_m",
    "backstress_pre_advance_Pa",
    "source_multiplicity_pre_advance",
    "cumulative_source_activations",
    "cumulative_line_content",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-registry", type=Path, required=True)
    parser.add_argument(
        "--candidate-id",
        default="v912_targeted_local_peak_013476_0083",
    )
    parser.add_argument(
        "--base-physics-json",
        type=Path,
        default=Path("mpz_v9_13_v10222_transfer_common_physics.json"),
    )
    parser.add_argument("--loading-map", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--current-peak-temperature-K", type=float, default=1000.0)
    parser.add_argument(
        "--target-peak-temperatures-K",
        type=float,
        nargs="+",
        default=(700.0, 1100.0),
    )
    parser.add_argument(
        "--identity-base-temperatures-K",
        type=float,
        nargs="+",
        default=(700.0, 800.0, 900.0, 1000.0, 1100.0, 1200.0),
    )
    parser.add_argument("--shelf-temperature-K", type=float, default=700.0)
    parser.add_argument("--anchor-temperature-K", type=float, default=1000.0)
    parser.add_argument(
        "--evaluation-temperatures-K",
        type=float,
        nargs="+",
        default=(
            300.0,
            400.0,
            500.0,
            600.0,
            700.0,
            800.0,
            900.0,
            1000.0,
            1100.0,
            1200.0,
        ),
    )
    parser.add_argument(
        "--shelf-energy-factors",
        type=float,
        nargs="+",
        default=(1.0, 0.95, 0.9),
    )
    parser.add_argument(
        "--shelf-stress-factors",
        type=float,
        nargs="+",
        default=(1.0, 0.8, 0.6),
    )
    parser.add_argument(
        "--global-cleavage-stress-scales",
        type=float,
        nargs="*",
        default=(0.8, 0.6, 0.4),
    )
    parser.add_argument("--minimum-zero-stress-energy-eV", type=float, default=0.05)
    parser.add_argument("--minimum-characteristic-stress-GPa", type=float, default=0.1)
    parser.add_argument("--target-extension-um", type=float, default=50.0)
    parser.add_argument("--translation-action-exponent", type=float, default=0.95)
    parser.add_argument("--max-hazard-increment", type=float, default=0.05)
    return parser.parse_args()


def read_candidate(path: Path, candidate_id: str) -> dict[str, Any]:
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    matches: list[dict[str, Any]] = []
    for source in rows:
        row: dict[str, Any] = dict(source)
        for field in ACTIVE_CANDIDATE_PARAMETER_FIELDS:
            raw_feature = row.get(f"x_raw__{field}")
            if row.get(field) in (None, "") and raw_feature not in (None, ""):
                row[field] = raw_feature
        if str(row.get("candidate_id")) == candidate_id:
            matches.append(row)
    if len(matches) != 1:
        raise RuntimeError(
            f"expected one row for {candidate_id!r} in {path}, found {len(matches)}"
        )
    return matches[0]


def run_curve(
    row: Mapping[str, Any],
    physics: Any,
    loading_map: RCurveLoadingMap,
    temperature_K: float,
    args: argparse.Namespace,
) -> RCurveResult:
    return run_autonomous_rcurve(
        candidate_from_registry_row(row),
        physics,
        loading_map,
        float(temperature_K),
        target_projected_extension_m=float(args.target_extension_um) * 1.0e-6,
        translation_action_exponent=float(args.translation_action_exponent),
        max_hazard_increment=float(args.max_hazard_increment),
    )


def max_event_differences(
    baseline: RCurveResult,
    transformed: RCurveResult,
) -> dict[str, float]:
    if len(baseline.events) != len(transformed.events):
        raise RuntimeError(
            "temperature scaling changed the event count: "
            f"{len(baseline.events)} versus {len(transformed.events)}"
        )
    output: dict[str, float] = {}
    for field in IDENTITY_FIELDS:
        values = (
            abs(float(getattr(left, field)) - float(getattr(right, field)))
            for left, right in zip(
                baseline.events,
                transformed.events,
                strict=True,
            )
        )
        output[f"max_abs_diff__{field}"] = max(values, default=0.0)
    return output


def event_rows(
    result: RCurveResult,
    variant_id: str,
    transform_kind: str,
) -> list[dict[str, Any]]:
    return [
        {
            "variant_id": variant_id,
            "transform_kind": transform_kind,
            "candidate_id": result.candidate_id,
            "temperature_K": result.temperature_K,
            "status": result.status,
            "seed": result.seed,
            **event.as_dict(),
        }
        for event in result.events
    ]


def create_shelf_variants(
    base_row: Mapping[str, Any],
    args: argparse.Namespace,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    variants = [
        {
            **dict(base_row),
            "candidate_id": f"{base_row['candidate_id']}__baseline",
            "variant_id": "baseline",
            "transform_kind": "baseline",
        }
    ]
    rejected: list[dict[str, Any]] = []
    temperatures = tuple(float(value) for value in args.evaluation_temperatures_K)

    for energy_factor, stress_factor in itertools.product(
        args.shelf_energy_factors,
        args.shelf_stress_factors,
    ):
        energy_factor = float(energy_factor)
        stress_factor = float(stress_factor)
        if np.isclose(energy_factor, 1.0) and np.isclose(stress_factor, 1.0):
            continue
        variant_id = f"pivot_G{energy_factor:.6g}_S{stress_factor:.6g}"
        row = anchored_cleavage_pivot_row(
            base_row,
            shelf_temperature_K=float(args.shelf_temperature_K),
            anchor_temperature_K=float(args.anchor_temperature_K),
            shelf_energy_factor=energy_factor,
            shelf_stress_factor=stress_factor,
            candidate_id=f"{base_row['candidate_id']}__{variant_id}",
        )
        row["variant_id"] = variant_id
        try:
            validate_positive_barrier_domain(
                row,
                temperatures,
                minimum_zero_stress_energy_eV=float(
                    args.minimum_zero_stress_energy_eV
                ),
                minimum_characteristic_stress_GPa=float(
                    args.minimum_characteristic_stress_GPa
                ),
            )
        except ValueError as exc:
            rejected.append(
                {
                    "variant_id": variant_id,
                    "transform_kind": row["transform_kind"],
                    "reason": str(exc),
                }
            )
        else:
            variants.append(row)

    for scale in args.global_cleavage_stress_scales:
        factor = float(scale)
        variant_id = f"global_cleavage_stress_{factor:.6g}"
        row = scale_cleavage_stress_axis_row(
            base_row,
            factor,
            candidate_id=f"{base_row['candidate_id']}__{variant_id}",
        )
        row["variant_id"] = variant_id
        validate_positive_barrier_domain(
            row,
            temperatures,
            minimum_zero_stress_energy_eV=float(args.minimum_zero_stress_energy_eV),
            minimum_characteristic_stress_GPa=float(
                args.minimum_characteristic_stress_GPa
            ),
        )
        variants.append(row)
    return variants, rejected


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    pd.DataFrame(rows).to_csv(path, index=False)


def main() -> int:
    args = parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    base_row = read_candidate(args.candidate_registry, args.candidate_id)
    physics, physics_metadata = load_physics(args.base_physics_json)
    loading_map = RCurveLoadingMap.from_dict(json.loads(args.loading_map.read_text()))
    coverage_m = float(np.sum(np.asarray(loading_map.projected_advances_m)))
    target_m = float(args.target_extension_um) * 1.0e-6
    if target_m > coverage_m + 1.0e-15:
        raise RuntimeError(
            f"requested target {target_m * 1.0e6:.9g} um exceeds "
            f"loading-map coverage {coverage_m * 1.0e6:.9g} um"
        )

    base_temperatures = tuple(
        float(value) for value in args.identity_base_temperatures_K
    )
    baseline_identity = {
        temperature: run_curve(base_row, physics, loading_map, temperature, args)
        for temperature in base_temperatures
    }
    identity_rows: list[dict[str, Any]] = []
    identity_events: list[dict[str, Any]] = []
    scaled_candidates: list[dict[str, Any]] = []
    for target_peak in args.target_peak_temperatures_K:
        scale = float(target_peak) / float(args.current_peak_temperature_K)
        variant_id = f"temperature_scale_{scale:.8g}"
        transformed = temperature_scale_candidate_row(
            base_row,
            scale,
            candidate_id=f"{base_row['candidate_id']}__{variant_id}",
        )
        scaled_candidates.append(active_parameter_row(transformed))
        for original_temperature in base_temperatures:
            scaled_temperature = scale * original_temperature
            result = run_curve(
                transformed,
                physics,
                loading_map,
                scaled_temperature,
                args,
            )
            differences = max_event_differences(
                baseline_identity[original_temperature],
                result,
            )
            identity_rows.append(
                {
                    "variant_id": variant_id,
                    "temperature_scale": scale,
                    "original_temperature_K": original_temperature,
                    "scaled_temperature_K": scaled_temperature,
                    "original_K50_MPa_sqrt_m": baseline_identity[
                        original_temperature
                    ].checkpoint_K(50.0e-6),
                    "scaled_K50_MPa_sqrt_m": result.checkpoint_K(50.0e-6),
                    "max_abs_event_K_difference_MPa_sqrt_m": differences[
                        "max_abs_diff__K_MPa_sqrt_m"
                    ],
                    **differences,
                }
            )
            identity_events.extend(
                event_rows(result, variant_id, "exact_temperature_axis_scale")
            )

    variants, rejected = create_shelf_variants(base_row, args)
    temperatures = tuple(float(value) for value in args.evaluation_temperatures_K)
    shelf_rows: list[dict[str, Any]] = []
    shelf_events: list[dict[str, Any]] = []
    aggregates: list[dict[str, Any]] = []
    shelf_candidates: list[dict[str, Any]] = []
    for variant in variants:
        variant_id = str(variant["variant_id"])
        transform_kind = str(variant["transform_kind"])
        shelf_candidates.append(
            {"variant_id": variant_id, **active_parameter_row(variant)}
        )
        local: dict[float, dict[str, Any]] = {}
        for temperature in temperatures:
            result = run_curve(variant, physics, loading_map, temperature, args)
            summary = result.summary_dict()
            summary["variant_id"] = variant_id
            summary["transform_kind"] = transform_kind
            summary["delta_K_50_minus_first_MPa_sqrt_m"] = (
                float(summary["K_50um_MPa_sqrt_m"])
                - float(summary["K_first_MPa_sqrt_m"])
            )
            shelf_rows.append(summary)
            local[temperature] = summary
            shelf_events.extend(event_rows(result, variant_id, transform_kind))
        shelf = local[float(args.shelf_temperature_K)]
        anchor = local[float(args.anchor_temperature_K)]
        aggregates.append(
            {
                "variant_id": variant_id,
                "transform_kind": transform_kind,
                "shelf_K_first_MPa_sqrt_m": shelf["K_first_MPa_sqrt_m"],
                "shelf_K50_MPa_sqrt_m": shelf["K_50um_MPa_sqrt_m"],
                "shelf_delta_K50_MPa_sqrt_m": shelf[
                    "delta_K_50_minus_first_MPa_sqrt_m"
                ],
                "anchor_K_first_MPa_sqrt_m": anchor["K_first_MPa_sqrt_m"],
                "anchor_K50_MPa_sqrt_m": anchor["K_50um_MPa_sqrt_m"],
                "anchor_delta_K50_MPa_sqrt_m": anchor[
                    "delta_K_50_minus_first_MPa_sqrt_m"
                ],
                "anchor_max_backstress_GPa": anchor["max_backstress_GPa"],
                "anchor_max_tip_radius_um": anchor["max_tip_radius_um"],
                "anchor_min_front_width_um": anchor["min_front_width_um"],
            }
        )

    baseline = next(row for row in aggregates if row["variant_id"] == "baseline")
    for row in aggregates:
        row["shelf_K50_change_MPa_sqrt_m"] = float(
            row["shelf_K50_MPa_sqrt_m"]
        ) - float(baseline["shelf_K50_MPa_sqrt_m"])
        row["anchor_delta_K50_fraction_of_baseline"] = float(
            row["anchor_delta_K50_MPa_sqrt_m"]
        ) / max(float(baseline["anchor_delta_K50_MPa_sqrt_m"]), 1.0e-30)
        row["anchor_backstress_fraction_of_baseline"] = float(
            row["anchor_max_backstress_GPa"]
        ) / max(float(baseline["anchor_max_backstress_GPa"]), 1.0e-30)

    write_csv(args.out / "temperature_scale_identity.csv", identity_rows)
    write_csv(args.out / "temperature_scale_events.csv", identity_events)
    write_csv(args.out / "temperature_scaled_candidates.csv", scaled_candidates)
    write_csv(args.out / "shelf_scan_aggregate.csv", aggregates)
    write_csv(args.out / "shelf_scan_temperature_summary.csv", shelf_rows)
    write_csv(args.out / "shelf_scan_events.csv", shelf_events)
    write_csv(
        args.out / "shelf_scan_candidates_and_rejections.csv",
        shelf_candidates + rejected,
    )
    manifest = {
        "schema": "v9.13_dbtt_temperature_shelf_diagnostic_v1",
        "candidate_id": args.candidate_id,
        "loading_map_coverage_um": coverage_m * 1.0e6,
        "target_extension_um": args.target_extension_um,
        "translation_action_exponent": args.translation_action_exponent,
        "max_hazard_increment": args.max_hazard_increment,
        "physics_metadata": physics_metadata,
        "temperature_identity_cases": len(identity_rows),
        "shelf_variants_completed": len(aggregates),
        "shelf_variants_rejected": len(rejected),
    }
    (args.out / "diagnostic_manifest.json").write_text(
        json.dumps(manifest, indent=2, allow_nan=False) + "\n"
    )
    maximum_error = max(
        (
            float(row["max_abs_event_K_difference_MPa_sqrt_m"])
            for row in identity_rows
        ),
        default=0.0,
    )
    print(
        "V913_DBTT_TRANSFORM_COMPLETE "
        f"identity_cases={len(identity_rows)} "
        f"max_event_K_difference={maximum_error:.9g} "
        f"shelf_variants={len(aggregates)} out={args.out}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
