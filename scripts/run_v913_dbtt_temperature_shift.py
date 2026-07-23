#!/usr/bin/env python3
"""Verify exact temperature-axis shifts of one calibrated v9.13 R-curve.

This entry point intentionally performs no cleavage-shelf modification.  It
runs the original candidate at a base temperature grid and the deterministically
scaled candidate at ``lambda*T`` using the same loading map and CRN thresholds.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys
from typing import Any, Mapping

# A script invoked as ``python scripts/name.py`` receives ``scripts/`` rather
# than the repository root on sys.path.  Add the root before package imports.
REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import numpy as np
import pandas as pd

from arrhenius_fracture.dbtt_transform_v913 import (
    active_parameter_row,
    temperature_scale_candidate_row,
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
from arrhenius_fracture.emergent_gnd_types_v913 import CommonPhysics


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
        "--base-temperatures-K",
        type=float,
        nargs="+",
        default=(700.0, 800.0, 900.0, 1000.0, 1100.0, 1200.0),
    )
    parser.add_argument("--target-extension-um", type=float, default=50.0)
    parser.add_argument("--translation-action-exponent", type=float, default=0.95)
    parser.add_argument("--max-hazard-increment", type=float, default=0.05)
    parser.add_argument(
        "--identity-tolerance-MPa-sqrt-m",
        type=float,
        default=1.0e-8,
    )
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
        "emission_geometry_extension_m",
    ):
        if key in common:
            common[key] = tuple(common[key])
    for key in (
        "forest_interaction_matrix",
        "gnd_stress_projection_matrix",
        "emission_geometry_factors",
    ):
        if key in common:
            common[key] = tuple(tuple(row) for row in common[key])
    physics = CommonPhysics(**common)
    physics.validate()
    return physics, metadata


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
    physics: CommonPhysics,
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
            "temperature scaling changed event count: "
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
    *,
    variant_id: str,
    curve_role: str,
    source_temperature_K: float,
) -> list[dict[str, Any]]:
    return [
        {
            "variant_id": variant_id,
            "curve_role": curve_role,
            "candidate_id": result.candidate_id,
            "source_temperature_K": source_temperature_K,
            "temperature_K": result.temperature_K,
            "status": result.status,
            "seed": result.seed,
            **event.as_dict(),
        }
        for event in result.events
    ]


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
    if float(args.current_peak_temperature_K) <= 0.0:
        raise ValueError("current peak temperature must be positive")
    if float(args.identity_tolerance_MPa_sqrt_m) < 0.0:
        raise ValueError("identity tolerance must be nonnegative")

    base_temperatures = tuple(float(value) for value in args.base_temperatures_K)
    target_peaks = tuple(float(value) for value in args.target_peak_temperatures_K)
    total_cases = len(base_temperatures) * (1 + len(target_peaks))
    completed = 0

    baseline: dict[float, RCurveResult] = {}
    all_events: list[dict[str, Any]] = []
    for temperature in base_temperatures:
        print(
            "V913_DBTT_SHIFT_CASE_START "
            f"role=baseline T={temperature:g} "
            f"case={completed + 1}/{total_cases}",
            flush=True,
        )
        result = run_curve(base_row, physics, loading_map, temperature, args)
        baseline[temperature] = result
        all_events.extend(
            event_rows(
                result,
                variant_id="baseline",
                curve_role="baseline",
                source_temperature_K=temperature,
            )
        )
        completed += 1
        print(
            "V913_DBTT_SHIFT_CASE_COMPLETE "
            f"role=baseline T={temperature:g} status={result.status} "
            f"K50={result.checkpoint_K(50.0e-6):.9g} "
            f"case={completed}/{total_cases}",
            flush=True,
        )

    identity_rows: list[dict[str, Any]] = []
    scaled_candidates: list[dict[str, Any]] = []
    maximum_event_K_error = 0.0

    for target_peak in target_peaks:
        scale = target_peak / float(args.current_peak_temperature_K)
        if scale <= 0.0 or not np.isfinite(scale):
            raise ValueError(f"invalid target peak temperature {target_peak!r}")
        variant_id = f"temperature_scale_{scale:.8g}"
        transformed = temperature_scale_candidate_row(
            base_row,
            scale,
            candidate_id=f"{base_row['candidate_id']}__{variant_id}",
        )
        scaled_candidates.append(active_parameter_row(transformed))

        for source_temperature in base_temperatures:
            shifted_temperature = scale * source_temperature
            print(
                "V913_DBTT_SHIFT_CASE_START "
                f"role={variant_id} source_T={source_temperature:g} "
                f"shifted_T={shifted_temperature:g} "
                f"case={completed + 1}/{total_cases}",
                flush=True,
            )
            result = run_curve(
                transformed,
                physics,
                loading_map,
                shifted_temperature,
                args,
            )
            differences = max_event_differences(baseline[source_temperature], result)
            event_K_error = differences["max_abs_diff__K_MPa_sqrt_m"]
            maximum_event_K_error = max(maximum_event_K_error, event_K_error)
            identity_rows.append(
                {
                    "variant_id": variant_id,
                    "temperature_scale": scale,
                    "target_peak_temperature_K": target_peak,
                    "source_temperature_K": source_temperature,
                    "shifted_temperature_K": shifted_temperature,
                    "baseline_status": baseline[source_temperature].status,
                    "shifted_status": result.status,
                    "baseline_K_first_MPa_sqrt_m": baseline[
                        source_temperature
                    ].checkpoint_K(0.0),
                    "shifted_K_first_MPa_sqrt_m": result.checkpoint_K(0.0),
                    "baseline_K50_MPa_sqrt_m": baseline[
                        source_temperature
                    ].checkpoint_K(50.0e-6),
                    "shifted_K50_MPa_sqrt_m": result.checkpoint_K(50.0e-6),
                    "max_abs_event_K_difference_MPa_sqrt_m": event_K_error,
                    **differences,
                }
            )
            all_events.extend(
                event_rows(
                    result,
                    variant_id=variant_id,
                    curve_role="temperature_scaled",
                    source_temperature_K=source_temperature,
                )
            )
            completed += 1
            print(
                "V913_DBTT_SHIFT_CASE_COMPLETE "
                f"role={variant_id} source_T={source_temperature:g} "
                f"shifted_T={shifted_temperature:g} status={result.status} "
                f"K50={result.checkpoint_K(50.0e-6):.9g} "
                f"event_K_error={event_K_error:.3e} "
                f"case={completed}/{total_cases}",
                flush=True,
            )

    write_csv(args.out / "temperature_scale_identity.csv", identity_rows)
    write_csv(args.out / "temperature_scale_events.csv", all_events)
    write_csv(args.out / "temperature_scaled_candidates.csv", scaled_candidates)

    passed = maximum_event_K_error <= float(args.identity_tolerance_MPa_sqrt_m)
    manifest = {
        "schema": "v9.13_dbtt_temperature_shift_only_v1",
        "candidate_id": args.candidate_id,
        "mode": "exact_temperature_axis_shift_only",
        "current_peak_temperature_K": args.current_peak_temperature_K,
        "target_peak_temperatures_K": list(target_peaks),
        "base_temperatures_K": list(base_temperatures),
        "loading_map_coverage_um": coverage_m * 1.0e6,
        "target_extension_um": args.target_extension_um,
        "translation_action_exponent": args.translation_action_exponent,
        "max_hazard_increment": args.max_hazard_increment,
        "identity_tolerance_MPa_sqrt_m": args.identity_tolerance_MPa_sqrt_m,
        "maximum_event_K_error_MPa_sqrt_m": maximum_event_K_error,
        "identity_passed": passed,
        "physics_metadata": physics_metadata,
        "cases_completed": completed,
    }
    (args.out / "temperature_shift_manifest.json").write_text(
        json.dumps(manifest, indent=2, allow_nan=False) + "\n"
    )

    print(
        "V913_DBTT_TEMPERATURE_SHIFT_COMPLETE "
        f"cases={completed} max_event_K_difference={maximum_event_K_error:.9g} "
        f"identity_passed={int(passed)} out={args.out}",
        flush=True,
    )
    return 0 if passed else 2


if __name__ == "__main__":
    raise SystemExit(main())
