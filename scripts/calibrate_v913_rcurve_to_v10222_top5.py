#!/usr/bin/env python3
"""Fit shared v9.13 constants to autonomous v10.2.22 R-curve targets.

Candidate registry rows are fingerprinted and never changed.  The optimizer
may alter only explicitly named shared one-dimensional constants or geometry
scales.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, dataclass, replace
import json
import math
from pathlib import Path
import time
from typing import Any, Mapping, Sequence

import numpy as np
from scipy.optimize import least_squares

from arrhenius_fracture.emergent_gnd_campaign_v913 import (
    candidate_from_registry_row,
)
from arrhenius_fracture.emergent_gnd_rcurve_v913 import (
    RCurveLoadingMap,
    RCurveResult,
    event_rows,
    run_autonomous_rcurve,
)
from arrhenius_fracture.emergent_gnd_types_v913 import (
    CandidateParameters,
    CommonPhysics,
)
from arrhenius_fracture.emergent_gnd_contract_v913 import (
    candidate_parameter_fingerprint,
)
from scripts.run_mpz_v9_13_persistent_top5 import load_physics


PARAMETER_BOUNDS = {
    "K_geometry_scale": (0.95, 1.05),
    "reference_source_area_scale": (0.5, 2.0),
    "emission_geometry_scale": (0.5, 2.0),
    "activation_to_line_content_scale": (0.5, 2.0),
    "persistent_backstress_scale": (0.5, 2.0),
    "blunting_slip_fraction_scale": (0.5, 2.0),
    "translation_action_exponent": (0.25, 2.0),
}
PARAMETER_PRIORS = {name: 1.0 for name in PARAMETER_BOUNDS}

DEFAULT_OPTIMIZE_PARAMETERS = ("translation_action_exponent",)

K_TARGET_FIELDS = (
    "K_first_MPa_sqrt_m",
    "K_10um_MPa_sqrt_m",
    "K_25um_MPa_sqrt_m",
    "K_50um_MPa_sqrt_m",
)


@dataclass(frozen=True)
class SharedCalibration:
    K_geometry_scale: float = 1.0
    reference_source_area_scale: float = 1.0
    emission_geometry_scale: float = 1.0
    activation_to_line_content_scale: float = 1.0
    persistent_backstress_scale: float = 1.0
    blunting_slip_fraction_scale: float = 1.0
    translation_action_exponent: float = 1.0

    def validate(self) -> None:
        for name, value in asdict(self).items():
            lo, hi = PARAMETER_BOUNDS[name]
            if not math.isfinite(value) or not lo <= value <= hi:
                raise ValueError(f"{name}={value} lies outside [{lo}, {hi}]")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mode",
        choices=("evaluate", "grid", "optimize"),
        default="evaluate",
    )
    parser.add_argument("--candidate-registry", type=Path, required=True)
    parser.add_argument("--base-physics-json", type=Path, required=True)
    parser.add_argument("--targets", type=Path, required=True)
    parser.add_argument(
        "--event-targets",
        type=Path,
        help="Optional event-resolved K-versus-extension target table.",
    )
    parser.add_argument("--target-manifest", type=Path, required=True)
    parser.add_argument("--loading-map", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--candidates", nargs="*", default=())
    parser.add_argument("--temperatures", nargs="*", type=float, default=())
    parser.add_argument(
        "--cases",
        nargs="*",
        default=(),
        metavar="CANDIDATE_ID:T_K",
        help=(
            "Optional exact candidate/temperature pairs. When supplied, "
            "these replace the candidates-by-temperatures cross product."
        ),
    )
    parser.add_argument("--train-candidates", nargs="*", default=())
    parser.add_argument(
        "--train-temperatures",
        nargs="*",
        type=float,
        default=(900.0, 1000.0, 1100.0, 1200.0),
    )
    parser.add_argument(
        "--train-cases",
        nargs="*",
        default=(),
        metavar="CANDIDATE_ID:T_K",
        help=(
            "Optional exact training pairs. When supplied, these replace "
            "the training candidates-by-temperatures cross product."
        ),
    )
    parser.add_argument(
        "--optimize-parameters",
        nargs="*",
        choices=tuple(PARAMETER_BOUNDS),
        default=DEFAULT_OPTIMIZE_PARAMETERS,
    )
    parser.add_argument("--max-nfev", type=int, default=12)
    parser.add_argument("--max-hazard-increment", type=float, default=0.05)
    parser.add_argument(
        "--translation-mode",
        choices=("hazard_coupled", "event_commit"),
        default="hazard_coupled",
    )
    parser.add_argument("--target-extension-um", type=float, default=50.0)
    parser.add_argument("--translation-action-exponent", type=float, default=1.0)
    parser.add_argument(
        "--translation-exponent-grid",
        nargs="*",
        type=float,
        default=(0.5, 0.65, 0.8, 0.9, 0.95, 1.0),
    )
    parser.add_argument("--K-residual-scale", type=float, default=2.0)
    parser.add_argument("--first-K-residual-scale", type=float, default=0.5)
    parser.add_argument("--event-residual-weight", type=float, default=0.25)
    parser.add_argument("--state-residual-weight", type=float, default=0.25)
    parser.add_argument("--prior-log10-scale", type=float, default=0.20)
    parser.add_argument("--write-case-json", action="store_true")
    return parser.parse_args()


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise RuntimeError(f"empty CSV: {path}")
    return rows


def _write_csv(path: Path, rows: list[Mapping[str, Any]]) -> None:
    if not rows:
        raise RuntimeError(f"refusing to write empty CSV: {path}")
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=fields,
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def _candidate_fingerprint(rows: list[Mapping[str, str]]) -> str:
    return candidate_parameter_fingerprint(rows)


def _load_candidates(
    registry_path: Path,
    target_manifest_path: Path,
) -> tuple[list[dict[str, str]], dict[str, CandidateParameters], str]:
    rows = _read_csv(registry_path)
    fingerprint = _candidate_fingerprint(rows)
    manifest = json.loads(target_manifest_path.read_text())
    expected = str(manifest["candidate_registry_sha256"])
    if fingerprint != expected:
        raise RuntimeError(
            "candidate registry changed after target extraction: "
            f"expected {expected}, got {fingerprint}"
        )
    candidates = {row["candidate_id"]: candidate_from_registry_row(row) for row in rows}
    return rows, candidates, fingerprint


def _load_targets(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for raw in _read_csv(path):
        row: dict[str, Any] = dict(raw)
        row["temperature_K"] = float(raw["temperature_K"])
        row["seed"] = int(raw["seed"])
        for field in K_TARGET_FIELDS + (
            "max_backstress_GPa",
            "min_front_width_um",
            "max_tip_radius_um",
        ):
            row[field] = float(raw[field])
        rows.append(row)
    return rows


def _load_event_targets(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for raw in _read_csv(path):
        row: dict[str, Any] = dict(raw)
        row["temperature_K"] = float(raw["temperature_K"])
        row["seed"] = int(raw["seed"])
        row["event_index"] = int(raw["event_index"])
        row["K_MPa_sqrt_m"] = float(raw["K_MPa_sqrt_m"])
        row["projected_extension_m"] = float(raw["projected_extension_m"])
        rows.append(row)
    return rows


def _select_targets(
    targets: Sequence[Mapping[str, Any]],
    *,
    candidate_ids: Sequence[str],
    temperatures: Sequence[float],
    case_specs: Sequence[str] = (),
) -> list[dict[str, Any]]:
    exact_cases: set[tuple[str, float]] = set()
    for spec in case_specs:
        try:
            candidate_id, temperature_text = spec.rsplit(":", 1)
            exact_cases.add((candidate_id, float(temperature_text)))
        except (ValueError, TypeError) as exc:
            raise ValueError(
                f"invalid case {spec!r}; expected CANDIDATE_ID:T_K"
            ) from exc
    candidate_filter = set(candidate_ids)
    temperature_filter = {float(value) for value in temperatures}

    def included(row: Mapping[str, Any]) -> bool:
        key = (
            str(row["candidate_id"]),
            float(row["temperature_K"]),
        )
        if exact_cases:
            return key in exact_cases
        return (not candidate_filter or row["candidate_id"] in candidate_filter) and (
            not temperature_filter or float(row["temperature_K"]) in temperature_filter
        )

    selected = [dict(row) for row in targets if included(row)]
    if not selected:
        raise RuntimeError("case selection produced no target rows")
    return sorted(
        selected,
        key=lambda row: (row["candidate_id"], row["temperature_K"]),
    )


def _scaled_inputs(
    base_physics: CommonPhysics,
    base_loading_map: RCurveLoadingMap,
    calibration: SharedCalibration,
) -> tuple[CommonPhysics, RCurveLoadingMap]:
    calibration.validate()
    geometry = tuple(
        tuple(calibration.emission_geometry_scale * float(value) for value in row)
        for row in base_physics.emission_geometry_factors
    )
    physics = replace(
        base_physics,
        reference_source_area_m2=(
            base_physics.reference_source_area_m2
            * calibration.reference_source_area_scale
        ),
        emission_geometry_factors=geometry,
        activation_to_line_content_per_system=tuple(
            calibration.activation_to_line_content_scale * float(value)
            for value in base_physics.activation_to_line_content_per_system
        ),
        persistent_backstress_scale=(
            base_physics.persistent_backstress_scale
            * calibration.persistent_backstress_scale
        ),
        blunting_slip_fraction=(
            base_physics.blunting_slip_fraction
            * calibration.blunting_slip_fraction_scale
        ),
    )
    loading = replace(
        base_loading_map,
        K_per_U_MPa_sqrt_m_per_m=tuple(
            calibration.K_geometry_scale * float(value)
            for value in base_loading_map.K_per_U_MPa_sqrt_m_per_m
        ),
    )
    physics.validate()
    loading.validate()
    return physics, loading


def _evaluate_cases(
    selected: Sequence[Mapping[str, Any]],
    candidates: Mapping[str, CandidateParameters],
    physics: CommonPhysics,
    loading_map: RCurveLoadingMap,
    *,
    target_extension_m: float,
    max_hazard_increment: float,
    translation_mode: str,
    translation_action_exponent: float,
    progress_prefix: str,
) -> list[RCurveResult]:
    results: list[RCurveResult] = []
    count = len(selected)
    for index, target in enumerate(selected, start=1):
        candidate_id = str(target["candidate_id"])
        temperature = float(target["temperature_K"])
        started = time.perf_counter()
        print(
            f"{progress_prefix}_CASE_START index={index}/{count} "
            f"candidate={candidate_id} T={temperature:g}",
            flush=True,
        )
        result = run_autonomous_rcurve(
            candidates[candidate_id],
            physics,
            loading_map,
            temperature,
            target_projected_extension_m=target_extension_m,
            max_hazard_increment=max_hazard_increment,
            translation_mode=translation_mode,
            translation_action_exponent=translation_action_exponent,
        )
        results.append(result)
        print(
            f"{progress_prefix}_CASE_COMPLETE candidate={candidate_id} "
            f"T={temperature:g} status={result.status} "
            f"Kfirst={result.summary_dict()['K_first_MPa_sqrt_m']:.8g} "
            f"K50={result.summary_dict()['K_50um_MPa_sqrt_m']:.8g} "
            f"wall_s={time.perf_counter() - started:.3f}",
            flush=True,
        )
    return results


def _comparison_rows(
    selected: Sequence[Mapping[str, Any]],
    results: Sequence[RCurveResult],
    calibration: SharedCalibration,
    *,
    dataset: str,
) -> list[dict[str, Any]]:
    by_key = {
        (result.candidate_id, float(result.temperature_K)): result for result in results
    }
    rows: list[dict[str, Any]] = []
    for target in selected:
        key = (str(target["candidate_id"]), float(target["temperature_K"]))
        result = by_key[key]
        predicted = result.summary_dict()
        row: dict[str, Any] = {
            "dataset": dataset,
            "candidate_id": key[0],
            "temperature_K": key[1],
            "status": result.status,
            **{f"scale_{name}": value for name, value in asdict(calibration).items()},
        }
        for field in K_TARGET_FIELDS:
            target_value = float(target[field])
            prediction = float(predicted[field])
            row[f"target_{field}"] = target_value
            row[f"predicted_{field}"] = prediction
            row[f"error_{field}"] = prediction - target_value
            row[f"absolute_error_{field}"] = abs(prediction - target_value)
        for field in (
            "max_backstress_GPa",
            "min_front_width_um",
            "max_tip_radius_um",
        ):
            target_value = float(target[field])
            prediction = float(predicted[field])
            row[f"target_{field}"] = target_value
            row[f"predicted_{field}"] = prediction
            row[f"error_{field}"] = prediction - target_value
        rows.append(row)
    return rows


def _event_comparison_rows(
    selected: Sequence[Mapping[str, Any]],
    results: Sequence[RCurveResult],
    calibration: SharedCalibration,
    *,
    dataset: str,
) -> list[dict[str, Any]]:
    by_key = {
        (result.candidate_id, float(result.temperature_K)): result for result in results
    }
    rows: list[dict[str, Any]] = []
    for target in selected:
        key = (str(target["candidate_id"]), float(target["temperature_K"]))
        result = by_key[key]
        event_index = int(target["event_index"])
        if event_index >= len(result.events):
            raise RuntimeError(f"{key} ended before target event {event_index}")
        prediction = result.events[event_index]
        target_extension = float(target["projected_extension_m"])
        if not np.isclose(
            prediction.cumulative_projected_extension_m,
            target_extension,
            rtol=1.0e-10,
            atol=1.0e-15,
        ):
            raise RuntimeError(
                f"{key} event {event_index} extension changed: "
                f"target={target_extension}, "
                f"prediction={prediction.cumulative_projected_extension_m}"
            )
        target_K = float(target["K_MPa_sqrt_m"])
        predicted_K = float(prediction.K_MPa_sqrt_m)
        rows.append(
            {
                "dataset": dataset,
                "candidate_id": key[0],
                "temperature_K": key[1],
                "event_index": event_index,
                "projected_extension_um": target_extension * 1.0e6,
                "target_K_MPa_sqrt_m": target_K,
                "predicted_K_MPa_sqrt_m": predicted_K,
                "error_K_MPa_sqrt_m": predicted_K - target_K,
                "absolute_error_K_MPa_sqrt_m": abs(predicted_K - target_K),
                **{
                    f"scale_{name}": value
                    for name, value in asdict(calibration).items()
                },
            }
        )
    return rows


def _metric_summary(rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    summary: dict[str, Any] = {"n_cases": len(rows)}
    all_errors: list[float] = []
    for field in K_TARGET_FIELDS:
        values = np.asarray(
            [float(row[f"error_{field}"]) for row in rows],
            dtype=float,
        )
        all_errors.extend(values.tolist())
        summary[field] = {
            "mean_error": float(np.mean(values)),
            "MAE": float(np.mean(np.abs(values))),
            "RMSE": float(np.sqrt(np.mean(values**2))),
            "max_absolute_error": float(np.max(np.abs(values))),
        }
    values = np.asarray(all_errors, dtype=float)
    summary["all_K_checkpoints"] = {
        "MAE": float(np.mean(np.abs(values))),
        "RMSE": float(np.sqrt(np.mean(values**2))),
        "max_absolute_error": float(np.max(np.abs(values))),
    }
    per_candidate: dict[str, Any] = {}
    for candidate_id in sorted({str(row["candidate_id"]) for row in rows}):
        candidate_errors = np.asarray(
            [
                float(row[f"error_{field}"])
                for row in rows
                if row["candidate_id"] == candidate_id
                for field in K_TARGET_FIELDS
            ]
        )
        per_candidate[candidate_id] = {
            "MAE": float(np.mean(np.abs(candidate_errors))),
            "RMSE": float(np.sqrt(np.mean(candidate_errors**2))),
            "max_absolute_error": float(np.max(np.abs(candidate_errors))),
        }
    summary["per_candidate"] = per_candidate
    per_temperature: dict[str, Any] = {}
    for temperature in sorted({float(row["temperature_K"]) for row in rows}):
        temperature_errors = np.asarray(
            [
                float(row[f"error_{field}"])
                for row in rows
                if float(row["temperature_K"]) == temperature
                for field in K_TARGET_FIELDS
            ]
        )
        per_temperature[f"{temperature:g}"] = {
            "MAE": float(np.mean(np.abs(temperature_errors))),
            "RMSE": float(np.sqrt(np.mean(temperature_errors**2))),
            "max_absolute_error": float(np.max(np.abs(temperature_errors))),
        }
    summary["per_temperature"] = per_temperature
    return summary


def _event_metric_summary(
    rows: Sequence[Mapping[str, Any]],
) -> dict[str, float | int]:
    if not rows:
        return {"n_events": 0}
    errors = np.asarray(
        [float(row["error_K_MPa_sqrt_m"]) for row in rows],
        dtype=float,
    )
    return {
        "n_events": len(rows),
        "MAE": float(np.mean(np.abs(errors))),
        "RMSE": float(np.sqrt(np.mean(errors**2))),
        "max_absolute_error": float(np.max(np.abs(errors))),
    }


def _residuals(
    selected: Sequence[Mapping[str, Any]],
    selected_events: Sequence[Mapping[str, Any]],
    results: Sequence[RCurveResult],
    calibration: SharedCalibration,
    *,
    K_scale: float,
    first_K_scale: float,
    event_weight: float,
    state_weight: float,
    optimized_parameters: Sequence[str],
    prior_log10_scale: float,
) -> np.ndarray:
    by_key = {
        (result.candidate_id, float(result.temperature_K)): result for result in results
    }
    residual: list[float] = []
    for target in selected:
        result = by_key[(str(target["candidate_id"]), float(target["temperature_K"]))]
        predicted = result.summary_dict()
        for field in K_TARGET_FIELDS:
            scale = first_K_scale if field.startswith("K_first") else K_scale
            residual.append((float(predicted[field]) - float(target[field])) / scale)
        if state_weight > 0.0:
            residual.extend(
                [
                    state_weight
                    * (
                        float(predicted["max_backstress_GPa"])
                        - float(target["max_backstress_GPa"])
                    )
                    / 0.5,
                    state_weight
                    * (
                        float(predicted["min_front_width_um"])
                        - float(target["min_front_width_um"])
                    )
                    / 0.25,
                    state_weight
                    * (
                        float(predicted["max_tip_radius_um"])
                        - float(target["max_tip_radius_um"])
                    )
                    / 0.25,
                ]
            )
    if event_weight > 0.0:
        for target in selected_events:
            result = by_key[
                (str(target["candidate_id"]), float(target["temperature_K"]))
            ]
            event_index = int(target["event_index"])
            if event_index >= len(result.events):
                residual.append(50.0 * event_weight)
                continue
            residual.append(
                event_weight
                * (
                    float(result.events[event_index].K_MPa_sqrt_m)
                    - float(target["K_MPa_sqrt_m"])
                )
                / K_scale
            )
    for name in optimized_parameters:
        residual.append(
            math.log10(getattr(calibration, name) / PARAMETER_PRIORS[name])
            / prior_log10_scale
        )
    return np.asarray(residual, dtype=float)


def _vector_to_calibration(
    vector: np.ndarray,
    optimized_parameters: Sequence[str],
    base: SharedCalibration | None = None,
) -> SharedCalibration:
    values = asdict(base if base is not None else SharedCalibration())
    values.update(
        {
            name: 10.0 ** float(value)
            for name, value in zip(optimized_parameters, vector)
        }
    )
    return SharedCalibration(**values)


def _calibration_to_vector(
    calibration: SharedCalibration,
    optimized_parameters: Sequence[str],
) -> np.ndarray:
    return np.asarray(
        [math.log10(getattr(calibration, name)) for name in optimized_parameters],
        dtype=float,
    )


def _bounds(optimized_parameters: Sequence[str]) -> tuple[np.ndarray, np.ndarray]:
    lo = [math.log10(PARAMETER_BOUNDS[name][0]) for name in optimized_parameters]
    hi = [math.log10(PARAMETER_BOUNDS[name][1]) for name in optimized_parameters]
    return np.asarray(lo), np.asarray(hi)


def _physics_payload(
    physics: CommonPhysics,
    base_payload: Mapping[str, Any],
    calibration: SharedCalibration,
    fingerprint: str,
) -> dict[str, Any]:
    common: dict[str, Any] = {}
    for name, value in vars(physics).items():
        if isinstance(value, tuple):
            common[name] = [
                list(item) if isinstance(item, tuple) else item for item in value
            ]
        else:
            common[name] = value
    return {
        "schema_version": 2,
        "model": "v9.13_persistent_site_autonomous_R_curve",
        "conversion_provenance": base_payload.get(
            "conversion_provenance",
            "v10.2.22_transfer",
        ),
        "common_physics": common,
        "shared_R_curve_calibration": asdict(calibration),
        "candidate_registry_sha256": fingerprint,
        "candidate_parameters_refit": False,
    }


def main() -> int:
    args = parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    registry_rows, candidates, fingerprint = _load_candidates(
        args.candidate_registry,
        args.target_manifest,
    )
    targets = _load_targets(args.targets)
    event_targets = (
        _load_event_targets(args.event_targets)
        if args.event_targets is not None
        else []
    )
    loading_map = RCurveLoadingMap.from_dict(json.loads(args.loading_map.read_text()))
    base_physics, base_payload = load_physics(args.base_physics_json)
    selected = _select_targets(
        targets,
        candidate_ids=args.candidates,
        temperatures=args.temperatures,
        case_specs=args.cases,
    )
    selected_events = (
        _select_targets(
            event_targets,
            candidate_ids=args.candidates,
            temperatures=args.temperatures,
            case_specs=args.cases,
        )
        if event_targets
        else []
    )
    train_candidate_ids = (
        args.train_candidates if args.train_candidates else sorted(candidates)
    )
    training = _select_targets(
        targets,
        candidate_ids=train_candidate_ids,
        temperatures=args.train_temperatures,
        case_specs=args.train_cases,
    )
    training_events = (
        _select_targets(
            event_targets,
            candidate_ids=train_candidate_ids,
            temperatures=args.train_temperatures,
            case_specs=args.train_cases,
        )
        if event_targets
        else []
    )
    optimized_parameters = tuple(args.optimize_parameters)
    history: list[dict[str, Any]] = []
    base_calibration = SharedCalibration(
        translation_action_exponent=args.translation_action_exponent,
    )
    calibration = base_calibration

    if args.mode == "grid":
        if len(optimized_parameters) != 1:
            raise RuntimeError(
                "--mode grid requires exactly one --optimize-parameters value"
            )
        grid_parameter = optimized_parameters[0]
        grid_values = (
            args.translation_exponent_grid
            if grid_parameter == "translation_action_exponent"
            else ()
        )
        if not grid_values:
            raise RuntimeError("--mode grid supports translation_action_exponent")
        best: tuple[float, SharedCalibration] | None = None
        for evaluation_index, grid_value in enumerate(
            grid_values,
            start=1,
        ):
            trial = replace(
                base_calibration,
                **{grid_parameter: float(grid_value)},
            )
            physics, loading = _scaled_inputs(
                base_physics,
                loading_map,
                trial,
            )
            started = time.perf_counter()
            results = _evaluate_cases(
                training,
                candidates,
                physics,
                loading,
                target_extension_m=args.target_extension_um * 1.0e-6,
                max_hazard_increment=args.max_hazard_increment,
                translation_mode=args.translation_mode,
                translation_action_exponent=(trial.translation_action_exponent),
                progress_prefix=f"V913_RCURVE_GRID_{evaluation_index}",
            )
            values = _residuals(
                training,
                training_events,
                results,
                trial,
                K_scale=args.K_residual_scale,
                first_K_scale=args.first_K_residual_scale,
                event_weight=args.event_residual_weight,
                state_weight=args.state_residual_weight,
                optimized_parameters=optimized_parameters,
                prior_log10_scale=args.prior_log10_scale,
            )
            objective_value = float(np.dot(values, values))
            history.append(
                {
                    "evaluation": evaluation_index,
                    "objective": objective_value,
                    "residual_RMS": float(np.sqrt(np.mean(values**2))),
                    "wall_s": time.perf_counter() - started,
                    **asdict(trial),
                }
            )
            if best is None or objective_value < best[0]:
                best = (objective_value, trial)
        if best is None:
            raise RuntimeError("calibration grid is empty")
        _write_csv(args.out / "objective_history.csv", history)
        calibration = best[1]
        optimizer_record = {
            "success": True,
            "status": 0,
            "message": "finite_grid_complete",
            "nfev": len(history),
            "best_objective": best[0],
        }
    elif args.mode == "optimize":
        if not optimized_parameters:
            raise RuntimeError("--mode optimize requires at least one parameter")
        evaluation_index = 0

        def objective(vector: np.ndarray) -> np.ndarray:
            nonlocal evaluation_index
            evaluation_index += 1
            trial = _vector_to_calibration(
                vector,
                optimized_parameters,
                base_calibration,
            )
            physics, loading = _scaled_inputs(
                base_physics,
                loading_map,
                trial,
            )
            started = time.perf_counter()
            print(
                "V913_RCURVE_OBJECTIVE_START "
                f"evaluation={evaluation_index} scales={asdict(trial)}",
                flush=True,
            )
            results = _evaluate_cases(
                training,
                candidates,
                physics,
                loading,
                target_extension_m=args.target_extension_um * 1.0e-6,
                max_hazard_increment=args.max_hazard_increment,
                translation_mode=args.translation_mode,
                translation_action_exponent=(trial.translation_action_exponent),
                progress_prefix=f"V913_RCURVE_OBJECTIVE_{evaluation_index}",
            )
            values = _residuals(
                training,
                training_events,
                results,
                trial,
                K_scale=args.K_residual_scale,
                first_K_scale=args.first_K_residual_scale,
                event_weight=args.event_residual_weight,
                state_weight=args.state_residual_weight,
                optimized_parameters=optimized_parameters,
                prior_log10_scale=args.prior_log10_scale,
            )
            objective_value = float(np.dot(values, values))
            history.append(
                {
                    "evaluation": evaluation_index,
                    "objective": objective_value,
                    "residual_RMS": float(np.sqrt(np.mean(values**2))),
                    "wall_s": time.perf_counter() - started,
                    **asdict(trial),
                }
            )
            _write_csv(args.out / "objective_history.csv", history)
            print(
                "V913_RCURVE_OBJECTIVE_COMPLETE "
                f"evaluation={evaluation_index} objective={objective_value:.9g}",
                flush=True,
            )
            return values

        x0 = _calibration_to_vector(base_calibration, optimized_parameters)
        lower, upper = _bounds(optimized_parameters)
        fit = least_squares(
            objective,
            x0,
            bounds=(lower, upper),
            max_nfev=args.max_nfev,
            verbose=2,
        )
        calibration = _vector_to_calibration(
            fit.x,
            optimized_parameters,
            base_calibration,
        )
        optimizer_record = {
            "success": bool(fit.success),
            "status": int(fit.status),
            "message": str(fit.message),
            "nfev": int(fit.nfev),
            "cost": float(fit.cost),
            "optimality": float(fit.optimality),
            "active_mask": fit.active_mask.tolist(),
        }
    else:
        optimizer_record = {
            "success": True,
            "status": 0,
            "message": "evaluation_only",
            "nfev": 0,
        }

    physics, calibrated_loading = _scaled_inputs(
        base_physics,
        loading_map,
        calibration,
    )
    results = _evaluate_cases(
        selected,
        candidates,
        physics,
        calibrated_loading,
        target_extension_m=args.target_extension_um * 1.0e-6,
        max_hazard_increment=args.max_hazard_increment,
        translation_mode=args.translation_mode,
        translation_action_exponent=(calibration.translation_action_exponent),
        progress_prefix="V913_RCURVE_FINAL",
    )
    comparisons = _comparison_rows(
        selected,
        results,
        calibration,
        dataset="final",
    )
    _write_csv(args.out / "R_curve_checkpoint_comparison.csv", comparisons)
    event_comparisons = (
        _event_comparison_rows(
            selected_events,
            results,
            calibration,
            dataset="final",
        )
        if selected_events
        else []
    )
    if event_comparisons:
        _write_csv(
            args.out / "R_curve_event_comparison.csv",
            event_comparisons,
        )
    _write_csv(
        args.out / "R_curve_event_predictions.csv",
        event_rows(results),
    )
    if args.write_case_json:
        case_root = args.out / "cases"
        case_root.mkdir(parents=True, exist_ok=True)
        for result in results:
            path = case_root / (
                f"{result.candidate_id}_T{int(round(result.temperature_K))}K.json"
            )
            path.write_text(json.dumps(result.as_dict(), indent=2))

    physics_record = _physics_payload(
        physics,
        base_payload,
        calibration,
        fingerprint,
    )
    (args.out / "v9_13_R_curve_calibrated_common_physics.json").write_text(
        json.dumps(physics_record, indent=2)
    )
    (args.out / "v9_13_R_curve_calibrated_loading_map.json").write_text(
        json.dumps(calibrated_loading.as_dict(), indent=2)
    )
    summary = {
        "schema": "v9.13_autonomous_R_curve_calibration",
        "mode": args.mode,
        "candidate_parameters_refit": False,
        "candidate_registry_sha256": fingerprint,
        "n_candidate_rows": len(registry_rows),
        "optimized_parameters": list(optimized_parameters),
        "shared_calibration": asdict(calibration),
        "translation_mode": args.translation_mode,
        "training_cases": [
            [row["candidate_id"], row["temperature_K"]] for row in training
        ],
        "final_cases": [
            [row["candidate_id"], row["temperature_K"]] for row in selected
        ],
        "optimizer": optimizer_record,
        "metrics": _metric_summary(comparisons),
        "event_metrics": _event_metric_summary(event_comparisons),
    }
    (args.out / "R_curve_calibration_summary.json").write_text(
        json.dumps(summary, indent=2)
    )
    print(
        "V913_RCURVE_CALIBRATION_COMPLETE "
        f"mode={args.mode} cases={len(results)} "
        f"RMSE={summary['metrics']['all_K_checkpoints']['RMSE']:.8g} "
        f"candidate_parameters_refit=false out={args.out}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
