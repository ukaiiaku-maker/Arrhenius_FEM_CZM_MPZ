#!/usr/bin/env python3
"""Transfer the v10.2.22 top-five 2-D geometry into the v9.13 1-D model.

The candidate rows are treated as immutable.  This script extracts only common
normalizations and geometry from a v10.2.22 archive:

* aggregate-activation to line-content conversion;
* encounter efficiency;
* the extension-resolved two-channel emission projection;
* moving-tip convection resolution; and
* the measured fact that direct signed shielding is negligible.

One plastic-free reference trajectory (top candidate at 300 K) supplies the
reduced geometry schedule.  Optional replay validation then holds that schedule
fixed and evaluates all other candidate/temperature histories as holdouts.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import asdict
import json
import math
from pathlib import Path, PurePosixPath
import re
import statistics
from typing import Any, Iterable, Mapping
import zipfile

import numpy as np

from arrhenius_fracture.emergent_gnd_campaign_v913 import (
    candidate_from_registry_row,
)
from arrhenius_fracture.emergent_gnd_state_v913 import EmergentGNDState
from arrhenius_fracture.emergent_gnd_types_v913 import CommonPhysics


STACK_NAME = "v10_2_17_final_signed_stochastic_stack.json"
REFERENCE_CANDIDATE = "v912_targeted_local_peak_013476_0368"
REFERENCE_TEMPERATURE_K = 300


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", required=True)
    parser.add_argument(
        "--base-physics",
        default="mpz_v9_13_persistent_sites_common_physics.json",
    )
    parser.add_argument(
        "--base-registry",
        default="candidates/v9_13_persistent_sites_top5_registry.csv",
    )
    parser.add_argument("--out", required=True)
    parser.add_argument(
        "--skip-replay",
        action="store_true",
        help="write transferred inputs without replaying 2-D accepted histories",
    )
    parser.add_argument(
        "--case-limit",
        type=int,
        default=0,
        help="limit replay cases for a smoke test; zero validates all holdouts",
    )
    return parser.parse_args()


def _read_json(archive: zipfile.ZipFile, member: str) -> Any:
    return json.loads(archive.read(member))


def _stack_members(archive: zipfile.ZipFile) -> list[str]:
    members = [
        name
        for name in archive.namelist()
        if PurePosixPath(name).name == STACK_NAME
    ]
    if not members:
        raise RuntimeError(f"archive contains no {STACK_NAME}")
    return sorted(members)


def _temperature_from_member(member: str) -> int:
    match = re.search(r"/T(\d+)K_", member)
    if match is None:
        raise RuntimeError(f"cannot parse temperature from {member}")
    return int(match.group(1))


def _reference_payload(
    archive: zipfile.ZipFile,
    members: Iterable[str],
) -> tuple[str, dict[str, Any]]:
    for member in members:
        payload = _read_json(archive, member)
        selected = payload["selected_option"]["candidate_id"]
        temperature = _temperature_from_member(member)
        if (
            selected == REFERENCE_CANDIDATE
            and temperature == REFERENCE_TEMPERATURE_K
        ):
            return member, payload
    raise RuntimeError(
        "plastic-free geometry reference case is missing: "
        f"{REFERENCE_CANDIDATE} at {REFERENCE_TEMPERATURE_K} K"
    )


def _same_numeric_or_text(left: Any, right: Any) -> bool:
    try:
        return math.isclose(
            float(left),
            float(right),
            rel_tol=1.0e-12,
            abs_tol=1.0e-15,
        )
    except (TypeError, ValueError):
        return str(left) == str(right)


def _exact_registry_rows(
    archive: zipfile.ZipFile,
    members: Iterable[str],
) -> dict[str, dict[str, str]]:
    exact: dict[str, dict[str, str]] = {}
    for member in members:
        payload = _read_json(archive, member)
        row = {
            str(key): str(value)
            for key, value in payload["selected_option"][
                "exact_registry_row"
            ].items()
        }
        candidate_id = row["candidate_id"]
        if candidate_id in exact:
            previous = exact[candidate_id]
            for key in set(previous) & set(row):
                if not _same_numeric_or_text(previous[key], row[key]):
                    raise RuntimeError(
                        f"candidate row changed across 2-D cases: "
                        f"{candidate_id} field={key}"
                    )
        else:
            exact[candidate_id] = row
    return exact


def _common_tuple(
    payloads: Iterable[Mapping[str, Any]],
    path: tuple[str, ...],
) -> tuple[float, ...]:
    values: list[tuple[float, ...]] = []
    for payload in payloads:
        value: Any = payload
        for key in path:
            value = value[key]
        values.append(tuple(float(item) for item in value))
    first = values[0]
    if any(
        len(value) != len(first)
        or not np.allclose(value, first, rtol=1.0e-12, atol=1.0e-15)
        for value in values[1:]
    ):
        raise RuntimeError(f"2-D common value is not invariant: {path}")
    return first


def _geometry_schedule(
    reference: Mapping[str, Any],
) -> tuple[list[float], list[list[float]]]:
    breakpoints: list[float] = []
    factors: list[list[float]] = []
    previous: np.ndarray | None = None
    for record in reference["records"]:
        current = np.asarray(
            record["anisotropic_drive_factors"],
            dtype=float,
        )
        if previous is not None and np.allclose(
            current,
            previous,
            rtol=1.0e-12,
            atol=1.0e-15,
        ):
            continue
        pre_extension = max(
            float(record["micro_advance_total_m"])
            - float(record["micro_advance_step_m"]),
            0.0,
        )
        if breakpoints and pre_extension <= breakpoints[-1]:
            pre_extension = float(np.nextafter(breakpoints[-1], math.inf))
        breakpoints.append(pre_extension)
        factors.append(current.tolist())
        previous = current
    if not breakpoints or breakpoints[0] > 1.0e-15:
        raise RuntimeError("geometry schedule does not begin at the virgin tip")
    breakpoints[0] = 0.0
    return breakpoints, factors


def _load_base_registry(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        rows = list(reader)
        fields = list(reader.fieldnames or [])
    if not rows:
        raise RuntimeError(f"empty base registry: {path}")
    return fields, rows


def _write_registry(
    path: Path,
    base_path: Path,
    exact_rows: Mapping[str, Mapping[str, str]],
) -> None:
    fields, rows = _load_base_registry(base_path)
    if "encounter_efficiency" not in fields:
        fields.append("encounter_efficiency")
    for row in rows:
        candidate_id = row["candidate_id"]
        exact = exact_rows.get(candidate_id)
        if exact is None:
            raise RuntimeError(
                f"base candidate is absent from 2-D archive: {candidate_id}"
            )
        for key, base_value in row.items():
            if key in exact and not _same_numeric_or_text(
                base_value,
                exact[key],
            ):
                raise RuntimeError(
                    "candidate parameter mismatch; refusing to tune a candidate "
                    f"row: candidate={candidate_id} field={key} "
                    f"base={base_value} archive={exact[key]}"
                )
        row["encounter_efficiency"] = exact["encounter_efficiency"]
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=fields,
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def _physics_from_payload(payload: Mapping[str, Any]) -> CommonPhysics:
    common = dict(payload["common_physics"])
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
    return physics


def _quantile(values: list[float], probability: float) -> float:
    if not values:
        return float("nan")
    return float(np.quantile(np.asarray(values, dtype=float), probability))


def _replay_case(
    payload: Mapping[str, Any],
    physics: CommonPhysics,
    temperature_K: float,
) -> dict[str, Any]:
    row = payload["selected_option"]["exact_registry_row"]
    state = EmergentGNDState(candidate_from_registry_row(row), physics)
    # Each archive record is already one accepted 2-D constitutive interval.
    # Do not subdivide stationary time beyond the accepted-record resolution;
    # moving-tip substeps remain controlled by ``moving_tip_cfl``.
    state.max_feedback_substep_s = float("inf")

    predicted_backstress: list[float] = []
    predicted_width: list[float] = []
    predicted_radius: list[float] = []
    target_records = payload["records"]
    for record in target_records:
        K = float(record["K_Pa_sqrt_m"]) / 1.0e6
        state.advance_coupled_segment(
            duration_s=float(record["kinetic_dt_consumed_s"]),
            da_m=float(record["micro_advance_step_m"]),
            K_start_MPa_sqrt_m=K,
            K_end_MPa_sqrt_m=K,
            T_K=float(temperature_K),
        )
        geometry = state.source_geometry()
        sigma_back = state.backstress_state()[2]
        predicted_backstress.append(float(np.mean(sigma_back)))
        predicted_width.append(float(geometry["front_width_m"]))
        predicted_radius.append(float(geometry["tip_radius_m"]))

    target_backstress = [
        float(record["persistent_sigma_back_Pa"])
        for record in target_records
    ]
    target_width = [
        float(record["persistent_site_front_width_m"])
        for record in target_records
    ]
    target_radius = [
        float(record["persistent_tip_radius_m"])
        for record in target_records
    ]

    predicted = {
        "max_backstress_GPa": max(predicted_backstress) / 1.0e9,
        "min_front_width_um": min(predicted_width) * 1.0e6,
        "max_tip_radius_um": max(predicted_radius) * 1.0e6,
    }
    target = {
        "max_backstress_GPa": max(target_backstress) / 1.0e9,
        "min_front_width_um": min(target_width) * 1.0e6,
        "max_tip_radius_um": max(target_radius) * 1.0e6,
    }
    row_out: dict[str, Any] = {
        "candidate_id": row["candidate_id"],
        "temperature_K": float(temperature_K),
        **{f"predicted_{key}": value for key, value in predicted.items()},
        **{f"target_{key}": value for key, value in target.items()},
    }
    for key in target:
        scale = max(abs(target[key]), 1.0e-12)
        row_out[f"relative_error_{key}"] = (
            predicted[key] - target[key]
        ) / scale
    return row_out


def _write_csv(path: Path, rows: list[Mapping[str, Any]]) -> None:
    if not rows:
        return
    fields = list(rows[0])
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=fields,
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    out = Path(args.out)
    if out.exists():
        raise RuntimeError(f"output already exists: {out}")
    out.mkdir(parents=True)

    archive_path = Path(args.archive)
    with zipfile.ZipFile(archive_path) as archive:
        members = _stack_members(archive)
        cases = [
            (member, _read_json(archive, member))
            for member in members
        ]
        payloads = [payload for _, payload in cases]
        reference_member, reference = _reference_payload(archive, members)
        exact_rows = _exact_registry_rows(archive, members)

    conversion = _common_tuple(
        payloads,
        (
            "signed_burgers_shared_physics",
            "kernel",
            "activation_to_line_content_by_system",
        ),
    )
    encounters = {
        float(
            payload["selected_option"]["exact_registry_row"][
                "encounter_efficiency"
            ]
        )
        for payload in payloads
    }
    if len(encounters) != 1:
        raise RuntimeError("encounter efficiency is not common across candidates")
    encounter_efficiency = encounters.pop()
    breakpoints, factors = _geometry_schedule(reference)

    base_physics = json.loads(Path(args.base_physics).read_text())
    common = base_physics.setdefault("common_physics", {})
    dx = float(common["mpz_length_m"]) / int(common["n_bins"])
    max_translation_substep = float(reference["config"][
        "max_translation_substep_m"
    ])
    common.update(
        {
            "activation_to_line_content_per_system": list(conversion),
            "encounter_efficiency": encounter_efficiency,
            "taylor_phi_max": 20.0,
            "mobile_transport_velocity_scale": 0.0,
            "emission_schmid_factors": list(factors[0]),
            "emission_geometry_extension_m": breakpoints,
            "emission_geometry_factors": factors,
            "shielding_orientation_factors": [0.0, 0.0],
            "minimum_front_width_m": float(common["b_m"]),
            "maximum_front_width_m": float(common["mpz_length_m"]),
            "coupled_moving_tip_enabled": True,
            "moving_tip_cfl": max_translation_substep / dx,
        }
    )
    base_physics["conversion_provenance"] = (
        "v10.2.22_archive_signed_kernel_activation_to_line_content"
    )
    base_physics["transfer_calibration"] = {
        "source_archive": archive_path.name,
        "candidate_parameters_refit": False,
        "geometry_reference_member": reference_member,
        "geometry_reference_candidate": REFERENCE_CANDIDATE,
        "geometry_reference_temperature_K": REFERENCE_TEMPERATURE_K,
        "geometry_reference_reason": "plastic_free_2d_trajectory",
        "emission_projection": (
            "piecewise_constant_in_cumulative_crack_extension"
        ),
        "projection_recomputed": "after_each_accepted_crack_event",
        "analytical_1d_signed_shielding_active": False,
        "shielding_reason": (
            "measured_2d_signed_kernel_is_negligible_relative_to_toughness"
        ),
        "holdout_definition": (
            "all candidate/temperature histories except the single 300 K "
            "geometry reference"
        ),
    }
    physics_path = out / "v9_13_v10222_transfer_common_physics.json"
    physics_path.write_text(json.dumps(base_physics, indent=2) + "\n")
    _write_registry(
        out / "v9_13_v10222_fixed_candidate_registry.csv",
        Path(args.base_registry),
        exact_rows,
    )

    manifest: dict[str, Any] = {
        "schema": "v9.13_to_v10.2.22_top5_transfer_calibration",
        "source_archive": archive_path.name,
        "n_cases": len(payloads),
        "n_candidates": len(exact_rows),
        "candidate_parameters_refit": False,
        "reference_member": reference_member,
        "line_conversion_per_system": list(conversion),
        "encounter_efficiency": encounter_efficiency,
        "moving_tip_cfl": common["moving_tip_cfl"],
        "mobile_transport_velocity_scale": common[
            "mobile_transport_velocity_scale"
        ],
        "geometry_breakpoints_m": breakpoints,
        "geometry_factors": factors,
        "common_physics": asdict(_physics_from_payload(base_physics)),
    }

    if not args.skip_replay:
        physics = _physics_from_payload(base_physics)
        validation_cases = [
            (member, payload)
            for member, payload in cases
            if not (
                payload["selected_option"]["candidate_id"]
                == REFERENCE_CANDIDATE
                and _temperature_from_member(member)
                == REFERENCE_TEMPERATURE_K
            )
        ]
        if args.case_limit > 0:
            validation_cases = validation_cases[: args.case_limit]
        rows: list[dict[str, Any]] = []
        for index, (member, payload) in enumerate(
            validation_cases,
            start=1,
        ):
            temperature = float(_temperature_from_member(member))
            row = _replay_case(payload, physics, temperature)
            rows.append(row)
            print(
                "V913_TRANSFER_REPLAY "
                f"{index}/{len(validation_cases)} "
                f"candidate={row['candidate_id']} T={temperature:g} "
                f"backstress_error={row['relative_error_max_backstress_GPa']:.4g} "
                f"width_error={row['relative_error_min_front_width_um']:.4g} "
                f"radius_error={row['relative_error_max_tip_radius_um']:.4g}",
                flush=True,
            )
        _write_csv(out / "v9_13_v10222_replay_validation.csv", rows)
        active = [
            row
            for row in rows
            if float(row["target_max_backstress_GPa"]) >= 0.05
        ]
        summary: dict[str, Any] = {
            "n_holdout_cases": len(rows),
            "n_active_holdout_cases": len(active),
            "active_case_threshold_backstress_GPa": 0.05,
        }
        for key in (
            "max_backstress_GPa",
            "min_front_width_um",
            "max_tip_radius_um",
        ):
            errors = [
                abs(float(row[f"relative_error_{key}"]))
                for row in active
            ]
            summary[f"median_absolute_relative_error_{key}"] = (
                statistics.median(errors) if errors else float("nan")
            )
            summary[f"p90_absolute_relative_error_{key}"] = _quantile(
                errors,
                0.90,
            )
        manifest["replay_validation"] = summary

    (out / "v9_13_v10222_transfer_manifest.json").write_text(
        json.dumps(manifest, indent=2, default=str) + "\n"
    )
    print(
        "V913_TRANSFER_CALIBRATION_COMPLETE "
        f"out={out} cases={len(payloads)} candidates={len(exact_rows)}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
