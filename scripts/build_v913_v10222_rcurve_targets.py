#!/usr/bin/env python3
"""Extract immutable v10.2.22 R-curve targets for v9.13 calibration."""

from __future__ import annotations

import argparse
import csv
import io
import json
from pathlib import Path
from typing import Any, Iterable, Mapping
import zipfile

import numpy as np

from arrhenius_fracture.emergent_gnd_rcurve_v913 import RCurveLoadingMap
from arrhenius_fracture.emergent_gnd_contract_v913 import (
    CANDIDATE_PARAMETER_FIELDS,
    candidate_parameter_fingerprint,
    effective_candidate_parameters,
)


REFERENCE_CANDIDATE_ID = "v912_targeted_local_peak_013476_0368"
REFERENCE_TEMPERATURE_K = 300

CALIBRATION_BOUNDARY_ROWS = (
    {
        "group": "candidate_parameterization",
        "name": "all_candidate_registry_fields",
        "treatment": "fixed_exactly",
        "nominal": "candidate_registry",
        "lower_bound": "",
        "upper_bound": "",
        "reason": "These are the parameters passed unchanged from 1-D to 2-D.",
    },
    {
        "group": "shared_loading_geometry",
        "name": "K_geometry_scale",
        "treatment": "calibratable_shared",
        "nominal": 1.0,
        "lower_bound": 0.95,
        "upper_bound": 1.05,
        "reason": "Common reduced displacement-to-K normalization.",
    },
    {
        "group": "shared_source_geometry",
        "name": "reference_source_area_scale",
        "treatment": "calibratable_shared",
        "nominal": 1.0,
        "lower_bound": 0.5,
        "upper_bound": 2.0,
        "reason": "One-dimensional active-arc representation of 2-D site area.",
    },
    {
        "group": "shared_emission_geometry",
        "name": "emission_geometry_scale",
        "treatment": "calibratable_shared",
        "nominal": 1.0,
        "lower_bound": 0.5,
        "upper_bound": 2.0,
        "reason": "Common reduction of the two 2-D slip-trace projections.",
    },
    {
        "group": "shared_line_normalization",
        "name": "activation_to_line_content_scale",
        "treatment": "calibratable_shared",
        "nominal": 1.0,
        "lower_bound": 0.5,
        "upper_bound": 2.0,
        "reason": "Shared aggregate-activation to represented-line conversion.",
    },
    {
        "group": "shared_backstress_mapping",
        "name": "persistent_backstress_scale",
        "treatment": "calibratable_shared",
        "nominal": 1.0,
        "lower_bound": 0.5,
        "upper_bound": 2.0,
        "reason": "One-dimensional homogenization of the 2-D Taylor backstress.",
    },
    {
        "group": "shared_blunting_mapping",
        "name": "blunting_slip_fraction_scale",
        "treatment": "calibratable_shared",
        "nominal": 1.0,
        "lower_bound": 0.5,
        "upper_bound": 2.0,
        "reason": "Shared conversion from accumulated line content to tip radius.",
    },
    {
        "group": "shared_event_geometry",
        "name": "translation_action_exponent",
        "treatment": "calibratable_shared",
        "nominal": 1.0,
        "lower_bound": 0.25,
        "upper_bound": 2.0,
        "reason": (
            "Distributes one accepted event advance over its autonomous "
            "cleavage-action interval; one value is shared by all cases."
        ),
    },
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--candidate-registry", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument(
        "--reference-candidate",
        default=REFERENCE_CANDIDATE_ID,
    )
    parser.add_argument(
        "--reference-temperature-K",
        type=int,
        default=REFERENCE_TEMPERATURE_K,
    )
    return parser.parse_args()


def _unique_member(
    archive: zipfile.ZipFile,
    *,
    suffix: str,
    contains: str | None = None,
) -> str:
    matches = [
        name
        for name in archive.namelist()
        if name.endswith(suffix) and (contains is None or contains in name)
    ]
    if len(matches) != 1:
        raise RuntimeError(
            f"expected one archive member ending in {suffix!r}, found {len(matches)}"
        )
    return matches[0]


def _read_csv(
    archive: zipfile.ZipFile,
    member: str,
) -> list[dict[str, str]]:
    with archive.open(member) as raw:
        return list(csv.DictReader(io.TextIOWrapper(raw, encoding="utf-8")))


def _read_json(archive: zipfile.ZipFile, member: str) -> Any:
    with archive.open(member) as raw:
        return json.load(io.TextIOWrapper(raw, encoding="utf-8"))


def _event_step_rows(rows: Iterable[Mapping[str, str]]) -> list[dict[str, float]]:
    events: list[dict[str, float]] = []
    previous_extension = 0.0
    for row in rows:
        extension = float(row["crack_extension_m"])
        if extension > previous_extension + 1.0e-15:
            events.append(
                {
                    "step": float(row["step"]),
                    "U_m": float(row["Uapp_m"]),
                    "K_MPa_sqrt_m": float(row["KJ_Pa_sqrtm"]) * 1.0e-6,
                    "projected_extension_m": extension,
                    "projected_advance_m": extension - previous_extension,
                    "projected_advance_reported_m": float(row["da_block_m"]),
                }
            )
        previous_extension = max(previous_extension, extension)
    if not events:
        raise RuntimeError("steps file contains no accepted crack events")
    return events


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


def _registry_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise RuntimeError(f"empty candidate registry: {path}")
    for row in rows:
        effective_candidate_parameters(row)
    return rows


def _candidate_fingerprint(rows: list[Mapping[str, str]]) -> str:
    return candidate_parameter_fingerprint(rows)


def _summary_member(archive: zipfile.ZipFile) -> str:
    return _unique_member(
        archive,
        suffix="/v10_2_22_dbtt_50um_screen_summary.csv",
    )


def _case_member(
    archive: zipfile.ZipFile,
    option_key: str,
    temperature_K: int,
    filename: str,
) -> str:
    token = f"/{option_key}/T{temperature_K:d}K_th45_seed3621/{filename}"
    return _unique_member(archive, suffix=token)


def _loading_map(
    archive: zipfile.ZipFile,
    summary_rows: list[Mapping[str, str]],
    reference_candidate: str,
    reference_temperature_K: int,
) -> tuple[RCurveLoadingMap, str]:
    reference_rows = [
        row
        for row in summary_rows
        if row["candidate_id"] == reference_candidate
        and int(round(float(row["temperature_K"]))) == reference_temperature_K
    ]
    if len(reference_rows) != 1:
        raise RuntimeError(
            "the loading-map reference must identify exactly one summary row"
        )
    option = reference_rows[0]["option_key"]
    case_token = f"/{option}/T{reference_temperature_K:d}K_th45_seed3621/"
    steps_member = _unique_member(
        archive,
        suffix=f"{case_token}steps_{reference_temperature_K:04d}K.csv",
    )
    event_member = _unique_member(
        archive,
        suffix=f"{case_token}stochastic_avalanche_geometry_events.json",
    )
    args_member = _unique_member(
        archive,
        suffix=f"{case_token}run_args.json",
    )
    step_events = _event_step_rows(_read_csv(archive, steps_member))
    geometry_events = _read_json(archive, event_member)
    run_args = _read_json(archive, args_member)
    if len(step_events) != len(geometry_events):
        raise RuntimeError(
            "reference steps and stochastic geometry event counts disagree"
        )

    K_per_U: list[float] = []
    thresholds: list[float] = []
    path_advances: list[float] = []
    projected_advances: list[float] = []
    reference_K: list[float] = []
    reference_U: list[float] = []
    for index, (step, event) in enumerate(zip(step_events, geometry_events)):
        if int(event["event_index"]) != index:
            raise RuntimeError("stochastic event indices are not contiguous")
        U = float(step["U_m"])
        K = float(step["K_MPa_sqrt_m"])
        if U <= 0.0:
            raise RuntimeError("reference event displacement must be positive")
        path = float(event["event_advance_m"])
        projected = float(event["x1"]) - float(event["x0"])
        if not np.isclose(
            projected,
            float(step["projected_advance_reported_m"]),
            rtol=1.0e-10,
            atol=1.0e-15,
        ):
            raise RuntimeError(
                f"reference event {index} projected advance disagrees with steps file"
            )
        K_per_U.append(K / U)
        thresholds.append(float(event["threshold_action"]))
        path_advances.append(path)
        projected_advances.append(projected)
        reference_K.append(K)
        reference_U.append(U)

    loading = RCurveLoadingMap(
        K_per_U_MPa_sqrt_m_per_m=tuple(K_per_U),
        threshold_actions=tuple(thresholds),
        path_advances_m=tuple(path_advances),
        projected_advances_m=tuple(projected_advances),
        nominal_dU_m=float(run_args["dU"]),
        nominal_dt_s=float(run_args["dt"]),
        seed=int(geometry_events[0]["hazard_seed"]),
        reference_candidate_id=reference_candidate,
        reference_temperature_K=float(reference_temperature_K),
        reference_event_K_MPa_sqrt_m=tuple(reference_K),
        reference_event_U_m=tuple(reference_U),
        provenance={
            "archive_member_steps": steps_member,
            "archive_member_events": event_member,
            "archive_member_run_args": args_member,
            "geometry_map": "event_preceding_KJ_over_applied_displacement",
            "state_translation": "path_advance",
            "R_curve_abscissa": "projected_x_advance",
        },
    )
    loading.validate()
    return loading, option


def main() -> int:
    args = parse_args()
    registry = _registry_rows(args.candidate_registry)
    registry_ids = {row["candidate_id"] for row in registry}
    args.out.mkdir(parents=True, exist_ok=True)

    with zipfile.ZipFile(args.archive) as archive:
        summary_member = _summary_member(archive)
        summary = _read_csv(archive, summary_member)
        summary_ids = {row["candidate_id"] for row in summary}
        if summary_ids != registry_ids:
            raise RuntimeError(
                "archive and candidate registry IDs differ: "
                f"archive_only={sorted(summary_ids - registry_ids)}, "
                f"registry_only={sorted(registry_ids - summary_ids)}"
            )
        if len(summary) != len(registry) * 10:
            raise RuntimeError(
                f"expected {len(registry) * 10} target cases, found {len(summary)}"
            )

        loading_map, reference_option = _loading_map(
            archive,
            summary,
            args.reference_candidate,
            args.reference_temperature_K,
        )
        event_targets: list[dict[str, Any]] = []
        for row in summary:
            option = row["option_key"]
            temperature = int(round(float(row["temperature_K"])))
            member = _case_member(
                archive,
                option,
                temperature,
                f"steps_{temperature:04d}K.csv",
            )
            events = _event_step_rows(_read_csv(archive, member))
            geometry_member = _case_member(
                archive,
                option,
                temperature,
                "stochastic_avalanche_geometry_events.json",
            )
            geometry_events = _read_json(archive, geometry_member)
            if len(events) != loading_map.n_events:
                raise RuntimeError(
                    f"{option} T={temperature} K has {len(events)} events; "
                    f"expected {loading_map.n_events}"
                )
            if len(geometry_events) != loading_map.n_events:
                raise RuntimeError(
                    f"{option} T={temperature} K has "
                    f"{len(geometry_events)} geometry events; "
                    f"expected {loading_map.n_events}"
                )
            for index, (event, geometry_event) in enumerate(
                zip(events, geometry_events)
            ):
                if int(geometry_event["event_index"]) != index:
                    raise RuntimeError(
                        f"{option} T={temperature} K has noncontiguous "
                        "geometry-event indices"
                    )
                common_values = (
                    (
                        float(geometry_event["threshold_action"]),
                        loading_map.threshold_actions[index],
                        "cleavage threshold",
                    ),
                    (
                        float(geometry_event["event_advance_m"]),
                        loading_map.path_advances_m[index],
                        "path advance",
                    ),
                    (
                        float(geometry_event["x1"]) - float(geometry_event["x0"]),
                        loading_map.projected_advances_m[index],
                        "projected advance",
                    ),
                )
                for actual, expected, label in common_values:
                    if not np.isclose(
                        actual,
                        expected,
                        rtol=1.0e-12,
                        atol=1.0e-15,
                    ):
                        raise RuntimeError(
                            f"{option} T={temperature} K event {index} "
                            f"{label} differs from the CRN loading map"
                        )
                event_targets.append(
                    {
                        "option_key": option,
                        "candidate_id": row["candidate_id"],
                        "role": row["role"],
                        "temperature_K": temperature,
                        "seed": int(row["seed"]),
                        "event_index": index,
                        "threshold_action": float(geometry_event["threshold_action"]),
                        "path_advance_m": float(geometry_event["event_advance_m"]),
                        **event,
                    }
                )

    target_fields = (
        "option_key",
        "candidate_id",
        "role",
        "temperature_K",
        "seed",
        "n_events",
        "achieved_extension_um",
        "K_first_MPa_sqrt_m",
        "K_10um_MPa_sqrt_m",
        "K_25um_MPa_sqrt_m",
        "K_50um_MPa_sqrt_m",
        "max_backstress_GPa",
        "max_shield_MPa_sqrt_m",
        "min_front_width_um",
        "max_tip_radius_um",
    )
    target_rows = [
        {field: row[field] for field in target_fields}
        for row in sorted(
            summary,
            key=lambda item: (
                item["candidate_id"],
                float(item["temperature_K"]),
            ),
        )
    ]

    _write_csv(args.out / "v10_2_22_rcurve_checkpoint_targets.csv", target_rows)
    _write_csv(args.out / "v10_2_22_rcurve_event_targets.csv", event_targets)
    _write_csv(
        args.out / "v9_13_rcurve_calibration_parameter_boundary.csv",
        list(CALIBRATION_BOUNDARY_ROWS),
    )
    (args.out / "v10_2_22_rcurve_loading_map.json").write_text(
        json.dumps(loading_map.as_dict(), indent=2)
    )
    manifest = {
        "schema": "v9.13_v10.2.22_autonomous_R_curve_targets",
        "source_archive": args.archive.name,
        "source_summary_member": summary_member,
        "candidate_registry": str(args.candidate_registry),
        "candidate_registry_sha256": _candidate_fingerprint(registry),
        "candidate_parameters_fixed": True,
        "candidate_parameter_fields": list(CANDIDATE_PARAMETER_FIELDS),
        "n_candidates": len(registry),
        "n_temperatures": 10,
        "n_cases": len(target_rows),
        "n_events_per_case": loading_map.n_events,
        "reference_candidate": args.reference_candidate,
        "reference_option": reference_option,
        "reference_temperature_K": args.reference_temperature_K,
        "acceptance_outputs": [
            "K_first_MPa_sqrt_m",
            "K_10um_MPa_sqrt_m",
            "K_25um_MPa_sqrt_m",
            "K_50um_MPa_sqrt_m",
            "event_resolved_K_MPa_sqrt_m",
            "max_backstress_GPa",
            "min_front_width_um",
            "max_tip_radius_um",
        ],
    }
    (args.out / "v9_13_v10_2_22_rcurve_target_manifest.json").write_text(
        json.dumps(manifest, indent=2)
    )
    print(
        "V913_RCURVE_TARGETS_COMPLETE "
        f"cases={len(target_rows)} events={len(event_targets)} "
        f"candidate_sha256={manifest['candidate_registry_sha256']} "
        f"out={args.out}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
