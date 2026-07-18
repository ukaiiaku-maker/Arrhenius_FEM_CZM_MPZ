"""Reporting-only normalization for the v10.0.5 reduced 2-D plastic channels.

The v10.0.5 mechanics and kinetics use two in-plane trace-resolved emission
channels. This module does not reinterpret them as a complete three-dimensional
BCC crystal-plasticity system. It reads completed v10.0.5 outputs, writes
terminology-correct long-form channel diagnostics, and proves that the source
outputs were not modified.

Per-channel diagnostics are written by the progressive formatter into
``kinetic_campaign_czm_progressive_2d_v10_0_3.json``. The legacy
``steps_*K.csv`` files contain outer-step summaries and are retained only as
immutable scalar compatibility outputs. CSV channel columns are accepted as a
fallback for older synthetic/test outputs, but are not required.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from pathlib import Path
from typing import Any, Iterable

SCHEMA = "reduced_2d_slip_trace_reporting_v10_0_5_1"
POINT_RELEASE = "10.0.5.1"
UPSTREAM_AUDIT = "parallel_opening_emission_v10_0_5_audit.json"
UPSTREAM_RESULTS = "mode_i_v10_0_5_results.json"
UPSTREAM_PROGRESSIVE = "kinetic_campaign_czm_progressive_2d_v10_0_3.json"
NORMALIZED_AUDIT = "slip_trace_reporting_v10_0_5_1.json"
NORMALIZED_RESULTS = "mode_i_v10_0_5_1_results.json"
CHANNEL_TABLE = "slip_trace_channels_v10_0_5_1.csv"

_CHANNEL_FIELDS = {
    "drive_factor": "slip_drive_factor_{i}",
    "tau_signed_Pa": "slip_tau_signed_Pa_{i}",
    "sigma_emission_effective_Pa": "sigma_emission_effective_Pa_{i}",
    "sigma_emission_backstress_Pa": "sigma_emission_backstress_Pa_{i}",
    "lambda_emit_s-1": "lambda_emit_s-1_{i}",
    "dN_emit": "dN_emit_{i}",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_json(path: Path) -> Any:
    return json.loads(path.read_text())


def _float_or_none(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _first_finite(*values: Any) -> float | None:
    for value in values:
        number = _float_or_none(value)
        if number is not None:
            return number
    return None


def _temperature_from_path(path: Path) -> float | None:
    match = re.search(r"steps_(\d+)K\.csv$", path.name)
    return float(match.group(1)) if match else None


def _channel_count(fieldnames: Iterable[str]) -> int:
    indices: set[int] = set()
    for name in fieldnames:
        match = re.fullmatch(r"slip_drive_factor_(\d+)", str(name))
        if match:
            indices.add(int(match.group(1)))
    if not indices:
        return 0
    expected = set(range(max(indices) + 1))
    if indices != expected:
        raise RuntimeError(
            f"non-contiguous reduced slip-trace channel indices: {sorted(indices)}"
        )
    return len(indices)


def _channel_names(audit: dict[str, Any], count: int) -> list[str]:
    latest = dict(audit.get("latest") or {})
    legacy = list(latest.get("slip_system_names") or [])
    names = []
    for index in range(count):
        label = str(legacy[index]) if index < len(legacy) else f"trace_{index}"
        names.append(f"2D_slip_trace_channel_{index}:{label}")
    return names


def _progressive_records(path: Path) -> list[dict[str, Any]]:
    if not path.is_file():
        return []
    payload = _read_json(path)
    if not isinstance(payload, dict):
        raise RuntimeError(f"{path.name} must contain a JSON object")
    records = payload.get("records")
    if records is None:
        nested = payload.get("v1002_runtime_after_run")
        if isinstance(nested, dict):
            records = nested.get("records")
    if records is None:
        return []
    if not isinstance(records, list) or any(not isinstance(row, dict) for row in records):
        raise RuntimeError(f"{path.name} records must be a list of objects")
    return [dict(row) for row in records]


def _csv_channel_records(step_paths: list[Path]) -> tuple[list[dict[str, Any]], int]:
    """Compatibility fallback for files that actually contain channel columns."""
    records: list[dict[str, Any]] = []
    channel_count: int | None = None
    for step_path in step_paths:
        with step_path.open(newline="") as handle:
            reader = csv.DictReader(handle)
            fields = list(reader.fieldnames or [])
            count = _channel_count(fields)
            if count < 1:
                continue
            if channel_count is None:
                channel_count = count
            elif count != channel_count:
                raise RuntimeError(
                    f"channel count changed across step files: {channel_count} -> {count}"
                )
            fallback_temperature = _temperature_from_path(step_path)
            for row_number, source in enumerate(reader):
                row = dict(source)
                row.setdefault("temperature_K", fallback_temperature)
                row.setdefault("step", row_number)
                row["_source_record_file"] = step_path.name
                row["_source_record_index"] = row_number
                records.append(row)
    return records, int(channel_count or 0)


def _record_channel_rows(
    records: list[dict[str, Any]],
    *,
    source_name: str,
    upstream_audit: dict[str, Any],
    expected_count: int | None = None,
) -> tuple[list[dict[str, Any]], int, list[str], int, bool]:
    if not records:
        return [], 0, [], 0, False

    first_count = _channel_count(records[0].keys())
    if first_count < 1:
        raise RuntimeError(
            f"{source_name} lacks v10.0.5 per-channel drive diagnostics"
        )
    if expected_count is not None and first_count != expected_count:
        raise RuntimeError(
            f"channel count mismatch: expected {expected_count}, found {first_count}"
        )
    names = _channel_names(upstream_audit, first_count)
    long_rows: list[dict[str, Any]] = []
    finite_drive_rows = 0
    emission_observed = False

    for record_index, row in enumerate(records):
        count = _channel_count(row.keys())
        if count != first_count:
            raise RuntimeError(
                f"channel count changed within {source_name}: {first_count} -> {count} "
                f"at record {record_index}"
            )
        source_file = str(row.get("_source_record_file") or source_name)
        source_index = int(
            _first_finite(row.get("_source_record_index"), record_index) or record_index
        )
        common = {
            "source_record_file": source_file,
            "source_record_index": source_index,
            "source_step_file": (
                source_file if source_file.startswith("steps_") else None
            ),
            "temperature_K": _first_finite(
                row.get("temperature_K"), row.get("T_K"), row.get("T")
            ),
            "accepted_step_index": _first_finite(
                row.get("step"), row.get("accepted_step"), record_index
            ),
            "carry_sequence_index": _float_or_none(row.get("carry_sequence_index")),
            "trial_event_id": _float_or_none(row.get("trial_event_id")),
            "KJ_Pa_sqrt_m": _first_finite(
                row.get("anisotropic_KJ_Pa_sqrt_m"),
                row.get("KJ_Pa_sqrt_m"),
                row.get("K_open_Pa_sqrt_m"),
            ),
            "crack_extension_m": _first_finite(
                row.get("checkpoint_committed_total_m"),
                row.get("crack_extension_m"),
                row.get("micro_advance_total_m"),
            ),
            "cleavage_clock_B": _first_finite(
                row.get("cleavage_clock_B"), row.get("B")
            ),
            "cumulative_emitted_count": _float_or_none(row.get("N_em")),
        }
        for index in range(first_count):
            values = {
                output_name: _float_or_none(row.get(template.format(i=index)))
                for output_name, template in _CHANNEL_FIELDS.items()
            }
            if values["drive_factor"] is not None:
                finite_drive_rows += 1
            if abs(values["dN_emit"] or 0.0) > 0.0:
                emission_observed = True
            long_rows.append({
                **common,
                "slip_trace_channel_index": index,
                "slip_trace_channel_name": names[index],
                **values,
            })
    return long_rows, first_count, names, finite_drive_rows, emission_observed


def normalize_output(root: str | Path) -> dict[str, Any]:
    """Normalize one completed v10.0.5 directory without recomputing physics."""

    root = Path(root)
    audit_path = root / UPSTREAM_AUDIT
    results_path = root / UPSTREAM_RESULTS
    progressive_path = root / UPSTREAM_PROGRESSIVE
    step_paths = sorted(root.glob("steps_*K.csv"))
    required = [audit_path, results_path]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing v10.0.5 source outputs: " + ", ".join(missing))
    if not step_paths:
        raise FileNotFoundError(f"no steps_*K.csv files found in {root}")

    source_paths = [audit_path, results_path, *step_paths]
    if progressive_path.is_file():
        source_paths.append(progressive_path)
    hashes_before = {path.name: _sha256(path) for path in source_paths}
    upstream_audit = dict(_read_json(audit_path))
    upstream_results = _read_json(results_path)
    if not isinstance(upstream_results, list):
        raise RuntimeError("mode_i_v10_0_5_results.json must contain a list")
    if upstream_audit.get("implementation_certified") is not True:
        raise RuntimeError("upstream v10.0.5 implementation audit is not certified")
    if int(upstream_audit.get("capture_count", 0)) < 1:
        raise RuntimeError("upstream v10.0.5 audit contains no tensor captures")
    for key in (
        "drive_factor_normalization_or_clipping_active",
        "directional_multiplier_applied_after_hazard",
        "fit_derived_shielding_cap_active",
    ):
        if upstream_audit.get(key) is not False:
            raise RuntimeError(f"upstream v10.0.5 audit has unexpected {key}")

    records = _progressive_records(progressive_path)
    diagnostic_source = "progressive_runtime_records"
    source_name = progressive_path.name
    if records:
        long_rows, channel_count, names, finite_drive_rows, emission_observed = (
            _record_channel_rows(
                records,
                source_name=source_name,
                upstream_audit=upstream_audit,
            )
        )
    else:
        records, csv_count = _csv_channel_records(step_paths)
        diagnostic_source = "legacy_step_csv_channel_columns"
        source_name = "steps_*K.csv"
        if not records or csv_count < 1:
            raise RuntimeError(
                f"{UPSTREAM_PROGRESSIVE} and steps_*K.csv lack v10.0.5 "
                "per-channel drive diagnostics"
            )
        long_rows, channel_count, names, finite_drive_rows, emission_observed = (
            _record_channel_rows(
                records,
                source_name=source_name,
                upstream_audit=upstream_audit,
                expected_count=csv_count,
            )
        )

    if not emission_observed:
        emission_observed = any(
            abs(_first_finite(row.get("N_em_final"), row.get("max_N_em")) or 0.0) > 0.0
            for row in upstream_results
            if isinstance(row, dict)
        )

    table_path = root / CHANNEL_TABLE
    table_fields = [
        "source_record_file",
        "source_record_index",
        "source_step_file",
        "temperature_K",
        "accepted_step_index",
        "carry_sequence_index",
        "trial_event_id",
        "KJ_Pa_sqrt_m",
        "crack_extension_m",
        "cleavage_clock_B",
        "cumulative_emitted_count",
        "slip_trace_channel_index",
        "slip_trace_channel_name",
        *_CHANNEL_FIELDS.keys(),
    ]
    with table_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=table_fields)
        writer.writeheader()
        writer.writerows(long_rows)

    normalized_results = []
    for source_row in upstream_results:
        row = dict(source_row)
        row.update({
            "point_release": POINT_RELEASE,
            "reporting_schema": SCHEMA,
            "plastic_channel_representation": "reduced_2d_slip_trace_channels",
            "slip_trace_channel_count": channel_count,
            "slip_trace_channel_names": names,
            "full_3d_bcc_slip_system_model_active": False,
            "legacy_slip_system_field_names_preserved": True,
            "emission_observation_required_for_implementation_certification": False,
            "emission_observed_in_this_run": emission_observed,
            "zero_emission_is_valid_implementation_outcome": True,
            "material_parameterization_assessment_performed": False,
            "response_classification_gate_active": False,
            "slip_trace_diagnostic_source": diagnostic_source,
            "slip_trace_channel_table": CHANNEL_TABLE,
        })
        normalized_results.append(row)
    normalized_results_path = root / NORMALIZED_RESULTS
    normalized_results_path.write_text(json.dumps(normalized_results, indent=2, default=str))

    expected_long_rows = len(records) * channel_count
    certified = bool(
        upstream_audit.get("implementation_certified") is True
        and len(long_rows) == expected_long_rows
        and len(long_rows) > 0
        and finite_drive_rows == len(long_rows)
    )
    normalized_audit = {
        "schema": SCHEMA,
        "point_release": POINT_RELEASE,
        "upstream_point_release": "10.0.5",
        "reporting_only": True,
        "physics_recomputed": False,
        "source_outputs_modified": False,
        "plastic_channel_representation": "reduced_2d_slip_trace_channels",
        "slip_trace_channel_count": channel_count,
        "slip_trace_channel_names": names,
        "full_3d_bcc_slip_system_model_active": False,
        "legacy_slip_system_field_names_preserved": True,
        "emission_observation_required_for_implementation_certification": False,
        "emission_observed_in_this_run": emission_observed,
        "zero_emission_is_valid_implementation_outcome": True,
        "material_parameterization_assessment_performed": False,
        "response_classification_gate_active": False,
        "slip_trace_diagnostic_source": diagnostic_source,
        "source_record_count": len(records),
        "channel_rows_written": len(long_rows),
        "finite_drive_rows": finite_drive_rows,
        "upstream_tensor_capture_count": int(upstream_audit.get("capture_count", 0)),
        "upstream_nonzero_drive_capture_count": int(
            upstream_audit.get("nonzero_emission_drive_capture_count", 0)
        ),
        "source_sha256_before": hashes_before,
        "channel_table": CHANNEL_TABLE,
        "normalized_results": NORMALIZED_RESULTS,
        "implementation_certified": certified,
    }
    normalized_audit_path = root / NORMALIZED_AUDIT
    normalized_audit_path.write_text(json.dumps(normalized_audit, indent=2, default=str))

    hashes_after = {path.name: _sha256(path) for path in source_paths}
    if hashes_after != hashes_before:
        raise RuntimeError("v10.0.5 source outputs changed during reporting normalization")
    normalized_audit["source_sha256_after"] = hashes_after
    normalized_audit_path.write_text(json.dumps(normalized_audit, indent=2, default=str))
    if not certified:
        raise RuntimeError(
            "v10.0.5.1 reporting certification failed: "
            + json.dumps(normalized_audit, default=str)
        )
    return normalized_audit


__all__ = [
    "SCHEMA",
    "POINT_RELEASE",
    "UPSTREAM_PROGRESSIVE",
    "NORMALIZED_AUDIT",
    "NORMALIZED_RESULTS",
    "CHANNEL_TABLE",
    "normalize_output",
]
