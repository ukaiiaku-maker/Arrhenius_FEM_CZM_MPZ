"""Reporting-only normalization for the v10.0.5 reduced 2-D plastic channels.

The v10.0.5 mechanics and kinetics use two in-plane trace-resolved emission
channels.  This module does not reinterpret them as a complete three-dimensional
BCC crystal-plasticity system.  It reads completed v10.0.5 outputs, writes
terminology-correct long-form channel diagnostics, and proves that the source
outputs were not modified.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from pathlib import Path
from typing import Any

SCHEMA = "reduced_2d_slip_trace_reporting_v10_0_5_1"
POINT_RELEASE = "10.0.5.1"
UPSTREAM_AUDIT = "parallel_opening_emission_v10_0_5_audit.json"
UPSTREAM_RESULTS = "mode_i_v10_0_5_results.json"
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


def _temperature_from_path(path: Path) -> float | None:
    match = re.search(r"steps_(\d+)K\.csv$", path.name)
    return float(match.group(1)) if match else None


def _channel_count(fieldnames: list[str]) -> int:
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


def normalize_output(root: str | Path) -> dict[str, Any]:
    """Normalize one completed v10.0.5 directory without recomputing physics."""

    root = Path(root)
    audit_path = root / UPSTREAM_AUDIT
    results_path = root / UPSTREAM_RESULTS
    step_paths = sorted(root.glob("steps_*K.csv"))
    required = [audit_path, results_path]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing v10.0.5 source outputs: " + ", ".join(missing))
    if not step_paths:
        raise FileNotFoundError(f"no steps_*K.csv files found in {root}")

    source_paths = [audit_path, results_path, *step_paths]
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

    long_rows: list[dict[str, Any]] = []
    channel_count: int | None = None
    names: list[str] = []
    finite_drive_rows = 0
    emission_observed = False

    for step_path in step_paths:
        with step_path.open(newline="") as handle:
            reader = csv.DictReader(handle)
            fields = list(reader.fieldnames or [])
            count = _channel_count(fields)
            if count < 1:
                raise RuntimeError(
                    f"{step_path.name} lacks v10.0.5 per-channel drive columns"
                )
            if channel_count is None:
                channel_count = count
                names = _channel_names(upstream_audit, count)
            elif count != channel_count:
                raise RuntimeError(
                    f"channel count changed across step files: {channel_count} -> {count}"
                )

            fallback_temperature = _temperature_from_path(step_path)
            for row_number, row in enumerate(reader):
                temperature = (
                    _float_or_none(row.get("T_K"))
                    or _float_or_none(row.get("T"))
                    or fallback_temperature
                )
                accepted_step = (
                    _float_or_none(row.get("step"))
                    or _float_or_none(row.get("accepted_step"))
                    or float(row_number)
                )
                common = {
                    "source_step_file": step_path.name,
                    "temperature_K": temperature,
                    "accepted_step_index": accepted_step,
                    "KJ_Pa_sqrt_m": (
                        _float_or_none(row.get("anisotropic_KJ_Pa_sqrt_m"))
                        or _float_or_none(row.get("KJ_Pa_sqrt_m"))
                    ),
                    "crack_extension_m": _float_or_none(row.get("crack_extension_m")),
                    "cleavage_clock_B": (
                        _float_or_none(row.get("cleavage_clock_B"))
                        if row.get("cleavage_clock_B") not in (None, "")
                        else _float_or_none(row.get("B"))
                    ),
                    "cumulative_emitted_count": _float_or_none(row.get("N_em")),
                }
                for index in range(count):
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

    assert channel_count is not None
    table_path = root / CHANNEL_TABLE
    table_fields = [
        "source_step_file",
        "temperature_K",
        "accepted_step_index",
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
            "slip_trace_channel_table": CHANNEL_TABLE,
        })
        normalized_results.append(row)
    normalized_results_path = root / NORMALIZED_RESULTS
    normalized_results_path.write_text(json.dumps(normalized_results, indent=2, default=str))

    expected_long_rows = sum(
        max(sum(1 for _ in csv.DictReader(path.open(newline=""))), 0)
        for path in step_paths
    ) * channel_count
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
    "NORMALIZED_AUDIT",
    "NORMALIZED_RESULTS",
    "CHANNEL_TABLE",
    "normalize_output",
]
