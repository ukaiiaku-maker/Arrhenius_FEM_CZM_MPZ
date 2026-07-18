#!/usr/bin/env python3
"""Fail-closed audit for v10.0.5.2 multicommit long-growth gates."""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path
import sys
from typing import Any


def _load(path: Path) -> Any:
    return json.loads(path.read_text())


def _records(payload: dict[str, Any]) -> list[dict[str, Any]]:
    rows = payload.get("records")
    if rows is None and isinstance(payload.get("v1002_runtime_after_run"), dict):
        rows = payload["v1002_runtime_after_run"].get("records")
    if not isinstance(rows, list):
        raise RuntimeError("progressive runtime contains no record list")
    return [dict(row) for row in rows]


def _channel_indices(row: dict[str, Any]) -> list[int]:
    prefix = "slip_drive_factor_"
    values = sorted({
        int(str(key)[len(prefix):])
        for key in row
        if str(key).startswith(prefix) and str(key)[len(prefix):].isdigit()
    })
    if values != list(range(len(values))):
        raise RuntimeError(f"non-contiguous slip-trace channel indices: {values}")
    return values


def _last_csv_row(path: Path) -> dict[str, str]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise RuntimeError(f"empty step file: {path}")
    return rows[-1]


def main(argv: list[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    if not args:
        raise SystemExit(
            "usage: audit_v10_0_5_2_long_growth.py OUTROOT "
            "[--target-um 100] [--expected-mpz-bins 200] [--da-um 5]"
        )
    root = Path(args[0])

    def option(name: str, default: float) -> float:
        if name in args:
            index = args.index(name)
            if index + 1 >= len(args):
                raise SystemExit(f"missing value after {name}")
            return float(args[index + 1])
        return float(default)

    target_um = option("--target-um", 100.0)
    expected_bins = int(option("--expected-mpz-bins", 200.0))
    da_um = option("--da-um", 5.0)
    expected_commits = int(round(target_um / da_um))

    required = {
        "completion": root / "run_completion_v10_0_5_2.json",
        "channel": root / "parallel_channel_diagnostics_v10_0_5_2.json",
        "progressive": root / "kinetic_campaign_czm_progressive_2d_v10_0_3.json",
        "quality": root / "explicit_quality_wrapper_chain_v91856.json",
        "args": root / "run_args.json",
        "steps": root / "steps_0700K.csv",
        "normalized": root / "slip_trace_reporting_v10_0_5_1.json",
    }
    missing = [str(path) for path in required.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing v10.0.5.2 audit inputs: " + ", ".join(missing))

    completion = dict(_load(required["completion"]))
    channel = dict(_load(required["channel"]))
    progressive = dict(_load(required["progressive"]))
    quality = dict(_load(required["quality"]))
    run_args = dict(_load(required["args"]))
    normalized = dict(_load(required["normalized"]))
    records = _records(progressive)
    final = _last_csv_row(required["steps"])

    final_extension_um = float(final["crack_extension_m"]) * 1.0e6
    committed = int(progressive.get("committed_events", 0))
    insertions = int(progressive.get("trial_insertions", 0))
    rejections = int(progressive.get("damage_rejections", 0))
    rollbacks = int(progressive.get("full_rollbacks", 0))
    mpz_bins = int(run_args.get("mpz_n_bins", 0))
    mpz_length_m = float(run_args.get("mpz_length_m", 0.0))

    channel_count = None
    finite_hazards = 0
    finite_increments = 0
    max_partition_residual = 0.0
    for record_index, row in enumerate(records):
        indices = _channel_indices(row)
        if not indices:
            raise RuntimeError(f"record {record_index} has no slip-trace channels")
        if channel_count is None:
            channel_count = len(indices)
        elif channel_count != len(indices):
            raise RuntimeError("slip-trace channel count changed during the run")
        if row.get("per_channel_emission_partition_verified") is not True:
            raise RuntimeError(f"record {record_index} lacks emission partition verification")
        residual = abs(float(row.get("per_channel_emission_partition_residual", math.inf)))
        if not math.isfinite(residual):
            raise RuntimeError(f"record {record_index} has nonfinite partition residual")
        max_partition_residual = max(max_partition_residual, residual)
        for index in indices:
            hazard = float(row[f"lambda_emit_s-1_{index}"])
            increment = float(row[f"dN_emit_{index}"])
            if not math.isfinite(hazard) or hazard < 0.0:
                raise RuntimeError(f"record {record_index} has invalid channel hazard")
            if not math.isfinite(increment) or increment < -1.0e-14:
                raise RuntimeError(f"record {record_index} has invalid channel increment")
            finite_hazards += 1
            finite_increments += 1

    tolerance_um = max(1.0e-6, 1.0e-9 * max(target_um, 1.0))
    checks = {
        "completion_manifest_complete": (
            completion.get("status") == "complete"
            and completion.get("run_completed_without_exception") is True
        ),
        "quality_audit_complete": quality.get("run_completed_without_exception") is True,
        "quality_vetoes_zero": len(quality.get("quality_vetoes") or []) == 0,
        "target_extension_exact": abs(final_extension_um - target_um) <= tolerance_um,
        "committed_event_count": committed == expected_commits,
        "trial_insertion_count": insertions == expected_commits,
        "damage_rejections_zero": rejections == 0,
        "full_rollbacks_zero": rollbacks == 0,
        "mpz_bin_count": mpz_bins == expected_bins,
        "mpz_length_100um": abs(mpz_length_m - 100.0e-6) <= 1.0e-12,
        "channel_audit_certified": channel.get("implementation_certified") is True,
        "channel_diagnostics_complete": (
            channel.get("per_channel_strang_diagnostics_complete") is True
        ),
        "normalized_reporting_certified": normalized.get("implementation_certified") is True,
        "finite_channel_hazards_complete": (
            finite_hazards == len(records) * int(channel_count or 0)
        ),
        "finite_channel_increments_complete": (
            finite_increments == len(records) * int(channel_count or 0)
        ),
        "no_response_classification_gate": (
            normalized.get("response_classification_gate_active") is False
        ),
    }
    failed = [name for name, passed in checks.items() if not passed]
    payload = {
        "schema": "v10_0_5_2_long_growth_gate_v1",
        "target_extension_um": target_um,
        "physical_increment_um": da_um,
        "expected_commits": expected_commits,
        "final_extension_um": final_extension_um,
        "committed_events": committed,
        "trial_insertions": insertions,
        "accepted_substep_records": len(records),
        "damage_rejections": rejections,
        "full_rollbacks": rollbacks,
        "mpz_n_bins": mpz_bins,
        "mpz_length_m": mpz_length_m,
        "slip_trace_channel_count": int(channel_count or 0),
        "finite_channel_hazards": finite_hazards,
        "finite_channel_increments": finite_increments,
        "maximum_emission_partition_residual": max_partition_residual,
        "constitutive_physics_changed_in_v10052": False,
        "material_parameterization_assessment_performed": False,
        "response_classification_gate_active": False,
        "checks": checks,
        "failed_checks": failed,
        "pass": not failed,
    }
    path = root / "long_growth_gate_v10_0_5_2.json"
    path.write_text(json.dumps(payload, indent=2, default=str))
    print(json.dumps(payload, indent=2, default=str))
    if failed:
        raise RuntimeError("v10.0.5.2 long-growth gate failed: " + ", ".join(failed))
    print(f"V10.0.5.2 LONG-GROWTH GATE PASSED: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
