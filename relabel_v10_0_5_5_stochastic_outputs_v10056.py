#!/usr/bin/env python3
"""Relabel completed v10.0.5.5 stochastic outputs without rerunning physics."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

from arrhenius_fracture.kj_audit_v10056 import (
    POINT_RELEASE,
    enrich_stochastic_block_rows,
)


BLOCK_INPUT = "fatigue_block_diagnostics_v10_0_5_4.csv"
CASE_INPUT = "K_vs_delta_sigma.csv"
BLOCK_OUTPUT = "fatigue_block_diagnostics_v10_0_5_6.csv"
CASE_OUTPUT = "K_vs_delta_sigma_v10_0_5_6.csv"
MANIFEST_OUTPUT = "diagnostic_relabel_v10_0_5_6.json"


def _read_csv(path: Path) -> list[dict[str, Any]]:
    with path.open(newline="") as handle:
        rows = []
        for raw in csv.DictReader(handle):
            row: dict[str, Any] = {}
            for key, value in raw.items():
                try:
                    row[key] = float(value)
                except (TypeError, ValueError):
                    token = str(value).strip().lower()
                    if token in {"true", "false"}:
                        row[key] = token == "true"
                    else:
                        row[key] = value
            rows.append(row)
        return rows


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def relabel(root: Path) -> dict[str, Any]:
    root = root.resolve()
    cases = _read_csv(root / CASE_INPUT)
    blocks = _read_csv(root / BLOCK_INPUT)
    output_cases = []
    output_blocks = []

    for case in cases:
        temperature = float(case["temperature_K"])
        delta_sigma = float(case["delta_sigma_requested_MPa"])
        case_dir = Path(str(case["run_directory"]))
        audit_path = case_dir / "stochastic_vhcf_v10_0_5_5.json"
        audit = json.loads(audit_path.read_text())
        records = list(audit.get("scheduler", {}).get("records", []))
        engines = list(audit.get("engines", []))
        if len(engines) != 1:
            raise RuntimeError(f"expected one engine in {audit_path}")
        engine = dict(engines[0])
        selected = [
            row
            for row in blocks
            if math.isclose(float(row["temperature_K"]), temperature)
            and math.isclose(
                float(row["delta_sigma_requested_MPa"]), delta_sigma
            )
        ]
        output_blocks.extend(enrich_stochastic_block_rows(selected, records))

        updated = dict(case)
        updated.update(
            {
                "diagnostic_schema": "v10.0.5.6",
                "cleavage_threshold": float(engine.get("cleavage_threshold", math.nan)),
                "cleavage_event_index": int(engine.get("cleavage_event_index", 0)),
                "source_budget_total": float(engine.get("source_budget_total", math.nan)),
                "source_consumed_final": float(engine.get("source_budget_consumed", math.nan)),
                "source_remaining_final": float(engine.get("source_budget_remaining", math.nan)),
                "mobile_count_final": float(engine.get("mobile_count", math.nan)),
                "retained_count_final": float(engine.get("retained_count", math.nan)),
                "active_count_final": float(engine.get("active_count", math.nan)),
                "cumulative_emitted": float(engine.get("cumulative_emitted", math.nan)),
                "stochastic_emission_channel_events": float(
                    engine.get("stochastic_emission_events", math.nan)
                ),
                "predictor_mean_field_calls": int(
                    engine.get("predictor_mean_field_calls", 0)
                ),
            }
        )
        output_cases.append(updated)

    _write_csv(root / BLOCK_OUTPUT, output_blocks)
    _write_csv(root / CASE_OUTPUT, output_cases)
    manifest = {
        "schema": "diagnostic_relabel_v10_0_5_6",
        "point_release": POINT_RELEASE,
        "source_release": "10.0.5.5",
        "root": str(root),
        "n_cases": len(output_cases),
        "n_blocks": len(output_blocks),
        "outputs": {"cases": CASE_OUTPUT, "blocks": BLOCK_OUTPUT},
        "simulation_rerun": False,
        "constitutive_physics_changed": False,
    }
    (root / MANIFEST_OUTPUT).write_text(json.dumps(manifest, indent=2))
    return manifest


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("roots", nargs="+", type=Path)
    args = parser.parse_args(argv)
    for root in args.roots:
        manifest = relabel(root)
        print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
