#!/usr/bin/env python3
"""Compact no-retuning regression for the frozen weak-T and ceramic rows."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_outputs" / "oneD_v2_peak_dbtt_rcurve_search"
SOURCE = ROOT / "analysis_outputs" / "oneD_v2_terminal_predictive_program"
sys.path.insert(0, str(ROOT))

from reduced_fracture_v2.rcurve_propensity import summarize_rcurve_propensity
from scripts.run_oneD_v2_predictive_campaign import inputs, run_case, summary


ROWS = {
    "weak-T": "oneD_v2_focused_weak_T_0016",
    "ceramic-like": "oneD_v2_focused_ceramic_like_0018",
}
SEEDS = {"weak-T": 2008666, "ceramic-like": 3008666}
TEMPERATURES = (300.0, 1000.0, 1200.0)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    registry = pd.read_csv(SOURCE / "oneD_v2_new_four_class_registry.csv")
    reference = pd.read_csv(SOURCE / "oneD_v2_final_four_class_results.csv")
    physics, _, providers = inputs()
    records = []
    for material, candidate_id in ROWS.items():
        row = registry[registry.candidate_id == candidate_id]
        if len(row) != 1:
            raise ValueError(f"missing unique frozen row {candidate_id}")
        row = row.iloc[0]
        for provider, (mechanics, drive) in providers.items():
            for temperature in TEMPERATURES:
                result = run_case(
                    row, material, temperature, provider, mechanics, drive,
                    physics, 100.0, seed=SEEDS[material], maximum_intervals=10_000,
                )
                record = summary(result)
                record.update(summarize_rcurve_propensity(result))
                record.update({
                    "hazard_seed": SEEDS[material],
                    "regression_status": "NEW_TEMPERATURE_POINT",
                    "reference_initial_onset_MPa_sqrt_m": np.nan,
                    "initial_onset_absolute_error_MPa_sqrt_m": np.nan,
                    "reference_physical_avalanche_count": np.nan,
                })
                old = reference[
                    (reference.material_class == material)
                    & (reference.provider == provider)
                    & (reference.temperature_K == temperature)
                ]
                if len(old) == 1:
                    expected = old.iloc[0]
                    error = float(record["initial_onset_native_KJ_MPa_sqrt_m"]) - float(
                        expected.first_event_native_KJ_MPa_sqrt_m
                    )
                    same_topology = int(record["physical_avalanche_count"]) == int(
                        expected.physical_avalanche_count
                    )
                    record.update({
                        "regression_status": "NUMERIC_REPRODUCTION_PASS"
                        if abs(error) <= 1.0e-12 and same_topology else "REGRESSION_FAILURE",
                        "reference_initial_onset_MPa_sqrt_m": float(
                            expected.first_event_native_KJ_MPa_sqrt_m
                        ),
                        "initial_onset_absolute_error_MPa_sqrt_m": error,
                        "reference_physical_avalanche_count": int(
                            expected.physical_avalanche_count
                        ),
                    })
                records.append(record)
    frame = pd.DataFrame(records)
    if (frame.regression_status == "REGRESSION_FAILURE").any():
        print(frame[frame.regression_status == "REGRESSION_FAILURE"][[
            "material_class", "provider", "temperature_K",
            "initial_onset_native_KJ_MPa_sqrt_m",
            "reference_initial_onset_MPa_sqrt_m",
            "initial_onset_absolute_error_MPa_sqrt_m",
            "physical_avalanche_count", "reference_physical_avalanche_count",
        ]].to_string(index=False))
        raise RuntimeError("weak-T/ceramic compact regression changed a frozen response")
    output = OUT / "oneD_v2_weak_ceramic_compact_regression.csv"
    frame.to_csv(output, index=False)
    manifest = {
        "schema": "oneD_v2_weak_ceramic_compact_regression_v1",
        "frozen_rows": ROWS,
        "temperatures_K": list(TEMPERATURES),
        "providers": list(providers),
        "case_count": len(frame),
        "numeric_reproduction_case_count": int(
            frame.regression_status.eq("NUMERIC_REPRODUCTION_PASS").sum()
        ),
        "new_temperature_case_count": int(
            frame.regression_status.eq("NEW_TEMPERATURE_POINT").sum()
        ),
        "output_sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
        "rows_retuned": False,
        "new_2D_PF_runs": 0,
        "new_2D_FEMCZM_runs": 0,
    }
    (OUT / "oneD_v2_weak_ceramic_compact_regression_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    print(frame[["material_class", "provider", "temperature_K", "status",
                 "physical_avalanche_count", "N_reinit", "regression_status"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
