#!/usr/bin/env python3
"""Validate focused Peak/DBTT R-propensity survivors across seeds and temperature.

This is a reduced-model calculation only.  It does not invoke either 2-D
production trajectory driver.  The 300 um continuations are deliberately
limited to the canonical seed and transition temperatures.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_outputs" / "oneD_v2_peak_dbtt_rcurve_search"
sys.path.insert(0, str(ROOT))

from arrhenius_fracture.emergent_gnd_contract_v913 import ACTIVE_CANDIDATE_PARAMETER_FIELDS
from reduced_fracture_v2.rcurve_propensity import summarize_rcurve_propensity
from scripts.run_oneD_v2_peak_dbtt_R_screen import (
    CLASS_SEEDS, CONTROL_IDS, _base_record,
)
from scripts.run_oneD_v2_predictive_campaign import inputs, run_case, summary


TEMPERATURES = (300.0, 600.0, 800.0, 900.0, 1000.0, 1100.0, 1200.0)
SEEDS = {
    "Peak": (8666, 8667, 8668),
    "DBTT": (1008666, 1008667, 1008668),
}
SURVIVORS = {
    "Peak": (
        "v913_zeroD_sobol_0242980",
        "oneD_v2_peak_R_41f8789bcbc1f097",
        "oneD_v2_peak_R_a2e923a7587543db",
        "oneD_v2_peak_R_7834115ae79cd559",
    ),
    "DBTT": (
        "v913_zeroD_sobol_0202500",
        "oneD_v2_dbtt_R_9f5160f509e713e2",
        "oneD_v2_dbtt_R_c4859a34963f15af",
        "oneD_v2_dbtt_R_7b3b5b45354e64ef",
        "oneD_v2_dbtt_R_3e168d9381eafa7b",
        "oneD_v2_dbtt_R_41a3f2096b3c4c54",
        "oneD_v2_dbtt_R_c2dabb42101cf812",
    ),
}
LONG_TEMPERATURES = {
    "Peak": (900.0, 1000.0),
    "DBTT": (900.0, 1000.0, 1100.0),
}


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _candidate_rows(material: str) -> pd.DataFrame:
    slug = material.lower()
    path = OUT / f"oneD_v2_{slug}_R_candidates.csv"
    rows = pd.read_csv(path)
    selected = rows[rows.candidate_id.isin(SURVIVORS[material])].copy()
    missing = sorted(set(SURVIVORS[material]) - set(selected.candidate_id.astype(str)))
    if missing:
        raise ValueError(f"missing {material} survivors: {missing}")
    return selected.set_index("candidate_id", drop=False).loc[list(SURVIVORS[material])].reset_index(drop=True)


def _record(row: pd.Series, material: str, provider: str, temperature: float,
            seed: int, target_um: float, mechanics, drive, physics) -> dict:
    record = _base_record(row, material, provider, temperature)
    record.update({
        "hazard_seed": int(seed),
        "target_um": float(target_um),
        "search_stage": "MULTISEED_DENSE_TEMPERATURE_VALIDATION"
        if target_um == 100.0 else "CANONICAL_SEED_300UM_CONTINUATION",
    })
    try:
        result = run_case(
            row, material, temperature, provider, mechanics, drive, physics,
            target_um, seed=seed, maximum_intervals=30_000,
        )
        record.update(summary(result))
        record.update(summarize_rcurve_propensity(result))
        onsets = json.loads(record["onset_candidates_json"])
        record.update({
            "maximum_onset_tip_radius_um": max(
                (x["tip_radius_m"] * 1.0e6 for x in onsets), default=np.nan
            ),
            "maximum_onset_backstress_GPa": max(
                (x["backstress_Pa"] * 1.0e-9 for x in onsets), default=np.nan
            ),
            "maximum_onset_mobile_density_m2": max(
                (x["mobile_density_m2"] for x in onsets), default=np.nan
            ),
            "maximum_onset_retained_density_m2": max(
                (x["retained_density_m2"] for x in onsets), default=np.nan
            ),
            "maximum_onset_source_multiplicity": max(
                (x["source_multiplicity"] for x in onsets), default=np.nan
            ),
            "map_oracle_fallback_count": int(
                result["status"] == "RIGHT_CENSORED_DRIVE_MAP_BOUND"
            ),
            "terminal_reason": result["status"],
            "case_error": "",
        })
    except Exception as exc:
        record.update({
            "status": "EVALUATION_FAILURE",
            "terminal_reason": type(exc).__name__,
            "case_error": str(exc),
            "N_reinit": 0,
            "physical_avalanche_count": 0,
            "largest_avalanche_fraction": 0.0,
            "map_oracle_fallback_count": 0,
        })
    return record


def planned_cases() -> list[tuple[str, str, float, int, float]]:
    cases = []
    for material, candidates in SURVIVORS.items():
        for candidate in candidates:
            for provider in ("PF", "FEMCZM"):
                for temperature in TEMPERATURES:
                    for seed in SEEDS[material]:
                        cases.append((material, candidate, temperature, seed, 100.0, provider))
                if candidate != CONTROL_IDS[material]:
                    for temperature in LONG_TEMPERATURES[material]:
                        cases.append((
                            material, candidate, temperature,
                            CLASS_SEEDS[material], 300.0, provider,
                        ))
    return cases


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--restart", action="store_true", help="discard an existing checkpoint")
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    checkpoint = OUT / "oneD_v2_peak_dbtt_R_multiseed_checkpoint.parquet"
    if args.restart:
        checkpoint.unlink(missing_ok=True)
    output = OUT / "oneD_v2_peak_dbtt_R_multiseed_validation.csv"
    if checkpoint.exists():
        records = pd.read_parquet(checkpoint).to_dict("records")
    elif output.exists() and not args.restart:
        records = pd.read_csv(output).to_dict("records")
    else:
        records = []
    done = {
        (str(x["target_response_class"]), str(x["candidate_id"]), str(x["provider"]),
         float(x["temperature_K"]), int(x["hazard_seed"]), float(x["target_um"]))
        for x in records
    }
    physics, _, providers = inputs()
    candidate_lookup = {
        material: _candidate_rows(material).set_index("candidate_id", drop=False)
        for material in SURVIVORS
    }
    cases = planned_cases()
    for index, (material, candidate, temperature, seed, target_um, provider) in enumerate(cases):
        key = (material, candidate, provider, temperature, seed, target_um)
        if key in done:
            continue
        mechanics, drive = providers[provider]
        row = candidate_lookup[material].loc[candidate]
        records.append(_record(
            row, material, provider, temperature, seed, target_um,
            mechanics, drive, physics,
        ))
        done.add(key)
        if len(records) % 12 == 0 or index + 1 == len(cases):
            pd.DataFrame(records).to_parquet(checkpoint, index=False)
            print(f"MULTISEED_PROGRESS cases={len(records)}/{len(cases)}", flush=True)

    frame = pd.DataFrame(records).sort_values(
        ["target_response_class", "candidate_id", "target_um", "provider",
         "temperature_K", "hazard_seed"], kind="stable",
    )
    frame.to_csv(output, index=False)
    checkpoint.unlink(missing_ok=True)
    manifest = {
        "schema": "oneD_v2_peak_dbtt_R_multiseed_validation_v1",
        "survivors": {key: list(value) for key, value in SURVIVORS.items()},
        "temperatures_K": list(TEMPERATURES),
        "seeds": {key: list(value) for key, value in SEEDS.items()},
        "canonical_seed_300um_temperatures_K": {
            key: list(value) for key, value in LONG_TEMPERATURES.items()
        },
        "case_count": int(len(frame)),
        "output_sha256": _sha(output),
        "new_2D_PF_runs": 0,
        "new_2D_FEMCZM_runs": 0,
        "material_fields": list(ACTIVE_CANDIDATE_PARAMETER_FIELDS),
    }
    (OUT / "oneD_v2_peak_dbtt_R_multiseed_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    print(f"MULTISEED_COMPLETE cases={len(frame)} output={output}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
