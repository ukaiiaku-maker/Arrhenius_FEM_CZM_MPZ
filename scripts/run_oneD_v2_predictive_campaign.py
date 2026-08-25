#!/usr/bin/env python3
"""Run the provider-native V2 baselines and current four-class matrix."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_outputs" / "oneD_v2_terminal_predictive_program"
DIAGNOSTIC = ROOT / "analysis_outputs" / "oneD_v2_predictive_model"
SOURCE_MAPS = OUT
MAPS = ROOT / "analysis_outputs" / "oneD_v2_mechanics_maps_and_baselines"
PF_DATA = Path("/Volumes/Data/Data/Nanopillar_calculation/PF-fracture-fatigue_v10_2_21_persistent_sites_top1")
REGISTRY = PF_DATA / "arrhenius_fracture/data/materials/v10_2_27_v913_four_class_paper_registry.csv"
sys.path.insert(0, str(ROOT))

from arrhenius_fracture.emergent_gnd_campaign_v913 import candidate_from_registry_row
from arrhenius_fracture.emergent_gnd_types_v913 import CommonPhysics
from reduced_fracture_v2.predictive import (
    ProviderMechanicsMap,
    SourceDriveMap,
    provider_loading_map,
    run_zero_d_predictive,
)

IDS = {
    "v913_zeroD_sobol_0242980": "Peak",
    "v913_zeroD_sobol_0202500": "DBTT",
    "v913_zeroD_sobol_0129902": "weak-T",
    "v913_zeroD_sobol_0077080": "ceramic-like",
}
SEEDS = {"Peak": 8666, "DBTT": 1008666, "weak-T": 2008666, "ceramic-like": 3008666}
TEMPERATURES = (300.0, 600.0, 900.0, 1100.0, 1200.0)


def _jsonable(value):
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    raise TypeError(type(value).__name__)


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def inputs():
    physics = CommonPhysics(**json.loads((ROOT / "mpz_v9_13_v10222_transfer_common_physics.json").read_text())["common_physics"])
    rows = pd.read_csv(REGISTRY)
    rows = rows[rows.candidate_id.isin(IDS)].copy()
    providers = {
        "PF": (
            ProviderMechanicsMap.from_csv("PF", MAPS / "oneD_v2_pf_mechanics_map.csv"),
            SourceDriveMap.from_csv(SOURCE_MAPS / "oneD_v2_pf_source_drive_map.csv"),
        ),
        "FEMCZM": (
            ProviderMechanicsMap.from_csv(
                "FEMCZM", MAPS / "oneD_v2_fem_native_mechanics_map.csv",
                MAPS / "oneD_v2_fem_qualified_G_map.csv",
            ),
            SourceDriveMap.from_csv(SOURCE_MAPS / "oneD_v2_fem_source_drive_map.csv"),
        ),
    }
    return physics, rows, providers


def run_case(
    row, material, temperature, backend, mechanics, drive, physics, target_um,
    *, seed=None,
):
    loading = provider_loading_map(
        mechanics,
        seed=SEEDS[material] if seed is None else int(seed),
        target_extension_m=target_um * 1.0e-6,
    )
    result = run_zero_d_predictive(
        candidate_from_registry_row(row), physics, mechanics, drive, loading,
        temperature, target_extension_m=target_um * 1.0e-6,
    )
    result.update({"material_class": material, "provider": backend, "target_um": target_um})
    return result


def summary(result):
    events = result["events"]
    onsets = [x for x in events if x["event_index"] == 0 or x["reload_separated"]]
    onset_values = [float(x["native_KJ_MPa_sqrt_m"]) for x in onsets]
    qualified = [float(x["qualified_G_J_per_m2"]) for x in onsets if x["qualified_G_J_per_m2"] is not None]
    sizes = [float(x["event_length_m"]) * 1.0e6 for x in events]
    return {
        "provider": result["provider"], "material_class": result["material_class"],
        "candidate_id": result["candidate_id"], "temperature_K": result["temperature_K"],
        "target_um": result["target_um"], "status": result["status"],
        "event_count": result["event_count"],
        "physical_avalanche_count": result["physical_avalanche_count"],
        "precursor_reinitiation_count": result["precursor_reinitiation_count"],
        "largest_avalanche_fraction": result["largest_avalanche_fraction"],
        "first_event_opening_um": None if result["first_event_opening_m"] is None else result["first_event_opening_m"] * 1.0e6,
        "first_event_native_KJ_MPa_sqrt_m": result["first_event_native_KJ_MPa_sqrt_m"],
        "onset_envelope_min_MPa_sqrt_m": min(onset_values, default=np.nan),
        "onset_envelope_max_MPa_sqrt_m": max(onset_values, default=np.nan),
        "qualified_G_onset_min_J_m2": min(qualified, default=np.nan),
        "qualified_G_onset_max_J_m2": max(qualified, default=np.nan),
        "mean_event_size_um": float(np.mean(sizes)) if sizes else np.nan,
        "median_event_size_um": float(np.median(sizes)) if sizes else np.nan,
        "terminal_extension_um": result["terminal_extension_m"] * 1.0e6,
        "terminal_opening_um": result["terminal_opening_m"] * 1.0e6,
        "max_tip_radius_um": result["max_tip_radius_m"] * 1.0e6,
        "max_backstress_GPa": result["max_backstress_Pa"] * 1.0e-9,
        "minimum_front_width_um": result["min_front_width_m"] * 1.0e6,
        "max_source_multiplicity": result["max_source_multiplicity"],
        "model_contract": result["model_contract"],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("baselines", "current", "long", "all"), default="all", nargs="?")
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    physics, rows, providers = inputs()
    payloads = []
    if args.mode in {"baselines", "all"}:
        for backend, (mechanics, drive) in providers.items():
            for _, row in rows[rows.candidate_id.isin(("v913_zeroD_sobol_0242980", "v913_zeroD_sobol_0202500"))].iterrows():
                material = IDS[row.candidate_id]
                payloads.append(run_case(row, material, 1000.0, backend, mechanics, drive, physics, 100.0))
        base = pd.DataFrame([summary(x) for x in payloads])
        pf2 = pd.read_csv(MAPS / "pf_2d_event_transactions_v2.csv")
        pf_first = pf2.groupby("material_class", as_index=False).first().set_index("material_class")
        for index, row in base.iterrows():
            if row.provider == "PF":
                target = float(pf_first.loc[row.material_class, "pre_event_native_KJ_MPa_sqrt_m"])
                base.loc[index, "authoritative_2D_initial_onset"] = target
                base.loc[index, "initial_onset_relative_error"] = (row.first_event_native_KJ_MPa_sqrt_m - target) / target
                base.loc[index, "authoritative_2D_physical_avalanche_count"] = 1 if row.material_class == "Peak" else 2
            else:
                target = 7454.99 if row.material_class == "Peak" else 8675.60
                base.loc[index, "authoritative_2D_initial_onset"] = target
                base.loc[index, "initial_onset_relative_error"] = (row.qualified_G_onset_min_J_m2 - target) / target
                base.loc[index, "authoritative_2D_physical_avalanche_count"] = 1 if row.material_class == "Peak" else 3
            base.loc[index, "topology_match"] = (
                int(row.physical_avalanche_count) == int(base.loc[index, "authoritative_2D_physical_avalanche_count"])
            )
        base.to_csv(OUT / "oneD_v2_native_baseline_results.csv", index=False)
        (OUT / "oneD_v2_native_baseline_cases.json").write_text(json.dumps(payloads, indent=2, sort_keys=True, default=_jsonable) + "\n")
    if args.mode in {"current", "all"}:
        payloads = []
        for backend, (mechanics, drive) in providers.items():
            for _, row in rows.iterrows():
                material = IDS[row.candidate_id]
                for temperature in TEMPERATURES:
                    payloads.append(run_case(row, material, temperature, backend, mechanics, drive, physics, 100.0))
        pd.DataFrame([summary(x) for x in payloads]).to_csv(
            OUT / "oneD_v2_current_four_class_results.csv", index=False
        )
        (OUT / "oneD_v2_current_four_class_cases.json").write_text(
            json.dumps(payloads, indent=2, sort_keys=True, default=_jsonable) + "\n"
        )
    if args.mode in {"long", "all"}:
        payloads = []
        for backend, (mechanics, drive) in providers.items():
            for _, row in rows[rows.candidate_id.isin(("v913_zeroD_sobol_0242980", "v913_zeroD_sobol_0202500"))].iterrows():
                payloads.append(run_case(row, IDS[row.candidate_id], 1000.0, backend, mechanics, drive, physics, 1000.0))
        pd.DataFrame([summary(x) for x in payloads]).to_csv(OUT / "oneD_v2_1000um_baseline_results.csv", index=False)
        (OUT / "oneD_v2_1000um_baseline_cases.json").write_text(json.dumps(payloads, indent=2, sort_keys=True, default=_jsonable) + "\n")
    manifest = {
        "schema": "oneD_v2_predictive_campaign_manifest_v1",
        "mode": args.mode, "canonical_parameters_changed": False,
        "new_2D_FEMCZM_runs": 0, "new_2D_PF_runs": 0,
        "physics_json_sha256": _sha(ROOT / "mpz_v9_13_v10222_transfer_common_physics.json"),
        "registry_sha256": _sha(REGISTRY),
        "source_drive_maps": {b: _sha(SOURCE_MAPS / f"oneD_v2_{b}_source_drive_map.csv") for b in ("pf", "fem")},
    }
    (OUT / "oneD_v2_predictive_campaign_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
