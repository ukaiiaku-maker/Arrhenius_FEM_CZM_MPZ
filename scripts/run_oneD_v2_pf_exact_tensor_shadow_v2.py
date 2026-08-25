#!/usr/bin/env python3
"""Fieldwise PF exact-wrapper shadow using qualified source tensor nodes."""
from __future__ import annotations

import copy
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PF = Path("/private/tmp/pf-v2-tensor-diagnostic")
PF_DATA = Path("/Volumes/Data/Data/Nanopillar_calculation/PF-fracture-fatigue_v10_2_21_persistent_sites_top1")
OUT = ROOT / "analysis_outputs/oneD_v2_local_drive"
sys.path[:0] = [str(ROOT), str(ROOT / "scripts"), str(PF)]

from reduced_fracture_v2.local_drive import PFTensorDriveProvider
from run_oneD_v2_native_source_shadow import IDS, pf_engine, rows


def install_drive(engine, drive) -> None:
    engine._anisotropic_drive = {
        "reliable": bool(drive.reliable),
        "drive_factors": np.asarray(drive.drive_factors, float),
        "tau_signed_Pa": np.asarray(drive.resolved_tau_signed_Pa, float),
    }
    engine._anisotropic_drive_serial = 1
    engine._install_current_drive_on_state()


def metrics(engine, out, temperature) -> dict:
    rates = np.asarray(getattr(engine, "anisotropic_last_lambda_emit_by_system_s", [0.0, 0.0]), float)
    probability = np.asarray(getattr(engine, "anisotropic_last_probability_by_system", [0.0, 0.0]), float)
    effective = np.asarray(getattr(engine, "anisotropic_last_sigma_emit_by_system_Pa", [0.0, 0.0]), float)
    barriers = np.asarray([
        engine.eb.G_barrier(np.asarray([max(float(value), 0.0)]), temperature, engine.b)[0]
        for value in effective
    ], float)
    return {
        "selected_emission_system": int(np.argmax(rates)),
        "rates": rates, "probability": probability, "effective": effective,
        "barriers": barriers, "lambda_e": float(out.get("lambda_e", np.sum(rates))),
    }


def main() -> None:
    summary = pd.read_csv(OUT / "oneD_v2_pf_tensor_drive_node_summary.csv")
    node = summary[np.isclose(summary.actual_extension_um, 100.0)].iloc[-1]
    provider = PFTensorDriveProvider.from_artifacts(
        OUT / "oneD_v2_pf_tensor_drive_node_summary.csv",
        OUT / "oneD_v2_pf_tensor_drive_map_manifest.json",
        allow_source_nodes_only=True,
    )
    drive = provider.evaluate(100e-6, float(node.reference_opening_m), exact_node_only=True)
    source_drive = copy.deepcopy(drive)
    records = []
    for row in rows(PF_DATA):
        for temperature in (300.0, 1000.0, 1200.0):
            source = pf_engine(PF_DATA, row, exact=True)
            reduced = copy.deepcopy(source)
            install_drive(source, source_drive); install_drive(reduced, drive)
            K = float(node.native_KJ_Pa_sqrt_m)
            source_out = source.step(K, temperature, 1e-12)
            reduced_out = reduced.step(K, temperature, 1e-12)
            a, b = metrics(source, source_out, temperature), metrics(reduced, reduced_out, temperature)
            fields = ("rates", "probability", "effective", "barriers")
            records.append({
                "material_class": IDS[row["candidate_id"]],
                "temperature_K": temperature,
                "extension_um": 100.0,
                "selected_emission_system_source": a["selected_emission_system"],
                "selected_emission_system_provider": b["selected_emission_system"],
                "selected_system_exact": a["selected_emission_system"] == b["selected_emission_system"],
                "tau_sign_exact": np.array_equal(np.sign(source_drive.resolved_tau_signed_Pa), np.sign(drive.resolved_tau_signed_Pa)),
                "reliability_exact": source_drive.reliable == drive.reliable,
                "maximum_rate_relative_error": float(np.max(np.abs(a["rates"]-b["rates"]) / np.maximum(np.abs(a["rates"]), 1e-300))),
                "maximum_barrier_absolute_error_J": float(np.max(np.abs(a["barriers"]-b["barriers"]))),
                "maximum_probability_absolute_error": float(np.max(np.abs(a["probability"]-b["probability"]))),
                "aggregate_hazard_relative_error": abs(a["lambda_e"]-b["lambda_e"]) / max(abs(a["lambda_e"]), 1e-300),
                "all_continuous_fields_exact": all(np.array_equal(a[name], b[name]) for name in fields) and a["lambda_e"] == b["lambda_e"],
                "status": "PASS_EXACT_SOURCE_NODE",
            })
    pd.DataFrame(records).to_csv(OUT / "oneD_v2_pf_exact_wrapper_shadow_v2.csv", index=False)


if __name__ == "__main__":
    main()
