#!/usr/bin/env python3
"""Audit historical DBTT candidates in the old and PF-equivalent reduced models."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

import optimize_mpz_v9_10_unified_global as old
from arrhenius_fracture.reduced_campaign_front_v9104 import (
    ReducedFrontSettings,
    best_adjacent_transition,
    simulate_reduced_response,
)

PARAMETER_COLUMNS = tuple(
    name
    for name in (
        "cleave_G00_eV", "cleave_gT_eV_per_K", "cleave_sigc0_GPa",
        "cleave_sT_GPa_per_K", "cleave_exp_a", "cleave_exp_n", "cleave_floor_frac",
        "emit_G00_eV", "emit_gT_eV_per_K", "emit_sigc0_GPa",
        "emit_sT_GPa_per_K", "emit_exp_a", "emit_exp_n", "emit_floor_frac",
        "peierls_H0_eV", "peierls_exp_a", "peierls_exp_n", "delta_H_PT_eV",
        "taylor_exp_a", "taylor_exp_n", "peierls_activation_entropy_kB",
        "taylor_activation_entropy_kB", "log10_taylor_corr_rho_c_m2",
        "log10_taylor_corr_scale", "log10_source_sites_per_system",
        "log10_encounter_efficiency", "log10_retained_recovery_rate_s",
        "log10_source_refresh_length_um", "c_blunt",
    )
)


def row_parameters(row: pd.Series) -> dict[str, float]:
    p = {name: float(row[name]) for name in PARAMETER_COLUMNS if name in row.index}
    p.update(
        {
            "taylor_H0_eV": float(row.get("taylor_H0_eV", p["peierls_H0_eV"] + p["delta_H_PT_eV"])),
            "taylor_corr_rho_c_m2": float(row.get("taylor_corr_rho_c_m2", 10.0 ** p["log10_taylor_corr_rho_c_m2"])),
            "taylor_corr_scale": float(row.get("taylor_corr_scale", 10.0 ** p["log10_taylor_corr_scale"])),
            "source_sites_per_system": float(row.get("source_sites_per_system", 10.0 ** p["log10_source_sites_per_system"])),
            "encounter_efficiency": float(row.get("encounter_efficiency", 10.0 ** p["log10_encounter_efficiency"])),
            "retained_recovery_rate_s": float(row.get("retained_recovery_rate_s", 10.0 ** p["log10_retained_recovery_rate_s"])),
            "source_refresh_length_um": float(row.get("source_refresh_length_um", 10.0 ** p["log10_source_refresh_length_um"])),
            "peierls_nu0_s": 1.0e12,
            "taylor_nu0_s": 1.0e11,
        }
    )
    return p


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", type=Path, required=True)
    ap.add_argument("--temperatures", default="300 400 500 600 700 800 900 1000 1100")
    ap.add_argument("--target-extension-um", type=float, default=5.0)
    ap.add_argument("--candidate-count", type=int, default=3)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    temperatures = np.asarray([float(x) for x in args.temperatures.split()], dtype=float)
    settings = ReducedFrontSettings(target_extension_um=float(args.target_extension_um))
    old_settings = old.ZeroDSettings(
        target_class="DBTT",
        temperatures=temperatures,
        targets=pd.DataFrame(),
        target_extension_um=float(args.target_extension_um),
    )
    manifest = pd.read_csv(args.manifest).head(args.candidate_count)
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    event_rows = []
    summaries = []
    for rank, row in manifest.iterrows():
        p = row_parameters(row)
        cid = str(row.get("candidate_id", f"candidate_{rank:02d}"))
        full_values = []
        off_values = []
        old_values = []
        for T in temperatures:
            new_full = simulate_reduced_response(p, float(T), settings, mode="full")
            new_off = simulate_reduced_response(p, float(T), settings, mode="plasticity_off")
            old_run = old.simulate_zero_d_rcurve(p, float(T), old_settings)
            full_values.append(float(new_full["K_init_proxy"]))
            off_values.append(float(new_off["K_init_proxy"]))
            old_values.append(float(old_run["K_init_proxy"]))
            rows.append(
                {
                    "candidate_id": cid,
                    "T_K": float(T),
                    "old_K_init": float(old_run["K_init_proxy"]),
                    "new_full_K_init": float(new_full["K_init_proxy"]),
                    "new_plasticity_off_K_init": float(new_off["K_init_proxy"]),
                    "new_plastic_increment": float(new_full["K_init_proxy"] - new_off["K_init_proxy"]),
                    "new_max_K_shield": float(new_full.get("max_K_shield_MPa_sqrt_m", np.nan)),
                    "new_max_sigma_back_Pa": float(new_full.get("max_sigma_back_Pa", np.nan)),
                }
            )
            event_rows.extend(
                {"candidate_id": cid, "T_K": float(T), **event}
                for event in new_full.get("events", [])
            )
        old_transition = best_adjacent_transition(temperatures, old_values)
        new_transition = best_adjacent_transition(
            temperatures,
            full_values,
            plasticity_off_toughness=off_values,
        )
        summaries.append(
            {
                "candidate_id": cid,
                "old_shelf_ratio": old_transition["shelf_ratio"],
                "old_jump_concentration": old_transition["jump_concentration"],
                "new_shelf_ratio": new_transition["shelf_ratio"],
                "new_jump_concentration": new_transition["jump_concentration"],
                "new_transition_low_K": new_transition["transition_low_K"],
                "new_transition_high_K": new_transition["transition_high_K"],
                "new_plasticity_off_ratio": new_transition["plasticity_off_ratio"],
                "cleave_gT_eV_per_K": p["cleave_gT_eV_per_K"],
                "cleave_sT_GPa_per_K": p["cleave_sT_GPa_per_K"],
            }
        )
    pd.DataFrame(rows).to_csv(out / "current_dbtt_old_vs_new_temperature_detail.csv", index=False)
    pd.DataFrame(event_rows).to_csv(out / "current_dbtt_new_event_detail.csv", index=False)
    pd.DataFrame(summaries).to_csv(out / "current_dbtt_old_vs_new_summary.csv", index=False)
    report = {
        "status": "V9_10_4_CURRENT_DBTT_AUDIT_COMPLETE",
        "manifest": str(args.manifest),
        "n_candidates": len(summaries),
        "target_extension_um": float(args.target_extension_um),
    }
    (out / "current_dbtt_audit.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
