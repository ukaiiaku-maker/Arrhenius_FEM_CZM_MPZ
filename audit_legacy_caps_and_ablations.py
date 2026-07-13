#!/usr/bin/env python3
"""Audit the frozen scalar first-passage model for active caps and compensations.

This script intentionally runs ``--front-state-model legacy_scalar``.  It is a
baseline/ablation tool, not the v9 constitutive model.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from arrhenius_fracture import sharp_front as sf
from arrhenius_fracture.config import make_emergent_config


def f(x):
    return str(float(x)) if np.isfinite(float(x)) else "inf"


def engine_args(row, dK, Kdot, ablation):
    argv = [
        "--mode", "1d", "--front-state-model", "legacy_scalar",
        "--temperatures", "300", "--Kdot", f(Kdot), "--Kmax", "60",
        "--dt", f(dK / Kdot), "--sigma-cap-GPa", "30",
        "--dN-cap", "50", "--v-emb-b3", "500", "--emb-sat-frac", "1",
        "--cleave-barrier-kind", "exp_floor", "--cleave-exp-T-mode", "linear",
        "--cleave-G00-eV", f(row.cleave_G00_eV),
        "--cleave-gT-eV-per-K", f(row.cleave_gT_eV_per_K),
        "--cleave-sigc0-GPa", f(row.cleave_sigc0_GPa),
        "--cleave-sT-GPa-per-K", f(row.cleave_sT_GPa_per_K),
        "--cleave-exp-a", f(row.cleave_exp_a), "--cleave-exp-n", f(row.cleave_exp_n),
        "--cleave-floor-frac", f(row.cleave_floor_frac),
        "--cleave-S-hs-kB", f(row.cleave_S_hs_kB),
        "--cleave-sigma-S-GPa", "6", "--cleave-S-hs-power", "2",
        "--emit-barrier-kind", "exp_floor",
        "--emit-G00-eV", f(row.emit_G00_eV_effective),
        "--emit-gT-eV-per-K", f(row.emit_gT_eV_per_K_effective),
        "--emit-sigc0-GPa", f(row.emit_sigc0_GPa),
        "--emit-sT-GPa-per-K", f(row.emit_sT_GPa_per_K),
        "--emit-exp-a", f(row.emit_exp_a), "--emit-exp-n", f(row.emit_exp_n),
        "--emit-floor-frac", f(row.emit_floor_frac),
        "--cleave-shield-chi", f(row.chi_shield), "--n-sat", f(row.N_sat),
    ]
    args = sf._build_parser().parse_args(argv)
    if ablation in ("no_dN_cap", "all_removed"):
        args.dN_cap = float("inf")
    if ablation in ("no_N_sat", "all_removed"):
        args.N_sat = float("inf")
    if ablation in ("no_chi_shield", "all_removed"):
        args.chi_shield = 0.0
    if ablation in ("no_stored_energy", "all_removed"):
        args.v_emb_b3 = 0.0
        args.emb_sat_frac = 0.0
    if ablation in ("no_sigma_cap", "all_removed"):
        args.sigma_cap_GPa = 0.0
    return args


def run_case(row, klass, T, dK, Kdot, Kmax, n_advances, ablation):
    args = engine_args(row, dK, Kdot, ablation)
    args.temperatures = [float(T)]
    args.Kmax = float(Kmax)
    mat = make_emergent_config().material
    eng = sf.build_engine(args, mat)
    dt = dK / Kdot
    records = []
    fire_K = []
    for i in range(int(math.ceil(Kmax / dK))):
        K = (i + 1) * dK * 1e6
        info = eng.step(K, T, dt)
        records.append({
            "step": i + 1, "K_MPa_sqrt_m": K / 1e6, "T_K": T,
            "B": info["B"], "N_em": info["N_em"],
            "sigma_tip_GPa": info["sigma_tip"] / 1e9,
            "sigma_back_GPa": info["sigma_back"] / 1e9,
            "r_eff_um": info["r_eff"] * 1e6,
            "lambda_e": info["lambda_e"], "lambda_c": info["lambda_c"],
            "dN_emit_raw": info.get("dN_emit_raw", 0.0),
            "dN_cap_active": int(info.get("dN_cap_active", False)),
            "N_sat_factor": info.get("N_sat_factor", 1.0),
            "N_sat_active": int(info.get("N_sat_active", False)),
            "sigma_cap_active": int(info.get("sigma_cap_active", False)),
            "dG_emb_eV": info.get("dG_emb_eV", 0.0),
            "G_cleave_raw_eV": info.get("G_cleave_raw_eV", np.nan),
            "G_cleave_eff_eV": info.get("G_cleave_eff_eV", np.nan),
            "n_fire": info.get("n_fire", 0), "a_adv_um": eng.a_adv * 1e6,
        })
        if info.get("n_fire", 0):
            fire_K.extend([K / 1e6] * int(info["n_fire"]))
        if eng.n_adv >= n_advances:
            break
    tr = pd.DataFrame(records)
    return tr, {
        "class": klass, "T_K": T, "dK_MPa_sqrt_m": dK,
        "dt_s": dt, "ablation": ablation,
        "first_K_MPa_sqrt_m": fire_K[0] if fire_K else np.nan,
        "last_K_MPa_sqrt_m": fire_K[-1] if fire_K else np.nan,
        "n_adv": int(eng.n_adv),
        "cap_step_fraction": float(tr.dN_cap_active.mean()) if len(tr) else 0.0,
        "N_sat_step_fraction": float(tr.N_sat_active.mean()) if len(tr) else 0.0,
        "sigma_cap_step_fraction": float(tr.sigma_cap_active.mean()) if len(tr) else 0.0,
        "max_N_em": float(tr.N_em.max()) if len(tr) else 0.0,
        "max_sigma_back_GPa": float(tr.sigma_back_GPa.max()) if len(tr) else 0.0,
        "max_dG_emb_eV": float(tr.dG_emb_eV.max()) if len(tr) else 0.0,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--parameters", default="legacy_four_class_v8_reference.csv")
    ap.add_argument("--out", default="runs/legacy_cap_ablation_audit_v9")
    ap.add_argument("--classes", default="ceramic peak weakT DBTT",
                    help="Subset of legacy target classes to audit.")
    ap.add_argument("--temperatures", default="300 500 700 900 1100")
    ap.add_argument("--dK", "--dK-values", dest="dK", default="0.25 0.05 0.02 0.005")
    ap.add_argument("--Kdot", type=float, default=0.005)
    ap.add_argument("--Kmax", type=float, default=60.0)
    ap.add_argument("--n-advances", type=int, default=20)
    ap.add_argument("--ablations", default="baseline no_dN_cap no_N_sat no_chi_shield no_stored_energy no_sigma_cap all_removed")
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    params = pd.read_csv(a.parameters)
    requested = set(a.classes.replace(",", " ").split())
    params = params[params.target_class.astype(str).isin(requested)]
    missing = requested.difference(set(params.target_class.astype(str)))
    if missing:
        raise SystemExit(f"unknown classes: {sorted(missing)}")
    temps = [float(x) for x in a.temperatures.replace(",", " ").split()]
    dks = [float(x) for x in a.dK.replace(",", " ").split()]
    ablations = a.ablations.replace(",", " ").split()
    summaries = []
    for _, row in params.iterrows():
        klass = str(row.target_class)
        for T in temps:
            for dk in dks:
                for abl in ablations:
                    tr, sm = run_case(row, klass, T, dk, a.Kdot, a.Kmax, a.n_advances, abl)
                    case = out / klass / f"T{int(T)}" / f"dK{dk:g}" / abl
                    case.mkdir(parents=True, exist_ok=True)
                    tr.to_csv(case / "trace.csv", index=False)
                    (case / "summary.json").write_text(json.dumps(sm, indent=2))
                    summaries.append(sm)
                    print(klass, T, dk, abl, "K1=", sm["first_K_MPa_sqrt_m"],
                          "cap=", sm["cap_step_fraction"], "Nsat=", sm["N_sat_step_fraction"])
    pd.DataFrame(summaries).to_csv(out / "audit_summary.csv", index=False)
    (out / "run_config.json").write_text(json.dumps(vars(a), indent=2))


if __name__ == "__main__":
    main()
