#!/usr/bin/env python3
"""Coarse V1 Paris-style sweep over plasticity barrier scales at fixed cleavage barriers.

This driver is intended to answer: once a plausible crack-opening/cleavage
EXP-floor barrier is selected, does the Paris curve remain cleavage-clock
controlled, or is it sensitive to emission/Peierls/Taylor plasticity barriers?

The sweep uses a compact parameterization:

  E_emit       = --emit-primary-scale
  S_emit       = E_emit * --plastic-entropy-mult
  E_peierls    = E_emit * --peierls-ratio
  S_peierls    = E_emit * --peierls-ratio * --plastic-entropy-mult
  E_taylor     = E_emit * --taylor-ratio
  S_taylor     = E_emit * --taylor-ratio * --plastic-entropy-mult

These are passed to arrhenius_fracture.fatigue_sharp_front as the existing
--emit/--peierls/--taylor energy and entropy scale arguments.  Stress scales are
also exposed, but are held at 1 by default.
"""
from __future__ import annotations
import argparse
import itertools
import json
import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import numpy as np
import pandas as pd


REQUIRED_CLEAVE = ["case", "G00_eV", "sigc0_GPa", "a", "n", "floor_frac"]
ALIASES = {
    "case": ["case", "case_id", "case_idx", "case_index", "idx", "id"],
    "G00_eV": ["G00_eV", "G0_eV", "cleave_G00_eV", "cleave-G00-eV"],
    "sigc0_GPa": ["sigc0_GPa", "sigc_GPa", "sigma_c_GPa", "cleave_sigc0_GPa", "cleave-sigc0-GPa"],
    "a": ["a", "exp_a", "cleave_exp_a", "cleave-exp-a"],
    "n": ["n", "exp_n", "cleave_exp_n", "cleave-exp-n"],
    "floor_frac": ["floor_frac", "floor", "cleave_floor_frac", "cleave-floor-frac"],
}


def _parse_float_token(tok: str) -> float:
    return float(tok.replace("p", "."))


def _try_parse_case_path(path: str) -> Optional[Dict[str, float]]:
    # supports case_0029_G1_sc2p5_a0p7_n0p6_ff0p02
    s = str(path)
    m = re.search(r"case[_-]?(?P<case>\d+).*?_G(?P<G>[0-9p.+-]+)_sc(?P<sc>[0-9p.+-]+)_a(?P<a>[0-9p.+-]+)_n(?P<n>[0-9p.+-]+)_ff(?P<ff>[0-9p.+-]+)", s)
    if not m:
        return None
    return {
        "case": int(m.group("case")),
        "G00_eV": _parse_float_token(m.group("G")),
        "sigc0_GPa": _parse_float_token(m.group("sc")),
        "a": _parse_float_token(m.group("a")),
        "n": _parse_float_token(m.group("n")),
        "floor_frac": _parse_float_token(m.group("ff")),
    }


def normalize_cleavage_summary(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for canon, aliases in ALIASES.items():
        if canon in out.columns:
            continue
        for a in aliases:
            if a in out.columns:
                out[canon] = out[a]
                break
    missing = [c for c in REQUIRED_CLEAVE if c not in out.columns]
    if missing:
        path_cols = [c for c in ["case_dir", "history_csv", "history_png", "_summary_path", "summary_path"] if c in out.columns]
        parsed_rows = []
        if path_cols:
            for _, row in out.iterrows():
                parsed = None
                for pc in path_cols:
                    parsed = _try_parse_case_path(str(row.get(pc, "")))
                    if parsed is not None:
                        break
                parsed_rows.append(parsed or {})
            parsed_df = pd.DataFrame(parsed_rows)
            for col in REQUIRED_CLEAVE:
                if col not in out.columns and col in parsed_df.columns:
                    out[col] = parsed_df[col]
        missing = [c for c in REQUIRED_CLEAVE if c not in out.columns]
    if missing:
        raise ValueError(
            "Cleavage sweep summary is missing required columns after alias/path parsing: "
            + ", ".join(missing)
            + "\nAvailable columns: " + ", ".join(map(str, df.columns))
        )
    out = out.copy()
    out["case"] = pd.to_numeric(out["case"]).astype(int)
    for c in ["G00_eV", "sigc0_GPa", "a", "n", "floor_frac"]:
        out[c] = pd.to_numeric(out[c], errors="coerce")
    return out.dropna(subset=["G00_eV", "sigc0_GPa", "a", "n", "floor_frac"])


def run(cmd: List[str], log: Path) -> None:
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("w") as fp:
        fp.write("$ " + " ".join(cmd) + "\n\n")
        fp.flush()
        p = subprocess.run(cmd, stdout=fp, stderr=subprocess.STDOUT)
    if p.returncode:
        raise RuntimeError(f"command failed; see {log}\n" + " ".join(cmd))


def parse_history(hist: pd.DataFrame, da_event_m: float, cycles_max: float, min_adv_measured: int) -> Dict[str, object]:
    last = hist.iloc[-1]
    n_adv = int(last.get("n_adv", 0))
    a_adv = float(last.get("a_adv_m", 0.0))
    cycles_total = float(last.get("cycles_total", np.nan))
    fired = hist[hist.get("n_fire", 0) > 0]
    first_fire = float(fired.iloc[0]["cycles_total"]) if len(fired) else np.nan
    if n_adv > 0 and cycles_total > 0:
        da_dN = a_adv / cycles_total
        bound = np.nan
        status = "measured_multi_event" if n_adv >= min_adv_measured else "single_event_estimate"
    else:
        da_dN = np.nan
        denom = cycles_total if np.isfinite(cycles_total) and cycles_total > 0 else cycles_max
        bound = da_event_m / denom
        status = "censored_upper_bound"
    out = {
        "cycles_total": cycles_total,
        "cycles_to_first_fire": first_fire,
        "n_adv": n_adv,
        "a_adv_m": a_adv,
        "da_dN_m_per_cycle": da_dN,
        "da_dN_upper_bound_m_per_cycle": bound,
        "log10_da_dN": math.log10(da_dN) if da_dN and da_dN > 0 else np.nan,
        "log10_da_dN_bound": math.log10(bound) if bound and bound > 0 else np.nan,
        "point_status": status,
    }
    for col in [
        "B", "N_em", "sigma_back_Pa", "dG_emb_eV",
        "G_cleave_eff_eV", "S_cleave_kB", "dGcleave_dsigma_eV_per_GPa", "vstar_cleave_b3",
        "mu_cleave_pred", "mu_emit", "store_per_cycle", "storage_fraction",
        "G_emit_eV", "S_emit_kB", "dGemit_dsigma_eV_per_GPa", "vstar_emit_b3",
        "G_peierls_eV", "S_peierls_kB", "G_taylor_eV", "S_taylor_kB",
    ]:
        if col in last.index:
            out[col] = float(last.get(col, np.nan))
    return out


def classify(points: pd.DataFrame) -> Dict[str, object]:
    p = points.sort_values("DeltaK_MPa_sqrtm")
    measured = p[p["da_dN_m_per_cycle"].notna()].copy()
    cens = p[p["point_status"].eq("censored_upper_bound")].copy()
    direct = measured[measured["cycles_to_first_fire"].fillna(np.inf) < 1.0]
    max_jump = np.nan
    slope = np.nan
    cls = "unclassified"
    if len(measured) == 0:
        cls = "inactive_or_below_growth_resolution"
    else:
        x = measured["DeltaK_MPa_sqrtm"].to_numpy(float)
        y = measured["log10_da_dN"].to_numpy(float)
        if len(y) >= 2:
            max_jump = float(np.nanmax(np.abs(np.diff(y))))
            slope = float(np.polyfit(x, y, 1)[0])
        low_censored = 0
        if len(cens):
            low_censored = int((cens["DeltaK_MPa_sqrtm"] < float(measured["DeltaK_MPa_sqrtm"].min())).sum())
        if len(direct) >= max(2, len(measured)//2):
            cls = "overdriven_direct_fracture_window"
        elif len(measured) >= 6 and max_jump <= 1.5 and low_censored <= 2:
            cls = "measured_paris_like_with_limited_censoring"
        elif len(measured) >= 4 and max_jump <= 2.5:
            cls = "partial_paris_like_sensitive_window"
        elif low_censored > 0:
            cls = "threshold_like_low_DeltaK_censored"
        else:
            cls = "steep_or_irregular_growth_curve"
    return {
        "paris_class": cls,
        "n_measured": int(len(measured)),
        "n_censored": int(len(cens)),
        "n_direct_lt_1_cycle": int(len(direct)),
        "max_adjacent_log10_jump": max_jump,
        "slope_log10_da_dN_per_DeltaK": slope,
        "min_measured_da_dN": float(measured["da_dN_m_per_cycle"].min()) if len(measured) else np.nan,
        "max_measured_da_dN": float(measured["da_dN_m_per_cycle"].max()) if len(measured) else np.nan,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cleavage-summary", required=True, help="case-level cleavage sweep_summary.csv containing selected crack-opening barriers")
    ap.add_argument("--out", default="runs/v1_plasticity_barrier_sweep")
    ap.add_argument("--cleavage-case-ids", nargs="+", type=int, default=[29, 64])
    ap.add_argument("--Kmax-MPa-sqrt-m", nargs="+", type=float, default=[4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0])
    ap.add_argument("--T", type=float, default=300.0)
    ap.add_argument("--R", type=float, default=0.1)
    ap.add_argument("--frequency-Hz", type=float, default=1000.0)
    ap.add_argument("--cycles-max", type=float, default=1e11)
    ap.add_argument("--max-blocks", type=int, default=260)
    ap.add_argument("--n-advances", type=int, default=5)
    ap.add_argument("--min-adv-measured", type=int, default=3)
    ap.add_argument("--da", type=float, default=2.0e-5)
    ap.add_argument("--target-dB", type=float, default=0.02)
    ap.add_argument("--target-dN-store", type=float, default=0.025)
    ap.add_argument("--target-dN-emit", type=float, default=0.25)
    ap.add_argument("--target-dN-mobile", type=float, default=0.25)
    ap.add_argument("--storage-model", choices=["fixed_fraction","all_retained","escape_limited"], default="fixed_fraction")
    ap.add_argument("--fixed-retained-fraction", type=float, default=0.1)
    ap.add_argument("--emit-primary-scale", nargs="+", type=float, default=[0.75, 1.0, 1.25], help="Primary emission/nucleation energy scale(s).")
    ap.add_argument("--plastic-entropy-mult", nargs="+", type=float, default=[0.0, 0.5, 1.0, 2.0], help="Multiplier applied to mechanism entropy scales relative to their energy scales.")
    ap.add_argument("--peierls-ratio", nargs="+", type=float, default=[0.01, 0.02, 0.05], help="Peierls barrier ratio relative to emission primary scale.")
    ap.add_argument("--taylor-ratio", nargs="+", type=float, default=[0.05, 0.10, 0.20], help="Taylor barrier ratio relative to emission primary scale.")
    ap.add_argument("--emit-stress-scale", nargs="+", type=float, default=[1.0])
    ap.add_argument("--peierls-stress-scale", nargs="+", type=float, default=[1.0])
    ap.add_argument("--taylor-stress-scale", nargs="+", type=float, default=[1.0])
    ap.add_argument("--exp-system", default="W[100]", choices=["W[100]", "Ta[111]", "Al0.7CoCrFeNi-BCC", "Al0.7CoCrFeNi-FCC", "Cu"])
    ap.add_argument("--exp-a", type=float, default=None)
    ap.add_argument("--exp-n", type=float, default=None)
    ap.add_argument("--cleave-exp-T-mode", choices=["linear","mu_scale"], default="mu_scale")
    ap.add_argument("--keep-existing", action="store_true")
    args = ap.parse_args()

    out = Path(args.out)
    if out.exists() and not args.keep_existing:
        shutil.rmtree(out)
    out.mkdir(parents=True, exist_ok=True)

    cleave = normalize_cleavage_summary(pd.read_csv(args.cleavage_summary))
    selected = cleave[cleave["case"].isin(args.cleavage_case_ids)].copy().sort_values("case")
    if len(selected) != len(set(args.cleavage_case_ids)):
        got = set(map(int, selected["case"].tolist()))
        missing = sorted(set(args.cleavage_case_ids) - got)
        raise ValueError(f"Missing cleavage case ids in summary: {missing}")
    selected.to_csv(out / "selected_cleavage_cases.csv", index=False)

    all_points: List[Dict[str, object]] = []
    case_rows: List[Dict[str, object]] = []
    combos = list(itertools.product(
        selected.to_dict("records"),
        args.emit_primary_scale,
        args.plastic_entropy_mult,
        args.peierls_ratio,
        args.taylor_ratio,
        args.emit_stress_scale,
        args.peierls_stress_scale,
        args.taylor_stress_scale,
    ))
    for idx, (c, Eemit, Sentropy, rP, rT, sE, sP, sT) in enumerate(combos):
        cleave_case = int(c["case"])
        emit_E = float(Eemit)
        emit_S = float(Eemit) * float(Sentropy)
        peierls_E = float(Eemit) * float(rP)
        peierls_S = float(Eemit) * float(rP) * float(Sentropy)
        taylor_E = float(Eemit) * float(rT)
        taylor_S = float(Eemit) * float(rT) * float(Sentropy)
        plast_label = f"pc{cleave_case:04d}_E{emit_E:g}_Sm{Sentropy:g}_Pr{rP:g}_Tr{rT:g}_sE{sE:g}_sP{sP:g}_sT{sT:g}".replace(".", "p")
        case_dir = out / plast_label
        points: List[Dict[str, object]] = []
        for K in args.Kmax_MPa_sqrt_m:
            klabel = f"K{K:g}".replace(".", "p")
            kout = case_dir / klabel
            cmd = [
                sys.executable, "-m", "arrhenius_fracture.fatigue_sharp_front",
                "--temperatures", str(args.T),
                "--Kmax-MPa-sqrt-m", str(K),
                "--R", str(args.R), "--frequency-Hz", str(args.frequency_Hz),
                "--cycles-max", str(args.cycles_max), "--max-blocks", str(args.max_blocks),
                "--n-advances", str(args.n_advances), "--da", str(args.da),
                "--block-cycles", "1e5", "--max-block-cycles", "inf", "--cycle-block-mode", "hazard_limited",
                "--target-dB", str(args.target_dB), "--target-dN-store", str(args.target_dN_store),
                "--target-dN-emit", str(args.target_dN_emit), "--target-dN-mobile", str(args.target_dN_mobile),
                "--storage-model", args.storage_model, "--dN-cap", "inf", "--sigma-cap-GPa", "0", "--no-plots",
                "--cleave-barrier-kind", "exp_floor", "--cleave-exp-T-mode", args.cleave_exp_T_mode,
                "--cleave-G00-eV", str(c["G00_eV"]), "--cleave-sigc0-GPa", str(c["sigc0_GPa"]),
                "--cleave-exp-a", str(c["a"]), "--cleave-exp-n", str(c["n"]), "--cleave-floor-frac", str(c["floor_frac"]),
                "--exp-system", args.exp_system,
                "--emit-energy-scale", str(emit_E), "--emit-entropy-scale", str(emit_S), "--emit-stress-scale", str(sE),
                "--peierls-energy-scale", str(peierls_E), "--peierls-entropy-scale", str(peierls_S), "--peierls-stress-scale", str(sP),
                "--taylor-energy-scale", str(taylor_E), "--taylor-entropy-scale", str(taylor_S), "--taylor-stress-scale", str(sT),
                "--out", str(kout),
            ]
            if args.storage_model == "fixed_fraction":
                cmd += ["--fixed-retained-fraction", str(args.fixed_retained_fraction)]
            if args.exp_a is not None:
                cmd += ["--exp-a", str(args.exp_a)]
            if args.exp_n is not None:
                cmd += ["--exp-n", str(args.exp_n)]
            run(cmd, case_dir / f"{klabel}.log")
            hist = pd.read_csv(kout / f"T{int(args.T)}K" / "fatigue_v1_history.csv")
            rec = parse_history(hist, args.da, args.cycles_max, args.min_adv_measured)
            rec.update({
                "plastic_case_index": idx,
                "cleavage_case": cleave_case,
                "Kmax_MPa_sqrtm": K,
                "DeltaK_MPa_sqrtm": (1.0 - args.R) * K,
                "G00_eV": c["G00_eV"], "sigc0_GPa": c["sigc0_GPa"], "cleave_a": c["a"], "cleave_n": c["n"], "floor_frac": c["floor_frac"],
                "emit_energy_scale": emit_E, "emit_entropy_scale": emit_S, "plastic_entropy_mult": Sentropy,
                "peierls_energy_scale": peierls_E, "peierls_entropy_scale": peierls_S, "peierls_ratio": rP,
                "taylor_energy_scale": taylor_E, "taylor_entropy_scale": taylor_S, "taylor_ratio": rT,
                "emit_stress_scale": sE, "peierls_stress_scale": sP, "taylor_stress_scale": sT,
                "case_dir": str(case_dir), "history_csv": str(kout / f"T{int(args.T)}K" / "fatigue_v1_history.csv"),
            })
            points.append(rec)
            all_points.append(rec)
        pdf = pd.DataFrame(points)
        pdf.to_csv(case_dir / "paris_points.csv", index=False)
        summ = classify(pdf)
        summ.update({
            "plastic_case_index": idx,
            "cleavage_case": cleave_case,
            "G00_eV": c["G00_eV"], "sigc0_GPa": c["sigc0_GPa"], "cleave_a": c["a"], "cleave_n": c["n"], "floor_frac": c["floor_frac"],
            "emit_energy_scale": emit_E, "emit_entropy_scale": emit_S, "plastic_entropy_mult": Sentropy,
            "peierls_energy_scale": peierls_E, "peierls_entropy_scale": peierls_S, "peierls_ratio": rP,
            "taylor_energy_scale": taylor_E, "taylor_entropy_scale": taylor_S, "taylor_ratio": rT,
            "emit_stress_scale": sE, "peierls_stress_scale": sP, "taylor_stress_scale": sT,
            "case_dir": str(case_dir),
        })
        case_rows.append(summ)
        pd.DataFrame(all_points).to_csv(out / "paris_points.csv", index=False)
        pd.DataFrame(case_rows).to_csv(out / "plasticity_sweep_summary.csv", index=False)
        print(f"[{idx+1}/{len(combos)}] cleave={cleave_case} E={emit_E:g} Sm={Sentropy:g} Pr={rP:g} Tr={rT:g}: {summ['paris_class']} ({summ['n_measured']} measured, {summ['n_censored']} censored)")
    with (out / "sweep_settings.json").open("w") as fp:
        json.dump(vars(args), fp, indent=2)
    print("Wrote", out / "plasticity_sweep_summary.csv")
    print("Wrote", out / "paris_points.csv")


if __name__ == "__main__":
    main()
