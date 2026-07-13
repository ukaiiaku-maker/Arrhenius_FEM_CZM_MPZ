#!/usr/bin/env python3
"""Run a detailed V1 Paris-style da/dN versus DeltaK analysis for a selected
subset of cleavage EXP-floor cases from a previous sweep_summary.csv.

The earlier sweep classified cases by first crack-advance/first-passage.  This
script reruns selected cases with --continue-after-fire so that the measured
quantity is crack growth per cycle at fixed K amplitude:

    da/dN = total_advance / total_cycles

If no advance occurs by cycles_max, the point is written as a censored upper
bound:

    da/dN < da_event / cycles_max

The model is still V1/K-controlled: Kmax is held fixed during repeated renewal
advances, matching the usual controlled-DeltaK interpretation for a crack-growth
curve.  v8 validation should be run later on only the most interesting curves.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Iterable, List, Tuple

import numpy as np
import pandas as pd


def _num(x, default=np.nan):
    try:
        if x is None or (isinstance(x, float) and np.isnan(x)):
            return default
        return float(x)
    except Exception:
        return default




def _canon_name(name: str) -> str:
    """Lowercase alphanumeric key for tolerant CSV column matching."""
    return "".join(ch for ch in str(name).lower() if ch.isalnum())


def _first_existing_column(df: pd.DataFrame, aliases: Iterable[str]):
    cmap = {_canon_name(c): c for c in df.columns}
    for a in aliases:
        key = _canon_name(a)
        if key in cmap:
            return cmap[key]
    return None


def _extract_case_from_dir(value) -> float:
    import re
    m = re.search(r"case[_-]?(\d+)", str(value))
    if m:
        return float(int(m.group(1)))
    return np.nan


def _float_from_path_token(tok: str) -> float:
    """Parse compact path tokens used by the sweep writer, e.g. 1p2 -> 1.2."""
    try:
        return float(str(tok).replace("p", "."))
    except Exception:
        return np.nan


def _extract_case_params_from_path(value) -> dict:
    """Extract case id and EXP-floor parameters from paths such as
    case_0024_G2_sc3p5_a1_n3_ff0p02/K10p0/...

    This is needed when a user accidentally passes a per-K summary CSV whose
    rows contain only history paths rather than the original sweep parameter
    table.
    """
    import re
    text = str(value)
    pat = re.compile(
        r"case[_-]?(?P<case>\d+).*?"
        r"G(?P<G00>[-+0-9p.eE]+)_"
        r"sc(?P<sigc>[-+0-9p.eE]+)_"
        r"a(?P<a>[-+0-9p.eE]+)_"
        r"n(?P<n>[-+0-9p.eE]+)_"
        r"ff(?P<ff>[-+0-9p.eE]+)"
    )
    m = pat.search(text)
    if not m:
        return {}
    return {
        "case": float(int(m.group("case"))),
        "G00_eV": _float_from_path_token(m.group("G00")),
        "sigc0_GPa": _float_from_path_token(m.group("sigc")),
        "a": _float_from_path_token(m.group("a")),
        "n": _float_from_path_token(m.group("n")),
        "floor_frac": _float_from_path_token(m.group("ff")),
    }


def _looks_like_per_k_result_table(df: pd.DataFrame) -> bool:
    cols = {_canon_name(c) for c in df.columns}
    return ("kmaxmpasqrtm" in cols or "kmaxmpasqrtm" in cols or "deltakmpasqrtm" in cols) and (
        "historycsv" in cols or "cyclesTotal".lower() in cols or "cycles_total" in df.columns
    )


def _find_alternate_sweep_summary(input_path: Path, case_ids: List[int]) -> Path | None:
    """Find an original EXP-floor sweep_summary.csv when the provided CSV is
    actually a per-K result table.  This keeps the workflow usable when a root
    sweep_summary.csv has been overwritten by a downstream summary.
    """
    search_roots = []
    for root in [input_path.parent, Path.cwd(), Path.cwd() / "runs"]:
        try:
            root = root.resolve()
        except Exception:
            pass
        if root.exists() and root not in search_roots:
            search_roots.append(root)
    candidates = []
    for root in search_roots:
        try:
            if root.is_file():
                continue
            # Keep search bounded enough for large run directories.
            patterns = ["sweep_summary.csv", "*/sweep_summary.csv", "*/*/sweep_summary.csv"]
            for pat in patterns:
                candidates.extend(root.glob(pat))
        except Exception:
            continue
    seen = set()
    for p in candidates:
        try:
            rp = p.resolve()
            if rp in seen or rp == input_path.resolve():
                continue
            seen.add(rp)
            df = pd.read_csv(rp, nrows=2000)
            norm = normalize_sweep_summary(df)  # recursion is safe if candidate has params or parsable paths
            if case_ids:
                have = set(norm["case"].astype(int))
                if not set(case_ids).issubset(have):
                    continue
            return rp
        except Exception:
            continue
    return None


def normalize_sweep_summary(summary: pd.DataFrame) -> pd.DataFrame:
    """Accept sweep-summary CSVs produced by older/newer scripts.

    Older versions used slightly different column names or, in some cases, only
    encoded the case id in case_dir.  The Paris subset runner needs a canonical
    set of columns; this adapter preserves all original columns and adds/renames
    the required ones.
    """
    out = summary.copy()

    alias_map = {
        "case": ["case", "case_id", "caseid", "case_idx", "case_index", "idx", "id"],
        "regime": ["regime", "classification", "class", "sweep_regime", "original_regime"],
        "G00_eV": ["G00_eV", "G0_eV", "G00", "G0", "cleave_G00_eV", "cleave_G0_eV", "cleave_g00_ev"],
        "sigc0_GPa": ["sigc0_GPa", "sigc_GPa", "sigma_c_GPa", "sigc0", "sigma_c", "cleave_sigc0_GPa", "cleave_sigc_GPa"],
        "a": ["a", "exp_a", "cleave_exp_a", "cleavage_exp_a"],
        "n": ["n", "exp_n", "cleave_exp_n", "cleavage_exp_n"],
        "floor_frac": ["floor_frac", "floor", "cleave_floor_frac", "cleavage_floor_frac", "floor_fraction"],
        "n_fire_conditions": ["n_fire_conditions", "num_fire_conditions", "n_fired", "n_fire", "fire_count"],
        "min_fire_cycles": ["min_fire_cycles", "min_cycles_to_fire", "min_cycles", "min_N_fire"],
        "case_dir": ["case_dir", "run_dir", "out_dir", "directory", "path", "_summary_path", "summary_path", "history_csv", "history_png"],
    }

    for canonical, aliases in alias_map.items():
        if canonical in out.columns:
            continue
        found = _first_existing_column(out, aliases)
        if found is not None:
            out[canonical] = out[found]

    # If only paths are present, recover case id and EXP-floor parameters from
    # the encoded case directory name.  This handles per-K result tables with
    # columns like history_csv or _summary_path.
    path_cols = [c for c in ["case_dir", "history_csv", "history_png", "_summary_path", "summary_path"] if c in out.columns]
    if path_cols:
        extracted_records = []
        for _, row in out.iterrows():
            rec = {}
            for pc in path_cols:
                rec = _extract_case_params_from_path(row.get(pc, ""))
                if rec:
                    break
            extracted_records.append(rec)
        if any(extracted_records):
            for key in ["case", "G00_eV", "sigc0_GPa", "a", "n", "floor_frac"]:
                if key not in out.columns or out[key].isna().all():
                    vals = [rec.get(key, np.nan) for rec in extracted_records]
                    if np.isfinite(pd.to_numeric(pd.Series(vals), errors="coerce")).any():
                        out[key] = vals

    if "case" not in out.columns:
        if "case_dir" in out.columns:
            extracted = out["case_dir"].map(_extract_case_from_dir)
            if extracted.notna().any():
                out["case"] = extracted
        if "case" not in out.columns:
            # Last-resort fallback: row index.  This keeps --auto-n usable and
            # also matches old sweeps whose rows were written in case order.
            out["case"] = np.arange(len(out), dtype=int)

    required = ["case", "G00_eV", "sigc0_GPa", "a", "n", "floor_frac"]
    missing = [c for c in required if c not in out.columns]
    if missing:
        raise ValueError(
            "sweep_summary.csv is missing required parameter columns after alias matching: "
            + ", ".join(missing)
            + "\nAvailable columns: " + ", ".join(map(str, summary.columns))
        )

    # Defaults for optional columns used by automatic selection/classification.
    if "regime" not in out.columns:
        out["regime"] = "unclassified_from_input"
    if "n_fire_conditions" not in out.columns:
        out["n_fire_conditions"] = np.nan
    if "min_fire_cycles" not in out.columns:
        out["min_fire_cycles"] = np.nan
    if "case_dir" not in out.columns:
        out["case_dir"] = ""

    # Numeric coercion, with explicit error for required columns.
    for c in ["case", "G00_eV", "sigc0_GPa", "a", "n", "floor_frac", "n_fire_conditions", "min_fire_cycles"]:
        out[c] = pd.to_numeric(out[c], errors="coerce")
    bad = [c for c in required if out[c].isna().any()]
    if bad:
        raise ValueError(
            "sweep_summary.csv has non-numeric or missing values in required columns: "
            + ", ".join(bad)
        )
    out["case"] = out["case"].astype(int)
    return out

def run(cmd: List[str], log: Path) -> None:
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("w") as f:
        f.write("$ " + " ".join(cmd) + "\n\n")
        f.flush()
        p = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)
    if p.returncode:
        raise RuntimeError(f"command failed with status {p.returncode}; see {log}")


def select_cases(summary: pd.DataFrame, case_ids: List[int], n_auto: int) -> pd.DataFrame:
    if case_ids:
        out = summary[summary["case"].isin(case_ids)].copy()
        missing = sorted(set(case_ids) - set(out["case"].astype(int)))
        if missing:
            raise ValueError(f"Requested case id(s) not found in summary: {missing}")
        # Per-K tables can contain several rows per case.  We only need one
        # parameter row for rerunning each selected case.
        return out.sort_values(["case"]).drop_duplicates("case", keep="first")

    picks = []
    # Candidate smooth/all-fire cases, preferring low n and low floor.
    all_fire = summary[summary["regime"].eq("all_fire_near_horizon")].copy()
    if len(all_fire):
        all_fire = all_fire.sort_values(["n", "floor_frac", "G00_eV", "sigc0_GPa", "a"])
        picks += list(all_fire.head(max(1, n_auto // 2))["case"].astype(int))

    # Threshold-like cases with several firing points but not all.  These are useful
    # controls for distinguishing a Paris-like curve from a cliff.
    th = summary[summary["regime"].eq("apparent_threshold_in_window")].copy()
    if len(th):
        th = th[(th["n_fire_conditions"] >= 3) & (th["n_fire_conditions"] <= 5)]
        th = th.sort_values(["n_fire_conditions", "min_fire_cycles"], ascending=[False, False])
        picks += list(th.head(max(1, n_auto // 3))["case"].astype(int))

    # One direct/overdriven all-fire case, as a negative control.
    direct = summary[summary["regime"].eq("direct_fracture_no_limit_in_window")].copy()
    if len(direct):
        direct = direct.sort_values(["n", "G00_eV", "sigc0_GPa", "a"])
        picks += list(direct.head(1)["case"].astype(int))

    # Preserve order and cap count.
    seen = set()
    pick_unique = []
    for c in picks:
        if c not in seen:
            seen.add(c)
            pick_unique.append(c)
        if len(pick_unique) >= n_auto:
            break
    return summary[summary["case"].isin(pick_unique)].copy().sort_values("case")


def parse_history(hist: pd.DataFrame, da_event_m: float, cycles_max: float, min_adv_measured: int) -> dict:
    last = hist.iloc[-1]
    fired = hist[hist.get("n_fire", 0) > 0]
    n_adv = int(last.get("n_adv", 0))
    a_adv = float(last.get("a_adv_m", 0.0))
    cycles_total = float(last.get("cycles_total", np.nan))
    first_fire = float(fired.iloc[0]["cycles_total"]) if len(fired) else np.nan
    last_fire = float(fired.iloc[-1]["cycles_total"]) if len(fired) else np.nan
    if n_adv > 0 and cycles_total > 0:
        da_dN = a_adv / cycles_total
        bound = np.nan
        if n_adv >= min_adv_measured:
            status = "measured_multi_event"
        else:
            status = "single_event_estimate"
    else:
        da_dN = np.nan
        # conservative upper bound: if one event would have been visible.
        denom = cycles_total if np.isfinite(cycles_total) and cycles_total > 0 else cycles_max
        bound = da_event_m / denom
        status = "censored_upper_bound"
    return {
        "cycles_total": cycles_total,
        "cycles_to_first_fire": first_fire,
        "cycles_to_last_fire": last_fire,
        "n_adv": n_adv,
        "a_adv_m": a_adv,
        "da_dN_m_per_cycle": da_dN,
        "da_dN_upper_bound_m_per_cycle": bound,
        "log10_da_dN": math.log10(da_dN) if da_dN and da_dN > 0 else np.nan,
        "log10_da_dN_bound": math.log10(bound) if bound and bound > 0 else np.nan,
        "point_status": status,
        "B_final": float(last.get("B", np.nan)),
        "N_em_final": float(last.get("N_em", np.nan)),
        "G_cleave_eff_eV": float(last.get("G_cleave_eff_eV", np.nan)),
        "G_cleave_raw_eV": float(last.get("G_cleave_raw_eV", np.nan)),
        "S_cleave_kB": float(last.get("S_cleave_kB", np.nan)),
        "dGcleave_dsigma_eV_per_GPa": float(last.get("dGcleave_dsigma_eV_per_GPa", np.nan)),
        "vstar_cleave_b3": float(last.get("vstar_cleave_b3", np.nan)),
        "mu_cleave_per_cycle": float(last.get("mu_cleave_pred", np.nan)),
        "mu_emit_per_cycle": float(last.get("mu_emit", np.nan)),
        "store_per_cycle": float(last.get("store_per_cycle", np.nan)),
        "storage_fraction": float(last.get("storage_fraction", np.nan)),
        "cycle_limiter_last": str(last.get("cycle_limiter", "")),
    }


def classify_paris_case(points: pd.DataFrame, cycles_max: float) -> dict:
    p = points.sort_values("DeltaK_MPa_sqrtm")
    measured = p[p["da_dN_m_per_cycle"].notna()].copy()
    censored = p[p["point_status"].eq("censored_upper_bound")].copy()
    directish = measured[measured["cycles_to_first_fire"].fillna(np.inf) < 1.0]
    cls = "unclassified"
    reason = ""
    max_adjacent_log_jump = np.nan
    slope_logdadn_per_DK = np.nan
    if len(measured) == 0:
        cls = "inactive_or_below_growth_resolution"
        reason = "no measured crack advance at any DeltaK"
    else:
        mm = measured.sort_values("DeltaK_MPa_sqrtm")
        y = mm["log10_da_dN"].to_numpy(float)
        x = mm["DeltaK_MPa_sqrtm"].to_numpy(float)
        if len(mm) >= 2:
            jumps = np.abs(np.diff(y))
            max_adjacent_log_jump = float(np.nanmax(jumps)) if len(jumps) else np.nan
            slope_logdadn_per_DK = float(np.polyfit(x, y, 1)[0]) if len(mm) >= 2 else np.nan
        low_censored = 0
        if len(censored):
            kth_min_meas = float(mm["DeltaK_MPa_sqrtm"].min())
            low_censored = int((censored["DeltaK_MPa_sqrtm"] < kth_min_meas).sum())
        if len(directish) >= max(2, len(measured)//2):
            cls = "overdriven_direct_fracture_window"
            reason = "many measured points first fire in less than one cycle"
        elif len(measured) >= 4 and low_censored == 0 and (not np.isfinite(max_adjacent_log_jump) or max_adjacent_log_jump <= 4.0):
            cls = "smooth_paris_like_no_threshold_in_window"
            reason = "multi-point finite growth curve without low-DeltaK censoring cliff"
        elif len(measured) >= 4 and low_censored <= 1 and (not np.isfinite(max_adjacent_log_jump) or max_adjacent_log_jump <= 5.0):
            cls = "mostly_smooth_with_low_DeltaK_bound"
            reason = "finite growth over most of window with limited low-DeltaK censoring"
        elif low_censored > 0:
            cls = "threshold_like_low_DeltaK_censored"
            reason = "lower-DeltaK points are censored while higher-DeltaK points grow"
        else:
            cls = "finite_growth_but_cliffy_or_sparse"
            reason = "finite growth exists but curve is too sparse or too steep for smooth Paris-like label"
    return {
        "paris_class": cls,
        "paris_reason": reason,
        "n_points": int(len(points)),
        "n_measured": int(len(measured)),
        "n_censored": int(len(censored)),
        "n_direct_lt_1_cycle": int(len(directish)),
        "max_adjacent_log_da_dN_jump": max_adjacent_log_jump,
        "slope_log_da_dN_per_DeltaK": slope_logdadn_per_DK,
    }


def make_plot(points: pd.DataFrame, case_row: pd.Series, out_png: Path) -> None:
    import matplotlib.pyplot as plt
    p = points.sort_values("DeltaK_MPa_sqrtm")
    fig, ax = plt.subplots(figsize=(6.2, 4.4))
    meas = p[p["da_dN_m_per_cycle"].notna()]
    cens = p[p["da_dN_upper_bound_m_per_cycle"].notna()]
    if len(meas):
        ax.plot(meas["DeltaK_MPa_sqrtm"], meas["da_dN_m_per_cycle"], marker="o", linestyle="-", label="measured")
    if len(cens):
        ax.scatter(cens["DeltaK_MPa_sqrtm"], cens["da_dN_upper_bound_m_per_cycle"], marker="v", label="upper bound")
    ax.set_yscale("log")
    ax.set_xlabel(r"$\Delta K$ (MPa $\sqrt{m}$)")
    ax.set_ylabel(r"$da/dN$ (m/cycle)")
    ax.set_title(f"case {int(case_row['case'])}: G00={case_row['G00_eV']:g}, sigc={case_row['sigc0_GPa']:g}, a={case_row['a']:g}, n={case_row['n']:g}, f={case_row['floor_frac']:g}")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="best")
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=200)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sweep-summary", required=True, help="Path to previous sweep_summary.csv")
    ap.add_argument("--out", default="runs/v1_paris_subset_from_sweep")
    ap.add_argument("--case-ids", nargs="*", type=int, default=[], help="Explicit case ids from sweep_summary.csv. If omitted, representative cases are selected automatically.")
    ap.add_argument("--auto-n", type=int, default=8, help="Number of cases for automatic selection.")
    ap.add_argument("--Kmax-MPa-sqrt-m", nargs="+", type=float, default=[5,6,7,8,9,10,11])
    ap.add_argument("--T", type=float, default=300.0)
    ap.add_argument("--R", type=float, default=0.1)
    ap.add_argument("--frequency-Hz", type=float, default=1000.0)
    ap.add_argument("--cycles-max", type=float, default=1e11)
    ap.add_argument("--max-blocks", type=int, default=400)
    ap.add_argument("--n-advances", type=int, default=10, help="Target repeated renewals for measured da/dN at fixed K.")
    ap.add_argument("--min-adv-measured", type=int, default=3, help="Fewer advances are labeled single_event_estimate.")
    ap.add_argument("--da", type=float, default=2.0e-5, help="V1 event advance in meters.")
    ap.add_argument("--block-cycles", type=float, default=1e5)
    ap.add_argument("--max-block-cycles", type=float, default=float("inf"))
    ap.add_argument("--cycle-block-mode", choices=["requested_cap", "hazard_limited"], default="hazard_limited")
    ap.add_argument("--target-dB", type=float, default=0.02)
    ap.add_argument("--target-dN-store", type=float, default=0.025)
    ap.add_argument("--target-dN-emit", type=float, default=0.25)
    ap.add_argument("--target-dN-mobile", type=float, default=0.25)
    ap.add_argument("--storage-model", choices=["fixed_fraction", "all_retained", "escape_limited"], default="fixed_fraction")
    ap.add_argument("--fixed-retained-fraction", type=float, default=0.1)
    ap.add_argument("--cleave-exp-T-mode", choices=["linear", "mu_scale"], default="mu_scale")
    ap.add_argument("--keep-existing", action="store_true")
    ap.add_argument("--no-plots", action="store_true")
    args = ap.parse_args()

    input_summary_path = Path(args.sweep_summary)
    raw_summary = pd.read_csv(input_summary_path)
    try:
        summary = normalize_sweep_summary(raw_summary)
    except ValueError as exc:
        alt = _find_alternate_sweep_summary(input_summary_path, args.case_ids)
        if alt is not None:
            print("Provided CSV does not contain case parameter columns; using alternate sweep summary:", alt)
            input_summary_path = alt
            raw_summary = pd.read_csv(input_summary_path)
            summary = normalize_sweep_summary(raw_summary)
        else:
            raise ValueError(
                str(exc)
                + "\n\nThis file looks like a per-K result table rather than the original EXP-floor case sweep summary. "
                + "Use the sweep_summary.csv written by run_v1_cleave_exp_floor_sweep.py, or pass a CSV whose paths "
                + "include case_####_G..._sc..._a..._n..._ff... so the parameters can be recovered."
            )
    print("Loaded sweep summary columns:", ", ".join(map(str, raw_summary.columns)))
    print("Using sweep summary file:", input_summary_path)
    print("Canonicalized required columns: case, G00_eV, sigc0_GPa, a, n, floor_frac")
    cases = select_cases(summary, args.case_ids, args.auto_n)
    if len(cases) == 0:
        raise SystemExit("No cases selected.")

    out = Path(args.out)
    if out.exists() and not args.keep_existing:
        shutil.rmtree(out)
    out.mkdir(parents=True, exist_ok=True)
    cases.to_csv(out / "selected_cases.csv", index=False)

    all_points = []
    case_rows = []
    for _, crow in cases.iterrows():
        cid = int(crow["case"])
        cdir = out / f"case_{cid:04d}_G{crow['G00_eV']:g}_sc{crow['sigc0_GPa']:g}_a{crow['a']:g}_n{crow['n']:g}_ff{crow['floor_frac']:g}".replace(".", "p")
        case_points = []
        for K in args.Kmax_MPa_sqrt_m:
            kout = cdir / f"K{K:g}".replace(".", "p")
            cmd = [sys.executable, "-m", "arrhenius_fracture.fatigue_sharp_front",
                   "--temperatures", str(args.T),
                   "--Kmax-MPa-sqrt-m", str(K),
                   "--R", str(args.R),
                   "--frequency-Hz", str(args.frequency_Hz),
                   "--cycles-max", str(args.cycles_max),
                   "--max-blocks", str(args.max_blocks),
                   "--n-advances", str(args.n_advances),
                   "--continue-after-fire",
                   "--da", str(args.da),
                   "--block-cycles", str(args.block_cycles),
                   "--max-block-cycles", str(args.max_block_cycles),
                   "--cycle-block-mode", args.cycle_block_mode,
                   "--target-dB", str(args.target_dB),
                   "--target-dN-store", str(args.target_dN_store),
                   "--target-dN-emit", str(args.target_dN_emit),
                   "--target-dN-mobile", str(args.target_dN_mobile),
                   "--storage-model", args.storage_model,
                   "--dN-cap", "inf", "--sigma-cap-GPa", "0", "--no-plots",
                   "--cleave-barrier-kind", "exp_floor",
                   "--cleave-exp-T-mode", args.cleave_exp_T_mode,
                   "--cleave-G00-eV", str(crow["G00_eV"]),
                   "--cleave-sigc0-GPa", str(crow["sigc0_GPa"]),
                   "--cleave-exp-a", str(crow["a"]),
                   "--cleave-exp-n", str(crow["n"]),
                   "--cleave-floor-frac", str(crow["floor_frac"]),
                   "--out", str(kout)]
            if args.storage_model == "fixed_fraction":
                cmd += ["--fixed-retained-fraction", str(args.fixed_retained_fraction)]
            run(cmd, cdir / f"K{K:g}.log".replace(".", "p"))
            hist_path = kout / f"T{int(round(args.T))}K" / "fatigue_v1_history.csv"
            hist = pd.read_csv(hist_path)
            pr = parse_history(hist, args.da, args.cycles_max, args.min_adv_measured)
            pr.update({
                "case": cid,
                "original_regime": str(crow.get("regime", "")),
                "G00_eV": float(crow["G00_eV"]),
                "sigc0_GPa": float(crow["sigc0_GPa"]),
                "cleave_exp_a": float(crow["a"]),
                "cleave_exp_n": float(crow["n"]),
                "floor_frac": float(crow["floor_frac"]),
                "Kmax_MPa_sqrtm": float(K),
                "DeltaK_MPa_sqrtm": float((1.0 - args.R) * K),
                "R": float(args.R),
                "T_K": float(args.T),
                "da_event_m": float(args.da),
                "run_dir": str(kout),
            })
            case_points.append(pr)
            all_points.append(pr)
        pdf = pd.DataFrame(case_points).sort_values("DeltaK_MPa_sqrtm")
        pdf.to_csv(cdir / "paris_points.csv", index=False)
        cls = classify_paris_case(pdf, args.cycles_max)
        cls.update({
            "case": cid,
            "original_regime": str(crow.get("regime", "")),
            "G00_eV": float(crow["G00_eV"]),
            "sigc0_GPa": float(crow["sigc0_GPa"]),
            "cleave_exp_a": float(crow["a"]),
            "cleave_exp_n": float(crow["n"]),
            "floor_frac": float(crow["floor_frac"]),
            "case_dir": str(cdir),
        })
        case_rows.append(cls)
        if not args.no_plots:
            try:
                make_plot(pdf, crow, cdir / "paris_da_dN_vs_DeltaK.png")
            except Exception as exc:
                print(f"WARNING: plotting failed for case {cid}: {exc}")
        pd.DataFrame(case_rows).to_csv(out / "paris_case_summary.csv", index=False)
        pd.DataFrame(all_points).to_csv(out / "paris_points.csv", index=False)
        print(f"case {cid}: {cls['paris_class']} ({cls['n_measured']} measured, {cls['n_censored']} censored)")

    with (out / "paris_settings.json").open("w") as f:
        json.dump(vars(args), f, indent=2, sort_keys=True)
    print(f"Wrote {out/'paris_points.csv'}")
    print(f"Wrote {out/'paris_case_summary.csv'}")


if __name__ == "__main__":
    main()
