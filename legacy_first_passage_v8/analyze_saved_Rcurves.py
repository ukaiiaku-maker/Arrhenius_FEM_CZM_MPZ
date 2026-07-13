#!/usr/bin/env python3
"""R-curve-style analysis for saved FEM/CZM seed simulations.

This script does not rerun simulations. It reads existing seed outputs and performs:
  - event-sampled KJ(Δa) extraction from R_curve_event_sampled.csv or steps_*K.csv
  - crack-extension binning using median K in each bin
  - early, late, AUC, slope, minimum, and amplitude metrics
  - saturating K_R(Δa) fit
  - power-law J_R(Δa) fit after converting K to J = K^2/E'
  - seed-level plots and class-level mean/binned plots

The analysis is "R-curve-like" rather than ASTM-valid; it is intended for comparing
simulation resistance histories across barrier classes and seeds.
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    from scipy.optimize import curve_fit
except Exception:
    curve_fit = None

CLASSES_DEFAULT = ["ceramic", "peak", "weakT", "DBTT"]


def first_existing(df: pd.DataFrame, names: Iterable[str]) -> str | None:
    lower = {c.lower(): c for c in df.columns}
    for n in names:
        if n in df.columns:
            return n
        if n.lower() in lower:
            return lower[n.lower()]
    return None


def parse_seed(path: Path) -> int | float:
    m = re.search(r"seed(\d+)", str(path))
    return int(m.group(1)) if m else np.nan


def parse_replicate(path: Path) -> int | float:
    m = re.search(r"replicate_(\d+)", str(path))
    return int(m.group(1)) if m else np.nan


def read_metadata(case_dir: Path, initial_a_mm: float, target_ext_um: float) -> dict:
    meta = {
        "case_dir": str(case_dir),
        "replicate": parse_replicate(case_dir),
        "seed": parse_seed(case_dir),
        "final_extension_um": np.nan,
        "target_hit": False,
        "complete": False,
        "status": "",
    }

    sj = case_dir / "summary.json"
    if sj.exists():
        try:
            s = json.loads(sj.read_text())
            if isinstance(s, list) and s:
                s = s[0]
            if s.get("status") is not None:
                meta["status"] = str(s.get("status"))
            if s.get("a_final_mm") is not None:
                meta["final_extension_um"] = (float(s["a_final_mm"]) - initial_a_mm) * 1000.0
            for key in ["final_crack_extension_um", "crack_extension_um", "extension_um"]:
                if s.get(key) is not None:
                    val = float(s[key])
                    meta["final_extension_um"] = max(meta["final_extension_um"], val) if math.isfinite(meta["final_extension_um"]) else val
        except Exception as e:
            meta["status"] += f" summary_read_error={e}"

    log = case_dir / "run.log"
    if log.exists():
        try:
            txt = log.read_text(errors="ignore")
            meta["target_hit"] = "reached target crack extension" in txt
            found = re.findall(r"reached target crack extension\s+([0-9.]+)\s+um", txt)
            if found:
                val = float(found[-1])
                meta["final_extension_um"] = max(meta["final_extension_um"], val) if math.isfinite(meta["final_extension_um"]) else val
        except Exception as e:
            meta["status"] += f" log_read_error={e}"

    meta["complete"] = bool(math.isfinite(meta["final_extension_um"]) and meta["final_extension_um"] >= 0.98 * target_ext_um)
    return meta


def read_curve(case_dir: Path) -> pd.DataFrame:
    """Return crack_extension_um and K_MPa_sqrt_m."""
    rc = case_dir / "R_curve_event_sampled.csv"
    if rc.exists():
        df = pd.read_csv(rc)
        if df.empty:
            return pd.DataFrame(columns=["crack_extension_um", "K_MPa_sqrt_m"])
        kcol = first_existing(df, ["KJ_MPa_sqrt_m", "K_MPa_sqrt_m", "Kc_MPa_sqrt_m", "K_MPa_sqrtm"])
        xcol = first_existing(df, ["crack_extension_um", "extension_um", "delta_a_um", "da_cumulative_um"])
        if kcol is None:
            return pd.DataFrame(columns=["crack_extension_um", "K_MPa_sqrt_m"])
        k = pd.to_numeric(df[kcol], errors="coerce")
        if xcol:
            x = pd.to_numeric(df[xcol], errors="coerce")
        else:
            x = np.arange(len(df), dtype=float)
        out = pd.DataFrame({"crack_extension_um": x, "K_MPa_sqrt_m": k})
        out = out.replace([np.inf, -np.inf], np.nan).dropna()
        out = out[out["K_MPa_sqrt_m"] > 0].copy()
        out = out.sort_values("crack_extension_um").drop_duplicates("crack_extension_um")
        return out.reset_index(drop=True)

    steps_files = sorted(case_dir.glob("steps_*K.csv"))
    if not steps_files:
        return pd.DataFrame(columns=["crack_extension_um", "K_MPa_sqrt_m"])
    df = pd.read_csv(steps_files[0])
    if df.empty:
        return pd.DataFrame(columns=["crack_extension_um", "K_MPa_sqrt_m"])

    keep = np.zeros(len(df), dtype=bool)
    da_col = first_existing(df, ["da_block_um", "da_block_m", "da_block"])
    if da_col:
        da = pd.to_numeric(df[da_col], errors="coerce").fillna(0).to_numpy(float)
        keep |= da > 0
    nfire_col = first_existing(df, ["n_fire", "nfire"])
    if nfire_col:
        nf = pd.to_numeric(df[nfire_col], errors="coerce").fillna(0).to_numpy(float)
        keep |= nf > 0

    ev = df.loc[keep].copy()
    if ev.empty:
        return pd.DataFrame(columns=["crack_extension_um", "K_MPa_sqrt_m"])

    kcol = first_existing(ev, ["KJ_MPa_sqrt_m", "K_MPa_sqrt_m", "Kc_MPa_sqrt_m", "KJ_Pa_sqrtm"])
    if kcol is None:
        return pd.DataFrame(columns=["crack_extension_um", "K_MPa_sqrt_m"])
    if kcol == "KJ_Pa_sqrtm":
        k = pd.to_numeric(ev[kcol], errors="coerce") / 1e6
    else:
        k = pd.to_numeric(ev[kcol], errors="coerce")

    xcol_um = first_existing(ev, ["crack_extension_um", "extension_um", "delta_a_um"])
    xcol_m = first_existing(ev, ["crack_extension_m"])
    if xcol_um:
        x = pd.to_numeric(ev[xcol_um], errors="coerce")
    elif xcol_m:
        x = pd.to_numeric(ev[xcol_m], errors="coerce") * 1e6
    else:
        x = np.arange(len(ev), dtype=float)

    out = pd.DataFrame({"crack_extension_um": x, "K_MPa_sqrt_m": k})
    out = out.replace([np.inf, -np.inf], np.nan).dropna()
    out = out[out["K_MPa_sqrt_m"] > 0].copy()
    out = out.sort_values("crack_extension_um").drop_duplicates("crack_extension_um")
    return out.reset_index(drop=True)


def binned_curve(curve: pd.DataFrame, bin_um: float, max_ext_um: float | None = None) -> pd.DataFrame:
    if curve.empty:
        return pd.DataFrame(columns=["bin_center_um", "K_median_MPa_sqrt_m", "K_mean_MPa_sqrt_m", "K_p10", "K_p90", "n"])
    x = curve["crack_extension_um"].to_numpy(float)
    y = curve["K_MPa_sqrt_m"].to_numpy(float)
    lo = 0.0
    hi = max_ext_um if max_ext_um is not None else max(float(np.nanmax(x)), bin_um)
    edges = np.arange(lo, hi + bin_um, bin_um)
    rows = []
    for a, b in zip(edges[:-1], edges[1:]):
        mask = (x >= a) & (x < b)
        vals = y[mask]
        vals = vals[np.isfinite(vals)]
        if len(vals) == 0:
            continue
        rows.append({
            "bin_start_um": a,
            "bin_end_um": b,
            "bin_center_um": 0.5 * (a + b),
            "K_median_MPa_sqrt_m": float(np.median(vals)),
            "K_mean_MPa_sqrt_m": float(np.mean(vals)),
            "K_p10_MPa_sqrt_m": float(np.quantile(vals, 0.10)),
            "K_p90_MPa_sqrt_m": float(np.quantile(vals, 0.90)),
            "n": int(len(vals)),
        })
    return pd.DataFrame(rows)


def window_median(curve: pd.DataFrame, a0: float, a1: float) -> float:
    if curve.empty:
        return np.nan
    m = (curve["crack_extension_um"] >= a0) & (curve["crack_extension_um"] <= a1)
    vals = curve.loc[m, "K_MPa_sqrt_m"].dropna()
    vals = vals[vals > 0]
    return float(vals.median()) if len(vals) else np.nan


def window_slope(curve: pd.DataFrame, a0: float, a1: float) -> float:
    if curve.empty:
        return np.nan
    m = (curve["crack_extension_um"] >= a0) & (curve["crack_extension_um"] <= a1)
    d = curve.loc[m, ["crack_extension_um", "K_MPa_sqrt_m"]].dropna()
    if len(d) < 3:
        return np.nan
    x = d["crack_extension_um"].to_numpy(float)
    y = d["K_MPa_sqrt_m"].to_numpy(float)
    if np.nanmax(x) <= np.nanmin(x):
        return np.nan
    return float(np.polyfit(x, y, 1)[0])


def normalized_auc(curve: pd.DataFrame, a0: float, a1: float) -> float:
    if curve.empty:
        return np.nan
    d = curve[(curve["crack_extension_um"] >= a0) & (curve["crack_extension_um"] <= a1)].copy()
    d = d.dropna()
    if len(d) < 2:
        return np.nan
    x = d["crack_extension_um"].to_numpy(float)
    y = d["K_MPa_sqrt_m"].to_numpy(float)
    order = np.argsort(x)
    x = x[order]
    y = y[order]
    if x[-1] <= x[0]:
        return np.nan
    return float(np.trapezoid(y, x) / (x[-1] - x[0]))


def sat_func(x, K0, dK, ell, p):
    x = np.asarray(x, float)
    ell = np.maximum(ell, 1e-9)
    p = np.maximum(p, 1e-6)
    return K0 + dK * (1.0 - np.exp(-np.power(np.maximum(x, 0) / ell, p)))


def power_j_func(x, J0, A, m):
    x = np.asarray(x, float)
    return J0 + A * np.power(np.maximum(x, 0), m)


def fit_saturating(binned: pd.DataFrame) -> dict:
    res = {
        "sat_success": False,
        "sat_K0": np.nan,
        "sat_dK": np.nan,
        "sat_Kss": np.nan,
        "sat_ell_um": np.nan,
        "sat_p": np.nan,
        "sat_rmse": np.nan,
    }
    if curve_fit is None or binned.empty or len(binned) < 5:
        return res
    x = binned["bin_center_um"].to_numpy(float)
    y = binned["K_median_MPa_sqrt_m"].to_numpy(float)
    ok = np.isfinite(x) & np.isfinite(y) & (y > 0)
    x, y = x[ok], y[ok]
    if len(x) < 5:
        return res
    K0_0 = float(np.nanmedian(y[:max(1, min(3, len(y)))]))
    Kss_0 = float(np.nanmedian(y[-max(1, min(3, len(y))):]))
    dK_0 = max(Kss_0 - K0_0, 0.1)
    ell_0 = max(0.25 * (np.nanmax(x) - np.nanmin(x)), 25.0)
    p_0 = 1.0
    try:
        popt, _ = curve_fit(
            sat_func, x, y,
            p0=[K0_0, dK_0, ell_0, p_0],
            bounds=([0.0, -200.0, 1.0, 0.2], [500.0, 500.0, 5000.0, 5.0]),
            maxfev=50000,
        )
        pred = sat_func(x, *popt)
        rmse = float(np.sqrt(np.mean((pred - y) ** 2)))
        K0, dK, ell, p = [float(v) for v in popt]
        res.update({
            "sat_success": True,
            "sat_K0": K0,
            "sat_dK": dK,
            "sat_Kss": K0 + dK,
            "sat_ell_um": ell,
            "sat_p": p,
            "sat_rmse": rmse,
        })
    except Exception:
        pass
    return res


def fit_power_j(binned: pd.DataFrame, Eprime_GPa: float) -> dict:
    res = {
        "J_power_success": False,
        "J0_kJ_m2": np.nan,
        "A_kJ_m2_per_um_m": np.nan,
        "m_power": np.nan,
        "J_power_rmse_kJ_m2": np.nan,
    }
    if curve_fit is None or binned.empty or len(binned) < 5:
        return res
    x = binned["bin_center_um"].to_numpy(float)
    K = binned["K_median_MPa_sqrt_m"].to_numpy(float)
    ok = np.isfinite(x) & np.isfinite(K) & (K > 0)
    x, K = x[ok], K[ok]
    if len(x) < 5:
        return res
    # K [MPa sqrt(m)] -> J [kJ/m^2], with E' in GPa:
    # K^2/E' = (MPa^2 m)/(GPa) = 1000 J/m^2 = 1 kJ/m^2.
    J = K ** 2 / Eprime_GPa
    J0_0 = max(float(np.nanmedian(J[:max(1, min(3, len(J)))])), 1e-6)
    A0 = max((float(np.nanmax(J)) - J0_0) / max(float(np.nanmax(x)) ** 0.5, 1), 1e-6)
    try:
        popt, _ = curve_fit(
            power_j_func, x, J,
            p0=[J0_0, A0, 0.5],
            bounds=([0.0, -1e4, 0.05], [1e5, 1e5, 3.0]),
            maxfev=50000,
        )
        pred = power_j_func(x, *popt)
        res.update({
            "J_power_success": True,
            "J0_kJ_m2": float(popt[0]),
            "A_kJ_m2_per_um_m": float(popt[1]),
            "m_power": float(popt[2]),
            "J_power_rmse_kJ_m2": float(np.sqrt(np.mean((pred - J) ** 2))),
        })
    except Exception:
        pass
    return res


def mean_curve_from_binned(binned_list: list[pd.DataFrame], bin_um: float, max_ext_um: float) -> pd.DataFrame:
    centers = np.arange(0.5 * bin_um, max_ext_um + 0.5 * bin_um, bin_um)
    rows = []
    for c in centers:
        vals = []
        for b in binned_list:
            if b.empty:
                continue
            idx = np.argmin(np.abs(b["bin_center_um"].to_numpy(float) - c))
            if abs(float(b["bin_center_um"].iloc[idx]) - c) < 0.51 * bin_um:
                vals.append(float(b["K_median_MPa_sqrt_m"].iloc[idx]))
        vals = np.array(vals, dtype=float)
        vals = vals[np.isfinite(vals)]
        if len(vals):
            rows.append({
                "bin_center_um": c,
                "K_mean_MPa_sqrt_m": float(np.mean(vals)),
                "K_median_across_seeds_MPa_sqrt_m": float(np.median(vals)),
                "K_std_MPa_sqrt_m": float(np.std(vals, ddof=1)) if len(vals) > 1 else 0.0,
                "K_p10_MPa_sqrt_m": float(np.quantile(vals, 0.10)),
                "K_p90_MPa_sqrt_m": float(np.quantile(vals, 0.90)),
                "n_seeds": int(len(vals)),
            })
    return pd.DataFrame(rows)


def publication_style():
    plt.rcParams.update({
        "font.size": 13,
        "axes.titlesize": 15,
        "axes.labelsize": 15,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "legend.fontsize": 10,
        "axes.linewidth": 1.2,
        "lines.linewidth": 2.2,
        "savefig.bbox": "tight",
    })


def plot_seed_analysis(klass: str, seed_curves: list[tuple[dict, pd.DataFrame, pd.DataFrame]],
                       mean_curve: pd.DataFrame, sat_fit: dict, out: Path, max_ext_um: float):
    publication_style()
    fig, ax = plt.subplots(figsize=(8.0, 5.6))
    cmap = plt.get_cmap("tab10")
    for i, (meta, raw, binned) in enumerate(seed_curves):
        if binned.empty:
            continue
        label = f"rep {int(meta['replicate']):02d}, seed {int(meta['seed'])}"
        if not meta["complete"]:
            label += " (incomplete)"
        ax.plot(
            binned["bin_center_um"], binned["K_median_MPa_sqrt_m"],
            color=cmap(i % 10),
            linestyle="-" if meta["complete"] else "--",
            alpha=0.85,
            marker="o",
            markersize=3.5,
            label=label,
        )

    if not mean_curve.empty:
        x = mean_curve["bin_center_um"].to_numpy(float)
        y = mean_curve["K_mean_MPa_sqrt_m"].to_numpy(float)
        s = mean_curve["K_std_MPa_sqrt_m"].to_numpy(float)
        ax.plot(x, y, color="black", linewidth=3.0, label="mean binned R-curve")
        ax.fill_between(x, y - s, y + s, color="black", alpha=0.12, linewidth=0, label="±1 SD")

    if sat_fit.get("sat_success") and not mean_curve.empty:
        xf = np.linspace(0, max_ext_um, 400)
        yf = sat_func(xf, sat_fit["sat_K0"], sat_fit["sat_dK"], sat_fit["sat_ell_um"], sat_fit["sat_p"])
        ax.plot(xf, yf, color="black", linewidth=2.0, linestyle=":", label="saturating fit")

    ax.set_xlim(0, max_ext_um)
    ax.set_xlabel("Crack extension, Δa (µm)")
    ax.set_ylabel(r"$K_J$ (MPa$\sqrt{m}$)")
    ax.set_title(f"{klass}: binned seed R-curves")
    ax.grid(False)
    ax.legend(frameon=False, ncol=1)
    fig.tight_layout()
    fig.savefig(out / f"{klass}_binned_seed_Rcurves_fit.png", dpi=320)
    fig.savefig(out / f"{klass}_binned_seed_Rcurves_fit.svg")
    plt.close(fig)


def plot_parameter_summary(seed_metrics: pd.DataFrame, out: Path):
    publication_style()
    complete = seed_metrics[seed_metrics["complete"]].copy()
    classes = [c for c in CLASSES_DEFAULT if c in set(complete["class"])]

    for metric, ylabel, fname in [
        ("K0_0_50_MPa_sqrt_m", r"$K_0$ 0–50 µm (MPa$\sqrt{m}$)", "summary_K0_early"),
        ("Kss_late_MPa_sqrt_m", r"$K_{ss}$ late-window (MPa$\sqrt{m}$)", "summary_Kss_late"),
        ("DeltaK_late_minus_early_MPa_sqrt_m", r"$\Delta K_R$ (MPa$\sqrt{m}$)", "summary_DeltaK"),
        ("AUC_0_max_MPa_sqrt_m", r"Normalized AUC (MPa$\sqrt{m}$)", "summary_AUC"),
        ("slope_100_400_MPa_sqrt_m_per_um", r"Slope 100–400 µm (MPa$\sqrt{m}$/µm)", "summary_slope"),
    ]:
        fig, ax = plt.subplots(figsize=(7.2, 5.0))
        x = np.arange(len(classes))
        means, stds, labels = [], [], []
        for klass in classes:
            vals = complete.loc[complete["class"] == klass, metric].dropna().to_numpy(float)
            means.append(float(np.nanmean(vals)) if len(vals) else np.nan)
            stds.append(float(np.nanstd(vals, ddof=1)) if len(vals) > 1 else 0.0)
            labels.append(f"{klass}\n(n={len(vals)})")
        ax.bar(x, means, yerr=stds, capsize=4, alpha=0.45)
        for i, klass in enumerate(classes):
            vals = complete.loc[complete["class"] == klass, metric].dropna().to_numpy(float)
            jitter = np.linspace(-0.10, 0.10, len(vals)) if len(vals) > 1 else np.zeros(len(vals))
            ax.scatter(np.full(len(vals), x[i]) + jitter, vals, s=48, zorder=3)
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        ax.set_ylabel(ylabel)
        ax.grid(False)
        fig.tight_layout()
        fig.savefig(out / f"{fname}.png", dpi=320)
        fig.savefig(out / f"{fname}.svg")
        plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="runs/four_class_exp_floor_CZM_500K_5rep_1000um_theta45")
    ap.add_argument("--out", default="runs/four_class_exp_floor_CZM_500K_5rep_1000um_theta45/Rcurve_analysis")
    ap.add_argument("--classes", default=" ".join(CLASSES_DEFAULT))
    ap.add_argument("--bin-um", type=float, default=25.0)
    ap.add_argument("--max-ext-um", type=float, default=1000.0)
    ap.add_argument("--initial-a-mm", type=float, default=0.5)
    ap.add_argument("--target-ext-um", type=float, default=1000.0)
    ap.add_argument("--E-GPa", type=float, default=410.0)
    ap.add_argument("--nu", type=float, default=0.28)
    ap.add_argument("--late-window", default="700 1000")
    args = ap.parse_args()

    root = Path(args.root)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    classes = args.classes.replace(",", " ").split()
    late0, late1 = [float(x) for x in args.late_window.replace(",", " ").split()]
    Eprime = args.E_GPa / (1.0 - args.nu ** 2)

    all_seed_metrics = []
    all_seed_binned = []
    class_fits = []

    for klass in classes:
        seed_curves = []
        binned_list = []

        for case in sorted((root / klass).glob("replicate_*_seed*/T500_th45")):
            if "geometry_veto" in str(case):
                continue
            meta = read_metadata(case, args.initial_a_mm, args.target_ext_um)
            raw = read_curve(case)
            binned = binned_curve(raw, args.bin_um, args.max_ext_um)
            meta["n_raw_points"] = len(raw)
            meta["n_bins"] = len(binned)

            if not raw.empty:
                K0 = window_median(raw, 0.0, 50.0)
                Klate = window_median(raw, late0, late1)
                Kmin = float(raw["K_MPa_sqrt_m"].min())
                amin = float(raw.loc[raw["K_MPa_sqrt_m"].idxmin(), "crack_extension_um"])
                AUC = normalized_auc(raw, 0.0, min(args.max_ext_um, float(raw["crack_extension_um"].max())))
                slope = window_slope(raw, 100.0, 400.0)
            else:
                K0 = Klate = Kmin = amin = AUC = slope = np.nan

            sat = fit_saturating(binned)
            power = fit_power_j(binned, Eprime_GPa=Eprime)

            row = {
                "class": klass,
                **meta,
                "K0_0_50_MPa_sqrt_m": K0,
                "Kss_late_MPa_sqrt_m": Klate,
                "DeltaK_late_minus_early_MPa_sqrt_m": Klate - K0 if np.isfinite(Klate) and np.isfinite(K0) else np.nan,
                "Kmin_MPa_sqrt_m": Kmin,
                "DeltaA_at_Kmin_um": amin,
                "Kss_minus_Kmin_MPa_sqrt_m": Klate - Kmin if np.isfinite(Klate) and np.isfinite(Kmin) else np.nan,
                "AUC_0_max_MPa_sqrt_m": AUC,
                "slope_100_400_MPa_sqrt_m_per_um": slope,
                "Eprime_GPa": Eprime,
                **sat,
                **power,
            }
            all_seed_metrics.append(row)

            if not binned.empty:
                bb = binned.copy()
                bb.insert(0, "class", klass)
                bb.insert(1, "replicate", meta["replicate"])
                bb.insert(2, "seed", meta["seed"])
                bb["complete"] = meta["complete"]
                all_seed_binned.append(bb)

            if not raw.empty:
                seed_curves.append((meta, raw, binned))
                binned_list.append(binned)

        mean_curve = mean_curve_from_binned(binned_list, args.bin_um, args.max_ext_um)
        mean_curve.to_csv(out / f"{klass}_mean_binned_Rcurve.csv", index=False)

        # Fit to class mean curve. Rename columns to fit_saturating expected schema.
        mean_for_fit = pd.DataFrame({
            "bin_center_um": mean_curve.get("bin_center_um", pd.Series(dtype=float)),
            "K_median_MPa_sqrt_m": mean_curve.get("K_mean_MPa_sqrt_m", pd.Series(dtype=float)),
        })
        class_sat = fit_saturating(mean_for_fit)
        class_power = fit_power_j(mean_for_fit, Eprime_GPa=Eprime)
        class_fits.append({"class": klass, **class_sat, **class_power})
        plot_seed_analysis(klass, seed_curves, mean_curve, class_sat, out, args.max_ext_um)

    seed_metrics = pd.DataFrame(all_seed_metrics).sort_values(["class", "replicate", "seed"])
    seed_metrics.to_csv(out / "seed_Rcurve_metrics_and_fits.csv", index=False)

    if all_seed_binned:
        pd.concat(all_seed_binned, ignore_index=True).to_csv(out / "seed_binned_Rcurves_long.csv", index=False)

    class_fit_df = pd.DataFrame(class_fits)
    class_fit_df.to_csv(out / "class_mean_Rcurve_fits.csv", index=False)

    complete = seed_metrics[seed_metrics["complete"]].copy()
    class_summary = (
        complete.groupby("class")
        .agg(
            n_complete=("seed", "count"),
            K0_mean=("K0_0_50_MPa_sqrt_m", "mean"),
            K0_std=("K0_0_50_MPa_sqrt_m", "std"),
            Kss_mean=("Kss_late_MPa_sqrt_m", "mean"),
            Kss_std=("Kss_late_MPa_sqrt_m", "std"),
            DeltaK_mean=("DeltaK_late_minus_early_MPa_sqrt_m", "mean"),
            DeltaK_std=("DeltaK_late_minus_early_MPa_sqrt_m", "std"),
            AUC_mean=("AUC_0_max_MPa_sqrt_m", "mean"),
            AUC_std=("AUC_0_max_MPa_sqrt_m", "std"),
            slope_100_400_mean=("slope_100_400_MPa_sqrt_m_per_um", "mean"),
            slope_100_400_std=("slope_100_400_MPa_sqrt_m_per_um", "std"),
            sat_K0_mean=("sat_K0", "mean"),
            sat_Kss_mean=("sat_Kss", "mean"),
            sat_ell_um_mean=("sat_ell_um", "mean"),
            sat_p_mean=("sat_p", "mean"),
            J0_kJ_m2_mean=("J0_kJ_m2", "mean"),
            J_power_m_mean=("m_power", "mean"),
        )
        .reset_index()
    )
    class_summary.to_csv(out / "class_Rcurve_metric_summary_complete_only.csv", index=False)
    plot_parameter_summary(seed_metrics, out)

    (out / "README_Rcurve_analysis.txt").write_text(
        "R-curve-like analysis of saved FEM/CZM seed outputs. No simulations are rerun.\\n\\n"
        "Key files:\\n"
        "  seed_Rcurve_metrics_and_fits.csv: per-seed window metrics and fit parameters.\\n"
        "  class_Rcurve_metric_summary_complete_only.csv: class mean/std over complete seeds.\\n"
        "  class_mean_Rcurve_fits.csv: saturating and J-power fits to the class mean binned R-curves.\\n"
        "  *_binned_seed_Rcurves_fit.png: seed binned R-curves, class mean, ±1 SD, and saturating fit.\\n\\n"
        "Metrics:\\n"
        "  K0_0_50: median KJ over 0-50 um.\\n"
        "  Kss_late: median KJ over late-window, default 700-1000 um.\\n"
        "  DeltaK: Kss_late - K0_0_50.\\n"
        "  AUC_0_max: normalized area under KJ(Δa).\\n"
        "  slope_100_400: linear slope of KJ over 100-400 um.\\n"
        "  sat_*: fit to K = K0 + dK*(1-exp(-(Δa/ell)^p)).\\n"
        "  J_power_*: fit to J = J0 + A*(Δa)^m after J=K^2/E'.\\n"
        "Incomplete seeds contribute to their own per-seed metrics where available, but class summary uses complete seeds only.\\n"
    )

    print(f"WROTE {out}")
    print(class_summary.to_string(index=False))


if __name__ == "__main__":
    main()
