#!/usr/bin/env python3
"""Compare PF sharp-front and FEM/CZM toughness summaries with analytical predictions.

This script now makes two plot families:
  1) first-passage fracture toughness vs. temperature;
  2) mean event-sampled propagation toughness over each run's available R-curve.

Incomplete runs are retained. They are shown as open markers and their shorter
final crack extension is recorded in the output tables.
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

CLASSES_DEFAULT = ["ceramic", "peak", "weakT", "DBTT"]
TEMPS_DEFAULT = list(range(300, 1201, 100))
A0_MM_DEFAULT = 0.5
TARGET_EXT_UM_DEFAULT = 500.0


def parse_list(text: str, cast=str):
    return [cast(x) for x in str(text).replace(",", " ").split() if x]


def first_existing(df: pd.DataFrame, candidates: Iterable[str]) -> str | None:
    lower = {c.lower(): c for c in df.columns}
    for c in candidates:
        if c in df.columns:
            return c
        if c.lower() in lower:
            return lower[c.lower()]
    return None


def normalize_class_name(x: object) -> str:
    s = str(x).strip()
    table = {
        "ceramic": "ceramic",
        "peak": "peak",
        "weakt": "weakT",
        "weakT": "weakT",
        "weak_t": "weakT",
        "weak-T": "weakT",
        "dbtt": "DBTT",
        "DBTT": "DBTT",
    }
    return table.get(s, table.get(s.lower(), s))


def read_summary_csv(root: Path, framework: str) -> pd.DataFrame | None:
    p = root / "four_class_temperature_summary.csv"
    if not p.exists():
        return None
    df = pd.read_csv(p)
    class_col = first_existing(df, ["class", "target_class", "fracture_class", "barrier_class"])
    T_col = first_existing(df, ["T_K", "temperature_K", "temperature", "T"])
    K_col = first_existing(df, [
        "Kc_first_MPa_sqrt_m", "Kc_first_MPa_sqrtm", "Kc_first",
        "Kinit_MPa_sqrt_m", "K_init_MPa_sqrt_m", "Kc_MPa_sqrtm",
        "Kc_MPa_sqrt_m", "KJ_first_MPa_sqrt_m",
    ])
    if class_col is None or T_col is None or K_col is None:
        raise SystemExit(
            f"{p} exists but could not identify required columns.\n"
            f"Found columns: {list(df.columns)}\n"
            f"Need class, temperature, and first-passage K columns."
        )
    out = pd.DataFrame({
        "framework": framework,
        "class": df[class_col].map(normalize_class_name),
        "T_K": pd.to_numeric(df[T_col], errors="coerce"),
        "Kc_first_MPa_sqrt_m": pd.to_numeric(df[K_col], errors="coerce"),
    })
    for c in ["status", "crack_extension_um", "extension_um", "n_growth_events", "n_r_curve_points", "case_dir"]:
        if c in df.columns:
            out[c] = df[c]
    if "crack_extension_um" not in out.columns and "extension_um" in out.columns:
        out["crack_extension_um"] = pd.to_numeric(out["extension_um"], errors="coerce")
    out["source_file"] = str(p)
    return out.dropna(subset=["T_K", "Kc_first_MPa_sqrt_m"]).reset_index(drop=True)


def parse_case_dir_temperature(name: str) -> int | None:
    m = re.search(r"T(\d{3,4})", name)
    return int(m.group(1)) if m else None


def summary_json_to_row(summary_path: Path, framework: str, klass: str,
                        a0_mm: float = A0_MM_DEFAULT) -> dict | None:
    try:
        data = json.loads(summary_path.read_text())
        if isinstance(data, list):
            data = data[0] if data else {}
    except Exception:
        return None
    T = data.get("T_K") or data.get("temperature_K") or parse_case_dir_temperature(summary_path.parent.name)
    K = None
    for key in ["Kc_first_MPa_sqrt_m", "Kc_first_MPa_sqrtm", "Kc_first", "Kinit_MPa_sqrt_m", "KJ_first_MPa_sqrt_m"]:
        if key in data and data[key] is not None:
            K = data[key]
            break
    if K is None:
        rcf = summary_path.parent / "R_curve_event_sampled.csv"
        if rcf.exists():
            try:
                rc = pd.read_csv(rcf)
                kcol = first_existing(rc, ["KJ_MPa_sqrt_m", "K_MPa_sqrt_m", "Kc_MPa_sqrt_m"])
                if len(rc) and kcol is not None:
                    K = float(rc[kcol].iloc[0])
            except Exception:
                pass
    if T is None or K is None:
        return None
    ext = None
    if "crack_extension_um" in data:
        ext = data.get("crack_extension_um")
    elif "a_final_mm" in data and data.get("a_final_mm") is not None:
        try:
            ext = (float(data["a_final_mm"]) - a0_mm) * 1000.0
        except Exception:
            ext = None
    return {
        "framework": framework,
        "class": normalize_class_name(klass),
        "T_K": float(T),
        "Kc_first_MPa_sqrt_m": float(K),
        "crack_extension_um": ext,
        "case_dir": str(summary_path.parent),
        "source_file": str(summary_path),
    }


def read_per_case_summaries(root: Path, framework: str, classes: list[str]) -> pd.DataFrame:
    rows = []
    for klass in classes:
        class_dir = root / klass
        if not class_dir.exists():
            continue
        for sf in sorted(class_dir.glob("T*_th*/summary.json")):
            row = summary_json_to_row(sf, framework, klass)
            if row is not None:
                rows.append(row)
    if not rows:
        return pd.DataFrame(columns=["framework", "class", "T_K", "Kc_first_MPa_sqrt_m"])
    return pd.DataFrame(rows)


def infer_case_dir(root: Path, klass: str, T_K: int | float) -> str | None:
    patt = root / klass
    if not patt.exists():
        return None
    cand = sorted(patt.glob(f"T{int(round(float(T_K)))}_th*"))
    if cand:
        return str(cand[0])
    return None


def load_framework(root: Path, framework: str, classes: list[str]) -> pd.DataFrame:
    if not root.exists():
        raise SystemExit(f"Missing {framework} root: {root}")
    df = read_summary_csv(root, framework)
    if df is None or df.empty:
        df = read_per_case_summaries(root, framework, classes)
    if df.empty:
        raise SystemExit(f"Could not find any first-passage data under {root}")
    if "case_dir" not in df.columns:
        df["case_dir"] = [infer_case_dir(root, k, t) for k, t in zip(df["class"], df["T_K"])]
    else:
        cd = df["case_dir"].astype(str)
        miss = (
            df["case_dir"].isna()
            | (cd.str.len() == 0)
            | (cd.str.lower() == "nan")
            | (~cd.map(lambda p: Path(p).exists()))
        )
        if miss.any():
            df.loc[miss, "case_dir"] = [
                infer_case_dir(root, k, t)
                for k, t in zip(df.loc[miss, "class"], df.loc[miss, "T_K"])
            ]
    return df


def load_analytic(path: Path, column: str) -> pd.DataFrame:
    if not path.exists():
        raise SystemExit(f"Missing analytical prediction CSV: {path}")
    a = pd.read_csv(path)
    class_col = first_existing(a, ["class", "target_class", "fracture_class", "barrier_class"])
    T_col = first_existing(a, ["T_K", "temperature_K", "temperature", "T"])
    K_col = first_existing(a, [column, "K_analytic_MPa_sqrt_m", "Kc_MPa_sqrtm", "Kc_MPa_sqrt_m", "target_Kc_MPa_sqrtm"])
    if class_col is None or T_col is None or K_col is None:
        raise SystemExit(
            f"Could not identify analytical columns in {path}.\n"
            f"Found columns: {list(a.columns)}\n"
            f"Requested analytical column: {column}"
        )
    out = pd.DataFrame({
        "class": a[class_col].map(normalize_class_name),
        "T_K": pd.to_numeric(a[T_col], errors="coerce"),
        "K_analytic_MPa_sqrt_m": pd.to_numeric(a[K_col], errors="coerce"),
    })
    for c in ["K_target_MPa_sqrt_m", "K_prior_weakT_MPa_sqrt_m", "analytic_source"]:
        if c in a.columns:
            out[c] = a[c]
    return out.dropna(subset=["T_K", "K_analytic_MPa_sqrt_m"]).reset_index(drop=True)


def nearest_interp(analytic: pd.DataFrame, klass: str, Tvals: np.ndarray) -> np.ndarray:
    g = analytic[analytic["class"] == klass].sort_values("T_K")
    if g.empty:
        return np.full(len(Tvals), np.nan)
    return np.interp(Tvals.astype(float), g.T_K.to_numpy(float), g.K_analytic_MPa_sqrt_m.to_numpy(float), left=np.nan, right=np.nan)


def extract_event_curve(case_dir: Path) -> pd.DataFrame:
    rcf = case_dir / "R_curve_event_sampled.csv"
    if rcf.exists():
        df = pd.read_csv(rcf)
    else:
        # fallback to raw steps file
        steps = sorted(case_dir.glob("steps_*K.csv"))
        if not steps:
            return pd.DataFrame()
        df = pd.read_csv(steps[0])
        da_col = first_existing(df, ["da_block_um", "da_block_uM", "da_block"])
        nfire_col = first_existing(df, ["n_fire", "nfire"])
        keep = np.zeros(len(df), dtype=bool)
        if da_col is not None:
            keep |= pd.to_numeric(df[da_col], errors="coerce").fillna(0).to_numpy(float) > 0.0
        if nfire_col is not None:
            keep |= pd.to_numeric(df[nfire_col], errors="coerce").fillna(0).to_numpy(float) > 0.0
        df = df.loc[keep].copy()
    if df.empty:
        return df
    kcol = first_existing(df, ["KJ_MPa_sqrt_m", "K_MPa_sqrt_m", "Kc_MPa_sqrt_m", "K_MPa_sqrtm"])
    xcol = first_existing(df, ["crack_extension_um", "extension_um", "delta_a_um", "da_cumulative_um"])
    if kcol is None:
        return pd.DataFrame()
    out = pd.DataFrame({"KJ_MPa_sqrt_m": pd.to_numeric(df[kcol], errors="coerce")})
    if xcol is not None:
        out["crack_extension_um"] = pd.to_numeric(df[xcol], errors="coerce")
    else:
        out["crack_extension_um"] = np.arange(1, len(out)+1, dtype=float)
    return out.dropna(subset=["KJ_MPa_sqrt_m", "crack_extension_um"]).reset_index(drop=True)


def add_rcurve_metrics(data: pd.DataFrame, target_ext_um: float) -> pd.DataFrame:
    rows = []
    for _, r in data.iterrows():
        d = r.to_dict()
        case_dir = d.get("case_dir")
        rc = pd.DataFrame()
        if case_dir and str(case_dir) != "nan":
            try:
                rc = extract_event_curve(Path(case_dir))
            except Exception:
                rc = pd.DataFrame()
        if not rc.empty:
            d["n_r_curve_points"] = int(len(rc))
            d["K_mean_all_points_MPa_sqrt_m"] = float(rc["KJ_MPa_sqrt_m"].mean())
            d["K_median_all_points_MPa_sqrt_m"] = float(rc["KJ_MPa_sqrt_m"].median())
            d["K_std_all_points_MPa_sqrt_m"] = float(rc["KJ_MPa_sqrt_m"].std(ddof=0))
            d["K_first_rcurve_MPa_sqrt_m"] = float(rc["KJ_MPa_sqrt_m"].iloc[0])
            d["K_last_rcurve_MPa_sqrt_m"] = float(rc["KJ_MPa_sqrt_m"].iloc[-1])
            d["crack_extension_um"] = float(np.nanmax(rc["crack_extension_um"].to_numpy(float)))
        else:
            d.setdefault("n_r_curve_points", np.nan)
            d["K_mean_all_points_MPa_sqrt_m"] = np.nan
            d["K_median_all_points_MPa_sqrt_m"] = np.nan
            d["K_std_all_points_MPa_sqrt_m"] = np.nan
            d["K_first_rcurve_MPa_sqrt_m"] = np.nan
            d["K_last_rcurve_MPa_sqrt_m"] = np.nan
        ext = d.get("crack_extension_um", np.nan)
        try:
            ext = float(ext)
        except Exception:
            ext = np.nan
        d["crack_extension_um"] = ext
        d["is_complete"] = bool(np.isfinite(ext) and ext >= 0.98 * float(target_ext_um))
        rows.append(d)
    return pd.DataFrame(rows)


def publication_style():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "font.size": 14,
        "axes.titlesize": 16,
        "axes.labelsize": 16,
        "xtick.labelsize": 13,
        "ytick.labelsize": 13,
        "legend.fontsize": 11,
        "axes.linewidth": 1.2,
        "lines.linewidth": 2.4,
        "lines.markersize": 7,
        "grid.linewidth": 0.7,
        "grid.alpha": 0.28,
        "savefig.bbox": "tight",
        "figure.dpi": 150,
    })
    return plt


FRAMEWORK_STYLES = {
    "PF sharp-front": {"color": "#0072B2", "marker": "o", "linestyle": (0, (5, 2))},
    "FEM/CZM": {"color": "#D55E00", "marker": "s", "linestyle": (0, (1, 1))},
}
ANALYTIC_STYLE = {"color": "#111111", "linestyle": "-", "linewidth": 2.8}


def make_overlay_plots(data: pd.DataFrame, analytic: pd.DataFrame, out: Path,
                       classes: list[str], metric_col: str, y_label: str,
                       fig_stub: str, class_stub_prefix: str,
                       analytic_label: str) -> None:
    plt = publication_style()
    from matplotlib.lines import Line2D

    def draw_panel(ax, klass: str):
        ga = analytic[analytic["class"] == klass].sort_values("T_K")
        if not ga.empty:
            ax.plot(ga.T_K, ga.K_analytic_MPa_sqrt_m, label=analytic_label, **ANALYTIC_STYLE)
        for fw in ["PF sharp-front", "FEM/CZM"]:
            style = FRAMEWORK_STYLES[fw]
            g = data[(data["class"] == klass) & (data["framework"] == fw)].sort_values("T_K")
            if g.empty:
                continue
            ym = pd.to_numeric(g[metric_col], errors="coerce")
            ok = np.isfinite(ym.to_numpy(float))
            g = g.loc[ok].copy()
            if g.empty:
                continue
            ym = ym.loc[ok]
            ax.plot(g.T_K, ym, color=style["color"], linestyle=style["linestyle"], linewidth=2.2, zorder=2)
            comp = g["is_complete"].fillna(False).to_numpy(bool)
            # completed: filled markers
            if comp.any():
                gg = g.loc[comp]
                yy = pd.to_numeric(gg[metric_col], errors="coerce")
                ax.scatter(gg.T_K, yy, marker=style["marker"], s=55,
                           facecolor=style["color"], edgecolor="white", linewidth=0.8, zorder=3)
            if (~comp).any():
                gg = g.loc[~comp]
                yy = pd.to_numeric(gg[metric_col], errors="coerce")
                ax.scatter(gg.T_K, yy, marker=style["marker"], s=55,
                           facecolor="white", edgecolor=style["color"], linewidth=1.6, zorder=3)
        ax.set_title(klass)
        ax.set_ylabel(y_label)
        ax.grid(False)
        ax.set_xlim(min(TEMPS_DEFAULT)-25, max(TEMPS_DEFAULT)+25)

    fig, axes = plt.subplots(2, 2, figsize=(12.8, 9.1), sharex=True)
    axes = axes.ravel()
    for ax, klass in zip(axes, classes):
        draw_panel(ax, klass)
    for ax in axes[-2:]:
        ax.set_xlabel("Temperature (K)")

    legend_items = [
        Line2D([0], [0], **ANALYTIC_STYLE, label=analytic_label),
        Line2D([0], [0], color=FRAMEWORK_STYLES["PF sharp-front"]["color"], linestyle=FRAMEWORK_STYLES["PF sharp-front"]["linestyle"], marker=FRAMEWORK_STYLES["PF sharp-front"]["marker"], markerfacecolor=FRAMEWORK_STYLES["PF sharp-front"]["color"], markeredgecolor="white", label="PF sharp-front (complete)"),
        Line2D([0], [0], color=FRAMEWORK_STYLES["FEM/CZM"]["color"], linestyle=FRAMEWORK_STYLES["FEM/CZM"]["linestyle"], marker=FRAMEWORK_STYLES["FEM/CZM"]["marker"], markerfacecolor=FRAMEWORK_STYLES["FEM/CZM"]["color"], markeredgecolor="white", label="FEM/CZM (complete)"),
        Line2D([0], [0], color="#666666", linestyle="None", marker="o", markerfacecolor="white", markeredgecolor="#666666", label="Open marker = incomplete run"),
    ]
    fig.legend(handles=legend_items, loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 1.02))
    fig.savefig(out / f"{fig_stub}.png", dpi=320)
    fig.savefig(out / f"{fig_stub}.svg")
    plt.close(fig)

    for klass in classes:
        fig, ax = plt.subplots(figsize=(7.4, 5.4))
        draw_panel(ax, klass)
        ax.set_xlabel("Temperature (K)")
        ax.legend(handles=legend_items, frameon=False, loc="best")
        safe = klass.replace("/", "_")
        fig.savefig(out / f"{class_stub_prefix}_{safe}.png", dpi=320)
        fig.savefig(out / f"{class_stub_prefix}_{safe}.svg")
        plt.close(fig)


def make_error_tables_first_passage(data: pd.DataFrame, analytic: pd.DataFrame, out: Path) -> pd.DataFrame:
    rows = []
    for (fw, klass), g in data.groupby(["framework", "class"]):
        T = g.T_K.to_numpy(float)
        Ka = nearest_interp(analytic, klass, T)
        for (_, r), ka in zip(g.iterrows(), Ka):
            err = float(r.Kc_first_MPa_sqrt_m) - float(ka) if np.isfinite(ka) else np.nan
            rel = err / float(ka) if np.isfinite(ka) and abs(float(ka)) > 1e-12 else np.nan
            d = r.to_dict()
            d.update({
                "K_analytic_interp_MPa_sqrt_m": ka,
                "error_vs_analytic_MPa_sqrt_m": err,
                "relative_error_vs_analytic": rel,
                "percent_error_vs_analytic": 100.0 * rel if np.isfinite(rel) else np.nan,
            })
            rows.append(d)
    outdf = pd.DataFrame(rows).sort_values(["class", "T_K", "framework"])
    outdf.to_csv(out / "first_passage_comparison_with_analytic.csv", index=False)
    metric = (outdf.groupby(["framework", "class"], dropna=False)
              .agg(n=("T_K", "count"),
                   RMSE_MPa_sqrt_m=("error_vs_analytic_MPa_sqrt_m", lambda x: float(np.sqrt(np.nanmean(np.asarray(x, float)**2))) if len(x) else np.nan),
                   mean_abs_percent_error=("percent_error_vs_analytic", lambda x: float(np.nanmean(np.abs(np.asarray(x, float)))) if len(x) else np.nan),
                   bias_MPa_sqrt_m=("error_vs_analytic_MPa_sqrt_m", lambda x: float(np.nanmean(np.asarray(x, float))) if len(x) else np.nan))
              .reset_index())
    metric.to_csv(out / "first_passage_error_metrics_by_class.csv", index=False)
    return outdf


def make_mean_rcurve_tables(data: pd.DataFrame, analytic: pd.DataFrame, out: Path) -> pd.DataFrame:
    rows = []
    valid = data[np.isfinite(pd.to_numeric(data["K_mean_all_points_MPa_sqrt_m"], errors="coerce"))].copy()
    for (fw, klass), g in valid.groupby(["framework", "class"]):
        T = g.T_K.to_numpy(float)
        Ka = nearest_interp(analytic, klass, T)
        for (_, r), ka in zip(g.iterrows(), Ka):
            kval = float(r.K_mean_all_points_MPa_sqrt_m)
            err = kval - float(ka) if np.isfinite(ka) else np.nan
            rel = err / float(ka) if np.isfinite(ka) and abs(float(ka)) > 1e-12 else np.nan
            d = r.to_dict()
            d.update({
                "K_mean_event_interp_reference_MPa_sqrt_m": ka,
                "meanK_minus_analytic_MPa_sqrt_m": err,
                "relative_meanK_vs_analytic": rel,
                "percent_meanK_vs_analytic": 100.0 * rel if np.isfinite(rel) else np.nan,
            })
            rows.append(d)
    if not rows:
        cols = [
            "framework", "class", "T_K", "K_mean_all_points_MPa_sqrt_m",
            "K_mean_event_interp_reference_MPa_sqrt_m",
            "meanK_minus_analytic_MPa_sqrt_m",
            "relative_meanK_vs_analytic", "percent_meanK_vs_analytic",
            "crack_extension_um", "n_r_curve_points", "is_complete", "case_dir",
        ]
        outdf = pd.DataFrame(columns=cols)
        outdf.to_csv(out / "mean_available_rcurve_comparison_with_analytic.csv", index=False)
        pd.DataFrame(columns=[
            "framework", "class", "n", "RMSE_MPa_sqrt_m",
            "mean_abs_percent_error", "bias_MPa_sqrt_m",
            "mean_final_extension_um",
        ]).to_csv(out / "mean_available_rcurve_error_metrics_by_class.csv", index=False)
        print("WARNING: no event-sampled R-curve files were found; mean-R-curve tables are empty.")
        return outdf

    outdf = pd.DataFrame(rows).sort_values(["class", "T_K", "framework"])
    outdf.to_csv(out / "mean_available_rcurve_comparison_with_analytic.csv", index=False)
    metric = (outdf.groupby(["framework", "class"], dropna=False)
              .agg(n=("T_K", "count"),
                   RMSE_MPa_sqrt_m=("meanK_minus_analytic_MPa_sqrt_m", lambda x: float(np.sqrt(np.nanmean(np.asarray(x, float)**2))) if len(x) else np.nan),
                   mean_abs_percent_error=("percent_meanK_vs_analytic", lambda x: float(np.nanmean(np.abs(np.asarray(x, float)))) if len(x) else np.nan),
                   bias_MPa_sqrt_m=("meanK_minus_analytic_MPa_sqrt_m", lambda x: float(np.nanmean(np.asarray(x, float))) if len(x) else np.nan),
                   mean_final_extension_um=("crack_extension_um", lambda x: float(np.nanmean(np.asarray(x, float))) if len(x) else np.nan))
              .reset_index())
    metric.to_csv(out / "mean_available_rcurve_error_metrics_by_class.csv", index=False)
    return outdf



def make_median_rcurve_tables(data: pd.DataFrame, analytic: pd.DataFrame, out: Path) -> pd.DataFrame:
    rows = []
    valid = data[np.isfinite(pd.to_numeric(data["K_median_all_points_MPa_sqrt_m"], errors="coerce"))].copy()
    for (fw, klass), g in valid.groupby(["framework", "class"]):
        T = g.T_K.to_numpy(float)
        Ka = nearest_interp(analytic, klass, T)
        for (_, r), ka in zip(g.iterrows(), Ka):
            kval = float(r.K_median_all_points_MPa_sqrt_m)
            err = kval - float(ka) if np.isfinite(ka) else np.nan
            rel = err / float(ka) if np.isfinite(ka) and abs(float(ka)) > 1e-12 else np.nan
            d = r.to_dict()
            d.update({
                "K_median_event_interp_reference_MPa_sqrt_m": ka,
                "medianK_minus_analytic_MPa_sqrt_m": err,
                "relative_medianK_vs_analytic": rel,
                "percent_medianK_vs_analytic": 100.0 * rel if np.isfinite(rel) else np.nan,
            })
            rows.append(d)
    if not rows:
        cols = [
            "framework", "class", "T_K", "K_median_all_points_MPa_sqrt_m",
            "K_median_event_interp_reference_MPa_sqrt_m",
            "medianK_minus_analytic_MPa_sqrt_m",
            "relative_medianK_vs_analytic", "percent_medianK_vs_analytic",
            "crack_extension_um", "n_r_curve_points", "is_complete", "case_dir",
        ]
        outdf = pd.DataFrame(columns=cols)
        outdf.to_csv(out / "median_available_rcurve_comparison_with_analytic.csv", index=False)
        pd.DataFrame(columns=[
            "framework", "class", "n", "RMSE_MPa_sqrt_m",
            "mean_abs_percent_error", "bias_MPa_sqrt_m",
            "mean_final_extension_um",
        ]).to_csv(out / "median_available_rcurve_error_metrics_by_class.csv", index=False)
        print("WARNING: no event-sampled R-curve files were found; median-R-curve tables are empty.")
        return outdf

    outdf = pd.DataFrame(rows).sort_values(["class", "T_K", "framework"])
    outdf.to_csv(out / "median_available_rcurve_comparison_with_analytic.csv", index=False)
    metric = (outdf.groupby(["framework", "class"], dropna=False)
              .agg(n=("T_K", "count"),
                   RMSE_MPa_sqrt_m=("medianK_minus_analytic_MPa_sqrt_m", lambda x: float(np.sqrt(np.nanmean(np.asarray(x, float)**2))) if len(x) else np.nan),
                   mean_abs_percent_error=("percent_medianK_vs_analytic", lambda x: float(np.nanmean(np.abs(np.asarray(x, float)))) if len(x) else np.nan),
                   bias_MPa_sqrt_m=("medianK_minus_analytic_MPa_sqrt_m", lambda x: float(np.nanmean(np.asarray(x, float))) if len(x) else np.nan),
                   mean_final_extension_um=("crack_extension_um", lambda x: float(np.nanmean(np.asarray(x, float))) if len(x) else np.nan))
              .reset_index())
    metric.to_csv(out / "median_available_rcurve_error_metrics_by_class.csv", index=False)
    return outdf

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pf-root", default="PF-four_class_exp_floor_PF_sharp_no_branch_500um_theta45")
    ap.add_argument("--czm-root", default="four_class_exp_floor_CZM_rates_no_branch_500um_theta45/rate_1x")
    ap.add_argument("--analytic-csv", default="four_class_analytical_prediction_final.csv")
    ap.add_argument("--analytic-column", default="K_analytic_MPa_sqrt_m",
                    help="Column to use from the analytical CSV. Use K_target_MPa_sqrt_m to overlay the original tuning target instead.")
    ap.add_argument("--out", default="PF_vs_CZM_first_passage_with_analytic")
    ap.add_argument("--classes", default=" ".join(CLASSES_DEFAULT))
    ap.add_argument("--temps", default=" ".join(map(str, TEMPS_DEFAULT)))
    ap.add_argument("--target-ext-um", type=float, default=TARGET_EXT_UM_DEFAULT)
    args = ap.parse_args()

    classes = parse_list(args.classes, str)
    temps = parse_list(args.temps, int)
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    pfroot = Path(args.pf_root)
    czroot = Path(args.czm_root)

    pf = load_framework(pfroot, "PF sharp-front", classes)
    cz = load_framework(czroot, "FEM/CZM", classes)
    data = pd.concat([pf, cz], ignore_index=True, sort=False)
    data = data[data["class"].isin(classes) & data["T_K"].isin(temps)].copy()
    data = add_rcurve_metrics(data, args.target_ext_um)
    if "n_r_curve_points" not in data.columns or not np.isfinite(pd.to_numeric(data["n_r_curve_points"], errors="coerce")).any():
        print("WARNING: no R_curve_event_sampled.csv files were found. "
              "First-passage plots will be produced; mean-R-curve plots may be empty.")
    data = data.sort_values(["class", "T_K", "framework"])
    data.to_csv(out / "combined_first_passage_summary.csv", index=False)
    coverage = data[["framework", "class", "T_K", "crack_extension_um", "n_r_curve_points", "is_complete", "case_dir"]].copy()
    coverage.to_csv(out / "r_curve_coverage_summary.csv", index=False)

    analytic = load_analytic(Path(args.analytic_csv), args.analytic_column)
    analytic = analytic[analytic["class"].isin(classes)].copy()
    analytic.to_csv(out / "analytical_prediction_used.csv", index=False)

    make_overlay_plots(
        data, analytic, out, classes,
        metric_col="Kc_first_MPa_sqrt_m",
        y_label=r"First-passage $K_c$ (MPa$\sqrt{m}$)",
        fig_stub="PF_CZM_first_passage_vs_analytic_publication",
        class_stub_prefix="first_passage_publication",
        analytic_label="V1 analytical prediction",
    )
    make_overlay_plots(
        data, analytic, out, classes,
        metric_col="K_mean_all_points_MPa_sqrt_m",
        y_label=r"Mean event-sampled $K_J$ (MPa$\sqrt{m}$)",
        fig_stub="PF_CZM_mean_available_rcurve_vs_analytic_publication",
        class_stub_prefix="mean_available_rcurve_publication",
        analytic_label="V1 analytical prediction (reference)",
    )
    make_overlay_plots(
        data, analytic, out, classes,
        metric_col="K_median_all_points_MPa_sqrt_m",
        y_label=r"Median event-sampled $K_J$ (MPa$\sqrt{m}$)",
        fig_stub="PF_CZM_median_available_rcurve_vs_analytic_publication",
        class_stub_prefix="median_available_rcurve_publication",
        analytic_label="V1 analytical prediction (reference)",
    )
    make_error_tables_first_passage(data, analytic, out)
    make_mean_rcurve_tables(data, analytic, out)
    make_median_rcurve_tables(data, analytic, out)

    cfg = {
        "pf_root": str(pfroot),
        "czm_root": str(czroot),
        "analytic_csv": str(args.analytic_csv),
        "analytic_column": args.analytic_column,
        "classes": classes,
        "temps_K": temps,
        "target_extension_um": args.target_ext_um,
        "outputs": [
            "PF_CZM_first_passage_vs_analytic_publication.png",
            "PF_CZM_mean_available_rcurve_vs_analytic_publication.png",
            "PF_CZM_median_available_rcurve_vs_analytic_publication.png",
            "first_passage_publication_<class>.png",
            "mean_available_rcurve_publication_<class>.png",
            "median_available_rcurve_publication_<class>.png",
            "combined_first_passage_summary.csv",
            "r_curve_coverage_summary.csv",
            "first_passage_comparison_with_analytic.csv",
            "mean_available_rcurve_comparison_with_analytic.csv",
        ],
    }
    (out / "comparison_config.json").write_text(json.dumps(cfg, indent=2))
    print(f"WROTE {out}")
    print(f"  first-passage plot: {out / 'PF_CZM_first_passage_vs_analytic_publication.png'}")
    print(f"  mean-R-curve plot: {out / 'PF_CZM_mean_available_rcurve_vs_analytic_publication.png'}")


if __name__ == "__main__":
    main()
