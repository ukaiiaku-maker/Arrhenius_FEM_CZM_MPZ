#!/usr/bin/env python3
"""Plot FEM/CZM four-class fracture toughness at three rates with V1 analytical predictions.

Outputs two figure families:
  1. First-passage Kc(T)
  2. Median event-sampled KJ(T) over the available crack-growth history

The script reads rate folders such as:
  ROOT/rate_1x
  ROOT/rate_10x
  ROOT/rate_100x

For analytical predictions, it first tries to compute rate-specific V1 curves by importing
run_v1_exp_floor_four_class_tuning.py and using four_class_exp_floor_exact_model_inputs.csv.
If that fails, it falls back to an analytical CSV. If the fallback CSV has no rate_factor
column, it is reused for all rates with a warning.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import math
import re
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

CLASSES_DEFAULT = ["ceramic", "peak", "weakT", "DBTT"]
RATES_DEFAULT = [1.0, 10.0, 100.0]
TEMPS_DEFAULT = list(range(300, 1201, 100))
A0_MM_DEFAULT = 0.5


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


def rate_label(rate: float) -> str:
    r = float(rate)
    if abs(r - round(r)) < 1e-12:
        return f"{int(round(r))}x"
    return f"{r:g}x"


def rate_folder(root: Path, rate: float) -> Path:
    candidates = [
        root / f"rate_{rate_label(rate)}",
        root / f"rate_{int(round(rate))}x",
        root / f"rate_{rate:g}x",
    ]
    for c in candidates:
        if c.exists():
            return c
    return candidates[0]


def parse_case_dir_temperature(name: str) -> int | None:
    m = re.search(r"T(\d{3,4})", name)
    return int(m.group(1)) if m else None


def infer_case_dir(root: Path, klass: str, T_K: int | float) -> str | None:
    class_dir = root / klass
    if not class_dir.exists():
        return None
    cand = sorted(class_dir.glob(f"T{int(round(float(T_K)))}_th*"))
    if cand:
        return str(cand[0])
    return None


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
    for key in [
        "Kc_first_MPa_sqrt_m", "Kc_first_MPa_sqrtm", "Kc_first",
        "Kinit_MPa_sqrt_m", "KJ_first_MPa_sqrt_m", "Kc_MPa_sqrtm",
    ]:
        if key in data and data[key] is not None:
            K = data[key]
            break
    if K is None:
        rc = summary_path.parent / "R_curve_event_sampled.csv"
        if rc.exists():
            try:
                dfr = pd.read_csv(rc)
                kcol = first_existing(dfr, ["KJ_MPa_sqrt_m", "K_MPa_sqrt_m", "Kc_MPa_sqrt_m"])
                if len(dfr) and kcol:
                    K = float(dfr[kcol].iloc[0])
            except Exception:
                pass
    if T is None or K is None:
        return None

    ext = None
    for key in ["crack_extension_um", "final_crack_extension_um", "extension_um"]:
        if key in data and data[key] is not None:
            ext = data[key]
            break
    if ext is None and "a_final_mm" in data and data.get("a_final_mm") is not None:
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


def read_summary_csv(rate_root: Path, framework: str, classes: list[str]) -> pd.DataFrame | None:
    p = rate_root / "four_class_temperature_summary.csv"
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
            f"{p} exists but required columns were not identified.\n"
            f"Found columns: {list(df.columns)}"
        )
    out = pd.DataFrame({
        "framework": framework,
        "class": df[class_col].map(normalize_class_name),
        "T_K": pd.to_numeric(df[T_col], errors="coerce"),
        "Kc_first_MPa_sqrt_m": pd.to_numeric(df[K_col], errors="coerce"),
    })
    for c in [
        "status", "crack_extension_um", "final_crack_extension_um", "extension_um",
        "n_growth_events", "n_r_curve_points", "case_dir",
    ]:
        if c in df.columns:
            out[c] = df[c]
    if "crack_extension_um" not in out.columns:
        for c in ["final_crack_extension_um", "extension_um"]:
            if c in out.columns:
                out["crack_extension_um"] = pd.to_numeric(out[c], errors="coerce")
                break
    if "case_dir" not in out.columns:
        out["case_dir"] = [infer_case_dir(rate_root, k, t) for k, t in zip(out["class"], out["T_K"])]
    else:
        cd = out["case_dir"].astype(str)
        miss = (
            out["case_dir"].isna()
            | (cd.str.len() == 0)
            | (cd.str.lower() == "nan")
            | (~cd.map(lambda pth: Path(pth).exists()))
        )
        if miss.any():
            out.loc[miss, "case_dir"] = [
                infer_case_dir(rate_root, k, t)
                for k, t in zip(out.loc[miss, "class"], out.loc[miss, "T_K"])
            ]
    return out.dropna(subset=["T_K", "Kc_first_MPa_sqrt_m"]).reset_index(drop=True)


def read_per_case(rate_root: Path, framework: str, classes: list[str]) -> pd.DataFrame:
    rows = []
    for klass in classes:
        cdir = rate_root / klass
        if not cdir.exists():
            continue
        for sf in sorted(cdir.glob("T*_th*/summary.json")):
            r = summary_json_to_row(sf, framework, klass)
            if r is not None:
                rows.append(r)
    return pd.DataFrame(rows) if rows else pd.DataFrame()


def extract_event_curve(case_dir: Path) -> pd.DataFrame:
    rcf = case_dir / "R_curve_event_sampled.csv"
    if rcf.exists():
        df = pd.read_csv(rcf)
    else:
        steps = sorted(case_dir.glob("steps_*K.csv"))
        if not steps:
            return pd.DataFrame()
        df = pd.read_csv(steps[0])
        da_col = first_existing(df, ["da_block_um", "da_block_uM", "da_block"])
        nfire_col = first_existing(df, ["n_fire", "nfire"])
        keep = np.zeros(len(df), dtype=bool)
        if da_col:
            keep |= pd.to_numeric(df[da_col], errors="coerce").fillna(0).to_numpy(float) > 0
        if nfire_col:
            keep |= pd.to_numeric(df[nfire_col], errors="coerce").fillna(0).to_numpy(float) > 0
        df = df.loc[keep].copy()
    if df.empty:
        return pd.DataFrame()
    kcol = first_existing(df, ["KJ_MPa_sqrt_m", "K_MPa_sqrt_m", "Kc_MPa_sqrt_m", "K_MPa_sqrtm"])
    xcol = first_existing(df, ["crack_extension_um", "extension_um", "delta_a_um", "da_cumulative_um"])
    if kcol is None:
        return pd.DataFrame()
    out = pd.DataFrame({"KJ_MPa_sqrt_m": pd.to_numeric(df[kcol], errors="coerce")})
    if xcol:
        out["crack_extension_um"] = pd.to_numeric(df[xcol], errors="coerce")
    else:
        out["crack_extension_um"] = np.arange(1, len(out) + 1, dtype=float)
    return out.dropna(subset=["KJ_MPa_sqrt_m", "crack_extension_um"]).reset_index(drop=True)


def load_rate_data(root: Path, rates: list[float], classes: list[str], target_ext_um: float) -> pd.DataFrame:
    """Load rate data by scanning both top-level summary CSVs and per-case folders.

    Earlier summary CSVs can be stale or incomplete after restart runs.  The
    per-case summaries are therefore always scanned and merged with the summary
    table.  When both exist for the same (class, T), the per-case record wins
    because it is closest to the actual case directory and R-curve file.
    """
    rows = []
    availability = []
    for r in rates:
        rr = rate_folder(root, r)
        if not rr.exists():
            print(f"WARNING: missing rate folder {rr}; skipping")
            continue

        frames = []
        sdf = read_summary_csv(rr, "FEM/CZM", classes)
        if sdf is not None and not sdf.empty:
            sdf = sdf.copy()
            sdf["record_source"] = "summary_csv"
            sdf["source_priority"] = 0
            frames.append(sdf)
        pdf = read_per_case(rr, "FEM/CZM", classes)
        if pdf is not None and not pdf.empty:
            pdf = pdf.copy()
            pdf["record_source"] = "per_case_summary"
            pdf["source_priority"] = 1
            frames.append(pdf)

        if not frames:
            print(f"WARNING: no summary or per-case data found under {rr}")
            continue

        df = pd.concat(frames, ignore_index=True, sort=False)
        df["class"] = df["class"].map(normalize_class_name)
        df["T_K"] = pd.to_numeric(df["T_K"], errors="coerce")
        df = df.dropna(subset=["class", "T_K", "Kc_first_MPa_sqrt_m"])
        # Prefer per-case rows over stale top-level summary rows for the same case.
        df = (df.sort_values(["class", "T_K", "source_priority"])
                .drop_duplicates(subset=["class", "T_K"], keep="last")
                .reset_index(drop=True))
        df["rate_factor"] = float(r)
        df["rate_label"] = rate_label(r)

        out_rows = []
        for _, row in df.iterrows():
            d = row.to_dict()
            case_dir = d.get("case_dir")
            # If a stale case_dir persisted from an old summary, infer the actual path.
            if not case_dir or str(case_dir).lower() == "nan" or not Path(str(case_dir)).exists():
                case_dir = infer_case_dir(rr, d.get("class"), d.get("T_K"))
                d["case_dir"] = case_dir
            rc = pd.DataFrame()
            if case_dir and str(case_dir).lower() != "nan":
                try:
                    rc = extract_event_curve(Path(case_dir))
                except Exception as e:
                    print(f"WARNING: failed reading R-curve for {case_dir}: {e}")
            if not rc.empty:
                d["n_r_curve_points"] = int(len(rc))
                d["K_mean_all_points_MPa_sqrt_m"] = float(rc["KJ_MPa_sqrt_m"].mean())
                d["K_median_all_points_MPa_sqrt_m"] = float(rc["KJ_MPa_sqrt_m"].median())
                d["K_p10_all_points_MPa_sqrt_m"] = float(rc["KJ_MPa_sqrt_m"].quantile(0.10))
                d["K_p90_all_points_MPa_sqrt_m"] = float(rc["KJ_MPa_sqrt_m"].quantile(0.90))
                d["crack_extension_um"] = float(np.nanmax(rc["crack_extension_um"].to_numpy(float)))
            else:
                d["K_mean_all_points_MPa_sqrt_m"] = np.nan
                d["K_median_all_points_MPa_sqrt_m"] = np.nan
                d["K_p10_all_points_MPa_sqrt_m"] = np.nan
                d["K_p90_all_points_MPa_sqrt_m"] = np.nan
            try:
                ext = float(d.get("crack_extension_um", np.nan))
            except Exception:
                ext = np.nan
            d["is_complete"] = bool(np.isfinite(ext) and ext >= 0.98 * float(target_ext_um))
            out_rows.append(d)
            availability.append({
                "rate_factor": float(r),
                "rate_label": rate_label(r),
                "class": d.get("class"),
                "T_K": d.get("T_K"),
                "record_source": d.get("record_source"),
                "case_dir": d.get("case_dir"),
                "case_dir_exists": bool(d.get("case_dir") and Path(str(d.get("case_dir"))).exists()),
                "n_r_curve_points": d.get("n_r_curve_points"),
                "crack_extension_um": d.get("crack_extension_um"),
                "is_complete": d.get("is_complete"),
            })
        print(f"Loaded rate {rate_label(r)}: {len(out_rows)} cases from {rr}")
        rows.extend(out_rows)
    out = pd.DataFrame(rows)
    if not out.empty:
        out.attrs["availability"] = pd.DataFrame(availability)
    return out


def import_v1_module(path: Path):
    spec = importlib.util.spec_from_file_location("v1_tuning_module", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not import V1 module from {path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def prepare_candidate_table(params_csv: Path, classes: list[str]) -> pd.DataFrame:
    df = pd.read_csv(params_csv)
    if "target_class" not in df.columns:
        if "class" in df.columns:
            df = df.rename(columns={"class": "target_class"})
        else:
            raise RuntimeError(f"{params_csv} must contain target_class or class")
    df["target_class"] = df["target_class"].map(normalize_class_name)
    df = df[df["target_class"].isin(classes)].copy()
    df = df.reset_index(drop=True)

    # The V1 tuning driver expects some metadata/default columns that are shared
    # constants in the manuscript table.  Add them here so the exact-parameter
    # CSV can stay compact.
    defaults = {
        "exp_Tref_K": 300.0,
        "surface_id": "final_exact",
        "surface_index": -1,
        "context_id": "final_exact",
        "thermal_shape_stratum": "final",
        "eta_G_Tref_over_G00": np.nan,
        "eta_sigc_Tref_over_sigc0": np.nan,
        "implied_S0_kB": np.nan,
    }
    for k, v in defaults.items():
        if k not in df.columns:
            df[k] = v

    # Older exact-input CSVs store the effective emission values already.
    # The V1 driver expects the values as columns named exp_* and does not
    # require any additional scaling here.
    required = [
        "exp_G00_eV", "exp_gT_eV_per_K", "exp_sigc0_GPa",
        "exp_sT_MPa_per_K", "exp_a", "exp_n", "exp_floor_frac",
        "exp_Tref_K",
        "cleave_G00_eV", "cleave_gT_eV_per_K", "cleave_sigc0_GPa",
        "cleave_sT_MPa_per_K", "cleave_exp_a", "cleave_exp_n",
        "cleave_floor_frac", "cleave_S_hs_kB", "chi_shield", "N_sat",
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise RuntimeError(f"{params_csv} is missing required model columns: {missing}")

    df.insert(0, "candidate_id", np.arange(len(df), dtype=int))

    # Normalize inf strings in N_sat.
    if "N_sat" in df.columns:
        df["N_sat"] = [np.inf if str(x).strip().lower() in {"inf", "infinity", "∞"} else float(x) for x in df["N_sat"]]
    return df


def compute_analytic_with_v1(v1_script: Path, params_csv: Path, classes: list[str],
                             rates: list[float], base_kdot: float, T_dense: np.ndarray,
                             Kmax: float, dK: float) -> pd.DataFrame:
    base = import_v1_module(v1_script)
    cand = prepare_candidate_table(params_csv, classes)
    from arrhenius_fracture.config import ElasticProperties
    mat = ElasticProperties()
    rows = []
    for rf in rates:
        K = base.simulate_candidates(
            cand, T_dense,
            Kmax=Kmax,
            dK=dK,
            Kdot=base_kdot * float(rf),
            G_Pa=mat.G,
            nu=mat.nu,
            b_m=mat.b,
        )
        for i, r in cand.iterrows():
            klass = normalize_class_name(r["target_class"])
            for T, kval in zip(T_dense, K[i]):
                rows.append({
                    "class": klass,
                    "T_K": float(T),
                    "rate_factor": float(rf),
                    "rate_label": rate_label(rf),
                    "K_analytic_MPa_sqrt_m": float(kval) if np.isfinite(kval) else np.nan,
                    "analytic_source": f"computed_V1_Kdot_{base_kdot:g}_times_{rf:g}",
                })
    return pd.DataFrame(rows)


def load_fallback_analytic(path: Path, rates: list[float], classes: list[str]) -> pd.DataFrame:
    if not path.exists():
        raise RuntimeError(f"Fallback analytical CSV does not exist: {path}")
    a = pd.read_csv(path)
    class_col = first_existing(a, ["class", "target_class", "fracture_class", "barrier_class"])
    T_col = first_existing(a, ["T_K", "temperature_K", "temperature", "T"])
    K_col = first_existing(a, ["K_analytic_MPa_sqrt_m", "Kc_MPa_sqrtm", "Kc_MPa_sqrt_m", "target_Kc_MPa_sqrtm"])
    rate_col = first_existing(a, ["rate_factor", "rate", "rate_multiplier"])
    if class_col is None or T_col is None or K_col is None:
        raise RuntimeError(f"Could not identify columns in {path}: {list(a.columns)}")
    out = pd.DataFrame({
        "class": a[class_col].map(normalize_class_name),
        "T_K": pd.to_numeric(a[T_col], errors="coerce"),
        "K_analytic_MPa_sqrt_m": pd.to_numeric(a[K_col], errors="coerce"),
    })
    if rate_col is not None:
        out["rate_factor"] = pd.to_numeric(a[rate_col], errors="coerce")
        out["rate_label"] = out["rate_factor"].map(rate_label)
    else:
        print("WARNING: analytical CSV has no rate_factor column; reusing it for all requested rates.")
        pieces = []
        for rf in rates:
            q = out.copy()
            q["rate_factor"] = float(rf)
            q["rate_label"] = rate_label(rf)
            pieces.append(q)
        out = pd.concat(pieces, ignore_index=True)
    out = out[out["class"].isin(classes)].dropna(subset=["T_K", "K_analytic_MPa_sqrt_m"])
    return out.reset_index(drop=True)


def get_analytic(args, classes: list[str], rates: list[float]) -> pd.DataFrame:
    T_dense = np.arange(args.T_min, args.T_max + 0.1, args.analytic_T_step)
    v1 = Path(args.v1_script)
    params = Path(args.params_csv)
    if args.compute_analytic and v1.exists() and params.exists():
        try:
            print(f"Computing analytical curves with {v1}")
            return compute_analytic_with_v1(
                v1, params, classes, rates,
                base_kdot=args.base_kdot,
                T_dense=T_dense,
                Kmax=args.analytic_Kmax,
                dK=args.analytic_dK,
            )
        except Exception as e:
            print(f"WARNING: V1 analytical computation failed: {e}")
            print("Falling back to analytical CSV.")
    return load_fallback_analytic(Path(args.analytic_csv), rates, classes)


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
        "savefig.bbox": "tight",
        "figure.dpi": 150,
    })
    return plt


RATE_STYLES = {
    1.0: {"color": "#0072B2", "marker": "o"},
    10.0: {"color": "#D55E00", "marker": "s"},
    100.0: {"color": "#009E73", "marker": "^"},
}


def style_for_rate(rate: float):
    if rate in RATE_STYLES:
        return RATE_STYLES[rate]
    colors = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9"]
    markers = ["o", "s", "^", "D", "v", "P"]
    idx = int(abs(hash(rate))) % len(colors)
    return {"color": colors[idx], "marker": markers[idx]}


def make_rate_panel_plot(data: pd.DataFrame, analytic: pd.DataFrame, out: Path,
                         classes: list[str], rates: list[float], metric_col: str,
                         ylabel: str, fig_name: str, title_prefix: str) -> None:
    plt = publication_style()
    from matplotlib.lines import Line2D

    fig, axes = plt.subplots(2, 2, figsize=(12.8, 9.0), sharex=True)
    axes = axes.ravel()

    for ax, klass in zip(axes, classes):
        for rf in rates:
            st = style_for_rate(float(rf))
            ga = analytic[(analytic["class"] == klass) & (np.isclose(analytic["rate_factor"], float(rf)))].sort_values("T_K")
            if not ga.empty:
                ax.plot(
                    ga["T_K"], ga["K_analytic_MPa_sqrt_m"],
                    color=st["color"], linestyle="-", linewidth=2.7,
                    alpha=0.90, zorder=1,
                )
            gd = data[(data["class"] == klass) & (np.isclose(data["rate_factor"], float(rf)))].sort_values("T_K")
            if not gd.empty:
                vals = pd.to_numeric(gd[metric_col], errors="coerce")
                keep = np.isfinite(vals.to_numpy(float))
                gd = gd.loc[keep].copy()
                vals = vals.loc[keep]
                if not gd.empty:
                    ax.plot(gd["T_K"], vals, color=st["color"], linestyle=(0, (4, 2)), linewidth=1.8, zorder=2)
                    comp = gd["is_complete"].fillna(False).to_numpy(bool)
                    if comp.any():
                        ax.scatter(
                            gd.loc[comp, "T_K"], vals.loc[comp],
                            s=62, marker=st["marker"], facecolor=st["color"],
                            edgecolor="white", linewidth=0.8, zorder=3,
                        )
                    if (~comp).any():
                        ax.scatter(
                            gd.loc[~comp, "T_K"], vals.loc[~comp],
                            s=62, marker=st["marker"], facecolor="white",
                            edgecolor=st["color"], linewidth=1.8, zorder=3,
                        )
        ax.set_title(klass)
        ax.set_ylabel(ylabel)
        ax.grid(False)
    for ax in axes[-2:]:
        ax.set_xlabel("Temperature (K)")

    legend = []
    for rf in rates:
        st = style_for_rate(float(rf))
        legend.append(Line2D([0], [0], color=st["color"], linestyle="-", linewidth=2.7, label=f"V1 {rate_label(rf)}"))
        legend.append(Line2D([0], [0], color=st["color"], linestyle=(0, (4, 2)), marker=st["marker"],
                             markerfacecolor=st["color"], markeredgecolor="white", linewidth=1.8,
                             label=f"CZM {rate_label(rf)}"))
    legend.append(Line2D([0], [0], color="#666666", linestyle="None", marker="o",
                         markerfacecolor="white", markeredgecolor="#666666",
                         label="open = incomplete"))
    fig.legend(handles=legend, loc="upper center", bbox_to_anchor=(0.5, 1.04), ncol=4, frameon=False)
    fig.savefig(out / f"{fig_name}.png", dpi=320)
    fig.savefig(out / f"{fig_name}.svg")
    plt.close(fig)

    for klass in classes:
        fig, ax = plt.subplots(figsize=(7.4, 5.4))
        for rf in rates:
            st = style_for_rate(float(rf))
            ga = analytic[(analytic["class"] == klass) & (np.isclose(analytic["rate_factor"], float(rf)))].sort_values("T_K")
            if not ga.empty:
                ax.plot(ga["T_K"], ga["K_analytic_MPa_sqrt_m"], color=st["color"], linestyle="-", linewidth=2.7, label=f"V1 {rate_label(rf)}")
            gd = data[(data["class"] == klass) & (np.isclose(data["rate_factor"], float(rf)))].sort_values("T_K")
            if not gd.empty:
                vals = pd.to_numeric(gd[metric_col], errors="coerce")
                keep = np.isfinite(vals.to_numpy(float))
                gd = gd.loc[keep].copy()
                vals = vals.loc[keep]
                if not gd.empty:
                    ax.plot(gd["T_K"], vals, color=st["color"], linestyle=(0, (4, 2)), linewidth=1.8)
                    comp = gd["is_complete"].fillna(False).to_numpy(bool)
                    if comp.any():
                        ax.scatter(gd.loc[comp, "T_K"], vals.loc[comp], s=62, marker=st["marker"],
                                   facecolor=st["color"], edgecolor="white", linewidth=0.8, zorder=3)
                    if (~comp).any():
                        ax.scatter(gd.loc[~comp, "T_K"], vals.loc[~comp], s=62, marker=st["marker"],
                                   facecolor="white", edgecolor=st["color"], linewidth=1.8, zorder=3)
        ax.set_title(f"{title_prefix}: {klass}")
        ax.set_xlabel("Temperature (K)")
        ax.set_ylabel(ylabel)
        ax.grid(False)
        ax.legend(frameon=False, ncol=2)
        safe = klass.replace("/", "_")
        fig.savefig(out / f"{fig_name}_{safe}.png", dpi=320)
        fig.savefig(out / f"{fig_name}_{safe}.svg")
        plt.close(fig)


def make_error_tables(data: pd.DataFrame, analytic: pd.DataFrame, metric_col: str,
                      value_name: str, out: Path) -> pd.DataFrame:
    rows = []
    for _, r in data.iterrows():
        klass = r["class"]
        rf = float(r["rate_factor"])
        ga = analytic[(analytic["class"] == klass) & (np.isclose(analytic["rate_factor"], rf))].sort_values("T_K")
        if ga.empty:
            ka = np.nan
        else:
            ka = float(np.interp(float(r["T_K"]), ga["T_K"].to_numpy(float), ga["K_analytic_MPa_sqrt_m"].to_numpy(float), left=np.nan, right=np.nan))
        kval = pd.to_numeric(pd.Series([r.get(metric_col, np.nan)]), errors="coerce").iloc[0]
        err = float(kval) - ka if np.isfinite(kval) and np.isfinite(ka) else np.nan
        rel = err / ka if np.isfinite(err) and abs(ka) > 1e-12 else np.nan
        d = r.to_dict()
        d.update({
            f"{value_name}_MPa_sqrt_m": kval,
            "K_analytic_interp_MPa_sqrt_m": ka,
            "error_vs_analytic_MPa_sqrt_m": err,
            "relative_error_vs_analytic": rel,
            "percent_error_vs_analytic": 100.0 * rel if np.isfinite(rel) else np.nan,
        })
        rows.append(d)
    outdf = pd.DataFrame(rows).sort_values(["class", "rate_factor", "T_K"])
    outdf.to_csv(out / f"{value_name}_comparison_with_analytic.csv", index=False)
    metrics = (
        outdf.groupby(["class", "rate_factor", "rate_label"], dropna=False)
        .agg(
            n=("T_K", "count"),
            RMSE_MPa_sqrt_m=("error_vs_analytic_MPa_sqrt_m", lambda x: float(np.sqrt(np.nanmean(np.asarray(x, float)**2))) if len(x) else np.nan),
            mean_abs_percent_error=("percent_error_vs_analytic", lambda x: float(np.nanmean(np.abs(np.asarray(x, float)))) if len(x) else np.nan),
            bias_MPa_sqrt_m=("error_vs_analytic_MPa_sqrt_m", lambda x: float(np.nanmean(np.asarray(x, float))) if len(x) else np.nan),
            complete_fraction=("is_complete", lambda x: float(np.mean(np.asarray(x, bool))) if len(x) else np.nan),
        )
        .reset_index()
    )
    metrics.to_csv(out / f"{value_name}_error_metrics_by_class_rate.csv", index=False)
    return outdf


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", default="runs/four_class_exp_floor_CZM_rates_no_branch_500um_theta45")
    ap.add_argument("--out", default="runs/CZM_three_rate_temperature_comparison")
    ap.add_argument("--rates", default="1 10 100")
    ap.add_argument("--classes", default=" ".join(CLASSES_DEFAULT))
    ap.add_argument("--temps", default=" ".join(map(str, TEMPS_DEFAULT)))
    ap.add_argument("--target-ext-um", type=float, default=500.0)

    ap.add_argument("--compute-analytic", action="store_true", default=True)
    ap.add_argument("--no-compute-analytic", dest="compute_analytic", action="store_false")
    ap.add_argument("--v1-script", default="run_v1_exp_floor_four_class_tuning.py")
    ap.add_argument("--params-csv", default="four_class_exp_floor_exact_model_inputs.csv")
    ap.add_argument("--analytic-csv", default="four_class_analytical_prediction_final.csv")
    ap.add_argument("--base-kdot", type=float, default=0.005)
    ap.add_argument("--analytic-dK", type=float, default=0.02)
    ap.add_argument("--analytic-Kmax", type=float, default=100.0)
    ap.add_argument("--analytic-T-step", type=float, default=5.0)
    ap.add_argument("--T-min", type=float, default=300.0)
    ap.add_argument("--T-max", type=float, default=1200.0)
    args = ap.parse_args()

    root = Path(args.root)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    rates = parse_list(args.rates, float)
    classes = parse_list(args.classes, str)
    temps = parse_list(args.temps, int)

    data = load_rate_data(root, rates, classes, args.target_ext_um)
    data = data[data["class"].isin(classes) & data["T_K"].isin(temps)].copy()
    if data.empty:
        raise SystemExit(f"No CZM rate data found under {root}")
    data.to_csv(out / "CZM_three_rate_toughness_summary.csv", index=False)
    avail = data.attrs.get("availability")
    if isinstance(avail, pd.DataFrame) and not avail.empty:
        avail.to_csv(out / "rate_case_availability.csv", index=False)
    print("Loaded case counts by rate:")
    print(data.groupby(["rate_label", "class"]).size().unstack(fill_value=0).to_string())

    analytic = get_analytic(args, classes, rates)
    analytic = analytic[analytic["class"].isin(classes)].copy()
    analytic.to_csv(out / "analytical_predictions_by_rate.csv", index=False)
    # Check whether analytical curves actually differ by rate.
    try:
        piv = analytic.pivot_table(index=["class", "T_K"], columns="rate_factor", values="K_analytic_MPa_sqrt_m", aggfunc="mean")
        if piv.shape[1] > 1 and float(piv.max(axis=1).sub(piv.min(axis=1)).max()) < 1e-8:
            print("WARNING: analytical curves are effectively identical across rates. "
                  "This usually means the fallback CSV was used rather than rate-specific V1 computation.")
    except Exception:
        pass

    make_rate_panel_plot(
        data, analytic, out, classes, rates,
        metric_col="Kc_first_MPa_sqrt_m",
        ylabel=r"$K_c$ (MPa$\sqrt{m}$)",
        fig_name="CZM_three_rate_first_passage_vs_analytic",
        title_prefix="First passage",
    )
    make_rate_panel_plot(
        data, analytic, out, classes, rates,
        metric_col="K_median_all_points_MPa_sqrt_m",
        ylabel=r"Median $K_J$ (MPa$\sqrt{m}$)",
        fig_name="CZM_three_rate_median_Rcurve_vs_analytic",
        title_prefix="Median propagation",
    )

    make_error_tables(data, analytic, "Kc_first_MPa_sqrt_m", "first_passage", out)
    make_error_tables(data, analytic, "K_median_all_points_MPa_sqrt_m", "median_Rcurve", out)

    cfg = vars(args).copy()
    cfg["rates"] = rates
    cfg["classes"] = classes
    cfg["temps"] = temps
    cfg["outputs"] = [
        "CZM_three_rate_first_passage_vs_analytic.png",
        "CZM_three_rate_median_Rcurve_vs_analytic.png",
        "CZM_three_rate_toughness_summary.csv",
        "analytical_predictions_by_rate.csv",
        "first_passage_comparison_with_analytic.csv",
        "median_Rcurve_comparison_with_analytic.csv",
    ]
    (out / "comparison_config.json").write_text(json.dumps(cfg, indent=2))

    print(f"WROTE {out}")
    print(f"  {out / 'CZM_three_rate_first_passage_vs_analytic.png'}")
    print(f"  {out / 'CZM_three_rate_median_Rcurve_vs_analytic.png'}")


if __name__ == "__main__":
    main()
