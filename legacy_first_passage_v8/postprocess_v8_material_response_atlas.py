#!/usr/bin/env python3
"""Postprocess 2-D v8 material-response atlas runs.

The script reads case subfolders produced by the atlas/long-growth/production
shell drivers, extracts v8_2D rows from compare_summary.csv, computes Paris-style
crack-growth points, and optionally extracts local block-by-block crack-growth
increments from steps_0300K.csv.

Important plotting convention:
  * No-growth points are retained in CSV as censored upper bounds but are not
    connected to the measured Paris curve by default.
  * For long-growth runs, the default x-axis is DeltaK_initial. DeltaK_final can
    be misleading because K increases as the crack extends.
  * Local block points give da_block/dN_block versus the local KJ of that block.
"""
from __future__ import annotations

import argparse
import math
import re
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import pandas as pd


def _safe_float(x, default=float("nan")):
    try:
        if pd.isna(x):
            return default
        return float(x)
    except Exception:
        return default


def _first_finite(*vals, default=float("nan")):
    for v in vals:
        f = _safe_float(v)
        if math.isfinite(f):
            return f
    return default


def _safe_int(x, default=0):
    try:
        if pd.isna(x):
            return default
        return int(float(x))
    except Exception:
        return default


def load_case_table(path: Path | None) -> pd.DataFrame:
    if path is None or not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path)


def _resolve_run_dir(run_dir_value, root: Path, case_dir: Path, k_label: str | None = None) -> Path | None:
    candidates = []
    if isinstance(run_dir_value, str) and run_dir_value.strip():
        p = Path(run_dir_value)
        candidates.append(p)
        if not p.is_absolute():
            candidates.append(root / p)
            candidates.append(case_dir / p)
    if k_label:
        candidates.append(case_dir / f"v8_2d_K{k_label}")
    for c in candidates:
        if c.exists():
            return c
    return None


def _k_label_from_row(r: pd.Series) -> str:
    k = _first_finite(r.get("target_Kmax_MPa_sqrtm"), r.get("actual_Kmax_MPa_sqrtm"), r.get("K_initial_MPa_sqrtm"))
    if not math.isfinite(k):
        return ""
    # The driver uses labels like K6, K6p5, etc.
    s = (f"{k:g}").replace(".", "p")
    return s


def compute_points(root: Path, case_table: pd.DataFrame, R: float) -> pd.DataFrame:
    rows = []
    if not case_table.empty and "case_label" in case_table.columns:
        labels = [str(x) for x in case_table["case_label"].tolist()]
    else:
        labels = [p.name for p in root.iterdir() if p.is_dir()]

    meta_by_label = {}
    if not case_table.empty and "case_label" in case_table.columns:
        for _, r in case_table.iterrows():
            meta_by_label[str(r["case_label"])] = r.to_dict()

    for label in labels:
        cdir = root / label
        csv_path = cdir / "compare_summary.csv"
        if not csv_path.exists():
            print(f"warning: missing {csv_path}")
            continue
        df = pd.read_csv(csv_path)
        if "model" not in df.columns:
            print(f"warning: no model column in {csv_path}")
            continue
        v8 = df[df["model"] == "v8_2D"].copy()
        meta = meta_by_label.get(label, {})
        for _, r in v8.iterrows():
            targetK = _safe_float(r.get("target_Kmax_MPa_sqrtm"))
            K_initial = _first_finite(r.get("K_initial_MPa_sqrtm"), targetK)
            K_first = _safe_float(r.get("K_first_fire_MPa_sqrtm"))
            K_final = _safe_float(r.get("K_final_MPa_sqrtm"))
            DeltaK_nominal = (1.0 - R) * targetK if math.isfinite(targetK) else float("nan")
            DeltaK_initial = (1.0 - R) * K_initial if math.isfinite(K_initial) else float("nan")
            DeltaK_first = (1.0 - R) * K_first if math.isfinite(K_first) else float("nan")
            DeltaK_final = (1.0 - R) * K_final if math.isfinite(K_final) else float("nan")
            DeltaK_summary = _first_finite(r.get("DeltaK_MPa_sqrtm"), DeltaK_initial, DeltaK_nominal)

            cycles = _safe_float(r.get("cycles_total"), 0.0)
            cycles_first = _safe_float(r.get("cycles_to_first_fire"))
            a_adv_m = _safe_float(r.get("a_adv_um"), 0.0) * 1e-6
            n_adv = _safe_int(r.get("n_adv_or_fire_total"), 0)
            da_dN = _safe_float(r.get("da_dN_m_per_cycle"))
            if not math.isfinite(da_dN):
                da_dN = a_adv_m / max(cycles, 1e-300)
            ub = _safe_float(r.get("da_dN_upper_bound_m_per_cycle"))
            if not math.isfinite(ub):
                # Default one-event resolution bound for a v8 5 micron advance quantum.
                ub = 5.0e-6 / max(cycles, 1e-300) if n_adv == 0 else float("nan")
            measured = (n_adv > 0 and da_dN > 0 and math.isfinite(da_dN))
            plot_da = da_dN if measured else ub

            k_label = _k_label_from_row(r)
            run_dir = _resolve_run_dir(r.get("run_dir", ""), root, cdir, k_label)
            run_dir_str = str(run_dir) if run_dir is not None else str(r.get("run_dir", cdir))

            rows.append({
                "case_label": label,
                "source_case": meta.get("source_case", ""),
                "material_response_class": meta.get("material_response_class", ""),
                "target_Kmax_MPa_sqrtm": targetK,
                "K_initial_MPa_sqrtm": K_initial,
                "K_first_fire_MPa_sqrtm": K_first,
                "K_final_MPa_sqrtm": K_final,
                "DeltaK_nominal_MPa_sqrtm": DeltaK_nominal,
                "DeltaK_initial_MPa_sqrtm": DeltaK_initial,
                "DeltaK_first_fire_MPa_sqrtm": DeltaK_first,
                "DeltaK_final_MPa_sqrtm": DeltaK_final,
                "DeltaK_MPa_sqrtm": DeltaK_summary,
                "cycles_total": cycles,
                "cycles_to_first_fire": cycles_first,
                "n_adv_or_fire_total": n_adv,
                "a_adv_m": a_adv_m,
                "da_dN_m_per_cycle": da_dN if measured else float("nan"),
                "da_dN_upper_bound_m_per_cycle": ub if not measured else float("nan"),
                "plot_da_dN_m_per_cycle": plot_da,
                "is_censored_upper_bound": (not measured),
                "direct_lt_1_cycle": (math.isfinite(cycles_first) and cycles_first < 1.0),
                "log10_da_dN_or_bound": math.log10(plot_da) if plot_da > 0 and math.isfinite(plot_da) else float("nan"),
                "B_final": _safe_float(r.get("B_final")),
                "N_em_final": _safe_float(r.get("N_em_final")),
                "dN_store_last_block": _safe_float(r.get("dN_store_last_block")),
                "dN_mobile_last_block": _safe_float(r.get("dN_mobile_last_block")),
                "dN_escape_last_block": _safe_float(r.get("dN_escape_last_block")),
                "cyclic_plastic_work_J": _safe_float(r.get("cyclic_plastic_work_J")),
                "G_cleave_eff_eV": _safe_float(r.get("G_cleave_eff_eV")),
                "S_cleave_kB": _safe_float(r.get("S_cleave_kB")),
                "K_calibration_rel_error": _safe_float(r.get("K_calibration_rel_error")),
                "cycle_limiter_code": _safe_float(r.get("cycle_limiter_code")),
                "blocks_completed": _safe_float(r.get("blocks_completed")),
                "run_dir": run_dir_str,
            })
    return pd.DataFrame(rows)


def extract_local_points(points: pd.DataFrame, R: float) -> pd.DataFrame:
    """Extract local block-level Paris points from steps_0300K.csv.

    For each accepted block with da_block_m > 0, compute da_block_m/fatigue_cycles
    and pair it with that block's KJ. This is better for long-growth runs where K
    evolves strongly during propagation.
    """
    rows = []
    if points.empty:
        return pd.DataFrame()
    for _, pr in points.iterrows():
        run_dir = Path(str(pr.get("run_dir", "")))
        if not run_dir.exists():
            continue
        step_files = sorted(run_dir.glob("steps_*K.csv"))
        if not step_files:
            continue
        # Use first temperature file; current atlas is one T per run.
        sf = step_files[0]
        try:
            df = pd.read_csv(sf)
        except Exception as e:
            print(f"warning: failed reading {sf}: {e}")
            continue
        needed = {"da_block_m", "fatigue_cycles", "KJ_Pa_sqrtm"}
        if not needed.issubset(df.columns):
            continue
        blocks = df[(df["da_block_m"].astype(float) > 0) & (df["fatigue_cycles"].astype(float) > 0)].copy()
        if blocks.empty:
            continue
        for _, r in blocks.iterrows():
            da = _safe_float(r.get("da_block_m"))
            dN = _safe_float(r.get("fatigue_cycles"))
            KJ = _safe_float(r.get("KJ_Pa_sqrtm")) / 1e6
            if not (da > 0 and dN > 0 and KJ > 0):
                continue
            rows.append({
                "case_label": pr["case_label"],
                "source_case": pr.get("source_case", ""),
                "material_response_class": pr.get("material_response_class", ""),
                "target_Kmax_MPa_sqrtm": pr.get("target_Kmax_MPa_sqrtm", float("nan")),
                "step": _safe_float(r.get("step")),
                "DeltaK_local_MPa_sqrtm": (1.0 - R) * KJ,
                "KJ_local_MPa_sqrtm": KJ,
                "dN_block": dN,
                "da_block_m": da,
                "da_dN_local_m_per_cycle": da / dN,
                "log10_da_dN_local": math.log10(da / dN) if da / dN > 0 else float("nan"),
                "crack_extension_m": _safe_float(r.get("crack_extension_m")),
                "B": _safe_float(r.get("B")),
                "N_em": _safe_float(r.get("N_em")),
                "pz_store_total": _safe_float(r.get("pz_store_total")),
                "pz_mobile_total": _safe_float(r.get("pz_mobile_total")),
                "cyclic_plastic_work_J": _safe_float(r.get("cyclic_plastic_work_J")),
                "G_cleave_eff_eV": _safe_float(r.get("G_cleave_eff_eV")),
                "S_cleave_kB": _safe_float(r.get("S_cleave_kB")),
                "run_dir": str(run_dir),
            })
    return pd.DataFrame(rows)


def summarize(points: pd.DataFrame, cycles_max: float = float("nan"), target_ext_um: float = float("nan")) -> pd.DataFrame:
    rows = []
    target_ext_m = target_ext_um * 1e-6 if math.isfinite(target_ext_um) else float("nan")
    for label, g in points.groupby("case_label", sort=False):
        measured = g[~g["is_censored_upper_bound"]].copy()
        censored = g[g["is_censored_upper_bound"]].copy()
        slope = float("nan")
        max_jump = float("nan")
        if len(measured) >= 2:
            measured = measured.sort_values("DeltaK_initial_MPa_sqrtm")
            x = measured["DeltaK_initial_MPa_sqrtm"].astype(float)
            y = measured["log10_da_dN_or_bound"].astype(float)
            dx = x.max() - x.min()
            slope = (y.iloc[-1] - y.iloc[0]) / dx if dx > 0 else float("nan")
            jumps = y.diff().abs().dropna()
            max_jump = jumps.max() if not jumps.empty else float("nan")
        cens_reached_horizon = float("nan")
        if math.isfinite(cycles_max) and cycles_max > 0 and not censored.empty:
            cens_reached_horizon = float((censored["cycles_total"] >= 0.999 * cycles_max).mean())
        reached_target = float("nan")
        if math.isfinite(target_ext_m) and not g.empty:
            reached_target = float((g["a_adv_m"] >= 0.999 * target_ext_m).mean())
        rows.append({
            "case_label": label,
            "source_case": g["source_case"].iloc[0] if "source_case" in g.columns else "",
            "material_response_class": g["material_response_class"].iloc[0] if "material_response_class" in g.columns else "",
            "n_points": len(g),
            "n_measured": int((~g["is_censored_upper_bound"]).sum()),
            "n_censored": int(g["is_censored_upper_bound"].sum()),
            "n_direct_lt_1_cycle": int(g["direct_lt_1_cycle"].sum()),
            "n_reached_target_extension": int((g["a_adv_m"] >= 0.999 * target_ext_m).sum()) if math.isfinite(target_ext_m) else "",
            "fraction_points_reaching_target_extension": reached_target,
            "min_measured_DeltaK_initial": measured["DeltaK_initial_MPa_sqrtm"].min() if not measured.empty else float("nan"),
            "max_measured_DeltaK_initial": measured["DeltaK_initial_MPa_sqrtm"].max() if not measured.empty else float("nan"),
            "min_log10_da_dN_or_bound": g["log10_da_dN_or_bound"].min(),
            "max_log10_da_dN_or_bound": g["log10_da_dN_or_bound"].max(),
            "slope_log10_da_dN_per_MPa_sqrtm_initial": slope,
            "max_adjacent_log10_jump_measured": max_jump,
            "min_censored_cycles_total": censored["cycles_total"].min() if not censored.empty else float("nan"),
            "max_censored_B_final": censored["B_final"].max() if not censored.empty else float("nan"),
            "fraction_censored_points_reaching_cycles_max": cens_reached_horizon,
            "max_abs_K_calibration_rel_error": g["K_calibration_rel_error"].abs().max(),
            "max_N_em_final": g["N_em_final"].max(),
            "max_cyclic_plastic_work_J": g["cyclic_plastic_work_J"].max(),
        })
    return pd.DataFrame(rows)


def _xcol_from_choice(choice: str) -> str:
    return {
        "nominal": "DeltaK_nominal_MPa_sqrtm",
        "initial": "DeltaK_initial_MPa_sqrtm",
        "first_fire": "DeltaK_first_fire_MPa_sqrtm",
        "final": "DeltaK_final_MPa_sqrtm",
        "summary": "DeltaK_MPa_sqrtm",
    }[choice]


def make_summary_plot(points: pd.DataFrame, out_png: Path, *, show_bounds: bool = False, loglog: bool = True, x_choice: str = "initial") -> None:
    fig, ax = plt.subplots(figsize=(8.4, 5.8))
    xcol = _xcol_from_choice(x_choice)
    for label, g in points.groupby("case_label", sort=False):
        g = g.sort_values(xcol)
        measured = g[~g["is_censored_upper_bound"]]
        censored = g[g["is_censored_upper_bound"]]
        if not measured.empty:
            ax.plot(measured[xcol], measured["da_dN_m_per_cycle"], marker="o", label=label)
        if show_bounds and not censored.empty:
            ax.scatter(censored[xcol], censored["da_dN_upper_bound_m_per_cycle"], marker="v", s=34, label=f"{label} upper bound")
    ax.set_yscale("log")
    if loglog:
        ax.set_xscale("log")
    xlab = {
        "nominal": r"Nominal $\Delta K=(1-R)K_{target}$ (MPa $\sqrt{m}$)",
        "initial": r"Initial $\Delta K=(1-R)K_J^0$ (MPa $\sqrt{m}$)",
        "first_fire": r"First-fire $\Delta K$ (MPa $\sqrt{m}$)",
        "final": r"Final $\Delta K$ (MPa $\sqrt{m}$)",
        "summary": r"$\Delta K$ (MPa $\sqrt{m}$)",
    }[x_choice]
    ax.set_xlabel(xlab)
    ax.set_ylabel(r"integrated $da/dN$ (m/cycle)")
    ax.set_title("2-D v8 material-response atlas")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_png, dpi=240)
    plt.close(fig)


def make_local_plot(local: pd.DataFrame, out_png: Path, *, loglog: bool = True) -> None:
    if local.empty:
        return
    fig, ax = plt.subplots(figsize=(8.4, 5.8))
    for label, g in local.groupby("case_label", sort=False):
        g = g.sort_values("DeltaK_local_MPa_sqrtm")
        ax.scatter(g["DeltaK_local_MPa_sqrtm"], g["da_dN_local_m_per_cycle"], s=18, alpha=0.72, label=label)
    ax.set_yscale("log")
    if loglog:
        ax.set_xscale("log")
    ax.set_xlabel(r"Local block $\Delta K=(1-R)K_J$ (MPa $\sqrt{m}$)")
    ax.set_ylabel(r"Local block $da/dN=\Delta a/\Delta N$ (m/cycle)")
    ax.set_title("2-D v8 local block Paris points")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_png, dpi=240)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", default="runs/v8_material_response_atlas_2d")
    ap.add_argument("--case-table", default="selected_2d_material_response_cases.csv")
    ap.add_argument("--R", type=float, default=0.1)
    ap.add_argument("--cycles-max", type=float, default=float("nan"), help="Physical cycle horizon used in the runs; used only to diagnose censored points.")
    ap.add_argument("--target-crack-extension-um", type=float, default=float("nan"), help="Target extension used by the driver; used only to label completed/incomplete points.")
    ap.add_argument("--show-censored-upper-bounds", action="store_true", help="Plot no-growth points as upper-bound markers. By default they remain in the CSV but are omitted from the plot.")
    ap.add_argument("--plot-scale", choices=["loglog", "semilogy"], default="loglog", help="Final Paris plot axis scaling.")
    ap.add_argument("--x-deltaK", choices=["nominal", "initial", "first_fire", "final", "summary"], default="initial", help="Which DeltaK value to use for the integrated Paris plot.")
    ap.add_argument("--extract-local-points", action="store_true", help="Read step histories and write/plot local block da/dN versus local KJ.")
    args = ap.parse_args()

    root = Path(args.root)
    root.mkdir(parents=True, exist_ok=True)
    case_table = load_case_table(Path(args.case_table))
    points = compute_points(root, case_table, args.R)
    points.to_csv(root / "atlas_2d_paris_points.csv", index=False)
    summary = summarize(points, cycles_max=args.cycles_max, target_ext_um=args.target_crack_extension_um) if not points.empty else pd.DataFrame()
    summary.to_csv(root / "atlas_2d_case_summary.csv", index=False)
    if not points.empty:
        make_summary_plot(points, root / "atlas_2d_da_dN_vs_DeltaK.png", show_bounds=args.show_censored_upper_bounds, loglog=(args.plot_scale == "loglog"), x_choice=args.x_deltaK)
        if not args.show_censored_upper_bounds:
            make_summary_plot(points, root / "atlas_2d_da_dN_vs_DeltaK_with_upper_bounds.png", show_bounds=True, loglog=(args.plot_scale == "loglog"), x_choice=args.x_deltaK)
    if args.extract_local_points and not points.empty:
        local = extract_local_points(points, args.R)
        local.to_csv(root / "atlas_2d_local_paris_points.csv", index=False)
        make_local_plot(local, root / "atlas_2d_local_da_dN_vs_local_DeltaK.png", loglog=(args.plot_scale == "loglog"))
        print(f"wrote {root / 'atlas_2d_local_paris_points.csv'}")
        if not local.empty:
            print(f"wrote {root / 'atlas_2d_local_da_dN_vs_local_DeltaK.png'}")
    print(f"wrote {root / 'atlas_2d_paris_points.csv'}")
    print(f"wrote {root / 'atlas_2d_case_summary.csv'}")
    if not points.empty:
        print(f"wrote {root / 'atlas_2d_da_dN_vs_DeltaK.png'}")
        if not args.show_censored_upper_bounds:
            print(f"wrote {root / 'atlas_2d_da_dN_vs_DeltaK_with_upper_bounds.png'}")


if __name__ == "__main__":
    main()
