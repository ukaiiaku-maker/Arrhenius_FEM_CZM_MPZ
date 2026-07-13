#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

plt.rcParams.update({
    "font.size": 13,
    "axes.labelsize": 15,
    "axes.titlesize": 16,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 10,
    "lines.linewidth": 2.0,
    "savefig.bbox": "tight",
})

CASE_ORDER = [
    "FCC_like_case29",
    "shifted_ductile_case64",
    "steep_cleavage_case35",
    "slow_threshold_case101",
    "higher_barrier_case171",
    "plastic_shielded_case64_M1",
]


def load_case_table(root: Path) -> pd.DataFrame:
    p = root / "fem_czm_fatigue_cases.csv"
    if p.exists():
        return pd.read_csv(p)
    return pd.DataFrame({"case_label": CASE_ORDER})


def find_summary_files(root: Path) -> list[Path]:
    return sorted(root.glob("*/compare_summary.csv")) + sorted(root.glob("*/case_summary.csv"))


def _to_num(s, default=np.nan):
    try:
        return pd.to_numeric(s, errors="coerce")
    except Exception:
        return default


def read_integrated_points(root: Path, R: float, cycles_max: float) -> pd.DataFrame:
    """Read one integrated da/dN point per K condition.

    Measured points use total accepted crack extension / total cycles.
    Runs with no accepted advance are plotted as one-advance upper bounds,
    da_phys / total cycles, matching the prior atlas convention.
    """
    rows = []
    for p in find_summary_files(root):
        try:
            df = pd.read_csv(p)
        except Exception:
            continue
        if df.empty:
            continue
        case = p.parent.name
        if "model" in df.columns:
            mask = df["model"].astype(str).str.contains("2D|v8", case=False, regex=True)
            if mask.any():
                df = df.loc[mask].copy()
        for _, r in df.iterrows():
            K0 = float(r.get("K_initial_MPa_sqrtm", r.get("actual_Kmax_MPa_sqrtm", r.get("target_Kmax_MPa_sqrtm", np.nan))))
            Ktarget = float(r.get("target_Kmax_MPa_sqrtm", K0))
            # The atlas x-axis is the initial cyclic range, not the geometry-evolved final K.
            deltaK = (1.0 - R) * K0
            cycles = float(r.get("cycles_total", np.nan))
            blocks_completed = float(r.get("blocks_completed", np.nan))
            nfire = float(r.get("n_adv_or_fire_total", 0.0))
            a_um = float(r.get("a_adv_um", np.nan))
            a_m = a_um * 1e-6 if np.isfinite(a_um) else float(r.get("a_adv_m", np.nan))
            measured_rate = float(r.get("da_dN_m_per_cycle", np.nan))
            upper = float(r.get("da_dN_upper_bound_m_per_cycle", np.nan))
            measured = bool(nfire > 0 and np.isfinite(a_m) and a_m > 0 and np.isfinite(cycles) and cycles > 0)
            if measured:
                rate = a_m / cycles if not (np.isfinite(measured_rate) and measured_rate > 0) else measured_rate
            else:
                if not (np.isfinite(upper) and upper > 0):
                    # Recover da_phys from a representative run_args.json when possible.
                    da_phys = np.nan
                    klabel = (f"{Ktarget:g}").replace(".", "p")
                    run_args = p.parent / f"v8_2d_K{klabel}" / "run_args.json"
                    if run_args.exists():
                        try:
                            da_phys = float(json.loads(run_args.read_text()).get("da_phys", np.nan))
                        except Exception:
                            pass
                    if not np.isfinite(da_phys) or da_phys <= 0:
                        da_phys = 5e-6
                    upper = da_phys / cycles if np.isfinite(cycles) and cycles > 0 else np.nan
                rate = upper
            rows.append({
                "case_label": case,
                "target_Kmax_MPa_sqrtm": Ktarget,
                "K_initial_MPa_sqrtm": K0,
                "DeltaK_initial_MPa_sqrtm": deltaK,
                "cycles_total": cycles,
                "blocks_completed": blocks_completed,
                "n_adv_or_fire_total": nfire,
                "a_adv_um": a_um,
                "da_dN_m_per_cycle": measured_rate,
                "da_dN_upper_bound_m_per_cycle": upper,
                "plot_da_dN_m_per_cycle": rate,
                "is_censored_upper_bound": not measured,
                "termination_is_cycle_horizon": bool(np.isfinite(cycles) and cycles >= 0.999 * cycles_max),
                "bound_kind": ("measured" if measured else ("cycle_horizon_upper_bound" if (np.isfinite(cycles) and cycles >= 0.999 * cycles_max) else "block_limited_bound")),
                "summary_file": str(p),
            })
    if not rows:
        return pd.DataFrame()
    out = pd.DataFrame(rows)
    out = out[np.isfinite(out["DeltaK_initial_MPa_sqrtm"]) & np.isfinite(out["plot_da_dN_m_per_cycle"]) & (out["plot_da_dN_m_per_cycle"] > 0)].copy()
    order = {c: i for i, c in enumerate(CASE_ORDER)}
    out["_order"] = out["case_label"].map(order).fillna(99)
    return out.sort_values(["_order", "DeltaK_initial_MPa_sqrtm"]).drop(columns="_order").reset_index(drop=True)


def _infer_K_target_from_dir(run_dir: Path) -> float:
    m = re.search(r"v8_2d_K([0-9p]+)$", run_dir.name)
    if not m:
        return np.nan
    try:
        return float(m.group(1).replace("p", "."))
    except Exception:
        return np.nan


def read_local_block_points(root: Path, R: float) -> pd.DataFrame:
    """Extract local block da/dN versus local J-derived DeltaK.

    Only accepted growth blocks are included: da_block_m > 0 and block cycles > 0.
    This is the same physical definition used in the prior PF/V8 atlas,
    (da/dN)_k = Delta a_k / Delta N_k.
    """
    rows = []
    for case_dir in sorted(root.iterdir()):
        if not case_dir.is_dir() or case_dir.name == "analysis" or case_dir.name.startswith("_"):
            continue
        case = case_dir.name
        for run_dir in sorted(case_dir.glob("v8_2d_K*")):
            if not run_dir.is_dir():
                continue
            steps_files = sorted(run_dir.glob("steps_*K.csv"))
            if not steps_files:
                continue
            try:
                st = pd.read_csv(steps_files[0])
            except Exception:
                continue
            required = {"KJ_Pa_sqrtm", "fatigue_cycles", "da_block_m"}
            if not required.issubset(st.columns):
                continue
            KJ = pd.to_numeric(st["KJ_Pa_sqrtm"], errors="coerce").to_numpy(float) / 1e6
            dN = pd.to_numeric(st["fatigue_cycles"], errors="coerce").to_numpy(float)
            da = pd.to_numeric(st["da_block_m"], errors="coerce").fillna(0.0).to_numpy(float)
            nfire = pd.to_numeric(st["n_fire"], errors="coerce").fillna(0.0).to_numpy(float) if "n_fire" in st.columns else np.zeros(len(st))
            local_dK = (1.0 - R) * KJ
            rate = da / np.maximum(dN, 1e-300)
            growth = np.isfinite(local_dK) & np.isfinite(rate) & (local_dK > 0) & (dN > 0) & ((da > 0) | (nfire > 0)) & (rate > 0)
            Ktarget = _infer_K_target_from_dir(run_dir)
            idx = np.flatnonzero(growth)
            for j in idx:
                rows.append({
                    "case_label": case,
                    "target_Kmax_MPa_sqrtm": Ktarget,
                    "step": float(st.iloc[j]["step"]) if "step" in st.columns else int(j),
                    "local_KJ_MPa_sqrtm": KJ[j],
                    "local_DeltaK_MPa_sqrtm": local_dK[j],
                    "da_block_m": da[j],
                    "block_cycles": dN[j],
                    "local_da_dN_m_per_cycle": rate[j],
                    "n_fire": nfire[j],
                    "crack_extension_um": float(st.iloc[j]["crack_extension_m"]) * 1e6 if "crack_extension_m" in st.columns else np.nan,
                    "run_dir": str(run_dir),
                })
    if not rows:
        return pd.DataFrame(columns=[
            "case_label", "target_Kmax_MPa_sqrtm", "step", "local_KJ_MPa_sqrtm",
            "local_DeltaK_MPa_sqrtm", "da_block_m", "block_cycles",
            "local_da_dN_m_per_cycle", "n_fire", "crack_extension_um", "run_dir"
        ])
    out = pd.DataFrame(rows)
    order = {c: i for i, c in enumerate(CASE_ORDER)}
    out["_order"] = out["case_label"].map(order).fillna(99)
    return out.sort_values(["_order", "local_DeltaK_MPa_sqrtm", "target_Kmax_MPa_sqrtm", "step"]).drop(columns="_order").reset_index(drop=True)


def summarize(integrated: pd.DataFrame, local: pd.DataFrame, cases: pd.DataFrame) -> pd.DataFrame:
    rows = []
    meta = cases.set_index("case_label") if "case_label" in cases.columns else pd.DataFrame()
    for case in CASE_ORDER:
        g = integrated[integrated.case_label == case] if not integrated.empty else pd.DataFrame()
        l = local[local.case_label == case] if not local.empty else pd.DataFrame()
        if g.empty and l.empty:
            continue
        measured = g[~g["is_censored_upper_bound"]] if not g.empty else pd.DataFrame()
        cens = g[g["is_censored_upper_bound"]] if not g.empty else pd.DataFrame()
        row = {
            "case_label": case,
            "n_integrated_points": len(g),
            "n_measured": len(measured),
            "n_upper_bounds": len(cens),
            "n_local_growth_blocks": len(l),
            "min_DeltaK_initial": float(g["DeltaK_initial_MPa_sqrtm"].min()) if len(g) else np.nan,
            "max_DeltaK_initial": float(g["DeltaK_initial_MPa_sqrtm"].max()) if len(g) else np.nan,
            "min_integrated_rate_or_bound": float(g["plot_da_dN_m_per_cycle"].min()) if len(g) else np.nan,
            "max_integrated_rate_or_bound": float(g["plot_da_dN_m_per_cycle"].max()) if len(g) else np.nan,
            "min_local_rate": float(l["local_da_dN_m_per_cycle"].min()) if len(l) else np.nan,
            "max_local_rate": float(l["local_da_dN_m_per_cycle"].max()) if len(l) else np.nan,
        }
        if not meta.empty and case in meta.index:
            for c in ["source_case", "material_response_class"]:
                if c in meta.columns:
                    row[c] = meta.loc[case, c]
        rows.append(row)
    return pd.DataFrame(rows)


def plot_integrated(points: pd.DataFrame, out: Path):
    fig, ax = plt.subplots(figsize=(10.2, 6.8))
    cmap = plt.get_cmap("tab10")
    for i, case in enumerate(CASE_ORDER):
        g = points[points.case_label == case].sort_values("DeltaK_initial_MPa_sqrtm")
        if g.empty:
            continue
        color = cmap(i % 10)
        m = ~g["is_censored_upper_bound"].astype(bool)
        c = ~m
        if m.any():
            gm = g.loc[m]
            ax.plot(gm["DeltaK_initial_MPa_sqrtm"], gm["plot_da_dN_m_per_cycle"], marker="o", markersize=7, color=color, label=case)
        elif c.any():
            # ensure the case still appears in the legend
            ax.plot([], [], marker="o", color=color, label=case)
        if c.any():
            gc = g.loc[c]
            gh = gc[gc["termination_is_cycle_horizon"].astype(bool)]
            gb = gc[~gc["termination_is_cycle_horizon"].astype(bool)]
            if not gh.empty:
                ax.scatter(gh["DeltaK_initial_MPa_sqrtm"], gh["plot_da_dN_m_per_cycle"], marker="v", s=75, color=color, label=f"{case} cycle-horizon upper bound")
            if not gb.empty:
                ax.scatter(gb["DeltaK_initial_MPa_sqrtm"], gb["plot_da_dN_m_per_cycle"], marker="v", s=75, facecolors="none", edgecolors=color, linewidths=1.5, label=f"{case} block-limited bound")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"Initial $\Delta K=(1-R)K_J^0$ (MPa $\sqrt{m}$)")
    ax.set_ylabel(r"Integrated $da/dN$ (m/cycle)")
    ax.set_title("2-D FEM/CZM material-response atlas")
    ax.grid(True, which="major", alpha=0.25)
    ax.legend(frameon=True, ncol=1)
    fig.tight_layout()
    for name in ["atlas_2d_da_dN_vs_DeltaK_with_upper_bounds", "fem_czm_da_dN_vs_DeltaK_with_upper_bounds"]:
        fig.savefig(out / f"{name}.png", dpi=320)
        fig.savefig(out / f"{name}.svg")
    plt.close(fig)


def plot_measured(points: pd.DataFrame, out: Path):
    fig, ax = plt.subplots(figsize=(9.2, 6.2))
    cmap = plt.get_cmap("tab10")
    for i, case in enumerate(CASE_ORDER):
        g = points[(points.case_label == case) & (~points["is_censored_upper_bound"].astype(bool))].sort_values("DeltaK_initial_MPa_sqrtm")
        if g.empty:
            continue
        ax.plot(g["DeltaK_initial_MPa_sqrtm"], g["plot_da_dN_m_per_cycle"], marker="o", color=cmap(i % 10), label=case)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"Initial $\Delta K=(1-R)K_J^0$ (MPa $\sqrt{m}$)")
    ax.set_ylabel(r"Integrated $da/dN$ (m/cycle)")
    ax.set_title("2-D FEM/CZM measured fatigue crack-growth rates")
    ax.grid(True, which="major", alpha=0.25)
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(frameon=True)
    fig.tight_layout()
    fig.savefig(out / "fem_czm_da_dN_vs_DeltaK_measured_only.png", dpi=320)
    fig.savefig(out / "fem_czm_da_dN_vs_DeltaK_measured_only.svg")
    plt.close(fig)


def plot_local(local: pd.DataFrame, out: Path):
    """Plot local Paris points with color representing accumulated crack extension.

    Material-response class is retained through marker shape.  A single global
    color normalization is used so that the same color corresponds to the same
    crack extension in every class.  This makes it possible to distinguish
    early-tip points from later propagation points without conflating crack
    extension with the material-class legend.
    """
    fig, ax = plt.subplots(figsize=(10.4, 7.0))

    markers = ["o", "s", "^", "D", "P", "X"]
    cmap = plt.get_cmap("viridis")

    if local.empty:
        ax.text(
            0.5, 0.5, "No accepted crack-growth blocks in this run",
            transform=ax.transAxes, ha="center", va="center"
        )
    else:
        ext_all = pd.to_numeric(local["crack_extension_um"], errors="coerce")
        finite_ext = ext_all[np.isfinite(ext_all)]
        if finite_ext.empty:
            vmin, vmax = 0.0, 1.0
        else:
            vmin = min(0.0, float(finite_ext.min()))
            vmax = float(finite_ext.max())
            if not np.isfinite(vmax) or vmax <= vmin:
                vmax = vmin + 1.0
        norm = Normalize(vmin=vmin, vmax=vmax)

        legend_handles = []
        for i, case in enumerate(CASE_ORDER):
            g = local[local.case_label == case].copy()
            if g.empty:
                continue
            g["crack_extension_um"] = pd.to_numeric(g["crack_extension_um"], errors="coerce")
            g = g.sort_values(["crack_extension_um", "target_Kmax_MPa_sqrtm", "step"])
            marker = markers[i % len(markers)]
            ax.scatter(
                g["local_DeltaK_MPa_sqrtm"],
                g["local_da_dN_m_per_cycle"],
                c=g["crack_extension_um"],
                cmap=cmap,
                norm=norm,
                marker=marker,
                s=38,
                alpha=0.86,
                linewidths=0.35,
                edgecolors="black",
                zorder=3,
            )
            legend_handles.append(
                Line2D(
                    [0], [0], marker=marker, linestyle="none",
                    markerfacecolor="0.55", markeredgecolor="black",
                    markeredgewidth=0.5, markersize=7.5, label=case,
                )
            )

        sm = ScalarMappable(norm=norm, cmap=cmap)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, pad=0.018)
        cbar.set_label(r"Crack extension $\Delta a$ ($\mu$m)")
        if legend_handles:
            ax.legend(handles=legend_handles, frameon=True, title="Material class")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"Local block $\Delta K=(1-R)K_J$ (MPa $\sqrt{m}$)")
    ax.set_ylabel(r"Local block $da/dN=\Delta a/\Delta N$ (m/cycle)")
    ax.set_title("2-D FEM/CZM local block Paris points")
    ax.grid(True, which="major", alpha=0.25)
    fig.tight_layout()

    # Preserve the established filenames, and also write an explicit alias that
    # records the color variable used in the plot.
    names = [
        "atlas_2d_local_da_dN_vs_local_DeltaK",
        "fem_czm_local_da_dN_vs_local_DeltaK",
        "atlas_2d_local_da_dN_vs_local_DeltaK_colored_by_extension",
    ]
    for name in names:
        fig.savefig(out / f"{name}.png", dpi=320)
        fig.savefig(out / f"{name}.svg")
    plt.close(fig)


def plot_cycles(points: pd.DataFrame, out: Path):
    fig, ax = plt.subplots(figsize=(8.5, 6.0))
    cmap = plt.get_cmap("tab10")
    for i, case in enumerate(CASE_ORDER):
        g = points[points.case_label == case].sort_values("target_Kmax_MPa_sqrtm")
        if g.empty:
            continue
        ax.plot(g["target_Kmax_MPa_sqrtm"], g["cycles_total"], marker="o", color=cmap(i % 10), label=case)
    ax.set_yscale("log")
    ax.set_xlabel(r"Target $K_{max}$ (MPa $\sqrt{m}$)")
    ax.set_ylabel("Accumulated cycles")
    ax.set_title("Cycles to target, arrest, or run limit")
    ax.grid(True, which="major", alpha=0.25)
    ax.legend(frameon=True)
    fig.tight_layout()
    fig.savefig(out / "fem_czm_cycles_vs_Kmax.png", dpi=320)
    fig.savefig(out / "fem_czm_cycles_vs_Kmax.svg")
    plt.close(fig)



def estimate_rate_thresholds(points: pd.DataFrame, threshold_rate: float) -> pd.DataFrame:
    """Estimate rate-defined DeltaK thresholds by log-linear interpolation.

    Only measured integrated rates are used for the interpolation. Bounds are
    retained separately in the atlas and are not promoted to measurements.
    """
    rows = []
    for case in CASE_ORDER:
        g = points[(points.case_label == case) & (~points["is_censored_upper_bound"].astype(bool))].copy()
        g = g[np.isfinite(g["DeltaK_initial_MPa_sqrtm"]) & np.isfinite(g["plot_da_dN_m_per_cycle"]) & (g["plot_da_dN_m_per_cycle"] > 0)]
        g = g.sort_values("DeltaK_initial_MPa_sqrtm")
        est = np.nan
        loK = hiK = loR = hiR = np.nan
        status = "not_bracketed"
        if len(g) >= 2:
            x = g["DeltaK_initial_MPa_sqrtm"].to_numpy(float)
            y = g["plot_da_dN_m_per_cycle"].to_numpy(float)
            for i in range(len(g)-1):
                y0, y1 = y[i], y[i+1]
                if (y0 <= threshold_rate <= y1) or (y1 <= threshold_rate <= y0):
                    if y0 > 0 and y1 > 0 and y0 != y1:
                        frac = (np.log10(threshold_rate)-np.log10(y0))/(np.log10(y1)-np.log10(y0))
                        est = x[i] + frac*(x[i+1]-x[i])
                    else:
                        est = 0.5*(x[i]+x[i+1])
                    loK, hiK, loR, hiR = x[i], x[i+1], y0, y1
                    status = "interpolated_measured"
                    break
        rows.append({
            "case_label": case,
            "threshold_rate_m_per_cycle": threshold_rate,
            "DeltaK_threshold_MPa_sqrtm": est,
            "status": status,
            "lower_DeltaK_MPa_sqrtm": loK,
            "upper_DeltaK_MPa_sqrtm": hiK,
            "lower_rate_m_per_cycle": loR,
            "upper_rate_m_per_cycle": hiR,
            "n_measured_points": len(g),
        })
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--R", type=float, default=0.1)
    ap.add_argument("--cycles-max", type=float, default=2e14)
    ap.add_argument("--target-crack-extension-um", type=float, default=250.0)
    ap.add_argument("--threshold-rate", type=float, default=1e-12, help="Rate criterion for DeltaK threshold interpolation [m/cycle]")
    args = ap.parse_args()

    root = Path(args.root)
    out = root / "analysis"
    out.mkdir(parents=True, exist_ok=True)
    cases = load_case_table(root)

    integrated = read_integrated_points(root, args.R, args.cycles_max)
    local = read_local_block_points(root, args.R)
    if integrated.empty:
        print(f"No compare_summary.csv/case_summary.csv files found under {root}")
        return

    integrated.to_csv(out / "atlas_2d_paris_points.csv", index=False)
    integrated.to_csv(out / "fem_czm_paris_points.csv", index=False)
    local.to_csv(out / "atlas_2d_local_paris_points.csv", index=False)
    local.to_csv(out / "fem_czm_local_paris_points.csv", index=False)

    summ = summarize(integrated, local, cases)
    summ.to_csv(out / "fem_czm_case_summary.csv", index=False)

    thresholds = estimate_rate_thresholds(integrated, args.threshold_rate)
    thresholds.to_csv(out / "fem_czm_rate_defined_thresholds.csv", index=False)

    plot_integrated(integrated, out)
    plot_measured(integrated, out)
    plot_local(local, out)
    plot_cycles(integrated, out)

    n_measured = int((~integrated["is_censored_upper_bound"].astype(bool)).sum())
    print(f"WROTE {out}")
    print("\nRate-defined threshold estimates:")
    print(thresholds.to_string(index=False))
    print(f"Integrated points: {len(integrated)}; measured: {n_measured}; upper bounds: {len(integrated)-n_measured}")
    print(f"Local accepted growth blocks: {len(local)}")
    if n_measured == 0:
        print("NOTE: this dataset contains no accepted crack advances; it cannot yet reveal a Paris slope or threshold knee.")
    print(summ.to_string(index=False))


if __name__ == "__main__":
    main()
