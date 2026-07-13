#!/usr/bin/env python3
"""Replot a Panel C V1 S-N sweep using Lambda_sh as the waterfall axis.

The input CSV is the raw output from build_panelC_SN_nucleation_waterfall_3d.py.
The script selects one fixed crack-entropy slice and plots

    x = log10 cycles to initiation
    y = Lambda_sh = Gshield / G0_emit(T)
    z = stress amplitude sigma_a [MPa]

Status conventions:
  failed          : filled circle and connected S-N curve
  right_censored  : open triangle at the achieved/censoring life
  block_limited   : open diamond at cycles_total (not treated as Nmax censoring)
  other           : x marker
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--raw-csv", type=Path, required=True)
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--entropy-mag-kB", type=float, default=0.0,
                   help="Fixed entropy-magnitude slice to plot (default: 0).")
    p.add_argument("--cycles-max", type=float, default=None,
                   help="Optional display horizon; inferred from data if omitted.")
    p.add_argument("--elev", type=float, default=24.0)
    p.add_argument("--azim", type=float, default=-58.0)
    p.add_argument("--dpi", type=int, default=300)
    return p.parse_args()


def _display_log10_cycles(df: pd.DataFrame, cycles_max: float) -> np.ndarray:
    if "log10_cycles_display" in df.columns:
        vals = pd.to_numeric(df["log10_cycles_display"], errors="coerce").to_numpy(float)
        if np.all(np.isfinite(vals)):
            return vals

    status = df["status"].astype(str).to_numpy()
    failed_cycles = pd.to_numeric(df.get("cycles_to_nucleation"), errors="coerce").to_numpy(float)
    total_cycles = pd.to_numeric(df.get("cycles_total"), errors="coerce").to_numpy(float)
    life = np.where(status == "failed", failed_cycles, total_cycles)
    life = np.where(np.isfinite(life) & (life > 0), life, cycles_max)
    return np.log10(np.maximum(life, 1e-300))


def select_slice(raw: pd.DataFrame, entropy_mag_kB: float) -> pd.DataFrame:
    if "Lambda_sh" not in raw.columns:
        raise ValueError("Input CSV does not contain required column 'Lambda_sh'.")
    if "entropy_mag_kB" not in raw.columns:
        raise ValueError("Input CSV does not contain required column 'entropy_mag_kB'.")

    available = np.array(sorted(pd.to_numeric(raw["entropy_mag_kB"], errors="coerce").dropna().unique()), float)
    if available.size == 0:
        raise ValueError("No finite entropy_mag_kB values were found.")
    idx = int(np.argmin(np.abs(available - float(entropy_mag_kB))))
    chosen = float(available[idx])
    if not math.isclose(chosen, float(entropy_mag_kB), rel_tol=0.0, abs_tol=1e-10):
        raise ValueError(
            f"Requested entropy_mag_kB={entropy_mag_kB:g} is not present. "
            f"Available values: {available.tolist()}"
        )
    out = raw[np.isclose(pd.to_numeric(raw["entropy_mag_kB"], errors="coerce").to_numpy(float), chosen,
                         atol=1e-10, rtol=0.0)].copy()
    if out.empty:
        raise ValueError("Selected entropy slice is empty.")
    return out


def stress_at_target_life(g: pd.DataFrame, target_N: float) -> float:
    f = g[g["status"].astype(str) == "failed"].copy()
    if len(f) < 2:
        return math.nan
    x = np.log10(pd.to_numeric(f["cycles_to_nucleation"], errors="coerce").to_numpy(float))
    y = pd.to_numeric(f["sigma_a_MPa"], errors="coerce").to_numpy(float)
    good = np.isfinite(x) & np.isfinite(y)
    x, y = x[good], y[good]
    if len(x) < 2:
        return math.nan
    order = np.argsort(x)
    x, y = x[order], y[order]
    xt = math.log10(target_N)
    if xt < np.nanmin(x) or xt > np.nanmax(x):
        return math.nan
    return float(np.interp(xt, x, y))


def make_summary(df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, float | int]] = []
    for lam, g in df.groupby("Lambda_sh", sort=True):
        status = g["status"].astype(str)
        row = {
            "Lambda_sh": float(lam),
            "Gshield_eV": float(pd.to_numeric(g["Gshield_eV"], errors="coerce").dropna().iloc[0])
                if "Gshield_eV" in g and pd.to_numeric(g["Gshield_eV"], errors="coerce").notna().any() else math.nan,
            "entropy_mag_kB": float(pd.to_numeric(g["entropy_mag_kB"], errors="coerce").dropna().iloc[0]),
            "n_points": int(len(g)),
            "n_failed": int((status == "failed").sum()),
            "n_right_censored": int((status == "right_censored").sum()),
            "n_block_limited": int((status == "block_limited").sum()),
            "sigma_at_1e6_cycles_MPa": stress_at_target_life(g, 1e6),
            "sigma_at_1e8_cycles_MPa": stress_at_target_life(g, 1e8),
            "sigma_at_1e10_cycles_MPa": stress_at_target_life(g, 1e10),
            "sigma_at_1e12_cycles_MPa": stress_at_target_life(g, 1e12),
        }
        rows.append(row)
    return pd.DataFrame(rows)


def plot_3d(df: pd.DataFrame, outdir: Path, cycles_max: float, elev: float, azim: float, dpi: int) -> None:
    fig = plt.figure(figsize=(10.2, 7.8))
    ax = fig.add_subplot(111, projection="3d")

    lamvals = np.array(sorted(pd.to_numeric(df["Lambda_sh"], errors="coerce").dropna().unique()), float)
    vmin, vmax = float(lamvals.min()), float(lamvals.max())
    norm = Normalize(vmin=vmin, vmax=(vmax if not math.isclose(vmin, vmax) else vmin + 1.0))
    cmap = plt.get_cmap("viridis")

    all_x: list[float] = []
    for lam in lamvals:
        g = df[np.isclose(pd.to_numeric(df["Lambda_sh"], errors="coerce").to_numpy(float), lam,
                          atol=1e-10, rtol=0.0)].copy()
        g = g.sort_values("sigma_a_MPa")
        color = cmap(norm(float(lam)))
        x = _display_log10_cycles(g, cycles_max)
        y = np.full(len(g), float(lam), float)
        z = pd.to_numeric(g["sigma_a_MPa"], errors="coerce").to_numpy(float)
        status = g["status"].astype(str).to_numpy()
        all_x.extend(x[np.isfinite(x)].tolist())

        failed = status == "failed"
        if np.any(failed):
            idx = np.where(failed)[0]
            ax.plot(x[idx], y[idx], z[idx], lw=2.0, color=color, alpha=0.95)
            ax.scatter(x[idx], y[idx], z[idx], s=18, color=[color], alpha=0.98)

        cens = status == "right_censored"
        if np.any(cens):
            idx = np.where(cens)[0]
            ax.scatter(x[idx], y[idx], z[idx], marker="^", s=42,
                       facecolors="none", edgecolors=[color], linewidths=1.1)

        block = status == "block_limited"
        if np.any(block):
            idx = np.where(block)[0]
            ax.scatter(x[idx], y[idx], z[idx], marker="D", s=34,
                       facecolors="none", edgecolors=[color], linewidths=1.1)

        other = ~(failed | cens | block)
        if np.any(other):
            idx = np.where(other)[0]
            ax.scatter(x[idx], y[idx], z[idx], marker="x", s=32,
                       color=[color], linewidths=1.0)

        if len(g):
            q = g.iloc[-1]
            qx = float(_display_log10_cycles(pd.DataFrame([q]), cycles_max)[0])
            ax.text(qx + 0.06, float(lam), float(q["sigma_a_MPa"]),
                    rf"$\Lambda_{{sh}}$={lam:.2f}", fontsize=8, color=color)

    ax.set_xlabel(r"$\log_{10}$ cycles to crack initiation, $N_i$", labelpad=10)
    ax.set_ylabel(r"$\Lambda_{sh}=G_{shield}/G_{0,e}(T)$", labelpad=10)
    ax.set_zlabel(r"Stress amplitude $\sigma_a$ [MPa]", labelpad=11)
    entropy = float(pd.to_numeric(df["entropy_mag_kB"], errors="coerce").dropna().iloc[0])
    ax.set_title(
        "Systematic V1 S–N crack-initiation waterfall\n"
        + rf"fixed $-S_{{crack}}^*/k_B={entropy:g}$; shielding axis $\Lambda_{{sh}}$",
        pad=18,
    )
    ax.view_init(elev=float(elev), azim=float(azim))

    xmax = math.log10(max(float(cycles_max), 1.0))
    xmin = max(0.0, min(all_x) - 0.25) if all_x else 0.0
    ax.set_xlim(xmin, xmax + 0.15)

    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.045, pad=0.10)
    cbar.set_label(r"$\Lambda_{sh}=G_{shield}/G_{0,e}(T)$")

    fig.text(0.07, 0.048,
             "Open triangle: right-censored; open diamond: block-limited (not treated as censoring)",
             fontsize=9.2)
    fig.tight_layout(rect=(0, 0.07, 0.93, 1))
    fig.savefig(outdir / "panelC_SN_lambda_waterfall_3d.png", dpi=int(dpi), bbox_inches="tight")
    fig.savefig(outdir / "panelC_SN_lambda_waterfall_3d.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_2d(df: pd.DataFrame, outdir: Path, dpi: int) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    lamvals = np.array(sorted(pd.to_numeric(df["Lambda_sh"], errors="coerce").dropna().unique()), float)
    vmin, vmax = float(lamvals.min()), float(lamvals.max())
    norm = Normalize(vmin=vmin, vmax=(vmax if not math.isclose(vmin, vmax) else vmin + 1.0))
    cmap = plt.get_cmap("viridis")

    for lam in lamvals:
        g = df[np.isclose(pd.to_numeric(df["Lambda_sh"], errors="coerce").to_numpy(float), lam,
                          atol=1e-10, rtol=0.0)].sort_values("sigma_a_MPa")
        color = cmap(norm(float(lam)))
        failed = g["status"].astype(str).to_numpy() == "failed"
        if np.any(failed):
            ax.plot(g.loc[failed, "cycles_to_nucleation"], g.loc[failed, "sigma_a_MPa"],
                    marker="o", ms=4, lw=1.7, color=color,
                    label=rf"$\Lambda_{{sh}}$={lam:.2f}")
        cens = g["status"].astype(str).to_numpy() == "right_censored"
        if np.any(cens):
            ax.scatter(g.loc[cens, "cycles_total"], g.loc[cens, "sigma_a_MPa"],
                       marker="^", s=38, facecolors="none", edgecolors=[color])
        block = g["status"].astype(str).to_numpy() == "block_limited"
        if np.any(block):
            ax.scatter(g.loc[block, "cycles_total"], g.loc[block, "sigma_a_MPa"],
                       marker="D", s=32, facecolors="none", edgecolors=[color])

    ax.set_xscale("log")
    ax.set_xlabel("Cycles to crack initiation")
    ax.set_ylabel(r"Stress amplitude $\sigma_a$ [MPa]")
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "panelC_SN_lambda_projection_diagnostic.png", dpi=int(dpi), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    raw = pd.read_csv(args.raw_csv)
    df = select_slice(raw, args.entropy_mag_kB)
    args.outdir.mkdir(parents=True, exist_ok=True)

    if args.cycles_max is not None:
        cycles_max = float(args.cycles_max)
    else:
        vals = pd.to_numeric(df.get("cycles_total"), errors="coerce")
        cycles_max = float(np.nanmax(vals.to_numpy(float))) if vals.notna().any() else 1e12
        cycles_max = max(cycles_max, 1.0)

    summary = make_summary(df)
    summary.to_csv(args.outdir / "panelC_SN_lambda_curve_summary.csv", index=False)
    df.to_csv(args.outdir / "panelC_SN_lambda_selected_slice.csv", index=False)
    plot_3d(df, args.outdir, cycles_max, args.elev, args.azim, args.dpi)
    plot_2d(df, args.outdir, args.dpi)

    print(f"Selected entropy_mag_kB={args.entropy_mag_kB:g}")
    print(f"Wrote: {args.outdir / 'panelC_SN_lambda_waterfall_3d.png'}")
    print(f"Wrote: {args.outdir / 'panelC_SN_lambda_waterfall_3d.pdf'}")
    print(f"Wrote: {args.outdir / 'panelC_SN_lambda_projection_diagnostic.png'}")
    print(f"Wrote: {args.outdir / 'panelC_SN_lambda_curve_summary.csv'}")


if __name__ == "__main__":
    main()
