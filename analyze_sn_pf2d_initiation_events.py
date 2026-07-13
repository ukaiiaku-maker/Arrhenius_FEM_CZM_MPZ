#!/usr/bin/env python3
"""Post-process 2-D Arrhenius fatigue histories into staged initiation metrics.

This script does not alter the physics solver.  It separates:
  1) embryo first-passage events from cumulative nucleation hazard B_nuc,
  2) phase-field localization thresholds from d_max,
  3) resolved crack formation from connected crack extent relative to ell.

It can be run directly on an existing production output directory.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def first_crossing(x: np.ndarray, y: np.ndarray, threshold: float) -> float:
    """Linearly interpolate x where monotone-ish y first crosses threshold."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if len(x) == 0:
        return math.nan
    hit = np.flatnonzero(y >= threshold)
    if len(hit) == 0:
        return math.nan
    i = int(hit[0])
    if i == 0:
        return float(x[0])
    x0, x1 = x[i - 1], x[i]
    y0, y1 = y[i - 1], y[i]
    if not np.isfinite(y0) or y1 == y0:
        return float(x1)
    f = (threshold - y0) / (y1 - y0)
    f = float(np.clip(f, 0.0, 1.0))
    return float(x0 + f * (x1 - x0))


def read_case(case_dir: Path) -> tuple[pd.DataFrame, dict]:
    hist = pd.read_csv(case_dir / "sn_pf2d_fullplastic_history.csv")
    with open(case_dir / "summary.json", "r", encoding="utf-8") as f:
        summary = json.load(f)
    return hist, summary


def iter_case_dirs(root: Path) -> Iterable[Path]:
    for p in sorted(root.glob("*/sigmaA_*MPa")):
        if (p / "sn_pf2d_fullplastic_history.csv").is_file() and (p / "summary.json").is_file():
            yield p


def analyze_case(case_dir: Path, embryo_probs: list[float], d_thresholds: list[float], extent_factors: list[float]) -> dict:
    h, s = read_case(case_dir)
    N = h["cycles_total"].to_numpy(float)
    ell = float(s["ell_m"])
    row = {
        "case": str(s["case"]),
        "sigma_a_MPa": float(s["sigma_a_MPa"]),
        "T_K": float(s.get("T_K", np.nan)),
        "cycles_horizon": float(N[-1]) if len(N) else np.nan,
        "ell_m": ell,
        "B_nuc_final_max": float(h["B_nuc_max"].iloc[-1]),
        "d_final_max": float(h["d_max"].iloc[-1]),
        "connected_extent_final_m": float(h["connected_crack_extent_m"].iloc[-1]),
        "connected_extent_final_over_ell": float(h["connected_crack_extent_m"].iloc[-1] / ell),
    }

    for q in embryo_probs:
        Bq = -math.log(max(1.0 - q, 1e-15))
        row[f"N_embryo_q{int(round(100*q)):02d}"] = first_crossing(N, h["B_nuc_max"].to_numpy(float), Bq)
        row[f"B_threshold_q{int(round(100*q)):02d}"] = Bq
    row["N_embryo_B1"] = first_crossing(N, h["B_nuc_max"].to_numpy(float), 1.0)

    for dcrit in d_thresholds:
        tag = str(dcrit).replace(".", "p")
        row[f"N_dmax_{tag}"] = first_crossing(N, h["d_max"].to_numpy(float), dcrit)

    for fac in extent_factors:
        tag = str(fac).replace(".", "p")
        row[f"N_extent_{tag}ell"] = first_crossing(
            N,
            h["connected_crack_extent_m"].to_numpy(float),
            fac * ell,
        )

    # Default staged event definitions for convenient plotting.
    row["N_embryo_default"] = row["N_embryo_B1"]
    row["N_localized_default"] = row.get("N_dmax_0p5", np.nan)
    row["N_resolved_default"] = row.get("N_extent_3p0ell", np.nan)
    if np.isfinite(row["N_embryo_default"]) and np.isfinite(row["N_resolved_default"]):
        row["small_crack_establishment_cycles"] = row["N_resolved_default"] - row["N_embryo_default"]
        row["resolved_over_embryo_life_ratio"] = row["N_resolved_default"] / max(row["N_embryo_default"], 1e-30)
    else:
        row["small_crack_establishment_cycles"] = np.nan
        row["resolved_over_embryo_life_ratio"] = np.nan
    return row


def plot_sn(events: pd.DataFrame, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 5.4))
    for case, d in events.groupby("case", sort=True):
        d = d.sort_values("sigma_a_MPa")
        ax.plot(d["N_embryo_default"], d["sigma_a_MPa"], marker="o", linestyle="--", label=f"{case}: embryo B=1")
        ax.plot(d["N_resolved_default"], d["sigma_a_MPa"], marker="s", linestyle="-", label=f"{case}: resolved 3ell")
    ax.set_xscale("log")
    ax.set_xlabel("Cycles")
    ax.set_ylabel("Stress amplitude (MPa)")
    ax.set_title("Staged 2-D fatigue initiation: embryo vs resolved crack")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "SN_embryo_vs_resolved.png", dpi=220)
    plt.close(fig)


def plot_extent_sensitivity(events: pd.DataFrame, factors: list[float], out: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 5.4))
    for case, d in events.groupby("case", sort=True):
        for fac in factors:
            tag = str(fac).replace(".", "p")
            col = f"N_extent_{tag}ell"
            if col not in d:
                continue
            dd = d.sort_values("sigma_a_MPa")
            ax.plot(dd[col], dd["sigma_a_MPa"], marker="o", label=f"{case}, {fac:g}ell")
    ax.set_xscale("log")
    ax.set_xlabel("Cycles to connected crack extent")
    ax.set_ylabel("Stress amplitude (MPa)")
    ax.set_title("Resolved-crack S-N sensitivity to extent criterion")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(out / "SN_extent_factor_sensitivity.png", dpi=220)
    plt.close(fig)


def plot_stage_bars(events: pd.DataFrame, out: Path) -> None:
    e = events.sort_values(["case", "sigma_a_MPa"]).reset_index(drop=True)
    labels = [f"{r.case}\n{r.sigma_a_MPa:g} MPa" for r in e.itertuples()]
    x = np.arange(len(e))
    fig, ax = plt.subplots(figsize=(10.0, 5.5))
    for col, marker, label in [
        ("N_embryo_default", "o", "embryo B=1"),
        ("N_localized_default", "^", "dmax=0.5"),
        ("N_resolved_default", "s", "connected extent=3ell"),
    ]:
        ax.scatter(x, e[col], marker=marker, s=55, label=label)
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=35, ha="right")
    ax.set_ylabel("Cycles to event")
    ax.set_title("Fatigue initiation stages by case")
    ax.grid(True, axis="y", which="both", alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / "initiation_stage_event_map.png", dpi=220)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("root", type=Path, help="Production output root containing case/sigmaA_*MPa directories")
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--embryo-probs", nargs="+", type=float, default=[0.5, 0.6321205588, 0.9])
    ap.add_argument("--d-thresholds", nargs="+", type=float, default=[0.4, 0.5, 0.6, 0.7])
    ap.add_argument("--extent-factors", nargs="+", type=float, default=[1.0, 2.0, 3.0, 4.0])
    args = ap.parse_args()

    root = args.root.resolve()
    out = args.out.resolve() if args.out else root / "initiation_event_analysis"
    out.mkdir(parents=True, exist_ok=True)

    rows = []
    for case_dir in iter_case_dirs(root):
        rows.append(analyze_case(case_dir, args.embryo_probs, args.d_thresholds, args.extent_factors))
    if not rows:
        raise SystemExit(f"No case histories found under {root}")

    events = pd.DataFrame(rows).sort_values(["case", "sigma_a_MPa"])
    events.to_csv(out / "staged_initiation_events.csv", index=False)

    # Long-form sensitivity table for easier plotting/statistics.
    sens_rows = []
    for r in events.to_dict("records"):
        for fac in args.extent_factors:
            tag = str(fac).replace(".", "p")
            sens_rows.append({
                "case": r["case"],
                "sigma_a_MPa": r["sigma_a_MPa"],
                "extent_factor_ell": fac,
                "cycles_to_extent": r.get(f"N_extent_{tag}ell", np.nan),
            })
    pd.DataFrame(sens_rows).to_csv(out / "extent_factor_sensitivity.csv", index=False)

    plot_sn(events, out)
    plot_extent_sensitivity(events, args.extent_factors, out)
    plot_stage_bars(events, out)

    print(events[[
        "case", "sigma_a_MPa", "N_embryo_default", "N_localized_default",
        "N_resolved_default", "connected_extent_final_over_ell",
    ]].to_string(index=False))
    print(f"\nWrote staged initiation analysis to: {out}")


if __name__ == "__main__":
    main()
