#!/usr/bin/env python3
"""Analyze staged fatigue-initiation events from 2-D full-plastic PF runs.

The solver physics is not modified.  Histories are post-processed into four
separate event families:
  * crack-opening/embryo hazard quantiles from B_nuc,
  * local PF localization from d_max,
  * connected crack formation from a_conn >= 2 ell,
  * established crack formation from a_conn >= 3 ell.

All absent events are retained as right-censored observations at the run horizon.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def first_crossing(x: np.ndarray, y: np.ndarray, threshold: float) -> float:
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
    if y1 == y0:
        return float(x1)
    f = float(np.clip((threshold - y0) / (y1 - y0), 0.0, 1.0))
    return float(x0 + f * (x1 - x0))


def iter_case_dirs(root: Path) -> Iterable[Path]:
    for p in sorted(root.glob("*/sigmaA_*MPa")):
        if (p / "sn_pf2d_fullplastic_history.csv").is_file() and (p / "summary.json").is_file():
            yield p


def read_case(case_dir: Path) -> tuple[pd.DataFrame, dict]:
    h = pd.read_csv(case_dir / "sn_pf2d_fullplastic_history.csv")
    with (case_dir / "summary.json").open("r", encoding="utf-8") as f:
        s = json.load(f)
    required = ["cycles_total", "B_nuc_max", "d_max", "connected_crack_extent_m"]
    missing = [c for c in required if c not in h.columns]
    if missing:
        raise ValueError(f"{case_dir}: missing history columns {missing}")
    return h, s


def status(n_event: float) -> str:
    return "event" if np.isfinite(n_event) else "right_censored"


def analyze_case(case_dir: Path, d_thresholds: list[float], extent_factors: list[float]) -> dict:
    h, s = read_case(case_dir)
    N = h["cycles_total"].to_numpy(float)
    B = h["B_nuc_max"].to_numpy(float)
    d = h["d_max"].to_numpy(float)
    a = h["connected_crack_extent_m"].to_numpy(float)
    ell = float(s["ell_m"])
    horizon = float(N[-1])

    row = {
        "case": str(s["case"]),
        "sigma_a_MPa": float(s["sigma_a_MPa"]),
        "T_K": float(s.get("T_K", np.nan)),
        "cycles_horizon": horizon,
        "ell_m": ell,
        "B_nuc_final_max": float(B[-1]),
        "d_final_max": float(d[-1]),
        "connected_extent_final_m": float(a[-1]),
        "connected_extent_final_over_ell": float(a[-1] / ell),
    }

    for q in (0.5, 0.6321205588, 0.9):
        Bq = -math.log(1.0 - q)
        tag = int(round(100 * q))
        row[f"N_hazard_q{tag:02d}"] = first_crossing(N, B, Bq)
        row[f"B_threshold_q{tag:02d}"] = Bq
    row["N_hazard_B1"] = first_crossing(N, B, 1.0)

    for dc in d_thresholds:
        tag = str(dc).replace(".", "p")
        row[f"N_dmax_{tag}"] = first_crossing(N, d, dc)

    for fac in extent_factors:
        tag = str(fac).replace(".", "p")
        row[f"N_extent_{tag}ell"] = first_crossing(N, a, fac * ell)

    # Staged definitions used for comparison, without changing solver physics.
    row["N_localization"] = row.get("N_dmax_0p5", np.nan)
    row["N_hazard_reference"] = row["N_hazard_B1"]
    row["N_connected_crack"] = row.get("N_extent_2p0ell", np.nan)
    row["N_established_crack"] = row.get("N_extent_3p0ell", np.nan)

    for name in ("localization", "hazard_reference", "connected_crack", "established_crack"):
        row[f"{name}_status"] = status(row[f"N_{name}"])

    if np.isfinite(row["N_hazard_reference"]) and np.isfinite(row["N_connected_crack"]):
        row["cycles_hazard_to_connected"] = row["N_connected_crack"] - row["N_hazard_reference"]
    else:
        row["cycles_hazard_to_connected"] = np.nan
    if np.isfinite(row["N_connected_crack"]) and np.isfinite(row["N_established_crack"]):
        row["cycles_connected_to_established"] = row["N_established_crack"] - row["N_connected_crack"]
    else:
        row["cycles_connected_to_established"] = np.nan
    return row


def plot_stage_sn(events: pd.DataFrame, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.0, 5.8))
    stages = [
        ("N_hazard_reference", "--", "o", "hazard B=1"),
        ("N_localization", ":", "^", "localization dmax=0.5"),
        ("N_connected_crack", "-", "s", "connected crack 2ell"),
        ("N_established_crack", "-.", "D", "established crack 3ell"),
    ]
    for case, dd in events.groupby("case", sort=True):
        dd = dd.sort_values("sigma_a_MPa")
        for col, ls, marker, label in stages:
            ax.plot(dd[col], dd["sigma_a_MPa"], linestyle=ls, marker=marker,
                    label=f"{case}: {label}")
    ax.set_xscale("log")
    ax.set_xlabel("Cycles to event")
    ax.set_ylabel("Stress amplitude (MPa)")
    ax.set_title("Staged 2-D fatigue events")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(out / "SN_staged_events.png", dpi=220)
    plt.close(fig)


def plot_crack_sn(events: pd.DataFrame, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.4, 5.5))
    for case, dd in events.groupby("case", sort=True):
        dd = dd.sort_values("sigma_a_MPa")
        ax.plot(dd["N_connected_crack"], dd["sigma_a_MPa"], marker="s", linestyle="-", label=f"{case}: 2ell")
        ax.plot(dd["N_established_crack"], dd["sigma_a_MPa"], marker="D", linestyle="--", label=f"{case}: 3ell")
        # Right-censored connected-crack cases plotted at their horizons with open x markers.
        cens = dd[dd["connected_crack_status"] == "right_censored"]
        if not cens.empty:
            ax.scatter(cens["cycles_horizon"], cens["sigma_a_MPa"], marker="x", s=45,
                       label=f"{case}: connected crack censored")
    ax.set_xscale("log")
    ax.set_xlabel("Cycles")
    ax.set_ylabel("Stress amplitude (MPa)")
    ax.set_title("Connected and established crack S-N response")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "SN_connected_and_established_crack.png", dpi=220)
    plt.close(fig)


def plot_final_state(events: pd.DataFrame, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    for case, dd in events.groupby("case", sort=True):
        dd = dd.sort_values("sigma_a_MPa")
        ax.plot(dd["sigma_a_MPa"], dd["connected_extent_final_over_ell"], marker="o", label=case)
    ax.axhline(2.0, linestyle="--", linewidth=1.0, label="connected crack: 2ell")
    ax.axhline(3.0, linestyle=":", linewidth=1.0, label="established crack: 3ell")
    ax.set_xlabel("Stress amplitude (MPa)")
    ax.set_ylabel("Final connected extent / ell")
    ax.set_title("Final connected localization extent")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "final_connected_extent_vs_stress.png", dpi=220)
    plt.close(fig)


def write_sensitivity(events: pd.DataFrame, extent_factors: list[float], out: Path) -> None:
    rows = []
    for r in events.to_dict("records"):
        for fac in extent_factors:
            tag = str(fac).replace(".", "p")
            n = r.get(f"N_extent_{tag}ell", np.nan)
            rows.append({
                "case": r["case"],
                "sigma_a_MPa": r["sigma_a_MPa"],
                "extent_factor_ell": fac,
                "cycles_to_extent": n,
                "status": status(n),
                "cycles_horizon": r["cycles_horizon"],
            })
    pd.DataFrame(rows).to_csv(out / "extent_factor_sensitivity.csv", index=False)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("root", type=Path)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--d-thresholds", nargs="+", type=float, default=[0.4, 0.5, 0.6, 0.7])
    ap.add_argument("--extent-factors", nargs="+", type=float, default=[1.0, 2.0, 3.0, 4.0])
    args = ap.parse_args()

    root = args.root.resolve()
    out = args.out.resolve() if args.out else root / "staged_initiation_analysis_v2"
    out.mkdir(parents=True, exist_ok=True)

    rows = [analyze_case(p, args.d_thresholds, args.extent_factors) for p in iter_case_dirs(root)]
    if not rows:
        raise SystemExit(f"No case histories found under {root}")
    events = pd.DataFrame(rows).sort_values(["case", "sigma_a_MPa"]).reset_index(drop=True)
    events.to_csv(out / "staged_initiation_events_v2.csv", index=False)
    write_sensitivity(events, args.extent_factors, out)
    plot_stage_sn(events, out)
    plot_crack_sn(events, out)
    plot_final_state(events, out)

    cols = ["case", "sigma_a_MPa", "N_hazard_reference", "N_localization",
            "N_connected_crack", "N_established_crack", "connected_extent_final_over_ell"]
    print(events[cols].to_string(index=False))
    print(f"\nWrote staged analysis to: {out}")


if __name__ == "__main__":
    main()
