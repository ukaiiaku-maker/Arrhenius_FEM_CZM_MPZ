#!/usr/bin/env python3
"""Censor-aware plotting for mixed-mode FEM/CZM v1 campaign results."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def parse_psi(tag: str) -> float:
    m = re.fullmatch(r"psi_([mp])(\d+)p(\d+)", tag)
    if not m:
        raise ValueError(tag)
    sign = -1.0 if m.group(1) == "m" else 1.0
    return sign * (int(m.group(2)) + int(m.group(3)) / 10.0)

def load(root: Path) -> pd.DataFrame:
    frames = []
    for p in sorted(root.glob("*/*/seed_*/mixed_mode_first_passage_summary.csv")):
        df = pd.read_csv(p)
        df["class"] = p.parents[2].name
        df["target_psi_deg"] = parse_psi(p.parents[1].name)
        df["seed"] = int(p.parent.name.split("_")[-1])
        frames.append(df)
    if not frames:
        raise RuntimeError(f"No summaries found under {root}")
    out = pd.concat(frames, ignore_index=True)
    out["event_observed"] = out["mode_classification"].eq("brittle")
    return out

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    data = load(args.root)
    data.to_csv(args.out / "mixed_mode_all_cases_censor_aware.csv", index=False)

    grouped = (
        data.groupby(["class", "target_psi_deg", "event_observed"], as_index=False)
        .agg(
            n=("seed", "count"),
            KI_mean=("KI_first_MPa_sqrt_m", "mean"),
            KII_mean=("KII_first_MPa_sqrt_m", "mean"),
            KJ_mean=("KJ_first_MPa_sqrt_m", "mean"),
            KJ_sd=("KJ_first_MPa_sqrt_m", "std"),
            Kopen_mean=("Kopen_maxhoop_first_MPa_sqrt_m", "mean"),
            Kopen_sd=("Kopen_maxhoop_first_MPa_sqrt_m", "std"),
            achieved_mean=("mode_phase_first_deg", "mean"),
            achieved_sd=("mode_phase_first_deg", "std"),
            kink_mean=("maxhoop_kink_first_deg", "mean"),
            kink_sd=("maxhoop_kink_first_deg", "std"),
        )
        .sort_values(["class", "target_psi_deg"])
    )
    grouped.to_csv(args.out / "mixed_mode_grouped_censor_aware.csv", index=False)

    event = data[data["event_observed"]]
    cens = data[~data["event_observed"]]
    eg = grouped[grouped["event_observed"]]
    cg = grouped[~grouped["event_observed"]]

    fig, ax = plt.subplots(figsize=(8, 7))
    ax.scatter(event["KI_first_MPa_sqrt_m"], event["KII_first_MPa_sqrt_m"], alpha=.45, label="First passage")
    ax.plot(eg["KI_mean"], eg["KII_mean"], marker="o", label="First-passage mean")
    ax.plot(cg["KI_mean"], cg["KII_mean"], marker="^", linestyle="--", label="Censored endpoint")
    ax.axhline(0, linewidth=.8); ax.axvline(0, linewidth=.8)
    ax.set_xlabel(r"$K_I$ [MPa$\sqrt{\mathrm{m}}$]")
    ax.set_ylabel(r"$K_{II}$ [MPa$\sqrt{\mathrm{m}}$]")
    ax.legend(); ax.grid(alpha=.25); fig.tight_layout()
    fig.savefig(args.out / "KI_KII_envelope_censor_aware.png", dpi=220)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.5, 5.7))
    ax.errorbar(eg["target_psi_deg"], eg["Kopen_mean"], yerr=eg["Kopen_sd"].fillna(0),
                marker="o", capsize=3, label="First-passage mean ± SD")
    ax.plot(cg["target_psi_deg"], cg["Kopen_mean"], marker="^", linestyle="--",
            label="Censored lower bound")
    ax.set_xlabel("Requested phase angle [deg]")
    ax.set_ylabel(r"$K_{\mathrm{open}}$ [MPa$\sqrt{\mathrm{m}}$]")
    ax.legend(); ax.grid(alpha=.25); fig.tight_layout()
    fig.savefig(args.out / "Kopen_vs_target_phase_censor_aware.png", dpi=220)
    plt.close(fig)

    print(f"Wrote plots to {args.out}")

if __name__ == "__main__":
    main()
