#!/usr/bin/env python3
"""Combine 1x/10x/100x four-class sweep summaries and make rate-comparison plots."""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def label(f: float) -> str:
    return f"rate_{int(f)}x" if float(f).is_integer() else f"rate_{f:g}x"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--rate-factors", default="1 10 100")
    args = ap.parse_args()
    root = Path(args.root)
    factors = [float(x) for x in args.rate_factors.replace(",", " ").split()]

    parts = []
    for f in factors:
        sf = root / label(f) / "four_class_temperature_summary.csv"
        if not sf.exists():
            print(f"WARNING missing {sf}")
            continue
        d = pd.read_csv(sf)
        d.insert(0, "rate_factor", f)
        parts.append(d)
    if not parts:
        raise SystemExit("No per-rate summary files found")
    df = pd.concat(parts, ignore_index=True, sort=False)
    df.to_csv(root / "rate_temperature_summary.csv", index=False)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    classes = list(dict.fromkeys(df["class"].astype(str)))
    nrows, ncols = 2, 2

    fig, axes = plt.subplots(nrows, ncols, figsize=(10.5, 7.6), sharex=True)
    for ax, klass in zip(axes.ravel(), classes):
        g0 = df[df["class"] == klass]
        for f in factors:
            g = g0[g0.rate_factor == f].sort_values("T_K")
            if len(g):
                ax.plot(g.T_K, g.Kc_first_MPa_sqrt_m, marker="o", label=f"{f:g}x")
        ax.set_title(klass)
        ax.set_ylabel(r"Initiation $K_c$ (MPa$\sqrt{m}$)")
        ax.grid(alpha=0.25)
    for ax in axes[-1, :]:
        ax.set_xlabel("Temperature (K)")
    axes[0, 0].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(root / "four_class_rate_effect_Kinit_vs_T.png", dpi=220)
    plt.close(fig)

    if "Kprop_200_500um_median" in df.columns:
        fig, axes = plt.subplots(nrows, ncols, figsize=(10.5, 7.6), sharex=True)
        for ax, klass in zip(axes.ravel(), classes):
            g0 = df[df["class"] == klass]
            for f in factors:
                g = g0[g0.rate_factor == f].sort_values("T_K")
                if len(g):
                    ax.plot(g.T_K, g.Kprop_200_500um_median, marker="o", label=f"{f:g}x")
            ax.set_title(klass)
            ax.set_ylabel(r"Propagation $K$ (MPa$\sqrt{m}$)")
            ax.grid(alpha=0.25)
        for ax in axes[-1, :]:
            ax.set_xlabel("Temperature (K)")
        axes[0, 0].legend(frameon=False)
        fig.tight_layout()
        fig.savefig(root / "four_class_rate_effect_Kprop_vs_T.png", dpi=220)
        plt.close(fig)

    print(f"WROTE {root / 'rate_temperature_summary.csv'}")


if __name__ == "__main__":
    main()
