#!/usr/bin/env python3
"""Link temperature-strength-anomaly metrics to excess fatigue resistance.

Inputs:
  1. A DBTT/fatigue link table containing Kc and rate-defined DeltaK_th.
  2. Strength-anomaly metrics generated from the same emission-barrier family.

The script first removes the common Kc -> DeltaK_th trend at each threshold
criterion, then asks whether the residual fatigue resistance correlates with
emission-controlled strength-anomaly metrics.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def prepare_link(path: Path, criterion: float, exclude_case: str):
    df = pd.read_csv(path)
    df = df[np.isclose(df["da_dN_criterion_m_per_cycle"], criterion)].copy()
    df = df[(df["threshold_class"] == "bracketed") & np.isfinite(df["Kc_first_MPa_sqrtm"]) & np.isfinite(df["DeltaK_threshold_estimate_MPa_sqrtm"])].copy()
    if df.empty:
        raise RuntimeError("no bracketed threshold points for requested criterion")
    train = df[df["case_label"] != exclude_case] if exclude_case else df
    if len(train) < 2:
        train = df
    slope, intercept = np.polyfit(train["Kc_first_MPa_sqrtm"], train["DeltaK_threshold_estimate_MPa_sqrtm"], 1)
    df["DeltaK_common_trend_MPa_sqrtm"] = slope*df["Kc_first_MPa_sqrtm"] + intercept
    df["fatigue_excess_residual_MPa_sqrtm"] = df["DeltaK_threshold_estimate_MPa_sqrtm"] - df["DeltaK_common_trend_MPa_sqrtm"]
    return df, float(slope), float(intercept)


def fatigue_summary(df: pd.DataFrame):
    rows = []
    for keys, g in df.groupby(["case_label", "S_emit_kB"], sort=False):
        g = g.sort_values("T_K")
        lowT = g.iloc[np.argmin(g["T_K"].to_numpy())]
        vals = g["DeltaK_threshold_estimate_MPa_sqrtm"].to_numpy(float)
        temps = g["T_K"].to_numpy(float)
        ref = float(lowT["DeltaK_threshold_estimate_MPa_sqrtm"])
        target = 0.5*ref
        T50 = np.nan
        for a, b, Ta, Tb in zip(vals[:-1], vals[1:], temps[:-1], temps[1:]):
            if (a-target)*(b-target) <= 0 and abs(b-a) > 1e-14:
                T50 = Ta + (Tb-Ta)*(target-a)/(b-a)
                break
        rows.append({
            "case_label": keys[0],
            "S_emit_kB": keys[1],
            "mean_fatigue_excess_MPa_sqrtm": float(g["fatigue_excess_residual_MPa_sqrtm"].mean()),
            "max_fatigue_excess_MPa_sqrtm": float(g["fatigue_excess_residual_MPa_sqrtm"].max()),
            "min_fatigue_excess_MPa_sqrtm": float(g["fatigue_excess_residual_MPa_sqrtm"].min()),
            "fatigue_threshold_lowT_MPa_sqrtm": ref,
            "fatigue_persistence_T50_K": T50,
            "n_temperature_points": len(g),
        })
    return pd.DataFrame(rows)


def correlations(merged: pd.DataFrame):
    metrics = ["anomaly_amplitude_GPa", "positive_slope_area_GPa", "plateau_width_K", "T_peak_K"]
    targets = ["mean_fatigue_excess_MPa_sqrtm", "max_fatigue_excess_MPa_sqrtm", "fatigue_persistence_T50_K"]
    rows = []
    for rate, gr in merged.groupby("strain_rate_s-1"):
        for x in metrics:
            for y in targets:
                g = gr[np.isfinite(gr[x]) & np.isfinite(gr[y])]
                if len(g) >= 3:
                    p = g[x].corr(g[y], method="pearson")
                    s = g[x].corr(g[y], method="spearman")
                else:
                    p = s = np.nan
                rows.append({"strain_rate_s-1": rate, "strength_metric": x, "fatigue_metric": y, "n": len(g), "pearson_r": p, "spearman_rho": s})
    return pd.DataFrame(rows)


def plot_amp_vs_excess(merged: pd.DataFrame, out: Path):
    rates = sorted(merged["strain_rate_s-1"].unique())
    fig, ax = plt.subplots(figsize=(8.3, 6.0))
    markers = ["o", "s", "^", "D", "v", "P"]
    for i, rate in enumerate(rates):
        g = merged[np.isclose(merged["strain_rate_s-1"], rate)]
        ax.scatter(g["anomaly_amplitude_GPa"], g["mean_fatigue_excess_MPa_sqrtm"],
                   marker=markers[i % len(markers)], s=58, label=f"rate={rate:g} s$^{{-1}}$")
    ax.axhline(0.0, linewidth=1.0)
    ax.set_xlabel("Strength-anomaly amplitude (GPa)")
    ax.set_ylabel(r"Mean excess fatigue resistance $\langle\delta_{fatigue}\rangle$ (MPa $\sqrt{m}$)")
    ax.set_title("Emission-barrier strength anomaly versus excess fatigue resistance")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "strength_anomaly_vs_fatigue_excess.png", dpi=260)
    plt.close(fig)


def plot_Tpeak_vs_T50(merged: pd.DataFrame, out: Path):
    g = merged[np.isfinite(merged["T_peak_K"]) & np.isfinite(merged["fatigue_persistence_T50_K"])].copy()
    if g.empty:
        return
    fig, ax = plt.subplots(figsize=(7.2, 5.8))
    for rate, gr in g.groupby("strain_rate_s-1"):
        ax.scatter(gr["T_peak_K"], gr["fatigue_persistence_T50_K"], s=58, label=f"rate={rate:g} s$^{{-1}}$")
    lo = min(g["T_peak_K"].min(), g["fatigue_persistence_T50_K"].min())
    hi = max(g["T_peak_K"].max(), g["fatigue_persistence_T50_K"].max())
    ax.plot([lo,hi],[lo,hi], linestyle="--", linewidth=1.0, alpha=0.6)
    ax.set_xlabel("Strength-anomaly peak temperature (K)")
    ax.set_ylabel("Fatigue persistence T50 (K)")
    ax.set_title("Temperature alignment of strength anomaly and fatigue persistence")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "strength_Tpeak_vs_fatigue_persistence_T50.png", dpi=250)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--link-csv", required=True)
    ap.add_argument("--strength-metrics", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--criterion", type=float, default=1e-10)
    ap.add_argument("--exclude-common-trend-case", default="plastic_shielded_case64_M1")
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    link, slope, intercept = prepare_link(Path(args.link_csv), args.criterion, args.exclude_common_trend_case)
    fsum = fatigue_summary(link)
    strength = pd.read_csv(args.strength_metrics)
    merged = strength.merge(fsum, on=["case_label", "S_emit_kB"], how="inner")
    corr = correlations(merged)

    link.to_csv(out / "fatigue_link_with_common_trend_residual.csv", index=False)
    fsum.to_csv(out / "fatigue_residual_summary_by_emission_barrier.csv", index=False)
    merged.to_csv(out / "strength_fatigue_merged_metrics.csv", index=False)
    corr.to_csv(out / "strength_fatigue_correlation_summary.csv", index=False)
    with open(out / "strength_fatigue_analysis_settings.json", "w") as f:
        json.dump({
            "criterion": args.criterion,
            "common_trend_slope": slope,
            "common_trend_intercept": intercept,
            "excluded_common_trend_case": args.exclude_common_trend_case,
        }, f, indent=2)
    if not merged.empty:
        plot_amp_vs_excess(merged, out)
        plot_Tpeak_vs_T50(merged, out)
    print(f"wrote {out / 'strength_fatigue_correlation_summary.csv'}")
    print(corr.sort_values("pearson_r", key=lambda s: s.abs(), ascending=False).head(12).to_string(index=False))


if __name__ == "__main__":
    main()
