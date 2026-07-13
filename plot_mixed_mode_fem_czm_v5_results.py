#!/usr/bin/env python3
"""Censor-aware plots for mixed-mode FEM/CZM v5."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def finish(fig, path):
    fig.tight_layout()
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--root", required=True, type=Path)
    p.add_argument("--out", type=Path)
    a = p.parse_args()
    out = a.out or (a.root/"plots_v5")
    out.mkdir(parents=True, exist_ok=True)
    csv_path = a.root/"mixed_mode_v5_anisotropic_all_cases.csv"
    if not csv_path.exists():
        raise SystemExit(f"missing {csv_path}")
    df = pd.read_csv(csv_path)
    df["event_observed"] = df["status"].eq("event") | df["status"].eq("event_phase_mismatch")
    df.to_csv(out/"mixed_mode_v5_plot_data.csv", index=False)

    # 1. Event KJ and censored endpoint KJ.
    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    for cls, g in df.groupby("class"):
        g = g.sort_values("target_psi_deg")
        ev = g[g["event_observed"]]
        ce = g[~g["event_observed"]]
        if len(ev):
            ax.plot(ev["target_psi_deg"], ev["KJ_reference_first_MPa_sqrt_m"],
                    marker="o", label=f"{cls}: first passage")
        if len(ce):
            ax.scatter(ce["target_psi_deg"], ce["KJ_reference_first_MPa_sqrt_m"],
                       marker="^", label=f"{cls}: censored endpoint")
    ax.set_xlabel("Target process-zone traction phase [deg]")
    ax.set_ylabel(r"$K_J$ [MPa$\sqrt{\mathrm{m}}$]")
    ax.set_title("Anisotropic first passage and right-censored endpoints")
    ax.grid(alpha=.25)
    if ax.get_legend_handles_labels()[0]: ax.legend(fontsize=8)
    finish(fig, out/"01_KJ_first_passage_censor_aware.png")

    # 2. Calibrated cleavage K drive.
    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    for cls, g in df.groupby("class"):
        g = g.sort_values("target_psi_deg")
        ev = g[g["event_observed"]]
        ce = g[~g["event_observed"]]
        if len(ev):
            ax.plot(ev["target_psi_deg"], ev["Kcleave_calibrated_first_MPa_sqrt_m"],
                    marker="o", label=f"{cls}: event")
        if len(ce):
            ax.scatter(ce["target_psi_deg"], ce["Kcleave_calibrated_first_MPa_sqrt_m"],
                       marker="^", label=f"{cls}: censored")
    ax.set_xlabel("Target process-zone traction phase [deg]")
    ax.set_ylabel(r"Directional cleavage drive [MPa$\sqrt{\mathrm{m}}$]")
    ax.set_title("Calibrated sharp-tip cleavage drive")
    ax.grid(alpha=.25)
    if ax.get_legend_handles_labels()[0]: ax.legend(fontsize=8)
    finish(fig, out/"02_calibrated_cleavage_drive.png")

    # 3. Directional factors.
    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    for cls, g in df.groupby("class"):
        g = g.sort_values("target_psi_deg")
        ax.plot(g["target_psi_deg"], g["cleavage_factor_first"], marker="o",
                label=f"{cls}: cleavage")
        ax.plot(g["target_psi_deg"], g["emission_factor_first"], marker="s",
                linestyle="--", label=f"{cls}: emission")
    ax.axhline(1.0, linewidth=.8, linestyle=":")
    ax.set_xlabel("Target process-zone traction phase [deg]")
    ax.set_ylabel("Dimensionless directional multiplier")
    ax.set_title("Anisotropic directional partition factors")
    ax.grid(alpha=.25)
    if ax.get_legend_handles_labels()[0]: ax.legend(fontsize=8, ncol=2)
    finish(fig, out/"03_directional_factors.png")

    # 4. Event-state phase audit.
    fig, ax = plt.subplots(figsize=(7.0, 6.2))
    lim = 65
    ax.plot([-lim, lim], [-lim, lim], linestyle=":", linewidth=1.0,
            label="target = achieved")
    for cls, g in df.groupby("class"):
        ax.scatter(g["target_psi_deg"], g["traction_phase_first_deg"], label=cls)
    ax.set_xlim(-lim, lim); ax.set_ylim(-95, 95)
    ax.set_xlabel("Target traction phase [deg]")
    ax.set_ylabel("Achieved event/endpoint phase [deg]")
    ax.set_title("Event-state mode-control audit")
    ax.grid(alpha=.25)
    if ax.get_legend_handles_labels()[0]: ax.legend(fontsize=8)
    finish(fig, out/"04_event_state_phase_audit.png")

    # 5. Candidate crack direction.
    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    for cls, g in df.groupby("class"):
        g = g.sort_values("target_psi_deg")
        ax.plot(g["target_psi_deg"], g["candidate_angle_first_deg"], marker="o", label=cls)
    ax.set_xlabel("Target traction phase [deg]")
    ax.set_ylabel("Selected anisotropic crack direction [deg]")
    ax.set_title("Crystallographic direction competition")
    ax.grid(alpha=.25)
    if ax.get_legend_handles_labels()[0]: ax.legend(fontsize=8)
    finish(fig, out/"05_candidate_direction.png")

    # 6. Class-specific hazard state: ensures different parameterizations are visible.
    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    plotted = False
    for cls, g in df.groupby("class"):
        if "B_final" in g and np.any(np.isfinite(pd.to_numeric(g["B_final"], errors="coerce"))):
            g = g.sort_values("target_psi_deg")
            y = pd.to_numeric(g["B_final"], errors="coerce")
            ax.semilogy(g["target_psi_deg"], np.maximum(y, 1e-300), marker="o", label=cls)
            plotted = True
    if plotted:
        ax.set_xlabel("Target traction phase [deg]")
        ax.set_ylabel("Final cleavage clock $B$")
        ax.set_title("Class-specific cumulative cleavage clock")
        ax.grid(alpha=.25)
        ax.legend(fontsize=8)
        finish(fig, out/"06_class_specific_cleavage_clock.png")
    else:
        plt.close(fig)

    print("wrote", out)


if __name__ == "__main__":
    main()
