#!/usr/bin/env python3
"""Matched Figure 1 Panels C and D using one entropy function.

Shared entropy function
-----------------------

    S*(sigma,T)/kB = clip[-A_T F_T(T) - A_sigma F_sigma(sigma), S_min, 0]

The stress-dependent entropy coordinate shown in the figure is derived from the
same function:

    DeltaG_S,sigma,ref = -T0 [S*(sigma_ref,T0)-S*(0,T0)]
    Lambda_S,ref       = DeltaG_S,sigma,ref / H0

Lambda_S,ref is therefore not an independent shielding variable.

Panel C
-------
    x = log10 N_i
    y = Lambda_S,ref
    z = sigma_a [MPa]
    color = A_T [kB]

Panel D
-------
    x = T [K]
    y = Lambda_S,ref
    z = sigma_y [MPa]
    color = A_T [kB]

All user-facing stresses and all plotted stresses are in MPa.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from scipy.optimize import brentq

from shared_entropy_family_v3 import build_parent_barrier, multihit_rate


def parse_float_list(text: str) -> list[float]:
    vals = [float(x) for x in str(text).replace(",", " ").split()]
    if not vals:
        raise argparse.ArgumentTypeError("expected at least one number")
    return vals


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", type=Path, default=Path("runs/panels_CD_entropy_family_v3"))

    # Shared entropy family.
    p.add_argument("--A-T-kB", type=parse_float_list,
                   default=parse_float_list("0 1 2 3"),
                   help="temperature-entropy magnitudes; limited below 4 by default and used as the color coordinate")
    p.add_argument("--A-sigma-kB", type=parse_float_list,
                   default=parse_float_list("0 4 8 12 16 24 32 48"),
                   help="stress-entropy amplitudes used to derive Lambda_S,ref")
    p.add_argument("--T0-K", type=float, default=300.0)
    p.add_argument("--T-S-K", type=float, default=400.0)
    p.add_argument("--T-gate-power", type=float, default=4.0)
    p.add_argument("--sigma-S-MPa", type=parse_float_list,
                   default=parse_float_list("400"),
                   help="stress-gate scale(s) in MPa; keep one value for the main figure")
    p.add_argument("--sigma-gate-power", type=float, default=6.0)
    p.add_argument("--S-min-kB", type=float, default=-100.0)

    # Derived Lambda reference state.
    p.add_argument("--sigma-ref-MPa", type=float, default=600.0,
                   help="physical reference stress used in Lambda_S,ref")
    p.add_argument("--lambda-T-ref-K", type=float, default=None,
                   help="reference temperature for Lambda; default is T0")

    # Common H - T*S - sigma*v barrier, all stress scales in MPa.
    p.add_argument("--H0-eV", type=float, default=0.8)
    p.add_argument("--sigma0-H-MPa", type=float, default=2500.0)
    p.add_argument("--v0-b3", type=float, default=0.6)
    p.add_argument("--sigma0-v-MPa", type=float, default=2500.0)
    p.add_argument("--b-m", type=float, default=2.74e-10)
    p.add_argument("--nu0-s", type=float, default=1e11)

    # Panel C.
    p.add_argument("--C-T-K", type=float, default=300.0)
    p.add_argument("--C-stresses-MPa", type=parse_float_list,
                   default=parse_float_list("100 150 200 250 300 350 400 450 500 600 700 800 900 1000 1150 1300"))
    p.add_argument("--R", type=float, default=0.1)
    p.add_argument("--frequency-Hz", type=float, default=1000.0)
    p.add_argument("--Kt", type=float, default=3.0)
    p.add_argument("--cycles-max", type=float, default=1e12)
    p.add_argument("--n-phase", type=int, default=128)
    p.add_argument("--multihit-m", type=float, default=3.0)
    p.add_argument("--multihit-tau-s", type=float, default=1e-6)

    # Panel D.
    p.add_argument("--D-temperatures-K", type=parse_float_list,
                   default=parse_float_list("250 300 350 400 450 500 550 600 650 700 750 800 850 900 950 1000 1050 1100 1150 1200"))
    p.add_argument("--yield-strain-rate-s", type=float, default=1e-4)
    p.add_argument("--yield-event-strain", type=float, default=1e-5)
    p.add_argument("--yield-sigma-max-MPa", type=float, default=100000.0)

    p.add_argument("--resume", action="store_true")
    p.add_argument("--plot-only", action="store_true")
    p.add_argument("--figure-dpi", type=int, default=300)
    p.add_argument("--show-censored-markers", action="store_true",
                   help="show open-triangle right-censored markers in Panel C; hidden by default")
    p.add_argument("--view-elev", type=float, default=24.0)
    p.add_argument("--view-azim", type=float, default=-62.0)
    return p.parse_args()


def barrier(args, A_T: float, A_sigma: float, sigma_S_MPa: float):
    return build_parent_barrier(
        A_T_kB=A_T,
        A_sigma_kB=A_sigma,
        T0_K=args.T0_K,
        T_S_K=args.T_S_K,
        T_gate_power=args.T_gate_power,
        sigma_S_MPa=sigma_S_MPa,
        sigma_gate_power=args.sigma_gate_power,
        S_min_kB=args.S_min_kB,
        H0_eV=args.H0_eV,
        sigma0_H_MPa=args.sigma0_H_MPa,
        v0_b3=args.v0_b3,
        sigma0_v_MPa=args.sigma0_v_MPa,
        b_m=args.b_m,
    )


def entropy_family_metrics(args, b) -> dict:
    Tref = args.T0_K if args.lambda_T_ref_K is None else float(args.lambda_T_ref_K)
    dG = b.entropy.deltaG_stress_entropy_ref_eV(args.sigma_ref_MPa, Tref)
    lam = dG / max(float(args.H0_eV), 1e-30)
    gate = float(b.entropy.stress_gate(np.array([args.sigma_ref_MPa * 1e6]))[0])
    return {
        "sigma_ref_MPa": float(args.sigma_ref_MPa),
        "lambda_T_ref_K": float(Tref),
        "stress_entropy_gate_ref": gate,
        "deltaG_S_sigma_ref_eV": float(dG),
        "Lambda_S_ref": float(lam),
    }


def waveform_stress_Pa(sigma_a_MPa: float, R: float, n_phase: int) -> np.ndarray:
    """Tension-tension waveform parameterized by alternating stress amplitude."""
    amp = float(sigma_a_MPa) * 1e6
    if abs(1.0 - float(R)) < 1e-12:
        raise ValueError("R must be less than 1 when sigma_a is the alternating amplitude")
    mean = amp * (1.0 + float(R)) / (1.0 - float(R))
    n = max(int(n_phase), 16)
    phase = 2.0 * np.pi * (np.arange(n, dtype=float) + 0.5) / n
    return mean + amp * np.sin(phase)


def panel_C_point(args, A_T: float, A_sigma: float, sigma_S_MPa: float, sigma_a_MPa: float) -> dict:
    b = barrier(args, A_T, A_sigma, sigma_S_MPa)
    metrics = entropy_family_metrics(args, b)
    nominal = waveform_stress_Pa(sigma_a_MPa, args.R, args.n_phase)
    local = float(args.Kt) * nominal
    raw = b.rate_s(local, args.C_T_K, args.nu0_s)
    effective = multihit_rate(raw, args.multihit_m, args.multihit_tau_s)
    period = 1.0 / max(float(args.frequency_Hz), 1e-300)
    mu_cycle = float(np.mean(effective) * period)
    N_i = 1.0 / max(mu_cycle, 1e-300)
    status = "failed" if N_i <= args.cycles_max else "right_censored"
    sigma_peak_local = float(args.Kt) * (2.0 * float(sigma_a_MPa) * 1e6 / max(1.0 - float(args.R), 1e-12))
    S_peak = float(b.entropy.S_kB(np.array([sigma_peak_local]), args.C_T_K)[0])
    return {
        "panel": "C",
        "A_T_kB": float(A_T),
        "A_sigma_kB": float(A_sigma),
        "sigma_S_MPa": float(sigma_S_MPa),
        **metrics,
        "sigma_a_MPa": float(sigma_a_MPa),
        "T_K": float(args.C_T_K),
        "cycles_to_initiation": float(N_i) if status == "failed" else float("nan"),
        "cycles_raw_first_passage": float(N_i),
        "status": status,
        "mu_crack_per_cycle": mu_cycle,
        "S_peak_kB": S_peak,
    }


def yield_strength_MPa(args, b, T_K: float) -> tuple[float, str]:
    target_rate = float(args.yield_strain_rate_s) / max(float(args.yield_event_strain), 1e-300)
    smax = float(args.yield_sigma_max_MPa) * 1e6

    def f(sig):
        rate = float(b.rate_s(np.array([sig]), T_K, args.nu0_s)[0])
        return math.log(max(rate, 1e-300)) - math.log(max(target_rate, 1e-300))

    f0 = f(0.0)
    f1 = f(smax)
    if f0 >= 0.0:
        return 0.0, "zero_stress_active"
    if f1 < 0.0:
        return float("nan"), "above_search_ceiling"
    root = brentq(f, 0.0, smax, xtol=1.0, rtol=1e-10, maxiter=200)
    return root / 1e6, "resolved"


def panel_D_point(args, A_T: float, A_sigma: float, sigma_S_MPa: float, T_K: float) -> dict:
    b = barrier(args, A_T, A_sigma, sigma_S_MPa)
    metrics = entropy_family_metrics(args, b)
    sy, status = yield_strength_MPa(args, b, T_K)
    S_y = float("nan")
    Ft = float(b.entropy.temperature_gate(T_K))
    Fs = float("nan")
    if np.isfinite(sy):
        S_y = float(b.entropy.S_kB(np.array([sy * 1e6]), T_K)[0])
        Fs = float(b.entropy.stress_gate(np.array([sy * 1e6]))[0])
    return {
        "panel": "D",
        "A_T_kB": float(A_T),
        "A_sigma_kB": float(A_sigma),
        "sigma_S_MPa": float(sigma_S_MPa),
        **metrics,
        "T_K": float(T_K),
        "sigma_y_MPa": float(sy),
        "status": status,
        "S_y_kB": S_y,
        "temperature_entropy_gate": Ft,
        "stress_entropy_gate_at_yield": Fs,
    }


def append_row(path: Path, row: dict) -> None:
    pd.DataFrame([row]).to_csv(path, mode="a", header=not path.exists(), index=False)


def key_C(r) -> tuple[float, float, float, float]:
    return (
        round(float(r["A_T_kB"]), 10),
        round(float(r["A_sigma_kB"]), 10),
        round(float(r["sigma_S_MPa"]), 8),
        round(float(r["sigma_a_MPa"]), 8),
    )


def key_D(r) -> tuple[float, float, float, float]:
    return (
        round(float(r["A_T_kB"]), 10),
        round(float(r["A_sigma_kB"]), 10),
        round(float(r["sigma_S_MPa"]), 8),
        round(float(r["T_K"]), 8),
    )


def run_sweeps(args, outdir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    cpath = outdir / "panelC_entropy_SN_raw_v3.csv"
    dpath = outdir / "panelD_entropy_strength_raw_v3.csv"
    c0 = pd.read_csv(cpath) if args.resume and cpath.exists() else pd.DataFrame()
    d0 = pd.read_csv(dpath) if args.resume and dpath.exists() else pd.DataFrame()
    cdone = set(key_C(r) for _, r in c0.iterrows()) if not c0.empty else set()
    ddone = set(key_D(r) for _, r in d0.iterrows()) if not d0.empty else set()

    for A_T in args.A_T_kB:
        for A_sigma in args.A_sigma_kB:
            for sigma_S in args.sigma_S_MPa:
                for stress in args.C_stresses_MPa:
                    k = (round(A_T,10), round(A_sigma,10), round(sigma_S,8), round(stress,8))
                    if k not in cdone:
                        print(f"C A_T={A_T:g} A_sigma={A_sigma:g} sigmaS={sigma_S:g} MPa sigma_a={stress:g} MPa")
                        append_row(cpath, panel_C_point(args, A_T, A_sigma, sigma_S, stress))
                        cdone.add(k)
                for T in args.D_temperatures_K:
                    k = (round(A_T,10), round(A_sigma,10), round(sigma_S,8), round(T,8))
                    if k not in ddone:
                        append_row(dpath, panel_D_point(args, A_T, A_sigma, sigma_S, T))
                        ddone.add(k)
    return pd.read_csv(cpath), pd.read_csv(dpath)


def norm_for(df: pd.DataFrame, col: str) -> Normalize:
    vals = pd.to_numeric(df[col], errors="coerce").to_numpy(float)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return Normalize(vmin=0.0, vmax=1.0)
    vmin, vmax = float(np.min(vals)), float(np.max(vals))
    if math.isclose(vmin, vmax):
        vmax = vmin + 1.0
    return Normalize(vmin=vmin, vmax=vmax)


def plot_C(args, c: pd.DataFrame, outdir: Path) -> None:
    norm = norm_for(c, "A_T_kB")
    cmap = plt.get_cmap("viridis")
    fig = plt.figure(figsize=(10.5, 8.2))
    ax = fig.add_subplot(111, projection="3d")
    Nmax_log = math.log10(args.cycles_max)

    group_cols = ["A_T_kB", "A_sigma_kB", "sigma_S_MPa", "Lambda_S_ref"]
    for (A_T, A_sigma, sigmaS, lam), g in c.groupby(group_cols, sort=True):
        g = g.sort_values("sigma_a_MPa")
        color = cmap(norm(float(A_T)))
        failed = g.status.astype(str).to_numpy() == "failed"
        cens = g.status.astype(str).to_numpy() == "right_censored"
        rawN = pd.to_numeric(g.cycles_raw_first_passage, errors="coerce").to_numpy(float)
        x = np.where(failed, np.log10(np.maximum(rawN, 1e-300)), Nmax_log)
        y = np.full(len(g), float(lam))
        z = g.sigma_a_MPa.to_numpy(float)
        ax.plot(x, y, z, color=color, lw=1.55, alpha=0.92)
        if failed.any():
            idx = np.where(failed)[0]
            ax.scatter(x[idx], y[idx], z[idx], color=[color], s=13)
        if args.show_censored_markers and cens.any():
            idx = np.where(cens)[0]
            ax.scatter(x[idx], y[idx], z[idx], marker="^", s=30,
                       facecolors="none", edgecolors=[color], linewidths=0.95)

    ax.set_xlabel(r"$\log_{10}$ cycles to crack initiation, $N_i$", labelpad=10)
    ax.set_ylabel(r"Entropic barrier ratio $\Lambda_{S,\mathrm{ref}}$", labelpad=10)
    ax.set_zlabel(r"Stress amplitude $\sigma_a$ [MPa]", labelpad=11)
    ax.view_init(elev=args.view_elev, azim=args.view_azim)
    failed_all = c.status.astype(str).to_numpy() == "failed"
    rawN = pd.to_numeric(c.cycles_raw_first_passage, errors="coerce").to_numpy(float)
    allx = np.where(failed_all, np.log10(np.maximum(rawN, 1e-300)), Nmax_log)
    if np.isfinite(allx).any():
        ax.set_xlim(max(0.0, float(np.nanmin(allx))-0.2), Nmax_log+0.12)
    sm = ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, fraction=0.045, pad=0.10)
    cb.set_label(r"Temperature-entropy magnitude $A_T=-S_T(T_0)/k_B$")
    fig.tight_layout(rect=(0,0.03,0.94,1))
    fig.savefig(outdir/"panelC_entropy_SN_LambdaS_waterfall_3d.png", dpi=args.figure_dpi, bbox_inches="tight")
    fig.savefig(outdir/"panelC_entropy_SN_LambdaS_waterfall_3d.pdf", bbox_inches="tight")
    plt.close(fig)

    # 2-D diagnostic grouped by A_T to verify endurance-knee development.
    ATs = sorted(pd.to_numeric(c.A_T_kB, errors="coerce").dropna().unique())
    fig, axes = plt.subplots(1, len(ATs), figsize=(4.6*len(ATs), 3.8), squeeze=False)
    lam_norm = norm_for(c, "Lambda_S_ref")
    for ax2, AT in zip(axes[0], ATs):
        q = c[np.isclose(c.A_T_kB, AT)]
        for lam, g in q.groupby("Lambda_S_ref", sort=True):
            g = g.sort_values("sigma_a_MPa")
            failed = g.status.astype(str).to_numpy() == "failed"
            Nplot = np.where(failed,
                             pd.to_numeric(g.cycles_raw_first_passage, errors="coerce").to_numpy(float),
                             args.cycles_max)
            ax2.plot(Nplot, g.sigma_a_MPa, lw=1.35, color=plt.get_cmap("plasma")(lam_norm(float(lam))))
        ax2.set_xscale("log")
        ax2.set_title(rf"$A_T={AT:g}$")
        ax2.set_xlabel(r"Cycles $N_i$")
        ax2.set_ylabel(r"$\sigma_a$ [MPa]")
        ax2.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(outdir/"panelC_SN_endurance_diagnostic.png", dpi=args.figure_dpi, bbox_inches="tight")
    plt.close(fig)


def plot_D(args, d: pd.DataFrame, outdir: Path) -> None:
    good = d[np.isfinite(pd.to_numeric(d.sigma_y_MPa, errors="coerce").to_numpy(float))].copy()
    norm = norm_for(good, "A_T_kB")
    cmap = plt.get_cmap("viridis")
    fig = plt.figure(figsize=(10.5, 8.2))
    ax = fig.add_subplot(111, projection="3d")

    group_cols = ["A_T_kB", "A_sigma_kB", "sigma_S_MPa", "Lambda_S_ref"]
    for (A_T, A_sigma, sigmaS, lam), g in good.groupby(group_cols, sort=True):
        g = g.sort_values("T_K")
        color = cmap(norm(float(A_T)))
        ax.plot(g.T_K, np.full(len(g), float(lam)), g.sigma_y_MPa,
                color=color, lw=1.65, alpha=0.92)
        ax.scatter(g.T_K, np.full(len(g), float(lam)), g.sigma_y_MPa,
                   color=[color], s=11)

    ax.set_xlabel("Temperature [K]", labelpad=10)
    ax.set_ylabel(r"Entropic barrier ratio $\Lambda_{S,\mathrm{ref}}$", labelpad=10)
    ax.set_zlabel(r"Yield strength $\sigma_y$ [MPa]", labelpad=11)
    ax.view_init(elev=args.view_elev, azim=args.view_azim)
    sm = ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, fraction=0.045, pad=0.10)
    cb.set_label(r"Temperature-entropy magnitude $A_T=-S_T(T_0)/k_B$")
    fig.tight_layout(rect=(0,0.03,0.94,1))
    fig.savefig(outdir/"panelD_entropy_strength_LambdaS_waterfall_3d.png", dpi=args.figure_dpi, bbox_inches="tight")
    fig.savefig(outdir/"panelD_entropy_strength_LambdaS_waterfall_3d.pdf", bbox_inches="tight")
    plt.close(fig)

    # 2-D diagnostic, grouped by the small A_T set and colored by Lambda_S_ref.
    ATs = sorted(pd.to_numeric(good.A_T_kB, errors="coerce").dropna().unique())
    lam_norm = norm_for(good, "Lambda_S_ref")
    lam_cmap = plt.get_cmap("plasma")
    fig, axes = plt.subplots(1, len(ATs), figsize=(4.6*len(ATs), 3.8), squeeze=False)
    for ax2, AT in zip(axes[0], ATs):
        q = good[np.isclose(good.A_T_kB, AT)]
        for lam, g in q.groupby("Lambda_S_ref", sort=True):
            g = g.sort_values("T_K")
            ax2.plot(g.T_K, g.sigma_y_MPa, lw=1.4, color=lam_cmap(lam_norm(float(lam))))
        ax2.set_title(rf"$A_T={AT:g}$")
        ax2.set_xlabel("T [K]")
        ax2.set_ylabel(r"$\sigma_y$ [MPa]")
        ax2.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(outdir/"panelD_strength_temperature_diagnostic_v3.png", dpi=args.figure_dpi, bbox_inches="tight")
    plt.close(fig)


def curve_summaries(args, c: pd.DataFrame, d: pd.DataFrame, outdir: Path) -> None:
    crows = []
    for keys, g in c.groupby(["A_T_kB","A_sigma_kB","sigma_S_MPa","Lambda_S_ref"], sort=True):
        AT, Asig, sS, lam = keys
        f = g[g.status.astype(str)=="failed"].sort_values("cycles_raw_first_passage")
        slope = float("nan")
        if len(f) >= 3:
            x = np.log10(f.cycles_raw_first_passage.to_numpy(float))
            y = f.sigma_a_MPa.to_numpy(float)
            use = np.arange(len(x)) >= max(0, len(x)//2 - 1)
            if use.sum() >= 2 and np.ptp(x[use]) > 0:
                slope = float(abs(np.polyfit(x[use], y[use], 1)[0]))
        crows.append({
            "A_T_kB":AT,"A_sigma_kB":Asig,"sigma_S_MPa":sS,"Lambda_S_ref":lam,
            "n_failed":int(len(f)),"n_censored":int((g.status.astype(str)=="right_censored").sum()),
            "high_cycle_abs_slope_MPa_per_decade":slope,
        })
    pd.DataFrame(crows).to_csv(outdir/"panelC_curve_summary_v3.csv", index=False)

    drows=[]
    for keys,g in d.groupby(["A_T_kB","A_sigma_kB","sigma_S_MPa","Lambda_S_ref"], sort=True):
        AT,Asig,sS,lam=keys
        q=g[np.isfinite(pd.to_numeric(g.sigma_y_MPa, errors="coerce"))].sort_values("T_K")
        if q.empty:
            continue
        sy=q.sigma_y_MPa.to_numpy(float); T=q.T_K.to_numpy(float)
        i=int(np.argmax(sy))
        drows.append({
            "A_T_kB":AT,"A_sigma_kB":Asig,"sigma_S_MPa":sS,"Lambda_S_ref":lam,
            "sigma_y_300_MPa":float(np.interp(300.0,T,sy)) if T.min()<=300<=T.max() else np.nan,
            "T_peak_K":float(T[i]),"sigma_peak_MPa":float(sy[i]),
            "peak_amp_frac_vs_300":float(sy[i]/max(np.interp(300.0,T,sy),1e-30)-1.0) if T.min()<=300<=T.max() else np.nan,
        })
    pd.DataFrame(drows).to_csv(outdir/"panelD_curve_summary_v3.csv", index=False)


def write_metadata(args, outdir: Path) -> None:
    rows = []
    for AT in args.A_T_kB:
        for Asig in args.A_sigma_kB:
            for sigmaS in args.sigma_S_MPa:
                b = barrier(args, AT, Asig, sigmaS)
                m = entropy_family_metrics(args, b)
                rows.append({
                    "A_T_kB": AT,
                    "A_sigma_kB": Asig,
                    "sigma_S_MPa": sigmaS,
                    **m,
                    "T0_K": args.T0_K,
                    "T_S_K": args.T_S_K,
                    "T_gate_power": args.T_gate_power,
                    "sigma_gate_power": args.sigma_gate_power,
                    "S_min_kB": args.S_min_kB,
                    "H0_eV": args.H0_eV,
                    "sigma0_H_MPa": args.sigma0_H_MPa,
                    "v0_b3": args.v0_b3,
                    "sigma0_v_MPa": args.sigma0_v_MPa,
                })
    pd.DataFrame(rows).to_csv(outdir/"panels_CD_entropy_family_design_v3.csv", index=False)
    manifest = {
        "shared_entropy_function": "S*/kB = clip[-A_T F_T(T) - A_sigma F_sigma(sigma), S_min, 0]",
        "temperature_gate": "F_T(T)=hill(T/T_S,p_T)/hill(T0/T_S,p_T)",
        "stress_gate": "F_sigma=x^p_sigma/(1+x^p_sigma), x=|sigma|/sigma_S",
        "derived_coordinate": "Lambda_S_ref = {-T_ref [S(sigma_ref,T_ref)-S(0,T_ref)]}/H0",
        "shared_barrier": "G*=H(sigma)-T*S*(sigma,T)-sigma*v(sigma)",
        "panel_C_axes": {"x":"log10 N_i","y":"Lambda_S_ref","z":"sigma_a_MPa","color":"A_T_kB"},
        "panel_D_axes": {"x":"T_K","y":"Lambda_S_ref","z":"sigma_y_MPa","color":"A_T_kB"},
        "stress_units": "MPa for all user-facing inputs, CSV stress fields, and plotted stress axes",
        "sigma_ref_MPa": args.sigma_ref_MPa,
        "lambda_T_ref_K": args.T0_K if args.lambda_T_ref_K is None else args.lambda_T_ref_K,
        "note": "Lambda_S_ref is derived from the shared entropy function. No independent chi_shield, Gshield, or Lambda_sh parameter is swept."
    }
    (outdir/"panels_CD_manifest_v3.json").write_text(json.dumps(manifest, indent=2))


def main():
    args = parse_args()
    outdir = args.out
    outdir.mkdir(parents=True, exist_ok=True)
    cpath = outdir/"panelC_entropy_SN_raw_v3.csv"
    dpath = outdir/"panelD_entropy_strength_raw_v3.csv"
    if args.plot_only:
        if not cpath.exists() or not dpath.exists():
            raise FileNotFoundError("--plot-only requires existing Panel C and D v2 raw CSV files")
        c = pd.read_csv(cpath); d = pd.read_csv(dpath)
    else:
        c, d = run_sweeps(args, outdir)
    write_metadata(args, outdir)
    curve_summaries(args, c, d, outdir)
    plot_C(args, c, outdir)
    plot_D(args, d, outdir)
    print("\nOutputs written to", outdir)
    for p in sorted(outdir.iterdir()):
        print(" ", p.name)


if __name__ == "__main__":
    main()
