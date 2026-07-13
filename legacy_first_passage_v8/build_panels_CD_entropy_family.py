#!/usr/bin/env python3
"""Calculate matched Figure 1 Panels C and D from one activation-entropy family.

Shared entropy function
-----------------------

    S*(sigma,T)/kB = clip[-A_T F_T(T) - A_sigma F_sigma(sigma), S_min, 0]

A_sigma is the common third axis in Panels C and D.  A_T is the common color
coordinate.  No chi_shield, Gshield, or Lambda_sh parameter is swept.

Panel C
-------
A cyclic first-passage calculation integrates the common Arrhenius barrier over
one stress cycle.  The cycle hazard is mu, and the deterministic first-passage
life is N_i = 1/mu.  An optional cooperative multihit renewal transform can be
applied to the instantaneous raw rate before cycle integration.

Panel D
-------
The same barrier and same entropy function are inverted at fixed event rate to
obtain sigma_y(T).
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

from shared_entropy_family import build_parent_barrier, multihit_rate


def parse_float_list(text: str) -> list[float]:
    vals = [float(x) for x in str(text).replace(",", " ").split()]
    if not vals:
        raise argparse.ArgumentTypeError("expected at least one number")
    return vals


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", type=Path, default=Path("runs/panels_CD_entropy_family"))

    # Shared entropy family.
    p.add_argument("--A-T-kB", type=parse_float_list, default=parse_float_list("0 2 4"),
                   help="temperature-entropy magnitude at T0; common color coordinate")
    p.add_argument("--A-sigma-kB", type=parse_float_list, default=parse_float_list("0 2 4 6 8"),
                   help="stress-entropy amplitude; common waterfall axis")
    p.add_argument("--T0-K", type=float, default=300.0)
    p.add_argument("--T-S-K", type=float, default=400.0)
    p.add_argument("--T-gate-power", type=float, default=4.0)
    p.add_argument("--sigma-S-GPa", type=float, default=3.0)
    p.add_argument("--sigma-gate-power", type=float, default=1.0)
    p.add_argument("--S-min-kB", type=float, default=-40.0)

    # Common H - T*S - sigma*v barrier.
    p.add_argument("--H0-eV", type=float, default=0.8)
    p.add_argument("--sigma0-H-GPa", type=float, default=2.5)
    p.add_argument("--v0-b3", type=float, default=0.6)
    p.add_argument("--sigma0-v-GPa", type=float, default=2.5)
    p.add_argument("--b-m", type=float, default=2.74e-10)
    p.add_argument("--nu0-s", type=float, default=1e11)

    # Panel C.
    p.add_argument("--C-T-K", type=float, default=300.0)
    p.add_argument("--C-stresses-MPa", type=parse_float_list,
                   default=parse_float_list("100 150 200 250 300 350 400 450 500 600 700 800"))
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
    p.add_argument("--yield-sigma-max-GPa", type=float, default=100.0)

    p.add_argument("--resume", action="store_true")
    p.add_argument("--plot-only", action="store_true")
    p.add_argument("--figure-dpi", type=int, default=300)
    p.add_argument("--view-elev", type=float, default=24.0)
    p.add_argument("--view-azim", type=float, default=-62.0)
    return p.parse_args()


def barrier(args, A_T: float, A_sigma: float):
    return build_parent_barrier(
        A_T_kB=A_T,
        A_sigma_kB=A_sigma,
        T0_K=args.T0_K,
        T_S_K=args.T_S_K,
        T_gate_power=args.T_gate_power,
        sigma_S_GPa=args.sigma_S_GPa,
        sigma_gate_power=args.sigma_gate_power,
        S_min_kB=args.S_min_kB,
        H0_eV=args.H0_eV,
        sigma0_H_GPa=args.sigma0_H_GPa,
        v0_b3=args.v0_b3,
        sigma0_v_GPa=args.sigma0_v_GPa,
        b_m=args.b_m,
    )


def waveform_stress_Pa(sigma_a_MPa: float, R: float, n_phase: int) -> np.ndarray:
    """Tension-tension waveform parameterized by alternating stress amplitude.

    With R=sigma_min/sigma_max, sigma_m/sigma_a=(1+R)/(1-R).
    """
    amp = float(sigma_a_MPa) * 1e6
    if abs(1.0 - float(R)) < 1e-12:
        raise ValueError("R must be less than 1 when sigma_a is the alternating amplitude")
    mean = amp * (1.0 + float(R)) / (1.0 - float(R))
    phase = 2.0 * np.pi * (np.arange(max(int(n_phase), 16), dtype=float) + 0.5) / max(int(n_phase), 16)
    return mean + amp * np.sin(phase)


def panel_C_point(args, A_T: float, A_sigma: float, sigma_a_MPa: float) -> dict:
    b = barrier(args, A_T, A_sigma)
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
    smax = float(args.yield_sigma_max_GPa) * 1e9

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


def panel_D_point(args, A_T: float, A_sigma: float, T_K: float) -> dict:
    b = barrier(args, A_T, A_sigma)
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
        "T_K": float(T_K),
        "sigma_y_MPa": float(sy),
        "status": status,
        "S_y_kB": S_y,
        "temperature_entropy_gate": Ft,
        "stress_entropy_gate_at_yield": Fs,
    }


def append_row(path: Path, row: dict) -> None:
    pd.DataFrame([row]).to_csv(path, mode="a", header=not path.exists(), index=False)


def key_C(r) -> tuple[float, float, float]:
    return (round(float(r["A_T_kB"]), 10), round(float(r["A_sigma_kB"]), 10), round(float(r["sigma_a_MPa"]), 8))


def key_D(r) -> tuple[float, float, float]:
    return (round(float(r["A_T_kB"]), 10), round(float(r["A_sigma_kB"]), 10), round(float(r["T_K"]), 8))


def run_sweeps(args, outdir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    cpath = outdir / "panelC_entropy_SN_raw.csv"
    dpath = outdir / "panelD_entropy_strength_raw.csv"
    c0 = pd.read_csv(cpath) if args.resume and cpath.exists() else pd.DataFrame()
    d0 = pd.read_csv(dpath) if args.resume and dpath.exists() else pd.DataFrame()
    cdone = set(key_C(r) for _, r in c0.iterrows()) if not c0.empty else set()
    ddone = set(key_D(r) for _, r in d0.iterrows()) if not d0.empty else set()

    for A_T in args.A_T_kB:
        for A_sigma in args.A_sigma_kB:
            for stress in args.C_stresses_MPa:
                k = (round(A_T,10), round(A_sigma,10), round(stress,8))
                if k not in cdone:
                    print(f"C A_T={A_T:g} A_sigma={A_sigma:g} sigma={stress:g} MPa")
                    append_row(cpath, panel_C_point(args, A_T, A_sigma, stress))
                    cdone.add(k)
            for T in args.D_temperatures_K:
                k = (round(A_T,10), round(A_sigma,10), round(T,8))
                if k not in ddone:
                    append_row(dpath, panel_D_point(args, A_T, A_sigma, T))
                    ddone.add(k)
    return pd.read_csv(cpath), pd.read_csv(dpath)


def norm_cmap(df: pd.DataFrame):
    vals = df["A_T_kB"].to_numpy(float)
    vmin, vmax = float(np.nanmin(vals)), float(np.nanmax(vals))
    if math.isclose(vmin, vmax):
        vmax = vmin + 1.0
    return Normalize(vmin=vmin, vmax=vmax), plt.get_cmap("viridis")


def plot_C(args, c: pd.DataFrame, outdir: Path) -> None:
    norm, cmap = norm_cmap(c)
    fig = plt.figure(figsize=(10.5, 8.2))
    ax = fig.add_subplot(111, projection="3d")
    Nmax_log = math.log10(args.cycles_max)
    for (A_T, A_sigma), g in c.groupby(["A_T_kB", "A_sigma_kB"], sort=True):
        g = g.sort_values("sigma_a_MPa")
        color = cmap(norm(float(A_T)))
        failed = g.status.astype(str).to_numpy() == "failed"
        cens = g.status.astype(str).to_numpy() == "right_censored"
        x = np.where(failed, np.log10(g.cycles_raw_first_passage.to_numpy(float)), Nmax_log)
        y = np.full(len(g), float(A_sigma))
        z = g.sigma_a_MPa.to_numpy(float)
        # Connect the displayed waterfall curve through both resolved and censored endpoints.
        ax.plot(x, y, z, color=color, lw=1.55, alpha=0.92)
        if failed.any():
            idx = np.where(failed)[0]
            ax.scatter(x[idx], y[idx], z[idx], color=[color], s=13)
        if cens.any():
            idx = np.where(cens)[0]
            ax.scatter(x[idx], y[idx], z[idx], marker="^", s=30,
                       facecolors="none", edgecolors=[color], linewidths=0.95)
    ax.set_xlabel(r"$\log_{10}$ cycles to crack initiation, $N_i$", labelpad=10)
    ax.set_ylabel(r"Stress-entropy amplitude $A_\sigma$ [$k_B$]", labelpad=10)
    ax.set_zlabel(r"Stress amplitude $\sigma_a$ [MPa]", labelpad=11)
    ax.view_init(elev=args.view_elev, azim=args.view_azim)
    allx = np.where(c.status.astype(str)=="failed", np.log10(c.cycles_raw_first_passage), Nmax_log)
    ax.set_xlim(max(0.0, float(np.nanmin(allx))-0.2), Nmax_log+0.12)
    sm = ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, fraction=0.045, pad=0.10)
    cb.set_label(r"Temperature-entropy magnitude $A_T=-S_T(T_0)/k_B$")
    fig.tight_layout(rect=(0,0.03,0.94,1))
    fig.savefig(outdir/"panelC_entropy_SN_waterfall_3d.png", dpi=args.figure_dpi, bbox_inches="tight")
    fig.savefig(outdir/"panelC_entropy_SN_waterfall_3d.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_D(args, d: pd.DataFrame, outdir: Path) -> None:
    good = d[np.isfinite(d.sigma_y_MPa.to_numpy(float))].copy()
    norm, cmap = norm_cmap(good)
    fig = plt.figure(figsize=(10.5, 8.2))
    ax = fig.add_subplot(111, projection="3d")
    for (A_T, A_sigma), g in good.groupby(["A_T_kB", "A_sigma_kB"], sort=True):
        g = g.sort_values("T_K")
        color = cmap(norm(float(A_T)))
        ax.plot(g.T_K, np.full(len(g), float(A_sigma)), g.sigma_y_MPa,
                color=color, lw=1.65, alpha=0.92)
        ax.scatter(g.T_K, np.full(len(g), float(A_sigma)), g.sigma_y_MPa,
                   color=[color], s=11)
    ax.set_xlabel("Temperature [K]", labelpad=10)
    ax.set_ylabel(r"Stress-entropy amplitude $A_\sigma$ [$k_B$]", labelpad=10)
    ax.set_zlabel(r"Yield strength $\sigma_y$ [MPa]", labelpad=11)
    ax.view_init(elev=args.view_elev, azim=args.view_azim)
    sm = ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, fraction=0.045, pad=0.10)
    cb.set_label(r"Temperature-entropy magnitude $A_T=-S_T(T_0)/k_B$")
    fig.tight_layout(rect=(0,0.03,0.94,1))
    fig.savefig(outdir/"panelD_entropy_strength_waterfall_3d.png", dpi=args.figure_dpi, bbox_inches="tight")
    fig.savefig(outdir/"panelD_entropy_strength_waterfall_3d.pdf", bbox_inches="tight")
    plt.close(fig)

    # Diagnostic 2-D strength-temperature curves for topology checking.
    ATs = sorted(good.A_T_kB.unique())
    fig, axes = plt.subplots(1, len(ATs), figsize=(4.5*len(ATs), 3.8), squeeze=False)
    for ax, AT in zip(axes[0], ATs):
        q = good[np.isclose(good.A_T_kB, AT)]
        for Asig, g in q.groupby("A_sigma_kB"):
            g = g.sort_values("T_K")
            ax.plot(g.T_K, g.sigma_y_MPa, lw=1.4, label=rf"$A_\sigma={Asig:g}$")
        ax.set_title(rf"$A_T={AT:g}$")
        ax.set_xlabel("T [K]")
        ax.set_ylabel(r"$\sigma_y$ [MPa]")
        ax.grid(alpha=0.25)
    axes[0,-1].legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(outdir/"panelD_strength_temperature_diagnostic.png", dpi=args.figure_dpi, bbox_inches="tight")
    plt.close(fig)


def write_metadata(args, outdir: Path) -> None:
    rows = []
    for AT in args.A_T_kB:
        for Asig in args.A_sigma_kB:
            rows.append({
                "A_T_kB": AT,
                "A_sigma_kB": Asig,
                "T0_K": args.T0_K,
                "T_S_K": args.T_S_K,
                "T_gate_power": args.T_gate_power,
                "sigma_S_GPa": args.sigma_S_GPa,
                "sigma_gate_power": args.sigma_gate_power,
                "S_min_kB": args.S_min_kB,
                "H0_eV": args.H0_eV,
                "sigma0_H_GPa": args.sigma0_H_GPa,
                "v0_b3": args.v0_b3,
                "sigma0_v_GPa": args.sigma0_v_GPa,
            })
    pd.DataFrame(rows).to_csv(outdir/"panels_CD_entropy_family_design.csv", index=False)
    manifest = {
        "shared_entropy_function": "S*/kB = clip[-A_T F_T(T) - A_sigma F_sigma(sigma), S_min, 0]",
        "temperature_gate": "F_T(T)=hill(T/T_S,p_T)/hill(T0/T_S,p_T)",
        "stress_gate": "F_sigma=x^p_sigma/(1+x^p_sigma), x=|sigma|/sigma_S",
        "shared_barrier": "G*=H(sigma)-T*S*(sigma,T)-sigma*v(sigma)",
        "panel_C": "cycle-integrated first-passage hazard; N_i=1/mu_cycle; optional cooperative multihit transform",
        "panel_D": "fixed-rate inversion of the same raw Arrhenius barrier",
        "panel_C_axes": {"x":"log10 N_i","y":"A_sigma_kB","z":"sigma_a_MPa","color":"A_T_kB"},
        "panel_D_axes": {"x":"T_K","y":"A_sigma_kB","z":"sigma_y_MPa","color":"A_T_kB"},
        "no_shielding_sweep": True,
        "note": "No chi_shield, Gshield, or Lambda_sh parameter is swept. Panels C and D share the exact same entropy helper and parent barrier parameters."
    }
    (outdir/"panels_CD_manifest.json").write_text(json.dumps(manifest, indent=2))


def main():
    args = parse_args()
    outdir = args.out
    outdir.mkdir(parents=True, exist_ok=True)
    cpath = outdir/"panelC_entropy_SN_raw.csv"
    dpath = outdir/"panelD_entropy_strength_raw.csv"
    if args.plot_only:
        if not cpath.exists() or not dpath.exists():
            raise FileNotFoundError("--plot-only requires existing Panel C and D raw CSV files")
        c = pd.read_csv(cpath); d = pd.read_csv(dpath)
    else:
        c, d = run_sweeps(args, outdir)
    write_metadata(args, outdir)
    plot_C(args, c, outdir)
    plot_D(args, d, outdir)
    print("\nOutputs written to", outdir)
    for p in sorted(outdir.iterdir()):
        print(" ", p.name)


if __name__ == "__main__":
    main()
