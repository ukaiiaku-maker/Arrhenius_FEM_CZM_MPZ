#!/usr/bin/env python3
"""Build Panel C: a systematic 3-D S-N crack-initiation waterfall from the reduced V1 model.

Plot geometry
-------------
    x = log10(cycles to crack initiation, N_i)
    y = -S_crack*/k_B  [dimensionless entropy magnitude]
    z = nominal stress amplitude sigma_a [MPa]

Line-family parameter
---------------------
    Lambda_sh = G_shield / G0_emit(T)

Why Lambda_sh instead of an artificial G0,c/G0,e knob?
-------------------------------------------------------
The current reduced V1 S-N implementation exposes the crack-opening entropy term
S_crack and the shielding barrier increment Gshield as independent native controls,
but does not expose a separate independent zero-stress crack-nucleation barrier
G0,c.  This driver therefore uses a dimensionless, model-native competition ratio,
Lambda_sh, rather than inventing a crack-barrier parameter that the solver would
not actually use.

The script is restartable. Each completed point is appended to the raw CSV and is
skipped on a subsequent run with --resume. Right-censored points are retained and
shown as open triangles at the cycle horizon.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from types import SimpleNamespace
from typing import Iterable

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401


# W[100] EXP-family fallback values used by the current project when no explicit
# custom EXP-floor override is present. These are used ONLY to normalize the
# dimensionless Lambda_sh family axis; the solver itself remains the source of
# all event rates and state evolution.
W100_G00_eV = 1.94
W100_gT_eV_per_K = 3.934e-3
W100_Tref_K = 481.33


def parse_float_list(text: str) -> list[float]:
    vals = [float(x) for x in text.replace(",", " ").split()]
    if not vals:
        raise argparse.ArgumentTypeError("expected at least one number")
    return vals


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--project-root", type=Path, default=Path.cwd())
    p.add_argument("--out", type=Path, default=Path("runs/panelC_SN_nucleation_waterfall"))

    p.add_argument("--entropy-mag-kB", type=parse_float_list,
                   default=parse_float_list("0 4 8 12 16 20 24"),
                   help="Quoted list for y=-S_crack*/kB; solver receives S_crack_kB=-y.")
    p.add_argument("--lambda-sh", type=parse_float_list,
                   default=parse_float_list("0 0.15 0.30 0.45 0.60"),
                   help="Quoted Lambda_sh=Gshield/G0_emit(T) line-family values.")
    p.add_argument("--stresses-MPa", type=parse_float_list,
                   default=parse_float_list(
                       "150 200 250 300 350 400 450 500 550 600 700 800 900 1050 1200 1400"
                   ))

    p.add_argument("--T-K", type=float, default=300.0)
    p.add_argument("--R", type=float, default=0.1)
    p.add_argument("--frequency-Hz", type=float, default=1000.0)
    p.add_argument("--chi-back", type=float, default=0.60,
                   help="Fixed back-stress coupling for all families.")

    p.add_argument("--cycles-max", type=float, default=1e12)
    p.add_argument("--max-blocks", type=int, default=20000)
    p.add_argument("--block-cycles", type=float, default=1e8)
    p.add_argument("--n-phase", type=int, default=64)
    p.add_argument("--target-dB-nuc", type=float, default=0.05)
    p.add_argument("--target-dep-eq-block", type=float, default=2e-4)
    p.add_argument("--target-rho-rel-block", type=float, default=0.05)

    p.add_argument("--emit-energy-scale", type=float, default=0.75)
    p.add_argument("--emit-entropy-scale", type=float, default=0.75)
    p.add_argument("--peierls-energy-scale", type=float, default=0.00375)
    p.add_argument("--peierls-entropy-scale", type=float, default=0.00375)
    p.add_argument("--taylor-energy-scale", type=float, default=0.015)
    p.add_argument("--taylor-entropy-scale", type=float, default=0.015)

    p.add_argument("--resume", action="store_true")
    p.add_argument("--plot-only", action="store_true")
    p.add_argument("--print-every-point", action="store_true")

    p.add_argument("--view-elev", type=float, default=24.0)
    p.add_argument("--view-azim", type=float, default=-62.0)
    p.add_argument("--figure-dpi", type=int, default=300)
    return p.parse_args()


def import_v1(project_root: Path):
    root = project_root.resolve()
    sys.path.insert(0, str(root))
    try:
        from arrhenius_fracture.sn_v1_arrhenius import (  # type: ignore
            build_parser as sn_parser,
            run_point,
            SNCase,
        )
    except Exception as exc:
        raise RuntimeError(
            "Could not import arrhenius_fracture.sn_v1_arrhenius from "
            f"{root}. Activate the fatigue-pf conda environment and run from the "
            "Fatigue-PF project root. Original error: " + repr(exc)
        ) from exc
    return sn_parser, run_point, SNCase


def make_base_args(sn_parser, args: argparse.Namespace):
    a = sn_parser().parse_args([])
    a.T = float(args.T_K)
    a.R = float(args.R)
    a.frequency_Hz = float(args.frequency_Hz)
    a.cycles_max = float(args.cycles_max)
    a.max_blocks = int(args.max_blocks)
    a.block_cycles = float(args.block_cycles)
    a.n_phase = int(args.n_phase)
    a.target_dB_nuc = float(args.target_dB_nuc)
    a.target_dep_eq_block = float(args.target_dep_eq_block)
    a.target_rho_rel_block = float(args.target_rho_rel_block)

    a.emit_energy_scale = float(args.emit_energy_scale)
    a.emit_entropy_scale = float(args.emit_entropy_scale)
    a.peierls_energy_scale = float(args.peierls_energy_scale)
    a.peierls_entropy_scale = float(args.peierls_entropy_scale)
    a.taylor_energy_scale = float(args.taylor_energy_scale)
    a.taylor_entropy_scale = float(args.taylor_entropy_scale)
    return a


def emission_zero_stress_G0_eV(a, T_K: float) -> float:
    """Reference zero-stress emission free energy used only to normalize Lambda_sh.

    Custom EXP-floor fields take precedence. Otherwise use the project's current
    W[100] family fallback. Energy and entropy scales are applied separately,
    matching the current chain parameterization.
    """
    # The V1 parser defines the custom EXP-floor override fields with default
    # None.  getattr(..., fallback) does not use the fallback when the
    # attribute exists but is None, so coalesce None explicitly.
    def _value_or_default(name: str, default: float) -> float:
        value = getattr(a, name, None)
        return float(default if value is None else value)

    G00 = _value_or_default("exp_G00_eV", W100_G00_eV)
    gT = _value_or_default("exp_gT_eV_per_K", W100_gT_eV_per_K)
    Tref = _value_or_default("exp_Tref_K", W100_Tref_K)
    e_scale = _value_or_default("emit_energy_scale", 1.0)
    s_scale = _value_or_default("emit_entropy_scale", e_scale)
    G = e_scale * G00 + s_scale * gT * (float(T_K) - Tref)
    return max(float(G), 1e-9)


def family_design(base_args, args: argparse.Namespace) -> pd.DataFrame:
    Gemit0 = emission_zero_stress_G0_eV(base_args, args.T_K)
    rows = []
    for i, lam in enumerate(sorted(set(float(x) for x in args.lambda_sh)), start=1):
        rows.append({
            "family_id": f"LAMBDA_{i:02d}",
            "family_index": i,
            "Lambda_sh": lam,
            "G0_emit_ref_eV": Gemit0,
            "Gshield_eV": lam * Gemit0,
            "chi_back": float(args.chi_back),
        })
    return pd.DataFrame(rows)


def read_existing_raw(path: Path) -> pd.DataFrame:
    if not path.is_file():
        return pd.DataFrame()
    df = pd.read_csv(path)
    return df


def point_key(row) -> tuple[float, float, float]:
    return (
        round(float(row["Lambda_sh"]), 12),
        round(float(row["entropy_mag_kB"]), 12),
        round(float(row["sigma_a_MPa"]), 9),
    )


def append_csv_row(path: Path, row: dict) -> None:
    frame = pd.DataFrame([row])
    frame.to_csv(path, mode="a", header=not path.exists(), index=False)


def run_sweep(args: argparse.Namespace, outdir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    sn_parser, run_point, SNCase = import_v1(args.project_root)
    base = make_base_args(sn_parser, args)
    design = family_design(base, args)
    design.to_csv(outdir / "panelC_family_design.csv", index=False)

    raw_path = outdir / "panelC_SN_nucleation_waterfall_raw.csv"
    old = read_existing_raw(raw_path) if args.resume else pd.DataFrame()
    done = set(point_key(r) for _, r in old.iterrows()) if not old.empty else set()
    rows = old.to_dict("records") if not old.empty else []

    entropy_vals = sorted(set(float(x) for x in args.entropy_mag_kB))
    stresses = sorted(set(float(x) for x in args.stresses_MPa))
    n_total = len(design) * len(entropy_vals) * len(stresses)
    n_counter = 0

    for _, fam in design.iterrows():
        case = SNCase(
            str(fam["family_id"]),
            float(fam["chi_back"]),
            float(fam["Gshield_eV"]),
        )
        print(
            f"=== {fam['family_id']}  Lambda_sh={fam['Lambda_sh']:.3f}  "
            f"Gshield={fam['Gshield_eV']:.4f} eV ==="
        )

        for Smag in entropy_vals:
            for sigma in stresses:
                n_counter += 1
                k = (round(float(fam["Lambda_sh"]), 12), round(Smag, 12), round(sigma, 9))
                if k in done:
                    continue

                q = SimpleNamespace(**vars(base))
                q.S_crack_kB = -float(Smag)
                q.T = float(args.T_K)

                rec, _hist = run_point(q, case, float(sigma))
                row = dict(rec)
                row.update({
                    "family_id": str(fam["family_id"]),
                    "family_index": int(fam["family_index"]),
                    "Lambda_sh": float(fam["Lambda_sh"]),
                    "G0_emit_ref_eV": float(fam["G0_emit_ref_eV"]),
                    "Gshield_eV": float(fam["Gshield_eV"]),
                    "chi_back": float(fam["chi_back"]),
                    "entropy_mag_kB": float(Smag),
                    "S_crack_kB": -float(Smag),
                    "sigma_a_MPa": float(sigma),
                    "T_K": float(args.T_K),
                })

                status = str(row.get("status", ""))
                Ntotal = float(row.get("cycles_total", math.nan))
                Nfail = row.get("cycles_to_nucleation", math.nan)
                try:
                    Nfail = float(Nfail)
                except Exception:
                    Nfail = math.nan
                fired = status == "failed" and math.isfinite(Nfail)
                Ndisplay = Nfail if fired else Ntotal
                if not math.isfinite(Ndisplay) or Ndisplay <= 0:
                    Ndisplay = float(args.cycles_max)
                row["fired"] = bool(fired)
                row["log10_cycles_display"] = math.log10(max(Ndisplay, 1.0))
                row["direct_lt_1_cycle"] = bool(fired and Nfail < 1.0)

                rows.append(row)
                append_csv_row(raw_path, row)
                done.add(k)

                if args.print_every_point:
                    print(
                        f"[{n_counter}/{n_total}] S_mag={Smag:g} kB  "
                        f"sigma={sigma:g} MPa  status={status}  N={Ndisplay:.3e}"
                    )

    raw = pd.DataFrame(rows)
    return raw, design


def stress_at_target_life(g: pd.DataFrame, target_N: float) -> float:
    """Interpolate sigma_a at target life from failed points only.

    This is descriptive and intentionally does not turn censored points into
    artificial failures.
    """
    f = g[(g["status"].astype(str) == "failed") & g["cycles_to_nucleation"].notna()].copy()
    if len(f) < 2:
        return math.nan
    x = np.log10(np.maximum(f["cycles_to_nucleation"].to_numpy(float), 1e-300))
    y = f["sigma_a_MPa"].to_numpy(float)
    order = np.argsort(x)
    x, y = x[order], y[order]
    xt = math.log10(target_N)
    if xt < np.nanmin(x) or xt > np.nanmax(x):
        return math.nan
    return float(np.interp(xt, x, y))


def summarize(raw: pd.DataFrame) -> pd.DataFrame:
    rows = []
    groups = raw.groupby(["family_id", "Lambda_sh", "Gshield_eV", "entropy_mag_kB"], sort=False)
    for keys, g in groups:
        family_id, lam, Gshield, Smag = keys
        failed = g[g["status"].astype(str) == "failed"]
        cens = g[g["status"].astype(str) == "right_censored"]
        rows.append({
            "family_id": family_id,
            "Lambda_sh": float(lam),
            "Gshield_eV": float(Gshield),
            "entropy_mag_kB": float(Smag),
            "S_crack_kB": -float(Smag),
            "n_points": int(len(g)),
            "n_failed": int(len(failed)),
            "n_right_censored": int(len(cens)),
            "sigma_at_1e6_cycles_MPa": stress_at_target_life(g, 1e6),
            "sigma_at_1e8_cycles_MPa": stress_at_target_life(g, 1e8),
            "sigma_at_1e10_cycles_MPa": stress_at_target_life(g, 1e10),
            "sigma_at_1e12_cycles_MPa": stress_at_target_life(g, 1e12),
        })
    return pd.DataFrame(rows)


def plot_waterfall(raw: pd.DataFrame, design: pd.DataFrame, outdir: Path,
                   cycles_max: float, elev: float, azim: float, dpi: int) -> None:
    if raw.empty:
        raise ValueError("No raw data available for plotting.")

    fig = plt.figure(figsize=(10.2, 7.8))
    ax = fig.add_subplot(111, projection="3d")

    lamvals = design["Lambda_sh"].to_numpy(float)
    vmin = float(np.nanmin(lamvals))
    vmax = float(np.nanmax(lamvals))
    if math.isclose(vmin, vmax):
        vmax = vmin + 1.0
    norm = Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.get_cmap("viridis")

    # Plot every S-N curve at fixed (Lambda_sh, entropy magnitude). Color encodes
    # Lambda_sh; the depth axis itself encodes entropy.
    for _, fam in design.sort_values("Lambda_sh").iterrows():
        color = cmap(norm(float(fam["Lambda_sh"])))
        fraw = raw[np.isclose(raw["Lambda_sh"].to_numpy(float), float(fam["Lambda_sh"]), rtol=0, atol=1e-10)]
        for Smag in sorted(fraw["entropy_mag_kB"].dropna().unique()):
            g = fraw[np.isclose(fraw["entropy_mag_kB"].to_numpy(float), float(Smag), rtol=0, atol=1e-10)].copy()
            if g.empty:
                continue
            g = g.sort_values("sigma_a_MPa")
            status = g["status"].astype(str).to_numpy()
            failed = status == "failed"
            cens = status == "right_censored"
            block = ~(failed | cens)

            x = g["log10_cycles_display"].to_numpy(float)
            y = np.full(len(g), float(Smag), float)
            z = g["sigma_a_MPa"].to_numpy(float)

            if np.any(failed):
                # Connecting in stress order preserves the S-N curve topology.
                idx = np.where(failed)[0]
                ax.plot(x[idx], y[idx], z[idx], lw=1.65, color=color, alpha=0.90)
                ax.scatter(x[idx], y[idx], z[idx], s=14, color=[color], alpha=0.95)
            if np.any(cens):
                idx = np.where(cens)[0]
                ax.scatter(x[idx], y[idx], z[idx], marker="^", s=30,
                           facecolors="none", edgecolors=[color], linewidths=0.95)
            if np.any(block):
                idx = np.where(block)[0]
                ax.scatter(x[idx], y[idx], z[idx], marker="x", s=26,
                           color=[color], linewidths=0.9)

        # One family label, placed at the highest entropy / highest-stress point.
        if not fraw.empty:
            q = fraw.sort_values(["entropy_mag_kB", "sigma_a_MPa"]).iloc[-1]
            ax.text(float(q["log10_cycles_display"]) + 0.08,
                    float(q["entropy_mag_kB"]),
                    float(q["sigma_a_MPa"]),
                    rf"$\Lambda_{{sh}}$={float(fam['Lambda_sh']):.2f}",
                    fontsize=7.5, color=color)

    ax.set_xlabel(r"$\log_{10}$ cycles to crack initiation, $N_i$", labelpad=10)
    ax.set_ylabel(r"$-S_{crack}^{*}/k_B$", labelpad=10)
    ax.set_zlabel(r"Stress amplitude $\sigma_a$ [MPa]", labelpad=11)
    ax.set_title(
        "Systematic V1 S–N crack-initiation waterfall\n"
        r"axes: $N_i$, $-S_{crack}^{*}/k_B$, $\sigma_a$; line family parameterized by $\Lambda_{sh}$",
        pad=18,
    )
    ax.view_init(elev=float(elev), azim=float(azim))

    # Keep the visible cycle axis bounded by the requested horizon.
    xmax = math.log10(max(float(cycles_max), 1.0))
    xmin = max(0.0, float(np.nanmin(raw["log10_cycles_display"].to_numpy(float))) - 0.25)
    ax.set_xlim(xmin, xmax + 0.15)

    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.045, pad=0.10)
    cbar.set_label(r"$\Lambda_{sh}=G_{shield}/G_{0,e}(T)$")

    fig.text(0.075, 0.055,
             rf"Open triangles: no initiation by $N_{{max}}={cycles_max:.0e}$ cycles",
             fontsize=9.5)
    fig.tight_layout(rect=(0, 0.075, 0.93, 1))

    png = outdir / "panelC_SN_nucleation_waterfall_3d.png"
    pdf = outdir / "panelC_SN_nucleation_waterfall_3d.pdf"
    fig.savefig(png, dpi=int(dpi), bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)


def plot_2d_projection(raw: pd.DataFrame, design: pd.DataFrame, outdir: Path, dpi: int) -> None:
    """Diagnostic 2-D S-N projection, useful for validating the 3-D waterfall."""
    ent = sorted(raw["entropy_mag_kB"].dropna().unique())
    if not ent:
        return
    # Choose low, middle, high entropy slices only, to avoid a dense diagnostic plot.
    picks = sorted(set([ent[0], ent[len(ent)//2], ent[-1]]))
    fig, axes = plt.subplots(1, len(picks), figsize=(4.6 * len(picks), 4.1), squeeze=False)
    axes = axes[0]
    norm = Normalize(vmin=float(design["Lambda_sh"].min()), vmax=float(design["Lambda_sh"].max()) or 1.0)
    cmap = plt.get_cmap("viridis")

    for ax, Smag in zip(axes, picks):
        for _, fam in design.sort_values("Lambda_sh").iterrows():
            g = raw[
                np.isclose(raw["Lambda_sh"].to_numpy(float), float(fam["Lambda_sh"]), atol=1e-10, rtol=0)
                & np.isclose(raw["entropy_mag_kB"].to_numpy(float), float(Smag), atol=1e-10, rtol=0)
            ].sort_values("sigma_a_MPa")
            if g.empty:
                continue
            color = cmap(norm(float(fam["Lambda_sh"])))
            failed = g["status"].astype(str).to_numpy() == "failed"
            if np.any(failed):
                ax.plot(g.loc[failed, "cycles_to_nucleation"], g.loc[failed, "sigma_a_MPa"],
                        marker="o", ms=3, lw=1.4, color=color,
                        label=rf"$\Lambda_{{sh}}$={float(fam['Lambda_sh']):.2f}")
            cens = g["status"].astype(str).to_numpy() == "right_censored"
            if np.any(cens):
                ax.scatter(g.loc[cens, "cycles_total"], g.loc[cens, "sigma_a_MPa"],
                           marker="^", s=28, facecolors="none", edgecolors=[color])
        ax.set_xscale("log")
        ax.set_xlabel("Cycles to initiation")
        ax.set_ylabel(r"$\sigma_a$ [MPa]")
        ax.set_title(rf"$-S_{{crack}}^*/k_B={Smag:g}$")
        ax.grid(True, alpha=0.25)
    handles, labels = axes[-1].get_legend_handles_labels()
    if handles:
        axes[-1].legend(handles, labels, fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(outdir / "panelC_SN_projection_diagnostic.png", dpi=int(dpi), bbox_inches="tight")
    plt.close(fig)


def write_manifest(args: argparse.Namespace, design: pd.DataFrame, outdir: Path) -> None:
    manifest = {
        "panel": "Figure 1C",
        "model": "reduced V1 fully Arrhenius S-N crack-initiation model",
        "axes": {
            "x": "log10 cycles to crack initiation N_i",
            "y": "-S_crack*/k_B",
            "z": "nominal stress amplitude sigma_a [MPa]",
            "color_family": "Lambda_sh = Gshield/G0_emit(T)",
        },
        "important_design_note": (
            "The entropy axis uses the explicit crack-opening entropy term S_crack_kB. "
            "The EXP-floor temperature slope gT is not used as the fixed-300-K entropy "
            "axis because at T=Tref it would not change the barrier and would produce a "
            "degenerate sweep."
        ),
        "T_K": float(args.T_K),
        "R": float(args.R),
        "frequency_Hz": float(args.frequency_Hz),
        "cycles_max": float(args.cycles_max),
        "chi_back": float(args.chi_back),
        "entropy_mag_kB": [float(x) for x in args.entropy_mag_kB],
        "lambda_sh": [float(x) for x in args.lambda_sh],
        "stresses_MPa": [float(x) for x in args.stresses_MPa],
        "family_design": design.to_dict("records"),
    }
    with (outdir / "panelC_run_manifest.json").open("w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)


def main() -> None:
    args = parse_args()
    outdir = args.out
    outdir.mkdir(parents=True, exist_ok=True)

    raw_path = outdir / "panelC_SN_nucleation_waterfall_raw.csv"
    design_path = outdir / "panelC_family_design.csv"

    if args.plot_only:
        if not raw_path.is_file() or not design_path.is_file():
            raise FileNotFoundError(
                f"--plot-only requires {raw_path} and {design_path}"
            )
        raw = pd.read_csv(raw_path)
        design = pd.read_csv(design_path)
    else:
        raw, design = run_sweep(args, outdir)

    summary = summarize(raw)
    summary.to_csv(outdir / "panelC_SN_curve_summary.csv", index=False)
    write_manifest(args, design, outdir)
    plot_waterfall(raw, design, outdir, args.cycles_max,
                   args.view_elev, args.view_azim, args.figure_dpi)
    plot_2d_projection(raw, design, outdir, args.figure_dpi)

    print("\nPanel C outputs:")
    for name in [
        "panelC_SN_nucleation_waterfall_raw.csv",
        "panelC_family_design.csv",
        "panelC_SN_curve_summary.csv",
        "panelC_run_manifest.json",
        "panelC_SN_nucleation_waterfall_3d.png",
        "panelC_SN_nucleation_waterfall_3d.pdf",
        "panelC_SN_projection_diagnostic.png",
    ]:
        print("  ", outdir / name)


if __name__ == "__main__":
    main()
