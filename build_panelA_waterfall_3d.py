#!/usr/bin/env python3
"""Build a systematic 3-D waterfall plot for V1 toughness response.

This workflow is intended for Figure 1 panel A.  It uses the 1-D
Arrhenius-hazard sharp-front engine to compute Kc(T) along a smooth line family
through the (H0,c, chi_shield, N_sat) control space.

Plot design
-----------
Axes:
    x = temperature T [K]
    y = chi_shield [-]
    z = Kc [MPa sqrt(m)]

Line family parameter:
    H0,c [eV]

The line family is built from a continuation path anchored to the four
canonical regimes:
    ceramic : H0=2.6, chi=0.00, N_sat=inf
    peak    : H0=3.6, chi=0.10, N_sat=inf
    weak-T  : H0=4.0, chi=0.20, N_sat=1500
    DBTT    : H0=6.0, chi=0.60, N_sat=2000

Because all four regimes cannot be traversed by varying a single physical input
alone, the workflow co-varies chi_shield and a hidden recovery/saturation path
N_sat(H0).  chi_shield is plotted explicitly; H0 labels each curve; N_sat is
saved to CSV so the exact path can be audited later.

Outputs
-------
- panel_A_waterfall_path.csv         : line definitions
- panel_A_waterfall_raw.csv          : replot-ready Kc(T) data for every line
- panel_A_waterfall_line_metrics.csv : per-line descriptors
- panel_A_waterfall_3d.png/.pdf      : 3-D waterfall plot
- panel_A_waterfall_H0_chi_path.png  : diagnostic H0-chi path plot
- README_PANELA_WATERFALL_3D.md      : usage and interpretation notes
"""
from __future__ import annotations

import argparse
import math
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401


@dataclass(frozen=True)
class Anchor:
    regime: str
    H0_eV: float
    chi_shield: float
    N_sat: float


DEFAULT_ANCHORS: tuple[Anchor, ...] = (
    Anchor("ceramic", 2.6, 0.00, math.inf),
    Anchor("peak",    3.6, 0.10, math.inf),
    Anchor("weakT",   4.0, 0.20, 1500.0),
    Anchor("dbtt",    6.0, 0.60, 2000.0),
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--project-root", type=Path, default=Path.cwd(),
                   help="Project root containing the arrhenius_fracture package.")
    p.add_argument("--out", type=Path, default=Path("runs/panelA_waterfall_3d"),
                   help="Output directory.")
    p.add_argument("--n-lines", type=int, default=10,
                   help="Number of waterfall curves along the continuation path.")
    p.add_argument("--H0-min-eV", type=float, default=float(DEFAULT_ANCHORS[0].H0_eV))
    p.add_argument("--H0-max-eV", type=float, default=float(DEFAULT_ANCHORS[-1].H0_eV))
    p.add_argument("--T-min", type=float, default=300.0)
    p.add_argument("--T-max", type=float, default=1200.0)
    p.add_argument("--T-step", type=float, default=50.0)
    p.add_argument("--Kdot-MPa-sqrtm-per-s", type=float, default=0.02,
                   help="Monotonic K ramp rate in MPa*sqrt(m)/s.")
    p.add_argument("--Kmax-MPa-sqrtm", type=float, default=80.0,
                   help="Maximum applied K in MPa*sqrt(m).")
    p.add_argument("--dt", type=float, default=1.0,
                   help="Monotonic ramp time step in seconds.")
    p.add_argument("--emit-S-T-c0-kB", type=float, default=-20.0)
    p.add_argument("--emit-S-T-c1", type=float, default=0.02)
    p.add_argument("--emit-S-sigma-max-kB", type=float, default=8.0)
    p.add_argument("--multihit-m", type=float, default=3.0)
    p.add_argument("--multihit-tau", type=float, default=1e-6)
    p.add_argument("--recover-k", type=float, default=0.0,
                   help="Optional explicit recovery term.  Defaults to 0, using only N_sat.")
    p.add_argument("--view-elev", type=float, default=23.0)
    p.add_argument("--view-azim", type=float, default=-63.0)
    p.add_argument("--hide-path-diagnostic", action="store_true",
                   help="Skip the auxiliary H0-chi path diagnostic plot.")
    return p.parse_args()


def _ensure_project_on_path(project_root: Path) -> None:
    root = project_root.resolve()
    if not root.exists():
        raise FileNotFoundError(f"project root not found: {root}")
    sys.path.insert(0, str(root))


def _build_engine_imports(project_root: Path):
    _ensure_project_on_path(project_root)
    from arrhenius_fracture.sharp_front import _build_parser, build_engine  # type: ignore
    from arrhenius_fracture.config import make_emergent_config  # type: ignore
    return _build_parser, build_engine, make_emergent_config


def _inverse_nsat(nsat: float) -> float:
    if nsat is None or not math.isfinite(float(nsat)) or float(nsat) <= 0.0:
        return 0.0
    return 1.0 / float(nsat)


def continuation_path(n_lines: int, H0_min: float, H0_max: float) -> pd.DataFrame:
    anchors = sorted(DEFAULT_ANCHORS, key=lambda a: a.H0_eV)
    H_anchor = np.array([a.H0_eV for a in anchors], dtype=float)
    chi_anchor = np.array([a.chi_shield for a in anchors], dtype=float)
    inv_n_anchor = np.array([_inverse_nsat(a.N_sat) for a in anchors], dtype=float)

    H_values = np.linspace(H0_min, H0_max, n_lines)
    chi_values = np.interp(H_values, H_anchor, chi_anchor)
    inv_n_values = np.interp(H_values, H_anchor, inv_n_anchor)

    nsat_values: list[float] = []
    for inv in inv_n_values:
        if inv <= 1e-12:
            nsat_values.append(math.inf)
        else:
            nsat_values.append(1.0 / float(inv))

    rows = []
    for i, (H0, chi, nsat) in enumerate(zip(H_values, chi_values, nsat_values), start=1):
        # loose regime hint from the nearest canonical anchor in H0-space
        j = int(np.argmin(np.abs(H_anchor - H0)))
        hint = anchors[j].regime
        rows.append({
            "line_id": f"line_{i:02d}",
            "line_index": i,
            "H0_eV": float(H0),
            "chi_shield": float(chi),
            "N_sat": float(nsat) if math.isfinite(nsat) else math.inf,
            "regime_hint": hint,
        })
    return pd.DataFrame(rows)


def build_base_args(_build_parser, args: argparse.Namespace):
    # The underlying parser expects at least a mode and temperatures token.
    base = _build_parser().parse_args(["--mode", "1d", "--temperatures", str(args.T_min)])
    base.emit_S_T_c0_kB = args.emit_S_T_c0_kB
    base.emit_S_T_c1 = args.emit_S_T_c1
    base.emit_S_sigma_max_kB = args.emit_S_sigma_max_kB
    base.multihit_m = args.multihit_m
    base.multihit_tau = args.multihit_tau
    base.emb_sat_frac = 1.0
    base.recover_k = args.recover_k
    return base


def kc_curve_for_line(base_args, build_engine, material, H0_eV: float, chi_shield: float,
                      N_sat: float, T_values: Sequence[float], Kdot_SI: float,
                      Kmax_SI: float, dt: float) -> pd.DataFrame:
    nstep = int(np.ceil(Kmax_SI / (Kdot_SI * dt)))
    rows: list[dict] = []

    for T in T_values:
        a = argparse.Namespace(**vars(base_args))
        a.chi_shield = float(chi_shield)
        a.cleave_H0_eV = float(H0_eV)
        a.N_sat = float(N_sat) if math.isfinite(N_sat) else float("inf")
        eng = build_engine(a, material)

        fired = False
        Kc = math.nan
        for i in range(nstep):
            K = Kdot_SI * (i + 1) * dt
            rec = eng.step(K, float(T), dt)
            if rec.get("fired", False):
                Kc = float(K) / 1e6
                fired = True
                break

        rows.append({
            "T_K": float(T),
            "Kc_pred_MPa_sqrtm": float(Kc) if fired else math.nan,
            "fired": bool(fired),
            "Kc_plot_MPa_sqrtm": float(Kc) if fired else float(Kmax_SI / 1e6),
            "censored_at_Kmax": not fired,
        })

    return pd.DataFrame(rows)


def classify_shape(T: np.ndarray, K: np.ndarray) -> str:
    mask = np.isfinite(T) & np.isfinite(K)
    T = np.asarray(T[mask], float)
    K = np.asarray(K[mask], float)
    if K.size < 3:
        return "indeterminate"
    d = np.gradient(K, T)
    if np.all(d < 0):
        return "monotonic_softening"
    if np.nanmax(K) > K[0] * 1.15 and K[-1] < np.nanmax(K) * 0.85:
        return "peak_like"
    if K[-1] > K[0] * 1.20:
        return "dbtt_like"
    if np.nanmax(np.abs(d)) < 0.01:
        return "weakT_like"
    return "mixed"


def summarize_lines(raw: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    for line_id, g in raw.groupby("line_id", sort=False):
        g = g.sort_values("T_K")
        T = g["T_K"].to_numpy(float)
        K = g["Kc_pred_MPa_sqrtm"].to_numpy(float)
        mask = np.isfinite(K)
        if mask.sum() >= 2:
            d = np.gradient(K[mask], T[mask])
            i = int(np.nanargmax(np.abs(d)))
            T_inflect = float(T[mask][i])
            Kmin = float(np.nanmin(K[mask]))
            Kmax = float(np.nanmax(K[mask]))
            amp = (Kmax - Kmin) / Kmin if Kmin > 0 else math.nan
        else:
            T_inflect = math.nan
            amp = math.nan
        meta = g.iloc[0]
        rows.append({
            "line_id": line_id,
            "line_index": int(meta["line_index"]),
            "H0_eV": float(meta["H0_eV"]),
            "chi_shield": float(meta["chi_shield"]),
            "N_sat": float(meta["N_sat"]) if math.isfinite(float(meta["N_sat"])) else math.inf,
            "n_temperatures": int(len(g)),
            "n_fired": int(np.isfinite(K).sum()),
            "n_censored": int(np.isnan(K).sum()),
            "DBTT_amplitude_index": float(amp) if math.isfinite(amp) else math.nan,
            "T_inflection_K": float(T_inflect) if math.isfinite(T_inflect) else math.nan,
            "shape_hint": classify_shape(T, K),
            "regime_hint": str(meta["regime_hint"]),
        })
    return pd.DataFrame(rows)


def write_readme(outdir: Path, args: argparse.Namespace) -> None:
    text = f"""# Panel A 3-D Waterfall Workflow

This workflow builds a systematic V1 toughness waterfall plot for Figure 1.

## Plot definition
- x-axis: temperature `T` [K]
- y-axis: shielding factor `chi_shield`
- z-axis: fracture toughness `K_c` [MPa sqrt(m)]
- line family parameter: cleavage barrier height `H0,c` [eV]

Each line is a full `K_c(T)` response at fixed `(H0,c, chi_shield, N_sat)`.
The displayed axis is `chi_shield`; `H0,c` varies systematically from line to
line and is also encoded in the colormap and in `panel_A_waterfall_path.csv`.

## Why `N_sat` is included
A single control parameter is not enough to traverse all four canonical
response classes. To move continuously from ceramic-like through peak and
weak-T to DBTT-like behavior, the workflow follows a continuation path through
`(H0,c, chi_shield, N_sat)` space.  `N_sat` is a hidden co-varied parameter and
is saved explicitly to CSV so the path remains auditable.

## Canonical anchors
- ceramic: `H0=2.6 eV`, `chi=0.00`, `N_sat=inf`
- peak: `H0=3.6 eV`, `chi=0.10`, `N_sat=inf`
- weak-T: `H0=4.0 eV`, `chi=0.20`, `N_sat=1500`
- DBTT: `H0=6.0 eV`, `chi=0.60`, `N_sat=2000`

## Run settings used here
- number of lines: {args.n_lines}
- H0 range: {args.H0_min_eV:.3f} to {args.H0_max_eV:.3f} eV
- temperature grid: {args.T_min:g} to {args.T_max:g} K in steps of {args.T_step:g} K
- K ramp rate: {args.Kdot_MPa_sqrtm_per_s:g} MPa sqrt(m)/s
- K max: {args.Kmax_MPa_sqrtm:g} MPa sqrt(m)
- dt: {args.dt:g} s

## Files
- `panel_A_waterfall_path.csv`: line definitions and hidden continuation path
- `panel_A_waterfall_raw.csv`: replot-ready raw Kc(T) values for every line
- `panel_A_waterfall_line_metrics.csv`: per-line descriptors
- `panel_A_waterfall_3d.png` and `.pdf`: main waterfall plot
- `panel_A_waterfall_H0_chi_path.png`: diagnostic path view in control space
"""
    (outdir / "README_PANELA_WATERFALL_3D.md").write_text(text)


def plot_path_diagnostic(path_df: pd.DataFrame, outdir: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    x = path_df["H0_eV"].to_numpy(float)
    y = path_df["chi_shield"].to_numpy(float)
    ns = path_df["N_sat"].to_numpy(float)
    sc = ax.scatter(x, y, c=np.where(np.isfinite(ns), np.log10(ns), np.nan), s=70)
    ax.plot(x, y, lw=1.5, alpha=0.8)
    for _, r in path_df.iterrows():
        ax.annotate(f"{r['line_index']}", (r['H0_eV'], r['chi_shield']),
                    xytext=(4, 4), textcoords="offset points", fontsize=8)
    ax.set_xlabel(r"$H_{{0,c}}$ [eV]")
    ax.set_ylabel(r"$\chi_{{shield}}$")
    ax.set_title(r"Continuation path in $(H_{{0,c}},\chi_{{shield}})$ space")
    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label(r"$\log_{{10}} N_{{sat}}$ for finite points")
    fig.tight_layout()
    fig.savefig(outdir / "panel_A_waterfall_H0_chi_path.png", dpi=220, bbox_inches="tight")
    fig.savefig(outdir / "panel_A_waterfall_H0_chi_path.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_waterfall_3d(raw: pd.DataFrame, path_df: pd.DataFrame, outdir: Path,
                      view_elev: float, view_azim: float, Kmax_plot: float) -> None:
    fig = plt.figure(figsize=(9.6, 7.4))
    ax = fig.add_subplot(111, projection="3d")

    H0_vals = path_df["H0_eV"].to_numpy(float)
    norm = Normalize(vmin=float(np.nanmin(H0_vals)), vmax=float(np.nanmax(H0_vals)))
    cmap = plt.get_cmap("viridis")

    for _, meta in path_df.iterrows():
        g = raw[raw["line_id"] == meta["line_id"]].sort_values("T_K")
        T = g["T_K"].to_numpy(float)
        chi = np.full_like(T, float(meta["chi_shield"]), dtype=float)
        K_plot = g["Kc_plot_MPa_sqrtm"].to_numpy(float)
        color = cmap(norm(float(meta["H0_eV"])))
        ax.plot(T, chi, K_plot, lw=2.0, color=color)

        cens = g["censored_at_Kmax"].to_numpy(bool)
        if np.any(cens):
            ax.scatter(T[cens], chi[cens], K_plot[cens], marker="^", s=28,
                       facecolors="none", edgecolors=[color], linewidths=1.0)
        # end label near highest T point
        ax.text(T[-1] + 8.0, chi[-1], K_plot[-1], f"H0={meta['H0_eV']:.2f}",
                color=color, fontsize=8)

    ax.set_xlabel("Temperature T [K]", labelpad=10)
    ax.set_ylabel(r"$\chi_{shield}$", labelpad=10)
    ax.set_zlabel(r"$K_c$ [MPa$\sqrt{m}$]", labelpad=8)
    ax.set_title("Systematic V1 toughness waterfall\naxes: T, chi_shield, Kc; line family parameterized by H0,c")
    ax.view_init(elev=view_elev, azim=view_azim)
    ax.set_zlim(0.0, max(Kmax_plot, 1.0))

    mappable = ScalarMappable(norm=norm, cmap=cmap)
    mappable.set_array(H0_vals)
    cbar = fig.colorbar(mappable, ax=ax, fraction=0.03, pad=0.08)
    cbar.set_label(r"$H_{0,c}$ [eV]")

    # small note for censored points
    ax.text2D(0.02, 0.04,
              f"Open triangles: no fire by Kmax = {Kmax_plot:g} MPa√m",
              transform=ax.transAxes, fontsize=9)
    fig.tight_layout()
    fig.savefig(outdir / "panel_A_waterfall_3d.png", dpi=240, bbox_inches="tight")
    fig.savefig(outdir / "panel_A_waterfall_3d.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    outdir = args.out.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    _build_parser, build_engine, make_emergent_config = _build_engine_imports(args.project_root)
    base_args = build_base_args(_build_parser, args)
    material = make_emergent_config().material

    path_df = continuation_path(args.n_lines, args.H0_min_eV, args.H0_max_eV)
    path_df.to_csv(outdir / "panel_A_waterfall_path.csv", index=False)

    T_values = np.arange(args.T_min, args.T_max + 0.5 * args.T_step, args.T_step, dtype=float)
    Kdot_SI = args.Kdot_MPa_sqrtm_per_s * 1e6
    Kmax_SI = args.Kmax_MPa_sqrtm * 1e6

    raw_frames: list[pd.DataFrame] = []
    for _, line in path_df.iterrows():
        g = kc_curve_for_line(base_args, build_engine, material,
                              H0_eV=float(line["H0_eV"]),
                              chi_shield=float(line["chi_shield"]),
                              N_sat=float(line["N_sat"]) if math.isfinite(float(line["N_sat"])) else math.inf,
                              T_values=T_values,
                              Kdot_SI=Kdot_SI,
                              Kmax_SI=Kmax_SI,
                              dt=float(args.dt))
        for col in ["line_id", "line_index", "H0_eV", "chi_shield", "N_sat", "regime_hint"]:
            g[col] = line[col]
        raw_frames.append(g)

    raw = pd.concat(raw_frames, ignore_index=True)
    raw.to_csv(outdir / "panel_A_waterfall_raw.csv", index=False)

    metrics = summarize_lines(raw)
    metrics.to_csv(outdir / "panel_A_waterfall_line_metrics.csv", index=False)

    plot_waterfall_3d(raw, path_df, outdir, args.view_elev, args.view_azim, args.Kmax_MPa_sqrtm)
    if not args.hide_path_diagnostic:
        plot_path_diagnostic(path_df, outdir)
    write_readme(outdir, args)

    print(f"Wrote {outdir / 'panel_A_waterfall_3d.png'}")
    print(f"Wrote {outdir / 'panel_A_waterfall_raw.csv'}")
    print(f"Wrote {outdir / 'panel_A_waterfall_path.csv'}")


if __name__ == "__main__":
    main()
