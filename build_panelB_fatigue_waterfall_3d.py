#!/usr/bin/env python3
"""Build a 3-D V1 fatigue-response waterfall analogous to Panel A.

Axes
----
    x = log10(cycles to first crack advance)
    y = chi_shield
    z = DeltaK [MPa sqrt(m)]

Line-family parameter
---------------------
    H0,c [eV]

The same continuation path used by the Panel-A toughness waterfall is used
here, so each fatigue curve has fixed (H0,c, chi_shield, N_sat) and the line
family traverses the ceramic -> peak -> weak-T -> DBTT control path.

The calculation uses the V1 cycle-integrated Arrhenius fatigue controller and
records first crack advance.  Points that do not fire by the cycle horizon are
kept as right-censored observations and plotted as open triangles at N_max.
Points that fire in less than one cycle are retained in the CSV and plotted at
N=1 only for display.
"""
from __future__ import annotations

import argparse
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401


@dataclass(frozen=True)
class Anchor:
    regime: str
    H0_eV: float
    chi_shield: float
    N_sat: float


DEFAULT_ANCHORS = (
    Anchor("ceramic", 2.6, 0.00, math.inf),
    Anchor("peak",    3.6, 0.10, math.inf),
    Anchor("weakT",   4.0, 0.20, 1500.0),
    Anchor("dbtt",    6.0, 0.60, 2000.0),
)


def parse_float_list(text: str) -> list[float]:
    vals = [float(x) for x in text.replace(",", " ").split()]
    if not vals:
        raise argparse.ArgumentTypeError("expected at least one numeric value")
    return vals


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--project-root", type=Path, default=Path.cwd())
    p.add_argument("--out", type=Path, default=Path("runs/panelB_fatigue_waterfall_3d"))
    p.add_argument("--n-lines", type=int, default=10)
    p.add_argument("--H0-min-eV", type=float, default=2.6)
    p.add_argument("--H0-max-eV", type=float, default=6.0)
    p.add_argument("--Kmax-grid", type=parse_float_list,
                   default=parse_float_list("6 8 10 12 14 16 18 20 22 24"),
                   help="Quoted list of Kmax values in MPa sqrt(m).")
    p.add_argument("--T-K", type=float, default=300.0)
    p.add_argument("--R", type=float, default=0.1)
    p.add_argument("--frequency-Hz", type=float, default=1000.0)
    p.add_argument("--cycles-max", type=float, default=1e12)
    p.add_argument("--max-blocks", type=int, default=20000)
    p.add_argument("--n-phase", type=int, default=96)
    p.add_argument("--target-dB", type=float, default=0.02)
    p.add_argument("--target-dN-store", type=float, default=0.01)
    p.add_argument("--min-block-cycles", type=float, default=1e-6)
    p.add_argument("--print-every-case", action="store_true")

    # Production fatigue-barrier scaling used in the current V1 family.
    p.add_argument("--emit-energy-scale", type=float, default=0.75)
    p.add_argument("--emit-entropy-scale", type=float, default=0.75)
    p.add_argument("--peierls-energy-scale", type=float, default=0.00375)
    p.add_argument("--peierls-entropy-scale", type=float, default=0.00375)
    p.add_argument("--taylor-energy-scale", type=float, default=0.015)
    p.add_argument("--taylor-entropy-scale", type=float, default=0.015)

    p.add_argument("--view-elev", type=float, default=24.0)
    p.add_argument("--view-azim", type=float, default=-62.0)
    return p.parse_args()


def continuation_path(n_lines: int, H0_min: float, H0_max: float) -> pd.DataFrame:
    anchors = sorted(DEFAULT_ANCHORS, key=lambda a: a.H0_eV)
    H_anchor = np.array([a.H0_eV for a in anchors], float)
    chi_anchor = np.array([a.chi_shield for a in anchors], float)
    inv_n_anchor = np.array([
        0.0 if not math.isfinite(a.N_sat) else 1.0 / a.N_sat
        for a in anchors
    ], float)

    H = np.linspace(H0_min, H0_max, n_lines)
    chi = np.interp(H, H_anchor, chi_anchor)
    inv_n = np.interp(H, H_anchor, inv_n_anchor)
    Nsat = [math.inf if x <= 1e-12 else 1.0 / float(x) for x in inv_n]

    rows = []
    for i, (h, c, n) in enumerate(zip(H, chi, Nsat), start=1):
        j = int(np.argmin(np.abs(H_anchor - h)))
        rows.append({
            "line_id": f"line_{i:02d}",
            "line_index": i,
            "H0_eV": float(h),
            "chi_shield": float(c),
            "N_sat": float(n) if math.isfinite(n) else math.inf,
            "regime_hint": anchors[j].regime,
        })
    return pd.DataFrame(rows)


def import_v1(project_root: Path):
    root = project_root.resolve()
    sys.path.insert(0, str(root))
    from arrhenius_fracture.fatigue_sharp_front import (  # type: ignore
        build_parser, _build_front, _build_controller,
    )
    from arrhenius_fracture.config import ElasticProperties  # type: ignore
    from arrhenius_fracture.fatigue_v1 import FatigueWaveform  # type: ignore
    return build_parser, _build_front, _build_controller, ElasticProperties, FatigueWaveform


def make_base_args(build_parser, args: argparse.Namespace):
    a = build_parser().parse_args([])
    a.R = float(args.R)
    a.frequency_Hz = float(args.frequency_Hz)
    a.cycles_max = float(args.cycles_max)
    # Let the adaptive controller take the full remaining horizon when hazards
    # are tiny; target increments reduce the block automatically near events.
    a.block_cycles = float(args.cycles_max)
    a.max_block_cycles = float(args.cycles_max)
    a.min_block_cycles = float(args.min_block_cycles)
    a.target_dB = float(args.target_dB)
    a.target_dN_store = float(args.target_dN_store)
    a.n_phase = int(args.n_phase)
    a.max_blocks = int(args.max_blocks)
    a.print_every = 0
    a.no_plots = True
    a.storage_model = "escape_limited"

    a.emit_energy_scale = float(args.emit_energy_scale)
    a.emit_entropy_scale = float(args.emit_entropy_scale)
    a.peierls_energy_scale = float(args.peierls_energy_scale)
    a.peierls_entropy_scale = float(args.peierls_entropy_scale)
    a.taylor_energy_scale = float(args.taylor_energy_scale)
    a.taylor_entropy_scale = float(args.taylor_entropy_scale)

    # Newer production variants may use these additional adaptive controls.
    # They are harmless attributes for older V1 implementations.
    a.cycle_block_mode = "hazard_limited"
    a.target_dN_emit = 0.2
    a.target_dN_mobile = 0.2
    a.target_dN_escape = math.inf
    a.target_dN_peierls = math.inf
    a.target_dN_taylor = math.inf
    return a


def run_one_case(base_args, _build_front, _build_controller, ElasticProperties,
                 FatigueWaveform, *, H0_eV: float, chi: float, N_sat: float,
                 Kmax_MPa: float, T_K: float, cycles_max: float, max_blocks: int) -> dict:
    a = argparse.Namespace(**vars(base_args))
    a.cleave_H0_eV = float(H0_eV)
    a.cleave_shield_chi = float(chi)
    a.N_sat = float(N_sat) if math.isfinite(N_sat) else math.inf
    a.Kmax_MPa_sqrt_m = float(Kmax_MPa)

    mat = ElasticProperties()
    front = _build_front(a, mat)
    controller = _build_controller(a)
    wave = FatigueWaveform(
        Kmax=float(Kmax_MPa) * 1e6,
        R=float(a.R),
        frequency_Hz=float(a.frequency_Hz),
        closure_clip=not bool(getattr(a, "no_closure_clip", False)),
    )

    cycles = 0.0
    fired = False
    n_blocks = 0
    for ib in range(max_blocks):
        if cycles >= cycles_max:
            break
        req = min(float(cycles_max - cycles), float(getattr(a, "block_cycles", cycles_max)))
        row = controller.cycle_step_front(front, wave, float(T_K), requested_cycles=req)
        dN = float(row["cycles"])
        if not math.isfinite(dN) or dN <= 0.0:
            raise RuntimeError(
                f"invalid cycle increment {dN} for H0={H0_eV}, chi={chi}, K={Kmax_MPa}"
            )
        cycles += dN
        n_blocks = ib + 1
        if bool(row.get("fired", False)):
            fired = True
            break

    tol = max(1e-8 * cycles_max, 1e-6)
    if fired:
        status = "failed"
    elif cycles >= cycles_max - tol:
        status = "right_censored"
    else:
        status = "block_limited"

    Kmax = float(Kmax_MPa)
    DeltaK = Kmax * (1.0 - float(a.R))
    return {
        "Kmax_MPa_sqrtm": Kmax,
        "DeltaK_MPa_sqrtm": DeltaK,
        "T_K": float(T_K),
        "R": float(a.R),
        "frequency_Hz": float(a.frequency_Hz),
        "cycles_total": float(cycles),
        "cycles_to_first_advance": float(cycles) if fired else math.nan,
        "status": status,
        "fired": bool(fired),
        "n_blocks": int(n_blocks),
        "direct_lt_1_cycle": bool(fired and cycles < 1.0),
        "log10_cycles_display": math.log10(max(float(cycles), 1.0)),
        "N_em_final": float(front.N_em),
        "B_final": float(front.B),
        "a_adv_m": float(front.a_adv),
        "n_adv": int(front.n_adv),
    }


def stress_at_life(g: pd.DataFrame, target_N: float) -> float:
    """Censoring-aware lower estimate of DeltaK at a target life.

    We use failed points bracketing the target when available and retain a
    simple midpoint in DeltaK. This is a descriptive atlas metric, not a fitted
    endurance limit.
    """
    g = g.sort_values("DeltaK_MPa_sqrtm")
    x = g["DeltaK_MPa_sqrtm"].to_numpy(float)
    N = g["cycles_total"].to_numpy(float)
    status = g["status"].astype(str).to_numpy()

    safe = [float(k) for k, n, s in zip(x, N, status)
            if (s == "right_censored" and n >= target_N) or (s == "failed" and n >= target_N)]
    fail = [float(k) for k, n, s in zip(x, N, status)
            if s == "failed" and n < target_N]
    if not safe or not fail:
        return math.nan
    lo = max(safe)
    hi_candidates = [k for k in fail if k >= lo]
    if not hi_candidates:
        return math.nan
    hi = min(hi_candidates)
    return 0.5 * (lo + hi)


def summarize_lines(raw: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for line_id, g in raw.groupby("line_id", sort=False):
        meta = g.iloc[0]
        failed = g[g["status"] == "failed"].copy()
        if len(failed) >= 2:
            xx = np.log10(np.maximum(failed["cycles_to_first_advance"].to_numpy(float), 1e-300))
            yy = failed["DeltaK_MPa_sqrtm"].to_numpy(float)
            order = np.argsort(xx)
            xx, yy = xx[order], yy[order]
            # high-cycle slope over the longest-life half of failed points
            keep = slice(max(0, len(xx) // 2), None)
            if len(xx[keep]) >= 2:
                slope = float(np.polyfit(xx[keep], yy[keep], 1)[0])
            else:
                slope = math.nan
        else:
            slope = math.nan
        rows.append({
            "line_id": line_id,
            "line_index": int(meta["line_index"]),
            "H0_eV": float(meta["H0_eV"]),
            "chi_shield": float(meta["chi_shield"]),
            "N_sat": float(meta["N_sat"]) if math.isfinite(float(meta["N_sat"])) else math.inf,
            "regime_hint": str(meta["regime_hint"]),
            "n_failed": int((g["status"] == "failed").sum()),
            "n_right_censored": int((g["status"] == "right_censored").sum()),
            "n_block_limited": int((g["status"] == "block_limited").sum()),
            "DeltaK_at_1e6_cycles": stress_at_life(g, 1e6),
            "DeltaK_at_1e8_cycles": stress_at_life(g, 1e8),
            "DeltaK_at_1e10_cycles": stress_at_life(g, 1e10),
            "DeltaK_at_1e12_cycles": stress_at_life(g, 1e12),
            "high_cycle_DeltaK_per_decade": slope,
        })
    return pd.DataFrame(rows)


def plot_waterfall(raw: pd.DataFrame, path: pd.DataFrame, outdir: Path,
                   cycles_max: float, elev: float, azim: float) -> None:
    fig = plt.figure(figsize=(9.8, 7.5))
    ax = fig.add_subplot(111, projection="3d")

    H0 = path["H0_eV"].to_numpy(float)
    norm = Normalize(vmin=float(np.min(H0)), vmax=float(np.max(H0)))
    cmap = plt.get_cmap("viridis")

    for _, meta in path.iterrows():
        g = raw[raw["line_id"] == meta["line_id"]].sort_values("DeltaK_MPa_sqrtm")
        color = cmap(norm(float(meta["H0_eV"])))
        x = g["log10_cycles_display"].to_numpy(float)
        y = np.full(len(g), float(meta["chi_shield"]), float)
        z = g["DeltaK_MPa_sqrtm"].to_numpy(float)

        failed = g["status"].to_numpy(str) == "failed"
        cens = g["status"].to_numpy(str) == "right_censored"
        block = g["status"].to_numpy(str) == "block_limited"

        # Connect only measured first-advance lives.  Censored points are shown
        # separately at Nmax so the horizon does not create a false vertical
        # segment that could be mistaken for an endurance plateau.
        if np.any(failed):
            of = np.argsort(x[failed])
            ax.plot(x[failed][of], y[failed][of], z[failed][of], lw=2.0, color=color)
            ax.scatter(x[failed], y[failed], z[failed], s=24, color=[color])
        if np.any(cens):
            ax.scatter(x[cens], y[cens], z[cens], marker="^", s=38,
                       facecolors="none", edgecolors=[color], linewidths=1.1)
            # One dashed bridge from the longest-life measured point to the
            # highest-DeltaK censored point gives a visual bracket only.
            if np.any(failed):
                jf = np.argmax(x[failed])
                cens_indices = np.where(cens)[0]
                jc = cens_indices[np.argmax(z[cens])]
                ax.plot([x[failed][jf], x[jc]],
                        [y[failed][jf], y[jc]],
                        [z[failed][jf], z[jc]],
                        lw=1.0, ls="--", color=color, alpha=0.75)
        if np.any(block):
            ax.scatter(x[block], y[block], z[block], marker="x", s=34,
                       color=[color], linewidths=1.0)

        # label next to the highest-DeltaK censored bracket if present, otherwise
        # at the longest-life measured point.
        if np.any(cens):
            candidates = np.where(cens)[0]
            j = candidates[np.argmax(z[cens])]
        else:
            j = int(np.argmax(x))
        ax.text(x[j] + 0.1, y[j], z[j], f"H0={meta['H0_eV']:.2f}",
                fontsize=8, color=color)

    ax.set_xlabel(r"$\log_{10}$ cycles to first advance", labelpad=10)
    ax.set_ylabel(r"$\chi_{shield}$", labelpad=10)
    ax.set_zlabel(r"$\Delta K$ [MPa$\sqrt{m}$]", labelpad=8)
    T = float(raw["T_K"].iloc[0])
    R = float(raw["R"].iloc[0])
    f = float(raw["frequency_Hz"].iloc[0])
    ax.set_title(
        f"Systematic V1 fatigue-response waterfall\nT={T:g} K, R={R:g}, f={f:g} Hz; line family parameterized by H0,c"
    )
    ax.view_init(elev=elev, azim=azim)
    ax.set_xlim(0.0, math.log10(cycles_max))

    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array(H0)
    cbar = fig.colorbar(sm, ax=ax, fraction=0.03, pad=0.08)
    cbar.set_label(r"$H_{0,c}$ [eV]")

    ax.text2D(0.02, 0.045,
              f"Open triangles: right-censored at Nmax={cycles_max:.0e}; x: block-limited",
              transform=ax.transAxes, fontsize=9)
    fig.tight_layout()
    fig.savefig(outdir / "panel_B_fatigue_waterfall_3d.png", dpi=240, bbox_inches="tight")
    fig.savefig(outdir / "panel_B_fatigue_waterfall_3d.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_2d_projection(raw: pd.DataFrame, path: pd.DataFrame, outdir: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    H0 = path["H0_eV"].to_numpy(float)
    norm = Normalize(vmin=float(np.min(H0)), vmax=float(np.max(H0)))
    cmap = plt.get_cmap("viridis")
    for _, meta in path.iterrows():
        g = raw[raw["line_id"] == meta["line_id"]].copy()
        color = cmap(norm(float(meta["H0_eV"])))
        x = np.maximum(g["cycles_total"].to_numpy(float), 1.0)
        z = g["DeltaK_MPa_sqrtm"].to_numpy(float)
        failed = g["status"].astype(str).to_numpy() == "failed"
        cens = g["status"].astype(str).to_numpy() == "right_censored"
        if np.any(failed):
            of = np.argsort(x[failed])
            ax.plot(x[failed][of], z[failed][of], lw=1.5, color=color, alpha=0.9)
        if np.any(cens):
            ax.scatter(x[cens], z[cens], marker="^", facecolors="none",
                       edgecolors=[color], s=34)
    ax.set_xscale("log")
    ax.set_xlabel("Cycles to first advance")
    ax.set_ylabel(r"$\Delta K$ [MPa$\sqrt{m}$]")
    ax.set_title("2-D projection of the fatigue waterfall")
    ax.grid(True, which="both", alpha=0.25)
    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array(H0)
    cbar = fig.colorbar(sm, ax=ax)
    cbar.set_label(r"$H_{0,c}$ [eV]")
    fig.tight_layout()
    fig.savefig(outdir / "panel_B_fatigue_waterfall_2d_projection.png", dpi=220, bbox_inches="tight")
    fig.savefig(outdir / "panel_B_fatigue_waterfall_2d_projection.pdf", bbox_inches="tight")
    plt.close(fig)


def write_readme(outdir: Path, args: argparse.Namespace) -> None:
    text = f"""# Panel B V1 Fatigue Waterfall

This is the fatigue-response analogue of the Panel-A toughness waterfall.

## Displayed axes
- x: `log10(cycles to first crack advance)`
- y: `chi_shield`
- z: `DeltaK` in MPa sqrt(m)

The line-family parameter is `H0,c`, shown by color and recorded for every line.
The same continuation path through `(H0,c, chi_shield, N_sat)` is used as in
Panel A, allowing the two panels to be interpreted together.

## Fixed fatigue conditions
- T = {args.T_K:g} K
- R = {args.R:g}
- frequency = {args.frequency_Hz:g} Hz
- cycle horizon = {args.cycles_max:.6g}
- Kmax grid = {' '.join(f'{x:g}' for x in args.Kmax_grid)} MPa sqrt(m)

## Censoring
- `failed`: first crack advance observed; life is measured.
- `right_censored`: no advance by Nmax; plotted as an open triangle at Nmax.
- `block_limited`: maximum block count reached before failure or Nmax; plotted as x.
- sub-cycle fires are preserved in CSV but displayed at N=1 in the waterfall.

## Outputs
- `panel_B_fatigue_waterfall_raw.csv`: all V1 fatigue points
- `panel_B_fatigue_waterfall_path.csv`: exact continuation path
- `panel_B_fatigue_waterfall_line_metrics.csv`: stress-at-life and slope descriptors
- `panel_B_fatigue_waterfall_3d.png/.pdf`: main panel
- `panel_B_fatigue_waterfall_2d_projection.png/.pdf`: conventional projection
"""
    (outdir / "README_PANELB_FATIGUE_WATERFALL.md").write_text(text)


def main() -> None:
    args = parse_args()
    outdir = args.out.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    build_parser, _build_front, _build_controller, ElasticProperties, FatigueWaveform = import_v1(args.project_root)
    base = make_base_args(build_parser, args)
    path = continuation_path(args.n_lines, args.H0_min_eV, args.H0_max_eV)
    path.to_csv(outdir / "panel_B_fatigue_waterfall_path.csv", index=False)

    rows = []
    total = len(path) * len(args.Kmax_grid)
    count = 0
    for _, line in path.iterrows():
        for Kmax in args.Kmax_grid:
            count += 1
            rec = run_one_case(
                base, _build_front, _build_controller, ElasticProperties, FatigueWaveform,
                H0_eV=float(line["H0_eV"]),
                chi=float(line["chi_shield"]),
                N_sat=float(line["N_sat"]) if math.isfinite(float(line["N_sat"])) else math.inf,
                Kmax_MPa=float(Kmax),
                T_K=float(args.T_K),
                cycles_max=float(args.cycles_max),
                max_blocks=int(args.max_blocks),
            )
            for col in ["line_id", "line_index", "H0_eV", "chi_shield", "N_sat", "regime_hint"]:
                rec[col] = line[col]
            rows.append(rec)
            if args.print_every_case:
                print(
                    f"[{count}/{total}] {line['line_id']} H0={line['H0_eV']:.3f} "
                    f"chi={line['chi_shield']:.3f} Kmax={Kmax:g} "
                    f"N={rec['cycles_total']:.3e} {rec['status']}"
                )
        # incremental checkpoint after each complete line
        pd.DataFrame(rows).to_csv(outdir / "panel_B_fatigue_waterfall_raw.partial.csv", index=False)
        print(f"completed {line['line_id']} ({count}/{total} cases)")

    raw = pd.DataFrame(rows)
    raw.to_csv(outdir / "panel_B_fatigue_waterfall_raw.csv", index=False)
    partial = outdir / "panel_B_fatigue_waterfall_raw.partial.csv"
    if partial.exists():
        partial.unlink()

    metrics = summarize_lines(raw)
    metrics.to_csv(outdir / "panel_B_fatigue_waterfall_line_metrics.csv", index=False)

    plot_waterfall(raw, path, outdir, args.cycles_max, args.view_elev, args.view_azim)
    plot_2d_projection(raw, path, outdir)
    write_readme(outdir, args)

    print(f"Wrote {outdir / 'panel_B_fatigue_waterfall_3d.png'}")
    print(f"Wrote {outdir / 'panel_B_fatigue_waterfall_raw.csv'}")


if __name__ == "__main__":
    main()
