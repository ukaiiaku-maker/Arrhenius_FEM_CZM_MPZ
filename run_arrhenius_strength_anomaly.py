#!/usr/bin/env python3
"""Compute temperature-strength-anomaly metrics from the emission barrier.

The calculation deliberately uses the *same emission free-energy surface* as the
corrected two-barrier fatigue map.  The existing map runner is imported, and its
`make_plastic_barriers(case, S_emit_kB, T_anchor_K)` function supplies the active
emission barrier.

For a fixed kinetic rate target, the strength proxy sigma_y is defined implicitly by

    target_rate = nu_emit exp[-DeltaG_emit(sigma_y,T)/(k_B T)].

By default target_rate is identified with the imposed strain-rate label.  An
optional `--strain-per-event` converts macroscopic strain rate to event rate:

    target_rate = strain_rate / strain_per_event.

Absolute stress magnitudes therefore depend on this projection prefactor, whereas
the temperature/rate trends and anomaly metrics are the main intended outputs.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import math
from pathlib import Path
from typing import Dict, List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def load_module(path: Path):
    spec = importlib.util.spec_from_file_location("two_barrier_map_runner", str(path))
    if spec is None or spec.loader is None:
        raise ImportError(path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    for name in ["read_cases", "make_plastic_barriers"]:
        if not hasattr(mod, name):
            raise AttributeError(f"{path} missing required function {name}")
    return mod


def choose_runner(explicit: str) -> Path:
    candidates = [Path(explicit)] if explicit else []
    candidates += [Path("run_v1_two_barrier_dbtt_fatigue_map_fixed.py"), Path("run_v1_two_barrier_dbtt_fatigue_map.py")]
    for p in candidates:
        if p.exists():
            return p.resolve()
    raise FileNotFoundError("Set --map-runner to the corrected two-barrier map script")


def solve_strength(barrier, T: float, target_rate: float, sigma_max_Pa: float, tol_rel: float = 1e-8):
    def rr(s):
        v = barrier.rate(float(s), float(T))
        return float(np.asarray(v))

    r0 = rr(0.0)
    if r0 >= target_rate:
        return 0.0, "athermal_zero_stress"
    rmax = rr(sigma_max_Pa)
    if rmax < target_rate:
        return np.nan, "above_sigma_search_limit"
    lo, hi = 0.0, float(sigma_max_Pa)
    for _ in range(140):
        mid = 0.5 * (lo + hi)
        if rr(mid) < target_rate:
            lo = mid
        else:
            hi = mid
        if (hi - lo) <= tol_rel * max(hi, 1.0):
            break
    return hi, "solved"


def largest_drawup(sig: np.ndarray, T: np.ndarray):
    best = -np.inf
    i0 = i1 = 0
    min_val = sig[0]
    min_i = 0
    for j in range(1, len(sig)):
        if sig[j] - min_val > best:
            best = sig[j] - min_val
            i0, i1 = min_i, j
        if sig[j] < min_val:
            min_val = sig[j]
            min_i = j
    if not math.isfinite(best):
        best = 0.0
    return float(max(best, 0.0)), float(T[i0]), float(T[i1])


def anomaly_metrics(g: pd.DataFrame, plateau_dlnsigma_dT: float) -> Dict:
    g = g.sort_values("T_K")
    ok = g[np.isfinite(g["sigma_y_GPa"])].copy()
    if len(ok) < 3:
        return {
            "n_solved": len(ok), "anomaly_amplitude_GPa": np.nan,
            "T_anomaly_onset_K": np.nan, "T_peak_K": np.nan,
            "positive_slope_area_GPa": np.nan, "plateau_width_K": np.nan,
            "sigma_min_GPa": np.nan, "sigma_max_GPa": np.nan,
        }
    T = ok["T_K"].to_numpy(float)
    s = ok["sigma_y_GPa"].to_numpy(float)
    amp, Ton, Tpk_draw = largest_drawup(s, T)
    ds_dT = np.gradient(s, T)
    positive_area = float(np.trapezoid(np.maximum(ds_dT, 0.0), T))
    with np.errstate(divide="ignore", invalid="ignore"):
        dln = np.gradient(np.log(np.maximum(s, 1e-12)), T)
    mask = np.abs(dln) <= plateau_dlnsigma_dT
    plateau_width = 0.0
    start = None
    for i, flag in enumerate(mask):
        if flag and start is None:
            start = i
        if (not flag or i == len(mask)-1) and start is not None:
            end = i if flag and i == len(mask)-1 else i-1
            plateau_width = max(plateau_width, T[end] - T[start])
            start = None
    return {
        "n_solved": int(len(ok)),
        "anomaly_amplitude_GPa": amp,
        "T_anomaly_onset_K": Ton,
        "T_peak_K": Tpk_draw,
        "positive_slope_area_GPa": positive_area,
        "plateau_width_K": float(plateau_width),
        "sigma_min_GPa": float(np.nanmin(s)),
        "sigma_max_GPa": float(np.nanmax(s)),
        "mean_positive_slope_GPa_per_K": float(np.mean(np.maximum(ds_dT, 0.0))),
    }


def plot_strength_curves(curves: pd.DataFrame, out: Path):
    rates = sorted(curves["strain_rate_s-1"].unique())
    ncols = 2
    nrows = int(math.ceil(len(rates) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(13, 4.8*nrows), squeeze=False)
    for ax, rate in zip(axes.flat, rates):
        gr = curves[np.isclose(curves["strain_rate_s-1"], rate)]
        for (case, Se), g in gr.groupby(["case_label", "S_emit_kB"], sort=False):
            g = g.sort_values("T_K")
            ax.plot(g["T_K"], g["sigma_y_GPa"], label=f"{case}, Se={Se:g}", linewidth=1.4)
        ax.set_title(rf"Kinetic rate label $\dot{{\epsilon}}={rate:g}$ s$^{{-1}}$")
        ax.set_xlabel("Temperature (K)")
        ax.set_ylabel(r"Arrhenius strength proxy $\sigma_y$ (GPa)")
        ax.grid(True, alpha=0.25)
        if len(gr.groupby(["case_label", "S_emit_kB"])) <= 14:
            ax.legend(fontsize=6.2, ncol=2)
    for ax in axes.flat[len(rates):]:
        ax.axis("off")
    fig.suptitle("Temperature-strength-anomaly response from the emission barrier", y=0.998)
    fig.tight_layout()
    fig.savefig(out / "strength_temperature_anomaly_curves.png", dpi=250)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--map-runner", default="")
    ap.add_argument("--case-table", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--case-filter", nargs="*", default=[])
    ap.add_argument("--emission-entropies-kB", nargs="+", type=float, required=True)
    ap.add_argument("--temperatures", nargs="+", type=float, required=True)
    ap.add_argument("--strain-rates", nargs="+", type=float, default=[1e-4, 1e-2, 1.0, 100.0])
    ap.add_argument("--strain-per-event", type=float, default=1.0,
                    help="Projection factor: target emission-event rate = strain_rate/strain_per_event.")
    ap.add_argument("--T-anchor-K", type=float, default=300.0)
    ap.add_argument("--sigma-max-GPa", type=float, default=250.0)
    ap.add_argument("--plateau-dlnsigma-dT", type=float, default=2.5e-4,
                    help="Plateau criterion |d ln sigma/dT| below this value [1/K].")
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    runner_path = choose_runner(args.map_runner)
    runner = load_module(runner_path)
    cases = runner.read_cases(Path(args.case_table))
    if args.case_filter:
        wanted = set(args.case_filter)
        cases = [c for c in cases if str(c["case_label"]) in wanted]
    rows: List[Dict] = []
    for case in cases:
        for Se in args.emission_entropies_kB:
            emit, _, _ = runner.make_plastic_barriers(case, Se, args.T_anchor_K)
            for edot in args.strain_rates:
                target_rate = float(edot) / max(args.strain_per_event, 1e-300)
                for T in args.temperatures:
                    sigma, status = solve_strength(emit, T, target_rate, args.sigma_max_GPa*1e9)
                    rows.append({
                        "case_label": case["case_label"],
                        "S_emit_kB": Se,
                        "strain_rate_s-1": edot,
                        "target_event_rate_s-1": target_rate,
                        "T_K": T,
                        "sigma_y_GPa": sigma/1e9 if math.isfinite(sigma) else np.nan,
                        "status": status,
                    })
    curves = pd.DataFrame(rows)
    metrics = []
    for keys, g in curves.groupby(["case_label", "S_emit_kB", "strain_rate_s-1"], sort=False):
        rec = {"case_label": keys[0], "S_emit_kB": keys[1], "strain_rate_s-1": keys[2]}
        rec.update(anomaly_metrics(g, args.plateau_dlnsigma_dT))
        metrics.append(rec)
    metrics_df = pd.DataFrame(metrics)
    curves.to_csv(out / "arrhenius_strength_curves.csv", index=False)
    metrics_df.to_csv(out / "strength_anomaly_metrics.csv", index=False)
    with open(out / "strength_model_settings.json", "w") as f:
        json.dump({
            "map_runner": str(runner_path),
            "case_table": str(Path(args.case_table).resolve()),
            "emission_entropies_kB": args.emission_entropies_kB,
            "temperatures": args.temperatures,
            "strain_rates": args.strain_rates,
            "strain_per_event": args.strain_per_event,
            "T_anchor_K": args.T_anchor_K,
            "sigma_max_GPa": args.sigma_max_GPa,
            "plateau_dlnsigma_dT": args.plateau_dlnsigma_dT,
            "interpretation": "kinetic strength proxy from the same emission barrier; absolute stress depends on rate-to-event projection",
        }, f, indent=2)
    plot_strength_curves(curves, out)
    print(f"wrote {out / 'arrhenius_strength_curves.csv'}")
    print(f"wrote {out / 'strength_anomaly_metrics.csv'}")


if __name__ == "__main__":
    main()
