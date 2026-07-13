#!/usr/bin/env python3
"""Adaptive V1 two-barrier study for Kc(T) and rate-defined fatigue thresholds.

This driver imports the existing corrected two-barrier map runner and reuses its
exact `run_fatigue` and `run_monotonic` functions.  It changes only the sampling
strategy: DeltaK is bracketed adaptively around one or more da/dN criteria instead
of evaluating a large uniform grid.

The existing map runner must expose:
  read_cases(path)
  run_fatigue(case, S_emit_kB, S_cleave_kB, T_K, Kmax_MPa_sqrtm, args)
  run_monotonic(case, S_emit_kB, S_cleave_kB, T_K, args)

The script is resumable.  Every completed V1 fatigue point is appended to
`fatigue_adaptive_points.csv`; monotonic points are written incrementally as well.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import math
from pathlib import Path
from types import SimpleNamespace
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

GROUP = ["case_label", "S_emit_kB", "S_cleave_kB", "T_K"]


def load_module(path: Path):
    spec = importlib.util.spec_from_file_location("two_barrier_map_runner", str(path))
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot import {path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    for name in ["read_cases", "run_fatigue", "run_monotonic"]:
        if not hasattr(mod, name):
            raise AttributeError(f"{path} does not define required function {name}")
    return mod


def choose_runner(explicit: str) -> Path:
    candidates = []
    if explicit:
        candidates.append(Path(explicit))
    candidates += [
        Path("run_v1_two_barrier_dbtt_fatigue_map_fixed.py"),
        Path("run_v1_two_barrier_dbtt_fatigue_map.py"),
    ]
    for p in candidates:
        if p.exists():
            return p.resolve()
    raise FileNotFoundError(
        "Could not find a two-barrier map runner. Set --map-runner to the corrected map script."
    )


def safe_num(v, default=np.nan):
    try:
        x = float(v)
        return x if math.isfinite(x) else default
    except Exception:
        return default


def effective_rate(rec: Dict) -> Tuple[float, str]:
    rate = safe_num(rec.get("da_dN_m_per_cycle"))
    if math.isfinite(rate) and rate > 0:
        return rate, "measured"
    ub = safe_num(rec.get("da_dN_upper_bound_m_per_cycle"))
    status = str(rec.get("status", ""))
    if math.isfinite(ub) and ub > 0 and "horizon" in status:
        return ub, "upper_bound"
    return np.nan, "unresolved"


def settings_namespace(args) -> SimpleNamespace:
    # Include both the fields used by the original corrected map runner and
    # fields added in subsequent adaptive-controller patches.  Extra attributes
    # are harmless and make the driver compatible across nearby code versions.
    return SimpleNamespace(
        T_anchor_K=args.T_anchor_K,
        da_m=args.da_m,
        n_phase=args.n_phase,
        target_dB=args.target_dB,
        target_dN_store=args.target_dN_store,
        R=args.R,
        frequency_Hz=args.frequency_Hz,
        max_blocks=args.max_blocks,
        cycles_max=args.cycles_max,
        n_advances=args.n_advances,
        Kdot_MPa_sqrtm_per_s=args.Kdot_MPa_sqrtm_per_s,
        monotonic_Kmax_MPa=args.monotonic_Kmax_MPa,
        monotonic_dK_MPa=args.monotonic_dK_MPa,
        max_block_cycles=args.max_block_cycles,
        min_block_cycles=args.min_block_cycles,
        target_state_fraction=args.target_state_fraction,
        saturation_tol_fraction=args.saturation_tol_fraction,
        target_dN_store_unbounded=args.target_dN_store_unbounded,
        target_dN_emit=args.target_dN_emit,
        target_dN_mobile=args.target_dN_mobile,
    )


def key_f(case, Se, Sc, T, Kmax):
    return (str(case), round(float(Se), 9), round(float(Sc), 9), round(float(T), 9), round(float(Kmax), 9))


def key_m(case, Se, Sc, T):
    return (str(case), round(float(Se), 9), round(float(Sc), 9), round(float(T), 9))


def save_df(rows: List[Dict], path: Path):
    pd.DataFrame(rows).to_csv(path, index=False)


def locate_bracket(points: List[Dict], criterion: float):
    vals = []
    for r in points:
        eff, source = effective_rate(r)
        if math.isfinite(eff) and eff > 0:
            vals.append((float(r["DeltaK_MPa_sqrtm"]), eff, source, r))
    vals.sort(key=lambda x: x[0])
    for a, b in zip(vals[:-1], vals[1:]):
        if a[1] < criterion <= b[1]:
            return a, b, vals
    return None, None, vals


def crossing_estimate(a, b, criterion: float) -> float:
    x0, y0 = a[0], math.log10(a[1])
    x1, y1 = b[0], math.log10(b[1])
    yc = math.log10(criterion)
    if abs(y1 - y0) < 1e-14:
        return 0.5 * (x0 + x1)
    x = x0 + (x1 - x0) * (yc - y0) / (y1 - y0)
    # Avoid pathological extrapolation and repeated endpoint evaluations.
    pad = 0.15 * (x1 - x0)
    return min(max(x, x0 + pad), x1 - pad)


def threshold_record(points: List[Dict], criterion: float, abs_tol: float, rel_tol: float):
    a, b, vals = locate_bracket(points, criterion)
    out = {
        "da_dN_criterion_m_per_cycle": criterion,
        "threshold_class": "unresolved",
        "DeltaK_threshold_lower_MPa_sqrtm": np.nan,
        "DeltaK_threshold_upper_MPa_sqrtm": np.nan,
        "DeltaK_threshold_estimate_MPa_sqrtm": np.nan,
        "threshold_width_MPa_sqrtm": np.nan,
        "n_evaluated_DeltaK": len(vals),
    }
    if a is not None:
        est = crossing_estimate(a, b, criterion)
        width = b[0] - a[0]
        out.update({
            "threshold_class": "bracketed",
            "DeltaK_threshold_lower_MPa_sqrtm": a[0],
            "DeltaK_threshold_upper_MPa_sqrtm": b[0],
            "DeltaK_threshold_estimate_MPa_sqrtm": est,
            "threshold_width_MPa_sqrtm": width,
            "lower_rate_m_per_cycle": a[1],
            "upper_rate_m_per_cycle": b[1],
            "lower_source": a[2],
            "upper_source": b[2],
            "converged": width <= max(abs_tol, rel_tol * max(est, 1e-12)),
        })
        return out
    if vals:
        if vals[0][1] >= criterion:
            out.update({
                "threshold_class": "below_search_range",
                "DeltaK_threshold_upper_MPa_sqrtm": vals[0][0],
                "upper_rate_m_per_cycle": vals[0][1],
                "upper_source": vals[0][2],
            })
        elif vals[-1][1] < criterion:
            out.update({
                "threshold_class": "above_search_range",
                "DeltaK_threshold_lower_MPa_sqrtm": vals[-1][0],
                "lower_rate_m_per_cycle": vals[-1][1],
                "lower_source": vals[-1][2],
            })
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--map-runner", default="")
    ap.add_argument("--case-table", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--case-filter", nargs="*", default=[])
    ap.add_argument("--temperatures", nargs="+", type=float, required=True)
    ap.add_argument("--emission-entropies-kB", nargs="+", type=float, required=True)
    ap.add_argument("--cleavage-entropies-kB", nargs="+", type=float, required=True)
    ap.add_argument("--criteria", nargs="+", type=float, default=[1e-12, 1e-10])
    ap.add_argument("--DeltaK-seeds", nargs="+", type=float,
                    default=[0.05, 0.10, 0.20, 0.40, 0.80, 1.60, 3.20, 6.40, 12.80])
    ap.add_argument("--DeltaK-min", type=float, default=0.025)
    ap.add_argument("--DeltaK-max", type=float, default=20.0)
    ap.add_argument("--threshold-abs-tol", type=float, default=0.05)
    ap.add_argument("--threshold-rel-tol", type=float, default=0.03)
    ap.add_argument("--max-refine-iters", type=int, default=10)
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--skip-monotonic", action="store_true")

    # Physics/controller settings copied from the corrected map defaults.
    ap.add_argument("--T-anchor-K", type=float, default=300.0)
    ap.add_argument("--R", type=float, default=0.1)
    ap.add_argument("--frequency-Hz", type=float, default=1000.0)
    ap.add_argument("--cycles-max", type=float, default=2e14)
    ap.add_argument("--max-blocks", type=int, default=10000)
    ap.add_argument("--max-block-cycles", type=float, default=float("inf"))
    ap.add_argument("--min-block-cycles", type=float, default=1e-10)
    ap.add_argument("--target-state-fraction", type=float, default=0.01)
    ap.add_argument("--saturation-tol-fraction", type=float, default=1e-4)
    ap.add_argument("--target-dN-store-unbounded", type=float, default=5.0)
    ap.add_argument("--n-advances", type=int, default=5)
    ap.add_argument("--da-m", type=float, default=20e-6)
    ap.add_argument("--n-phase", type=int, default=96)
    ap.add_argument("--target-dB", type=float, default=0.02)
    ap.add_argument("--target-dN-store", type=float, default=0.01)
    ap.add_argument("--target-dN-emit", type=float, default=0.20)
    ap.add_argument("--target-dN-mobile", type=float, default=0.20)
    ap.add_argument("--Kdot-MPa-sqrtm-per-s", type=float, default=0.005)
    ap.add_argument("--monotonic-Kmax-MPa", type=float, default=40.0)
    ap.add_argument("--monotonic-dK-MPa", type=float, default=0.025)
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    runner_path = choose_runner(args.map_runner)
    runner = load_module(runner_path)
    cases = runner.read_cases(Path(args.case_table))
    if args.case_filter:
        wanted = set(args.case_filter)
        cases = [c for c in cases if str(c["case_label"]) in wanted]
    if not cases:
        raise RuntimeError("no cases selected")

    ns = settings_namespace(args)
    fatigue_path = out / "fatigue_adaptive_points.csv"
    mono_path = out / "monotonic_adaptive_points.csv"
    fatigue_rows = pd.read_csv(fatigue_path).to_dict("records") if args.resume and fatigue_path.exists() else []
    mono_rows = pd.read_csv(mono_path).to_dict("records") if args.resume and mono_path.exists() else []
    fcache = {key_f(r["case_label"], r["S_emit_kB"], r["S_cleave_kB"], r["T_K"], r["Kmax_MPa_sqrtm"]): r for r in fatigue_rows}
    mcache = {key_m(r["case_label"], r["S_emit_kB"], r["S_cleave_kB"], r["T_K"]): r for r in mono_rows}

    def evaluate(case: Dict, Se: float, Sc: float, T: float, DeltaK: float) -> Dict:
        DeltaK = min(max(float(DeltaK), args.DeltaK_min), args.DeltaK_max)
        Kmax = DeltaK / max(1.0 - args.R, 1e-12)
        k = key_f(case["case_label"], Se, Sc, T, Kmax)
        if k in fcache:
            return fcache[k]
        rec = runner.run_fatigue(case, Se, Sc, T, Kmax, ns)
        rec.setdefault("DeltaK_MPa_sqrtm", DeltaK)
        rec["scenario"] = f"adaptive_Se{Se:g}_Sc{Sc:g}"
        fatigue_rows.append(rec)
        fcache[k] = rec
        save_df(fatigue_rows, fatigue_path)
        eff, src = effective_rate(rec)
        print(f"  F {case['case_label']} Se={Se:g} Sc={Sc:g} T={T:g} DK={DeltaK:.5g} rate={eff:.3e} ({src})")
        return rec

    thresholds: List[Dict] = []
    total_groups = len(cases) * len(args.emission_entropies_kB) * len(args.cleavage_entropies_kB) * len(args.temperatures)
    ig = 0
    for case in cases:
        for Se in args.emission_entropies_kB:
            for Sc in args.cleavage_entropies_kB:
                for T in args.temperatures:
                    ig += 1
                    print(f"\n=== group {ig}/{total_groups}: {case['case_label']} Se={Se:g} Sc={Sc:g} T={T:g} K ===")
                    if not args.skip_monotonic:
                        km = key_m(case["case_label"], Se, Sc, T)
                        if km not in mcache:
                            mrec = runner.run_monotonic(case, Se, Sc, T, ns)
                            mrec["scenario"] = f"adaptive_Se{Se:g}_Sc{Sc:g}"
                            mono_rows.append(mrec)
                            mcache[km] = mrec
                            save_df(mono_rows, mono_path)
                            print(f"  M Kc={mrec.get('Kc_first_MPa_sqrtm')}")

                    # Evaluate the coarse geometric seed grid once.  All criteria
                    # reuse these points, so sensitivity thresholds cost little extra.
                    group_points = []
                    for dk in sorted(set(args.DeltaK_seeds)):
                        if args.DeltaK_min <= dk <= args.DeltaK_max:
                            group_points.append(evaluate(case, Se, Sc, T, dk))

                    for crit in sorted(args.criteria):
                        for _ in range(args.max_refine_iters):
                            a, b, vals = locate_bracket(group_points, crit)
                            if a is None:
                                # Expand toward the missing side if the seed grid did
                                # not cover the crossing.
                                if vals and vals[0][1] >= crit and vals[0][0] > args.DeltaK_min * 1.001:
                                    new_dk = max(args.DeltaK_min, vals[0][0] / 2.0)
                                elif vals and vals[-1][1] < crit and vals[-1][0] < args.DeltaK_max / 1.001:
                                    new_dk = min(args.DeltaK_max, vals[-1][0] * 1.5)
                                else:
                                    break
                            else:
                                est = crossing_estimate(a, b, crit)
                                width = b[0] - a[0]
                                if width <= max(args.threshold_abs_tol, args.threshold_rel_tol * max(est, 1e-12)):
                                    break
                                new_dk = est
                            # Avoid duplicate evaluations due to floating-point
                            # roundoff or clipping at search limits.
                            existing = np.array([float(r["DeltaK_MPa_sqrtm"]) for r in group_points])
                            if len(existing) and np.min(np.abs(existing - new_dk)) < 1e-7:
                                break
                            group_points.append(evaluate(case, Se, Sc, T, new_dk))

                        tr = threshold_record(group_points, crit, args.threshold_abs_tol, args.threshold_rel_tol)
                        tr.update({
                            "case_label": case["case_label"],
                            "S_emit_kB": Se,
                            "S_cleave_kB": Sc,
                            "T_K": T,
                        })
                        thresholds.append(tr)

    thr = pd.DataFrame(thresholds)
    thr.to_csv(out / "adaptive_rate_thresholds.csv", index=False)
    mono_df = pd.DataFrame(mono_rows)
    if not mono_df.empty:
        link = thr.merge(mono_df, on=GROUP, how="left")
        link.to_csv(out / "adaptive_DBTT_fatigue_link.csv", index=False)
    with open(out / "adaptive_study_settings.json", "w") as f:
        json.dump({
            "map_runner": str(runner_path),
            "case_table": str(Path(args.case_table).resolve()),
            "temperatures": args.temperatures,
            "emission_entropies_kB": args.emission_entropies_kB,
            "cleavage_entropies_kB": args.cleavage_entropies_kB,
            "criteria": args.criteria,
            "DeltaK_seeds": args.DeltaK_seeds,
            "DeltaK_bounds": [args.DeltaK_min, args.DeltaK_max],
            "threshold_abs_tol": args.threshold_abs_tol,
            "threshold_rel_tol": args.threshold_rel_tol,
            "cycles_max": args.cycles_max,
            "max_blocks": args.max_blocks,
        }, f, indent=2)
    print(f"\nwrote {out / 'adaptive_rate_thresholds.csv'}")
    if (out / "adaptive_DBTT_fatigue_link.csv").exists():
        print(f"wrote {out / 'adaptive_DBTT_fatigue_link.csv'}")


if __name__ == "__main__":
    main()
