#!/usr/bin/env python3
"""Focused V1 EXP-floor search for a weakly temperature-dependent convex Kc(T) response.

The old weak-T response is used only to set the toughness scale and endpoint anchors.
The desired calibration response is a smooth, weakly temperature-dependent, convex
reference curve.  The optimization separately penalizes:
  * deviation from the smooth weak-T reference,
  * excessive total Kc range,
  * abrupt adjacent-temperature jumps,
  * violations of convexity (non-monotone slope), and
  * excessive roughness.

The search reuses both the original four-class broad screen and, when available,
the continuous global population from the expanded peak search.  Both cleavage
and emission EXP-floor parameters are perturbed continuously during refinement.
The selected expanded-search peak recommendation is copied into the output as a
locked/saved regime recommendation.
"""
from __future__ import annotations

import argparse
import json
import math
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

import run_v1_exp_floor_four_class_tuning as base
import run_v1_exp_floor_peak_expanded_search as expanded


def weak_reference(targets: pd.DataFrame, T: np.ndarray, tau_K: float) -> np.ndarray:
    """Smooth decreasing convex reference anchored to old weakT endpoints.

    K_ref(T) = K_inf + A exp(-(T-T0)/tau), with A and K_inf chosen so the
    300 and 1200 K values equal the old weakT response endpoints.  This keeps
    the prior toughness scale while deliberately removing the old 800-900 K
    discontinuity.
    """
    g = base.target_curve(targets, "weakT").sort_values("T_K")
    T = np.asarray(T, float)
    T0 = float(g.T_K.iloc[0]); T1 = float(g.T_K.iloc[-1])
    K0 = float(g.target_Kc_MPa_sqrtm.iloc[0])
    K1 = float(g.target_Kc_MPa_sqrtm.iloc[-1])
    tau = max(float(tau_K), 1.0)
    e = math.exp(-(T1 - T0) / tau)
    A = (K0 - K1) / max(1.0 - e, 1e-12)
    Kinf = K0 - A
    return Kinf + A * np.exp(-(T - T0) / tau)


def _curve_metrics_one(k: np.ndarray, T: np.ndarray, ref: np.ndarray,
                       Kmax: float, convex_tol_per_100K: float) -> dict[str, float]:
    T = np.asarray(T, float)
    k = np.asarray(k, float)
    cens = ~np.isfinite(k)
    p = np.where(cens, 1.15 * Kmax, k)

    # Whole-grid reference error retains the desired weak overall scale.
    ref_rmse = float(np.sqrt(np.mean(((p - ref) / np.maximum(ref, 2.0)) ** 2)))

    # Topology/convexity metrics are evaluated on the 100 K grid used by the
    # PF comparison.  This makes them robust to dK first-passage quantization
    # at the finer intermediate search grids.
    Ta = np.arange(300.0, 1200.1, 100.0)
    pa = np.interp(Ta, T, p)
    ra = np.interp(Ta, T, ref)
    meanK = max(float(np.mean(pa)), 1.0)
    range_rel = float((np.max(pa) - np.min(pa)) / meanK)
    std_rel = float(np.std(pa) / meanK)

    dp = np.diff(pa)
    jump_rel = float(np.max(np.abs(dp)) / meanK) if len(dp) else 0.0
    jump_rms_rel = float(np.sqrt(np.mean(dp**2)) / meanK) if len(dp) else 0.0

    # Fractional slope per 100 K.  Convexity requires the slope sequence to
    # be nondecreasing.  Positive slope is also softly penalized because the
    # desired weak-T reference is gently decreasing and convex.
    slope = np.diff(pa) / meanK
    ds = np.diff(slope)
    tol = max(float(convex_tol_per_100K), 0.0)
    viol = np.maximum(-(ds) - tol, 0.0)
    convex_rms = float(np.sqrt(np.mean(viol**2))) if len(viol) else 0.0
    convex_max = float(np.max(viol)) if len(viol) else 0.0
    convex_fraction = float(np.mean(viol > 0.0)) if len(viol) else 0.0
    positive = np.maximum(slope - tol, 0.0)
    monotone_violation_rms = float(np.sqrt(np.mean(positive**2))) if len(positive) else 0.0
    monotone_violation_fraction = float(np.mean(positive > 0.0)) if len(positive) else 0.0
    roughness = float(np.sqrt(np.mean(ds**2))) if len(ds) else 0.0

    anchor = 0.5 * (
        abs(float(pa[0] - ra[0])) / max(float(ra[0]), 2.0)
        + abs(float(pa[-1] - ra[-1])) / max(float(ra[-1]), 2.0)
    )
    cpen = float(np.mean(cens))

    # A total range <=12% of the mean is treated as already weakly
    # temperature dependent; only excess range is penalized.  Convexity is
    # enforced more strongly through both magnitude and violation fraction,
    # so a nearly flat but non-convex curve cannot beat a smooth convex one.
    excess_range = max(range_rel - 0.12, 0.0)
    score = (
        0.45 * ref_rmse
        + 0.10 * excess_range
        + 0.06 * std_rel
        + 0.10 * jump_rel
        + 1.25 * convex_rms
        + 0.08 * convex_fraction
        + 0.35 * monotone_violation_rms
        + 0.03 * monotone_violation_fraction
        + 0.02 * roughness
        + 0.05 * anchor
        + 2.0 * cpen
    )
    return {
        "weak_convex_score": score,
        "reference_rmse": ref_rmse,
        "range_rel": range_rel,
        "std_rel": std_rel,
        "max_jump_rel": jump_rel,
        "jump_rms_rel": jump_rms_rel,
        "convex_violation_rms": convex_rms,
        "convex_violation_max": convex_max,
        "convex_violation_fraction": convex_fraction,
        "monotone_violation_rms": monotone_violation_rms,
        "monotone_violation_fraction": monotone_violation_fraction,
        "roughness": roughness,
        "anchor_penalty": anchor,
        "K_mean": float(np.mean(pa)),
        "K_min": float(np.min(pa)),
        "K_max": float(np.max(pa)),
        "n_censored": int(np.sum(cens)),
    }

def weak_scores(Kc: np.ndarray, T: np.ndarray, targets: pd.DataFrame,
                Kmax: float, tau_K: float, convex_tol_per_100K: float) -> pd.DataFrame:
    ref = weak_reference(targets, np.asarray(T, float), tau_K)
    rows = [_curve_metrics_one(k, T, ref, Kmax, convex_tol_per_100K) for k in Kc]
    return pd.DataFrame(rows)


def add_scores(cand: pd.DataFrame, K: np.ndarray, T: np.ndarray, targets: pd.DataFrame,
               Kmax: float, tau_K: float, convex_tol_per_100K: float) -> pd.DataFrame:
    s = weak_scores(K, T, targets, Kmax, tau_K, convex_tol_per_100K)
    return pd.concat([cand.reset_index(drop=True), s.reset_index(drop=True)], axis=1)


def run_stage(cand: pd.DataFrame, T: np.ndarray, dK: float, Kdot: float, Kmax: float,
              targets: pd.DataFrame, out_npz: Path, out_csv: Path, resume: bool,
              tau_K: float, convex_tol_per_100K: float) -> tuple[np.ndarray, pd.DataFrame]:
    if resume and out_npz.exists() and out_csv.exists():
        z = np.load(out_npz)
        return z["Kc"], pd.read_csv(out_csv)
    from arrhenius_fracture.config import ElasticProperties
    mat = ElasticProperties()
    K = base.simulate_candidates(cand, T, Kmax=Kmax, dK=dK, Kdot=Kdot,
                                 G_Pa=mat.G, nu=mat.nu, b_m=mat.b)
    s = add_scores(cand, K, T, targets, Kmax, tau_K, convex_tol_per_100K)
    np.savez_compressed(out_npz, Kc=K, T=T)
    s.to_csv(out_csv, index=False)
    return K, s


def _param_subset(df: pd.DataFrame) -> list[str]:
    return [c for c in expanded.CLEAVE_KEYS + expanded.EMIT_KEYS if c in df.columns]


def diverse_seed_take(scored: pd.DataFrame, n: int) -> pd.DataFrame:
    if n <= 0 or scored.empty:
        return scored.iloc[:0].copy()
    s = scored[scored.n_censored == 0].copy()
    if s.empty:
        s = scored.copy()

    q = max(2, n // 4)
    parts = [
        s.nsmallest(q, "weak_convex_score"),
        s.nsmallest(q, "range_rel"),
        s.nsmallest(q, "convex_violation_rms"),
        s.nsmallest(q, "reference_rmse"),
    ]
    out = pd.concat(parts, ignore_index=True, sort=False)
    subset = _param_subset(out)
    if subset:
        out = out.drop_duplicates(subset=subset)
    # Fill any remaining slots by total score while limiting collapse onto a
    # single original surface/context where those metadata exist.
    selected = [r for _, r in out.sort_values("weak_convex_score").iterrows()]
    used_ids = {int(r.get("candidate_id", -10**9)) for r in selected}
    surf_count: dict[str, int] = {}
    ctx_count: dict[str, int] = {}
    for r in selected:
        surf = str(r.get("surface_id", "")); ctx = str(r.get("context_id", ""))
        surf_count[surf] = surf_count.get(surf, 0) + 1
        ctx_count[ctx] = ctx_count.get(ctx, 0) + 1
    for _, r in s.sort_values("weak_convex_score").iterrows():
        if len(selected) >= n:
            break
        cid = int(r.get("candidate_id", -1))
        if cid in used_ids:
            continue
        surf = str(r.get("surface_id", "")); ctx = str(r.get("context_id", ""))
        if surf and surf_count.get(surf, 0) >= max(3, n // 4):
            continue
        if ctx and ctx_count.get(ctx, 0) >= 3:
            continue
        selected.append(r)
        used_ids.add(cid)
        surf_count[surf] = surf_count.get(surf, 0) + 1
        ctx_count[ctx] = ctx_count.get(ctx, 0) + 1
    return pd.DataFrame(selected).head(n).reset_index(drop=True)


def combined_seed_pool(source_scored: pd.DataFrame, global_scored: pd.DataFrame | None,
                       source_run: Path, n_seeds: int) -> pd.DataFrame:
    n_source = n_seeds if global_scored is None else max(1, n_seeds // 2)
    pieces = [diverse_seed_take(source_scored, n_source)]
    if global_scored is not None:
        pieces.append(diverse_seed_take(global_scored, n_seeds - n_source))

    # Always preserve the original weakT recommendation as an explicit seed.
    rec = source_run / "recommended_exp_floor_four_class.csv"
    if rec.exists():
        r = pd.read_csv(rec)
        rw = r[r.target_class.astype(str) == "weakT"].copy()
        if not rw.empty:
            pieces.append(rw)

    seeds = pd.concat(pieces, ignore_index=True, sort=False)
    subset = _param_subset(seeds)
    if subset:
        seeds = seeds.drop_duplicates(subset=subset)
    return seeds.head(n_seeds).reset_index(drop=True)


def select_generation_seeds(scored: pd.DataFrame, n: int) -> pd.DataFrame:
    s = scored[scored.n_censored == 0].copy()
    if s.empty:
        s = scored.copy()
    q = max(1, n // 4)
    parts = [
        s.nsmallest(max(q, n // 2), "weak_convex_score"),
        s.nsmallest(q, "reference_rmse"),
        s.nsmallest(q, "convex_violation_rms"),
        s.nsmallest(q, "range_rel"),
    ]
    out = pd.concat(parts, ignore_index=True, sort=False)
    subset = _param_subset(out)
    if subset:
        out = out.drop_duplicates(subset=subset)
    return out.sort_values("weak_convex_score").head(n).reset_index(drop=True)


def make_plot(curve: pd.DataFrame, targets: pd.DataFrame, tau_K: float, out: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    old = base.target_curve(targets, "weakT").sort_values("T_K")
    T = curve.T_K.to_numpy(float)
    ref = weak_reference(targets, T, tau_K)
    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    ax.plot(old.T_K, old.target_Kc_MPa_sqrtm, label="prior weakT response")
    ax.plot(T, ref, label="smooth convex weak-T reference")
    ax.plot(T, curve.Kc_MPa_sqrtm, label="tuned EXP-floor V1")
    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel(r"$K_c$ (MPa$\sqrt{m}$)")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / "weakT_convex_fit.png", dpi=220)
    plt.close(fig)


def save_locked_peak(peak_run: Path, out: Path) -> None:
    here = Path(__file__).resolve().parent
    items = [
        ("recommended_peak_expanded.csv", "locked_peak_expanded_recommendation.csv", "saved_peak_recommendation.csv"),
        ("recommended_peak_curve_dense.csv", "locked_peak_expanded_curve_dense.csv", "saved_peak_curve_dense.csv"),
    ]
    for run_name, fallback_name, dst_name in items:
        src = peak_run / run_name
        if not src.exists():
            src = here / fallback_name
        if src.exists():
            shutil.copy2(src, out / dst_name)


def build_regime_manifest(source_run: Path, peak_run: Path, weak_row: pd.Series,
                          out: Path) -> None:
    rec = source_run / "recommended_exp_floor_four_class.csv"
    if not rec.exists():
        return
    base_rec = pd.read_csv(rec)
    rows: list[pd.DataFrame] = []
    rows.append(base_rec[base_rec.target_class.isin(["ceramic", "DBTT"])].copy())

    pfile = peak_run / "recommended_peak_expanded.csv"
    if pfile.exists():
        p = pd.read_csv(pfile).copy()
        p.insert(0, "target_class", "peak")
        p.insert(1, "recommendation_source", "expanded_peak_search_locked")
        rows.append(p)
    w = weak_row.to_frame().T.copy()
    w.insert(0, "target_class", "weakT")
    w.insert(1, "recommendation_source", "weakT_convex_refinement")
    rows.append(w)
    pd.concat(rows, ignore_index=True, sort=False).to_csv(out / "saved_regime_recommendations.csv", index=False)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source-run", default="runs/v1_exp_floor_four_class_tuning")
    ap.add_argument("--peak-run", default="runs/v1_exp_floor_peak_expanded_search")
    ap.add_argument("--targets", default="exp_floor_four_class_target_curves.csv")
    ap.add_argument("--out", default="runs/v1_exp_floor_weakT_convex_refinement")
    ap.add_argument("--Kdot", type=float, default=0.005)
    ap.add_argument("--Kmax", type=float, default=80.0)
    ap.add_argument("--tau-K", type=float, default=450.0,
                    help="Decay scale of the smooth convex weak-T reference [K]")
    ap.add_argument("--convex-tol", type=float, default=0.002,
                    help="Tolerance in fractional slope change per 100 K")
    ap.add_argument("--n-seeds", type=int, default=64)
    ap.add_argument("--gen1-perturb", type=int, default=64)
    ap.add_argument("--gen2-seeds", type=int, default=36)
    ap.add_argument("--gen2-perturb", type=int, default=72)
    ap.add_argument("--gen3-seeds", type=int, default=20)
    ap.add_argument("--gen3-perturb", type=int, default=64)
    ap.add_argument("--finalists", type=int, default=12)
    ap.add_argument("--gen1-dK", type=float, default=0.075)
    ap.add_argument("--gen2-dK", type=float, default=0.04)
    ap.add_argument("--gen3-dK", type=float, default=0.025)
    ap.add_argument("--final-dK", type=float, default=0.02)
    ap.add_argument("--seed", type=int, default=20260713)
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()

    if args.smoke:
        args.n_seeds = min(args.n_seeds, 8)
        args.gen1_perturb = min(args.gen1_perturb, 3)
        args.gen2_seeds = min(args.gen2_seeds, 6)
        args.gen2_perturb = min(args.gen2_perturb, 3)
        args.gen3_seeds = min(args.gen3_seeds, 4)
        args.gen3_perturb = min(args.gen3_perturb, 2)
        args.finalists = min(args.finalists, 3)
        args.gen1_dK = max(args.gen1_dK, 0.5)
        args.gen2_dK = max(args.gen2_dK, 0.25)
        args.gen3_dK = max(args.gen3_dK, 0.15)
        args.final_dK = max(args.final_dK, 0.1)

    src = Path(args.source_run)
    peak = Path(args.peak_run)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    targets = base.load_targets(args.targets)
    save_locked_peak(peak, out)

    # Rescore the original broad population for the weak-convex objective.
    broad = pd.read_csv(src / "broad_candidates.csv")
    bz = np.load(src / "broad_Kc.npz")
    sb = add_scores(broad, bz["Kc"], bz["T"], targets, args.Kmax,
                    args.tau_K, args.convex_tol)
    sb.to_csv(out / "source_broad_weak_convex_scores.csv", index=False)

    # Reuse the continuous global population generated during the expanded
    # peak search.  This gives broad continuous coverage of both barrier
    # channels at essentially no additional global-search cost.
    sg = None
    if (peak / "global_candidates.csv").exists() and (peak / "global_Kc.npz").exists():
        gc = pd.read_csv(peak / "global_candidates.csv")
        gz = np.load(peak / "global_Kc.npz")
        sg = add_scores(gc, gz["Kc"], gz["T"], targets, args.Kmax,
                        args.tau_K, args.convex_tol)
        sg.to_csv(out / "source_continuous_global_weak_convex_scores.csv", index=False)

    seeds = combined_seed_pool(sb, sg, src, args.n_seeds)
    seeds.to_csv(out / "multibasin_weakT_seeds.csv", index=False)
    print(f"selected {len(seeds)} weak-convex seeds", flush=True)

    T1 = np.arange(300.0, 1200.1, 50.0)
    T2 = np.arange(300.0, 1200.1, 25.0)
    T3 = np.arange(300.0, 1200.1, 10.0)
    Td = base.make_temperature_grid(targets, "dense")

    g1 = expanded.generate_population(seeds, args.gen1_perturb, args.seed,
                                      shrink=0.80, allow_surface_mutation=True)
    g1.to_csv(out / "generation1_candidates.csv", index=False)
    print(f"generation 1: {len(g1)} candidates", flush=True)
    K1, s1 = run_stage(g1, T1, args.gen1_dK, args.Kdot, args.Kmax, targets,
                       out / "generation1_Kc.npz", out / "generation1_scores.csv",
                       args.resume, args.tau_K, args.convex_tol)

    seed2 = select_generation_seeds(s1, args.gen2_seeds)
    seed2.to_csv(out / "generation2_seeds.csv", index=False)
    g2 = expanded.generate_population(seed2, args.gen2_perturb, args.seed + 1,
                                      shrink=0.42, allow_surface_mutation=True)
    g2.to_csv(out / "generation2_candidates.csv", index=False)
    print(f"generation 2: {len(g2)} candidates", flush=True)
    K2, s2 = run_stage(g2, T2, args.gen2_dK, args.Kdot, args.Kmax, targets,
                       out / "generation2_Kc.npz", out / "generation2_scores.csv",
                       args.resume, args.tau_K, args.convex_tol)

    seed3 = select_generation_seeds(s2, args.gen3_seeds)
    seed3.to_csv(out / "generation3_seeds.csv", index=False)
    g3 = expanded.generate_population(seed3, args.gen3_perturb, args.seed + 2,
                                      shrink=0.20, allow_surface_mutation=True)
    g3.to_csv(out / "generation3_candidates.csv", index=False)
    print(f"generation 3: {len(g3)} candidates", flush=True)
    K3, s3 = run_stage(g3, T3, args.gen3_dK, args.Kdot, args.Kmax, targets,
                       out / "generation3_Kc.npz", out / "generation3_scores.csv",
                       args.resume, args.tau_K, args.convex_tol)

    finalists = select_generation_seeds(s3, args.finalists)
    ids = set(finalists.candidate_id.astype(int))
    fc = g3[g3.candidate_id.astype(int).isin(ids)].copy().reset_index(drop=True)
    print(f"final: {len(fc)} candidates, {len(Td)} temperatures", flush=True)
    Kf, sf = run_stage(fc, Td, args.final_dK, args.Kdot, args.Kmax, targets,
                       out / "final_Kc.npz", out / "final_scores.csv",
                       args.resume, args.tau_K, args.convex_tol)

    best = sf.sort_values("weak_convex_score").iloc[0]
    bid = int(best.candidate_id)
    j = int(np.flatnonzero(fc.candidate_id.to_numpy(int) == bid)[0])
    best.to_frame().T.to_csv(out / "recommended_weakT_convex.csv", index=False)
    curve = pd.DataFrame({
        "T_K": Td,
        "Kc_MPa_sqrtm": Kf[j],
        "smooth_convex_reference_MPa_sqrtm": weak_reference(targets, Td, args.tau_K),
        "prior_weakT_Kc_MPa_sqrtm": base.interp_target(targets, "weakT", Td),
    })
    curve.to_csv(out / "recommended_weakT_curve_dense.csv", index=False)
    make_plot(curve, targets, args.tau_K, out)
    build_regime_manifest(src, peak, best, out)

    cfg = vars(args).copy()
    old = base.target_curve(targets, "weakT").sort_values("T_K")
    cfg["reference_definition"] = {
        "form": "Kinf + A*exp(-(T-T0)/tau_K)",
        "T0_K": float(old.T_K.iloc[0]),
        "T1_K": float(old.T_K.iloc[-1]),
        "K0_MPa_sqrtm": float(old.target_Kc_MPa_sqrtm.iloc[0]),
        "K1_MPa_sqrtm": float(old.target_Kc_MPa_sqrtm.iloc[-1]),
        "tau_K": args.tau_K,
    }
    cfg["objective_notes"] = [
        "old weakT interior discontinuity is not a fitting target",
        "old weakT endpoints set the toughness scale",
        "smooth weak-dependence and convexity are explicit optimization terms",
        "expanded-search peak result is copied and locked, not re-optimized",
    ]
    with (out / "run_config.json").open("w") as f:
        json.dump(cfg, f, indent=2)

    cols = [c for c in [
        "candidate_id", "surface_id", "cleave_G00_eV", "cleave_gT_eV_per_K",
        "cleave_sigc0_GPa", "cleave_sT_MPa_per_K", "cleave_exp_a",
        "cleave_exp_n", "cleave_floor_frac", "cleave_S_hs_kB",
        "chi_shield", "N_sat", "exp_G00_eV", "exp_sigc0_GPa", "exp_a",
        "exp_n", "exp_floor_frac", "exp_gT_eV_per_K", "exp_sT_MPa_per_K",
        "weak_convex_score", "reference_rmse", "range_rel", "max_jump_rel",
        "convex_violation_rms", "K_mean", "K_min", "K_max",
    ] if c in best.index]
    print("\nBest weak-T convex candidate:")
    print(best[cols].to_frame().T.to_string(index=False))
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
