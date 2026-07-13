#!/usr/bin/env python3
"""Focused V1 refinement for the narrow peak-like fracture regime.

This script reuses a completed four-class broad screen, but *does not* select
peak seeds by the generic four-class RMS score.  Instead it identifies actual
interior-peak phenotypes, refines them in several shrinking local generations,
and finally evaluates rate sensitivity around the calibration rate.

The prior 1-D Kc(T) peak curve is the observable target.  Kdot is held fixed at
the target value during calibration and is varied only in the post-fit rate
sensitivity diagnostic.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

import run_v1_exp_floor_four_class_tuning as base


def _interp_target(targets: pd.DataFrame, T: np.ndarray) -> np.ndarray:
    return base.interp_target(targets, "peak", np.asarray(T, float))


def _nearest_idx(T: np.ndarray, x: float) -> int:
    return int(np.argmin(np.abs(np.asarray(T, float) - float(x))))


def target_features(targets: pd.DataFrame) -> dict[str, float]:
    g = base.target_curve(targets, "peak").sort_values("T_K")
    T = g.T_K.to_numpy(float)
    K = g.target_Kc_MPa_sqrtm.to_numpy(float)
    mask = (T >= 820.0) & (T <= 980.0)
    ii = np.where(mask)[0]
    ip = ii[int(np.argmax(K[ii]))]
    Tpk = float(T[ip]); Kpk = float(K[ip])
    Tpre = Tpk - 55.0
    Tpost = Tpk + 45.0
    Kpre = float(np.interp(Tpre, T, K))
    Kpost = float(np.interp(Tpost, T, K))
    return {
        "T_peak_target_K": Tpk,
        "K_peak_target": Kpk,
        "T_pre_K": Tpre,
        "T_post_K": Tpost,
        "K_pre_target": Kpre,
        "K_post_target": Kpost,
        "rise_target": Kpk - Kpre,
        "fall_target": Kpk - Kpost,
    }


def peak_scores(Kc: np.ndarray, Tgrid: np.ndarray, targets: pd.DataFrame,
                min_prominence: float = 1.0) -> pd.DataFrame:
    """Phenotype-aware score for narrow peak-like curves."""
    T = np.asarray(Tgrid, float)
    K = np.asarray(Kc, float)
    target = _interp_target(targets, T)
    feat = target_features(targets)

    pred = np.where(np.isfinite(K), K, 1.15 * np.nanmax(np.where(np.isfinite(K), K, 0.0), axis=1, keepdims=True) + 10.0)
    # Weight the narrow transition zone heavily while retaining global scale.
    w = np.ones_like(T, float)
    w[(T >= 825.0) & (T <= 975.0)] = 5.0
    w[(T >= feat["T_peak_target_K"] - 20.0) & (T <= feat["T_peak_target_K"] + 20.0)] = 12.0
    w /= np.sum(w)

    scale = np.maximum(target, 2.0)
    rel = np.sqrt(np.sum(w[None, :] * ((pred - target[None, :]) / scale[None, :]) ** 2, axis=1))
    log = np.sqrt(np.sum(w[None, :] * (np.log10(pred + 0.25) - np.log10(target[None, :] + 0.25)) ** 2, axis=1))

    win = np.where((T >= 820.0) & (T <= 980.0))[0]
    local_arg = np.argmax(pred[:, win], axis=1)
    peak_idx = win[local_arg]
    row = np.arange(len(pred))
    Tpk = T[peak_idx]
    Kpk = pred[row, peak_idx]
    i_pre = _nearest_idx(T, feat["T_pre_K"])
    i_post = _nearest_idx(T, feat["T_post_K"])
    Kpre = pred[:, i_pre]
    Kpost = pred[:, i_post]
    rise = Kpk - Kpre
    fall = Kpk - Kpost
    prominence = np.minimum(rise, fall)

    peak_T_pen = np.abs(Tpk - feat["T_peak_target_K"]) / 25.0
    peak_K_pen = np.abs(Kpk - feat["K_peak_target"]) / max(feat["K_peak_target"], 1.0)
    pre_pen = np.abs(Kpre - feat["K_pre_target"]) / max(feat["K_pre_target"], 2.0)
    post_pen = np.abs(Kpost - feat["K_post_target"]) / max(feat["K_post_target"], 2.0)
    rise_pen = np.abs(rise - feat["rise_target"]) / max(feat["rise_target"], 1.0)
    fall_pen = np.abs(fall - feat["fall_target"]) / max(feat["fall_target"], 1.0)

    topo_pen = np.zeros(len(pred), float)
    topo_pen[prominence < min_prominence] += 2.0
    topo_pen[rise <= 0.0] += 2.0
    topo_pen[fall <= 0.0] += 2.0

    score = (
        0.30 * rel + 0.12 * log + 0.14 * peak_T_pen + 0.10 * peak_K_pen
        + 0.07 * pre_pen + 0.07 * post_pen + 0.08 * rise_pen + 0.08 * fall_pen
        + topo_pen
    )
    return pd.DataFrame({
        "peak_score": score,
        "weighted_rel_rmse": rel,
        "weighted_log_rmse": log,
        "peak_T_K": Tpk,
        "peak_Kc": Kpk,
        "K_pre": Kpre,
        "K_post": Kpost,
        "peak_rise": rise,
        "peak_fall": fall,
        "peak_prominence": prominence,
        "peak_valid": (prominence >= min_prominence) & (rise > 0.0) & (fall > 0.0),
    })


def add_peak_scores(cand: pd.DataFrame, Kc: np.ndarray, T: np.ndarray,
                    targets: pd.DataFrame, min_prominence: float) -> pd.DataFrame:
    q = cand.copy().reset_index(drop=True)
    f = peak_scores(Kc, T, targets, min_prominence=min_prominence)
    return pd.concat([q, f], axis=1)


def choose_diverse_seeds(scored: pd.DataFrame, n: int) -> pd.DataFrame:
    """Retain real peaks while preventing one context/surface from monopolizing seeds."""
    v = scored[scored.peak_valid].sort_values("peak_score").copy()
    if v.empty:
        raise RuntimeError("No actual interior-peak candidates found in source broad screen")
    rows = []
    ctx_count: dict[str, int] = {}
    surf_count: dict[str, int] = {}
    for _, r in v.iterrows():
        ctx = str(r.get("context_id", "")); surf = str(r.get("surface_id", ""))
        if ctx_count.get(ctx, 0) >= 3:
            continue
        if surf_count.get(surf, 0) >= max(6, n // 3):
            continue
        rows.append(r)
        ctx_count[ctx] = ctx_count.get(ctx, 0) + 1
        surf_count[surf] = surf_count.get(surf, 0) + 1
        if len(rows) >= n:
            break
    if len(rows) < min(n, len(v)):
        used = {int(r.candidate_id) for r in rows}
        for _, r in v.iterrows():
            if int(r.candidate_id) in used:
                continue
            rows.append(r)
            if len(rows) >= n:
                break
    return pd.DataFrame(rows).reset_index(drop=True)


def _clip_scalar(x: float, name: str) -> float:
    return float(base._clip(np.array([x], float), name)[0])


def perturb_peak_candidates(seeds: pd.DataFrame, n_perturb: int, seed: int,
                            shrink: float = 1.0) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows: list[dict] = []
    for _, r in seeds.iterrows():
        b = r.to_dict()
        # Remove score-only columns while preserving model parameters and metadata.
        drop = {c for c in b if c.startswith("peak_") or c.startswith("weighted_") or c in {"K_pre", "K_post"}}
        b = {k: v for k, v in b.items() if k not in drop}
        rows.append(dict(b))
        for _ in range(n_perturb):
            x = dict(b)
            x["cleave_G00_eV"] = _clip_scalar(b["cleave_G00_eV"] * math.exp(rng.normal(0, 0.08 * shrink)), "cleave_G00_eV")
            x["cleave_gT_eV_per_K"] = _clip_scalar(b["cleave_gT_eV_per_K"] + rng.normal(0, 0.00045 * shrink), "cleave_gT_eV_per_K")
            x["cleave_sigc0_GPa"] = _clip_scalar(b["cleave_sigc0_GPa"] * math.exp(rng.normal(0, 0.08 * shrink)), "cleave_sigc0_GPa")
            x["cleave_sT_MPa_per_K"] = _clip_scalar(b["cleave_sT_MPa_per_K"] + rng.normal(0, 0.35 * shrink), "cleave_sT_MPa_per_K")
            x["cleave_exp_a"] = _clip_scalar(b["cleave_exp_a"] * math.exp(rng.normal(0, 0.10 * shrink)), "cleave_exp_a")
            x["cleave_exp_n"] = _clip_scalar(b["cleave_exp_n"] * math.exp(rng.normal(0, 0.08 * shrink)), "cleave_exp_n")
            x["cleave_floor_frac"] = _clip_scalar(b["cleave_floor_frac"] * math.exp(rng.normal(0, 0.16 * shrink)), "cleave_floor_frac")
            x["cleave_S_hs_kB"] = _clip_scalar(b["cleave_S_hs_kB"] + rng.normal(0, 3.0 * shrink), "cleave_S_hs_kB")
            x["chi_shield"] = _clip_scalar(b["chi_shield"] + rng.normal(0, 0.025 * shrink), "chi_shield")
            if math.isfinite(float(b["N_sat"])):
                x["N_sat"] = float(np.clip(float(b["N_sat"]) * math.exp(rng.normal(0, 0.16 * shrink)), *base.RANGES["N_sat_finite"]))
            rows.append(x)
    out = pd.DataFrame(rows)
    # Stage-local identity; preserve parent ID.
    if "candidate_id" in out.columns:
        out = out.rename(columns={"candidate_id": "parent_candidate_id"})
    out.insert(0, "candidate_id", np.arange(len(out), dtype=int))
    return out


def focused_grid() -> np.ndarray:
    return np.array(sorted(set([
        300, 400, 500, 600, 700, 750, 800, 825, 840, 850, 860, 870, 875, 880,
        885, 890, 895, 900, 905, 910, 915, 920, 925, 930, 935, 940, 945, 950,
        955, 960, 970, 980, 1000, 1050, 1100, 1200
    ])), float)


def save_plot(curve: pd.DataFrame, rate: pd.DataFrame, targets: pd.DataFrame, out: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    tg = base.target_curve(targets, "peak").sort_values("T_K")
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    ax.plot(tg.T_K, tg.target_Kc_MPa_sqrtm, label="prior 1-D target")
    ax.plot(curve.T_K, curve.Kc_MPa_sqrtm, label="focused EXP-floor fit")
    ax.set_xlabel("Temperature (K)"); ax.set_ylabel(r"$K_c$ (MPa$\sqrt{m}$)")
    ax.grid(alpha=0.25); ax.legend(); fig.tight_layout()
    fig.savefig(out / "peak_focused_fit.png", dpi=220); plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for kd, g in rate.groupby("Kdot"):
        ax.plot(g.T_K, g.Kc_MPa_sqrtm, marker="o", ms=3, label=f"Kdot={kd:g}")
    ax.set_xlabel("Temperature (K)"); ax.set_ylabel(r"$K_c$ (MPa$\sqrt{m}$)")
    ax.grid(alpha=0.25); ax.legend(fontsize=8); fig.tight_layout()
    fig.savefig(out / "peak_rate_sensitivity.png", dpi=220); plt.close(fig)


def run_stage(cand: pd.DataFrame, T: np.ndarray, dK: float, Kdot: float,
              Kmax: float, targets: pd.DataFrame, out_npz: Path, out_csv: Path,
              min_prominence: float, resume: bool) -> tuple[np.ndarray, pd.DataFrame]:
    if resume and out_npz.exists() and out_csv.exists():
        z = np.load(out_npz); K = z["Kc"]; score = pd.read_csv(out_csv)
        return K, score
    from arrhenius_fracture.config import ElasticProperties
    mat = ElasticProperties()
    K = base.simulate_candidates(cand, T, Kmax=Kmax, dK=dK, Kdot=Kdot,
                                 G_Pa=mat.G, nu=mat.nu, b_m=mat.b)
    score = add_peak_scores(cand, K, T, targets, min_prominence)
    np.savez_compressed(out_npz, Kc=K, T=T)
    score.to_csv(out_csv, index=False)
    return K, score


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source-run", default="runs/v1_exp_floor_four_class_tuning")
    ap.add_argument("--targets", default="exp_floor_four_class_target_curves.csv")
    ap.add_argument("--out", default="runs/v1_exp_floor_peak_focused_refinement")
    ap.add_argument("--Kdot", type=float, default=0.005)
    ap.add_argument("--Kmax", type=float, default=80.0)
    ap.add_argument("--n-seeds", type=int, default=24)
    ap.add_argument("--gen1-perturb", type=int, default=48)
    ap.add_argument("--gen2-seeds", type=int, default=24)
    ap.add_argument("--gen2-perturb", type=int, default=64)
    ap.add_argument("--finalists", type=int, default=8)
    ap.add_argument("--gen1-dK", type=float, default=0.05)
    ap.add_argument("--gen2-dK", type=float, default=0.025)
    ap.add_argument("--final-dK", type=float, default=0.02)
    ap.add_argument("--min-prominence", type=float, default=1.0)
    ap.add_argument("--rates", default="0.00125 0.0025 0.005 0.01 0.02")
    ap.add_argument("--seed", type=int, default=20260711)
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()

    if args.smoke:
        args.n_seeds = min(args.n_seeds, 4)
        args.gen1_perturb = min(args.gen1_perturb, 3)
        args.gen2_seeds = min(args.gen2_seeds, 4)
        args.gen2_perturb = min(args.gen2_perturb, 3)
        args.finalists = min(args.finalists, 2)
        args.gen1_dK = max(args.gen1_dK, 0.5)
        args.gen2_dK = max(args.gen2_dK, 0.25)
        args.final_dK = max(args.final_dK, 0.1)

    src = Path(args.source_run); out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    targets = base.load_targets(args.targets)
    broad = pd.read_csv(src / "broad_candidates.csv")
    z = np.load(src / "broad_Kc.npz"); Kb = z["Kc"]; Tb = z["T"]
    broad_scored = add_peak_scores(broad, Kb, Tb, targets, args.min_prominence)
    broad_scored.to_csv(out / "source_broad_peak_scores.csv", index=False)
    seeds = choose_diverse_seeds(broad_scored, args.n_seeds)
    seeds.to_csv(out / "peak_seed_candidates.csv", index=False)
    print(f"selected {len(seeds)} actual-peak seeds from broad screen; best source peak score={seeds.peak_score.min():.4g}", flush=True)

    Tfocus = focused_grid()
    gen1 = perturb_peak_candidates(seeds, args.gen1_perturb, args.seed, shrink=1.0)
    gen1.to_csv(out / "generation1_candidates.csv", index=False)
    print(f"generation 1: {len(gen1)} candidates, {len(Tfocus)} temperatures, dK={args.gen1_dK}", flush=True)
    K1, s1 = run_stage(gen1, Tfocus, args.gen1_dK, args.Kdot, args.Kmax, targets,
                       out / "generation1_Kc.npz", out / "generation1_scores.csv",
                       args.min_prominence, args.resume)

    seed2 = s1[s1.peak_valid].nsmallest(args.gen2_seeds, "peak_score")
    if seed2.empty:
        seed2 = s1.nsmallest(args.gen2_seeds, "peak_score")
    gen2 = perturb_peak_candidates(seed2, args.gen2_perturb, args.seed + 1, shrink=0.50)
    gen2.to_csv(out / "generation2_candidates.csv", index=False)
    print(f"generation 2: {len(gen2)} candidates, dK={args.gen2_dK}", flush=True)
    K2, s2 = run_stage(gen2, Tfocus, args.gen2_dK, args.Kdot, args.Kmax, targets,
                       out / "generation2_Kc.npz", out / "generation2_scores.csv",
                       args.min_prominence, args.resume)

    finalists = s2[s2.peak_valid].nsmallest(args.finalists, "peak_score")
    if finalists.empty:
        finalists = s2.nsmallest(args.finalists, "peak_score")
    # Strip score columns and keep candidate model fields.
    model_ids = finalists.candidate_id.astype(int).to_numpy()
    final_cand = gen2[gen2.candidate_id.isin(model_ids)].copy().reset_index(drop=True)
    Td = base.make_temperature_grid(targets, "dense")
    print(f"final: {len(final_cand)} candidates, dense {len(Td)} temperatures, dK={args.final_dK}", flush=True)
    Kf, sf = run_stage(final_cand, Td, args.final_dK, args.Kdot, args.Kmax, targets,
                       out / "final_Kc.npz", out / "final_scores.csv",
                       args.min_prominence, args.resume)
    best = sf.sort_values("peak_score").iloc[0]
    cid = int(best.candidate_id)
    j = int(np.flatnonzero(final_cand.candidate_id.to_numpy(int) == cid)[0])
    rec = best.to_frame().T
    rec.to_csv(out / "recommended_peak_exp_floor.csv", index=False)
    curve = pd.DataFrame({
        "T_K": Td,
        "Kc_MPa_sqrtm": Kf[j],
        "target_Kc_MPa_sqrtm": _interp_target(targets, Td),
    })
    curve.to_csv(out / "recommended_peak_curve_dense.csv", index=False)

    # Rate sensitivity is diagnostic only; Kdot is not a fitted material parameter.
    rates = [float(x) for x in args.rates.split()]
    Trate = np.array(sorted(set([300, 500, 700, 800, 825, 850, 875, 890, 900, 905, 910, 920, 935, 950, 975, 1000, 1050, 1100, 1200])), float)
    rate_rows = []
    from arrhenius_fracture.config import ElasticProperties
    mat = ElasticProperties()
    one = final_cand.iloc[[j]].copy().reset_index(drop=True)
    for kd in rates:
        Kr = base.simulate_candidates(one, Trate, Kmax=args.Kmax, dK=args.final_dK,
                                      Kdot=kd, G_Pa=mat.G, nu=mat.nu, b_m=mat.b)[0]
        pf = peak_scores(Kr[None, :], Trate, targets, min_prominence=0.0).iloc[0]
        for T, kval in zip(Trate, Kr):
            rate_rows.append({"Kdot": kd, "T_K": T, "Kc_MPa_sqrtm": kval,
                              "peak_T_K": pf.peak_T_K, "peak_Kc": pf.peak_Kc,
                              "peak_prominence": pf.peak_prominence})
    rate = pd.DataFrame(rate_rows)
    rate.to_csv(out / "peak_rate_sensitivity.csv", index=False)
    save_plot(curve, rate, targets, out)

    cfg = vars(args).copy()
    cfg.update({"focused_temperatures_K": Tfocus.tolist(), "dense_temperatures_K": Td.tolist(),
                "target_features": target_features(targets),
                "calibration_note": "Kdot fixed during fit; rates varied only after fitting"})
    with (out / "run_config.json").open("w") as f:
        json.dump(cfg, f, indent=2)
    print("\nBest focused peak candidate:")
    cols = [c for c in ["candidate_id", "surface_id", "cleave_G00_eV", "cleave_gT_eV_per_K",
                        "cleave_sigc0_GPa", "cleave_sT_MPa_per_K", "cleave_exp_a", "cleave_exp_n",
                        "cleave_floor_frac", "cleave_S_hs_kB", "chi_shield", "N_sat", "peak_score",
                        "peak_T_K", "peak_Kc", "peak_prominence"] if c in rec.columns]
    print(rec[cols].to_string(index=False))
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
