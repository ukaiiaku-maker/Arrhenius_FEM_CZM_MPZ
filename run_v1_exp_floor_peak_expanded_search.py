#!/usr/bin/env python3
"""Expanded multi-basin V1 search for a strong narrow peak-like Kc(T) response.

The old prior 1-D Kc(T) curve is the observable target.  This search addresses
three restrictions in the first focused refinement:

1. seeds are retained from several peak-producing basins rather than only the
   lowest scalar-score basin;
2. both cleavage and emission EXP-floor parameters are varied continuously;
3. N_sat may switch between finite and effectively unsaturated states.

The calibration loading rate is fixed.  Rate sensitivity is evaluated only
post-fit.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

import run_v1_exp_floor_four_class_tuning as base
import run_v1_exp_floor_peak_focused_refinement as focused


# Deliberately wider than the original generic four-class screen.
EXPANDED_RANGES = {
    "cleave_G00_eV": (0.4, 4.5),
    "cleave_gT_eV_per_K": (-0.005, 0.012),
    "cleave_sigc0_GPa": (0.8, 10.0),
    "cleave_sT_MPa_per_K": (-8.0, 6.0),
    "cleave_exp_a": (0.08, 2.5),
    "cleave_exp_n": (0.35, 2.5),
    "cleave_floor_frac": (5e-4, 0.12),
    "cleave_S_hs_kB": (-70.0, 70.0),
    "chi_shield": (0.0, 1.0),
    "N_sat_finite": (300.0, 50000.0),
    "exp_G00_eV": (0.5, 3.5),
    "exp_sigc0_GPa": (0.35, 6.0),
    "exp_a": (0.03, 1.4),
    "exp_n": (0.35, 2.2),
    "exp_floor_frac": (0.002, 0.10),
    "exp_gT_eV_per_K": (-0.004, 0.028),
    "exp_sT_MPa_per_K": (-9.0, 6.0),
}

CLEAVE_KEYS = [
    "cleave_G00_eV", "cleave_gT_eV_per_K", "cleave_sigc0_GPa",
    "cleave_sT_MPa_per_K", "cleave_exp_a", "cleave_exp_n",
    "cleave_floor_frac", "cleave_S_hs_kB", "chi_shield", "N_sat",
]
EMIT_KEYS = [
    "exp_G00_eV", "exp_sigc0_GPa", "exp_a", "exp_n", "exp_floor_frac",
    "exp_Tref_K", "exp_gT_eV_per_K", "exp_sT_MPa_per_K",
]
META_KEYS = [
    "context_id", "surface_id", "surface_index", "thermal_shape_stratum",
    "eta_G_Tref_over_G00", "eta_sigc_Tref_over_sigc0", "implied_S0_kB",
]


def _clip(x: float, key: str) -> float:
    lo, hi = EXPANDED_RANGES[key]
    return float(np.clip(float(x), lo, hi))


def _logmut(x: float, sd: float, rng: np.random.Generator, key: str) -> float:
    return _clip(float(x) * math.exp(rng.normal(0.0, sd)), key)


def model_row(r: pd.Series | dict) -> dict:
    d = r.to_dict() if hasattr(r, "to_dict") else dict(r)
    keep = [k for k in CLEAVE_KEYS + EMIT_KEYS + META_KEYS if k in d]
    out = {k: d[k] for k in keep}
    return out


def target_feature_distance(scored: pd.DataFrame, feat: dict[str, float]) -> np.ndarray:
    # Dimensionless feature-space metric.  Scales intentionally reflect the
    # narrowness and amplitude of the target peak rather than generic RMS.
    d = (
        ((scored.peak_T_K.to_numpy(float) - feat["T_peak_target_K"]) / 20.0) ** 2
        + ((scored.peak_Kc.to_numpy(float) - feat["K_peak_target"]) / 5.0) ** 2
        + ((scored.K_pre.to_numpy(float) - feat["K_pre_target"]) / 3.0) ** 2
        + ((scored.K_post.to_numpy(float) - feat["K_post_target"]) / 2.0) ** 2
        + ((scored.peak_rise.to_numpy(float) - feat["rise_target"]) / 5.0) ** 2
        + ((scored.peak_fall.to_numpy(float) - feat["fall_target"]) / 5.0) ** 2
    )
    return d


def expanded_peak_scores(Kc: np.ndarray, T: np.ndarray, targets: pd.DataFrame,
                         min_prominence: float = 2.0) -> pd.DataFrame:
    q = focused.peak_scores(Kc, T, targets, min_prominence=min_prominence)
    feat = focused.target_features(targets)
    fd = target_feature_distance(q, feat)
    # Stronger emphasis on amplitude and shoulders than in the first focused
    # run; retain a small global-shape term to suppress pathological spikes.
    amp_pen = np.abs(q.peak_Kc - feat["K_peak_target"]) / max(feat["K_peak_target"], 1.0)
    pre_pen = np.abs(q.K_pre - feat["K_pre_target"]) / max(feat["K_pre_target"], 2.0)
    post_pen = np.abs(q.K_post - feat["K_post_target"]) / max(feat["K_post_target"], 2.0)
    rise_pen = np.abs(q.peak_rise - feat["rise_target"]) / max(feat["rise_target"], 1.0)
    fall_pen = np.abs(q.peak_fall - feat["fall_target"]) / max(feat["fall_target"], 1.0)
    T_pen = np.abs(q.peak_T_K - feat["T_peak_target_K"]) / 25.0
    topo = np.where(q.peak_valid.to_numpy(bool), 0.0, 3.0)
    score = (
        0.14 * q.weighted_rel_rmse.to_numpy(float)
        + 0.05 * q.weighted_log_rmse.to_numpy(float)
        + 0.20 * amp_pen.to_numpy(float)
        + 0.14 * pre_pen.to_numpy(float)
        + 0.14 * post_pen.to_numpy(float)
        + 0.10 * rise_pen.to_numpy(float)
        + 0.10 * fall_pen.to_numpy(float)
        + 0.08 * T_pen.to_numpy(float)
        + 0.05 * np.sqrt(fd)
        + topo
    )
    q = q.copy()
    q["expanded_peak_score"] = score
    q["feature_distance"] = fd
    return q


def add_scores(cand: pd.DataFrame, K: np.ndarray, T: np.ndarray,
               targets: pd.DataFrame, min_prominence: float) -> pd.DataFrame:
    s = expanded_peak_scores(K, T, targets, min_prominence=min_prominence)
    return pd.concat([cand.reset_index(drop=True), s.reset_index(drop=True)], axis=1)


def diverse_take(df: pd.DataFrame, n: int, key: str) -> pd.DataFrame:
    if len(df) == 0 or n <= 0:
        return df.iloc[:0].copy()
    work = df.sort_values(key).copy()
    rows = []
    surf_count: dict[str, int] = {}
    ctx_count: dict[str, int] = {}
    for _, r in work.iterrows():
        surf = str(r.get("surface_id", "")); ctx = str(r.get("context_id", ""))
        if surf_count.get(surf, 0) >= max(3, n // 4):
            continue
        if ctx_count.get(ctx, 0) >= 2:
            continue
        rows.append(r)
        surf_count[surf] = surf_count.get(surf, 0) + 1
        ctx_count[ctx] = ctx_count.get(ctx, 0) + 1
        if len(rows) >= n:
            break
    if len(rows) < n:
        used = {int(r.candidate_id) for r in rows if "candidate_id" in r}
        for _, r in work.iterrows():
            if int(r.candidate_id) in used:
                continue
            rows.append(r)
            if len(rows) >= n:
                break
    return pd.DataFrame(rows).reset_index(drop=True)



def _lin(u: np.ndarray, lo: float, hi: float) -> np.ndarray:
    return lo + np.asarray(u, float) * (hi - lo)


def _logmap(u: np.ndarray, lo: float, hi: float) -> np.ndarray:
    return np.exp(np.log(lo) + np.asarray(u, float) * (np.log(hi) - np.log(lo)))


def generate_global_sobol(n: int, seed: int) -> pd.DataFrame:
    """Continuous expanded-space Sobol design over both barrier channels.

    This is intentionally independent of the 96-surface discrete emission
    design used by the first broad screen.  It is the safeguard against a
    narrow peak basin falling between those discrete emission surfaces.
    """
    if n <= 0:
        return pd.DataFrame()
    from scipy.stats import qmc
    m = int(math.ceil(math.log2(max(n, 2))))
    U = qmc.Sobol(d=17, scramble=True, seed=seed).random_base2(m)[:n]
    R = EXPANDED_RANGES
    d = {
        "cleave_G00_eV": _lin(U[:,0], *R["cleave_G00_eV"]),
        "cleave_gT_eV_per_K": _lin(U[:,1], *R["cleave_gT_eV_per_K"]),
        "cleave_sigc0_GPa": _lin(U[:,2], *R["cleave_sigc0_GPa"]),
        "cleave_sT_MPa_per_K": _lin(U[:,3], *R["cleave_sT_MPa_per_K"]),
        "cleave_exp_a": _logmap(U[:,4], *R["cleave_exp_a"]),
        "cleave_exp_n": _logmap(U[:,5], *R["cleave_exp_n"]),
        "cleave_floor_frac": _logmap(U[:,6], *R["cleave_floor_frac"]),
        "cleave_S_hs_kB": _lin(U[:,7], *R["cleave_S_hs_kB"]),
        "chi_shield": _lin(U[:,8], *R["chi_shield"]),
        "exp_G00_eV": _lin(U[:,10], *R["exp_G00_eV"]),
        "exp_sigc0_GPa": _lin(U[:,11], *R["exp_sigc0_GPa"]),
        "exp_a": _logmap(U[:,12], *R["exp_a"]),
        "exp_n": _logmap(U[:,13], *R["exp_n"]),
        "exp_floor_frac": _logmap(U[:,14], *R["exp_floor_frac"]),
        "exp_gT_eV_per_K": _lin(U[:,15], *R["exp_gT_eV_per_K"]),
        "exp_sT_MPa_per_K": _lin(U[:,16], *R["exp_sT_MPa_per_K"]),
    }
    ns = np.full(n, np.inf)
    finite = U[:,9] >= 0.25
    uf = np.clip((U[finite,9] - 0.25) / 0.75, 0.0, 1.0)
    ns[finite] = _logmap(uf, *R["N_sat_finite"])
    d["N_sat"] = ns
    out = pd.DataFrame(d)
    out["exp_Tref_K"] = 300.0
    out["context_id"] = [f"GLOBAL_CTX_{i:06d}" for i in range(n)]
    out["surface_id"] = [f"GLOBAL_EXP_{i:06d}" for i in range(n)]
    out["surface_index"] = -1
    out["thermal_shape_stratum"] = "continuous_global"
    out["eta_G_Tref_over_G00"] = np.nan
    out["eta_sigc_Tref_over_sigc0"] = np.nan
    out["implied_S0_kB"] = np.nan
    out.insert(0, "candidate_id", np.arange(n, dtype=int))
    return out


def global_grid() -> np.ndarray:
    return np.array([300, 500, 700, 800, 825, 850, 875, 890, 900, 905,
                     910, 920, 935, 950, 975, 1000, 1100, 1200], float)

def build_multibasin_seed_pool(broad_scored: pd.DataFrame, focused_scored: pd.DataFrame | None,
                               targets: pd.DataFrame, n_seeds: int) -> pd.DataFrame:
    feat = focused.target_features(targets)
    b = broad_scored.copy()
    b["feature_distance"] = target_feature_distance(b, feat)
    valid = b[b.peak_valid].copy()
    pools = []
    q = max(4, n_seeds // 4)
    pools.append(diverse_take(valid, q, "expanded_peak_score"))
    pools.append(diverse_take(valid, q, "feature_distance"))
    # Strong-peak pool near the desired transition temperature, ranked by
    # amplitude mismatch plus shoulder mismatch rather than scalar score.
    nearT = valid[(valid.peak_T_K >= feat["T_peak_target_K"] - 40)
                  & (valid.peak_T_K <= feat["T_peak_target_K"] + 40)].copy()
    nearT["amp_shoulder_score"] = (
        np.abs(nearT.peak_Kc - feat["K_peak_target"]) / 5.0
        + np.abs(nearT.K_pre - feat["K_pre_target"]) / 3.0
        + np.abs(nearT.K_post - feat["K_post_target"]) / 2.0
    )
    pools.append(diverse_take(nearT, q, "amp_shoulder_score"))
    prom = valid.sort_values("peak_prominence", ascending=False).head(max(4 * q, q)).copy()
    prom["neg_prom"] = -prom.peak_prominence
    pools.append(diverse_take(prom, q, "neg_prom"))

    if focused_scored is not None and len(focused_scored):
        f = focused_scored.copy()
        if "feature_distance" not in f:
            f["feature_distance"] = target_feature_distance(f, feat)
        if "expanded_peak_score" not in f:
            # Reconstruct a lightweight ranking from feature distance and prior score.
            f["expanded_peak_score"] = 0.2 * np.sqrt(f.feature_distance) + 0.8 * f.peak_score
        fv = f[f.peak_valid].copy()
        pools.append(diverse_take(fv, q, "expanded_peak_score"))
        pools.append(diverse_take(fv, q, "feature_distance"))

    seeds = pd.concat(pools, ignore_index=True, sort=False)
    # Deduplicate on parameter vector rather than stage-local ID.
    subset = [c for c in CLEAVE_KEYS + EMIT_KEYS if c in seeds.columns]
    seeds = seeds.drop_duplicates(subset=subset).reset_index(drop=True)
    # Final balanced ranking.
    seeds["seed_rank"] = np.minimum(
        seeds.get("expanded_peak_score", pd.Series(np.inf, index=seeds.index)).to_numpy(float),
        0.2 * np.sqrt(seeds.get("feature_distance", pd.Series(np.inf, index=seeds.index)).to_numpy(float))
        + 0.8 * seeds.get("peak_score", pd.Series(np.inf, index=seeds.index)).to_numpy(float),
    )
    return seeds.sort_values("seed_rank").head(n_seeds).reset_index(drop=True)


def mutate_candidate(b: dict, rng: np.random.Generator, shrink: float,
                     allow_surface_mutation: bool = True) -> dict:
    x = dict(b)
    # Cleavage EXP-floor + thermal/state coupling.
    x["cleave_G00_eV"] = _logmut(b["cleave_G00_eV"], 0.16 * shrink, rng, "cleave_G00_eV")
    x["cleave_gT_eV_per_K"] = _clip(b["cleave_gT_eV_per_K"] + rng.normal(0, 0.0012 * shrink), "cleave_gT_eV_per_K")
    x["cleave_sigc0_GPa"] = _logmut(b["cleave_sigc0_GPa"], 0.18 * shrink, rng, "cleave_sigc0_GPa")
    x["cleave_sT_MPa_per_K"] = _clip(b["cleave_sT_MPa_per_K"] + rng.normal(0, 0.9 * shrink), "cleave_sT_MPa_per_K")
    x["cleave_exp_a"] = _logmut(b["cleave_exp_a"], 0.20 * shrink, rng, "cleave_exp_a")
    x["cleave_exp_n"] = _logmut(b["cleave_exp_n"], 0.16 * shrink, rng, "cleave_exp_n")
    x["cleave_floor_frac"] = _logmut(b["cleave_floor_frac"], 0.30 * shrink, rng, "cleave_floor_frac")
    x["cleave_S_hs_kB"] = _clip(b["cleave_S_hs_kB"] + rng.normal(0, 7.0 * shrink), "cleave_S_hs_kB")
    x["chi_shield"] = _clip(b["chi_shield"] + rng.normal(0, 0.07 * shrink), "chi_shield")

    # Permit transitions between finite saturation and effectively unsaturated.
    ns = float(b["N_sat"])
    r = rng.random()
    if math.isfinite(ns):
        if r < 0.10 * shrink:
            x["N_sat"] = np.inf
        else:
            x["N_sat"] = float(np.clip(ns * math.exp(rng.normal(0, 0.35 * shrink)), *EXPANDED_RANGES["N_sat_finite"]))
    else:
        if r < 0.25 * shrink:
            lo, hi = EXPANDED_RANGES["N_sat_finite"]
            x["N_sat"] = float(math.exp(math.log(lo) + rng.random() * (math.log(hi) - math.log(lo))))
        else:
            x["N_sat"] = np.inf

    if allow_surface_mutation:
        x["exp_G00_eV"] = _logmut(b["exp_G00_eV"], 0.14 * shrink, rng, "exp_G00_eV")
        x["exp_sigc0_GPa"] = _logmut(b["exp_sigc0_GPa"], 0.18 * shrink, rng, "exp_sigc0_GPa")
        x["exp_a"] = _logmut(b["exp_a"], 0.20 * shrink, rng, "exp_a")
        x["exp_n"] = _logmut(b["exp_n"], 0.16 * shrink, rng, "exp_n")
        x["exp_floor_frac"] = _logmut(b["exp_floor_frac"], 0.25 * shrink, rng, "exp_floor_frac")
        x["exp_gT_eV_per_K"] = _clip(b["exp_gT_eV_per_K"] + rng.normal(0, 0.0018 * shrink), "exp_gT_eV_per_K")
        x["exp_sT_MPa_per_K"] = _clip(b["exp_sT_MPa_per_K"] + rng.normal(0, 0.9 * shrink), "exp_sT_MPa_per_K")
        x["surface_id"] = f"CONT_{b.get('surface_id', 'seed')}"
        x["surface_index"] = -1
    return x


def generate_population(seeds: pd.DataFrame, n_perturb: int, seed: int,
                        shrink: float, allow_surface_mutation: bool = True) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows: list[dict] = []
    for _, r in seeds.iterrows():
        b = model_row(r)
        b["parent_source_id"] = int(r.get("candidate_id", -1))
        rows.append(dict(b))
        for _ in range(n_perturb):
            rows.append(mutate_candidate(b, rng, shrink, allow_surface_mutation=allow_surface_mutation))
    out = pd.DataFrame(rows)
    out.insert(0, "candidate_id", np.arange(len(out), dtype=int))
    return out


def run_stage(cand: pd.DataFrame, T: np.ndarray, dK: float, Kdot: float, Kmax: float,
              targets: pd.DataFrame, out_npz: Path, out_csv: Path,
              min_prominence: float, resume: bool) -> tuple[np.ndarray, pd.DataFrame]:
    if resume and out_npz.exists() and out_csv.exists():
        z = np.load(out_npz)
        return z["Kc"], pd.read_csv(out_csv)
    from arrhenius_fracture.config import ElasticProperties
    mat = ElasticProperties()
    K = base.simulate_candidates(cand, T, Kmax=Kmax, dK=dK, Kdot=Kdot,
                                 G_Pa=mat.G, nu=mat.nu, b_m=mat.b)
    s = add_scores(cand, K, T, targets, min_prominence)
    np.savez_compressed(out_npz, Kc=K, T=T)
    s.to_csv(out_csv, index=False)
    return K, s


def select_generation_seeds(scored: pd.DataFrame, n: int, targets: pd.DataFrame) -> pd.DataFrame:
    feat = focused.target_features(targets)
    s = scored[scored.peak_valid].copy()
    if s.empty:
        s = scored.copy()
    s["feature_distance"] = target_feature_distance(s, feat)
    parts = [
        s.nsmallest(max(1, n // 2), "expanded_peak_score"),
        s.nsmallest(max(1, n // 3), "feature_distance"),
    ]
    # Preserve strong peak basin even if shoulders remain imperfect.
    near = s[(s.peak_T_K >= feat["T_peak_target_K"] - 35)
             & (s.peak_T_K <= feat["T_peak_target_K"] + 35)]
    if len(near):
        near = near.assign(amp_error=np.abs(near.peak_Kc - feat["K_peak_target"]))
        parts.append(near.nsmallest(max(1, n // 3), "amp_error"))
    out = pd.concat(parts, ignore_index=True).drop_duplicates(subset=["candidate_id"])
    return out.nsmallest(n, "expanded_peak_score").reset_index(drop=True) if len(out) > n else out.reset_index(drop=True)


def make_plot(curve: pd.DataFrame, targets: pd.DataFrame, out: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    tg = base.target_curve(targets, "peak").sort_values("T_K")
    fig, ax = plt.subplots(figsize=(7.4, 4.9))
    ax.plot(tg.T_K, tg.target_Kc_MPa_sqrtm, label="prior 1-D target")
    ax.plot(curve.T_K, curve.Kc_MPa_sqrtm, label="expanded EXP-floor fit")
    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel(r"$K_c$ (MPa$\sqrt{m}$)")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / "peak_expanded_fit.png", dpi=220)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source-run", default="runs/v1_exp_floor_four_class_tuning")
    ap.add_argument("--focused-run", default="runs/v1_exp_floor_peak_focused_refinement")
    ap.add_argument("--targets", default="exp_floor_four_class_target_curves.csv")
    ap.add_argument("--out", default="runs/v1_exp_floor_peak_expanded_search")
    ap.add_argument("--Kdot", type=float, default=0.005)
    ap.add_argument("--Kmax", type=float, default=100.0)
    ap.add_argument("--global-n", type=int, default=16384)
    ap.add_argument("--global-dK", type=float, default=0.5)
    ap.add_argument("--n-seeds", type=int, default=64)
    ap.add_argument("--gen1-perturb", type=int, default=72)
    ap.add_argument("--gen2-seeds", type=int, default=40)
    ap.add_argument("--gen2-perturb", type=int, default=80)
    ap.add_argument("--gen3-seeds", type=int, default=24)
    ap.add_argument("--gen3-perturb", type=int, default=64)
    ap.add_argument("--finalists", type=int, default=12)
    ap.add_argument("--gen1-dK", type=float, default=0.075)
    ap.add_argument("--gen2-dK", type=float, default=0.04)
    ap.add_argument("--gen3-dK", type=float, default=0.025)
    ap.add_argument("--final-dK", type=float, default=0.02)
    ap.add_argument("--min-prominence", type=float, default=2.0)
    ap.add_argument("--seed", type=int, default=20260712)
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()

    if args.smoke:
        args.global_n = min(args.global_n, 64)
        args.global_dK = max(args.global_dK, 1.0)
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
    frun = Path(args.focused_run)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    targets = base.load_targets(args.targets)

    broad = pd.read_csv(src / "broad_candidates.csv")
    bz = np.load(src / "broad_Kc.npz")
    bs = add_scores(broad, bz["Kc"], bz["T"], targets, args.min_prominence)
    bs.to_csv(out / "source_broad_expanded_scores.csv", index=False)

    fs = None
    if (frun / "generation2_scores.csv").exists():
        fs = pd.read_csv(frun / "generation2_scores.csv")
        # Re-score the focused generation using the expanded objective if K is available.
        if (frun / "generation2_Kc.npz").exists() and (frun / "generation2_candidates.csv").exists():
            fz = np.load(frun / "generation2_Kc.npz")
            fc = pd.read_csv(frun / "generation2_candidates.csv")
            fs = add_scores(fc, fz["Kc"], fz["T"], targets, args.min_prominence)
            fs.to_csv(out / "source_focused_expanded_scores.csv", index=False)

    # Independent continuous global screen in the expanded parameter box.
    gglobal = generate_global_sobol(args.global_n, args.seed - 1)
    Tg = global_grid()
    print(f"global continuous screen: {len(gglobal)} candidates; dK={args.global_dK}", flush=True)
    Kg, sg = run_stage(gglobal, Tg, args.global_dK, args.Kdot, args.Kmax, targets,
                       out / "global_Kc.npz", out / "global_scores.csv",
                       args.min_prominence, args.resume)
    gglobal.to_csv(out / "global_candidates.csv", index=False)

    # Half the seed budget comes from existing known peak basins and half from
    # the new continuous global screen.  This balances exploitation/exploration.
    n_old = max(1, args.n_seeds // 2)
    old_seeds = build_multibasin_seed_pool(bs, fs, targets, n_old)
    global_seeds = select_generation_seeds(sg, args.n_seeds - n_old, targets)
    seeds = pd.concat([old_seeds, global_seeds], ignore_index=True, sort=False)
    subset = [c for c in CLEAVE_KEYS + EMIT_KEYS if c in seeds.columns]
    seeds = seeds.drop_duplicates(subset=subset).head(args.n_seeds).reset_index(drop=True)
    seeds.to_csv(out / "multibasin_seeds.csv", index=False)
    print(f"selected {len(seeds)} combined old-basin + global seeds", flush=True)

    Tfocus = focused.focused_grid()

    g1 = generate_population(seeds, args.gen1_perturb, args.seed, shrink=1.0, allow_surface_mutation=True)
    g1.to_csv(out / "generation1_candidates.csv", index=False)
    print(f"generation 1: {len(g1)} candidates; dK={args.gen1_dK}", flush=True)
    K1, s1 = run_stage(g1, Tfocus, args.gen1_dK, args.Kdot, args.Kmax, targets,
                       out / "generation1_Kc.npz", out / "generation1_scores.csv",
                       args.min_prominence, args.resume)

    seed2 = select_generation_seeds(s1, args.gen2_seeds, targets)
    seed2.to_csv(out / "generation2_seeds.csv", index=False)
    g2 = generate_population(seed2, args.gen2_perturb, args.seed + 1, shrink=0.55, allow_surface_mutation=True)
    g2.to_csv(out / "generation2_candidates.csv", index=False)
    print(f"generation 2: {len(g2)} candidates; dK={args.gen2_dK}", flush=True)
    K2, s2 = run_stage(g2, Tfocus, args.gen2_dK, args.Kdot, args.Kmax, targets,
                       out / "generation2_Kc.npz", out / "generation2_scores.csv",
                       args.min_prominence, args.resume)

    seed3 = select_generation_seeds(s2, args.gen3_seeds, targets)
    seed3.to_csv(out / "generation3_seeds.csv", index=False)
    g3 = generate_population(seed3, args.gen3_perturb, args.seed + 2, shrink=0.28, allow_surface_mutation=True)
    g3.to_csv(out / "generation3_candidates.csv", index=False)
    print(f"generation 3: {len(g3)} candidates; dK={args.gen3_dK}", flush=True)
    K3, s3 = run_stage(g3, Tfocus, args.gen3_dK, args.Kdot, args.Kmax, targets,
                       out / "generation3_Kc.npz", out / "generation3_scores.csv",
                       args.min_prominence, args.resume)

    finalists = select_generation_seeds(s3, args.finalists, targets)
    ids = set(finalists.candidate_id.astype(int))
    fcand = g3[g3.candidate_id.astype(int).isin(ids)].copy().reset_index(drop=True)
    Td = base.make_temperature_grid(targets, "dense")
    print(f"final: {len(fcand)} candidates; dense {len(Td)} temperatures; dK={args.final_dK}", flush=True)
    Kf, sf = run_stage(fcand, Td, args.final_dK, args.Kdot, args.Kmax, targets,
                       out / "final_Kc.npz", out / "final_scores.csv",
                       args.min_prominence, args.resume)
    best = sf.sort_values("expanded_peak_score").iloc[0]
    bid = int(best.candidate_id)
    jj = int(np.flatnonzero(fcand.candidate_id.to_numpy(int) == bid)[0])
    best.to_frame().T.to_csv(out / "recommended_peak_expanded.csv", index=False)
    curve = pd.DataFrame({
        "T_K": Td,
        "Kc_MPa_sqrtm": Kf[jj],
        "target_Kc_MPa_sqrtm": base.interp_target(targets, "peak", Td),
    })
    curve.to_csv(out / "recommended_peak_curve_dense.csv", index=False)
    make_plot(curve, targets, out)

    cfg = vars(args).copy()
    cfg["target_features"] = focused.target_features(targets)
    cfg["expanded_ranges"] = {k: list(v) for k, v in EXPANDED_RANGES.items()}
    cfg["notes"] = [
        "Kdot fixed during calibration",
        "continuous mutation of both cleavage and emission EXP-floor parameters",
        "finite/infinite N_sat switching enabled",
        "multi-basin seed preservation used to avoid collapse onto weak-peak basin",
    ]
    with (out / "run_config.json").open("w") as f:
        json.dump(cfg, f, indent=2)

    cols = [c for c in [
        "candidate_id", "surface_id", "cleave_G00_eV", "cleave_gT_eV_per_K",
        "cleave_sigc0_GPa", "cleave_sT_MPa_per_K", "cleave_exp_a", "cleave_exp_n",
        "cleave_floor_frac", "cleave_S_hs_kB", "chi_shield", "N_sat",
        "exp_G00_eV", "exp_sigc0_GPa", "exp_a", "exp_n", "exp_floor_frac",
        "exp_gT_eV_per_K", "exp_sT_MPa_per_K", "expanded_peak_score",
        "peak_T_K", "peak_Kc", "K_pre", "K_post", "peak_prominence",
    ] if c in best.index]
    print("\nBest expanded peak candidate:")
    print(best[cols].to_frame().T.to_string(index=False))
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
