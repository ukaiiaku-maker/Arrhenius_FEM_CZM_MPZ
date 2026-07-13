#!/usr/bin/env python3
"""Search fully EXP-floor V1 barrier parameterizations for four Kc(T) targets.

The calibration target is the *observable response* from the prior reduced
model (the curves in forward_input_parameter_KT_overlay), not its old
linearized cleavage-barrier parameters.

Physics retained from the later V5.7 monotonic V1 architecture:
  * EXP-floor crack-tip emission surface,
  * EXP-floor crack-opening surface,
  * temperature-dependent G0(T) and sigmac(T) for both channels,
  * optional high-stress entropy crossover on cleavage,
  * emission accumulation, saturation, back-stress shielding and blunting,
  * stored-energy embrittlement,
  * cooperative m-hit cleavage renewal,
  * monotonic K ramp and first-passage failure.

The search is staged so it is practical on a laptop:
  1. Sobol broad cleavage/state contexts x the supplied 96 emission surfaces.
  2. Local perturbation around the best candidates for each target class.
  3. Fine re-evaluation of the finalists and dense target-grid curves.

Primary outputs:
  broad_candidate_scores.csv
  local_candidate_scores.csv
  finalist_scores.csv
  recommended_exp_floor_four_class.csv
  recommended_curves_dense.csv
  exp_floor_four_class_tuning.png
  run_config.json
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
from scipy.special import gammainc
from scipy.stats import qmc

KBEV = 8.617333262145e-5
TARGET_CLASSES = ("ceramic", "peak", "weakT", "DBTT")
CLASS_ALIASES = {"dbtt": "DBTT", "DBTT": "DBTT", "weakT": "weakT", "peak": "peak", "ceramic": "ceramic"}

# Broad ranges are intentionally wider than the five later V8 cleavage anchors,
# while remaining in the same EXP-floor family used by the production code.
RANGES = {
    "cleave_G00_eV": (0.6, 3.2),
    "cleave_gT_eV_per_K": (-0.0030, 0.0100),
    "cleave_sigc0_GPa": (1.2, 7.0),
    "cleave_sT_MPa_per_K": (-5.0, 4.0),
    "cleave_exp_a": (0.12, 1.8),
    "cleave_exp_n": (0.45, 2.0),
    "cleave_floor_frac": (0.0015, 0.08),
    "cleave_S_hs_kB": (-45.0, 45.0),
    "chi_shield": (0.0, 0.8),
    "N_sat_finite": (500.0, 12000.0),
}


def _loguniform(u: np.ndarray, lo: float, hi: float) -> np.ndarray:
    return np.exp(np.log(lo) + u * (np.log(hi) - np.log(lo)))


def _clip(v: np.ndarray, key: str) -> np.ndarray:
    lo, hi = RANGES[key]
    return np.clip(v, lo, hi)


def load_targets(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    # Accept either the packaged target schema or the original overlay prediction.
    if {"regime_key", "T_K", "Kc_pred_MPa_sqrtm"}.issubset(df.columns):
        out = df[["regime_key", "T_K", "Kc_pred_MPa_sqrtm"]].copy()
        out["target_class"] = out["regime_key"].map(lambda x: CLASS_ALIASES.get(str(x), str(x)))
        out = out.rename(columns={"Kc_pred_MPa_sqrtm": "target_Kc_MPa_sqrtm"})
        out["source"] = "prior_1D_forward_model"
        return out[["target_class", "T_K", "target_Kc_MPa_sqrtm", "source"]].sort_values(["target_class", "T_K"])
    req = {"target_class", "T_K", "target_Kc_MPa_sqrtm"}
    missing = req - set(df.columns)
    if missing:
        raise SystemExit(f"target CSV missing columns: {sorted(missing)}")
    df = df.copy()
    df["target_class"] = df["target_class"].map(lambda x: CLASS_ALIASES.get(str(x), str(x)))
    if "source" not in df.columns:
        df["source"] = "user_target"
    return df.sort_values(["target_class", "T_K"])


def target_curve(targets: pd.DataFrame, klass: str) -> pd.DataFrame:
    g = targets[targets.target_class.astype(str) == klass].sort_values("T_K")
    if g.empty:
        raise ValueError(f"no target curve for class {klass}")
    return g


def make_temperature_grid(targets: pd.DataFrame, stage: str) -> np.ndarray:
    allT = np.array(sorted(targets.T_K.unique()), float)
    if stage == "broad":
        wanted = set(np.arange(300.0, 1200.1, 100.0))
        wanted.update(np.arange(850.0, 950.1, 25.0))
    elif stage == "local":
        wanted = set(np.arange(300.0, 1200.1, 50.0))
        wanted.update(np.arange(850.0, 950.1, 10.0))
    elif stage == "dense":
        return allT
    else:
        raise ValueError(stage)
    return np.array([T for T in allT if any(abs(T - w) < 1e-9 for w in wanted)], float)


def make_broad_contexts(n_contexts: int, seed: int) -> pd.DataFrame:
    d = 10
    m = int(math.ceil(math.log2(max(2, n_contexts))))
    U = qmc.Sobol(d=d, scramble=True, seed=seed).random_base2(m=m)[:n_contexts]
    out = pd.DataFrame({
        "cleave_G00_eV": RANGES["cleave_G00_eV"][0] + U[:, 0] * np.diff(RANGES["cleave_G00_eV"])[0],
        "cleave_gT_eV_per_K": RANGES["cleave_gT_eV_per_K"][0] + U[:, 1] * np.diff(RANGES["cleave_gT_eV_per_K"])[0],
        "cleave_sigc0_GPa": RANGES["cleave_sigc0_GPa"][0] + U[:, 2] * np.diff(RANGES["cleave_sigc0_GPa"])[0],
        "cleave_sT_MPa_per_K": RANGES["cleave_sT_MPa_per_K"][0] + U[:, 3] * np.diff(RANGES["cleave_sT_MPa_per_K"])[0],
        "cleave_exp_a": RANGES["cleave_exp_a"][0] + U[:, 4] * np.diff(RANGES["cleave_exp_a"])[0],
        "cleave_exp_n": RANGES["cleave_exp_n"][0] + U[:, 5] * np.diff(RANGES["cleave_exp_n"])[0],
        "cleave_floor_frac": _loguniform(U[:, 6], *RANGES["cleave_floor_frac"]),
        "cleave_S_hs_kB": RANGES["cleave_S_hs_kB"][0] + U[:, 7] * np.diff(RANGES["cleave_S_hs_kB"])[0],
        "chi_shield": RANGES["chi_shield"][0] + U[:, 8] * np.diff(RANGES["chi_shield"])[0],
    })
    finite = U[:, 9] >= 0.20
    nsat = np.full(n_contexts, np.inf)
    uf = np.clip((U[finite, 9] - 0.20) / 0.80, 0.0, 1.0)
    nsat[finite] = _loguniform(uf, *RANGES["N_sat_finite"])
    out["N_sat"] = nsat
    out.insert(0, "context_id", [f"CTX_{i:05d}" for i in range(n_contexts)])
    return out


def cross_emission_contexts(design: pd.DataFrame, contexts: pd.DataFrame) -> pd.DataFrame:
    nctx, nem = len(contexts), len(design)
    c = contexts.loc[contexts.index.repeat(nem)].reset_index(drop=True)
    e = pd.concat([design.reset_index(drop=True)] * nctx, ignore_index=True)
    out = pd.concat([c, e], axis=1)
    out.insert(0, "candidate_id", np.arange(len(out), dtype=int))
    return out


def _emit_terms(cand: pd.DataFrame, T: float):
    G0 = np.maximum(
        0.75 * cand["exp_G00_eV"].to_numpy(float)
        + 0.75 * cand["exp_gT_eV_per_K"].to_numpy(float) * (T - cand["exp_Tref_K"].to_numpy(float)),
        1e-10,
    )
    sigc = np.maximum(
        cand["exp_sigc0_GPa"].to_numpy(float) * 1e9
        + cand["exp_sT_MPa_per_K"].to_numpy(float) * 1e6 * (T - cand["exp_Tref_K"].to_numpy(float)),
        1e6,
    )
    ff = cand["exp_floor_frac"].to_numpy(float)
    floor = np.minimum(0.95 * G0, np.maximum(1e-4 * 0.75, ff * G0))
    return G0, sigc, floor, cand["exp_a"].to_numpy(float), cand["exp_n"].to_numpy(float)


def _exp_floor(G0, floor, sigc, a, n, sigma):
    x = np.maximum(np.asarray(sigma, float), 0.0) / np.maximum(sigc, 1e6)
    return np.maximum(floor + (G0 - floor) * np.exp(-a * np.power(x, n)), 0.0)


def _cleave_G_eV(cand: pd.DataFrame, idx: np.ndarray, sigma: np.ndarray, T: float) -> np.ndarray:
    G00 = cand["cleave_G00_eV"].to_numpy(float)[idx]
    gT = cand["cleave_gT_eV_per_K"].to_numpy(float)[idx]
    sig0 = cand["cleave_sigc0_GPa"].to_numpy(float)[idx] * 1e9
    sT = cand["cleave_sT_MPa_per_K"].to_numpy(float)[idx] * 1e6
    a = cand["cleave_exp_a"].to_numpy(float)[idx]
    n = cand["cleave_exp_n"].to_numpy(float)[idx]
    ff = cand["cleave_floor_frac"].to_numpy(float)[idx]
    Shs = cand["cleave_S_hs_kB"].to_numpy(float)[idx]

    dT = T - 300.0
    G0 = np.maximum(G00 + gT * dT, 1e-9)
    sigc = np.maximum(sig0 + sT * dT, 1e6)
    floor = np.minimum(0.95 * G0, np.maximum(1e-4, ff * G0))
    x = np.maximum(sigma, 0.0) / sigc
    G = floor + (G0 - floor) * np.exp(-a * np.power(x, n))
    # Same optional high-stress entropy crossover used by FractureBarrier exp_floor.
    xs = np.power(np.maximum(sigma, 0.0) / 6.0e9, 2.0)
    gate = xs / (1.0 + xs)
    G = G - dT * Shs * gate * KBEV
    return np.maximum(G, 0.0)


def _logmean(lo: np.ndarray, hi: np.ndarray) -> np.ndarray:
    lo = np.maximum(np.asarray(lo, float), 0.0)
    hi = np.maximum(np.asarray(hi, float), 0.0)
    a = np.minimum(lo, hi); b = np.maximum(lo, hi)
    out = np.empty_like(a)
    zero = a <= 0.0
    close = (~zero) & (np.abs(b-a) <= 1e-12*np.maximum(b,1.0))
    reg = (~zero) & (~close)
    out[zero] = 0.5*b[zero]
    out[close] = b[close]
    out[reg] = (b[reg]-a[reg]) / np.log(b[reg]/a[reg])
    return out


def simulate_candidates(cand: pd.DataFrame, temperatures: Sequence[float], *, Kmax: float, dK: float,
                        Kdot: float, G_Pa: float, nu: float, b_m: float) -> np.ndarray:
    """Vectorized monotonic V1 first-passage calculation."""
    nC = len(cand); temps = np.asarray(temperatures, float)
    Kc = np.full((nC, len(temps)), np.nan)

    r0 = 1e-6; sigma_cap = 30e9; nu0_c = 1e12; nu0_e = 1e11
    m_hits = 3.0; tau_c = 1e-6; dN_cap = 50.0
    beta_back = 1.0; L_pz = 1e-6; rho0 = 5e12
    c_blunt = 1.0; v_emb_b3 = 500.0; emb_sat_frac = 1.0
    back_coeff = beta_back * G_Pa * b_m / (2*np.pi*(1-nu)*L_pz)
    stored_coeff_eV = (0.5*G_Pa*b_m**2) * (v_emb_b3*b_m**3) / 1.602176634e-19
    chi_all = cand["chi_shield"].to_numpy(float)
    Nsat_all = cand["N_sat"].to_numpy(float)

    dK_SI = dK*1e6; Kmax_SI=Kmax*1e6; Kdot_SI=Kdot*1e6
    dt = dK_SI/max(Kdot_SI,1e-300); nsteps=int(math.ceil(Kmax_SI/dK_SI))

    for it,T in enumerate(temps):
        G0e,sigce,Gfe,ae,ne = _emit_terms(cand,T)
        Nem=np.zeros(nC); B=np.zeros(nC); Kprev=np.zeros(nC); active=np.ones(nC,bool)
        for istep in range(1,nsteps+1):
            if not active.any(): break
            idx=np.flatnonzero(active); K=min(istep*dK_SI,Kmax_SI); Ni=Nem[idx]
            r_eff=r0+c_blunt*b_m*Ni
            sig_tip=np.minimum(K/np.sqrt(2*np.pi*r_eff),sigma_cap)
            sig_back=back_coeff*Ni
            sig_emit=np.maximum(sig_tip-sig_back,0.0)
            Ge=_exp_floor(G0e[idx],Gfe[idx],sigce[idx],ae[idx],ne[idx],sig_emit)
            lam_e=nu0_e*np.exp(np.clip(-Ge/max(KBEV*T,1e-30),-700,0))
            prod=np.minimum(lam_e*dt,dN_cap)
            ns=Nsat_all[idx]; finite=np.isfinite(ns)&(ns>0)
            if finite.any(): prod[finite]*=np.maximum(1.0-Ni[finite]/ns[finite],0.0)
            Ni2=np.maximum(Ni+prod,0.0); Nem[idx]=Ni2

            r_eff2=r0+c_blunt*b_m*Ni2
            sig_tip2=np.minimum(K/np.sqrt(2*np.pi*r_eff2),sigma_cap)
            sig_back2=back_coeff*Ni2
            sig_c=np.maximum(sig_tip2-chi_all[idx]*sig_back2,0.0)
            Gc=_cleave_G_eV(cand,idx,sig_c,T)
            rho=rho0+Ni2/(L_pz**2)
            dGemb=stored_coeff_eV*rho
            Geff=np.maximum(Gc-np.minimum(dGemb,emb_sat_frac*Gc),0.0)
            lam_raw=nu0_c*np.exp(np.clip(-Geff/max(KBEV*T,1e-30),-700,0))
            lam_c=gammainc(m_hits,np.minimum(lam_raw*tau_c,1e12))/tau_c

            kp=Kprev[idx]
            sig_prev=np.zeros_like(sig_tip2); pos=kp>0
            sig_prev[pos]=np.minimum(kp[pos]/np.sqrt(2*np.pi*r_eff2[pos]),sigma_cap)
            sig_cp=np.maximum(sig_prev-chi_all[idx]*sig_back2,0.0)
            Gcp=_cleave_G_eV(cand,idx,sig_cp,T)
            Geffp=np.maximum(Gcp-np.minimum(dGemb,emb_sat_frac*Gcp),0.0)
            lam_raw_p=nu0_c*np.exp(np.clip(-Geffp/max(KBEV*T,1e-30),-700,0))
            lam_p=gammainc(m_hits,np.minimum(lam_raw_p*tau_c,1e12))/tau_c
            lam_eff=lam_c.copy(); lam_eff[pos]=_logmean(lam_p[pos],lam_c[pos])
            Bnew=B[idx]+lam_eff*dt; B[idx]=Bnew; Kprev[idx]=K
            fire=Bnew>=1.0
            if fire.any():
                fi=idx[fire]; Kc[fi,it]=K/1e6; active[fi]=False
    return Kc


def interp_target(targets: pd.DataFrame, klass: str, Tgrid: np.ndarray) -> np.ndarray:
    g=target_curve(targets,klass)
    return np.interp(Tgrid,g.T_K.to_numpy(float),g.target_Kc_MPa_sqrtm.to_numpy(float))


def score_curves(Kc: np.ndarray, Tgrid: np.ndarray, targets: pd.DataFrame, klass: str, Kmax: float):
    target=interp_target(targets,klass,Tgrid)
    cens=~np.isfinite(Kc); p=np.where(cens,1.15*Kmax,Kc)
    scale=np.maximum(target,2.0)
    rel=np.sqrt(np.mean(((p-target[None,:])/scale[None,:])**2,axis=1))
    log=np.sqrt(np.mean((np.log10(p+0.25)-np.log10(target[None,:]+0.25))**2,axis=1))
    pnorm=p/np.maximum(p[:,[0]],1e-6); tnorm=target/max(target[0],1e-6)
    shape=np.sqrt(np.mean((pnorm-tnorm[None,:])**2,axis=1))
    dp=np.diff(p,axis=1)/np.maximum(p[:,[0]],1e-6)
    dt=np.diff(target)/max(target[0],1e-6)
    slope=np.sqrt(np.mean((dp-dt[None,:])**2,axis=1))
    cpen=cens.mean(axis=1)
    score=0.35*rel+0.25*log+0.30*shape+0.10*slope+2.0*cpen
    return score,rel,log,shape,slope,cens.sum(axis=1)


def append_scores(cand: pd.DataFrame,Kc: np.ndarray,Tgrid: np.ndarray,targets: pd.DataFrame,stage: str,Kmax: float)->pd.DataFrame:
    blocks=[]
    for klass in TARGET_CLASSES:
        score,rel,log,shape,slope,nc=score_curves(Kc,Tgrid,targets,klass,Kmax)
        q=cand.copy(); q.insert(1,"target_class",klass); q.insert(2,"stage",stage)
        q["score"]=score; q["relative_rmse"]=rel; q["log_rmse"]=log; q["shape_rmse"]=shape; q["slope_rmse"]=slope; q["n_censored"]=nc
        blocks.append(q)
    return pd.concat(blocks,ignore_index=True)


def top_ids(scores: pd.DataFrame,n_per_class:int)->np.ndarray:
    ids=set()
    for k in TARGET_CLASSES:
        ids.update(int(x) for x in scores[scores.target_class==k].nsmallest(n_per_class,"score").candidate_id)
    return np.array(sorted(ids),int)


def perturb_candidates(seeds: pd.DataFrame,n_perturb:int,seed:int)->pd.DataFrame:
    rng=np.random.default_rng(seed); rows=[]
    for _,r in seeds.iterrows():
        base=r.to_dict(); rows.append(base.copy())
        for _ in range(n_perturb):
            x=base.copy()
            x["cleave_G00_eV"]=float(_clip(np.array([base["cleave_G00_eV"]*math.exp(rng.normal(0,0.16))]),"cleave_G00_eV")[0])
            x["cleave_gT_eV_per_K"]=float(_clip(np.array([base["cleave_gT_eV_per_K"]+rng.normal(0,0.0012)]),"cleave_gT_eV_per_K")[0])
            x["cleave_sigc0_GPa"]=float(_clip(np.array([base["cleave_sigc0_GPa"]*math.exp(rng.normal(0,0.16))]),"cleave_sigc0_GPa")[0])
            x["cleave_sT_MPa_per_K"]=float(_clip(np.array([base["cleave_sT_MPa_per_K"]+rng.normal(0,0.8)]),"cleave_sT_MPa_per_K")[0])
            x["cleave_exp_a"]=float(_clip(np.array([base["cleave_exp_a"]*math.exp(rng.normal(0,0.20))]),"cleave_exp_a")[0])
            x["cleave_exp_n"]=float(_clip(np.array([base["cleave_exp_n"]*math.exp(rng.normal(0,0.16))]),"cleave_exp_n")[0])
            x["cleave_floor_frac"]=float(_clip(np.array([base["cleave_floor_frac"]*math.exp(rng.normal(0,0.30))]),"cleave_floor_frac")[0])
            x["cleave_S_hs_kB"]=float(_clip(np.array([base["cleave_S_hs_kB"]+rng.normal(0,7.5)]),"cleave_S_hs_kB")[0])
            x["chi_shield"]=float(_clip(np.array([base["chi_shield"]+rng.normal(0,0.06)]),"chi_shield")[0])
            if math.isfinite(float(base["N_sat"])):
                x["N_sat"]=float(np.clip(base["N_sat"]*math.exp(rng.normal(0,0.30)),*RANGES["N_sat_finite"]))
            elif rng.random()<0.15:
                x["N_sat"]=float(_loguniform(np.array([rng.random()]),*RANGES["N_sat_finite"])[0])
            rows.append(x)
    out=pd.DataFrame(rows)
    # IDs are stage-local; preserve source broad candidate for traceability.
    out=out.drop(columns=[c for c in ["target_class","stage","score","relative_rmse","log_rmse","shape_rmse","slope_rmse","n_censored"] if c in out.columns],errors="ignore")
    if "candidate_id" in out.columns: out=out.rename(columns={"candidate_id":"parent_candidate_id"})
    out.insert(0,"candidate_id",np.arange(len(out),dtype=int))
    return out


def make_plot(curves:pd.DataFrame,targets:pd.DataFrame,validation:pd.DataFrame|None,out:Path):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(12,8.5),constrained_layout=True)
    for ax,k in zip(axes.ravel(),TARGET_CLASSES):
        tg=target_curve(targets,k); cg=curves[curves.target_class==k].sort_values("T_K")
        ax.plot(tg.T_K,tg.target_Kc_MPa_sqrtm,label="prior 1-D target")
        ax.plot(cg.T_K,cg.Kc_MPa_sqrtm,label="tuned EXP-floor V1")
        if validation is not None and not validation.empty:
            vk="dbtt" if k=="DBTT" else k
            vg=validation[validation.regime_key.astype(str)==vk]
            if not vg.empty: ax.scatter(vg.T_K,vg.Kc_2D_MPa_sqrtm,s=18,label="prior 2-D validation")
        ax.set_title(k); ax.set_xlabel("Temperature (K)"); ax.set_ylabel(r"$K_c$ (MPa$\sqrt{m}$)"); ax.grid(alpha=.25); ax.legend(fontsize=8)
    fig.savefig(out/"exp_floor_four_class_tuning.png",dpi=220); plt.close(fig)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--design",default="exp_floor_emission_design_seed.csv")
    ap.add_argument("--targets",default="exp_floor_four_class_target_curves.csv")
    ap.add_argument("--validation-2d",default="exp_floor_four_class_2D_validation.csv")
    ap.add_argument("--out",default="runs/v1_exp_floor_four_class_tuning")
    ap.add_argument("--n-contexts",type=int,default=256)
    ap.add_argument("--seed",type=int,default=20260710)
    ap.add_argument("--Kmax",type=float,default=80.0)
    ap.add_argument("--Kdot",type=float,default=0.005)
    ap.add_argument("--broad-dK",type=float,default=0.25)
    ap.add_argument("--local-dK",type=float,default=0.05)
    ap.add_argument("--final-dK",type=float,default=0.02)
    ap.add_argument("--broad-top-per-class",type=int,default=24)
    ap.add_argument("--local-perturb-per-seed",type=int,default=16)
    ap.add_argument("--final-top-per-class",type=int,default=6)
    ap.add_argument("--resume",action="store_true", help="reuse completed stage tables/arrays in the same OUT directory")
    ap.add_argument("--smoke",action="store_true")
    args=ap.parse_args()

    if args.smoke:
        args.n_contexts=min(args.n_contexts,4); args.broad_top_per_class=min(args.broad_top_per_class,2); args.local_perturb_per_seed=min(args.local_perturb_per_seed,2); args.final_top_per_class=min(args.final_top_per_class,2)
        args.broad_dK=max(args.broad_dK,1.0); args.local_dK=max(args.local_dK,0.5); args.final_dK=max(args.final_dK,0.25)

    out=Path(args.out); out.mkdir(parents=True,exist_ok=True)
    design=pd.read_csv(args.design).copy().reset_index(drop=True)
    design=design.rename(columns={"scenario":"surface_id","design_index":"surface_index"})
    targets=load_targets(args.targets)
    for k in TARGET_CLASSES: target_curve(targets,k)
    validation=None
    vp=Path(args.validation_2d)
    if vp.exists(): validation=pd.read_csv(vp)

    from arrhenius_fracture.config import ElasticProperties
    mat=ElasticProperties()

    Tb=make_temperature_grid(targets,"broad")
    Tl=make_temperature_grid(targets,"local")
    Td=make_temperature_grid(targets,"dense")
    if args.smoke: Tb=np.array([300.,600.,900.,1200.]); Tl=Tb.copy(); Td=Tb.copy()

    contexts=make_broad_contexts(args.n_contexts,args.seed)
    broad=cross_emission_contexts(design,contexts)
    broad.to_csv(out/"broad_candidates.csv",index=False)
    broad_npz=out/"broad_Kc.npz"; broad_scores_path=out/"broad_candidate_scores.csv"
    if args.resume and broad_npz.exists() and broad_scores_path.exists():
        z=np.load(broad_npz); Kb=z["Kc"]; bs=pd.read_csv(broad_scores_path)
        print(f"resume broad stage: {len(broad)} candidates",flush=True)
    else:
        print(f"broad stage: {len(broad)} candidates; {len(Tb)} temperatures; dK={args.broad_dK}",flush=True)
        Kb=simulate_candidates(broad,Tb,Kmax=args.Kmax,dK=args.broad_dK,Kdot=args.Kdot,G_Pa=mat.G,nu=mat.nu,b_m=mat.b)
        np.savez_compressed(broad_npz,Kc=Kb,T=Tb)
        bs=append_scores(broad,Kb,Tb,targets,"broad",args.Kmax); bs.to_csv(broad_scores_path,index=False)

    bids=top_ids(bs,args.broad_top_per_class)
    seeds=broad[broad.candidate_id.isin(bids)].copy().reset_index(drop=True)
    local=perturb_candidates(seeds,args.local_perturb_per_seed,args.seed+1)
    local.to_csv(out/"local_candidates.csv",index=False)
    local_npz=out/"local_Kc.npz"; local_scores_path=out/"local_candidate_scores.csv"
    if args.resume and local_npz.exists() and local_scores_path.exists():
        z=np.load(local_npz); Kl=z["Kc"]; ls=pd.read_csv(local_scores_path)
        print(f"resume local stage: {len(local)} candidates",flush=True)
    else:
        print(f"local stage: {len(local)} candidates; {len(Tl)} temperatures; dK={args.local_dK}",flush=True)
        Kl=simulate_candidates(local,Tl,Kmax=args.Kmax,dK=args.local_dK,Kdot=args.Kdot,G_Pa=mat.G,nu=mat.nu,b_m=mat.b)
        np.savez_compressed(local_npz,Kc=Kl,T=Tl)
        ls=append_scores(local,Kl,Tl,targets,"local",args.Kmax); ls.to_csv(local_scores_path,index=False)

    lids=top_ids(ls,args.final_top_per_class)
    final=local[local.candidate_id.isin(lids)].copy().reset_index(drop=True)
    final_npz=out/"final_Kc.npz"; final_scores_path=out/"finalist_scores.csv"
    if args.resume and final_npz.exists() and final_scores_path.exists():
        z=np.load(final_npz); Kf=z["Kc"]; fs=pd.read_csv(final_scores_path)
        print(f"resume final stage: {len(final)} candidates",flush=True)
    else:
        print(f"final stage: {len(final)} candidates; {len(Td)} temperatures; dK={args.final_dK}",flush=True)
        Kf=simulate_candidates(final,Td,Kmax=args.Kmax,dK=args.final_dK,Kdot=args.Kdot,G_Pa=mat.G,nu=mat.nu,b_m=mat.b)
        np.savez_compressed(final_npz,Kc=Kf,T=Td)
        fs=append_scores(final,Kf,Td,targets,"final",args.Kmax); fs.to_csv(final_scores_path,index=False)

    rec_rows=[]; curve_rows=[]
    for klass in TARGET_CLASSES:
        g=fs[fs.target_class==klass].sort_values("score")
        if g.empty: raise RuntimeError(f"no finalist for {klass}")
        best=g.iloc[0].copy(); rec_rows.append(best)
        cid=int(best.candidate_id); j=int(np.flatnonzero(final.candidate_id.to_numpy(int)==cid)[0])
        for i,T in enumerate(Td):
            curve_rows.append({"target_class":klass,"candidate_id":cid,"T_K":float(T),"Kc_MPa_sqrtm":float(Kf[j,i]) if np.isfinite(Kf[j,i]) else np.nan,"target_Kc_MPa_sqrtm":float(interp_target(targets,klass,np.array([T]))[0])})
    rec=pd.DataFrame(rec_rows)
    front=["target_class","candidate_id","surface_id","surface_index","exp_G00_eV","exp_gT_eV_per_K","exp_sigc0_GPa","exp_sT_MPa_per_K","exp_a","exp_n","exp_floor_frac","cleave_G00_eV","cleave_gT_eV_per_K","cleave_sigc0_GPa","cleave_sT_MPa_per_K","cleave_exp_a","cleave_exp_n","cleave_floor_frac","cleave_S_hs_kB","chi_shield","N_sat","score","relative_rmse","log_rmse","shape_rmse","slope_rmse","n_censored"]
    rec=rec[[c for c in front if c in rec.columns]+[c for c in rec.columns if c not in front]]
    rec.to_csv(out/"recommended_exp_floor_four_class.csv",index=False)
    curves=pd.DataFrame(curve_rows); curves.to_csv(out/"recommended_curves_dense.csv",index=False)
    targets.to_csv(out/"target_curves_used.csv",index=False)
    make_plot(curves,targets,validation,out)

    cfg=vars(args).copy(); cfg.update({"broad_temperatures_K":Tb.tolist(),"local_temperatures_K":Tl.tolist(),"dense_temperatures_K":Td.tolist(),"n_emission_surfaces":len(design),"n_broad_candidates":len(broad),"n_local_candidates":len(local),"n_final_candidates":len(final),"target_role":"observable prior 1-D Kc(T) response; old barrier parameters are not fitted targets"})
    with (out/"run_config.json").open("w") as f: json.dump(cfg,f,indent=2)
    print("\nRecommended fully EXP-floor seeds:")
    print(rec[["target_class","surface_id","cleave_G00_eV","cleave_gT_eV_per_K","cleave_sigc0_GPa","cleave_sT_MPa_per_K","cleave_exp_a","cleave_exp_n","cleave_floor_frac","cleave_S_hs_kB","chi_shield","N_sat","score"]].to_string(index=False))
    print(f"\nWrote {out}")


if __name__=="__main__":
    main()
