#!/usr/bin/env python3
"""Test categorical strength-anomaly/non-softening vs S-N endurance association.

Unlike the earlier pre-cracked DeltaK study, this calculation uses the blunt-site
S-N initiation model with evolving Arrhenius plastic strain and rho.  For each
barrier scenario it computes:

1. monotonic fixed-rate strength sigma_y(T) by inverting the same fully
   Arrhenius emission -> Peierls -> Taylor chain at rho0;
2. a 300 K S-N initiation curve from the same barrier scenario;
3. categorical labels for an above-room-temperature anomaly/non-softening
   regime and for an endurance-like high-cycle S-N plateau;
4. a contingency table, odds ratio, and Fisher exact test.

The point is deliberately categorical: the hypothesis is that a non-softening
or anomalous strength regime above the fatigue-test temperature is associated
with the *existence* of endurance behavior, not that anomaly amplitude must be
linearly proportional to endurance stress.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.optimize import brentq
from scipy.stats import fisher_exact

from arrhenius_fracture.sn_arrhenius_chain import build_chain_from_namespace
from arrhenius_fracture.sn_v1_arrhenius import build_parser as sn_parser, run_point, SNCase


def write_csv(path, rows):
    if not rows: return
    keys=sorted(set().union(*(r.keys() for r in rows)))
    with open(path,"w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=keys); w.writeheader(); w.writerows(rows)


def strength_at_T(ns, T, epsdot_target, sigma_hi_GPa=20.0):
    chain=build_chain_from_namespace(ns, ns.b_m)
    rho=np.array([ns.rho0],float)
    def rate(sig):
        return float(chain.rates(np.array([sig]),rho,T)["dot_ep"][0])
    r0=rate(0.0)
    if r0>=epsdot_target:
        return 0.0
    hi=sigma_hi_GPa*1e9
    rhi=rate(hi)
    if rhi<epsdot_target:
        return float("nan")
    root=brentq(lambda s: math.log(max(rate(s),1e-300))-math.log(epsdot_target),0.0,hi,maxiter=200)
    return root/1e6


def classify_strength(T, sig, plateau_tol_frac_per_100K=0.02, min_window_K=150.0, anomaly_amp_frac=0.02):
    T=np.asarray(T,float); sig=np.asarray(sig,float)
    good=np.isfinite(sig); T=T[good]; sig=sig[good]
    if len(T)<3:
        return {"has_anomaly":False,"has_plateau":False,"has_nonsoftening":False,
                "anomaly_amp_frac":float("nan"),"max_slope_MPa_per_K":float("nan")}
    order=np.argsort(T); T=T[order]; sig=sig[order]
    i0=int(np.argmin(np.abs(T-300.0))); s300=max(sig[i0],1e-12)
    amp=(np.nanmax(sig[T>=T[i0]])-s300)/s300
    slopes=np.diff(sig)/np.maximum(np.diff(T),1e-30)
    has_anom=bool(amp>=anomaly_amp_frac or np.any(slopes[T[:-1]>=T[i0]]>0.0))
    has_plateau=False
    for i in range(i0,len(T)-1):
        for j in range(i+1,len(T)):
            if T[j]-T[i] < min_window_K: continue
            frac_drop=(sig[i]-sig[j])/max(sig[i],1e-12)
            allowed=plateau_tol_frac_per_100K*((T[j]-T[i])/100.0)
            if frac_drop <= allowed:
                has_plateau=True; break
        if has_plateau: break
    return {"has_anomaly":has_anom,"has_plateau":has_plateau,
            "has_nonsoftening":bool(has_anom or has_plateau),
            "anomaly_amp_frac":float(amp),"max_slope_MPa_per_K":float(np.nanmax(slopes))}


def classify_endurance(points, plateau_tol_frac_per_decade=0.05, min_logN=6.0):
    # Treat right-censored life as a lower bound at cycles_total.  The plateau
    # label requires a slow stress change over at least one full high-cycle decade.
    rows=[]
    for r in points:
        N=r["cycles_to_nucleation"] if r["cycles_to_nucleation"] is not None else r["cycles_total"]
        if N is None or not np.isfinite(N) or N<=0: continue
        rows.append((math.log10(float(N)),float(r["sigma_a_MPa"]),r["status"]))
    rows=sorted(rows)
    best=float("inf"); has=False
    for i in range(len(rows)):
        for j in range(i+1,len(rows)):
            dlog=rows[j][0]-rows[i][0]
            if dlog<1.0 or rows[i][0]<min_logN: continue
            slope=abs(rows[j][1]-rows[i][1])/max(0.5*(rows[j][1]+rows[i][1]),1e-12)/dlog
            best=min(best,slope)
            if slope<=plateau_tol_frac_per_decade:
                has=True
    # A bracketed censored/failed transition is an additional endurance-like
    # signal, but not sufficient alone unless it is in the high-cycle regime.
    cens=[r for r in points if r["status"]=="right_censored"]
    fail=[r for r in points if r["status"]=="failed" and r["cycles_to_nucleation"] and r["cycles_to_nucleation"]>=10**min_logN]
    bracket=bool(cens and fail and max(x["sigma_a_MPa"] for x in cens)<max(x["sigma_a_MPa"] for x in fail))
    return {"endurance_like":bool(has or bracket),
            "best_highcycle_frac_stress_change_per_decade":best if np.isfinite(best) else float("nan"),
            "has_censored_failure_bracket":bracket}


def plot_results(out, strength_rows, sn_rows, summary_rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig,ax=plt.subplots(figsize=(8.5,5.8))
    for key,g in _group(strength_rows,["scenario"]):
        ax.plot([r["T_K"] for r in g],[r["strength_MPa"] for r in g],marker="o",lw=1,label=key[0])
    ax.set_xlabel("Temperature (K)"); ax.set_ylabel("Fixed-rate strength (MPa)")
    ax.set_title("Fully Arrhenius plastic-chain strength curves")
    ax.grid(True,alpha=.3); ax.legend(fontsize=6,ncol=2); fig.tight_layout()
    fig.savefig(Path(out)/"strength_temperature_curves.png",dpi=220); plt.close(fig)

    fig,ax=plt.subplots(figsize=(8.5,5.8))
    for key,g in _group(sn_rows,["scenario"]):
        g=sorted(g,key=lambda r:r["sigma_a_MPa"])
        N=[r["cycles_to_nucleation"] if r["cycles_to_nucleation"] is not None else r["cycles_total"] for r in g]
        ax.plot(N,[r["sigma_a_MPa"] for r in g],marker="o",lw=1,label=key[0])
    ax.set_xscale("log"); ax.set_xlabel("Cycles to nucleation / censoring horizon"); ax.set_ylabel("Stress amplitude (MPa)")
    ax.set_title("Blunt-site S-N initiation curves")
    ax.grid(True,which="both",alpha=.3); ax.legend(fontsize=6,ncol=2); fig.tight_layout()
    fig.savefig(Path(out)/"sn_curves_by_barrier_scenario.png",dpi=220); plt.close(fig)

    fig,ax=plt.subplots(figsize=(7.0,5.5))
    for r in summary_rows:
        marker="o" if r["endurance_like"] else "x"
        ax.scatter(r["anomaly_amp_frac"],r["best_highcycle_frac_stress_change_per_decade"],marker=marker,s=55)
        ax.annotate(r["scenario"],(r["anomaly_amp_frac"],r["best_highcycle_frac_stress_change_per_decade"]),fontsize=6)
    ax.set_xlabel("Strength anomaly amplitude / sigma(300 K)")
    ax.set_ylabel("Best high-cycle fractional stress change per decade")
    ax.set_title("Continuous metrics (categorical test reported separately)")
    ax.grid(True,alpha=.3); fig.tight_layout(); fig.savefig(Path(out)/"anomaly_vs_SN_metric.png",dpi=220); plt.close(fig)


def _group(rows,keys):
    d={}
    for r in rows:
        k=tuple(r[x] for x in keys); d.setdefault(k,[]).append(r)
    return d.items()


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out",default="runs/sn_v1_temperature_anomaly_correlation")
    ap.add_argument("--systems",nargs="+",default=["W[100]","Ta[111]","Cu","Al0.7CoCrFeNi-BCC","Al0.7CoCrFeNi-FCC"])
    ap.add_argument("--entropy-multipliers",nargs="+",type=float,default=[0.0,0.5,1.0,1.5])
    ap.add_argument("--temperatures",nargs="+",type=float,default=[300,350,400,450,500,600,700,800,900,1000])
    ap.add_argument("--sn-stresses",nargs="+",type=float,default=[200,250,300,350,400,450,500,550,600,650,700,750,800,900])
    ap.add_argument("--strength-epsdot",type=float,default=1e-4)
    ap.add_argument("--cycles-max",type=float,default=1e10)
    ap.add_argument("--max-blocks",type=int,default=5000)
    ap.add_argument("--sn-n-phase",type=int,default=48)
    ap.add_argument("--sn-block-cycles",type=float,default=1e8)
    ap.add_argument("--sn-target-dep",type=float,default=5e-4)
    ap.add_argument("--sn-target-rho-rel",type=float,default=0.10)
    ap.add_argument("--sn-target-dB",type=float,default=0.10)
    ap.add_argument("--case",choices=["shielded","no_shield"],default="shielded")
    args=ap.parse_args(); out=Path(args.out); out.mkdir(parents=True,exist_ok=True)

    base_ns=sn_parser().parse_args([])
    strength_rows=[]; sn_rows=[]; summary=[]
    for system in args.systems:
        for M in args.entropy_multipliers:
            ns=SimpleNamespace(**vars(base_ns))
            ns.exp_system=system
            ns.emit_entropy_scale=ns.emit_energy_scale*M
            # Preserve the implemented mechanism scaling ratios for entropy.
            ns.peierls_entropy_scale=ns.peierls_energy_scale*M
            ns.taylor_entropy_scale=ns.taylor_energy_scale*M
            ns.cycles_max=args.cycles_max; ns.max_blocks=args.max_blocks; ns.T=300.0
            ns.n_phase=args.sn_n_phase; ns.block_cycles=args.sn_block_cycles
            ns.target_dep_eq_block=args.sn_target_dep
            ns.target_rho_rel_block=args.sn_target_rho_rel
            ns.target_dB_nuc=args.sn_target_dB
            scenario=f"{system}_M{M:g}"
            for T in args.temperatures:
                strength=strength_at_T(ns,float(T),args.strength_epsdot)
                strength_rows.append({"scenario":scenario,"exp_system":system,"entropy_multiplier":M,"T_K":T,"strength_MPa":strength})
            sr=[r for r in strength_rows if r["scenario"]==scenario]
            sc=classify_strength([r["T_K"] for r in sr],[r["strength_MPa"] for r in sr])

            case=SNCase(args.case, ns.shield_chi if args.case=="shielded" else 0.0,
                        ns.Gshield_eV if args.case=="shielded" else 0.0)
            pts=[]
            for s in args.sn_stresses:
                q=SimpleNamespace(**vars(ns)); q.T=300.0
                rec,_=run_point(q,case,float(s)); rec["scenario"]=scenario; rec["exp_system"]=system; rec["entropy_multiplier"]=M
                sn_rows.append(rec); pts.append(rec)
            ec=classify_endurance(pts)
            summary.append({"scenario":scenario,"exp_system":system,"entropy_multiplier":M,**sc,**ec})
            print(scenario, sc, ec)

    # Categorical contingency table: rows nonsoftening yes/no, columns endurance yes/no.
    a=sum(r["has_nonsoftening"] and r["endurance_like"] for r in summary)
    b=sum(r["has_nonsoftening"] and not r["endurance_like"] for r in summary)
    c=sum((not r["has_nonsoftening"]) and r["endurance_like"] for r in summary)
    d=sum((not r["has_nonsoftening"]) and (not r["endurance_like"]) for r in summary)
    odds,p=fisher_exact([[a,b],[c,d]])
    corr=[{"nonsoftening_and_endurance":a,"nonsoftening_no_endurance":b,
           "softening_and_endurance":c,"softening_no_endurance":d,
           "odds_ratio":odds,"fisher_p":p,"n_scenarios":len(summary)}]

    write_csv(out/"strength_temperature_curves.csv",strength_rows)
    write_csv(out/"sn_initiation_points.csv",sn_rows)
    write_csv(out/"scenario_classification.csv",summary)
    write_csv(out/"categorical_correlation.csv",corr)
    plot_results(out,strength_rows,sn_rows,summary)
    with (out/"run_config.json").open("w") as f: json.dump(vars(args),f,indent=2)
    print(corr[0])

if __name__=="__main__": main()
