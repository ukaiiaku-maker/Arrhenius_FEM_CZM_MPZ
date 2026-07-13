"""Reduced 1-D S-N initiation model using the fully Arrhenius plastic-event chain.

The model starts from a finite elastic stress concentration Kt and no crack.
Plasticity is governed by the same scaled EXP-floor barriers used in the fatigue
framework:

    surface emission -> Peierls glide -> Taylor junction depinning.

All three steps are evaluated as Arrhenius rates.  Peierls and Taylor combine by
series residence times, and the emission step is in series with that mobility
branch.  No quasi-static Peierls stress, Taylor stress, athermal floor, or hard
yield gate is used.

The completed chain produces equivalent plastic strain; rho evolves from the
accepted Arrhenius strain through a Kocks-Mecking state-transition law.  Rho
feeds back only through the Taylor Arrhenius amplification phi_T(rho).  The
accumulated plastic state then alters the crack-nucleation hazard, allowing a
clean comparison between shielded and unshielded S-N initiation.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.special import gammainc

from .sn_arrhenius_chain import build_chain_from_namespace, ArrheniusPlasticChain
from .sn_v1 import make_barriers, waveform_sigma, KBEV


@dataclass
class SNCase:
    name: str
    chi_back: float
    Gshield_eV: float


@dataclass
class State:
    cycles: float = 0.0
    epsp_acc: float = 0.0
    rho_m2: float = 1.0e12
    B_nuc: float = 0.0


def _write_csv(path: Path, rows):
    if not rows:
        return
    keys = sorted(set().union(*(r.keys() for r in rows)))
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader(); w.writerows(rows)


def _effective_multihit_rate(lam_raw, m_hits: float, tau_c_s: float):
    lam_raw = np.asarray(lam_raw, float)
    if m_hits <= 1.0 + 1e-12:
        return lam_raw
    return gammainc(m_hits, np.minimum(lam_raw*max(tau_c_s,1e-30), 1e12)) / max(tau_c_s,1e-30)


def _bounded_states(epsp_acc: float, shield_scale: float, damage_scale: float):
    P = 1.0 - math.exp(-max(epsp_acc,0.0)/max(shield_scale,1e-30))
    D = 1.0 - math.exp(-max(epsp_acc,0.0)/max(damage_scale,1e-30))
    return P, D


def cycle_quantities(args, case: SNCase, state: State,
                     sigma_a_MPa: float, chain: ArrheniusPlasticChain, crack):
    sig_nom, sigmax, sigmin = waveform_sigma(sigma_a_MPa*1e6, args.R, args.n_phase)
    sig_local = args.Kt * sig_nom
    P, Dloc = _bounded_states(state.epsp_acc, args.epsp_shield_scale, args.epsp_damage_scale)
    sigma_back = args.sigma_back_max_GPa*1e9*P

    # Plastic chain uses equivalent amplitude.  Back stress is an evolving
    # internal stress state, not a quasi-static yield threshold.
    sig_pl = np.maximum(np.abs(sig_local)-sigma_back, 0.0)
    sig_hist = sig_pl[:,None]
    cyc = chain.cycle_integrals(sig_hist, np.array([state.rho_m2]), args.T, args.frequency_Hz)
    dep_cycle = float(cyc["dep_eq_per_cycle"][0])

    # Kocks-Mecking rho transition driven by the Arrhenius plastic increment.
    rho = max(state.rho_m2, args.rho_floor)
    drho_cycle = args.k_store*math.sqrt(rho)/chain.b_m*dep_cycle - args.k_dyn*rho*dep_cycle

    sig_nuc = np.maximum(sig_local - case.chi_back*sigma_back, 0.0)
    G0 = crack.deltaG_eV(sig_nuc, args.T)
    Geff = np.maximum(G0 + case.Gshield_eV*P - args.Gstored_eV*Dloc, 1e-12)
    lam_raw = crack.rate_prefactor*np.exp(np.clip(-Geff/max(KBEV*args.T,1e-30), -700.0, 0.0))
    lam_eff = _effective_multihit_rate(lam_raw, args.multihit_m, args.multihit_tau_s)
    mu_nuc = float(np.mean(lam_eff)/max(args.frequency_Hz,1e-300))

    bd = chain.barrier_diagnostics(np.array([np.max(sig_pl)]), np.array([rho]), args.T)
    return {
        "dep_eq_per_cycle": dep_cycle,
        "drho_per_cycle": drho_cycle,
        "mu_emit": float(cyc["mu_emit"][0]),
        "mu_peierls": float(cyc["mu_peierls"][0]),
        "mu_taylor": float(cyc["mu_taylor"][0]),
        "mu_escape": float(cyc["mu_escape"][0]),
        "mu_flow": float(cyc["mu_flow"][0]),
        "phi_taylor": float(cyc["phi_taylor_mean"][0]),
        "mu_nuc": mu_nuc,
        "P": P, "Dloc": Dloc, "sigma_back_Pa": sigma_back,
        "sigma_max_nom_Pa": sigmax, "sigma_min_nom_Pa": sigmin,
        "G_nuc_min_eV": float(np.min(Geff)),
        "G_emit_eV": float(bd["G_emit_eV"][0]),
        "G_peierls_eV": float(bd["G_peierls_eV"][0]),
        "G_taylor_eV": float(bd["G_taylor_eV"][0]),
    }


def run_point(args, case: SNCase, sigma_a_MPa: float):
    chain = build_chain_from_namespace(args, args.b_m)
    _, crack = make_barriers(-40.0, args.S_crack_kB, args.emit_energy_scale)
    st = State(rho_m2=args.rho0)
    rows=[]; failure=None

    for ib in range(args.max_blocks):
        if st.cycles >= args.cycles_max or st.B_nuc >= 1.0:
            break
        q0 = cycle_quantities(args, case, st, sigma_a_MPa, chain, crack)
        remaining = args.cycles_max-st.cycles
        dN = min(args.block_cycles, remaining)
        if q0["dep_eq_per_cycle"] > 0:
            dN=min(dN, args.target_dep_eq_block/q0["dep_eq_per_cycle"])
        if abs(q0["drho_per_cycle"])>0:
            dN=min(dN, args.target_rho_rel_block*max(st.rho_m2,args.rho0)/abs(q0["drho_per_cycle"]))
        if q0["mu_nuc"]>0:
            dN=min(dN, args.target_dB_nuc/q0["mu_nuc"])
        dN=max(min(dN,remaining),args.min_block_cycles); dN=min(dN,remaining)
        if dN<=0: break

        epsp_old=st.epsp_acc; rho_old=st.rho_m2; B_old=st.B_nuc
        st.epsp_acc += q0["dep_eq_per_cycle"]*dN
        st.rho_m2 = float(np.clip(st.rho_m2 + q0["drho_per_cycle"]*dN,
                                  args.rho_floor,args.rho_cap))
        q1 = cycle_quantities(args, case, st, sigma_a_MPa, chain, crack)
        mu_nuc_eff=0.5*(q0["mu_nuc"]+q1["mu_nuc"])
        st.B_nuc += mu_nuc_eff*dN
        st.cycles += dN
        if failure is None and st.B_nuc>=1.0:
            frac=(1.0-B_old)/max(st.B_nuc-B_old,1e-300)
            failure=st.cycles-dN+float(np.clip(frac,0.0,1.0))*dN

        rows.append({
            "block":ib,"case":case.name,"sigma_a_MPa":sigma_a_MPa,
            "T_K":args.T,"cycles_total":st.cycles,"dN":dN,
            "epsp_acc":st.epsp_acc,"rho_m2":st.rho_m2,"B_nuc":st.B_nuc,
            "depsp_block":st.epsp_acc-epsp_old,"drho_block":st.rho_m2-rho_old,
            **q1,
        })

    status="failed" if failure is not None else ("right_censored" if st.cycles>=0.999999*args.cycles_max else "block_limited")
    P,D=_bounded_states(st.epsp_acc,args.epsp_shield_scale,args.epsp_damage_scale)
    return {
        "model":"SN_V1_fully_Arrhenius_emission_Peierls_Taylor",
        "case":case.name,"sigma_a_MPa":sigma_a_MPa,"T_K":args.T,
        "R":args.R,"frequency_Hz":args.frequency_Hz,"Kt":args.Kt,
        "cycles_to_nucleation":failure,"cycles_total":st.cycles,"status":status,
        "epsp_acc_final":st.epsp_acc,"rho_final_m2":st.rho_m2,
        "P_final":P,"Dloc_final":D,"B_nuc_final":st.B_nuc,
        "chi_back":case.chi_back,"Gshield_eV":case.Gshield_eV,
    }, rows


def run_sweep(args):
    out=Path(args.out); out.mkdir(parents=True,exist_ok=True)
    cases=[SNCase("no_shield",0.0,0.0),SNCase("shielded",args.shield_chi,args.Gshield_eV)]
    summaries=[]
    for case in cases:
        cdir=out/case.name; cdir.mkdir(exist_ok=True)
        for s in args.sigma_a_MPa:
            print(f"ARRH-V1 case={case.name} T={args.T:g}K sigma_a={s:g} MPa")
            summary,hist=run_point(args,case,float(s)); summaries.append(summary)
            tag=f"sigmaA_{float(s):g}MPa".replace(".","p")
            _write_csv(cdir/f"history_{tag}.csv",hist)
    _write_csv(out/"sn_v1_arrhenius_summary.csv",summaries)
    with (out/"run_args.json").open("w") as f: json.dump(vars(args),f,indent=2,sort_keys=True)
    return summaries


def add_common_args(p):
    p.add_argument("--T",type=float,default=300.0)
    p.add_argument("--R",type=float,default=0.1)
    p.add_argument("--frequency-Hz",type=float,default=1000.0,dest="frequency_Hz")
    p.add_argument("--Kt",type=float,default=3.0)
    p.add_argument("--b-m",type=float,default=2.74e-10,dest="b_m")
    p.add_argument("--rho0",type=float,default=1e12)
    p.add_argument("--rho-floor",type=float,default=1e8,dest="rho_floor")
    p.add_argument("--rho-cap",type=float,default=1e17,dest="rho_cap")
    p.add_argument("--k-store",type=float,default=np.sqrt(2.0),dest="k_store")
    p.add_argument("--k-dyn",type=float,default=1.0,dest="k_dyn")
    p.add_argument("--exp-system",default="W[100]",dest="exp_system")
    p.add_argument("--exp-a",type=float,default=None,dest="exp_a")
    p.add_argument("--exp-n",type=float,default=None,dest="exp_n")
    p.add_argument("--emit-energy-scale",type=float,default=0.75,dest="emit_energy_scale")
    p.add_argument("--emit-entropy-scale",type=float,default=0.75,dest="emit_entropy_scale")
    p.add_argument("--emit-stress-scale",type=float,default=1.0,dest="emit_stress_scale")
    p.add_argument("--peierls-energy-scale",type=float,default=0.00375,dest="peierls_energy_scale")
    p.add_argument("--peierls-entropy-scale",type=float,default=0.00375,dest="peierls_entropy_scale")
    p.add_argument("--peierls-stress-scale",type=float,default=1.0,dest="peierls_stress_scale")
    p.add_argument("--taylor-energy-scale",type=float,default=0.015,dest="taylor_energy_scale")
    p.add_argument("--taylor-entropy-scale",type=float,default=0.015,dest="taylor_entropy_scale")
    p.add_argument("--taylor-stress-scale",type=float,default=1.0,dest="taylor_stress_scale")
    p.add_argument("--nu0-emit-pz",type=float,default=1e11,dest="nu0_emit_pz")
    p.add_argument("--nu0-peierls",type=float,default=1e11,dest="nu0_peierls")
    p.add_argument("--nu0-taylor",type=float,default=1e11,dest="nu0_taylor")
    p.add_argument("--plastic-event-strain",type=float,default=1e-5,dest="plastic_event_strain")
    p.add_argument("--phi-taylor-max",type=float,default=20.0,dest="phi_taylor_max")
    p.add_argument("--epsp-shield-scale",type=float,default=5e-3,dest="epsp_shield_scale")
    p.add_argument("--epsp-damage-scale",type=float,default=2e-2,dest="epsp_damage_scale")
    p.add_argument("--sigma-back-max-GPa",type=float,default=1.0,dest="sigma_back_max_GPa")
    p.add_argument("--S-crack-kB",type=float,default=0.0,dest="S_crack_kB")
    p.add_argument("--Gstored-eV",type=float,default=0.25,dest="Gstored_eV")
    p.add_argument("--Gshield-eV",type=float,default=0.35,dest="Gshield_eV")
    p.add_argument("--shield-chi",type=float,default=0.6,dest="shield_chi")
    p.add_argument("--multihit-m",type=float,default=3.0,dest="multihit_m")
    p.add_argument("--multihit-tau-s",type=float,default=1e-6,dest="multihit_tau_s")
    p.add_argument("--n-phase",type=int,default=128,dest="n_phase")
    return p


def build_parser():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out",default="runs/sn_v1_arrhenius_two_case")
    p.add_argument("--sigma-a-MPa",nargs="+",type=float,
                   default=[250,300,350,400,450,500,550,600,650,700,750,800,900],dest="sigma_a_MPa")
    p.add_argument("--cycles-max",type=float,default=1e10,dest="cycles_max")
    p.add_argument("--block-cycles",type=float,default=1e7,dest="block_cycles")
    p.add_argument("--min-block-cycles",type=float,default=1e-6,dest="min_block_cycles")
    p.add_argument("--max-blocks",type=int,default=5000,dest="max_blocks")
    p.add_argument("--target-dep-eq-block",type=float,default=2e-4,dest="target_dep_eq_block")
    p.add_argument("--target-rho-rel-block",type=float,default=0.05,dest="target_rho_rel_block")
    p.add_argument("--target-dB-nuc",type=float,default=0.05,dest="target_dB_nuc")
    add_common_args(p)
    return p


def main(argv=None):
    args=build_parser().parse_args(argv); run_sweep(args)


if __name__=="__main__": main()
