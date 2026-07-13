"""Fast batched reduced S-N solver for barrier/phenomenon correlation maps.

This module is a numerical acceleration of :mod:`sn_v1_arrhenius`, not a new
constitutive model.  It advances all stress amplitudes for one barrier surface,
fatigue temperature, and shielding case in one vectorized state batch.

Physics retained from the scalar solver
---------------------------------------
* fully Arrhenius emission -> Peierls -> Taylor series chain;
* Kocks-Mecking rho evolution driven by accepted Arrhenius plastic strain;
* rho feedback through the Taylor barrier amplification, not an athermal stress;
* accumulated plastic strain mapped to bounded shielding/localization states;
* fixed crack-opening barrier sampled after the plastic-state update;
* trapezoidal crack-hazard integration in evolving state;
* adaptive cycle blocks and first passage at B_nuc = 1.

The speedup comes from vectorizing phase quadrature and stress amplitudes and
from constructing the barrier chain only once per surface/temperature batch.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from types import SimpleNamespace
from typing import Dict, Iterable, List, Sequence

import numpy as np
from scipy.special import gammainc

from .sn_arrhenius_chain import ArrheniusPlasticChain, build_chain_from_namespace
from .sn_v1 import KBEV, make_barriers
from .sn_v1_arrhenius import SNCase


@dataclass
class BatchState:
    cycles: np.ndarray
    epsp_acc: np.ndarray
    rho_m2: np.ndarray
    B_nuc: np.ndarray
    failure_cycles: np.ndarray
    blocks_used: np.ndarray


def _effective_multihit_rate(lam_raw, m_hits: float, tau_c_s: float):
    lam_raw = np.asarray(lam_raw, float)
    if m_hits <= 1.0 + 1e-12:
        return lam_raw
    return gammainc(m_hits, np.minimum(lam_raw * max(tau_c_s, 1e-30), 1e12)) / max(tau_c_s, 1e-30)


def _bounded_states(epsp_acc, shield_scale: float, damage_scale: float):
    e = np.maximum(np.asarray(epsp_acc, float), 0.0)
    P = 1.0 - np.exp(-e / max(shield_scale, 1e-30))
    D = 1.0 - np.exp(-e / max(damage_scale, 1e-30))
    return P, D


def _waveform_matrix(sigma_a_MPa: np.ndarray, R: float, n_phase: int):
    s = np.asarray(sigma_a_MPa, float) * 1e6
    sigma_max = 2.0 * s / max(1.0 - R, 1e-30)
    sigma_min = R * sigma_max
    mean = 0.5 * (sigma_max + sigma_min)
    amp = 0.5 * (sigma_max - sigma_min)
    phase = np.linspace(0.0, 2.0 * np.pi, int(n_phase), endpoint=False)
    hist = mean[None, :] + np.cos(phase)[:, None] * amp[None, :]
    return hist, sigma_max, sigma_min


def cycle_quantities_batch(
    args,
    case: SNCase,
    epsp_acc: np.ndarray,
    rho_m2: np.ndarray,
    sigma_a_MPa: np.ndarray,
    chain: ArrheniusPlasticChain,
    crack,
) -> Dict[str, np.ndarray]:
    sigma_a_MPa = np.asarray(sigma_a_MPa, float)
    epsp_acc = np.asarray(epsp_acc, float)
    rho_m2 = np.asarray(rho_m2, float)

    sig_nom, sigmax, sigmin = _waveform_matrix(sigma_a_MPa, args.R, args.n_phase)
    sig_local = args.Kt * sig_nom
    P, Dloc = _bounded_states(epsp_acc, args.epsp_shield_scale, args.epsp_damage_scale)
    sigma_back = args.sigma_back_max_GPa * 1e9 * P

    sig_pl = np.maximum(np.abs(sig_local) - sigma_back[None, :], 0.0)
    cyc = chain.cycle_integrals(sig_pl, rho_m2, args.T, args.frequency_Hz)
    dep_cycle = np.asarray(cyc["dep_eq_per_cycle"], float)

    rho = np.maximum(rho_m2, args.rho_floor)
    drho_cycle = args.k_store * np.sqrt(rho) / chain.b_m * dep_cycle - args.k_dyn * rho * dep_cycle

    sig_nuc = np.maximum(sig_local - case.chi_back * sigma_back[None, :], 0.0)
    G0 = crack.deltaG_eV(sig_nuc, args.T)
    Geff = np.maximum(G0 + case.Gshield_eV * P[None, :] - args.Gstored_eV * Dloc[None, :], 1e-12)
    lam_raw = crack.rate_prefactor * np.exp(np.clip(-Geff / max(KBEV * args.T, 1e-30), -700.0, 0.0))
    lam_eff = _effective_multihit_rate(lam_raw, args.multihit_m, args.multihit_tau_s)
    mu_nuc = np.mean(lam_eff, axis=0) / max(args.frequency_Hz, 1e-300)

    return {
        "dep_eq_per_cycle": dep_cycle,
        "drho_per_cycle": drho_cycle,
        "mu_emit": np.asarray(cyc["mu_emit"], float),
        "mu_peierls": np.asarray(cyc["mu_peierls"], float),
        "mu_taylor": np.asarray(cyc["mu_taylor"], float),
        "mu_escape": np.asarray(cyc["mu_escape"], float),
        "mu_flow": np.asarray(cyc["mu_flow"], float),
        "phi_taylor": np.asarray(cyc["phi_taylor_mean"], float),
        "mu_nuc": np.asarray(mu_nuc, float),
        "P": P,
        "Dloc": Dloc,
        "sigma_back_Pa": sigma_back,
        "sigma_max_nom_Pa": sigmax,
        "sigma_min_nom_Pa": sigmin,
        "G_nuc_min_eV": np.min(Geff, axis=0),
    }


def run_stress_grid(
    args,
    case: SNCase,
    sigma_a_MPa: Sequence[float],
    *,
    chain: ArrheniusPlasticChain | None = None,
) -> List[dict]:
    """Run one case/temperature stress grid in a vectorized adaptive batch."""
    sigma = np.asarray(sorted(set(float(x) for x in sigma_a_MPa)), float)
    n = len(sigma)
    if n == 0:
        return []
    chain = chain or build_chain_from_namespace(args, args.b_m)
    _, crack = make_barriers(-40.0, args.S_crack_kB, args.emit_energy_scale)

    st = BatchState(
        cycles=np.zeros(n, float),
        epsp_acc=np.zeros(n, float),
        rho_m2=np.full(n, float(args.rho0), float),
        B_nuc=np.zeros(n, float),
        failure_cycles=np.full(n, np.nan, float),
        blocks_used=np.zeros(n, int),
    )

    last_q: Dict[str, np.ndarray] | None = None
    # Global iterations are vectorized; each amplitude tracks its own block count.
    while True:
        active = (
            (st.cycles < args.cycles_max * (1.0 - 1e-14))
            & (st.B_nuc < 1.0)
            & (st.blocks_used < args.max_blocks)
        )
        if not np.any(active):
            break
        ia = np.flatnonzero(active)
        q0 = cycle_quantities_batch(
            args, case, st.epsp_acc[ia], st.rho_m2[ia], sigma[ia], chain, crack
        )
        remaining = np.maximum(args.cycles_max - st.cycles[ia], 0.0)
        dN = np.minimum(np.full(len(ia), float(args.block_cycles)), remaining)

        dep = q0["dep_eq_per_cycle"]
        # Limit changes in the bounded crack-coupled plastic states rather than
        # unbounded accumulated strain.  Once P and Dloc saturate, raw epsp can
        # advance in large physical cycle blocks without forcing thousands of
        # numerically unnecessary updates.
        P0 = q0["P"]
        D0 = q0["Dloc"]
        dP_dN = dep * (1.0 - P0) / max(args.epsp_shield_scale, 1e-30)
        dD_dN = dep * (1.0 - D0) / max(args.epsp_damage_scale, 1e-30)
        target_dP = float(getattr(args, "target_dP", 0.03))
        target_dD = float(getattr(args, "target_dD", 0.03))
        good = dP_dN > 0
        dN[good] = np.minimum(dN[good], target_dP / dP_dN[good])
        good = dD_dN > 0
        dN[good] = np.minimum(dN[good], target_dD / dD_dN[good])

        drho_signed = q0["drho_per_cycle"]
        drho = np.abs(drho_signed)
        rho_now = st.rho_m2[ia]
        # A clipped state at rho_cap/floor must not continue to restrict dN in
        # the outward direction; doing so caused the earlier long runs to spend
        # tens of thousands of blocks advancing an already saturated rho state.
        free_rho = ~(((rho_now >= args.rho_cap * (1.0 - 1e-12)) & (drho_signed > 0)) |
                     ((rho_now <= args.rho_floor * (1.0 + 1e-12)) & (drho_signed < 0)))
        good = (drho > 0) & free_rho
        if np.any(good):
            rho_scale = np.maximum(rho_now[good], args.rho0)
            dN[good] = np.minimum(dN[good], args.target_rho_rel_block * rho_scale / drho[good])

        mun = q0["mu_nuc"]
        good = mun > 0
        dN[good] = np.minimum(dN[good], args.target_dB_nuc / mun[good])

        dN = np.maximum(np.minimum(dN, remaining), args.min_block_cycles)
        dN = np.minimum(dN, remaining)
        valid = dN > 0
        if not np.any(valid):
            # Defensive exit for zero remaining-cycle roundoff.
            st.cycles[ia] = np.minimum(st.cycles[ia], args.cycles_max)
            break
        if not np.all(valid):
            ia = ia[valid]
            dN = dN[valid]
            q0 = {k: np.asarray(v)[valid] for k, v in q0.items()}

        eps_old = st.epsp_acc[ia].copy()
        rho_old = st.rho_m2[ia].copy()
        B_old = st.B_nuc[ia].copy()

        def rho_after_deps(rho_old, deps):
            y0 = np.sqrt(np.maximum(rho_old, args.rho_floor))
            if args.k_dyn > 0:
                yss = args.k_store / (chain.b_m * args.k_dyn)
                y = yss + (y0 - yss) * np.exp(-0.5 * args.k_dyn * np.maximum(deps, 0.0))
            else:
                y = y0 + 0.5 * args.k_store / chain.b_m * np.maximum(deps, 0.0)
            return np.clip(y * y, args.rho_floor, args.rho_cap)

        # Predictor-corrector in accumulated Arrhenius plastic strain.  The
        # Kocks-Mecking rho transition is integrated analytically in depsp.
        deps_pred = q0["dep_eq_per_cycle"] * dN
        eps_pred = eps_old + deps_pred
        rho_pred = rho_after_deps(rho_old, deps_pred)
        qpred = cycle_quantities_batch(
            args, case, eps_pred, rho_pred, sigma[ia], chain, crack
        )
        deps_eff = 0.5 * (q0["dep_eq_per_cycle"] + qpred["dep_eq_per_cycle"]) * dN
        st.epsp_acc[ia] = eps_old + deps_eff
        st.rho_m2[ia] = rho_after_deps(rho_old, deps_eff)

        q1 = cycle_quantities_batch(
            args, case, st.epsp_acc[ia], st.rho_m2[ia], sigma[ia], chain, crack
        )
        mu_eff = 0.5 * (q0["mu_nuc"] + q1["mu_nuc"])
        st.B_nuc[ia] += mu_eff * dN
        st.cycles[ia] += dN
        st.blocks_used[ia] += 1

        crossed = np.isnan(st.failure_cycles[ia]) & (st.B_nuc[ia] >= 1.0)
        if np.any(crossed):
            ii = ia[crossed]
            frac = (1.0 - B_old[crossed]) / np.maximum(st.B_nuc[ii] - B_old[crossed], 1e-300)
            st.failure_cycles[ii] = st.cycles[ii] - dN[crossed] + np.clip(frac, 0.0, 1.0) * dN[crossed]
        last_q = q1

    P, D = _bounded_states(st.epsp_acc, args.epsp_shield_scale, args.epsp_damage_scale)
    out: List[dict] = []
    for i, s in enumerate(sigma):
        if np.isfinite(st.failure_cycles[i]):
            status = "failed"
        elif st.cycles[i] >= 0.999999 * args.cycles_max:
            status = "right_censored"
        else:
            status = "block_limited"
        out.append({
            "model": "SN_V1_fully_Arrhenius_batch",
            "case": case.name,
            "sigma_a_MPa": float(s),
            "T_K": float(args.T),
            "R": float(args.R),
            "frequency_Hz": float(args.frequency_Hz),
            "Kt": float(args.Kt),
            "cycles_to_nucleation": float(st.failure_cycles[i]) if np.isfinite(st.failure_cycles[i]) else np.nan,
            "cycles_total": float(st.cycles[i]),
            "status": status,
            "epsp_acc_final": float(st.epsp_acc[i]),
            "rho_final_m2": float(st.rho_m2[i]),
            "P_final": float(P[i]),
            "Dloc_final": float(D[i]),
            "B_nuc_final": float(st.B_nuc[i]),
            "blocks_used": int(st.blocks_used[i]),
            "chi_back": float(case.chi_back),
            "Gshield_eV": float(case.Gshield_eV),
        })
    return out
