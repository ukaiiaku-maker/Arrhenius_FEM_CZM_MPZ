#!/usr/bin/env python3
"""Two-barrier temperature map linking monotonic DBTT and fatigue thresholds.

This driver uses the existing 1-D sharp-front fatigue architecture but imposes
clean, 300 K-anchored activation entropies on the two competing crack-tip
channels:

    G_e(sigma,T) = G_e(sigma,T0) - (T-T0) S_e
    G_c(sigma,T) = G_c(sigma,T0) - (T-T0) S_c

with T0=300 K by default.  The 300 K free-energy surfaces are therefore
preserved exactly while only their temperature slopes are changed.

The same barrier pair is then used in two loading histories:

  1. Monotonic K ramp -> Kc(T), a DBTT-like response measure.
  2. Cyclic K(t)      -> da/dN(DeltaK,T) and threshold brackets.

Default "core" scope:
  * all six material-response cases at S_emit=-40 kB and
    S_cleave={-5,0,+5} kB;
  * the plastic-shielded case additionally at S_emit={-30,-50} kB over the
    same cleavage-entropy values, yielding the full 3x3 entropy map.

Place this script in the root of the current fatigue-PF code tree and run it
there.  It imports arrhenius_fracture from that tree.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import asdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from arrhenius_fracture.config import ElasticProperties, FractureBarrier, KB, EV_TO_J
from arrhenius_fracture.sharp_front import FrontConfig, FrontEngine
from arrhenius_fracture.fatigue_v1 import (
    ExpFloorBarrierParams,
    ScaledExpFloorBarrier,
    FatigueWaveform,
    FatigueControllerConfig,
    FatigueCycleHazardController,
)

KB_EV_PER_K = KB / EV_TO_J


def _compat_barrier_diagnostics(barrier, sigma_Pa, T_K: float, b: float = 2.74e-10):
    """Return the barrier-diagnostic dictionary expected by newer sharp_front trees.

    The current fatigue-PF FrontEngine calls ``barrier.diagnostics(...)`` for
    audit fields such as G*, S*, dG*/dsigma and v*.  The anchored wrappers in
    the first map package implemented only G_barrier/S/rate, so the solver
    failed before the first monotonic step.  This helper synthesizes a
    forward-compatible diagnostic dictionary directly from the anchored free
    energy surface.

    The finite-difference stress derivative is taken at fixed T.  Because the
    anchored entropy shift is stress independent, this derivative is exactly
    the derivative of the 300 K base surface, apart from the numerical finite
    difference.
    """
    sigma = np.asarray(sigma_Pa, dtype=float)
    G_J = np.asarray(barrier.G_barrier(sigma, T_K, b), dtype=float)
    S_J_per_K = np.asarray(barrier.S(sigma, T_K), dtype=float)

    # Stable centered finite difference in stress.  Use a relative perturbation
    # with a small absolute floor so sigma=0 is also well-defined.
    ds = np.maximum(np.abs(sigma) * 1.0e-5, 1.0e4)
    sp = sigma + ds
    sm = np.maximum(sigma - ds, 0.0)
    Gp = np.asarray(barrier.G_barrier(sp, T_K, b), dtype=float)
    Gm = np.asarray(barrier.G_barrier(sm, T_K, b), dtype=float)
    denom = np.maximum(sp - sm, 1.0)
    dG_dsigma_J_per_Pa = (Gp - Gm) / denom
    vstar_m3 = np.maximum(-dG_dsigma_J_per_Pa, 0.0)

    G_eV = G_J / EV_TO_J
    S_kB = S_J_per_K / KB
    dGdT_eV_per_K = -S_J_per_K / EV_TO_J
    dG_dsigma_eV_per_GPa = dG_dsigma_J_per_Pa * 1.0e9 / EV_TO_J
    vstar_b3 = vstar_m3 / max(float(b) ** 3, 1.0e-300)
    H_eff_eV = (G_J + float(T_K) * S_J_per_K) / EV_TO_J

    # Include both canonical names and aliases used by several recent audit
    # branches.  Extra keys are harmless; missing keys are not.
    return {
        "G_J": G_J,
        "G_eV": G_eV,
        "DeltaG_J": G_J,
        "DeltaG_eV": G_eV,
        "S_J_per_K": S_J_per_K,
        "S_kB": S_kB,
        "entropy_kB": S_kB,
        "dGdT_eV_per_K": dGdT_eV_per_K,
        "dG_dsigma_J_per_Pa": dG_dsigma_J_per_Pa,
        "dG_dsigma_eV_per_GPa": dG_dsigma_eV_per_GPa,
        "dGcleave_dsigma_eV_per_GPa": dG_dsigma_eV_per_GPa,
        "dGemit_dsigma_eV_per_GPa": dG_dsigma_eV_per_GPa,
        "vstar_m3": vstar_m3,
        "vstar_b3": vstar_b3,
        "H_eff_eV": H_eff_eV,
    }


class AnchoredScaledBarrier:
    """Freeze a ScaledExpFloorBarrier at T_anchor, then add an anchored entropy.

    The base barrier can retain any calibrated 300 K value.  The imposed entropy
    changes only the T-slope:

        DeltaG(T) = DeltaG(T_anchor) - (T-T_anchor) S*.

    S_kB may be negative.  For S_kB<0, the free-energy barrier rises with T.
    """

    def __init__(self, base: ScaledExpFloorBarrier, S_kB: float, T_anchor_K: float = 300.0):
        self.base = base
        self.S_kB = float(S_kB)
        self.T_anchor_K = float(T_anchor_K)
        self.mechanism = getattr(base, "mechanism", "anchored_scaled_barrier")
        self.rate_prefactor = float(getattr(base, "rate_prefactor", 1e11))

    def deltaG_eV(self, sigma_Pa, T_K: float):
        baseG = np.asarray(self.base.deltaG_eV(sigma_Pa, self.T_anchor_K), dtype=float)
        shift_eV = -(float(T_K) - self.T_anchor_K) * self.S_kB * KB_EV_PER_K
        return np.maximum(baseG + shift_eV, 0.0)

    def rate(self, sigma_Pa, T_K: float):
        DG_J = self.deltaG_eV(sigma_Pa, T_K) * EV_TO_J
        exponent = -DG_J / max(KB * float(T_K), 1e-30)
        return self.rate_prefactor * np.exp(np.clip(exponent, -700.0, 0.0))

    def entropy_over_kB_numeric(self, sigma_Pa: float, T_K: float, dT: float = 1.0) -> float:
        return self.S_kB

    def G_barrier(self, sigma, T: float = 0.0, b: float = 2.74e-10):
        return self.deltaG_eV(sigma, T) * EV_TO_J

    # Compatibility aliases used by newer barrier/front audit code.
    def G_J(self, sigma, T: float = 0.0):
        return self.G_barrier(sigma, T)

    def S(self, sigma, T: float = 0.0):
        return np.ones_like(np.asarray(sigma, dtype=float)) * self.S_kB * KB

    def diagnostics(self, sigma_Pa, T_K: float, b: float = 2.74e-10):
        return _compat_barrier_diagnostics(self, sigma_Pa, T_K, b)

    def __getattr__(self, name):
        # Preserve access to calibrated EXP-floor metadata expected by some
        # diagnostics while keeping the anchored wrapper as the active surface.
        return getattr(self.base, name)

    def as_dict(self):
        return {
            "kind": "anchored_scaled_exp_floor",
            "mechanism": self.mechanism,
            "S_kB": self.S_kB,
            "T_anchor_K": self.T_anchor_K,
            "base": self.base.as_dict() if hasattr(self.base, "as_dict") else {},
        }


class AnchoredCleavageBarrier:
    """300 K-anchored cleavage free-energy surface with imposed S_cleave*."""

    def __init__(self, base: FractureBarrier, S_kB: float, T_anchor_K: float = 300.0):
        self.base = base
        self.S_kB = float(S_kB)
        self.T_anchor_K = float(T_anchor_K)

    def G_barrier(self, sigma, T: float = 0.0, b: float = 2.74e-10):
        baseG = np.asarray(self.base.G_barrier(sigma, self.T_anchor_K, b), dtype=float)
        shift_J = -(float(T) - self.T_anchor_K) * self.S_kB * KB
        return np.maximum(baseG + shift_J, 0.0)

    def S(self, sigma, T: float = 0.0):
        return np.ones_like(np.asarray(sigma, dtype=float)) * self.S_kB * KB

    def deltaG_eV(self, sigma_Pa, T_K: float):
        return np.asarray(self.G_barrier(sigma_Pa, T_K), dtype=float) / EV_TO_J

    def G_J(self, sigma, T: float = 0.0):
        return self.G_barrier(sigma, T)

    def rate(self, sigma_Pa, T_K: float, nu0: float = 1.0e12):
        G_J = np.asarray(self.G_barrier(sigma_Pa, T_K), dtype=float)
        exponent = -G_J / max(KB * float(T_K), 1.0e-30)
        return float(nu0) * np.exp(np.clip(exponent, -700.0, 0.0))

    def entropy_over_kB_numeric(self, sigma_Pa: float, T_K: float, dT: float = 1.0) -> float:
        return self.S_kB

    def diagnostics(self, sigma_Pa, T_K: float, b: float = 2.74e-10):
        return _compat_barrier_diagnostics(self, sigma_Pa, T_K, b)

    def __getattr__(self, name):
        return getattr(self.base, name)

    def as_dict(self):
        payload = asdict(self.base) if hasattr(self.base, "__dataclass_fields__") else {}
        return {
            "kind": "anchored_cleavage_exp_floor",
            "S_kB": self.S_kB,
            "T_anchor_K": self.T_anchor_K,
            "base": payload,
        }


def _parse_float_or_inf(value, default=float("inf")) -> float:
    if value is None:
        return float(default)
    if isinstance(value, str):
        t = value.strip().lower()
        if t in {"inf", "+inf", "infinity", "+infinity"}:
            return float("inf")
    try:
        return float(value)
    except Exception:
        return float(default)


def _front_config(case: Dict, da_m: float = 20e-6) -> FrontConfig:
    """Build the native sharp-front state with case-dependent shielding restored.

    The first two-barrier map accidentally hard-coded chi_shield=0 and N_sat=inf
    for every case.  This corrected map reads the native shielding/saturation
    hierarchy from the case table so emission history can actually feed back on
    cleavage through blunting, back-stress shielding, and the retained ledger.
    """
    f = FrontConfig()
    f.r0 = 1e-6
    f.sigma_cap = 30e9
    f.m_hits = 3.0
    f.tau_c = 1e-6
    f.nu0_c = 1e12
    f.nu0_e = 1e11
    f.beta_back = 1.0
    f.c_blunt = float(case.get("c_blunt", 1.0))
    f.L_pz = 1e-6
    f.v_emb_b3 = float(case.get("v_emb_b3", 500.0))
    f.wake_retain = float(case.get("wake_retain", 0.3))
    f.chi_shield = float(case.get("chi_shield", 0.0))
    f.emb_sat_frac = float(case.get("emb_sat_frac", 1.0))
    f.N_sat = _parse_float_or_inf(case.get("N_sat", float("inf")))
    f.recover_k = float(case.get("recover_k", 0.0))
    f.rho0 = 5e12
    f.da = float(da_m)
    return f


def make_cleavage_barrier(case: Dict, S_c_kB: float, T_anchor: float) -> AnchoredCleavageBarrier:
    cb = FractureBarrier()
    cb.barrier_kind = "exp_floor"
    cb.ef_G00_eV = float(case["cleave_G00_eV"])
    cb.ef_gT_eV_per_K = 0.0
    cb.ef_sigc0_Pa = float(case["cleave_sigc0_GPa"]) * 1e9
    cb.ef_sT_Pa_per_K = 0.0
    cb.ef_a = float(case["cleave_exp_a"])
    cb.ef_n = float(case["cleave_exp_n"])
    cb.ef_Tref_K = float(T_anchor)
    cb.ef_floor_frac = float(case["cleave_floor_frac"])
    cb.ef_floor_min_eV = 1e-4
    cb.ef_floor_max_frac = 0.95
    cb.ef_S_hs_kB = 0.0
    return AnchoredCleavageBarrier(cb, S_kB=S_c_kB, T_anchor_K=T_anchor)


def make_plastic_barriers(case: Dict, S_e_kB: float, T_anchor: float):
    base = ExpFloorBarrierParams.preset("W[100]")
    emit0 = ScaledExpFloorBarrier(
        base=base,
        mechanism="crack_tip_dislocation_emission",
        energy_scale=float(case.get("emit_energy_scale", 0.75)),
        entropy_scale=float(case.get("emit_entropy_scale_calibrated", 0.0)),
        stress_scale=1.0,
        rate_prefactor=1e11,
    )
    peierls0 = ScaledExpFloorBarrier(
        base=base,
        mechanism="peierls_glide_escape",
        energy_scale=0.00375,
        entropy_scale=0.0,
        stress_scale=1.0,
        rate_prefactor=1e12,
    )
    taylor0 = ScaledExpFloorBarrier(
        base=base,
        mechanism="taylor_junction_depinning_escape",
        energy_scale=0.015,
        entropy_scale=0.0,
        stress_scale=1.0,
        rate_prefactor=1e11,
    )
    emit = AnchoredScaledBarrier(emit0, S_kB=S_e_kB, T_anchor_K=T_anchor)
    peierls = AnchoredScaledBarrier(peierls0, S_kB=0.0, T_anchor_K=T_anchor)
    taylor = AnchoredScaledBarrier(taylor0, S_kB=0.0, T_anchor_K=T_anchor)
    return emit, peierls, taylor


def make_front_and_controller(case: Dict, S_e_kB: float, S_c_kB: float, T_anchor: float,
                              da_m: float, n_phase: int, target_dB: float,
                              target_dN_store: float):
    mat = ElasticProperties()
    emit, peierls, taylor = make_plastic_barriers(case, S_e_kB, T_anchor)
    cb = make_cleavage_barrier(case, S_c_kB, T_anchor)
    front = FrontEngine(_front_config(case, da_m), cb, emit, mat.G, mat.nu, mat.b)
    cfg = FatigueControllerConfig(
        n_phase=int(n_phase),
        block_cycles=1e5,
        adaptive_cycles=True,
        max_block_cycles=float("inf"),
        min_block_cycles=1e-6,
        target_dB=float(target_dB),
        target_dN_store=float(target_dN_store),
        recovery_per_s=float(case.get("recover_k", 0.0)),
        N_sat=_parse_float_or_inf(case.get("N_sat", float("inf"))),
        storage_model=str(case.get("map_storage_model", "all_retained")),
        fixed_retained_fraction=1.0,
    )
    # Newer controller versions contain extra adaptive limits.  Set them when
    # present without breaking older trees.
    for name, value in [
        ("target_dN_emit", 0.20),
        ("target_dN_mobile", 0.20),
        ("target_dN_escape", float("inf")),
        ("target_dN_peierls", float("inf")),
        ("target_dN_taylor", float("inf")),
    ]:
        if hasattr(cfg, name):
            setattr(cfg, name, value)
    controller = FatigueCycleHazardController(cfg, emit, peierls, taylor)
    return front, controller


def _cycle_cleave_hazard(front, controller, waveform: FatigueWaveform, T_K: float) -> float:
    """Cycle-integrated cleavage hazard at the current front state."""
    phase = controller._phases()
    K = waveform.K_phase(phase)
    dtw = waveform.period_s / len(phase)
    vals = []
    for k in K:
        sig = front.sigma_tip(float(k))
        lam_c, _, _ = front.lambda_cleave(sig, T_K)
        vals.append(max(float(lam_c), 0.0))
    return float(np.sum(vals) * dtw)


def _renew_from_clock(front, dt_block_s: float) -> Dict:
    """Apply the same renewal semantics as the fatigue controller."""
    N_pre = float(front.N_em)
    n_fire = int(np.floor(min(max(float(front.B), 0.0), 1.0e7)))
    fired = n_fire >= 1
    if fired:
        front.B -= n_fire
        retain = float(np.clip(front.f.wake_retain, 0.0, 1.0)) ** n_fire
        front.N_em = N_pre * retain
        front.a_adv += front.f.da * n_fire
        front.n_adv += n_fire
    return {
        "fired": fired,
        "n_fire": n_fire,
        "N_em_pre_renewal": N_pre,
        "N_em_post_renewal": float(front.N_em),
    }


def _exact_saturating_state_update(N0: float, store_per_cycle: float, cycles: float,
                                    N_sat: float) -> float:
    """Exact block update for dN/dNcyc = store*(1-N/Nsat)."""
    N0 = max(float(N0), 0.0)
    q = max(float(store_per_cycle), 0.0)
    cycles = max(float(cycles), 0.0)
    if q <= 0.0 or cycles <= 0.0:
        return N0
    if math.isfinite(N_sat) and N_sat > 0.0:
        N0 = min(N0, N_sat)
        x = q * cycles / N_sat
        if x > 700.0:
            return float(N_sat)
        return float(N_sat - (N_sat - N0) * math.exp(-x))
    return float(N0 + q * cycles)


def _state_limited_cycles(front, store_per_cycle: float, args) -> float:
    """Limit a block by change in normalized shielding state, not raw event count.

    The original map used target_dN_store=0.01 even when N_sat~2000.  That can
    require hundreds of thousands of blocks before the bounded shielding state
    equilibrates.  For finite N_sat we instead limit Delta(N/Nsat), which is the
    physically relevant state resolution.  For unbounded ledgers the original
    absolute-event limit is retained through target_dN_store_unbounded.
    """
    q = max(float(store_per_cycle), 0.0)
    if q <= 0.0:
        return float("inf")
    Nsat = float(front.f.N_sat)
    if math.isfinite(Nsat) and Nsat > 0.0:
        p = float(np.clip(front.N_em / Nsat, 0.0, 1.0))
        rem = max(1.0 - p, 0.0)
        if rem <= args.saturation_tol_fraction:
            return float("inf")
        dp = min(float(args.target_state_fraction), 0.5 * rem)
        frac = min(max(dp / rem, 0.0), 1.0 - 1.0e-15)
        return float(-Nsat * math.log1p(-frac) / q)
    return float(args.target_dN_store_unbounded / q)


def _map_cycle_step(front, controller, waveform: FatigueWaveform, T_K: float,
                    remaining_cycles: float, args) -> Dict:
    """One saturation-aware adaptive fatigue block for the corrected map."""
    pred_pre = controller.integrate_one_cycle(front, waveform, T_K)
    dN = min(float(remaining_cycles), float(args.max_block_cycles))

    if pred_pre.mu_cleave > 0.0:
        dN = min(dN, args.target_dB / pred_pre.mu_cleave)
    dN = min(dN, _state_limited_cycles(front, pred_pre.store_per_cycle, args))
    dN = max(min(dN, float(remaining_cycles)), float(args.min_block_cycles))

    # Trial the bounded plastic-state update and use both pre- and post-state
    # cleavage hazards when enforcing the dB limit.  This catches either
    # shielding-dominated or embrittlement-dominated state evolution.
    N0 = float(front.N_em)
    Nsat = float(front.f.N_sat)
    N_trial = _exact_saturating_state_update(N0, pred_pre.store_per_cycle, dN, Nsat)
    front.N_em = N_trial
    mu_c_post = _cycle_cleave_hazard(front, controller, waveform, T_K)
    pred_post = controller.integrate_one_cycle(front, waveform, T_K)
    front.N_em = N0

    mu_c_lim = max(float(pred_pre.mu_cleave), float(mu_c_post), 0.0)
    if mu_c_lim > 0.0:
        dN = min(dN, args.target_dB / mu_c_lim)
    dN = max(min(dN, float(remaining_cycles)), float(args.min_block_cycles))

    # Recompute the trial endpoint for the accepted cycle count.  Averaging the
    # pre/post storage rates is second-order accurate for the small normalized
    # state increments imposed above.
    N_trial_pre = _exact_saturating_state_update(N0, pred_pre.store_per_cycle, dN, Nsat)
    front.N_em = N_trial_pre
    pred_end = controller.integrate_one_cycle(front, waveform, T_K)
    front.N_em = N0
    q_store = 0.5 * (max(pred_pre.store_per_cycle, 0.0) + max(pred_end.store_per_cycle, 0.0))
    N_new = _exact_saturating_state_update(N0, q_store, dN, Nsat)

    # Optional recovery.  Current corrected cases use zero, but keep the
    # implementation explicit for reproducibility.
    dt_block = dN * waveform.period_s
    krec = max(float(front.f.recover_k), 0.0)
    if krec > 0.0 and dt_block > 0.0:
        N_new *= math.exp(-min(krec * dt_block, 700.0))

    front.N_em = N_new
    mu_emit_avg = 0.5 * (max(pred_pre.mu_emit, 0.0) + max(pred_end.mu_emit, 0.0))
    sig_emit_avg = 0.5 * (pred_pre.avg_sigma_emit_eff + pred_end.avg_sigma_emit_eff)
    front.W_emit += sig_emit_avg * front.b * front.f.L_pz * mu_emit_avg * dN

    mu_c_commit = _cycle_cleave_hazard(front, controller, waveform, T_K)
    dB = mu_c_commit * dN
    front.B += dB
    front.t += dt_block
    front.K_prev = waveform.Kmax
    renew = _renew_from_clock(front, dt_block)

    return {
        "cycles": float(dN),
        "dB_block": float(dB),
        "mu_emit": float(pred_pre.mu_emit),
        "mu_cleave_pre": float(pred_pre.mu_cleave),
        "mu_cleave_post": float(mu_c_commit),
        "store_per_cycle_pre": float(pred_pre.store_per_cycle),
        "N_em": float(front.N_em),
        "B": float(front.B),
        **renew,
    }


def run_fatigue(case: Dict, S_e_kB: float, S_c_kB: float, T: float, Kmax_MPa: float,
                args) -> Dict:
    front, controller = make_front_and_controller(
        case, S_e_kB, S_c_kB, args.T_anchor_K, args.da_m,
        args.n_phase, args.target_dB, args.target_dN_store,
    )
    wave = FatigueWaveform(
        Kmax=float(Kmax_MPa) * 1e6,
        R=args.R,
        frequency_Hz=args.frequency_Hz,
        closure_clip=True,
    )
    pred0 = controller.integrate_one_cycle(front, wave, T)
    log_ratio0 = math.log10(max(pred0.mu_emit, 1e-300) / max(pred0.mu_cleave, 1e-300))

    cycles_done = 0.0
    cycles_first = float("nan")
    blocks = 0
    last = {}
    for ib in range(args.max_blocks):
        if cycles_done >= args.cycles_max or front.n_adv >= args.n_advances:
            break
        remaining = args.cycles_max - cycles_done
        last = _map_cycle_step(front, controller, wave, T, remaining, args)
        cycles = float(last.get("cycles", 0.0))
        if not math.isfinite(cycles) or cycles <= 0.0:
            break
        cycles_done += cycles
        blocks = ib + 1
        if math.isnan(cycles_first) and int(last.get("n_fire", 0)) > 0:
            cycles_first = cycles_done

    measured = front.n_adv > 0 and cycles_done > 0.0
    da_dN = front.a_adv / cycles_done if measured else float("nan")
    ub = args.da_m / cycles_done if (not measured and cycles_done > 0.0) else float("nan")
    if measured:
        status = "measured"
    elif cycles_done >= 0.999 * args.cycles_max:
        status = "censored_cycle_horizon"
    else:
        status = "unresolved_block_limited"

    return {
        "case_label": case["case_label"],
        "response_regime": case.get("response_regime", ""),
        "chi_shield": float(case.get("chi_shield", 0.0)),
        "N_sat": _parse_float_or_inf(case.get("N_sat", float("inf"))),
        "S_emit_kB": S_e_kB,
        "S_cleave_kB": S_c_kB,
        "T_K": T,
        "Kmax_MPa_sqrtm": Kmax_MPa,
        "DeltaK_MPa_sqrtm": (1.0 - args.R) * Kmax_MPa,
        "cycles_total": cycles_done,
        "cycles_to_first_fire": cycles_first,
        "n_adv": int(front.n_adv),
        "a_adv_m": float(front.a_adv),
        "da_dN_m_per_cycle": da_dN,
        "da_dN_upper_bound_m_per_cycle": ub,
        "status": status,
        "direct_lt_1_cycle": math.isfinite(cycles_first) and cycles_first < 1.0,
        "N_em_final": float(front.N_em),
        "sigma_back_Pa": float(front.sigma_back()),
        "dG_emb_eV": float(front.dG_emb() / EV_TO_J),
        "B_final": float(front.B),
        "blocks_completed": blocks,
        "mu_emit_initial_per_cycle": float(pred0.mu_emit),
        "mu_cleave_initial_per_cycle": float(pred0.mu_cleave),
        "log10_mu_emit_over_cleave_initial": log_ratio0,
        "last_dB_block": float(last.get("dB_block", float("nan"))) if last else float("nan"),
    }


def run_monotonic(case: Dict, S_e_kB: float, S_c_kB: float, T: float, args) -> Dict:
    mat = ElasticProperties()
    emit, _, _ = make_plastic_barriers(case, S_e_kB, args.T_anchor_K)
    cb = make_cleavage_barrier(case, S_c_kB, args.T_anchor_K)
    front = FrontEngine(_front_config(case, args.da_m), cb, emit, mat.G, mat.nu, mat.b)

    dK = args.monotonic_dK_MPa * 1e6
    Kdot = args.Kdot_MPa_sqrtm_per_s * 1e6
    dt = dK / max(Kdot, 1e-300)
    Kc = float("nan")
    n_steps = int(math.ceil(args.monotonic_Kmax_MPa / args.monotonic_dK_MPa))
    last = {}
    for i in range(1, n_steps + 1):
        K = min(i * dK, args.monotonic_Kmax_MPa * 1e6)
        last = front.step(K, T, dt)
        if bool(last.get("fired", False)):
            Kc = K / 1e6
            break
    return {
        "case_label": case["case_label"],
        "S_emit_kB": S_e_kB,
        "S_cleave_kB": S_c_kB,
        "T_K": T,
        "Kc_first_MPa_sqrtm": Kc,
        "reached_monotonic_Kmax": not math.isfinite(Kc),
        "N_em_at_end": float(front.N_em),
        "sigma_back_Pa_at_end": float(front.sigma_back()),
        "dG_emb_eV_at_end": float(front.dG_emb() / EV_TO_J),
        "B_at_end": float(front.B),
        "n_adv_at_end": int(front.n_adv),
    }


def read_cases(path: Path) -> List[Dict]:
    df = pd.read_csv(path)
    return df.to_dict("records")


def scenario_pairs(scope: str, emit_common: List[float], cleave_vals: List[float],
                   shield_emit: List[float]) -> List[Tuple[str, float, float, str]]:
    out = []
    # Common map across all six classes.
    for Se in emit_common:
        for Sc in cleave_vals:
            out.append((scenario_label(Se, Sc), Se, Sc, "all6"))
    # Full shielded-case entropy map adds Se values not already covered.
    existing = {(Se, Sc) for _, Se, Sc, group in out if group == "all6"}
    if scope in {"core", "shielded_full", "all6_full"}:
        for Se in shield_emit:
            for Sc in cleave_vals:
                if (Se, Sc) not in existing:
                    out.append((scenario_label(Se, Sc), Se, Sc, "shielded_only"))
    if scope == "all6_full":
        out = []
        vals = sorted(set(emit_common + shield_emit))
        for Se in vals:
            for Sc in cleave_vals:
                out.append((scenario_label(Se, Sc), Se, Sc, "all6"))
    return out


def fmt_entropy(v: float) -> str:
    if v < 0:
        return f"m{abs(v):g}".replace(".", "p")
    if v > 0:
        return f"p{v:g}".replace(".", "p")
    return "0"


def scenario_label(Se: float, Sc: float) -> str:
    return f"Se_{fmt_entropy(Se)}_Sc_{fmt_entropy(Sc)}"


def save_incremental(rows: List[Dict], path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, index=False)


def threshold_summary(points: pd.DataFrame) -> pd.DataFrame:
    """Summarize threshold brackets without treating unresolved points as no-growth.

    Only ``censored_cycle_horizon`` is accepted as a confirmed no-growth upper
    bound.  ``unresolved_block_limited`` is reported separately.  If growth is
    already measured at the lowest tested DeltaK, the threshold is recorded as
    ``below_grid`` and the first measured DeltaK is an upper bound, not a value.
    """
    rows = []
    if points.empty:
        return pd.DataFrame()
    group_cols = ["case_label", "S_emit_kB", "S_cleave_kB", "T_K"]
    for keys, g in points.groupby(group_cols):
        g = g.sort_values("DeltaK_MPa_sqrtm")
        measured = g[g["status"] == "measured"]
        horizon = g[g["status"] == "censored_cycle_horizon"]
        unresolved = g[g["status"] == "unresolved_block_limited"]
        direct = g[g["direct_lt_1_cycle"]]

        grid_min = float(g["DeltaK_MPa_sqrtm"].min())
        grid_max = float(g["DeltaK_MPa_sqrtm"].max())
        lo = float(horizon["DeltaK_MPa_sqrtm"].max()) if len(horizon) else float("nan")
        hi = float(measured["DeltaK_MPa_sqrtm"].min()) if len(measured) else float("nan")

        threshold_class = "unresolved"
        midpoint = float("nan")
        threshold_upper = float("nan")
        threshold_lower = float("nan")

        if math.isfinite(lo) and math.isfinite(hi) and lo < hi:
            # A valid bracket exists only if no unresolved point lies between the
            # confirmed no-growth lower bound and first measured-growth point.
            between = unresolved[(unresolved["DeltaK_MPa_sqrtm"] > lo) &
                                 (unresolved["DeltaK_MPa_sqrtm"] < hi)]
            if len(between) == 0:
                threshold_class = "bracketed"
                threshold_lower = lo
                threshold_upper = hi
                midpoint = 0.5 * (lo + hi)
            else:
                threshold_class = "unresolved_between_bounds"
                threshold_lower = lo
                threshold_upper = hi
        elif math.isfinite(hi) and math.isclose(hi, grid_min, rel_tol=0.0, abs_tol=1e-12) and not math.isfinite(lo):
            threshold_class = "below_grid"
            threshold_upper = hi
        elif math.isfinite(lo) and not math.isfinite(hi):
            threshold_class = "above_grid_or_no_growth"
            threshold_lower = lo
        elif len(unresolved) > 0:
            threshold_class = "unresolved_block_limited"

        rows.append({
            "case_label": keys[0],
            "S_emit_kB": keys[1],
            "S_cleave_kB": keys[2],
            "T_K": keys[3],
            "threshold_class": threshold_class,
            "DeltaK_threshold_lower_MPa_sqrtm": threshold_lower,
            "DeltaK_threshold_upper_MPa_sqrtm": threshold_upper,
            "DeltaK_threshold_midpoint_MPa_sqrtm": midpoint,
            "DeltaK_highest_confirmed_no_growth_MPa_sqrtm": lo,
            "DeltaK_first_measured_MPa_sqrtm": hi,
            "DeltaK_grid_min_MPa_sqrtm": grid_min,
            "DeltaK_grid_max_MPa_sqrtm": grid_max,
            "DeltaK_first_direct_lt_1_cycle_MPa_sqrtm": float(direct.iloc[0]["DeltaK_MPa_sqrtm"]) if len(direct) else float("nan"),
            "n_measured": len(measured),
            "n_horizon_censored": len(horizon),
            "n_unresolved": len(unresolved),
        })
    return pd.DataFrame(rows)


def make_plots(out: Path, fatigue: pd.DataFrame, thr: pd.DataFrame, mono: pd.DataFrame):
    if fatigue.empty:
        return

    # Common Se=-40: bracketed thresholds only.  Below-grid cases are shown as
    # open downward triangles at their upper bound rather than as fake values.
    common = thr[np.isclose(thr["S_emit_kB"], -40.0)]
    for case, g in common.groupby("case_label", sort=False):
        fig, ax = plt.subplots(figsize=(7.4, 5.2))
        for Sc, ss in g.groupby("S_cleave_kB"):
            ss = ss.sort_values("T_K")
            br = ss[ss["threshold_class"] == "bracketed"]
            if not br.empty:
                y = br["DeltaK_threshold_midpoint_MPa_sqrtm"]
                ylo = y - br["DeltaK_threshold_lower_MPa_sqrtm"]
                yhi = br["DeltaK_threshold_upper_MPa_sqrtm"] - y
                ax.errorbar(br["T_K"], y, yerr=np.vstack([ylo, yhi]), marker="o", capsize=3,
                            label=f"S_c={Sc:g} kB")
            below = ss[ss["threshold_class"] == "below_grid"]
            if not below.empty:
                ax.scatter(below["T_K"], below["DeltaK_threshold_upper_MPa_sqrtm"],
                           marker="v", facecolors="none", edgecolors="black", s=52)
            unres = ss[ss["threshold_class"].str.startswith("unresolved", na=False)]
            if not unres.empty:
                ax.scatter(unres["T_K"], np.full(len(unres), ss["DeltaK_grid_min_MPa_sqrtm"].min()),
                           marker="x", s=42)
        ax.set_xlabel("Temperature (K)")
        ax.set_ylabel(r"Fatigue threshold bracket $\Delta K_{th}$ (MPa $\sqrt{m}$)")
        ax.set_title(f"{case}: corrected fatigue threshold vs temperature")
        ax.grid(True, alpha=0.3)
        ax.legend()
        fig.tight_layout()
        fig.savefig(out / f"{case}_fatigue_threshold_vs_T.png", dpi=220)
        plt.close(fig)

    # Monotonic Kc(T) counterpart.
    if not mono.empty:
        common_m = mono[np.isclose(mono["S_emit_kB"], -40.0)]
        for case, g in common_m.groupby("case_label", sort=False):
            fig, ax = plt.subplots(figsize=(7.4, 5.2))
            for Sc, ss in g.groupby("S_cleave_kB"):
                ss = ss.sort_values("T_K")
                ax.plot(ss["T_K"], ss["Kc_first_MPa_sqrtm"], marker="o", label=f"S_c={Sc:g} kB")
            ax.set_xlabel("Temperature (K)")
            ax.set_ylabel(r"Monotonic first-passage $K_c$ (MPa $\sqrt{m}$)")
            ax.set_title(f"{case}: corrected monotonic DBTT-like response")
            ax.grid(True, alpha=0.3)
            ax.legend()
            fig.tight_layout()
            fig.savefig(out / f"{case}_monotonic_Kc_vs_T.png", dpi=220)
            plt.close(fig)

    # Shielded-case entropy maps use only genuinely bracketed values.
    sh = thr[thr["case_label"] == "plastic_shielded_case64_M1"]
    for T, g in sh.groupby("T_K"):
        gb = g[g["threshold_class"] == "bracketed"]
        piv = gb.pivot_table(index="S_emit_kB", columns="S_cleave_kB",
                             values="DeltaK_threshold_midpoint_MPa_sqrtm", aggfunc="first")
        if piv.empty:
            continue
        fig, ax = plt.subplots(figsize=(6.2, 4.8))
        im = ax.imshow(piv.values, aspect="auto", origin="lower")
        ax.set_xticks(range(len(piv.columns)))
        ax.set_xticklabels([f"{x:g}" for x in piv.columns])
        ax.set_yticks(range(len(piv.index)))
        ax.set_yticklabels([f"{x:g}" for x in piv.index])
        ax.set_xlabel(r"Cleavage entropy $S_c^*/k_B$")
        ax.set_ylabel(r"Emission entropy $S_e^*/k_B$")
        ax.set_title(f"Shielded-case bracketed threshold map at {T:g} K")
        fig.colorbar(im, ax=ax, label=r"$\Delta K_{th}$ midpoint (MPa $\sqrt{m}$)")
        fig.tight_layout()
        fig.savefig(out / f"shielded_entropy_threshold_map_T{int(T)}K.png", dpi=220)
        plt.close(fig)

    # Common-hazard link: correlation uses bracketed fatigue thresholds only.
    if mono.empty:
        return
    merged = thr.merge(mono, on=["case_label", "S_emit_kB", "S_cleave_kB", "T_K"], how="inner")
    merged.to_csv(out / "DBTT_fatigue_link_points.csv", index=False)
    m = merged[(merged["threshold_class"] == "bracketed") &
               np.isfinite(merged["Kc_first_MPa_sqrtm"]) &
               np.isfinite(merged["DeltaK_threshold_midpoint_MPa_sqrtm"])]
    if not m.empty:
        fig, ax = plt.subplots(figsize=(6.6, 5.3))
        for case, g in m.groupby("case_label", sort=False):
            ax.scatter(g["Kc_first_MPa_sqrtm"], g["DeltaK_threshold_midpoint_MPa_sqrtm"],
                       label=case, alpha=0.75)
        ax.set_xlabel(r"Monotonic $K_c$ (MPa $\sqrt{m}$)")
        ax.set_ylabel(r"Bracketed fatigue $\Delta K_{th}$ midpoint (MPa $\sqrt{m}$)")
        ax.set_title("Corrected common-hazard DBTT–fatigue threshold link")
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(out / "DBTT_vs_fatigue_threshold_map.png", dpi=220)
        plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--case-table", default="selected_v1_temperature_cases_corrected.csv")
    ap.add_argument("--out", default="runs/v1_two_barrier_dbtt_fatigue_map_corrected")
    ap.add_argument("--scope", choices=["core", "shielded_full", "all6_full"], default="core")
    ap.add_argument("--temperatures", nargs="+", type=float, default=[300, 400, 500, 600, 700, 900])
    ap.add_argument("--Kmax-MPa-sqrt-m", nargs="+", type=float, default=[0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 10, 12, 14])
    ap.add_argument("--common-emission-entropies-kB", nargs="+", type=float, default=[-40.0])
    ap.add_argument("--cleavage-entropies-kB", nargs="+", type=float, default=[-5.0, 0.0, 5.0])
    ap.add_argument("--shielded-emission-entropies-kB", nargs="+", type=float, default=[-30.0, -40.0, -50.0])
    ap.add_argument("--T-anchor-K", type=float, default=300.0)
    ap.add_argument("--R", type=float, default=0.1)
    ap.add_argument("--frequency-Hz", type=float, default=1000.0)
    ap.add_argument("--cycles-max", type=float, default=2e14)
    ap.add_argument("--max-blocks", type=int, default=10000)
    ap.add_argument("--max-block-cycles", type=float, default=float("inf"))
    ap.add_argument("--min-block-cycles", type=float, default=1e-10)
    ap.add_argument("--target-state-fraction", type=float, default=0.01, help="Max change in N_em/N_sat per block for bounded shielding states.")
    ap.add_argument("--saturation-tol-fraction", type=float, default=1e-4, help="Stop state-limiting blocks once remaining capacity is below this fraction.")
    ap.add_argument("--target-dN-store-unbounded", type=float, default=5.0, help="Absolute retained-state increment limit for cases with N_sat=inf.")
    ap.add_argument("--n-advances", type=int, default=5)
    ap.add_argument("--da-m", type=float, default=20e-6)
    ap.add_argument("--n-phase", type=int, default=96)
    ap.add_argument("--target-dB", type=float, default=0.02)
    ap.add_argument("--target-dN-store", type=float, default=0.01)
    ap.add_argument("--Kdot-MPa-sqrtm-per-s", type=float, default=0.005)
    ap.add_argument("--monotonic-Kmax-MPa", type=float, default=40.0)
    ap.add_argument("--monotonic-dK-MPa", type=float, default=0.025)
    ap.add_argument("--resume", action="store_true", help="Resume from existing CORRECTED output CSVs and skip completed task keys.")
    ap.add_argument("--case-filter", nargs="*", default=[], help="Optional exact case labels to run; comma-separated entries are also accepted.")
    ap.add_argument("--skip-monotonic", action="store_true")
    ap.add_argument("--skip-fatigue", action="store_true")
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    cases = read_cases(Path(args.case_table))
    if args.case_filter:
        wanted = set()
        for token in args.case_filter:
            wanted.update(x.strip() for x in str(token).split(",") if x.strip())
        cases = [c for c in cases if str(c["case_label"]) in wanted]
        missing = wanted - {str(c["case_label"]) for c in cases}
        if missing:
            raise SystemExit(f"Unknown case label(s): {sorted(missing)}")
    case_by_label = {c["case_label"]: c for c in cases}
    scenarios = scenario_pairs(args.scope, args.common_emission_entropies_kB,
                               args.cleavage_entropies_kB,
                               args.shielded_emission_entropies_kB)

    fatigue_path = out / "fatigue_paris_points.csv"
    mono_path = out / "monotonic_DBTT_points.csv"
    fatigue_rows: List[Dict] = []
    mono_rows: List[Dict] = []
    if args.resume and fatigue_path.exists():
        fatigue_rows = pd.read_csv(fatigue_path).to_dict("records")
        print(f"resume: loaded {len(fatigue_rows)} fatigue rows")
    if args.resume and mono_path.exists():
        mono_rows = pd.read_csv(mono_path).to_dict("records")
        print(f"resume: loaded {len(mono_rows)} monotonic rows")

    def fkey(rec):
        return (str(rec["case_label"]), round(float(rec["S_emit_kB"]), 9),
                round(float(rec["S_cleave_kB"]), 9), round(float(rec["T_K"]), 9),
                round(float(rec["Kmax_MPa_sqrtm"]), 9))

    def mkey(rec):
        return (str(rec["case_label"]), round(float(rec["S_emit_kB"]), 9),
                round(float(rec["S_cleave_kB"]), 9), round(float(rec["T_K"]), 9))

    done_f = {fkey(r) for r in fatigue_rows}
    done_m = {mkey(r) for r in mono_rows}
    total_fatigue = 0
    total_mono = 0
    for _, _, _, group in scenarios:
        nc = (1 if "plastic_shielded_case64_M1" in case_by_label else 0) if group == "shielded_only" else len(cases)
        total_fatigue += 0 if args.skip_fatigue else nc * len(args.temperatures) * len(args.Kmax_MPa_sqrt_m)
        total_mono += 0 if args.skip_monotonic else nc * len(args.temperatures)
    print(f"planned fatigue tasks: {total_fatigue}; monotonic tasks: {total_mono}")

    i_f = 0
    i_m = 0
    for label, Se, Sc, group in scenarios:
        use_cases = ([case_by_label["plastic_shielded_case64_M1"]] if "plastic_shielded_case64_M1" in case_by_label else []) if group == "shielded_only" else cases
        print(f"\n=== scenario {label}: Se={Se:g} kB, Sc={Sc:g} kB, group={group} ===")
        for case in use_cases:
            for T in args.temperatures:
                if not args.skip_monotonic:
                    i_m += 1
                    key_m = (str(case["case_label"]), round(float(Se), 9), round(float(Sc), 9), round(float(T), 9))
                    if key_m not in done_m:
                        rec = run_monotonic(case, Se, Sc, T, args)
                        rec["scenario"] = label
                        mono_rows.append(rec)
                        done_m.add(key_m)
                        save_incremental(mono_rows, mono_path)
                        print(f"  M {i_m}/{total_mono}: {case['case_label']} T={T:g} K Kc={rec['Kc_first_MPa_sqrtm']}")
                    elif args.resume:
                        print(f"  M {i_m}/{total_mono}: resume skip {case['case_label']} T={T:g}")
                if not args.skip_fatigue:
                    for K in args.Kmax_MPa_sqrt_m:
                        i_f += 1
                        key_f = (str(case["case_label"]), round(float(Se), 9), round(float(Sc), 9), round(float(T), 9), round(float(K), 9))
                        if key_f in done_f:
                            if args.resume and i_f % 100 == 0:
                                print(f"  F {i_f}/{total_fatigue}: resume skip through {case['case_label']} T={T:g} K={K:g}")
                            continue
                        rec = run_fatigue(case, Se, Sc, T, K, args)
                        rec["scenario"] = label
                        fatigue_rows.append(rec)
                        done_f.add(key_f)
                        save_incremental(fatigue_rows, fatigue_path)
                        if i_f % 25 == 0 or rec["direct_lt_1_cycle"]:
                            print(f"  F {i_f}/{total_fatigue}: {case['case_label']} T={T:g} K={K:g} status={rec['status']} N={rec['cycles_total']:.3g}")

    fatigue = pd.DataFrame(fatigue_rows)
    mono = pd.DataFrame(mono_rows)
    if not fatigue.empty:
        thr = threshold_summary(fatigue)
        thr.to_csv(out / "fatigue_threshold_map.csv", index=False)
    else:
        thr = pd.DataFrame()
    make_plots(out, fatigue, thr, mono)

    settings = vars(args).copy()
    settings["scenarios"] = [dict(label=l, S_emit_kB=Se, S_cleave_kB=Sc, group=g) for l, Se, Sc, g in scenarios]
    with (out / "map_settings.json").open("w") as fp:
        json.dump(settings, fp, indent=2)
    print(f"\nWrote map to {out}")


if __name__ == "__main__":
    main()
