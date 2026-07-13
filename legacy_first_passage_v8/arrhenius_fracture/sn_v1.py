"""Reduced 1-D S-N crack-initiation model for the Arrhenius fatigue framework.

The model starts from a finite stress concentration Kt representing a blunt
surface feature.  There is no pre-existing crack and no imposed Paris law.

Common topology
---------------
  cyclic nominal stress
    -> cycle-integrated plastic/emission hazard
    -> cumulative local plastic state P and localization/stored-energy state D
    -> crack-nucleation hazard modified by P and D
    -> cumulative first-passage clock B_nuc

Two pilot cases use identical barriers and differ only in the coupling of the
plastic state to crack opening:
  no_shield : G_shield = 0
  shielded  : G_shield > 0 and chi_back > 0

This makes the two-case comparison a clean test of whether a plastic shielding
feedback can create an S-N endurance-like response from the same local event
physics.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Iterable

import numpy as np
from scipy.special import gammainc

from .config import KB, EV_TO_J
from .fatigue_v1 import ExpFloorBarrierParams, ScaledExpFloorBarrier

KBEV = KB / EV_TO_J


@dataclass
class AnchoredExpFloorBarrier:
    """EXP-floor barrier anchored at T0 with an explicit activation entropy."""
    base_scaled: ScaledExpFloorBarrier
    S_kB: float = 0.0
    T0_K: float = 300.0
    rate_prefactor: float = 1.0e11

    def deltaG_eV(self, sigma_Pa, T_K: float):
        G300 = np.asarray(self.base_scaled.deltaG_eV(sigma_Pa, self.T0_K), float)
        G = G300 - (float(T_K) - self.T0_K) * self.S_kB * KBEV
        return np.maximum(G, 1.0e-12)

    def rate(self, sigma_Pa, T_K: float):
        G = self.deltaG_eV(sigma_Pa, T_K)
        return self.rate_prefactor * np.exp(np.clip(-G / max(KBEV * T_K, 1e-30), -700.0, 0.0))


@dataclass
class SNCase:
    name: str
    chi_back: float
    Gshield_eV: float


@dataclass
class SNState:
    cycles: float = 0.0
    B_emit: float = 0.0
    B_nuc: float = 0.0
    P: float = 0.0
    Dloc: float = 0.0


def make_barriers(S_emit_kB=-40.0, S_crack_kB=0.0, emit_energy_scale=0.75):
    # Emission baseline: the case-64-M1 family at 300 K, then re-anchored.
    eb = ExpFloorBarrierParams.preset("W[100]")
    emit_scaled = ScaledExpFloorBarrier(
        base=eb,
        mechanism="emission",
        energy_scale=emit_energy_scale,
        entropy_scale=0.0,
        stress_scale=1.0,
        rate_prefactor=1e11,
    )
    emit = AnchoredExpFloorBarrier(emit_scaled, S_kB=S_emit_kB, rate_prefactor=1e11)

    cb = ExpFloorBarrierParams(
        name="case64_crack_nucleation",
        G00_eV=1.0,
        gT_eV_per_K=0.0,
        sigc0_Pa=3.0e9,
        sT_Pa_per_K=0.0,
        Tref_K=300.0,
        a=0.70,
        n=0.60,
        Gfloor_fraction=0.010,
        Gfloor_min_eV=1e-4,
        Gfloor_max_fraction=0.95,
    )
    crack_scaled = ScaledExpFloorBarrier(
        base=cb,
        mechanism="crack_nucleation",
        energy_scale=1.0,
        entropy_scale=0.0,
        stress_scale=1.0,
        rate_prefactor=1e11,
    )
    crack = AnchoredExpFloorBarrier(crack_scaled, S_kB=S_crack_kB, rate_prefactor=1e11)
    return emit, crack


def waveform_sigma(sigma_a_Pa: float, R: float, n_phase: int):
    # R = sigma_min/sigma_max and sigma_a=(sigma_max-sigma_min)/2.
    sigma_max = 2.0 * sigma_a_Pa / max(1.0 - R, 1e-30)
    sigma_min = R * sigma_max
    mean = 0.5 * (sigma_max + sigma_min)
    amp = 0.5 * (sigma_max - sigma_min)
    phase = np.linspace(0.0, 2.0 * np.pi, n_phase, endpoint=False)
    return mean + amp * np.cos(phase), sigma_max, sigma_min


def effective_multihit_rate(lam_raw, m_hits: float, tau_c_s: float):
    lam_raw = np.asarray(lam_raw, float)
    if m_hits <= 1.0 + 1e-12:
        return lam_raw
    return gammainc(m_hits, np.minimum(lam_raw * max(tau_c_s, 1e-30), 1e12)) / max(tau_c_s, 1e-30)


def cycle_mean_per_cycle(rate_s, frequency_Hz: float):
    return float(np.mean(np.asarray(rate_s, float)) / max(frequency_Hz, 1e-300))


def state_from_B(B_emit: float, B_shield_events: float, B_damage_events: float):
    P = 1.0 - math.exp(-max(B_emit, 0.0) / max(B_shield_events, 1e-30))
    D = 1.0 - math.exp(-max(B_emit, 0.0) / max(B_damage_events, 1e-30))
    return P, D


def cycle_hazards(
    state: SNState,
    case: SNCase,
    sigma_a_MPa: float,
    T_K: float,
    R: float,
    frequency_Hz: float,
    Kt: float,
    emit_barrier,
    crack_barrier,
    n_phase: int,
    sigma_back_max_GPa: float,
    Gstored_eV: float,
    emit_site_multiplicity: float,
    m_hits: float,
    tau_c_s: float,
):
    sig_nom, sigmax, sigmin = waveform_sigma(sigma_a_MPa * 1e6, R, n_phase)
    sig_local = Kt * sig_nom
    sigma_back = sigma_back_max_GPa * 1e9 * state.P

    # Local hardening opposes repeated emission in both cases.
    sig_emit = np.maximum(np.abs(sig_local) - sigma_back, 0.0)
    mu_emit = emit_site_multiplicity * cycle_mean_per_cycle(emit_barrier.rate(sig_emit, T_K), frequency_Hz)

    # Crack opening sees tensile stress only.  Shielded case additionally sees
    # a fraction of the internal back stress and a positive barrier shift.
    sig_nuc = np.maximum(sig_local - case.chi_back * sigma_back, 0.0)
    G0 = crack_barrier.deltaG_eV(sig_nuc, T_K)
    Geff = np.maximum(G0 + case.Gshield_eV * state.P - Gstored_eV * state.Dloc, 1e-12)
    lam_raw = crack_barrier.rate_prefactor * np.exp(
        np.clip(-Geff / max(KBEV * T_K, 1e-30), -700.0, 0.0)
    )
    lam_eff = effective_multihit_rate(lam_raw, m_hits, tau_c_s)
    mu_nuc = cycle_mean_per_cycle(lam_eff, frequency_Hz)
    return {
        "mu_emit": mu_emit,
        "mu_nuc": mu_nuc,
        "sigma_max_nom_Pa": sigmax,
        "sigma_min_nom_Pa": sigmin,
        "sigma_back_Pa": sigma_back,
        "G_nuc_min_eV": float(np.min(Geff)),
    }


def run_point(args, case: SNCase, sigma_a_MPa: float):
    emit, crack = make_barriers(args.S_emit_kB, args.S_crack_kB, args.emit_energy_scale)
    st = SNState()
    rows = []
    failure_cycles = None

    for ib in range(args.max_blocks):
        if st.cycles >= args.cycles_max or st.B_nuc >= 1.0:
            break
        hz = cycle_hazards(
            st, case, sigma_a_MPa, args.T, args.R, args.frequency_Hz,
            args.Kt, emit, crack, args.n_phase, args.sigma_back_max_GPa,
            args.Gstored_eV, args.emit_site_multiplicity, args.multihit_m, args.multihit_tau_s,
        )
        remaining = args.cycles_max - st.cycles
        dN = min(args.block_cycles, remaining)
        if hz["mu_emit"] > 0:
            # Limit changes in the *bounded state variables*, not raw event count.
            # Once P or Dloc saturates, large cycle jumps are safe.
            dP_dB = (1.0 - st.P) / max(args.B_shield_events, 1e-30)
            dD_dB = (1.0 - st.Dloc) / max(args.B_damage_events, 1e-30)
            state_sens = max(dP_dB / max(args.target_dP, 1e-30),
                             dD_dB / max(args.target_dD, 1e-30))
            if state_sens > 0:
                dN = min(dN, 1.0 / (hz["mu_emit"] * state_sens))
        if hz["mu_nuc"] > 0:
            dN = min(dN, args.target_dB_nuc / hz["mu_nuc"])
        dN = max(min(dN, remaining), args.min_block_cycles)
        dN = min(dN, remaining)
        if dN <= 0:
            break

        st.B_emit += hz["mu_emit"] * dN
        st.P, st.Dloc = state_from_B(st.B_emit, args.B_shield_events, args.B_damage_events)
        # Re-evaluate crack hazard at the post-plastic state; trapezoid in state.
        hz2 = cycle_hazards(
            st, case, sigma_a_MPa, args.T, args.R, args.frequency_Hz,
            args.Kt, emit, crack, args.n_phase, args.sigma_back_max_GPa,
            args.Gstored_eV, args.emit_site_multiplicity, args.multihit_m, args.multihit_tau_s,
        )
        mu_nuc_eff = 0.5 * (hz["mu_nuc"] + hz2["mu_nuc"])
        Bprev = st.B_nuc
        st.B_nuc += mu_nuc_eff * dN
        st.cycles += dN
        if failure_cycles is None and st.B_nuc >= 1.0:
            frac = (1.0 - Bprev) / max(st.B_nuc - Bprev, 1e-300)
            failure_cycles = st.cycles - dN + frac * dN

        rows.append({
            "block": ib,
            "case": case.name,
            "sigma_a_MPa": sigma_a_MPa,
            "T_K": args.T,
            "cycles_total": st.cycles,
            "dN": dN,
            "B_emit": st.B_emit,
            "P_shield": st.P,
            "D_localization": st.Dloc,
            "B_nuc": st.B_nuc,
            **hz2,
        })

    if failure_cycles is not None:
        status = "failed"
    elif st.cycles >= 0.999999 * args.cycles_max:
        status = "right_censored"
    else:
        status = "block_limited"
    summary = {
        "model": "SN_V1_blunt_feature",
        "case": case.name,
        "sigma_a_MPa": sigma_a_MPa,
        "T_K": args.T,
        "R": args.R,
        "frequency_Hz": args.frequency_Hz,
        "Kt": args.Kt,
        "cycles_to_nucleation": failure_cycles,
        "cycles_total": st.cycles,
        "status": status,
        "B_emit_final": st.B_emit,
        "P_shield_final": st.P,
        "D_localization_final": st.Dloc,
        "B_nuc_final": st.B_nuc,
        "chi_back": case.chi_back,
        "Gshield_eV": case.Gshield_eV,
    }
    return summary, rows


def write_csv(path: Path, rows):
    if not rows:
        return
    keys = sorted(set().union(*(r.keys() for r in rows)))
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader(); w.writerows(rows)


def run_sweep(args):
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    cases = [
        SNCase("no_shield", chi_back=0.0, Gshield_eV=0.0),
        SNCase("shielded", chi_back=args.shield_chi, Gshield_eV=args.Gshield_eV),
    ]
    all_summary = []
    for case in cases:
        cdir = out / case.name; cdir.mkdir(exist_ok=True)
        for s in args.sigma_a_MPa:
            print(f"V1 case={case.name} sigma_a={s:g} MPa")
            summary, hist = run_point(args, case, float(s))
            all_summary.append(summary)
            tag = f"sigmaA_{float(s):g}MPa".replace(".", "p")
            write_csv(cdir / f"history_{tag}.csv", hist)
    write_csv(out / "sn_v1_summary.csv", all_summary)
    with (out / "run_args.json").open("w") as f:
        json.dump(vars(args), f, indent=2, sort_keys=True)
    return all_summary


def build_parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", default="runs/sn_v1_two_case")
    p.add_argument("--T", type=float, default=300.0)
    p.add_argument("--sigma-a-MPa", nargs="+", type=float,
                   default=[250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 900], dest="sigma_a_MPa")
    p.add_argument("--R", type=float, default=0.1)
    p.add_argument("--frequency-Hz", type=float, default=1000.0, dest="frequency_Hz")
    p.add_argument("--Kt", type=float, default=3.0)
    p.add_argument("--cycles-max", type=float, default=1e10, dest="cycles_max")
    p.add_argument("--block-cycles", type=float, default=1e7, dest="block_cycles")
    p.add_argument("--min-block-cycles", type=float, default=1e-6, dest="min_block_cycles")
    p.add_argument("--max-blocks", type=int, default=5000, dest="max_blocks")
    p.add_argument("--n-phase", type=int, default=128, dest="n_phase")
    p.add_argument("--target-dP", type=float, default=0.02, dest="target_dP")
    p.add_argument("--target-dD", type=float, default=0.02, dest="target_dD")
    p.add_argument("--target-dB-nuc", type=float, default=0.05, dest="target_dB_nuc")
    p.add_argument("--S-emit-kB", type=float, default=-40.0, dest="S_emit_kB")
    p.add_argument("--emit-energy-scale", type=float, default=0.75, dest="emit_energy_scale")
    p.add_argument("--S-crack-kB", type=float, default=0.0, dest="S_crack_kB")
    p.add_argument("--emit-site-multiplicity", type=float, default=5e8, dest="emit_site_multiplicity")
    p.add_argument("--B-shield-events", type=float, default=50.0, dest="B_shield_events")
    p.add_argument("--B-damage-events", type=float, default=500.0, dest="B_damage_events")
    p.add_argument("--sigma-back-max-GPa", type=float, default=1.0, dest="sigma_back_max_GPa")
    p.add_argument("--Gstored-eV", type=float, default=0.25, dest="Gstored_eV")
    p.add_argument("--Gshield-eV", type=float, default=0.35, dest="Gshield_eV")
    p.add_argument("--shield-chi", type=float, default=0.6, dest="shield_chi")
    p.add_argument("--multihit-m", type=float, default=3.0, dest="multihit_m")
    p.add_argument("--multihit-tau-s", type=float, default=1e-6, dest="multihit_tau_s")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    run_sweep(args)


if __name__ == "__main__":
    main()
