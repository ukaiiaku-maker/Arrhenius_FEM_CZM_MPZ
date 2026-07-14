#!/usr/bin/env python3
"""Broad global search using unified Peierls transport/Taylor retention physics.

Unlike v9.8/v9.9, every v9.10 search starts from a full Sobol population over
one common 25-dimensional domain.  It does not start from the previous
shortlist or from a class-specific down-selection.  Cleavage, emission,
Peierls, Taylor, and process-zone state parameters are optimized together.

The fidelity-0 state uses the same two-population closure as the spatial v9.10
MPZ:

    v_P = jump * lambda_P
    k_enc = eta * v_P * sqrt(rho_f)
    k_release = lambda_T_completion
    k_escape = v_P / L_pz.

Peierls and Taylor microscopic attempt frequencies remain fixed at 1e12 and
1e11 s^-1.  Independent activation entropies provide the effective-prefactor
variation without reopening an eleven-decade prefactor/entropy degeneracy.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.optimize import differential_evolution, minimize

import optimize_mpz_v9_8_joint_response as core
from arrhenius_fracture.config import EV_TO_J, KB
from arrhenius_fracture.emission_derived_plasticity import (
    CorrelatedTaylorConfig,
    EmissionDerivedPeierlsTaylorConfig,
    ExpFloorSurface,
)
from arrhenius_fracture.emission_derived_plasticity_v97 import (
    EmissionDerivedPeierlsTaylorModel,
    IndependentEntropyMechanismScale,
)
from arrhenius_fracture.moving_process_zone_v910 import MovingProcessZoneState


PARAMETER_NAMES = (
    "cleave_G00_eV", "cleave_gT_eV_per_K", "cleave_sigc0_GPa",
    "cleave_sT_GPa_per_K", "cleave_exp_a", "cleave_exp_n", "cleave_floor_frac",
    "emit_G00_eV", "emit_gT_eV_per_K", "emit_sigc0_GPa",
    "emit_sT_GPa_per_K", "emit_exp_a", "emit_exp_n", "emit_floor_frac",
    "peierls_H0_eV", "delta_H_PT_eV",
    "peierls_activation_entropy_kB", "taylor_activation_entropy_kB",
    "log10_taylor_corr_rho_c_m2", "log10_taylor_corr_scale",
    "log10_source_sites_per_system", "log10_encounter_efficiency",
    "log10_retained_recovery_rate_s", "log10_source_refresh_length_um",
    "c_blunt",
)

BOUNDS = {
    "cleave_G00_eV": (0.4, 5.5),
    "cleave_gT_eV_per_K": (-0.004, 0.010),
    "cleave_sigc0_GPa": (0.5, 10.0),
    "cleave_sT_GPa_per_K": (-0.005, 0.005),
    "cleave_exp_a": (0.02, 2.0),
    "cleave_exp_n": (0.45, 2.5),
    "cleave_floor_frac": (0.001, 0.20),
    "emit_G00_eV": (0.4, 5.5),
    "emit_gT_eV_per_K": (-0.004, 0.012),
    "emit_sigc0_GPa": (0.3, 10.0),
    "emit_sT_GPa_per_K": (-0.005, 0.005),
    "emit_exp_a": (0.02, 2.0),
    "emit_exp_n": (0.45, 2.5),
    "emit_floor_frac": (0.001, 0.20),
    "peierls_H0_eV": (0.02, 8.0),
    "delta_H_PT_eV": (0.001, 12.0),
    "peierls_activation_entropy_kB": (-40.0, 40.0),
    "taylor_activation_entropy_kB": (-40.0, 40.0),
    "log10_taylor_corr_rho_c_m2": (10.0, 17.0),
    "log10_taylor_corr_scale": (-1.0, 1.0),
    "log10_source_sites_per_system": (0.0, 4.0),
    "log10_encounter_efficiency": (-3.0, 2.0),
    "log10_retained_recovery_rate_s": (-10.0, 1.0),
    "log10_source_refresh_length_um": (-1.0, 3.0),
    "c_blunt": (0.0, 3.0),
}


def bounds_list() -> list[tuple[float, float]]:
    return [BOUNDS[name] for name in PARAMETER_NAMES]


def decode(x: np.ndarray) -> dict[str, float]:
    p = {name: float(value) for name, value in zip(PARAMETER_NAMES, x)}
    p.update({
        "taylor_H0_eV": p["peierls_H0_eV"] + p["delta_H_PT_eV"],
        "taylor_corr_rho_c_m2": 10.0 ** p["log10_taylor_corr_rho_c_m2"],
        "taylor_corr_scale": 10.0 ** p["log10_taylor_corr_scale"],
        "source_sites_per_system": 10.0 ** p["log10_source_sites_per_system"],
        "encounter_efficiency": 10.0 ** p["log10_encounter_efficiency"],
        "retained_recovery_rate_s": 10.0 ** p["log10_retained_recovery_rate_s"],
        "source_refresh_length_um": 10.0 ** p["log10_source_refresh_length_um"],
        "peierls_nu0_s": 1.0e12,
        "taylor_nu0_s": 1.0e11,
    })
    return p


def build_model(p: dict[str, float], Tref_K: float) -> EmissionDerivedPeierlsTaylorModel:
    emit0 = max(p["emit_G00_eV"], 1.0e-12)
    parent = ExpFloorSurface(
        G00_eV=p["emit_G00_eV"],
        gT_eV_per_K=p["emit_gT_eV_per_K"],
        sigc0_Pa=p["emit_sigc0_GPa"] * 1.0e9,
        sT_Pa_per_K=p["emit_sT_GPa_per_K"] * 1.0e9,
        Tref_K=Tref_K,
        a=p["emit_exp_a"],
        n=p["emit_exp_n"],
        floor_fraction=p["emit_floor_frac"],
        floor_min_eV=1.0e-4,
        floor_max_fraction=0.95,
    )
    return EmissionDerivedPeierlsTaylorModel(
        EmissionDerivedPeierlsTaylorConfig(
            parent=parent,
            peierls=IndependentEntropyMechanismScale(
                p["peierls_H0_eV"] / emit0,
                p["peierls_activation_entropy_kB"],
                1.0,
                p["peierls_nu0_s"],
            ),
            taylor=IndependentEntropyMechanismScale(
                p["taylor_H0_eV"] / emit0,
                p["taylor_activation_entropy_kB"],
                1.0,
                p["taylor_nu0_s"],
            ),
            correlated_taylor=CorrelatedTaylorConfig(
                rho_c_m2=p["taylor_corr_rho_c_m2"],
                renewal_time_s=1.0,
                m_exponent=1.0,
                m_scale=p["taylor_corr_scale"],
                m_cap=float("inf"),
            ),
            mobile_saturation_density_m2=float("inf"),
            mobile_density_floor_m2=0.0,
            jump_length_min_m=0.0,
            taylor_phi_max=float("inf"),
            rate_cap_s=float("inf"),
        )
    )


def barrier_eV(p: dict[str, float], mechanism: str, sigma_Pa: float, T_K: float, Tref_K: float) -> float:
    prefix = "cleave" if mechanism == "cleave" else "emit"
    return float(core.exp_floor_barrier_eV(
        sigma_Pa, T_K,
        p[f"{prefix}_G00_eV"], p[f"{prefix}_gT_eV_per_K"],
        p[f"{prefix}_sigc0_GPa"], p[f"{prefix}_sT_GPa_per_K"],
        p[f"{prefix}_exp_a"], p[f"{prefix}_exp_n"],
        p[f"{prefix}_floor_frac"], Tref_K,
    ))


@dataclass
class ZeroDSettings:
    target_class: str
    temperatures: np.ndarray
    targets: pd.DataFrame
    Kdot: float = 0.005
    dK: float = 0.5
    Kmax: float = 80.0
    target_extension_um: float = 500.0
    da_um: float = 5.0
    r_pz_m: float = 1.0e-6
    L_pz_m: float = 50.0e-6
    b_m: float = 2.74e-10
    G_shear_Pa: float = 160.0e9
    poisson: float = 0.28
    rho0_m2: float = 5.0e12
    nu0_c: float = 1.0e12
    nu0_e: float = 1.0e11
    cleavage_hits: float = 3.0
    cleavage_tau_s: float = 1.0e-6
    Tref_K: float = 481.33
    mobile_shield_fraction: float = 0.0


def simulate_zero_d_rcurve(p: dict[str, float], T_K: float, s: ZeroDSettings) -> dict[str, Any]:
    model = build_model(p, s.Tref_K)
    capacity = 2.0 * p["source_sites_per_system"]
    available = capacity
    mobile = 0.0
    retained = 0.0
    slip = 0.0
    B = 0.0
    K = 0.0
    a_um = 0.0
    events: list[dict[str, float]] = []
    max_events = max(int(round(s.target_extension_um / s.da_um)) + 1, 2)
    k_line = (
        s.G_shear_Pa * s.b_m
        / max(1.0 - s.poisson, 1.0e-9)
        / math.sqrt(2.0 * math.pi * max(0.5 * s.r_pz_m, s.b_m))
        / 1.0e6
    )
    min_pi = float("inf")
    max_shield = 0.0
    last_diag: dict[str, float] = {}

    for _ in range(int(math.ceil(s.Kmax / s.dK))):
        K = min(K + s.dK, s.Kmax)
        dt = s.dK / max(s.Kdot, 1.0e-30)
        Kshield = k_line * (retained + s.mobile_shield_fraction * mobile)
        max_shield = max(max_shield, Kshield)
        r_eff = s.r_pz_m + p["c_blunt"] * s.b_m * max(slip, 0.0)
        sigma = max(K - Kshield, 0.0) * 1.0e6 / math.sqrt(
            2.0 * math.pi * max(r_eff, 1.0e-30)
        )

        Ge = barrier_eV(p, "emit", sigma, T_K, s.Tref_K)
        emit_rate = s.nu0_e * math.exp(
            float(np.clip(-Ge * EV_TO_J / (KB * T_K), -700.0, 0.0))
        )
        emitted = available * (1.0 - math.exp(-min(emit_rate * dt, 700.0)))
        available = max(available - emitted, 0.0)
        mobile += emitted
        slip += emitted

        rho = max(s.rho0_m2 + retained / max(math.pi * s.r_pz_m**2, 1.0e-30), 1.0)
        rates = model.rates(sigma, rho, T_K, s.b_m)
        p_rate = float(np.asarray(rates["peierls_rate_s"]))
        t_rate = float(np.asarray(rates["taylor_completion_rate_s"]))
        jump = float(np.asarray(rates["jump_length_m"]))
        velocity = max(jump * p_rate, 0.0)
        encounter = float(MovingProcessZoneState.encounter_rate_s(
            p_rate, jump, rho, p["encounter_efficiency"]
        ))
        total = mobile + retained
        exchange = encounter + t_rate
        if exchange > 0.0 and total > 0.0:
            r_eq = encounter / exchange * total
            retained = r_eq + (retained - r_eq) * math.exp(-min(exchange * dt, 700.0))
            retained = min(max(retained, 0.0), total)
            mobile = total - retained
        frec = 1.0 - math.exp(-min(p["retained_recovery_rate_s"] * dt, 700.0))
        retained *= 1.0 - frec
        escape_rate = velocity / max(s.L_pz_m, 1.0e-30)
        mobile *= math.exp(-min(escape_rate * dt, 700.0))

        pi_p = velocity * max(K / max(s.Kdot, 1.0e-30), dt) / max(s.r_pz_m, 1.0e-30)
        min_pi = min(min_pi, pi_p)
        Kshield = k_line * (retained + s.mobile_shield_fraction * mobile)
        r_eff = s.r_pz_m + p["c_blunt"] * s.b_m * max(slip, 0.0)
        sigma_c = max(K - Kshield, 0.0) * 1.0e6 / math.sqrt(
            2.0 * math.pi * max(r_eff, 1.0e-30)
        )
        Gc = barrier_eV(p, "cleave", sigma_c, T_K, s.Tref_K)
        lam_c = float(core.cleavage_effective_rate(
            np.asarray([Gc]), np.asarray([T_K]), s.nu0_c,
            s.cleavage_hits, s.cleavage_tau_s,
        )[0])
        B += lam_c * dt
        nfire = min(int(math.floor(max(B, 0.0))), max_events - len(events))
        for _event in range(nfire):
            B -= 1.0
            events.append({
                "a_um": a_um,
                "K_MPa_sqrt_m": K,
                "K_shield_MPa_sqrt_m": Kshield,
                "mobile_count": mobile,
                "retained_count": retained,
                "peierls_rate_s": p_rate,
                "taylor_completion_rate_s": t_rate,
                "encounter_rate_s": encounter,
                "peierls_traverse_number": pi_p,
            })
            a_um += s.da_um
            refresh = 1.0 - math.exp(-s.da_um / max(p["source_refresh_length_um"], 1.0e-12))
            available += (capacity - available) * refresh
            keep = math.exp(-s.da_um * 1.0e-6 / max(s.L_pz_m, 1.0e-30))
            mobile *= keep
            retained *= keep
            slip *= keep
        last_diag = {
            "peierls_rate_s": p_rate,
            "taylor_completion_rate_s": t_rate,
            "encounter_rate_s": encounter,
            "peierls_traverse_number": pi_p,
            "mobile_count": mobile,
            "retained_count": retained,
            "K_shield_MPa_sqrt_m": Kshield,
        }
        if len(events) >= max_events:
            break

    if not events:
        return {
            "completed": False, "n_events": 0,
            "K_init_proxy": float("nan"), "K_plateau_proxy": float("nan"),
            "delta_KR_proxy": float("nan"),
            "early_rise_per_100um_proxy": float("nan"),
            "plateau_rise_per_100um_proxy": float("nan"),
            "min_peierls_traverse_number": min_pi,
            "max_K_shield_MPa_sqrt_m": max_shield,
            **last_diag,
        }
    df = pd.DataFrame(events)
    Kinit = float(df.iloc[0].K_MPa_sqrt_m)
    plateau_start = 0.70 * s.target_extension_um
    plateau = df[df.a_um >= plateau_start]
    Kplateau = float(plateau.K_MPa_sqrt_m.median()) if not plateau.empty else float(df.K_MPa_sqrt_m.iloc[-1])
    early_a = min(100.0, 0.5 * s.target_extension_um)
    early = df[df.a_um <= early_a]
    early_rise = float(early.K_MPa_sqrt_m.iloc[-1] - Kinit) if len(early) >= 2 else 0.0
    late0 = 0.70 * s.target_extension_um
    late1 = 0.90 * s.target_extension_um
    q0 = df.iloc[(df.a_um - late0).abs().argsort()[:1]].K_MPa_sqrt_m.iloc[0]
    q1 = df.iloc[(df.a_um - late1).abs().argsort()[:1]].K_MPa_sqrt_m.iloc[0]
    late_rise = float(q1 - q0)
    return {
        "completed": len(events) >= max_events,
        "n_events": len(events),
        "K_init_proxy": Kinit,
        "K_plateau_proxy": Kplateau,
        "delta_KR_proxy": max(Kplateau - Kinit, 0.0),
        "early_rise_per_100um_proxy": early_rise,
        "plateau_rise_per_100um_proxy": late_rise,
        "min_peierls_traverse_number": float(df.peierls_traverse_number.min()),
        "max_K_shield_MPa_sqrt_m": max(float(df.K_shield_MPa_sqrt_m.max()), max_shield),
        "final_mobile_count": float(df.iloc[-1].mobile_count),
        "final_retained_count": float(df.iloc[-1].retained_count),
        **last_diag,
        "events": events,
    }


class UnifiedObjective:
    def __init__(self, settings: ZeroDSettings):
        self.s = settings
        self.bounds = np.asarray(bounds_list(), dtype=float)

    def evaluate(self, x: np.ndarray, details: bool = False) -> dict[str, Any]:
        x = np.asarray(x, dtype=float)
        if not np.all(np.isfinite(x)):
            return {"objective": 1.0e12, "nonfinite_parameter_vector": True}
        outside = np.maximum(self.bounds[:, 0] - x, 0.0) + np.maximum(x - self.bounds[:, 1], 0.0)
        if np.any(outside > 0.0):
            return {"objective": 1.0e10 + 1.0e7 * float(np.sum(outside**2))}
        p = decode(x)
        model = build_model(p, self.s.Tref_K)
        stress_grid = np.linspace(0.0, 30.0e9, 31)
        order_margin = core.barrier_order_margin_eV(model, self.s.temperatures, stress_grid)
        raw = [model.raw_zero_stress_barrier_eV(m, T) for m in ("peierls", "taylor") for T in self.s.temperatures]
        min_raw = float(np.min(raw))
        if order_margin < -1.0e-9 or min_raw <= 0.0:
            return {
                "objective": 1.0e8 + 1.0e6 * max(-order_margin, 0.0) + 1.0e6 * max(-min_raw, 0.0),
                "barrier_order_margin_eV": order_margin,
                "min_raw_barrier_eV": min_raw,
            }

        rows = []
        event_detail = []
        for T in self.s.temperatures:
            r = simulate_zero_d_rcurve(p, float(T), self.s)
            events = r.pop("events", [])
            rows.append({"T_K": float(T), **r})
            if details:
                for event in events:
                    event_detail.append({"T_K": float(T), **event})
        pred = pd.DataFrame(rows)
        merged = self.s.targets.merge(pred, on="T_K", how="left")
        incomplete = int(np.sum(~pred.completed.astype(bool)))
        components: dict[str, float] = {"completion_loss": 500.0 * incomplete}

        if self.s.target_class in {"ceramic", "weakT"}:
            for name in ("K_init", "K_plateau", "early_rise", "plateau_rise"):
                components[f"{name}_loss"] = 0.0
            components["delta_window_loss"] = 0.0
            for _, row in merged.iterrows():
                w = float(row.weight)
                components["K_init_loss"] += w * core.huber((row.K_init_proxy - row.K_init_target) / max(row.K_init_scale, 1e-9))
                components["K_plateau_loss"] += w * core.huber((row.K_plateau_proxy - row.K_plateau_target) / max(row.K_plateau_scale, 1e-9))
                components["early_rise_loss"] += w * core.huber((row.early_rise_per_100um_proxy - row.early_rise_per_100um_target) / max(row.early_rise_scale, 1e-9))
                components["plateau_rise_loss"] += w * core.huber((row.plateau_rise_per_100um_proxy - row.plateau_rise_per_100um_target) / max(row.plateau_rise_scale, 1e-9))
                if row.delta_KR_proxy < row.delta_KR_min:
                    components["delta_window_loss"] += w * core.huber(row.delta_KR_min - row.delta_KR_proxy)
                elif row.delta_KR_proxy > row.delta_KR_max:
                    components["delta_window_loss"] += w * core.huber(row.delta_KR_proxy - row.delta_KR_max)
        else:
            low = pred[pred.T_K <= 700.0]
            high = pred[pred.T_K >= 900.0]
            Kp_low = float(low.K_plateau_proxy.median())
            Kp_high = float(high.K_plateau_proxy.median())
            Ki_low = float(low.K_init_proxy.median())
            Ki_high = float(high.K_init_proxy.median())
            high_dKR = float(high.delta_KR_proxy.median())
            low_dKR = float(low.delta_KR_proxy.median())
            components.update({
                "DBTT_plateau_rise_loss": core.huber(max(15.0 - (Kp_high - Kp_low), 0.0) / 3.0),
                "DBTT_initiation_rise_loss": core.huber(max(10.0 - (Ki_high - Ki_low), 0.0) / 3.0),
                "DBTT_high_Rcurve_loss": core.huber(max(3.0 - high_dKR, 0.0) / 1.5),
                "DBTT_low_Rcurve_loss": core.huber(max(low_dKR - 3.0, 0.0) / 1.5),
            })

        if self.s.target_class == "weakT":
            min_pi = float(pred.min_peierls_traverse_number.min())
            ratio = p["taylor_H0_eV"] / max(p["peierls_H0_eV"], 1.0e-12)
            components["fast_Peierls_loss"] = 20.0 * core.huber(max(-math.log10(max(min_pi, 1e-300)), 0.0))
            components["FCC_barrier_separation_loss"] = 5.0 * core.huber(max(2.0 - ratio, 0.0))
        objective = float(sum(components.values()))
        if not np.isfinite(objective):
            return {"objective": 1.0e12, "nonfinite_objective_replaced": True}
        result: dict[str, Any] = {
            "objective": objective,
            **components,
            "barrier_order_margin_eV": order_margin,
            "min_raw_barrier_eV": min_raw,
            "min_peierls_traverse_number": float(pred.min_peierls_traverse_number.min()),
            "plateau_temperature_rise": float(pred.sort_values("T_K").K_plateau_proxy.iloc[-1] - pred.sort_values("T_K").K_plateau_proxy.iloc[0]),
            "max_K_shield_MPa_sqrt_m": float(pred.max_K_shield_MPa_sqrt_m.max()),
            "parameters": p,
        }
        if details:
            result["temperature_detail"] = merged.to_dict(orient="records")
            result["event_detail"] = event_detail
        return result

    def __call__(self, x: np.ndarray) -> float:
        return float(self.evaluate(x, details=False)["objective"])


def acceptance(target_class: str, detail: pd.DataFrame, summary: dict[str, Any]) -> tuple[bool, str]:
    if not detail.completed.astype(bool).all():
        return False, "incomplete_zero_d_growth"
    if target_class == "ceramic":
        ok = float(detail.delta_KR_proxy.max()) <= 3.0
        return ok, "ceramic_small_Rcurve" if ok else "ceramic_Rcurve_too_large"
    if target_class == "weakT":
        pi_ok = float(detail.min_peierls_traverse_number.min()) >= 1.0
        flat = float(detail.K_plateau_proxy.max() - detail.K_plateau_proxy.min()) <= 6.0
        dkr = float(detail.delta_KR_proxy.median())
        ratio = float(summary["taylor_H0_eV"]) / max(float(summary["peierls_H0_eV"]), 1e-12)
        ok = pi_ok and flat and 2.0 <= dkr <= 9.0 and ratio >= 2.0
        return ok, "weakT_FCC_fastP_highTbarrier" if ok else "weakT_acceptance_failed"
    low = detail[detail.T_K <= 700.0]
    high = detail[detail.T_K >= 900.0]
    rise = float(high.K_plateau_proxy.median() - low.K_plateau_proxy.median())
    high_dkr = float(high.delta_KR_proxy.median())
    ok = rise >= 15.0 and high_dkr >= 3.0
    return ok, "DBTT_persists_in_unified_zeroD_growth" if ok else "DBTT_unified_trend_insufficient"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target-class", choices=["ceramic", "weakT", "DBTT"], required=True)
    ap.add_argument("--temperatures", default="300 700 900 1200")
    ap.add_argument("--targets", type=Path, default=Path("mpz_three_class_design_targets.csv"))
    ap.add_argument("--restarts", type=int, default=3)
    ap.add_argument("--de-maxiter", type=int, default=60)
    ap.add_argument("--de-popsize", type=int, default=8)
    ap.add_argument("--local-maxiter", type=int, default=250)
    ap.add_argument("--seed", type=int, default=910017)
    ap.add_argument("--Kdot", type=float, default=0.005)
    ap.add_argument("--dK", type=float, default=0.5)
    ap.add_argument("--Kmax", type=float, default=80.0)
    ap.add_argument("--target-extension-um", type=float, default=500.0)
    ap.add_argument("--da-um", type=float, default=5.0)
    ap.add_argument("--population-keep", type=int, default=12)
    ap.add_argument("--shortlist-count", type=int, default=20)
    ap.add_argument("--out", type=Path, default=Path("runs/mpz_v9_10_unified_global_search_v1"))
    a = ap.parse_args()

    temperatures = np.asarray(core.parse_floats(a.temperatures), dtype=float)
    targets = core.load_targets(a.targets, a.target_class, temperatures.tolist())
    settings = ZeroDSettings(
        target_class=a.target_class,
        temperatures=temperatures,
        targets=targets,
        Kdot=a.Kdot,
        dK=a.dK,
        Kmax=a.Kmax,
        target_extension_um=a.target_extension_um,
        da_um=a.da_um,
    )
    objective = UnifiedObjective(settings)
    out = (a.out / a.target_class).resolve()
    checkpoints = out / "checkpoints"
    checkpoints.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    temp_rows: list[dict[str, Any]] = []
    event_rows: list[dict[str, Any]] = []
    history_rows: list[dict[str, Any]] = []

    for restart in range(a.restarts):
        checkpoint = checkpoints / f"restart_{restart:03d}.json"
        if checkpoint.exists():
            saved = json.loads(checkpoint.read_text())
            if saved.get("status") == "COMPLETE":
                rows.extend(saved.get("candidates", []))
                temp_rows.extend(saved.get("temperature_detail", []))
                event_rows.extend(saved.get("event_detail", []))
                history_rows.extend(saved.get("history", []))
                print(f"resumed restart {restart}", flush=True)
                continue
        history: list[dict[str, Any]] = []
        def callback(xk: np.ndarray, convergence: float) -> bool:
            rec = {
                "restart": restart,
                "generation": len(history),
                "objective": objective(xk),
                "convergence": float(convergence),
            }
            history.append(rec)
            print(
                f"class={a.target_class} restart={restart} generation={len(history)} "
                f"objective={rec['objective']:.6g}", flush=True,
            )
            return False
        de = differential_evolution(
            objective,
            bounds_list(),
            maxiter=a.de_maxiter,
            popsize=a.de_popsize,
            seed=a.seed + 1009 * restart,
            init="sobol",
            polish=False,
            updating="immediate",
            workers=1,
            callback=callback,
            tol=1.0e-3,
        )
        local = minimize(
            objective, de.x, method="Powell", bounds=bounds_list(),
            options={"maxiter": a.local_maxiter, "xtol": 1e-4, "ftol": 1e-5},
        )
        candidates = [(float(de.fun), np.asarray(de.x, dtype=float), "de_best")]
        if np.isfinite(local.fun):
            candidates.append((float(local.fun), np.asarray(local.x, dtype=float), "local_best"))
        order = np.argsort(np.asarray(de.population_energies, dtype=float))[:a.population_keep]
        for rank in order:
            candidates.append((float(de.population_energies[rank]), np.asarray(de.population[rank], dtype=float), f"population_{rank}"))
        seen: set[tuple[float, ...]] = set()
        restart_rows = []
        restart_temp = []
        restart_events = []
        for rank, (_, x, source) in enumerate(sorted(candidates, key=lambda z: z[0])):
            key = tuple(np.round(x, 10))
            if key in seen:
                continue
            seen.add(key)
            detail = objective.evaluate(x, details=True)
            p = detail.pop("parameters")
            tdf = pd.DataFrame(detail.pop("temperature_detail"))
            edf = pd.DataFrame(detail.pop("event_detail"))
            candidate_id = f"{a.target_class}_restart{restart:02d}_candidate{rank:02d}"
            summary = {
                "candidate_id": candidate_id,
                "target_class": a.target_class,
                "restart": restart,
                "candidate_source": source,
                "objective": float(detail["objective"]),
                "de_success": bool(de.success),
                "local_success": bool(local.success),
                **{name: float(x[i]) for i, name in enumerate(PARAMETER_NAMES)},
                **{k: float(v) for k, v in p.items() if k not in PARAMETER_NAMES},
                **{k: v for k, v in detail.items() if k != "objective"},
                "search_initialization": "FULL_SOBOL_NO_PRIOR_SHORTLIST",
                "status": "V9_10_UNIFIED_ZEROD_REQUIRES_SPATIAL_PROMOTION",
            }
            accepted, reason = acceptance(a.target_class, tdf, summary)
            summary["accepted_for_spatial_promotion"] = bool(accepted)
            summary["acceptance_reason"] = reason
            tdf["candidate_id"] = candidate_id
            tdf["target_class"] = a.target_class
            edf["candidate_id"] = candidate_id
            edf["target_class"] = a.target_class
            restart_rows.append(summary)
            restart_temp.extend(tdf.to_dict(orient="records"))
            restart_events.extend(edf.to_dict(orient="records"))
        payload = {
            "status": "COMPLETE", "restart": restart,
            "candidates": restart_rows,
            "temperature_detail": restart_temp,
            "event_detail": restart_events,
            "history": history,
        }
        checkpoint.write_text(json.dumps(payload, indent=2, allow_nan=True))
        rows.extend(restart_rows)
        temp_rows.extend(restart_temp)
        event_rows.extend(restart_events)
        history_rows.extend(history)

    results = pd.DataFrame(rows).sort_values("objective").drop_duplicates(PARAMETER_NAMES).reset_index(drop=True)
    accepted = results[results.accepted_for_spatial_promotion.astype(bool)].copy()
    pool = accepted if not accepted.empty else results
    shortlist = pool.head(a.shortlist_count).reset_index(drop=True)
    pd.DataFrame(temp_rows).to_csv(out / "unified_global_temperature_detail.csv", index=False)
    pd.DataFrame(event_rows).to_csv(out / "unified_global_event_detail.csv", index=False)
    pd.DataFrame(history_rows).to_csv(out / "unified_global_generation_history.csv", index=False)
    results.to_csv(out / "unified_global_all_candidates.csv", index=False)
    accepted.to_csv(out / "unified_global_accepted.csv", index=False)
    shortlist.to_csv(out / "unified_global_shortlist.csv", index=False)
    shortlist.to_csv(out / "spatial_promotion_manifest.csv", index=False)
    report = {
        "target_class": a.target_class,
        "n_candidates": int(len(results)),
        "n_accepted": int(len(accepted)),
        "best_objective": float(results.iloc[0].objective),
        "full_search_space": True,
        "prior_shortlist_used": False,
        "parameter_count": len(PARAMETER_NAMES),
        "output": str(out),
        "status": "V9_10_UNIFIED_GLOBAL_SEARCH_COMPLETE",
    }
    (out / "unified_global_summary.json").write_text(json.dumps(report, indent=2))
    config = vars(a).copy()
    for key, value in list(config.items()):
        if isinstance(value, Path):
            config[key] = str(value)
    config.update({"parameter_names": PARAMETER_NAMES, "bounds": BOUNDS})
    (out / "unified_global_config.json").write_text(json.dumps(config, indent=2))
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
