#!/usr/bin/env python3
"""Continue v9.8.1 basins toward lower absolute PT barriers.

The successful v9.8.1 basins are used as warm starts.  Peierls and Taylor
reference enthalpies are scaled together so their ratio and ordering are
preserved, while activation entropies, effective prefactors, correlation
parameters, source inventory, recovery, development length, and blunting are
locally reoptimized.

Class handling is deliberately limited:

* ceramic: no Peierls-mobility requirement;
* weakT: Peierls motion must traverse a process-zone scale during loading;
* DBTT: fit the shelf/transition trends rather than the exact provisional curve.

No assumption is made about the controlling mechanism in ceramic or DBTT.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.optimize import minimize

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


LOCAL_NAMES = (
    "peierls_activation_entropy_kB",
    "taylor_activation_entropy_kB",
    "log10_peierls_nu0_s",
    "log10_taylor_nu0_s",
    "log10_taylor_corr_rho_c_m2",
    "log10_taylor_corr_scale",
    "log10_mobile_fraction",
    "log10_source_sites_per_system",
    "log10_recovery_rate_s",
    "log10_source_refresh_length_um",
    "c_blunt",
)

LOCAL_BOUNDS = {
    "peierls_activation_entropy_kB": (-80.0, 40.0),
    "taylor_activation_entropy_kB": (-80.0, 40.0),
    "log10_peierls_nu0_s": (3.0, 14.0),
    "log10_taylor_nu0_s": (3.0, 14.0),
    "log10_taylor_corr_rho_c_m2": (10.0, 17.0),
    "log10_taylor_corr_scale": (-1.0, 1.0),
    "log10_mobile_fraction": (-5.0, -0.3),
    "log10_source_sites_per_system": (0.0, 4.0),
    "log10_recovery_rate_s": (-10.0, 1.0),
    "log10_source_refresh_length_um": (-1.0, 3.0),
    "c_blunt": (0.0, 3.0),
}


def parse_floats(text: str) -> list[float]:
    return [float(x) for x in str(text).replace(",", " ").split() if x]


def finite_value(row: pd.Series, name: str, default: float) -> float:
    try:
        value = float(row.get(name, default))
    except (TypeError, ValueError):
        return float(default)
    return value if np.isfinite(value) else float(default)


def shape_from_row(row: pd.Series) -> dict[str, float]:
    defaults = {
        "cleave_sT_GPa_per_K": 0.0,
        "cleave_exp_a": 0.2,
        "cleave_exp_n": 1.0,
        "cleave_floor_frac": 0.02,
        "emit_sT_GPa_per_K": 0.0,
        "emit_exp_a": 0.2,
        "emit_exp_n": 1.0,
        "emit_floor_frac": 0.02,
    }
    out = {}
    for name, default in defaults.items():
        out[name] = finite_value(row, f"shape_{name}", finite_value(row, name, default))
    return out


def fixed_parameters(row: pd.Series, scale: float) -> dict[str, float]:
    hp = finite_value(row, "peierls_H0_eV", 1.0)
    ht = finite_value(row, "taylor_H0_eV", hp + finite_value(row, "delta_H_PT_eV", 1.0))
    if ht < hp:
        ht = hp
    return {
        "cleave_G00_eV": finite_value(row, "cleave_G00_eV", 2.0),
        "cleave_gT_eV_per_K": finite_value(row, "cleave_gT_eV_per_K", 0.0),
        "cleave_sigc0_GPa": finite_value(row, "cleave_sigc0_GPa", 4.0),
        "emit_G00_eV": finite_value(row, "emit_G00_eV", 1.5),
        "emit_gT_eV_per_K": finite_value(row, "emit_gT_eV_per_K", 0.0),
        "emit_sigc0_GPa": finite_value(row, "emit_sigc0_GPa", 2.5),
        "peierls_H0_eV": max(float(scale) * hp, 1.0e-6),
        "taylor_H0_eV": max(float(scale) * ht, float(scale) * hp),
    }


def local_to_parameters(x: np.ndarray) -> dict[str, float]:
    raw = {name: float(value) for name, value in zip(LOCAL_NAMES, x)}
    raw.update(
        {
            "peierls_nu0_s": 10.0 ** raw["log10_peierls_nu0_s"],
            "taylor_nu0_s": 10.0 ** raw["log10_taylor_nu0_s"],
            "taylor_corr_rho_c_m2": 10.0 ** raw["log10_taylor_corr_rho_c_m2"],
            "taylor_corr_scale": 10.0 ** raw["log10_taylor_corr_scale"],
            "mobile_fraction": 10.0 ** raw["log10_mobile_fraction"],
            "source_sites_per_system": 10.0 ** raw["log10_source_sites_per_system"],
            "recovery_rate_s": 10.0 ** raw["log10_recovery_rate_s"],
            "source_refresh_length_um": 10.0 ** raw["log10_source_refresh_length_um"],
        }
    )
    return raw


def initial_local_vector(row: pd.Series, scale: float, reference_T_K: float) -> np.ndarray:
    hp = finite_value(row, "peierls_H0_eV", 1.0)
    ht = finite_value(row, "taylor_H0_eV", hp + finite_value(row, "delta_H_PT_eV", 1.0))
    logp = math.log10(max(finite_value(row, "peierls_nu0_s", 1.0e12), 1.0))
    logt = math.log10(max(finite_value(row, "taylor_nu0_s", 1.0e11), 1.0))
    # Warm-start transformation that approximately preserves each zero-stress
    # rate at the reference temperature after the common enthalpy scaling.
    denom = KB * float(reference_T_K) * math.log(10.0)
    logp += (float(scale) - 1.0) * hp * EV_TO_J / denom
    logt += (float(scale) - 1.0) * ht * EV_TO_J / denom
    values = {
        "peierls_activation_entropy_kB": finite_value(
            row, "peierls_activation_entropy_kB", -20.0
        ),
        "taylor_activation_entropy_kB": finite_value(
            row, "taylor_activation_entropy_kB", -20.0
        ),
        "log10_peierls_nu0_s": logp,
        "log10_taylor_nu0_s": logt,
        "log10_taylor_corr_rho_c_m2": finite_value(
            row, "log10_taylor_corr_rho_c_m2", 14.0
        ),
        "log10_taylor_corr_scale": finite_value(
            row, "log10_taylor_corr_scale", 0.0
        ),
        "log10_mobile_fraction": finite_value(row, "log10_mobile_fraction", -2.0),
        "log10_source_sites_per_system": finite_value(
            row, "log10_source_sites_per_system", math.log10(200.0)
        ),
        "log10_recovery_rate_s": finite_value(row, "log10_recovery_rate_s", -5.0),
        "log10_source_refresh_length_um": finite_value(
            row, "log10_source_refresh_length_um", math.log10(0.25)
        ),
        "c_blunt": finite_value(row, "c_blunt", 1.0),
    }
    x = np.array([values[name] for name in LOCAL_NAMES], dtype=float)
    bounds = np.asarray([LOCAL_BOUNDS[name] for name in LOCAL_NAMES], dtype=float)
    return np.clip(x, bounds[:, 0], bounds[:, 1])


def build_model(
    fixed: dict[str, float],
    local: dict[str, float],
    shape: dict[str, float],
    Tref_K: float,
) -> EmissionDerivedPeierlsTaylorModel:
    emit0 = max(fixed["emit_G00_eV"], 1.0e-12)
    parent = ExpFloorSurface(
        G00_eV=fixed["emit_G00_eV"],
        gT_eV_per_K=fixed["emit_gT_eV_per_K"],
        sigc0_Pa=fixed["emit_sigc0_GPa"] * 1.0e9,
        sT_Pa_per_K=shape["emit_sT_GPa_per_K"] * 1.0e9,
        Tref_K=Tref_K,
        a=shape["emit_exp_a"],
        n=shape["emit_exp_n"],
        floor_fraction=shape["emit_floor_frac"],
        floor_min_eV=1.0e-4,
        floor_max_fraction=0.95,
    )
    return EmissionDerivedPeierlsTaylorModel(
        EmissionDerivedPeierlsTaylorConfig(
            parent=parent,
            peierls=IndependentEntropyMechanismScale(
                fixed["peierls_H0_eV"] / emit0,
                local["peierls_activation_entropy_kB"],
                1.0,
                local["peierls_nu0_s"],
            ),
            taylor=IndependentEntropyMechanismScale(
                fixed["taylor_H0_eV"] / emit0,
                local["taylor_activation_entropy_kB"],
                1.0,
                local["taylor_nu0_s"],
            ),
            correlated_taylor=CorrelatedTaylorConfig(
                rho_c_m2=local["taylor_corr_rho_c_m2"],
                renewal_time_s=1.0,
                m_exponent=1.0,
                m_scale=local["taylor_corr_scale"],
                m_cap=float("inf"),
            ),
            mobile_fraction_low_density=local["mobile_fraction"],
            mobile_saturation_density_m2=float("inf"),
            mobile_density_floor_m2=0.0,
            jump_fraction_of_forest_spacing=1.0,
            jump_length_min_m=0.0,
            taylor_phi_max=float("inf"),
            rate_cap_s=float("inf"),
        )
    )


def interval_loss(value: float, low: float, high: float, scale: float) -> float:
    if not np.isfinite(value):
        return 100.0
    if value < low:
        return core.huber((value - low) / max(scale, 1.0e-12))
    if value > high:
        return core.huber((value - high) / max(scale, 1.0e-12))
    return 0.0


class ContinuationObjective:
    def __init__(
        self,
        target_class: str,
        fixed: dict[str, float],
        shape: dict[str, float],
        temperatures: np.ndarray,
        targets: pd.DataFrame,
        *,
        Kdot: float,
        dK: float,
        Kmax: float,
        r_pz_m: float,
        b_m: float = 2.74e-10,
        G_shear_Pa: float = 160.0e9,
        poisson: float = 0.28,
        rho0_m2: float = 5.0e12,
        Tref_K: float = 481.33,
    ):
        self.target_class = target_class
        self.fixed = fixed
        self.shape = shape
        self.temperatures = np.asarray(temperatures, dtype=float)
        self.targets = targets
        self.Kdot = float(Kdot)
        self.dK = float(dK)
        self.Kmax = float(Kmax)
        self.r_pz_m = float(r_pz_m)
        self.b_m = float(b_m)
        self.G_shear_Pa = float(G_shear_Pa)
        self.poisson = float(poisson)
        self.rho0_m2 = float(rho0_m2)
        self.Tref_K = float(Tref_K)
        self.bounds = np.asarray([LOCAL_BOUNDS[name] for name in LOCAL_NAMES], dtype=float)

    def evaluate(self, x: np.ndarray, details: bool = False) -> dict[str, Any]:
        x = np.asarray(x, dtype=float)
        if not np.all(np.isfinite(x)):
            return {"objective": 1.0e9, "nonfinite_parameter_vector": True}
        outside = np.maximum(self.bounds[:, 0] - x, 0.0) + np.maximum(
            x - self.bounds[:, 1], 0.0
        )
        if np.any(outside > 0.0):
            return {"objective": 1.0e8 + 1.0e6 * float(np.sum(outside**2))}
        local = local_to_parameters(x)
        model = build_model(self.fixed, local, self.shape, self.Tref_K)
        stress_grid = np.linspace(0.0, 30.0e9, 31)
        order_margin = core.barrier_order_margin_eV(
            model, self.temperatures, stress_grid
        )
        raw = [
            model.raw_zero_stress_barrier_eV(mech, T)
            for mech in ("peierls", "taylor")
            for T in self.temperatures
        ]
        min_raw = float(np.min(raw))
        if order_margin < -1.0e-9 or min_raw <= 0.0:
            penalty = 1.0e5 + 1.0e5 * max(-order_margin, 0.0) + 1.0e5 * max(-min_raw, 0.0)
            return {
                "objective": penalty,
                "barrier_order_margin_eV": order_margin,
                "min_raw_barrier_eV": min_raw,
            }

        fp = core.first_passage_proxy(
            self.fixed,
            self.shape,
            self.temperatures,
            Kdot=self.Kdot,
            dK=self.dK,
            Kmax=self.Kmax,
            r_pz_m=self.r_pz_m,
            nu0_c=1.0e12,
            nu0_e=1.0e11,
            cleavage_hits=3.0,
            cleavage_correlation_time_s=1.0e-6,
            Tref_K=self.Tref_K,
        )
        K0 = fp["K_intrinsic"]
        Hc = fp["H_emit_at_Kc"]
        if not np.all(np.isfinite(K0)) or not np.all(np.isfinite(Hc)):
            return {"objective": 5.0e4, "unresolved_first_passage": True}

        n_sites = 2.0 * local["source_sites_per_system"]
        K_per_line = (
            self.G_shear_Pa
            * self.b_m
            / max(1.0 - self.poisson, 1.0e-6)
            / math.sqrt(2.0 * math.pi * max(0.5 * self.r_pz_m, self.b_m))
            / 1.0e6
        )
        rows = []
        for j, T in enumerate(self.temperatures):
            sigma_tip = K0[j] * 1.0e6 / math.sqrt(2.0 * math.pi * self.r_pz_m)
            G_emit = float(model.barrier_eV("emission", sigma_tip, T))
            emit_rate = float(model._arrhenius_rate(G_emit, T, 1.0e11))
            N = 0.0
            rates: dict[str, Any] = {}
            retention = 0.0
            escape_rate = 0.0
            for _ in range(32):
                rho = self.rho0_m2 + N / max(math.pi * self.r_pz_m**2, 1.0e-30)
                rates = model.rates(sigma_tip, rho, T, self.b_m)
                escape_rate = float(np.asarray(rates["series_rate_s"]))
                denom = emit_rate + escape_rate + local["recovery_rate_s"]
                retention = emit_rate / denom if denom > 0.0 else 0.0
                Nnew = n_sites * retention
                if abs(Nnew - N) <= 1.0e-9 * max(Nnew, 1.0):
                    N = Nnew
                    break
                N = 0.5 * (N + Nnew)
            initial_emitted = n_sites * (1.0 - math.exp(-min(float(Hc[j]), 700.0)))
            N_init = min(initial_emitted * retention, N)
            N_dev = max(N, N_init)

            def resistance(count: float) -> tuple[float, float]:
                shield = K_per_line * max(count, 0.0)
                r_eff = self.r_pz_m + local["c_blunt"] * self.b_m * max(count, 0.0)
                blunt = K0[j] * (math.sqrt(r_eff / self.r_pz_m) - 1.0)
                return shield, blunt

            shield_i, blunt_i = resistance(N_init)
            shield_d, blunt_d = resistance(N_dev)
            Kinit = K0[j] + shield_i + blunt_i
            Kplateau = K0[j] + shield_d + blunt_d
            dKR = max(Kplateau - Kinit, 0.0)
            L = max(local["source_refresh_length_um"], 1.0e-9)
            early = dKR * (math.exp(-20.0 / L) - math.exp(-120.0 / L))
            late = dKR * (math.exp(-300.0 / L) - math.exp(-400.0 / L))
            p_rate = float(np.asarray(rates["peierls_rate_s"]))
            t_rate = float(np.asarray(rates["taylor_completion_rate_s"]))
            jump = float(np.asarray(rates["jump_length_m"]))
            load_time = max(K0[j] / max(self.Kdot, 1.0e-30), 1.0e-12)
            pi_p = jump * p_rate * load_time / max(self.r_pz_m, 1.0e-30)
            rows.append(
                {
                    "T_K": float(T),
                    "K_intrinsic": float(K0[j]),
                    "K_init_proxy": Kinit,
                    "K_plateau_proxy": Kplateau,
                    "delta_KR_proxy": dKR,
                    "early_rise_per_100um_proxy": early,
                    "plateau_rise_per_100um_proxy": late,
                    "N_initial_proxy": N_init,
                    "N_developed_proxy": N_dev,
                    "retained_fraction_proxy": retention,
                    "emission_rate_s": emit_rate,
                    "escape_rate_s": escape_rate,
                    "peierls_rate_s": p_rate,
                    "taylor_completion_rate_s": t_rate,
                    "jump_length_m": jump,
                    "peierls_traverse_number": pi_p,
                    "barrier_bottleneck_log10_P_over_T": math.log10(
                        max(p_rate, 1.0e-300) / max(t_rate, 1.0e-300)
                    ),
                }
            )
        pred = pd.DataFrame(rows)
        merged = self.targets.merge(pred, on="T_K", how="left")

        components: dict[str, float] = {}
        if self.target_class in {"ceramic", "weakT"}:
            components = {
                "K_init_loss": 0.0,
                "K_plateau_loss": 0.0,
                "early_rise_loss": 0.0,
                "plateau_rise_loss": 0.0,
                "delta_window_loss": 0.0,
            }
            for _, row in merged.iterrows():
                w = float(row.weight)
                components["K_init_loss"] += w * core.huber(
                    (row.K_init_proxy - row.K_init_target) / max(row.K_init_scale, 1.0e-9)
                )
                components["K_plateau_loss"] += w * core.huber(
                    (row.K_plateau_proxy - row.K_plateau_target)
                    / max(row.K_plateau_scale, 1.0e-9)
                )
                components["early_rise_loss"] += w * core.huber(
                    (row.early_rise_per_100um_proxy - row.early_rise_per_100um_target)
                    / max(row.early_rise_scale, 1.0e-9)
                )
                components["plateau_rise_loss"] += w * core.huber(
                    (row.plateau_rise_per_100um_proxy - row.plateau_rise_per_100um_target)
                    / max(row.plateau_rise_scale, 1.0e-9)
                )
                if row.delta_KR_proxy < row.delta_KR_min:
                    components["delta_window_loss"] += w * core.huber(
                        row.delta_KR_min - row.delta_KR_proxy
                    )
                elif row.delta_KR_proxy > row.delta_KR_max:
                    components["delta_window_loss"] += w * core.huber(
                        row.delta_KR_proxy - row.delta_KR_max
                    )
        else:
            low = pred[pred.T_K <= 700.0]
            high = pred[pred.T_K >= 900.0]
            Kp_low = float(low.K_plateau_proxy.median())
            Kp_high = float(high.K_plateau_proxy.median())
            Ki_low = float(low.K_init_proxy.median())
            Ki_high = float(high.K_init_proxy.median())
            high_dKR = float(high.delta_KR_proxy.median())
            low_dKR = float(low.delta_KR_proxy.median())
            monotonic_drop = float(
                np.sum(np.maximum(-np.diff(pred.sort_values("T_K").K_plateau_proxy.to_numpy()) - 2.0, 0.0))
            )
            components = {
                "DBTT_plateau_rise_loss": core.huber(max(20.0 - (Kp_high - Kp_low), 0.0) / 4.0),
                "DBTT_initiation_rise_loss": core.huber(max(10.0 - (Ki_high - Ki_low), 0.0) / 3.0),
                "DBTT_high_Rcurve_loss": core.huber(max(5.0 - high_dKR, 0.0) / 2.0),
                "DBTT_low_Rcurve_loss": core.huber(max(low_dKR - 3.0, 0.0) / 1.5),
                "DBTT_low_shelf_loss": interval_loss(Kp_low, 8.0, 30.0, 4.0),
                "DBTT_high_branch_loss": interval_loss(Kp_high, 30.0, 65.0, 5.0),
                "DBTT_monotonicity_loss": core.huber(monotonic_drop / 3.0),
            }

        if self.target_class == "weakT":
            min_pi = float(np.min(pred.peierls_traverse_number))
            components["fast_Peierls_loss"] = 20.0 * core.huber(
                max(-math.log10(max(min_pi, 1.0e-300)), 0.0)
            )
        objective = float(sum(components.values()))
        if not np.isfinite(objective):
            return {"objective": 1.0e9, "nonfinite_objective_replaced": True}
        result: dict[str, Any] = {
            "objective": objective,
            **components,
            "barrier_order_margin_eV": order_margin,
            "min_raw_barrier_eV": min_raw,
            "min_peierls_traverse_number": float(pred.peierls_traverse_number.min()),
            "plateau_temperature_rise": float(
                pred.sort_values("T_K").K_plateau_proxy.iloc[-1]
                - pred.sort_values("T_K").K_plateau_proxy.iloc[0]
            ),
            "local_parameters": local,
        }
        if details:
            result["temperature_detail"] = merged.to_dict(orient="records")
        return result

    def __call__(self, x: np.ndarray) -> float:
        return float(self.evaluate(x, details=False)["objective"])


def acceptance(target_class: str, detail: pd.DataFrame, summary: dict[str, Any]) -> tuple[bool, str]:
    if not np.isfinite(float(summary.get("objective", np.nan))):
        return False, "nonfinite_objective"
    if float(summary.get("barrier_order_margin_eV", -1.0)) < -1.0e-8:
        return False, "barrier_order_violation"
    if target_class == "ceramic":
        ok = float(detail.delta_KR_proxy.max()) <= 3.0
        return ok, "ceramic_small_Rcurve" if ok else "ceramic_Rcurve_too_large"
    if target_class == "weakT":
        pi_ok = float(detail.peierls_traverse_number.min()) >= 1.0
        flat = float(detail.K_plateau_proxy.max() - detail.K_plateau_proxy.min()) <= 6.0
        rise = float(detail.delta_KR_proxy.median())
        ok = pi_ok and flat and 2.0 <= rise <= 9.0
        return ok, "weakT_fast_Peierls" if ok else "weakT_acceptance_failed"
    low = detail[detail.T_K <= 700.0]
    high = detail[detail.T_K >= 900.0]
    plateau_rise = float(high.K_plateau_proxy.median() - low.K_plateau_proxy.median())
    high_rcurve = float(high.delta_KR_proxy.median())
    ok = plateau_rise >= 15.0 and high_rcurve >= 3.0
    return ok, "DBTT_trend_present" if ok else "DBTT_trend_insufficient"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--input-root",
        type=Path,
        default=Path("runs/mpz_v9_8_1_joint_response_optimization_v1"),
    )
    ap.add_argument("--classes", default="ceramic weakT DBTT")
    ap.add_argument("--scales", default="1.0 0.8 0.6 0.4 0.3")
    ap.add_argument("--temperatures", default="300 700 900 1200")
    ap.add_argument("--candidates-per-class", type=int, default=3)
    ap.add_argument("--local-maxiter", type=int, default=400)
    ap.add_argument("--Kdot", type=float, default=0.005)
    ap.add_argument("--dK", type=float, default=0.5)
    ap.add_argument("--Kmax", type=float, default=65.0)
    ap.add_argument("--r-pz-m", type=float, default=1.0e-6)
    ap.add_argument("--reference-temperature-K", type=float, default=700.0)
    ap.add_argument(
        "--targets", type=Path, default=Path("mpz_three_class_design_targets.csv")
    )
    ap.add_argument(
        "--out", type=Path, default=Path("runs/mpz_v9_9_barrier_continuation_v1")
    )
    a = ap.parse_args()

    classes = str(a.classes).replace(",", " ").split()
    scales = parse_floats(a.scales)
    temperatures = np.asarray(parse_floats(a.temperatures), dtype=float)
    out = a.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    summary_rows: list[dict[str, Any]] = []
    detail_rows: list[dict[str, Any]] = []

    for target_class in classes:
        path = a.input_root / target_class / "joint_response_shortlist.csv"
        if not path.exists():
            raise SystemExit(f"v9.8.1 shortlist not found: {path}")
        candidates = pd.read_csv(path).sort_values("objective").head(a.candidates_per_class)
        targets = core.load_targets(a.targets, target_class, temperatures.tolist())
        for rank, (_, row) in enumerate(candidates.iterrows()):
            shape = shape_from_row(row)
            for scale in scales:
                fixed = fixed_parameters(row, scale)
                objective = ContinuationObjective(
                    target_class,
                    fixed,
                    shape,
                    temperatures,
                    targets,
                    Kdot=a.Kdot,
                    dK=a.dK,
                    Kmax=a.Kmax,
                    r_pz_m=a.r_pz_m,
                )
                x0 = initial_local_vector(row, scale, a.reference_temperature_K)
                initial = objective.evaluate(x0, details=False)
                result = minimize(
                    objective,
                    x0,
                    method="Powell",
                    bounds=[LOCAL_BOUNDS[name] for name in LOCAL_NAMES],
                    options={"maxiter": a.local_maxiter, "xtol": 1.0e-4, "ftol": 1.0e-5},
                )
                best_x = result.x if np.isfinite(result.fun) and result.fun <= initial["objective"] else x0
                best = objective.evaluate(best_x, details=True)
                local = best.pop("local_parameters")
                detail = pd.DataFrame(best.pop("temperature_detail"))
                accepted, reason = acceptance(target_class, detail, best)
                candidate_id = f"{target_class}_rank{rank:02d}_scale{scale:g}"
                rec = {
                    "continuation_candidate_id": candidate_id,
                    "target_class": target_class,
                    "source_rank": rank,
                    "source_seed_candidate_id": str(row.get("seed_candidate_id", "")),
                    "source_basin_index": int(finite_value(row, "basin_index", rank)),
                    "barrier_scale": float(scale),
                    "source_objective": finite_value(row, "objective", np.nan),
                    "initial_objective": float(initial["objective"]),
                    "objective": float(best["objective"]),
                    "local_success": bool(result.success),
                    "local_message": str(result.message),
                    "accepted_for_spatial_promotion": bool(accepted),
                    "acceptance_reason": reason,
                    **fixed,
                    **{name: float(best_x[i]) for i, name in enumerate(LOCAL_NAMES)},
                    **{
                        k: float(v)
                        for k, v in local.items()
                        if k not in LOCAL_NAMES
                    },
                    **{k: v for k, v in best.items() if k != "objective"},
                    **{f"shape_{k}": float(v) for k, v in shape.items()},
                    "status": "ANALYTICAL_CONTINUATION_REQUIRES_SPATIAL_MPZ",
                }
                summary_rows.append(rec)
                detail.insert(0, "continuation_candidate_id", candidate_id)
                detail.insert(1, "target_class", target_class)
                detail.insert(2, "barrier_scale", float(scale))
                detail_rows.extend(detail.to_dict(orient="records"))
                print(
                    f"class={target_class} rank={rank} scale={scale:g} "
                    f"objective={best['objective']:.6g} accepted={accepted} "
                    f"minPiP={rec['min_peierls_traverse_number']:.3g}",
                    flush=True,
                )

    summary = pd.DataFrame(summary_rows).sort_values(
        ["target_class", "accepted_for_spatial_promotion", "objective"],
        ascending=[True, False, True],
    )
    detail = pd.DataFrame(detail_rows)
    accepted = summary[summary.accepted_for_spatial_promotion.astype(bool)].copy()
    promotion = (
        accepted.sort_values(["target_class", "objective"])
        .groupby("target_class", as_index=False)
        .head(5)
        .reset_index(drop=True)
    )
    summary.to_csv(out / "barrier_continuation_all.csv", index=False)
    accepted.to_csv(out / "barrier_continuation_accepted.csv", index=False)
    promotion.to_csv(out / "spatial_promotion_manifest.csv", index=False)
    detail.to_csv(out / "barrier_continuation_temperature_detail.csv", index=False)
    report = {
        "n_evaluated": int(len(summary)),
        "n_accepted": int(len(accepted)),
        "accepted_by_class": accepted.target_class.value_counts().to_dict(),
        "scales": scales,
        "classes": classes,
        "output": str(out),
        "status": "V9_9_ANALYTICAL_CONTINUATION_COMPLETE",
    }
    (out / "barrier_continuation_summary.json").write_text(json.dumps(report, indent=2))
    config = vars(a).copy()
    for key, value in list(config.items()):
        if isinstance(value, Path):
            config[key] = str(value)
    config.update(
        {
            "local_parameter_names": LOCAL_NAMES,
            "local_bounds": LOCAL_BOUNDS,
            "common_barrier_scale_preserves_Ht_over_Hp": True,
            "mechanism_rules": {"ceramic": [], "weakT": ["Peierls traverse number >= 1"], "DBTT": []},
            "production_status": "NOT_ACTIVATED",
        }
    )
    (out / "barrier_continuation_config.json").write_text(json.dumps(config, indent=2))
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
