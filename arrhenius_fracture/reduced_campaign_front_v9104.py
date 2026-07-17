"""PF-equivalent reduced moving-tip front model for narrow-DBTT searches.

This is a deliberately reduced, prescribed-K counterpart of the finalized
sharp-front campaign model. It preserves the constitutive channel separation
and moving-tip state evolution needed for calibration while omitting the 2-D
FEM solve, branching, and mesh topology.

The three stress channels are intentionally distinct:

* opening:   sigma_open = K / sqrt(2*pi*r_eff)
* cleavage:  sigma_c = (K - K_shield) / sqrt(2*pi*r_eff)
* emission:  sigma_e,s = max(w_s*sigma_open - sigma_back,s, 0)

Elastic shielding therefore acts on cleavage only. The local Taylor back
stress acts on emission only. Source recovery is driven only by crack advance.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
import math
from typing import Any, Iterable

import numpy as np
from scipy.special import gammainc

from .config import EV_TO_J, KB
from .emission_derived_plasticity import (
    CorrelatedTaylorConfig,
    EmissionDerivedPeierlsTaylorConfig,
    ExpFloorSurface,
)
from .emission_derived_plasticity_v9102 import (
    EmissionDerivedPeierlsTaylorModel,
    IndependentShapeEntropyMechanismScale,
)


@dataclass(frozen=True)
class ReducedFrontSettings:
    Kdot_MPa_sqrt_m_s: float = 0.005
    Kmax_MPa_sqrt_m: float = 80.0
    max_dK_substep_MPa_sqrt_m: float = 0.05
    target_extension_um: float = 5.0
    checkpoint_da_um: float = 5.0
    r0_m: float = 1.0e-6
    L_pz_m: float = 50.0e-6
    b_m: float = 2.74e-10
    G_Pa: float = 160.0e9
    poisson: float = 0.28
    rho0_m2: float = 5.0e12
    nu0_c_s: float = 1.0e12
    nu0_e_s: float = 1.0e11
    cleavage_hits: float = 3.0
    cleavage_tau_s: float = 1.0e-6
    Tref_K: float = 481.33
    n_systems: int = 2
    resolved_stress_fraction: float = 1.0
    backstress_scale: float = 1.0
    max_K_shield_MPa_sqrt_m: float = 1.0
    mobile_shield_fraction: float = 0.0
    max_action_substep: float = 0.01
    max_emit_fraction_substep: float = 0.02
    max_exchange_fraction_substep: float = 0.10
    min_dt_s: float = 1.0e-12
    max_internal_steps: int = 500_000
    event_tolerance: float = 1.0e-10

    def with_extension(self, extension_um: float) -> "ReducedFrontSettings":
        return replace(self, target_extension_um=float(extension_um))


@dataclass(frozen=True)
class TransitionRequirements:
    low_shelf_min: float = 8.0
    low_shelf_max: float = 25.0
    high_shelf_max: float = 70.0
    min_ratio: float = 2.0
    robust_ratio: float = 1.8
    max_low_span_fraction: float = 0.15
    max_high_span_fraction: float = 0.20
    min_jump_concentration: float = 0.75
    max_plasticity_off_ratio: float = 1.25
    min_points_per_shelf: int = 2


def exp_floor_barrier_eV(
    sigma_Pa: float | np.ndarray,
    T_K: float | np.ndarray,
    G00_eV: float,
    gT_eV_per_K: float,
    sigc0_GPa: float,
    sT_GPa_per_K: float,
    exp_a: float,
    exp_n: float,
    floor_fraction: float,
    Tref_K: float,
) -> np.ndarray:
    sigma = np.maximum(np.asarray(sigma_Pa, dtype=float), 0.0)
    T = np.asarray(T_K, dtype=float)
    G0 = np.maximum(float(G00_eV) + float(gT_eV_per_K) * (T - Tref_K), 1.0e-12)
    sigc = np.maximum(
        (float(sigc0_GPa) + float(sT_GPa_per_K) * (T - Tref_K)) * 1.0e9,
        1.0,
    )
    floor = np.minimum(
        0.95 * G0,
        np.maximum(1.0e-4, float(floor_fraction) * G0),
    )
    return floor + (G0 - floor) * np.exp(
        -max(float(exp_a), 0.0)
        * np.power(sigma / sigc, max(float(exp_n), 1.0e-9))
    )


def cleavage_effective_rate_s(
    G_eV: float,
    T_K: float,
    nu0_s: float,
    m_hits: float,
    correlation_time_s: float,
) -> float:
    raw = float(nu0_s) * math.exp(
        float(np.clip(-float(G_eV) * EV_TO_J / (KB * float(T_K)), -700.0, 0.0))
    )
    if m_hits <= 1.0 + 1.0e-12:
        return raw
    return float(gammainc(float(m_hits), min(raw * float(correlation_time_s), 1.0e12))) / float(
        correlation_time_s
    )


def exact_depletion(available: np.ndarray, rate_s: np.ndarray, dt_s: float) -> np.ndarray:
    """Finite-source depletion, exactly invariant under timestep partition."""
    avail = np.maximum(np.asarray(available, dtype=float), 0.0)
    rate = np.maximum(np.asarray(rate_s, dtype=float), 0.0)
    fraction = 1.0 - np.exp(-np.minimum(rate * max(float(dt_s), 0.0), 700.0))
    return np.minimum(avail * fraction, avail)


def exact_refresh(
    available: np.ndarray,
    capacity: np.ndarray,
    da_m: float,
    refresh_length_m: float,
) -> np.ndarray:
    """Crack-advance-only source recovery, invariant under advance partition."""
    avail = np.maximum(np.asarray(available, dtype=float), 0.0)
    cap = np.maximum(np.asarray(capacity, dtype=float), 0.0)
    if da_m <= 0.0:
        return avail.copy()
    fraction = 1.0 - math.exp(-max(float(da_m), 0.0) / max(float(refresh_length_m), 1.0e-30))
    return np.minimum(avail + (cap - avail) * fraction, cap)


def _build_pt_model(p: dict[str, float], Tref_K: float) -> EmissionDerivedPeierlsTaylorModel:
    emit0 = max(float(p["emit_G00_eV"]), 1.0e-12)
    parent = ExpFloorSurface(
        G00_eV=float(p["emit_G00_eV"]),
        gT_eV_per_K=float(p["emit_gT_eV_per_K"]),
        sigc0_Pa=float(p["emit_sigc0_GPa"]) * 1.0e9,
        sT_Pa_per_K=float(p["emit_sT_GPa_per_K"]) * 1.0e9,
        Tref_K=float(Tref_K),
        a=float(p["emit_exp_a"]),
        n=float(p["emit_exp_n"]),
        floor_fraction=float(p["emit_floor_frac"]),
        floor_min_eV=1.0e-4,
        floor_max_fraction=0.95,
    )
    return EmissionDerivedPeierlsTaylorModel(
        EmissionDerivedPeierlsTaylorConfig(
            parent=parent,
            peierls=IndependentShapeEntropyMechanismScale(
                energy_ratio=float(p["peierls_H0_eV"]) / emit0,
                activation_entropy_kB=float(p["peierls_activation_entropy_kB"]),
                exp_a=float(p["peierls_exp_a"]),
                exp_n=float(p["peierls_exp_n"]),
                stress_ratio=1.0,
                rate_prefactor_s=float(p.get("peierls_nu0_s", 1.0e12)),
            ),
            taylor=IndependentShapeEntropyMechanismScale(
                energy_ratio=float(p["taylor_H0_eV"]) / emit0,
                activation_entropy_kB=float(p["taylor_activation_entropy_kB"]),
                exp_a=float(p["taylor_exp_a"]),
                exp_n=float(p["taylor_exp_n"]),
                stress_ratio=1.0,
                rate_prefactor_s=float(p.get("taylor_nu0_s", 1.0e11)),
            ),
            correlated_taylor=CorrelatedTaylorConfig(
                rho_c_m2=float(p["taylor_corr_rho_c_m2"]),
                renewal_time_s=1.0,
                m_exponent=1.0,
                m_scale=float(p["taylor_corr_scale"]),
                m_cap=float("inf"),
            ),
            mobile_saturation_density_m2=float("inf"),
            mobile_density_floor_m2=0.0,
            jump_length_min_m=0.0,
            taylor_phi_max=float("inf"),
            rate_cap_s=float("inf"),
        )
    )


class ReducedCampaignFront:
    """Prescribed-K moving-tip state with PF-equivalent channel separation."""

    def __init__(
        self,
        parameters: dict[str, float],
        settings: ReducedFrontSettings,
        *,
        mode: str = "full",
    ) -> None:
        allowed = {
            "full",
            "plasticity_off",
            "backstress_off",
            "shielding_off",
            "blunting_off",
        }
        if mode not in allowed:
            raise ValueError(f"unknown reduced-front mode: {mode}")
        self.p = dict(parameters)
        self.s = settings
        self.mode = mode
        self.pt_model = _build_pt_model(self.p, settings.Tref_K)
        n = int(settings.n_systems)
        capacity_per_system = max(float(self.p["source_sites_per_system"]), 0.0)
        self.capacity = np.full(n, capacity_per_system, dtype=float)
        self.available = self.capacity.copy()
        self.mobile = np.zeros(n, dtype=float)
        self.retained = np.zeros(n, dtype=float)
        self.slip = np.zeros(n, dtype=float)
        self.system_weights = np.ones(n, dtype=float)
        self.B = 0.0
        self.K = 0.0
        self.time_s = 0.0
        self.a_um = 0.0
        self.cumulative_emitted = 0.0
        self.cumulative_refreshed = 0.0
        self.cumulative_trapped = 0.0
        self.cumulative_released = 0.0
        self.cumulative_recovered = 0.0
        self.cumulative_escaped = 0.0
        self.events: list[dict[str, float]] = []

    @property
    def plasticity_active(self) -> bool:
        return self.mode != "plasticity_off"

    def r_eff_m(self) -> float:
        if self.mode == "blunting_off" or not self.plasticity_active:
            return float(self.s.r0_m)
        return float(
            self.s.r0_m
            + max(float(self.p.get("c_blunt", 0.0)), 0.0)
            * self.s.b_m
            * float(np.sum(np.maximum(self.slip, 0.0)))
        )

    def active_density_m2(self) -> np.ndarray:
        area = math.pi * max(self.s.r0_m, self.s.b_m) ** 2
        return self.s.rho0_m2 + (np.maximum(self.mobile, 0.0) + np.maximum(self.retained, 0.0)) / area

    def emission_backstress_Pa(self) -> np.ndarray:
        if self.mode == "backstress_off" or not self.plasticity_active:
            return np.zeros_like(self.mobile)
        rho = np.maximum(self.active_density_m2(), 0.0)
        tau = (
            max(float(self.s.backstress_scale), 0.0)
            * self.s.G_Pa
            * self.s.b_m
            * np.sqrt(rho)
        )
        return tau / max(float(self.s.resolved_stress_fraction), 1.0e-12)

    def K_shield_MPa_sqrt_m(self) -> float:
        if self.mode == "shielding_off" or not self.plasticity_active:
            return 0.0
        coefficient = (
            self.s.G_Pa
            * self.s.b_m
            / max(1.0 - self.s.poisson, 1.0e-12)
            / math.sqrt(2.0 * math.pi * max(0.5 * self.s.r0_m, self.s.b_m))
            / 1.0e6
        )
        count = float(
            np.sum(
                np.maximum(self.retained, 0.0)
                + max(self.s.mobile_shield_fraction, 0.0) * np.maximum(self.mobile, 0.0)
            )
        )
        raw = coefficient * count
        return float(np.clip(raw, 0.0, max(self.s.max_K_shield_MPa_sqrt_m, 0.0)))

    def stress_channels(self, K_MPa_sqrt_m: float) -> dict[str, Any]:
        r_eff = self.r_eff_m()
        denominator = math.sqrt(2.0 * math.pi * max(r_eff, 1.0e-30))
        sigma_open = max(float(K_MPa_sqrt_m), 0.0) * 1.0e6 / denominator
        K_shield = self.K_shield_MPa_sqrt_m()
        sigma_cleave = max(float(K_MPa_sqrt_m) - K_shield, 0.0) * 1.0e6 / denominator
        back = self.emission_backstress_Pa()
        sigma_emit = np.maximum(self.system_weights * sigma_open - back, 0.0)
        return {
            "r_eff_m": r_eff,
            "K_shield_MPa_sqrt_m": K_shield,
            "sigma_open_Pa": sigma_open,
            "sigma_cleave_Pa": sigma_cleave,
            "sigma_back_Pa": back,
            "sigma_emit_Pa": sigma_emit,
        }

    def _barrier(self, mechanism: str, sigma_Pa: float | np.ndarray, T_K: float) -> np.ndarray:
        prefix = "cleave" if mechanism == "cleave" else "emit"
        return exp_floor_barrier_eV(
            sigma_Pa,
            T_K,
            self.p[f"{prefix}_G00_eV"],
            self.p[f"{prefix}_gT_eV_per_K"],
            self.p[f"{prefix}_sigc0_GPa"],
            self.p[f"{prefix}_sT_GPa_per_K"],
            self.p[f"{prefix}_exp_a"],
            self.p[f"{prefix}_exp_n"],
            self.p[f"{prefix}_floor_frac"],
            self.s.Tref_K,
        )

    def instantaneous_rates(self, K_MPa_sqrt_m: float, T_K: float) -> dict[str, Any]:
        channels = self.stress_channels(K_MPa_sqrt_m)
        Gc = float(self._barrier("cleave", channels["sigma_cleave_Pa"], T_K))
        lambda_c = cleavage_effective_rate_s(
            Gc,
            T_K,
            self.s.nu0_c_s,
            self.s.cleavage_hits,
            self.s.cleavage_tau_s,
        )
        if self.plasticity_active:
            Ge = np.asarray(self._barrier("emit", channels["sigma_emit_Pa"], T_K), dtype=float)
            lambda_e = self.s.nu0_e_s * np.exp(
                np.clip(-Ge * EV_TO_J / (KB * float(T_K)), -700.0, 0.0)
            )
            rho = np.maximum(self.active_density_m2(), 1.0)
            pt = self.pt_model.rates(
                channels["sigma_open_Pa"], rho, T_K, self.s.b_m
            )
            p_rate = np.asarray(pt["peierls_rate_s"], dtype=float).reshape(-1)
            t_rate = np.asarray(pt["taylor_completion_rate_s"], dtype=float).reshape(-1)
            jump = np.asarray(pt["jump_length_m"], dtype=float).reshape(-1)
            if p_rate.size == 1:
                p_rate = np.full(self.s.n_systems, float(p_rate[0]))
                t_rate = np.full(self.s.n_systems, float(t_rate[0]))
                jump = np.full(self.s.n_systems, float(jump[0]))
            encounter = (
                max(float(self.p["encounter_efficiency"]), 0.0)
                * jump
                * p_rate
                * np.sqrt(rho)
            )
            velocity = jump * p_rate
        else:
            Ge = np.full(self.s.n_systems, np.nan)
            lambda_e = np.zeros(self.s.n_systems)
            p_rate = np.zeros(self.s.n_systems)
            t_rate = np.zeros(self.s.n_systems)
            encounter = np.zeros(self.s.n_systems)
            velocity = np.zeros(self.s.n_systems)
        return {
            **channels,
            "G_cleave_eV": Gc,
            "G_emit_eV": Ge,
            "lambda_c_s": lambda_c,
            "lambda_e_s": lambda_e,
            "peierls_rate_s": p_rate,
            "taylor_rate_s": t_rate,
            "encounter_rate_s": encounter,
            "velocity_m_s": velocity,
        }

    def _choose_dt(self, rates: dict[str, Any]) -> float:
        h = self.s.max_dK_substep_MPa_sqrt_m / max(self.s.Kdot_MPa_sqrt_m_s, 1.0e-30)
        lam_c = max(float(rates["lambda_c_s"]), 0.0)
        if lam_c > 0.0:
            h = min(h, self.s.max_action_substep / lam_c)
            h = min(h, max(1.0 - self.B, self.s.event_tolerance) / lam_c)
        if self.plasticity_active:
            lam_e = float(np.max(np.asarray(rates["lambda_e_s"], dtype=float)))
            if lam_e > 0.0 and self.s.max_emit_fraction_substep < 1.0:
                h = min(
                    h,
                    -math.log(max(1.0 - self.s.max_emit_fraction_substep, 1.0e-12)) / lam_e,
                )
            kinetic = float(
                np.max(
                    np.asarray(rates["encounter_rate_s"], dtype=float)
                    + np.asarray(rates["taylor_rate_s"], dtype=float)
                    + np.asarray(rates["velocity_m_s"], dtype=float) / max(self.s.L_pz_m, 1.0e-30)
                )
            )
            if kinetic > 0.0:
                h = min(h, self.s.max_exchange_fraction_substep / kinetic)
        return max(float(h), self.s.min_dt_s)

    def _plastic_step(self, dt_s: float, K_MPa_sqrt_m: float, T_K: float) -> dict[str, float]:
        if not self.plasticity_active or dt_s <= 0.0:
            return {
                "emitted": 0.0,
                "trapped": 0.0,
                "released": 0.0,
                "recovered": 0.0,
                "escaped": 0.0,
            }
        rates = self.instantaneous_rates(K_MPa_sqrt_m, T_K)
        emitted = exact_depletion(self.available, rates["lambda_e_s"], dt_s)
        self.available -= emitted
        self.mobile += emitted
        self.slip += emitted

        total = np.maximum(self.mobile, 0.0) + np.maximum(self.retained, 0.0)
        ke = np.maximum(np.asarray(rates["encounter_rate_s"], dtype=float), 0.0)
        kt = np.maximum(np.asarray(rates["taylor_rate_s"], dtype=float), 0.0)
        exchange = ke + kt
        frac_eq = np.divide(ke, exchange, out=np.zeros_like(exchange), where=exchange > 0.0)
        retained_eq = frac_eq * total
        decay = np.exp(-np.minimum(exchange * dt_s, 700.0))
        new_retained = retained_eq + (self.retained - retained_eq) * decay
        new_retained = np.clip(new_retained, 0.0, total)
        trapped = float(np.sum(np.maximum(new_retained - self.retained, 0.0)))
        released = float(np.sum(np.maximum(self.retained - new_retained, 0.0)))
        self.retained = new_retained
        self.mobile = total - self.retained

        recovery_rate = max(float(self.p["retained_recovery_rate_s"]), 0.0)
        recovered_vec = self.retained * (1.0 - math.exp(-min(recovery_rate * dt_s, 700.0)))
        self.retained -= recovered_vec

        escape_rate = np.maximum(np.asarray(rates["velocity_m_s"], dtype=float), 0.0) / max(
            self.s.L_pz_m, 1.0e-30
        )
        mobile_before = self.mobile.copy()
        self.mobile *= np.exp(-np.minimum(escape_rate * dt_s, 700.0))
        escaped = float(np.sum(np.maximum(mobile_before - self.mobile, 0.0)))

        values = {
            "emitted": float(np.sum(emitted)),
            "trapped": trapped,
            "released": released,
            "recovered": float(np.sum(recovered_vec)),
            "escaped": escaped,
        }
        self.cumulative_emitted += values["emitted"]
        self.cumulative_trapped += values["trapped"]
        self.cumulative_released += values["released"]
        self.cumulative_recovered += values["recovered"]
        self.cumulative_escaped += values["escaped"]
        return values

    def _translate_and_refresh(self, da_um: float) -> float:
        da_m = max(float(da_um), 0.0) * 1.0e-6
        if da_m <= 0.0:
            return 0.0
        keep = math.exp(-da_m / max(self.s.L_pz_m, 1.0e-30))
        self.mobile *= keep
        self.retained *= keep
        self.slip *= keep
        before = self.available.copy()
        self.available = exact_refresh(
            self.available,
            self.capacity,
            da_m,
            max(float(self.p["source_refresh_length_um"]), 1.0e-12) * 1.0e-6,
        )
        refreshed = float(np.sum(self.available - before))
        self.cumulative_refreshed += refreshed
        return refreshed

    def _record_event(self, T_K: float) -> None:
        rates = self.instantaneous_rates(self.K, T_K)
        self.events.append(
            {
                "event_index": float(len(self.events) + 1),
                "a_um": float(self.a_um),
                "K_MPa_sqrt_m": float(self.K),
                "K_shield_MPa_sqrt_m": float(rates["K_shield_MPa_sqrt_m"]),
                "sigma_open_Pa": float(rates["sigma_open_Pa"]),
                "sigma_cleave_Pa": float(rates["sigma_cleave_Pa"]),
                "sigma_back_max_Pa": float(np.max(rates["sigma_back_Pa"])),
                "sigma_emit_max_Pa": float(np.max(rates["sigma_emit_Pa"])),
                "mobile_count": float(np.sum(self.mobile)),
                "retained_count": float(np.sum(self.retained)),
                "available_sites": float(np.sum(self.available)),
                "cumulative_emitted": float(self.cumulative_emitted),
                "cumulative_refreshed": float(self.cumulative_refreshed),
                "r_eff_m": float(rates["r_eff_m"]),
                "time_s": float(self.time_s),
            }
        )

    def run(self, T_K: float) -> dict[str, Any]:
        target_events = max(int(math.ceil(self.s.target_extension_um / self.s.checkpoint_da_um)), 1)
        internal_steps = 0
        while (
            self.K < self.s.Kmax_MPa_sqrt_m
            and len(self.events) < target_events
            and internal_steps < self.s.max_internal_steps
        ):
            internal_steps += 1
            rates0 = self.instantaneous_rates(self.K, T_K)
            h = self._choose_dt(rates0)
            h = min(
                h,
                (self.s.Kmax_MPa_sqrt_m - self.K)
                / max(self.s.Kdot_MPa_sqrt_m_s, 1.0e-30),
            )
            if h < self.s.min_dt_s:
                break
            K_mid = self.K + 0.5 * self.s.Kdot_MPa_sqrt_m_s * h

            snapshot = (
                self.available.copy(),
                self.mobile.copy(),
                self.retained.copy(),
                self.slip.copy(),
                self.cumulative_emitted,
                self.cumulative_trapped,
                self.cumulative_released,
                self.cumulative_recovered,
                self.cumulative_escaped,
            )
            self._plastic_step(0.5 * h, K_mid, T_K)
            mid = self.instantaneous_rates(K_mid, T_K)
            lam_c = max(float(mid["lambda_c_s"]), 0.0)
            allowed_action = min(self.s.max_action_substep, max(1.0 - self.B, 0.0))
            if lam_c > 0.0 and lam_c * h > allowed_action * (1.0 + 1.0e-10):
                (
                    self.available,
                    self.mobile,
                    self.retained,
                    self.slip,
                    self.cumulative_emitted,
                    self.cumulative_trapped,
                    self.cumulative_released,
                    self.cumulative_recovered,
                    self.cumulative_escaped,
                ) = snapshot
                h = max(allowed_action / lam_c, self.s.min_dt_s)
                K_mid = self.K + 0.5 * self.s.Kdot_MPa_sqrt_m_s * h
                self._plastic_step(0.5 * h, K_mid, T_K)
                mid = self.instantaneous_rates(K_mid, T_K)
                lam_c = max(float(mid["lambda_c_s"]), 0.0)

            dB = min(lam_c * h, self.s.max_action_substep, max(1.0 - self.B, 0.0))
            da_um = self.s.checkpoint_da_um * dB
            self.B += dB
            self.a_um += da_um
            self._translate_and_refresh(da_um)
            self._plastic_step(0.5 * h, K_mid, T_K)
            self.time_s += h
            self.K += self.s.Kdot_MPa_sqrt_m_s * h

            if self.B >= 1.0 - self.s.event_tolerance:
                self.B = 0.0
                self.a_um = len(self.events) * self.s.checkpoint_da_um + self.s.checkpoint_da_um
                self._record_event(T_K)

        return summarize_run(self, T_K, target_events, internal_steps)


def summarize_run(
    front: ReducedCampaignFront,
    T_K: float,
    target_events: int,
    internal_steps: int,
) -> dict[str, Any]:
    events = front.events
    if not events:
        return {
            "completed": False,
            "T_K": float(T_K),
            "mode": front.mode,
            "n_events": 0,
            "internal_steps": int(internal_steps),
            "K_init_proxy": float("nan"),
            "K_plateau_proxy": float("nan"),
            "delta_KR_proxy": float("nan"),
            "events": [],
        }
    K = np.asarray([e["K_MPa_sqrt_m"] for e in events], dtype=float)
    a = np.asarray([e["a_um"] for e in events], dtype=float)
    late = K[a >= 0.70 * front.s.target_extension_um]
    Kplateau = float(np.median(late)) if late.size else float(K[-1])
    return {
        "completed": len(events) >= target_events,
        "T_K": float(T_K),
        "mode": front.mode,
        "n_events": len(events),
        "internal_steps": int(internal_steps),
        "K_init_proxy": float(K[0]),
        "K_plateau_proxy": Kplateau,
        "K_peak_proxy": float(np.max(K)),
        "K_final_proxy": float(K[-1]),
        "delta_KR_proxy": max(Kplateau - float(K[0]), 0.0),
        "final_mobile_count": float(np.sum(front.mobile)),
        "final_retained_count": float(np.sum(front.retained)),
        "cumulative_emitted": float(front.cumulative_emitted),
        "cumulative_refreshed": float(front.cumulative_refreshed),
        "max_K_shield_MPa_sqrt_m": float(
            max(e["K_shield_MPa_sqrt_m"] for e in events)
        ),
        "max_sigma_back_Pa": float(max(e["sigma_back_max_Pa"] for e in events)),
        "events": events,
    }


def simulate_reduced_response(
    parameters: dict[str, float],
    T_K: float,
    settings: ReducedFrontSettings,
    *,
    mode: str = "full",
) -> dict[str, Any]:
    return ReducedCampaignFront(parameters, settings, mode=mode).run(float(T_K))


def _safe_span_fraction(values: np.ndarray, center: float) -> float:
    if values.size <= 1:
        return 0.0
    return float((np.max(values) - np.min(values)) / max(abs(center), 1.0e-12))


def best_adjacent_transition(
    temperatures_K: Iterable[float],
    toughness: Iterable[float],
    *,
    plasticity_off_toughness: Iterable[float] | None = None,
    requirements: TransitionRequirements = TransitionRequirements(),
) -> dict[str, Any]:
    """Find the best free adjacent DBTT split and return feature losses."""
    T = np.asarray(list(temperatures_K), dtype=float)
    K = np.asarray(list(toughness), dtype=float)
    order = np.argsort(T)
    T = T[order]
    K = K[order]
    if T.size < 2 * requirements.min_points_per_shelf:
        raise ValueError("insufficient temperatures for two DBTT shelves")
    if not np.all(np.isfinite(K)):
        return {"valid": False, "loss": 1.0e12, "reason": "nonfinite_toughness"}
    Koff = None
    if plasticity_off_toughness is not None:
        Koff = np.asarray(list(plasticity_off_toughness), dtype=float)[order]

    candidates: list[dict[str, Any]] = []
    lo_n = requirements.min_points_per_shelf
    hi_n = requirements.min_points_per_shelf
    for j in range(lo_n - 1, T.size - hi_n):
        low = K[: j + 1]
        high = K[j + 1 :]
        KL = float(np.median(low))
        KH = float(np.median(high))
        delta = KH - KL
        ratio = KH / max(KL, 1.0e-12)
        robust_ratio = float(np.min(high) / max(np.max(low), 1.0e-12))
        main_jump = float(K[j + 1] - K[j])
        concentration = main_jump / max(delta, 1.0e-12)
        low_span = _safe_span_fraction(low, KL)
        high_span = _safe_span_fraction(high, KH)
        outside_jumps = np.diff(K).copy()
        outside_jumps = np.delete(outside_jumps, j)
        max_secondary = float(np.max(np.abs(outside_jumps))) if outside_jumps.size else 0.0
        secondary_ratio = max_secondary / max(abs(main_jump), 1.0e-12)
        off_ratio = 1.0
        if Koff is not None and np.all(np.isfinite(Koff)):
            off_ratio = float(np.median(Koff[j + 1 :]) / max(np.median(Koff[: j + 1]), 1.0e-12))

        penalties = {
            "low_floor": max(requirements.low_shelf_min - KL, 0.0) / 2.0,
            "low_ceiling": max(KL - requirements.low_shelf_max, 0.0) / 3.0,
            "high_ceiling": max(KH - requirements.high_shelf_max, 0.0) / 5.0,
            "ratio": max(requirements.min_ratio - ratio, 0.0) / 0.25,
            "robust_ratio": max(requirements.robust_ratio - robust_ratio, 0.0) / 0.25,
            "low_flatness": max(low_span - requirements.max_low_span_fraction, 0.0) / 0.05,
            "high_flatness": max(high_span - requirements.max_high_span_fraction, 0.0) / 0.05,
            "jump_concentration": max(requirements.min_jump_concentration - concentration, 0.0) / 0.10,
            "secondary_jump": max(secondary_ratio - 0.50, 0.0) / 0.20,
            "plasticity_off_ratio": max(off_ratio - requirements.max_plasticity_off_ratio, 0.0) / 0.10,
        }
        loss = float(sum(value * value for value in penalties.values()))
        candidates.append(
            {
                "valid": True,
                "loss": loss,
                "split_index": int(j),
                "transition_low_K": float(T[j]),
                "transition_high_K": float(T[j + 1]),
                "low_shelf": KL,
                "high_shelf": KH,
                "shelf_ratio": ratio,
                "robust_shelf_ratio": robust_ratio,
                "main_jump": main_jump,
                "jump_concentration": concentration,
                "low_span_fraction": low_span,
                "high_span_fraction": high_span,
                "secondary_jump_ratio": secondary_ratio,
                "plasticity_off_ratio": off_ratio,
                "penalties": penalties,
            }
        )
    return min(candidates, key=lambda row: row["loss"])


__all__ = [
    "ReducedCampaignFront",
    "ReducedFrontSettings",
    "TransitionRequirements",
    "best_adjacent_transition",
    "cleavage_effective_rate_s",
    "exact_depletion",
    "exact_refresh",
    "exp_floor_barrier_eV",
    "simulate_reduced_response",
]
