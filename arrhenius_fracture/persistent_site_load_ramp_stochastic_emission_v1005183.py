"""Exact stochastic emission integrated over the accepted FEM load ramp.

v10.0.5.18.3 removes the deterministic state-change predictor introduced in
v10.0.5.18.2.  The global FEM solver supplies the start and end K values for an
accepted mechanics interval.  Cleavage, stochastic signed emission, transport,
backstress, blunting, source area, and site multiplicity are then integrated by
cheap front-local substeps along the linear K(t) ramp.

The stochastic event law and four-class parameters are unchanged.  No expected
emission increment is used to reject a global FEM trial.
"""
from __future__ import annotations

import copy
import math
import os
from typing import Any

import numpy as np

from .persistent_site_stochastic_emission_v100518 import (
    EMISSION_HAZARD_SCHEMA,
    EMISSION_TRANSPORT_SCHEMA,
    PersistentSiteStochasticEmissionMovingTipFrontEngineV100518,
    _emission_hazards,
)

MODEL_ID = "FEM_CZM_exact_stochastic_emission_load_ramp_v10_0_5_18_3"
LOAD_RAMP_SCHEMA = "v10.0.5.18.3_linear_K_ramp_local_exact_stochastic_kinetics"
PROBE_FALLBACK_SCHEMA = "v10.0.5.18.3_bounded_last_reliable_two_channel_probe"


def _linear_value(start: float, end: float, alpha: float) -> float:
    a = min(max(float(alpha), 0.0), 1.0)
    return (1.0 - a) * float(start) + a * float(end)


def _maximum_log_change(initial: np.ndarray, final: np.ndarray) -> float:
    a = np.asarray(initial, dtype=float).reshape(-1)
    b = np.asarray(final, dtype=float).reshape(-1)
    if a.shape != b.shape:
        raise ValueError("rate arrays must have the same shape")
    active = np.maximum(a, b) > 1.0e-300
    if not np.any(active):
        return 0.0
    la = np.log(np.maximum(a[active], 1.0e-300))
    lb = np.log(np.maximum(b[active], 1.0e-300))
    return float(np.max(np.abs(lb - la)))


def two_channel_absolute_opening_drives_v1005183(self, fallback_K: float):
    """Apply the resolved signed-channel anisotropy exactly once.

    A temporarily unavailable tensor probe may use the engine's bounded cached
    reliable metadata.  The fallback is numerical continuity only; the signed
    event hazards and material parameters are not modified.
    """

    KJ = max(float(fallback_K), 0.0)
    latest = dict(getattr(getattr(self, "_mm", None), "latest", {}) or {})
    reliable = bool(latest.get("two_channel_drive_reliable", False))
    if not reliable:
        cached = dict(getattr(self, "_last_reliable_two_channel_metadata", {}) or {})
        if cached:
            latest = {**latest, **cached}
            latest["two_channel_probe_fallback_active"] = True
            latest["two_channel_probe_fallback_schema"] = PROBE_FALLBACK_SCHEMA
    fc = float(latest.get("cleavage_factor", 1.0))
    fe = float(latest.get("emission_factor", 1.0))
    resolved = bool(latest.get("two_channel_drive_reliable", False))
    use_absolute = bool(
        getattr(self, "persistent_site_source_active", False) and resolved
    )
    K_cleave = KJ * fc
    K_emit = KJ if use_absolute else KJ * fe
    metadata = {
        "KJ": KJ,
        "fc": fc,
        "fe": fe,
        **latest,
        "persistent_two_channel_absolute_opening_drive": use_absolute,
        "scalar_emission_factor_applied_to_Kemit": not use_absolute,
        "two_channel_factor_applied_inside_emission_hazard": use_absolute,
    }
    return K_cleave, K_emit, metadata


class PersistentSiteLoadRampStochasticEmissionFrontEngineV1005183(
    PersistentSiteStochasticEmissionMovingTipFrontEngineV100518
):
    """Exact stochastic events with local integration along K(t)."""

    load_ramp_local_kinetics_active = True
    _ramp_max_log_rate_change_default = 0.5
    _ramp_action_floor_default = 1.0e-6
    _probe_fallback_max_calls_default = 8

    @classmethod
    def configure_stochastic(cls, **kwargs) -> None:
        super().configure_stochastic(**kwargs)
        cls._ramp_max_log_rate_change_default = max(
            float(os.environ.get("EMISSION_RAMP_MAX_LOG_RATE_CHANGE", "0.5")),
            1.0e-6,
        )
        cls._ramp_action_floor_default = max(
            float(os.environ.get("EMISSION_RAMP_ACTION_FLOOR", "1e-6")),
            1.0e-15,
        )
        cls._probe_fallback_max_calls_default = max(
            int(os.environ.get("TWO_CHANNEL_PROBE_FALLBACK_MAX_CALLS", "8")),
            0,
        )

    def reset(self):
        super().reset()
        self.load_ramp_K_cleave_prev_Pa_sqrt_m = 0.0
        self.load_ramp_K_emit_prev_Pa_sqrt_m = 0.0
        self.load_ramp_substeps_total = 0
        self.load_ramp_refinements_total = 0
        self.load_ramp_last_diagnostics: dict[str, Any] = {}
        self.ramp_max_log_rate_change = float(
            type(self)._ramp_max_log_rate_change_default
        )
        self.ramp_action_floor = float(type(self)._ramp_action_floor_default)
        self.probe_fallback_max_calls = int(
            type(self)._probe_fallback_max_calls_default
        )
        self._last_reliable_two_channel_factors: np.ndarray | None = None
        self._last_reliable_two_channel_tau: np.ndarray | None = None
        self._last_reliable_two_channel_metadata: dict[str, Any] = {}
        self.two_channel_probe_fallback_calls = 0
        self.two_channel_probe_fallback_total = 0

    def _capture_state(self) -> dict[str, Any]:
        snapshot = super()._capture_state()
        if hasattr(self, "load_ramp_K_cleave_prev_Pa_sqrt_m"):
            snapshot.update(
                {
                    "load_ramp_K_cleave_prev_Pa_sqrt_m": float(
                        self.load_ramp_K_cleave_prev_Pa_sqrt_m
                    ),
                    "load_ramp_K_emit_prev_Pa_sqrt_m": float(
                        self.load_ramp_K_emit_prev_Pa_sqrt_m
                    ),
                    "load_ramp_substeps_total": int(self.load_ramp_substeps_total),
                    "load_ramp_refinements_total": int(
                        self.load_ramp_refinements_total
                    ),
                    "load_ramp_last_diagnostics": copy.deepcopy(
                        self.load_ramp_last_diagnostics
                    ),
                    "last_reliable_two_channel_factors": (
                        None
                        if self._last_reliable_two_channel_factors is None
                        else self._last_reliable_two_channel_factors.copy()
                    ),
                    "last_reliable_two_channel_tau": (
                        None
                        if self._last_reliable_two_channel_tau is None
                        else self._last_reliable_two_channel_tau.copy()
                    ),
                    "last_reliable_two_channel_metadata": copy.deepcopy(
                        self._last_reliable_two_channel_metadata
                    ),
                    "two_channel_probe_fallback_calls": int(
                        self.two_channel_probe_fallback_calls
                    ),
                    "two_channel_probe_fallback_total": int(
                        self.two_channel_probe_fallback_total
                    ),
                }
            )
        return snapshot

    def _restore_state(self, snapshot: dict[str, Any]) -> None:
        super()._restore_state(snapshot)
        if "load_ramp_K_cleave_prev_Pa_sqrt_m" not in snapshot:
            return
        self.load_ramp_K_cleave_prev_Pa_sqrt_m = float(
            snapshot["load_ramp_K_cleave_prev_Pa_sqrt_m"]
        )
        self.load_ramp_K_emit_prev_Pa_sqrt_m = float(
            snapshot["load_ramp_K_emit_prev_Pa_sqrt_m"]
        )
        self.load_ramp_substeps_total = int(snapshot["load_ramp_substeps_total"])
        self.load_ramp_refinements_total = int(
            snapshot["load_ramp_refinements_total"]
        )
        self.load_ramp_last_diagnostics = copy.deepcopy(
            snapshot["load_ramp_last_diagnostics"]
        )
        factors = snapshot["last_reliable_two_channel_factors"]
        tau = snapshot["last_reliable_two_channel_tau"]
        self._last_reliable_two_channel_factors = (
            None if factors is None else np.asarray(factors, dtype=float).copy()
        )
        self._last_reliable_two_channel_tau = (
            None if tau is None else np.asarray(tau, dtype=float).copy()
        )
        self._last_reliable_two_channel_metadata = copy.deepcopy(
            snapshot["last_reliable_two_channel_metadata"]
        )
        self.two_channel_probe_fallback_calls = int(
            snapshot["two_channel_probe_fallback_calls"]
        )
        self.two_channel_probe_fallback_total = int(
            snapshot["two_channel_probe_fallback_total"]
        )

    def _two_channel_drive(self) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
        context = getattr(self, "_mm", None)
        latest = {} if context is None else dict(getattr(context, "latest", {}) or {})
        reliable = bool(latest.get("two_channel_drive_reliable", False))
        factors = np.asarray(
            latest.get("two_channel_drive_factors", []), dtype=float
        ).reshape(-1)
        tau = np.asarray(
            latest.get("two_channel_tau_signed_Pa", []), dtype=float
        ).reshape(-1)
        valid = bool(
            reliable
            and factors.shape == (2,)
            and tau.shape == (2,)
            and np.all(np.isfinite(factors))
            and np.all(np.isfinite(tau))
        )
        if valid:
            self._last_reliable_two_channel_factors = factors.copy()
            self._last_reliable_two_channel_tau = tau.copy()
            self._last_reliable_two_channel_metadata = copy.deepcopy(latest)
            self.two_channel_probe_fallback_calls = 0
            return factors, tau, latest

        cached = bool(
            self._last_reliable_two_channel_factors is not None
            and self._last_reliable_two_channel_tau is not None
        )
        if cached and self.two_channel_probe_fallback_calls < self.probe_fallback_max_calls:
            self.two_channel_probe_fallback_calls += 1
            self.two_channel_probe_fallback_total += 1
            metadata = {
                **latest,
                **copy.deepcopy(self._last_reliable_two_channel_metadata),
                "two_channel_probe_fallback_active": True,
                "two_channel_probe_fallback_calls": int(
                    self.two_channel_probe_fallback_calls
                ),
                "two_channel_probe_fallback_schema": PROBE_FALLBACK_SCHEMA,
            }
            return (
                self._last_reliable_two_channel_factors.copy(),
                self._last_reliable_two_channel_tau.copy(),
                metadata,
            )
        raise RuntimeError(
            "persistent signed emission could not obtain a reliable two-channel "
            "FEM tensor drive and the bounded numerical fallback was exhausted"
        )

    def _ramp_values(
        self,
        *,
        consumed_s: float,
        total_s: float,
        offset_s: float,
        K_cleave_end: float,
        K_emit_end: float,
    ) -> tuple[float, float]:
        if total_s <= 0.0:
            return float(K_cleave_end), float(K_emit_end)
        alpha = (float(consumed_s) + float(offset_s)) / float(total_s)
        return (
            _linear_value(
                self.load_ramp_K_cleave_prev_Pa_sqrt_m,
                K_cleave_end,
                alpha,
            ),
            _linear_value(
                self.load_ramp_K_emit_prev_Pa_sqrt_m,
                K_emit_end,
                alpha,
            ),
        )

    def _local_ramp_substep(
        self,
        *,
        remaining_s: float,
        consumed_s: float,
        total_s: float,
        K_cleave_end: float,
        K_emit_end: float,
        T_K: float,
        drive_factors: np.ndarray,
        tau_signed_Pa: np.ndarray,
    ) -> tuple[float, dict[str, float]]:
        remaining = max(float(remaining_s), 0.0)
        Kc0, Ke0 = self._ramp_values(
            consumed_s=consumed_s,
            total_s=total_s,
            offset_s=0.0,
            K_cleave_end=K_cleave_end,
            K_emit_end=K_emit_end,
        )
        threshold = max(float(self.hazard_threshold_action), 1.0e-300)
        lam0, _, _ = self.lambda_cleave(self.sigma_tip(Kc0), T_K)
        lam0 = max(float(lam0), 0.0) if math.isfinite(lam0) else 0.0
        progress0 = lam0 / threshold
        h = self._substep_limit_stochastic(remaining, progress0)
        refinements = 0
        diagnostics = {
            "cleavage_action": 0.0,
            "emission_action": 0.0,
            "cleavage_log_change": 0.0,
            "emission_log_change": 0.0,
        }
        for _ in range(16):
            h = min(max(float(h), min(remaining, self.min_substep_s)), remaining)
            Kc1, Ke1 = self._ramp_values(
                consumed_s=consumed_s,
                total_s=total_s,
                offset_s=h,
                K_cleave_end=K_cleave_end,
                K_emit_end=K_emit_end,
            )
            lam1, _, _ = self.lambda_cleave(self.sigma_tip(Kc1), T_K)
            lam1 = max(float(lam1), 0.0) if math.isfinite(lam1) else 0.0
            progress1 = lam1 / threshold
            opening0 = self._opening_stress(Ke0)[0]
            opening1 = self._opening_stress(Ke1)[0]
            emission0 = _emission_hazards(
                self.mpz_state,
                T_K=T_K,
                opening_stress_Pa=opening0,
                drive_factors=drive_factors,
                tau_signed_Pa=tau_signed_Pa,
                rate_function=self._emission_rate_per_site,
            )["aggregate_hazard_s"]
            emission1 = _emission_hazards(
                self.mpz_state,
                T_K=T_K,
                opening_stress_Pa=opening1,
                drive_factors=drive_factors,
                tau_signed_Pa=tau_signed_Pa,
                rate_function=self._emission_rate_per_site,
            )["aggregate_hazard_s"]
            cleavage_action = 0.5 * (progress0 + progress1) * h
            emission_action = float(
                np.max(0.5 * (np.asarray(emission0) + np.asarray(emission1)) * h)
            )
            cleavage_log = _maximum_log_change(
                np.asarray([progress0]), np.asarray([progress1])
            )
            emission_log = _maximum_log_change(emission0, emission1)
            diagnostics = {
                "cleavage_action": float(cleavage_action),
                "emission_action": float(emission_action),
                "cleavage_log_change": float(cleavage_log),
                "emission_log_change": float(emission_log),
            }
            action_ratio = max(
                cleavage_action / max(float(self.max_action_substep), 1.0e-12),
                emission_action
                / max(float(self.emission_max_action_substep), 1.0e-12),
            )
            log_ratio = 0.0
            if cleavage_action > self.ramp_action_floor:
                log_ratio = max(
                    log_ratio,
                    cleavage_log / self.ramp_max_log_rate_change,
                )
            if emission_action > self.ramp_action_floor:
                log_ratio = max(
                    log_ratio,
                    emission_log / self.ramp_max_log_rate_change,
                )
            ratio = max(action_ratio, log_ratio)
            if ratio <= 1.0 + 1.0e-10 or h <= self.min_substep_s:
                break
            refinements += 1
            factor = min(0.5, max(0.05, 0.8 / max(ratio, 1.0)))
            h *= factor
        diagnostics["refinements"] = float(refinements)
        return h, diagnostics

    def _integrate_coupled(
        self,
        *,
        K_cleave: float,
        K_emit: float,
        T_K: float,
        dt_s: float,
        drive_factors: np.ndarray,
        tau_signed_Pa: np.ndarray,
    ) -> dict[str, Any]:
        self._synchronize_driver_checkpoint_length()
        dt_requested = max(float(dt_s), 0.0)
        remaining = dt_requested
        consumed = 0.0
        dB_total = 0.0
        dH_total = 0.0
        da_total = 0.0
        packet_mean = 0.0
        packet_variance = 0.0
        plastic_totals: dict[str, float] = {}
        advance_totals: dict[str, float] = {}
        fired = False
        microsteps = 0
        ramp_refinements = 0
        max_cleavage_action = 0.0
        max_emission_action = 0.0
        max_cleavage_log = 0.0
        max_emission_log = 0.0
        last_lambda = 0.0
        last_lambda_raw = 0.0
        last_Gc = 0.0
        last_sigma = self.sigma_tip(K_cleave)
        last_opening = self._opening_stress(K_emit)[0]
        completed_threshold = 0.0
        completed_action = 0.0
        completed_length = 0.0
        completed_factor = 0.0

        while remaining > 0.0:
            microsteps += 1
            if microsteps > int(self.max_internal_steps):
                raise RuntimeError(
                    "load-ramp stochastic moving-tip cell exceeded max_internal_steps; "
                    "the local kinetic ramp requires refinement"
                )

            threshold = max(float(self.hazard_threshold_action), 1.0e-300)
            event_length = float(self.stochastic_event_advance_m)
            event_factor = float(self.stochastic_event_length_factor)
            h, ramp_diag = self._local_ramp_substep(
                remaining_s=remaining,
                consumed_s=consumed,
                total_s=dt_requested,
                K_cleave_end=float(K_cleave),
                K_emit_end=float(K_emit),
                T_K=float(T_K),
                drive_factors=drive_factors,
                tau_signed_Pa=tau_signed_Pa,
            )
            ramp_refinements += int(ramp_diag["refinements"])
            max_cleavage_action = max(
                max_cleavage_action, ramp_diag["cleavage_action"]
            )
            max_emission_action = max(
                max_emission_action, ramp_diag["emission_action"]
            )
            max_cleavage_log = max(
                max_cleavage_log, ramp_diag["cleavage_log_change"]
            )
            max_emission_log = max(
                max_emission_log, ramp_diag["emission_log_change"]
            )
            microstep_start = self._capture_state()

            Kc_quarter, Ke_quarter = self._ramp_values(
                consumed_s=consumed,
                total_s=dt_requested,
                offset_s=0.25 * h,
                K_cleave_end=K_cleave,
                K_emit_end=K_emit,
            )
            Kc_mid, _ = self._ramp_values(
                consumed_s=consumed,
                total_s=dt_requested,
                offset_s=0.5 * h,
                K_cleave_end=K_cleave,
                K_emit_end=K_emit,
            )
            first, _ = self._plastic_half_step(
                dt_s=0.5 * h,
                T_K=T_K,
                K_emit=Ke_quarter,
                drive_factors=drive_factors,
                tau_signed_Pa=tau_signed_Pa,
            )
            sigma_mid = self.sigma_tip(Kc_mid)
            lambda_mid, raw_mid, Gc_mid = self.lambda_cleave(sigma_mid, T_K)
            lambda_mid = (
                max(float(lambda_mid), 0.0)
                if math.isfinite(lambda_mid)
                else 0.0
            )
            progress_mid = lambda_mid / threshold
            action_remaining = max(1.0 - float(self.B), 0.0)

            if (
                progress_mid > 0.0
                and progress_mid * h > action_remaining + 1.0e-12
            ):
                self._restore_state(microstep_start)
                threshold = max(float(self.hazard_threshold_action), 1.0e-300)
                event_length = float(self.stochastic_event_advance_m)
                event_factor = float(self.stochastic_event_length_factor)
                h = max(action_remaining / progress_mid, float(self.min_substep_s))
                h = min(h, remaining)
                Kc_quarter, Ke_quarter = self._ramp_values(
                    consumed_s=consumed,
                    total_s=dt_requested,
                    offset_s=0.25 * h,
                    K_cleave_end=K_cleave,
                    K_emit_end=K_emit,
                )
                Kc_mid, _ = self._ramp_values(
                    consumed_s=consumed,
                    total_s=dt_requested,
                    offset_s=0.5 * h,
                    K_cleave_end=K_cleave,
                    K_emit_end=K_emit,
                )
                first, _ = self._plastic_half_step(
                    dt_s=0.5 * h,
                    T_K=T_K,
                    K_emit=Ke_quarter,
                    drive_factors=drive_factors,
                    tau_signed_Pa=tau_signed_Pa,
                )
                sigma_mid = self.sigma_tip(Kc_mid)
                lambda_mid, raw_mid, Gc_mid = self.lambda_cleave(sigma_mid, T_K)
                lambda_mid = (
                    max(float(lambda_mid), 0.0)
                    if math.isfinite(lambda_mid)
                    else 0.0
                )
                progress_mid = lambda_mid / threshold

            action_remaining = max(1.0 - float(self.B), 0.0)
            dB = min(progress_mid * h, action_remaining)
            dH = dB * threshold
            da = event_length * dB
            advance = self.mpz_state.advance(da) if da > 0.0 else {}

            _, Ke_three_quarter = self._ramp_values(
                consumed_s=consumed,
                total_s=dt_requested,
                offset_s=0.75 * h,
                K_cleave_end=K_cleave,
                K_emit_end=K_emit,
            )
            second, last_opening = self._plastic_half_step(
                dt_s=0.5 * h,
                T_K=T_K,
                K_emit=Ke_three_quarter,
                drive_factors=drive_factors,
                tau_signed_Pa=tau_signed_Pa,
            )
            self._sum_numeric(plastic_totals, first)
            self._sum_numeric(plastic_totals, second)
            self._sum_numeric(advance_totals, advance)

            packet_rate = (
                event_length / max(float(self.b), 1.0e-30) * progress_mid
            )
            packet_n = packet_rate * h
            packet_var = float(self.b) ** 2 * packet_n
            self.B += dB
            self.hazard_action_current += dH
            self.micro_advance_total_m += da
            self.kinetic_packet_count_mean_total += packet_n
            self.kinetic_packet_variance_total_m2 += packet_var
            dB_total += dB
            dH_total += dH
            da_total += da
            packet_mean += packet_n
            packet_variance += packet_var
            consumed += h
            remaining = max(remaining - h, 0.0)
            self.t += h
            last_lambda = lambda_mid
            last_lambda_raw = float(raw_mid)
            last_Gc = float(Gc_mid)
            last_sigma = self.sigma_tip(Kc_mid)

            if self.B >= 1.0 - 1.0e-10:
                self.B = max(self.B - 1.0, 0.0)
                self.a_adv += event_length
                self.checkpoint_advance_total_m += event_length
                self.n_adv += 1
                fired = True
                completed_threshold = threshold
                completed_action = float(self.hazard_action_current)
                completed_length = event_length
                completed_factor = event_factor
                self.hazard_last_completed_threshold = completed_threshold
                self.hazard_last_completed_action = completed_action
                self.hazard_threshold_history.append(completed_threshold)
                self.stochastic_last_completed_advance_m = completed_length
                self.stochastic_last_completed_factor = completed_factor
                self.stochastic_event_length_history.append(completed_length)
                self.hazard_event_index += 1
                self.hazard_action_current = 0.0
                self.hazard_threshold_action = self._draw_threshold()
                self._set_current_event_length()
                break
            if h <= 0.0:
                break

        self.kinetic_internal_substeps_total += microsteps
        self.load_ramp_substeps_total += microsteps
        self.load_ramp_refinements_total += ramp_refinements
        self.load_ramp_last_diagnostics = {
            "local_substeps": int(microsteps),
            "local_refinements": int(ramp_refinements),
            "max_cleavage_action": float(max_cleavage_action),
            "max_emission_action": float(max_emission_action),
            "max_cleavage_log_rate_change": float(max_cleavage_log),
            "max_emission_log_rate_change": float(max_emission_log),
            "K_cleave_start_Pa_sqrt_m": float(
                self.load_ramp_K_cleave_prev_Pa_sqrt_m
            ),
            "K_cleave_end_Pa_sqrt_m": float(K_cleave),
            "K_emit_start_Pa_sqrt_m": float(
                self.load_ramp_K_emit_prev_Pa_sqrt_m
            ),
            "K_emit_end_Pa_sqrt_m": float(K_emit),
        }
        return {
            "fired": fired,
            "n_fire": 1 if fired else 0,
            "v_crack": da_total / consumed if consumed > 0.0 else 0.0,
            "dB": dB_total,
            "physical_hazard_action_step": dH_total,
            "da": da_total,
            "dt_consumed": consumed,
            "dt_unused": max(dt_requested - consumed, 0.0),
            "packet_mean": packet_mean,
            "packet_variance_m2": packet_variance,
            "lambda_c": last_lambda,
            "lambda_c_raw": last_lambda_raw,
            "Gc_J": last_Gc,
            "sigma_tip": last_sigma,
            "sigma_emit_tip": last_opening,
            "plastic": plastic_totals,
            "advance": advance_totals,
            "microsteps": microsteps,
            "load_ramp_refinements": ramp_refinements,
            "load_ramp_diagnostics": copy.deepcopy(
                self.load_ramp_last_diagnostics
            ),
            "hazard_threshold_completed_action": completed_threshold,
            "hazard_action_completed": completed_action,
            "hazard_threshold_next_action": float(self.hazard_threshold_action),
            "geometry_event_advance_m": completed_length,
            "stochastic_event_advance_m": completed_length,
            "stochastic_event_length_factor": completed_factor,
            "stochastic_current_event_advance_m": float(
                self.stochastic_event_advance_m
            ),
            "stochastic_current_event_length_factor": float(
                self.stochastic_event_length_factor
            ),
        }

    def step_drives(self, K_cleave, K_emit, T, dt, metadata=None):
        out = super().step_drives(
            K_cleave,
            K_emit,
            T,
            dt,
            metadata=metadata,
        )
        self.load_ramp_K_cleave_prev_Pa_sqrt_m = float(K_cleave)
        self.load_ramp_K_emit_prev_Pa_sqrt_m = float(K_emit)
        out.update(
            {
                "load_ramp_local_kinetics_active": True,
                "load_ramp_schema": LOAD_RAMP_SCHEMA,
                "load_ramp_deterministic_state_predictor_active": False,
                "load_ramp_local_substeps_total": int(
                    self.load_ramp_substeps_total
                ),
                "load_ramp_local_refinements_total": int(
                    self.load_ramp_refinements_total
                ),
                "load_ramp_last_diagnostics": copy.deepcopy(
                    self.load_ramp_last_diagnostics
                ),
                "two_channel_probe_fallback_total": int(
                    self.two_channel_probe_fallback_total
                ),
            }
        )
        return out

    @classmethod
    def audit_payload(cls) -> dict[str, Any]:
        payload = super().audit_payload()
        payload["model_id"] = MODEL_ID
        emission = dict(payload.get("stochastic_emission", {}) or {})
        emission.update(
            {
                "load_ramp_schema": LOAD_RAMP_SCHEMA,
                "global_FEM_state_change_predictor": False,
                "local_linear_K_ramp": True,
                "accepted_update": "exact_event_localized_one_activation_packets",
                "post_onset_mean_field_burst": False,
                "backstress_recomputed_after_each_event": True,
                "tip_radius_recomputed_after_each_event": True,
                "front_width_recomputed_after_each_event": True,
                "site_multiplicity_recomputed_after_each_event": True,
                "signed_channel_factors_applied_once": True,
                "bounded_probe_fallback_schema": PROBE_FALLBACK_SCHEMA,
            }
        )
        payload["stochastic_emission"] = emission
        return payload


__all__ = [
    "EMISSION_HAZARD_SCHEMA",
    "EMISSION_TRANSPORT_SCHEMA",
    "LOAD_RAMP_SCHEMA",
    "MODEL_ID",
    "PROBE_FALLBACK_SCHEMA",
    "PersistentSiteLoadRampStochasticEmissionFrontEngineV1005183",
    "two_channel_absolute_opening_drives_v1005183",
]
