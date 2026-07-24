"""PF stochastic cleavage threshold and event-length parity for FEM/CZM.

This release retains the v10.0.5.15 PF local update map and continuous moving-tip
Strang coupling.  It adds the exact PF v10.1.7.2/v10.1.7.3 renewal semantics:

* Xi ~ Exponential(1) is the integrated-hazard threshold for one event;
* normalized progress evolves at lambda_c / Xi;
* the same Xi sets the event reward through a bounded mean-preserving map;
* no noise is added to K, J, barriers, shielding, backstress, or material data.
"""
from __future__ import annotations

import copy
import math
from typing import Any

import numpy as np

from .persistent_site_moving_tip_v100515 import (
    COUPLING_SCHEME,
    PersistentSitePFMovingTipFrontEngineV100515,
)
from .stochastic_geometry_v100516 import enqueue_stochastic_geometry_event

MODEL_ID = "FEM_CZM_PF_stochastic_hazard_event_length_moving_tip_v10_0_5_16"
HAZARD_SCHEMA = "PF_v10_1_7_2_exponential_integrated_hazard_threshold"
EVENT_LENGTH_SCHEMA = "PF_v10_1_7_3_threshold_correlated_mean_preserving_event_length"


def draw_hazard_threshold(
    mode: str,
    rng: np.random.Generator,
    minimum_threshold: float = 1.0e-12,
) -> float:
    selected = str(mode).strip().lower()
    floor = max(float(minimum_threshold), 1.0e-300)
    if selected == "deterministic":
        return 1.0
    if selected != "exponential":
        raise ValueError("hazard threshold mode must be deterministic or exponential")
    return max(float(rng.exponential(1.0)), floor)


def clipped_exponential_mean(minimum_factor: float, maximum_factor: float) -> float:
    """Return E[clip(X,a,b)] for X~Exponential(1), matching PF v10.1.7.3."""

    a = max(float(minimum_factor), 0.0)
    b = max(float(maximum_factor), a)
    return max(a + math.exp(-a) - math.exp(-b), 1.0e-300)


def threshold_event_length_factor(
    threshold_action: float,
    *,
    mode: str = "threshold_scaled",
    minimum_factor: float = 0.5,
    maximum_factor: float = 4.0,
    deterministic_threshold: bool = False,
) -> float:
    selected = str(mode).strip().lower()
    if selected == "fixed" or deterministic_threshold:
        return 1.0
    if selected != "threshold_scaled":
        raise ValueError("event length mode must be fixed or threshold_scaled")
    a = max(float(minimum_factor), 1.0e-12)
    b = max(float(maximum_factor), a)
    clipped = min(max(float(threshold_action), a), b)
    return clipped / clipped_exponential_mean(a, b)


class PersistentSitePFStochasticMovingTipFrontEngineV100516(
    PersistentSitePFMovingTipFrontEngineV100515
):
    """v10.0.5.15 engine with PF-correlated stochastic renewal rewards."""

    stochastic_hazard_threshold_active = True
    stochastic_event_length_active = True
    _hazard_mode_default = "exponential"
    _hazard_seed_default = 0
    _hazard_minimum_threshold_default = 1.0e-12
    _event_length_mode_default = "threshold_scaled"
    _event_minimum_factor_default = 0.5
    _event_maximum_factor_default = 4.0

    @classmethod
    def configure_stochastic(
        cls,
        *,
        hazard_mode: str = "exponential",
        hazard_seed: int = 0,
        hazard_minimum_threshold: float = 1.0e-12,
        event_length_mode: str = "threshold_scaled",
        event_minimum_factor: float = 0.5,
        event_maximum_factor: float = 4.0,
    ) -> None:
        hazard = str(hazard_mode).strip().lower()
        length = str(event_length_mode).strip().lower()
        if hazard not in {"deterministic", "exponential"}:
            raise ValueError("hazard mode must be deterministic or exponential")
        if length not in {"fixed", "threshold_scaled"}:
            raise ValueError("event length mode must be fixed or threshold_scaled")
        if int(hazard_seed) < 0:
            raise ValueError("hazard seed must be nonnegative")
        minimum = max(float(event_minimum_factor), 1.0e-12)
        maximum = max(float(event_maximum_factor), minimum)
        cls._hazard_mode_default = hazard
        cls._hazard_seed_default = int(hazard_seed)
        cls._hazard_minimum_threshold_default = max(
            float(hazard_minimum_threshold), 1.0e-300
        )
        cls._event_length_mode_default = length
        cls._event_minimum_factor_default = minimum
        cls._event_maximum_factor_default = maximum

    def reset(self):
        super().reset()
        self.hazard_mode = str(type(self)._hazard_mode_default)
        self.hazard_seed = int(type(self)._hazard_seed_default)
        self.hazard_minimum_threshold = float(
            type(self)._hazard_minimum_threshold_default
        )
        self.event_length_mode = str(type(self)._event_length_mode_default)
        self.event_minimum_factor = float(type(self)._event_minimum_factor_default)
        self.event_maximum_factor = float(type(self)._event_maximum_factor_default)
        sequence = np.random.SeedSequence(
            [self.hazard_seed, int(getattr(self, "_engine_id", 0))]
        )
        self._hazard_rng = np.random.default_rng(sequence)
        self.hazard_threshold_action = self._draw_threshold()
        self.hazard_action_current = 0.0
        self.hazard_event_index = 0
        self.hazard_last_completed_threshold = 0.0
        self.hazard_last_completed_action = 0.0
        self.hazard_threshold_history: list[float] = []
        self.stochastic_base_checkpoint_m = max(float(self.f.da), 1.0e-30)
        self.stochastic_checkpoint_synchronized = False
        self.stochastic_event_length_factor = 1.0
        self.stochastic_event_advance_m = self.stochastic_base_checkpoint_m
        self.stochastic_last_completed_factor = 0.0
        self.stochastic_last_completed_advance_m = 0.0
        self.stochastic_event_length_history: list[float] = []
        self._set_current_event_length()

    def _draw_threshold(self) -> float:
        return draw_hazard_threshold(
            self.hazard_mode,
            self._hazard_rng,
            self.hazard_minimum_threshold,
        )

    def _set_current_event_length(self) -> None:
        factor = threshold_event_length_factor(
            self.hazard_threshold_action,
            mode=self.event_length_mode,
            minimum_factor=self.event_minimum_factor,
            maximum_factor=self.event_maximum_factor,
            deterministic_threshold=(self.hazard_mode == "deterministic"),
        )
        self.stochastic_event_length_factor = float(factor)
        self.stochastic_event_advance_m = (
            float(self.stochastic_base_checkpoint_m) * float(factor)
        )

    def _synchronize_driver_checkpoint_length(self) -> None:
        configured = max(float(self.f.da), 1.0e-30)
        same = math.isclose(
            configured,
            float(self.stochastic_base_checkpoint_m),
            rel_tol=1.0e-12,
            abs_tol=1.0e-18,
        )
        if self.stochastic_checkpoint_synchronized and same:
            return
        event_started = (
            bool(self.stochastic_event_length_history)
            or int(self.hazard_event_index) > 0
            or abs(float(self.B)) > 1.0e-12
            or abs(float(self.hazard_action_current)) > 1.0e-12
        )
        if event_started and not same:
            raise RuntimeError(
                "nominal checkpoint length changed after stochastic event evolution "
                f"began: old={self.stochastic_base_checkpoint_m:.9e} m, "
                f"new={configured:.9e} m"
            )
        self.stochastic_base_checkpoint_m = configured
        self.stochastic_checkpoint_synchronized = True
        self._set_current_event_length()

    def _capture_state(self) -> dict[str, Any]:
        snapshot = super()._capture_state()
        if hasattr(self, "_hazard_rng"):
            snapshot.update(
                {
                    "hazard_rng_state": copy.deepcopy(
                        self._hazard_rng.bit_generator.state
                    ),
                    "hazard_threshold_action": float(self.hazard_threshold_action),
                    "hazard_action_current": float(self.hazard_action_current),
                    "hazard_event_index": int(self.hazard_event_index),
                    "hazard_last_completed_threshold": float(
                        self.hazard_last_completed_threshold
                    ),
                    "hazard_last_completed_action": float(
                        self.hazard_last_completed_action
                    ),
                    "hazard_threshold_history": list(self.hazard_threshold_history),
                    "stochastic_base_checkpoint_m": float(
                        self.stochastic_base_checkpoint_m
                    ),
                    "stochastic_checkpoint_synchronized": bool(
                        self.stochastic_checkpoint_synchronized
                    ),
                    "stochastic_event_length_factor": float(
                        self.stochastic_event_length_factor
                    ),
                    "stochastic_event_advance_m": float(
                        self.stochastic_event_advance_m
                    ),
                    "stochastic_last_completed_factor": float(
                        self.stochastic_last_completed_factor
                    ),
                    "stochastic_last_completed_advance_m": float(
                        self.stochastic_last_completed_advance_m
                    ),
                    "stochastic_event_length_history": list(
                        self.stochastic_event_length_history
                    ),
                }
            )
        return snapshot

    def _restore_state(self, snapshot: dict[str, Any]) -> None:
        super()._restore_state(snapshot)
        if "hazard_rng_state" not in snapshot:
            return
        rng = np.random.default_rng()
        rng.bit_generator.state = copy.deepcopy(snapshot["hazard_rng_state"])
        self._hazard_rng = rng
        for name in (
            "hazard_threshold_action",
            "hazard_action_current",
            "hazard_last_completed_threshold",
            "hazard_last_completed_action",
            "stochastic_base_checkpoint_m",
            "stochastic_event_length_factor",
            "stochastic_event_advance_m",
            "stochastic_last_completed_factor",
            "stochastic_last_completed_advance_m",
        ):
            setattr(self, name, float(snapshot[name]))
        self.hazard_event_index = int(snapshot["hazard_event_index"])
        self.hazard_threshold_history = list(snapshot["hazard_threshold_history"])
        self.stochastic_checkpoint_synchronized = bool(
            snapshot["stochastic_checkpoint_synchronized"]
        )
        self.stochastic_event_length_history = list(
            snapshot["stochastic_event_length_history"]
        )

    def _substep_limit_stochastic(
        self, dt_remaining: float, progress_rate_s: float
    ) -> float:
        remaining = max(float(dt_remaining), 0.0)
        h = remaining
        rate = max(float(progress_rate_s), 0.0)
        if rate > 0.0 and math.isfinite(rate):
            h = min(h, float(self.max_action_substep) / rate)
            translation_rate = max(
                float(self.stochastic_event_advance_m) * rate, 1.0e-300
            )
            h = min(h, float(self.max_translation_substep_m) / translation_rate)
            action_remaining = max(1.0 - float(self.B), 0.0)
            if action_remaining > 0.0:
                h = min(h, action_remaining / rate)
        return max(min(h, remaining), min(float(self.min_substep_s), remaining))

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
                    "PF stochastic moving-tip cell exceeded max_internal_steps; "
                    "reduce the outer timestep without changing the hazard law"
                )

            threshold = max(float(self.hazard_threshold_action), 1.0e-300)
            event_length = float(self.stochastic_event_advance_m)
            event_factor = float(self.stochastic_event_length_factor)
            sigma0 = self.sigma_tip(K_cleave)
            lambda0, _, _ = self.lambda_cleave(sigma0, T_K)
            lambda0 = max(float(lambda0), 0.0) if math.isfinite(lambda0) else 0.0
            progress0 = lambda0 / threshold
            h = self._substep_limit_stochastic(remaining, progress0)
            microstep_start = self._capture_state()

            first, _ = self._plastic_half_step(
                dt_s=0.5 * h,
                T_K=T_K,
                K_emit=K_emit,
                drive_factors=drive_factors,
                tau_signed_Pa=tau_signed_Pa,
            )
            sigma_mid = self.sigma_tip(K_cleave)
            lambda_mid, raw_mid, Gc_mid = self.lambda_cleave(sigma_mid, T_K)
            lambda_mid = (
                max(float(lambda_mid), 0.0) if math.isfinite(lambda_mid) else 0.0
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
                first, _ = self._plastic_half_step(
                    dt_s=0.5 * h,
                    T_K=T_K,
                    K_emit=K_emit,
                    drive_factors=drive_factors,
                    tau_signed_Pa=tau_signed_Pa,
                )
                sigma_mid = self.sigma_tip(K_cleave)
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

            second, last_opening = self._plastic_half_step(
                dt_s=0.5 * h,
                T_K=T_K,
                K_emit=K_emit,
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
            last_sigma = self.sigma_tip(K_cleave)

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

    def predict_clock_increment_drives(self, K_cleave, K_emit, T, dt):
        factors, tau_signed, _ = self._two_channel_drive()
        trial = copy.copy(self)
        trial.mpz_state = self.mpz_state.copy()
        trial._checkpoint_origin_snapshot = None
        trial._geometry_veto_snapshot = None
        trial._checkpoint_fired_last_step = False
        rng = np.random.default_rng()
        rng.bit_generator.state = copy.deepcopy(self._hazard_rng.bit_generator.state)
        trial._hazard_rng = rng
        predicted = trial._integrate_coupled(
            K_cleave=float(K_cleave),
            K_emit=float(K_emit),
            T_K=float(T),
            dt_s=max(float(dt), 0.0),
            drive_factors=factors,
            tau_signed_Pa=tau_signed,
        )
        self.kinetic_prediction_calls += 1
        return float(max(predicted["dB"], 0.0))

    def step_drives(self, K_cleave, K_emit, T, dt, metadata=None):
        out = super().step_drives(K_cleave, K_emit, T, dt, metadata=metadata)
        fired = bool(out.get("fired", False))
        event_length = (
            float(self.stochastic_last_completed_advance_m) if fired else 0.0
        )
        event_factor = (
            float(self.stochastic_last_completed_factor) if fired else 0.0
        )
        completed_threshold = (
            float(self.hazard_last_completed_threshold) if fired else 0.0
        )
        completed_action = (
            float(self.hazard_last_completed_action) if fired else 0.0
        )
        out.update(
            {
                "stochastic_hazard_schema": HAZARD_SCHEMA,
                "stochastic_event_length_schema": EVENT_LENGTH_SCHEMA,
                "stochastic_hazard_enabled": self.hazard_mode == "exponential",
                "hazard_threshold_mode": self.hazard_mode,
                "hazard_seed": int(self.hazard_seed),
                "hazard_event_index": int(self.hazard_event_index),
                "hazard_threshold_completed_action": completed_threshold,
                "hazard_action_completed": completed_action,
                "hazard_threshold_next_action": float(
                    self.hazard_threshold_action
                ),
                "stochastic_event_length_enabled": (
                    self.event_length_mode == "threshold_scaled"
                ),
                "stochastic_event_length_mode": self.event_length_mode,
                "stochastic_event_minimum_factor": float(
                    self.event_minimum_factor
                ),
                "stochastic_event_maximum_factor": float(
                    self.event_maximum_factor
                ),
                "stochastic_base_checkpoint_m": float(
                    self.stochastic_base_checkpoint_m
                ),
                "stochastic_event_advance_m": event_length,
                "geometry_event_advance_m": event_length,
                "stochastic_event_length_factor": event_factor,
                "stochastic_current_event_advance_m": float(
                    self.stochastic_event_advance_m
                ),
                "stochastic_current_event_length_factor": float(
                    self.stochastic_event_length_factor
                ),
                "kinetic_checkpoint_progress_m": float(
                    self.B * self.stochastic_event_advance_m
                ),
                "mean_event_length_preserved": True,
                "event_length_uses_same_hazard_threshold": True,
                "noise_added_to_K": False,
                "noise_added_to_barriers": False,
            }
        )
        if fired:
            enqueue_stochastic_geometry_event(
                {
                    "event_advance_m": event_length,
                    "event_length_factor": event_factor,
                    "threshold_action": completed_threshold,
                    "hazard_action_completed": completed_action,
                    "hazard_seed": int(self.hazard_seed),
                    "hazard_event_index": int(self.hazard_event_index - 1),
                    "nominal_checkpoint_m": float(
                        self.stochastic_base_checkpoint_m
                    ),
                }
            )
        return out

    @classmethod
    def audit_payload(cls) -> dict[str, Any]:
        payload = super().audit_payload()
        payload.update(
            {
                "model_id": MODEL_ID,
                "coupling_scheme": COUPLING_SCHEME,
                "stochastic_hazard": {
                    "schema": HAZARD_SCHEMA,
                    "mode": cls._hazard_mode_default,
                    "seed": cls._hazard_seed_default,
                    "minimum_threshold": cls._hazard_minimum_threshold_default,
                    "distribution": "exponential_unit_mean",
                },
                "stochastic_event_length": {
                    "schema": EVENT_LENGTH_SCHEMA,
                    "mode": cls._event_length_mode_default,
                    "minimum_factor": cls._event_minimum_factor_default,
                    "maximum_factor": cls._event_maximum_factor_default,
                    "same_threshold_as_waiting_time": True,
                    "mean_length_preserved": True,
                    "geometry_realization": "single_checked_adaptive_CZM_commit",
                },
                "noise_added_to_K": False,
                "noise_added_to_J": False,
                "noise_added_to_barriers": False,
                "noise_added_to_material_parameters": False,
            }
        )
        return payload


__all__ = [
    "MODEL_ID",
    "HAZARD_SCHEMA",
    "EVENT_LENGTH_SCHEMA",
    "PersistentSitePFStochasticMovingTipFrontEngineV100516",
    "clipped_exponential_mean",
    "draw_hazard_threshold",
    "threshold_event_length_factor",
]
