"""Joint linear-K ramp and exact event-driven signed emission for v10.0.5.18.3.2.

The outer moving-tip cell controls only cleavage action and fractional crack-tip
translation.  The signed-emission solver receives the actual start and end
``K_emit`` values for each plastic half-step, evaluates the stochastic hazards
at the corresponding ramp endpoints, advances PF transport at the midpoint
opening stress, and localizes every exponential-threshold crossing inside that
same ramp interval.

No event batching, expected emission increment, constitutive refit, or physical
event clipping is introduced.
"""
from __future__ import annotations

import copy
import math
from typing import Any

import numpy as np

from .persistent_site_event_driven_emission_v10051831 import (
    EMISSION_HAZARD_SCHEMA,
    EMISSION_TRANSPORT_SCHEMA,
    PROBE_FALLBACK_SCHEMA,
    PersistentSiteEventDrivenSingleLayerLoadRampFrontEngineV10051831,
    two_channel_absolute_opening_drives_v1005183,
)
from .persistent_site_stochastic_emission_v100518 import (
    PF_UPDATE_MAP,
    _deposit_one_emission_event,
    _draw_exponential,
    _emission_hazards,
    _linear_hazard_crossing_time,
    _sum_numeric,
    advance_pf_transport_only_v10222,
)

MODEL_ID = "FEM_CZM_joint_K_ramp_event_driven_emission_v10_0_5_18_3_2"
LOAD_RAMP_SCHEMA = (
    "v10.0.5.18.3.2_outer_cleavage_cell_inner_joint_linear_K_event_horizon"
)
EVENT_DRIVEN_SCHEMA = (
    "v10.0.5.18.3.2_joint_K_state_event_horizon_individual_threshold_localization"
)
OUTER_EMISSION_LIMITER_SCHEMA = (
    "v10.0.5.18.3.2_no_outer_emission_action_or_rate_variation_limiter"
)


def _linear_value(start: float, end: float, alpha: float) -> float:
    a = min(max(float(alpha), 0.0), 1.0)
    return (1.0 - a) * float(start) + a * float(end)


def _active_log_hazard_change(initial: np.ndarray, final: np.ndarray) -> float:
    a = np.asarray(initial, dtype=float).reshape(-1)
    b = np.asarray(final, dtype=float).reshape(-1)
    active = np.maximum(a, b) > 1.0e-300
    if not np.any(active):
        return 0.0
    la = np.log(np.maximum(a[active], 1.0e-300))
    lb = np.log(np.maximum(b[active], 1.0e-300))
    return float(np.max(np.abs(lb - la)))


class PersistentSiteJointKRampEventDrivenFrontEngineV10051832(
    PersistentSiteEventDrivenSingleLayerLoadRampFrontEngineV10051831
):
    """Exact emission thresholds integrated with the accepted linear K ramp."""

    joint_K_ramp_inside_event_horizon_active = True
    outer_emission_rate_variation_limiter_active = False

    def reset(self):
        super().reset()
        self._joint_K_emit_context_active = False
        self._joint_K_emit_rate_Pa_sqrt_m_per_s = 0.0
        self.joint_K_ramp_half_steps_total = 0
        self.joint_K_ramp_events_total = 0

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
        """Limit the outer cell only by cleavage and moving-tip requirements."""

        remaining = max(float(remaining_s), 0.0)
        Kc0, _ = self._ramp_values(
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
            "outer_emission_action_limiter_active": 0.0,
            "outer_emission_rate_variation_limiter_active": 0.0,
        }

        for _ in range(16):
            h = min(max(float(h), min(remaining, self.min_substep_s)), remaining)
            Kc1, _ = self._ramp_values(
                consumed_s=consumed_s,
                total_s=total_s,
                offset_s=h,
                K_cleave_end=K_cleave_end,
                K_emit_end=K_emit_end,
            )
            lam1, _, _ = self.lambda_cleave(self.sigma_tip(Kc1), T_K)
            lam1 = max(float(lam1), 0.0) if math.isfinite(lam1) else 0.0
            progress1 = lam1 / threshold
            cleavage_action = 0.5 * (progress0 + progress1) * h
            active = max(progress0, progress1) > 1.0e-300
            cleavage_log = 0.0
            if active:
                cleavage_log = abs(
                    math.log(max(progress1, 1.0e-300))
                    - math.log(max(progress0, 1.0e-300))
                )
            diagnostics = {
                "cleavage_action": float(cleavage_action),
                "emission_action": 0.0,
                "cleavage_log_change": float(cleavage_log),
                "emission_log_change": 0.0,
                "outer_emission_action_limiter_active": 0.0,
                "outer_emission_rate_variation_limiter_active": 0.0,
            }

            ratio = cleavage_action / max(float(self.max_action_substep), 1.0e-12)
            if cleavage_action > self.ramp_action_floor:
                ratio = max(
                    ratio,
                    cleavage_log / max(self.ramp_max_log_rate_change, 1.0e-12),
                )
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
        dt = max(float(dt_s), 0.0)
        saved_active = bool(getattr(self, "_joint_K_emit_context_active", False))
        saved_rate = float(
            getattr(self, "_joint_K_emit_rate_Pa_sqrt_m_per_s", 0.0)
        )
        self._joint_K_emit_context_active = True
        self._joint_K_emit_rate_Pa_sqrt_m_per_s = (
            (float(K_emit) - float(self.load_ramp_K_emit_prev_Pa_sqrt_m)) / dt
            if dt > 0.0
            else 0.0
        )
        try:
            return super()._integrate_coupled(
                K_cleave=K_cleave,
                K_emit=K_emit,
                T_K=T_K,
                dt_s=dt,
                drive_factors=drive_factors,
                tau_signed_Pa=tau_signed_Pa,
            )
        finally:
            self._joint_K_emit_context_active = saved_active
            self._joint_K_emit_rate_Pa_sqrt_m_per_s = saved_rate

    def _plastic_half_step(
        self,
        *,
        dt_s: float,
        T_K: float,
        K_emit: float,
        drive_factors: np.ndarray,
        tau_signed_Pa: np.ndarray,
    ) -> tuple[dict[str, Any], float]:
        dt = max(float(dt_s), 0.0)
        if dt <= 0.0:
            return {
                "dN_emit": 0.0,
                "dN_trapped": 0.0,
                "dN_released": 0.0,
                "dN_escaped": 0.0,
                "dN_recovered": 0.0,
                "stochastic_emission_event_count": 0,
                "transport_integrator": PF_UPDATE_MAP,
                "joint_K_ramp_inside_event_horizon": True,
            }, 0.0

        rate = (
            float(self._joint_K_emit_rate_Pa_sqrt_m_per_s)
            if self._joint_K_emit_context_active
            else 0.0
        )
        K_mid = float(K_emit)
        K_start = K_mid - 0.5 * dt * rate
        K_end = K_mid + 0.5 * dt * rate
        result = self._stochastic_emission_transport_ramp(
            dt_s=dt,
            T_K=float(T_K),
            K_emit_start=K_start,
            K_emit_end=K_end,
            drive_factors=drive_factors,
            tau_signed_Pa=tau_signed_Pa,
        )
        work_stress_line_sum = float(
            result.get("emission_work_stress_line_sum_Pa", 0.0)
        )
        self.W_emit += (
            work_stress_line_sum
            * float(self.b)
            * float(self.f.L_pz)
        )
        self.joint_K_ramp_half_steps_total += 1
        self.joint_K_ramp_events_total += int(
            result.get("stochastic_emission_event_count", 0)
        )
        self._sync_compat()
        opening_end = self._opening_stress(K_end)[0]
        return result, opening_end

    def _stochastic_emission_transport_ramp(
        self,
        *,
        dt_s: float,
        T_K: float,
        K_emit_start: float,
        K_emit_end: float,
        drive_factors: np.ndarray,
        tau_signed_Pa: np.ndarray,
    ) -> dict[str, Any]:
        total = max(float(dt_s), 0.0)
        remaining = total
        diagnostics: dict[str, float] = {}
        events_this_step: list[dict[str, Any]] = []
        emitted_lines = 0.0
        work_stress_line_sum = 0.0
        proposals = 0
        refinements = 0
        elapsed = 0.0
        last_hazard = None
        max_log_change = 0.0
        max_proposed_action = 0.0

        def K_at(time_s: float) -> float:
            alpha = 1.0 if total <= 0.0 else float(time_s) / total
            return _linear_value(K_emit_start, K_emit_end, alpha)

        while remaining > 0.0:
            proposals += 1
            if proposals > 4 * self.emission_max_events_per_half_step + 2000:
                raise RuntimeError(
                    "joint-K event-driven stochastic emission exceeded the proposal "
                    "budget; the accepted local interval requires numerical refinement"
                )

            K0 = K_at(elapsed)
            opening0 = self._opening_stress(K0)[0]
            start = _emission_hazards(
                self.mpz_state,
                T_K=T_K,
                opening_stress_Pa=opening0,
                drive_factors=drive_factors,
                tau_signed_Pa=tau_signed_Pa,
                rate_function=self._emission_rate_per_site,
            )
            hazard0 = np.asarray(start["aggregate_hazard_s"], dtype=float)
            residual = np.maximum(
                self.emission_threshold_action - self.emission_action_current,
                0.0,
            )
            positive = hazard0 > 0.0
            predicted = np.divide(
                residual,
                hazard0,
                out=np.full_like(residual, np.inf),
                where=positive,
            )
            next_predicted = float(np.min(predicted))
            h = remaining
            if math.isfinite(next_predicted):
                h = min(
                    h,
                    max(
                        self.inner_event_horizon_factor * next_predicted,
                        min(remaining, 1.0e-15),
                    ),
                )
            h = min(max(float(h), min(remaining, 1.0e-15)), remaining)

            state_start = self.mpz_state.copy()
            proposal = None
            end = None
            hazard1 = None
            crossing_times = None
            K1 = K_at(elapsed + h)
            opening1 = self._opening_stress(K1)[0]

            for _ in range(self.inner_max_refinements):
                self.mpz_state = state_start.copy()
                K_mid = K_at(elapsed + 0.5 * h)
                opening_mid = self._opening_stress(K_mid)[0]
                proposal = advance_pf_transport_only_v10222(
                    self.mpz_state,
                    dt_s=h,
                    T_K=T_K,
                    opening_stress_Pa=opening_mid,
                )
                K1 = K_at(elapsed + h)
                opening1 = self._opening_stress(K1)[0]
                end = _emission_hazards(
                    self.mpz_state,
                    T_K=T_K,
                    opening_stress_Pa=opening1,
                    drive_factors=drive_factors,
                    tau_signed_Pa=tau_signed_Pa,
                    rate_function=self._emission_rate_per_site,
                )
                hazard1 = np.asarray(end["aggregate_hazard_s"], dtype=float)
                action_vector = 0.5 * (hazard0 + hazard1) * h
                proposed_action = float(np.max(action_vector))
                log_change = _active_log_hazard_change(hazard0, hazard1)
                max_log_change = max(max_log_change, log_change)
                max_proposed_action = max(max_proposed_action, proposed_action)

                crossing_times = np.full(self.mpz_state.n_systems, np.inf)
                for system in range(self.mpz_state.n_systems):
                    crossing = _linear_hazard_crossing_time(
                        hazard0[system],
                        hazard1[system],
                        h,
                        residual[system],
                    )
                    if crossing is not None:
                        crossing_times[system] = crossing

                needs_refinement = bool(
                    proposed_action > self.inner_action_floor
                    and log_change > self.inner_max_log_hazard_change
                    and h > 1.0e-15
                )
                if not needs_refinement:
                    break

                refinements += 1
                factor = min(
                    0.5,
                    max(
                        0.05,
                        0.8
                        * self.inner_max_log_hazard_change
                        / max(log_change, self.inner_max_log_hazard_change),
                    ),
                )
                h = max(min(h * factor, remaining), min(remaining, 1.0e-15))
            else:
                self.mpz_state = state_start
                raise RuntimeError(
                    "joint-K event-driven stochastic emission could not resolve "
                    "the state-and-load-dependent hazard variation"
                )

            if proposal is None or end is None or hazard1 is None or crossing_times is None:
                self.mpz_state = state_start
                raise RuntimeError(
                    "joint-K event-driven stochastic emission proposal was incomplete"
                )

            if not np.any(np.isfinite(crossing_times)):
                self.emission_action_current += 0.5 * (hazard0 + hazard1) * h
                _sum_numeric(diagnostics, proposal)
                elapsed += h
                remaining = max(remaining - h, 0.0)
                last_hazard = end
                continue

            system = int(np.argmin(crossing_times))
            tau = min(max(float(crossing_times[system]), 0.0), h)
            self.mpz_state = state_start
            K_partial_mid = K_at(elapsed + 0.5 * tau)
            opening_partial_mid = self._opening_stress(K_partial_mid)[0]
            partial = advance_pf_transport_only_v10222(
                self.mpz_state,
                dt_s=tau,
                T_K=T_K,
                opening_stress_Pa=opening_partial_mid,
            )
            K_event = K_at(elapsed + tau)
            opening_event = self._opening_stress(K_event)[0]
            at_event = _emission_hazards(
                self.mpz_state,
                T_K=T_K,
                opening_stress_Pa=opening_event,
                drive_factors=drive_factors,
                tau_signed_Pa=tau_signed_Pa,
                rate_function=self._emission_rate_per_site,
            )
            hazard_event = np.asarray(at_event["aggregate_hazard_s"], dtype=float)
            self.emission_action_current += 0.5 * (hazard0 + hazard_event) * tau
            self.emission_action_current[system] = self.emission_threshold_action[system]
            threshold_completed = float(self.emission_threshold_action[system])
            sign = float(np.asarray(at_event["signs"])[system])
            line_content = _deposit_one_emission_event(
                self.mpz_state,
                system,
                sign,
            )
            emitted_lines += line_content
            work_stress_line_sum += opening_event * line_content
            event_index = int(self.emission_event_index[system])
            event_record = {
                "system": system,
                "sign": sign,
                "event_index": event_index,
                "time_within_half_step_s": elapsed + tau,
                "K_emit_at_event_Pa_sqrt_m": float(K_event),
                "opening_stress_at_event_Pa": float(opening_event),
                "threshold_action": threshold_completed,
                "aggregate_hazard_at_event_s": float(hazard_event[system]),
                "rate_per_site_at_event_s": float(
                    np.asarray(at_event["rate_per_site_s"])[system]
                ),
                "multiplicity_at_event": float(at_event["multiplicity"]),
                "line_content_added": line_content,
                "tip_radius_before_event_m": float(
                    start["geometry"]["tip_radius_m"]
                ),
                "tip_radius_after_event_m": float(self.mpz_state.blunted_radius()),
                "front_width_at_event_m": float(
                    at_event["geometry"]["front_width_m"]
                ),
            }
            events_this_step.append(event_record)
            self.emission_event_history.append(copy.deepcopy(event_record))
            self.emission_event_index[system] += 1
            self.emission_event_count_total += 1
            self.emission_action_current[system] = 0.0
            self.emission_threshold_action[system] = _draw_exponential(
                self._emission_rngs[system],
                self.emission_minimum_threshold,
            )
            _sum_numeric(diagnostics, partial)
            elapsed += tau
            remaining = max(remaining - tau, 0.0)
            last_hazard = at_event

            if len(events_this_step) > self.emission_max_events_per_half_step:
                raise RuntimeError(
                    "stochastic emission exceeded EMISSION_MAX_EVENTS_PER_HALF_STEP; "
                    "the exact event density requires a smaller accepted local interval"
                )
            if tau <= 0.0 and remaining > 0.0:
                remaining = max(remaining - min(remaining, 1.0e-15), 0.0)

        self.emission_internal_substeps_total += proposals
        self.emission_transport_time_total_s += elapsed
        self.event_driven_proposals_total += proposals
        self.event_driven_refinements_total += refinements
        geometry = self.mpz_state.source_geometry()
        if last_hazard is None:
            K_final = K_at(total)
            opening_final = self._opening_stress(K_final)[0]
            last_hazard = _emission_hazards(
                self.mpz_state,
                T_K=T_K,
                opening_stress_Pa=opening_final,
                drive_factors=drive_factors,
                tau_signed_Pa=tau_signed_Pa,
                rate_function=self._emission_rate_per_site,
            )

        rho_back = np.asarray(last_hazard["rho_back_m2"], dtype=float).copy()
        sigma_back = np.asarray(last_hazard["sigma_back_Pa"], dtype=float).copy()
        result: dict[str, Any] = {
            **diagnostics,
            "dN_emit": emitted_lines,
            "source_activations": float(len(events_this_step)),
            "stochastic_emission_events": events_this_step,
            "stochastic_emission_event_count": len(events_this_step),
            "stochastic_emission_event_count_total": int(
                self.emission_event_count_total
            ),
            "stochastic_emission_threshold_action": self.emission_threshold_action.copy(),
            "stochastic_emission_action_current": self.emission_action_current.copy(),
            "stochastic_emission_event_index": self.emission_event_index.copy(),
            "stochastic_emission_stream_seeds": self.emission_stream_seeds.copy(),
            "aggregate_hazard_final_by_system_s": np.asarray(
                last_hazard["aggregate_hazard_s"], dtype=float
            ),
            "rate_final_by_system_s": np.asarray(
                last_hazard["rate_per_site_s"], dtype=float
            ),
            "persistent_site_geometry": geometry,
            "persistent_site_multiplicity_per_system": float(
                geometry["multiplicity_per_system"]
            ),
            "rho_back_final_by_system_m2": rho_back,
            "rho_back_final_max_m2": float(np.max(rho_back)),
            "sigma_back_final_by_system_Pa": sigma_back,
            "sigma_back_final_max_Pa": float(np.max(sigma_back)),
            "stochastic_emission_schema": EMISSION_HAZARD_SCHEMA,
            "stochastic_emission_transport_schema": EMISSION_TRANSPORT_SCHEMA,
            "transport_integrator": PF_UPDATE_MAP,
            "transport_drive_quadrature": "midpoint_opening_stress_on_each_adaptive_horizon",
            "hazard_drive_quadrature": "linear_K_endpoint_trapezoid",
            "joint_K_ramp_inside_event_horizon": True,
            "K_emit_start_Pa_sqrt_m": float(K_emit_start),
            "K_emit_end_Pa_sqrt_m": float(K_emit_end),
            "emission_work_stress_line_sum_Pa": float(work_stress_line_sum),
            "emission_internal_substeps": proposals,
            "event_driven_refinements": refinements,
            "event_driven_schema": EVENT_DRIVEN_SCHEMA,
            "event_driven_max_log_hazard_change": float(max_log_change),
            "event_driven_max_proposed_action": float(max_proposed_action),
            "emission_transport_time_s": elapsed,
            "finite_source_inventory_active": False,
            "source_depletion_active": False,
            "source_refresh_active": False,
            "explicit_recovery_active": False,
        }
        self.mpz_state.last_emission = copy.deepcopy(result)
        return result

    @classmethod
    def audit_payload(cls) -> dict[str, Any]:
        payload = super().audit_payload()
        payload["model_id"] = MODEL_ID
        emission = dict(payload.get("stochastic_emission", {}) or {})
        emission.update(
            {
                "load_ramp_schema": LOAD_RAMP_SCHEMA,
                "event_driven_transport_schema": EVENT_DRIVEN_SCHEMA,
                "outer_emission_action_limiter_active": False,
                "outer_emission_rate_variation_limiter_active": False,
                "K_ramp_evaluated_inside_event_horizon": True,
                "hazard_drive_quadrature": "linear_K_endpoint_trapezoid",
                "transport_drive_quadrature": (
                    "midpoint_opening_stress_on_each_adaptive_horizon"
                ),
                "threshold_crossings_localized_individually": True,
                "events_batched": False,
                "accepted_update": "exact_event_localized_one_activation_packets",
            }
        )
        payload["stochastic_emission"] = emission
        return payload


__all__ = [
    "EMISSION_HAZARD_SCHEMA",
    "EMISSION_TRANSPORT_SCHEMA",
    "EVENT_DRIVEN_SCHEMA",
    "LOAD_RAMP_SCHEMA",
    "MODEL_ID",
    "OUTER_EMISSION_LIMITER_SCHEMA",
    "PROBE_FALLBACK_SCHEMA",
    "PersistentSiteJointKRampEventDrivenFrontEngineV10051832",
    "two_channel_absolute_opening_drives_v1005183",
]
