"""Action-error controlled joint-K stochastic emission for v10.0.5.18.3.2.

The exact stochastic threshold and one-packet event laws are unchanged.  The
adaptive event horizon compares a one-interval transport/hazard estimate with
two sequential half-interval estimates.  The two-half state is accepted when
the integrated signed-channel hazard actions agree within a prescribed absolute
action tolerance.  Refinement is therefore driven by the quantity that advances
the stochastic clocks, rather than by a potentially huge logarithmic rate
change after backstress begins to evolve.
"""
from __future__ import annotations

import copy
import math
import os
from typing import Any

import numpy as np

from .persistent_site_joint_K_ramp_emission_v10051832 import (
    EMISSION_HAZARD_SCHEMA,
    EMISSION_TRANSPORT_SCHEMA,
    LOAD_RAMP_SCHEMA,
    OUTER_EMISSION_LIMITER_SCHEMA,
    PROBE_FALLBACK_SCHEMA,
    PersistentSiteJointKRampEventDrivenFrontEngineV10051832,
    _active_log_hazard_change,
    _linear_value,
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

MODEL_ID = "FEM_CZM_joint_K_action_error_event_driven_emission_v10_0_5_18_3_2"
EVENT_DRIVEN_SCHEMA = (
    "v10.0.5.18.3.2_joint_K_two_half_action_error_individual_threshold_localization"
)
ACTION_ERROR_SCHEMA = (
    "v10.0.5.18.3.2_full_vs_two_half_integrated_hazard_action_error"
)


def _sum_transport_diagnostics(*records: dict[str, Any]) -> dict[str, float]:
    out: dict[str, float] = {}
    for record in records:
        _sum_numeric(out, record)
    return out


class PersistentSiteJointKActionErrorFrontEngineV10051832(
    PersistentSiteJointKRampEventDrivenFrontEngineV10051832
):
    """Exact event localization with action-error controlled adaptive horizons."""

    action_error_control_active = True
    _inner_max_action_error_default = 0.01

    @classmethod
    def configure_stochastic(cls, **kwargs) -> None:
        super().configure_stochastic(**kwargs)
        cls._inner_max_action_error_default = max(
            float(os.environ.get("EMISSION_INNER_MAX_ACTION_ERROR", "0.01")),
            1.0e-8,
        )

    def reset(self):
        super().reset()
        self.inner_max_action_error = float(
            type(self)._inner_max_action_error_default
        )
        self.action_error_refinements_total = 0
        self.action_error_max_observed = 0.0

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
        max_action_error = 0.0

        def K_at(time_s: float) -> float:
            alpha = 1.0 if total <= 0.0 else float(time_s) / total
            return _linear_value(K_emit_start, K_emit_end, alpha)

        while remaining > 0.0:
            proposals += 1
            if proposals > 4 * self.emission_max_events_per_half_step + 2000:
                raise RuntimeError(
                    "joint-K action-error stochastic emission exceeded the proposal "
                    "budget; the accepted local interval requires refinement"
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
            accepted = False
            first_transport = None
            second_transport = None
            state_mid = None
            hazard_mid = None
            hazard_end = None
            action_first = None
            action_second = None
            action_two = None

            for _ in range(self.inner_max_refinements):
                half = 0.5 * h

                self.mpz_state = state_start.copy()
                K_full_mid = K_at(elapsed + 0.5 * h)
                opening_full_mid = self._opening_stress(K_full_mid)[0]
                advance_pf_transport_only_v10222(
                    self.mpz_state,
                    dt_s=h,
                    T_K=T_K,
                    opening_stress_Pa=opening_full_mid,
                )
                K1 = K_at(elapsed + h)
                opening1 = self._opening_stress(K1)[0]
                full_end = _emission_hazards(
                    self.mpz_state,
                    T_K=T_K,
                    opening_stress_Pa=opening1,
                    drive_factors=drive_factors,
                    tau_signed_Pa=tau_signed_Pa,
                    rate_function=self._emission_rate_per_site,
                )
                hazard_full_end = np.asarray(
                    full_end["aggregate_hazard_s"], dtype=float
                )
                action_full = 0.5 * (hazard0 + hazard_full_end) * h

                self.mpz_state = state_start.copy()
                K_first_mid = K_at(elapsed + 0.25 * h)
                opening_first_mid = self._opening_stress(K_first_mid)[0]
                first_transport = advance_pf_transport_only_v10222(
                    self.mpz_state,
                    dt_s=half,
                    T_K=T_K,
                    opening_stress_Pa=opening_first_mid,
                )
                Kmid = K_at(elapsed + half)
                opening_mid = self._opening_stress(Kmid)[0]
                mid = _emission_hazards(
                    self.mpz_state,
                    T_K=T_K,
                    opening_stress_Pa=opening_mid,
                    drive_factors=drive_factors,
                    tau_signed_Pa=tau_signed_Pa,
                    rate_function=self._emission_rate_per_site,
                )
                hazard_mid = np.asarray(mid["aggregate_hazard_s"], dtype=float)
                state_mid = self.mpz_state.copy()

                K_second_mid = K_at(elapsed + 0.75 * h)
                opening_second_mid = self._opening_stress(K_second_mid)[0]
                second_transport = advance_pf_transport_only_v10222(
                    self.mpz_state,
                    dt_s=half,
                    T_K=T_K,
                    opening_stress_Pa=opening_second_mid,
                )
                end = _emission_hazards(
                    self.mpz_state,
                    T_K=T_K,
                    opening_stress_Pa=opening1,
                    drive_factors=drive_factors,
                    tau_signed_Pa=tau_signed_Pa,
                    rate_function=self._emission_rate_per_site,
                )
                hazard_end = np.asarray(end["aggregate_hazard_s"], dtype=float)

                action_first = 0.5 * (hazard0 + hazard_mid) * half
                action_second = 0.5 * (hazard_mid + hazard_end) * half
                action_two = action_first + action_second
                action_error = float(np.max(np.abs(action_two - action_full)))
                max_action_error = max(max_action_error, action_error)
                max_log_change = max(
                    max_log_change,
                    _active_log_hazard_change(hazard0, hazard_mid),
                    _active_log_hazard_change(hazard_mid, hazard_end),
                )
                max_proposed_action = max(
                    max_proposed_action,
                    float(np.max(action_two)),
                )

                if (
                    action_error <= self.inner_max_action_error
                    or h <= 1.0e-15
                ):
                    last_hazard = end
                    accepted = True
                    break

                refinements += 1
                h = max(0.5 * h, min(remaining, 1.0e-15))

            if not accepted:
                self.mpz_state = state_start
                raise RuntimeError(
                    "joint-K stochastic emission could not resolve integrated "
                    "hazard action with the full-vs-two-half estimator"
                )

            if (
                first_transport is None
                or second_transport is None
                or state_mid is None
                or hazard_mid is None
                or hazard_end is None
                or action_first is None
                or action_second is None
                or action_two is None
            ):
                self.mpz_state = state_start
                raise RuntimeError(
                    "joint-K action-error stochastic emission proposal was incomplete"
                )

            half = 0.5 * h
            crossing_times = np.full(self.mpz_state.n_systems, np.inf)
            crossing_segments = np.full(self.mpz_state.n_systems, -1, dtype=int)
            for system in range(self.mpz_state.n_systems):
                first_crossing = _linear_hazard_crossing_time(
                    hazard0[system],
                    hazard_mid[system],
                    half,
                    residual[system],
                )
                if first_crossing is not None:
                    crossing_times[system] = first_crossing
                    crossing_segments[system] = 0
                    continue
                residual_second = max(
                    residual[system] - action_first[system],
                    0.0,
                )
                second_crossing = _linear_hazard_crossing_time(
                    hazard_mid[system],
                    hazard_end[system],
                    half,
                    residual_second,
                )
                if second_crossing is not None:
                    crossing_times[system] = half + second_crossing
                    crossing_segments[system] = 1

            if not np.any(np.isfinite(crossing_times)):
                self.emission_action_current += action_two
                committed = _sum_transport_diagnostics(
                    first_transport,
                    second_transport,
                )
                _sum_numeric(diagnostics, committed)
                elapsed += h
                remaining = max(remaining - h, 0.0)
                continue

            system = int(np.argmin(crossing_times))
            tau = min(max(float(crossing_times[system]), 0.0), h)
            segment = int(crossing_segments[system])

            if segment == 0:
                self.mpz_state = state_start.copy()
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
                hazard_event = np.asarray(
                    at_event["aggregate_hazard_s"], dtype=float
                )
                action_to_event = 0.5 * (hazard0 + hazard_event) * tau
                committed = partial
            else:
                self.mpz_state = state_mid.copy()
                tau_second = max(tau - half, 0.0)
                K_partial_mid = K_at(elapsed + half + 0.5 * tau_second)
                opening_partial_mid = self._opening_stress(K_partial_mid)[0]
                partial_second = advance_pf_transport_only_v10222(
                    self.mpz_state,
                    dt_s=tau_second,
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
                hazard_event = np.asarray(
                    at_event["aggregate_hazard_s"], dtype=float
                )
                action_to_event = action_first + 0.5 * (
                    hazard_mid + hazard_event
                ) * tau_second
                committed = _sum_transport_diagnostics(
                    first_transport,
                    partial_second,
                )

            self.emission_action_current += action_to_event
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
                "action_error_estimator": ACTION_ERROR_SCHEMA,
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
            _sum_numeric(diagnostics, committed)
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
        self.action_error_refinements_total += refinements
        self.action_error_max_observed = max(
            self.action_error_max_observed,
            max_action_error,
        )
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
            "transport_drive_quadrature": (
                "two_sequential_half_horizons_with_midpoint_opening_stress"
            ),
            "hazard_drive_quadrature": (
                "piecewise_linear_K_endpoint_trapezoid_two_half_horizons"
            ),
            "joint_K_ramp_inside_event_horizon": True,
            "action_error_control_active": True,
            "action_error_schema": ACTION_ERROR_SCHEMA,
            "action_error_tolerance": float(self.inner_max_action_error),
            "action_error_max_observed": float(max_action_error),
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
                "event_driven_transport_schema": EVENT_DRIVEN_SCHEMA,
                "state_dependent_hazard_error_control": "integrated_action",
                "action_error_control_active": True,
                "action_error_schema": ACTION_ERROR_SCHEMA,
                "action_error_default_tolerance": float(
                    cls._inner_max_action_error_default
                ),
                "log_rate_change_is_refinement_criterion": False,
                "threshold_crossings_localized_individually": True,
                "events_batched": False,
                "accepted_update": "exact_event_localized_one_activation_packets",
            }
        )
        payload["stochastic_emission"] = emission
        return payload


__all__ = [
    "ACTION_ERROR_SCHEMA",
    "EMISSION_HAZARD_SCHEMA",
    "EMISSION_TRANSPORT_SCHEMA",
    "EVENT_DRIVEN_SCHEMA",
    "LOAD_RAMP_SCHEMA",
    "MODEL_ID",
    "OUTER_EMISSION_LIMITER_SCHEMA",
    "PROBE_FALLBACK_SCHEMA",
    "PersistentSiteJointKActionErrorFrontEngineV10051832",
    "two_channel_absolute_opening_drives_v1005183",
]
