"""Event-driven exact signed-emission transport for v10.0.5.18.3.1.

The stochastic threshold law is unchanged.  Instead of advancing a fixed 0.05
hazard action before every event, the integrator proposes a horizon slightly
beyond the next threshold predicted from the current aggregate hazard.  The PF
transport map is evaluated transactionally over that horizon, the hazard
variation is checked, and the horizon is refined only when the state-dependent
hazard changes too strongly.  A threshold crossing is then localized with the
same linearly integrated hazard solve used by v10.0.5.18.

This removes the approximately twenty transport proposals per event without
batching events or replacing them by an expected increment.
"""
from __future__ import annotations

import copy
import math
import os
from typing import Any

import numpy as np

from .persistent_site_load_ramp_stochastic_emission_v10051831 import (
    LOAD_RAMP_SCHEMA,
    OUTER_EMISSION_LIMITER_SCHEMA,
    PROBE_FALLBACK_SCHEMA,
    PersistentSiteSingleLayerLoadRampStochasticEmissionFrontEngineV10051831,
    two_channel_absolute_opening_drives_v1005183,
)
from .persistent_site_stochastic_emission_v100518 import (
    EMISSION_HAZARD_SCHEMA,
    EMISSION_TRANSPORT_SCHEMA,
    PF_UPDATE_MAP,
    _deposit_one_emission_event,
    _draw_exponential,
    _emission_hazards,
    _linear_hazard_crossing_time,
    _sum_numeric,
    advance_pf_transport_only_v10222,
)

MODEL_ID = "FEM_CZM_event_driven_single_layer_load_ramp_v10_0_5_18_3_1"
EVENT_DRIVEN_SCHEMA = (
    "v10.0.5.18.3.1_next_threshold_event_horizon_adaptive_hazard_variation"
)


def _active_log_hazard_change(initial: np.ndarray, final: np.ndarray) -> float:
    a = np.asarray(initial, dtype=float).reshape(-1)
    b = np.asarray(final, dtype=float).reshape(-1)
    active = np.maximum(a, b) > 1.0e-300
    if not np.any(active):
        return 0.0
    la = np.log(np.maximum(a[active], 1.0e-300))
    lb = np.log(np.maximum(b[active], 1.0e-300))
    return float(np.max(np.abs(lb - la)))


class PersistentSiteEventDrivenSingleLayerLoadRampFrontEngineV10051831(
    PersistentSiteSingleLayerLoadRampStochasticEmissionFrontEngineV10051831
):
    """One exact event per localized threshold, with adaptive event horizons."""

    event_driven_emission_transport_active = True
    _inner_max_log_hazard_change_default = 0.5
    _inner_action_floor_default = 1.0e-6
    _inner_event_horizon_factor_default = 1.25
    _inner_max_refinements_default = 24

    @classmethod
    def configure_stochastic(cls, **kwargs) -> None:
        super().configure_stochastic(**kwargs)
        cls._inner_max_log_hazard_change_default = max(
            float(os.environ.get("EMISSION_INNER_MAX_LOG_HAZARD_CHANGE", "0.5")),
            1.0e-6,
        )
        cls._inner_action_floor_default = max(
            float(os.environ.get("EMISSION_INNER_ACTION_FLOOR", "1e-6")),
            1.0e-15,
        )
        cls._inner_event_horizon_factor_default = max(
            float(os.environ.get("EMISSION_EVENT_HORIZON_FACTOR", "1.25")),
            1.0,
        )
        cls._inner_max_refinements_default = max(
            int(os.environ.get("EMISSION_INNER_MAX_REFINEMENTS", "24")),
            1,
        )

    def reset(self):
        super().reset()
        self.inner_max_log_hazard_change = float(
            type(self)._inner_max_log_hazard_change_default
        )
        self.inner_action_floor = float(type(self)._inner_action_floor_default)
        self.inner_event_horizon_factor = float(
            type(self)._inner_event_horizon_factor_default
        )
        self.inner_max_refinements = int(type(self)._inner_max_refinements_default)
        self.event_driven_proposals_total = 0
        self.event_driven_refinements_total = 0

    def _stochastic_emission_transport(
        self,
        *,
        dt_s: float,
        T_K: float,
        opening_stress_Pa: float,
        drive_factors: np.ndarray,
        tau_signed_Pa: np.ndarray,
    ) -> dict[str, Any]:
        remaining = max(float(dt_s), 0.0)
        diagnostics: dict[str, float] = {}
        events_this_step: list[dict[str, Any]] = []
        emitted_lines = 0.0
        proposals = 0
        refinements = 0
        elapsed = 0.0
        last_hazard = None
        max_log_change = 0.0
        max_proposed_action = 0.0

        while remaining > 0.0:
            proposals += 1
            if proposals > 4 * self.emission_max_events_per_half_step + 2000:
                raise RuntimeError(
                    "event-driven stochastic emission exceeded the proposal budget; "
                    "the accepted local interval requires numerical refinement"
                )

            start = _emission_hazards(
                self.mpz_state,
                T_K=T_K,
                opening_stress_Pa=opening_stress_Pa,
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
            if not np.any(positive):
                transport = advance_pf_transport_only_v10222(
                    self.mpz_state,
                    dt_s=remaining,
                    T_K=T_K,
                    opening_stress_Pa=opening_stress_Pa,
                )
                _sum_numeric(diagnostics, transport)
                elapsed += remaining
                remaining = 0.0
                last_hazard = start
                break

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

            for _ in range(self.inner_max_refinements):
                self.mpz_state = state_start.copy()
                proposal = advance_pf_transport_only_v10222(
                    self.mpz_state,
                    dt_s=h,
                    T_K=T_K,
                    opening_stress_Pa=opening_stress_Pa,
                )
                end = _emission_hazards(
                    self.mpz_state,
                    T_K=T_K,
                    opening_stress_Pa=opening_stress_Pa,
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
                    "event-driven stochastic emission could not resolve the "
                    "state-dependent hazard variation"
                )

            if proposal is None or end is None or hazard1 is None or crossing_times is None:
                self.mpz_state = state_start
                raise RuntimeError("event-driven stochastic emission proposal was incomplete")

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
            partial = advance_pf_transport_only_v10222(
                self.mpz_state,
                dt_s=tau,
                T_K=T_K,
                opening_stress_Pa=opening_stress_Pa,
            )
            at_event = _emission_hazards(
                self.mpz_state,
                T_K=T_K,
                opening_stress_Pa=opening_stress_Pa,
                drive_factors=drive_factors,
                tau_signed_Pa=tau_signed_Pa,
                rate_function=self._emission_rate_per_site,
            )
            hazard_event = np.asarray(
                at_event["aggregate_hazard_s"], dtype=float
            )
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
            event_index = int(self.emission_event_index[system])
            event_record = {
                "system": system,
                "sign": sign,
                "event_index": event_index,
                "time_within_half_step_s": elapsed + tau,
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
            last_hazard = _emission_hazards(
                self.mpz_state,
                T_K=T_K,
                opening_stress_Pa=opening_stress_Pa,
                drive_factors=drive_factors,
                tau_signed_Pa=tau_signed_Pa,
                rate_function=self._emission_rate_per_site,
            )

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
            "multiplicity_final": float(last_hazard["multiplicity"]),
            "rho_back_final_m2": float(last_hazard["rho_back_m2"]),
            "sigma_back_final_Pa": float(last_hazard["sigma_back_Pa"]),
            "source_area_final_m2": float(geometry["source_area_m2"]),
            "front_width_final_m": float(geometry["front_width_m"]),
            "tip_radius_final_m": float(geometry["tip_radius_m"]),
            "transport_integrator": PF_UPDATE_MAP,
            "emission_internal_substeps": proposals,
            "event_driven_refinements": refinements,
            "event_driven_schema": EVENT_DRIVEN_SCHEMA,
            "event_driven_max_log_hazard_change": float(max_log_change),
            "event_driven_max_proposed_action": float(max_proposed_action),
            "emission_transport_time_s": elapsed,
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
                "event_driven_transport_active": True,
                "event_driven_transport_schema": EVENT_DRIVEN_SCHEMA,
                "fixed_inner_action_substep_per_event": False,
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
    "PersistentSiteEventDrivenSingleLayerLoadRampFrontEngineV10051831",
    "two_channel_absolute_opening_drives_v1005183",
]
