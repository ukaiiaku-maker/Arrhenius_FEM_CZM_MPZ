"""Hybrid stochastic-onset/mean-field emission bursts for FEM/CZM v10.0.5.18.1.

The exact microscopic signed-emission SSA introduced in v10.0.5.18 becomes
stiff when hundreds of persistent sites operate near their attempt frequency:
tens of thousands of individual activations can occur before Taylor
backstress arrests the episode.  The physically meaningful stochastic
quantity requested for production is the onset of each signed emission
episode.  This module therefore retains one exact exponential integrated-
hazard onset clock per signed channel, deposits the triggering activation at
the localized onset time, and resolves the subsequent high-rate burst with
the pre-existing backstress-limited mean-field activation law.  Transport,
blunting, source-area evolution, signed populations, and trial rollback are
unchanged.
"""
from __future__ import annotations

import copy
import math
from typing import Any

import numpy as np

from .persistent_site_pf_update_v100515 import PF_REFERENCE_COMMIT, PF_UPDATE_MAP
from .persistent_site_stochastic_emission_v100518 import (
    EMISSION_HAZARD_SCHEMA as MICROSCOPIC_EMISSION_HAZARD_SCHEMA,
    PersistentSiteStochasticEmissionMovingTipFrontEngineV100518,
    _deposit_one_emission_event,
    _draw_exponential,
    _emission_hazards,
    _linear_hazard_crossing_time,
    _sum_numeric,
    advance_pf_transport_only_v10222,
)

MODEL_ID = "FEM_CZM_stochastic_emission_onset_mean_field_burst_v10_0_5_18_1"
EMISSION_HAZARD_SCHEMA = (
    "v10.0.5.18.1_independent_signed_channel_exponential_onset_hazard"
)
EMISSION_TRANSPORT_SCHEMA = (
    "v10.0.5.18.1_exact_onset_then_backstress_limited_burst_"
    "interleaved_PF_v10_2_22_transport"
)


def _masked_channel(values: np.ndarray, channel: int) -> np.ndarray:
    out = np.zeros_like(np.asarray(values, dtype=float).reshape(-1))
    out[int(channel)] = float(np.asarray(values, dtype=float).reshape(-1)[int(channel)])
    return out


class PersistentSiteStochasticOnsetBurstMovingTipFrontEngineV1005181(
    PersistentSiteStochasticEmissionMovingTipFrontEngineV100518
):
    """Exact stochastic episode onset followed by mean-field burst relaxation."""

    stochastic_emission_active = True
    stochastic_emission_onset_exact = True
    post_onset_burst_mean_field = True

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
        elapsed = 0.0
        diagnostics: dict[str, float] = {}
        onset_events: list[dict[str, Any]] = []
        n_systems = int(self.mpz_state.n_systems)
        resolved = np.zeros(n_systems, dtype=bool)
        line_content_by_system = np.zeros(n_systems, dtype=float)
        trigger_activations_by_system = np.zeros(n_systems, dtype=float)
        burst_activations_by_system = np.zeros(n_systems, dtype=float)
        first_hazard = None
        substeps = 0

        while remaining > 0.0:
            substeps += 1
            if substeps > 10000:
                raise RuntimeError(
                    "stochastic-onset emission exceeded the internal transport budget; "
                    "reduce the outer timestep"
                )

            start = _emission_hazards(
                self.mpz_state,
                T_K=T_K,
                opening_stress_Pa=opening_stress_Pa,
                drive_factors=drive_factors,
                tau_signed_Pa=tau_signed_Pa,
                rate_function=self._emission_rate_per_site,
            )
            if first_hazard is None:
                first_hazard = copy.deepcopy(start)
            hazard0 = np.asarray(start["aggregate_hazard_s"], dtype=float).copy()
            hazard0[resolved] = 0.0
            active = hazard0 > 0.0

            if not np.any(active):
                transport = advance_pf_transport_only_v10222(
                    self.mpz_state,
                    dt_s=remaining,
                    T_K=T_K,
                    opening_stress_Pa=opening_stress_Pa,
                )
                _sum_numeric(diagnostics, transport)
                elapsed += remaining
                remaining = 0.0
                break

            residual = np.maximum(
                self.emission_threshold_action - self.emission_action_current,
                0.0,
            )
            max_hazard = float(np.max(hazard0))
            h = min(
                remaining,
                self.emission_max_action_substep / max(max_hazard, 1.0e-300),
            )
            predicted = np.divide(
                residual,
                hazard0,
                out=np.full_like(residual, np.inf),
                where=active,
            )
            next_predicted = float(np.min(predicted))
            if math.isfinite(next_predicted):
                h = min(h, max(1.25 * next_predicted, 1.0e-15))
            h = min(max(h, min(remaining, 1.0e-15)), remaining)

            state_start = self.mpz_state.copy()
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
            hazard1 = np.asarray(end["aggregate_hazard_s"], dtype=float).copy()
            hazard1[resolved] = 0.0
            crossing_times = np.full(n_systems, np.inf)
            for system in range(n_systems):
                if resolved[system]:
                    continue
                crossing = _linear_hazard_crossing_time(
                    hazard0[system],
                    hazard1[system],
                    h,
                    residual[system],
                )
                if crossing is not None:
                    crossing_times[system] = crossing

            if not np.any(np.isfinite(crossing_times)):
                increment = 0.5 * (hazard0 + hazard1) * h
                self.emission_action_current[~resolved] += increment[~resolved]
                _sum_numeric(diagnostics, proposal)
                elapsed += h
                remaining = max(remaining - h, 0.0)
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
            increment = 0.5 * (hazard0 + hazard_event) * tau
            self.emission_action_current[~resolved] += increment[~resolved]
            self.emission_action_current[system] = self.emission_threshold_action[
                system
            ]

            threshold_completed = float(self.emission_threshold_action[system])
            sign = float(np.asarray(at_event["signs"], dtype=float)[system])
            radius_before = float(self.mpz_state.blunted_radius())
            trigger_line_content = _deposit_one_emission_event(
                self.mpz_state,
                system,
                sign,
            )
            trigger_activations_by_system[system] += 1.0
            line_content_by_system[system] += trigger_line_content

            remaining_after_onset = max(remaining - tau, 0.0)
            burst = {
                "dN_emit": 0.0,
                "source_activations": 0.0,
                "activations_by_system": np.zeros(n_systems),
            }
            if remaining_after_onset > 0.0:
                burst = self.mpz_state.emit_persistent(
                    dt_s=remaining_after_onset,
                    T_K=float(T_K),
                    opening_stress_Pa=float(opening_stress_Pa),
                    drive_factors=_masked_channel(drive_factors, system),
                    tau_signed_Pa=_masked_channel(tau_signed_Pa, system),
                    rate_function=self._emission_rate_per_site,
                )
                burst_activations = np.asarray(
                    burst.get("activations_by_system", np.zeros(n_systems)),
                    dtype=float,
                )
                burst_lines = np.asarray(
                    burst.get("line_content_by_system", np.zeros(n_systems)),
                    dtype=float,
                )
                burst_activations_by_system += burst_activations
                line_content_by_system += burst_lines

            event_index = int(self.emission_event_index[system])
            event_record = {
                "system": system,
                "sign": sign,
                "event_index": event_index,
                "time_within_half_step_s": elapsed + tau,
                "threshold_action": threshold_completed,
                "aggregate_onset_hazard_s": float(hazard_event[system]),
                "rate_per_site_at_onset_s": float(
                    np.asarray(at_event["rate_per_site_s"], dtype=float)[system]
                ),
                "multiplicity_at_onset": float(at_event["multiplicity"]),
                "trigger_line_content_added": trigger_line_content,
                "post_onset_burst_activations": float(
                    burst_activations_by_system[system]
                ),
                "post_onset_burst_line_content": float(
                    line_content_by_system[system] - trigger_line_content
                ),
                "tip_radius_before_onset_m": radius_before,
                "tip_radius_after_burst_m": float(
                    self.mpz_state.blunted_radius()
                ),
                "front_width_at_onset_m": float(
                    at_event["geometry"]["front_width_m"]
                ),
            }
            onset_events.append(event_record)
            self.emission_event_history.append(copy.deepcopy(event_record))
            self.emission_event_index[system] += 1
            self.emission_event_count_total += 1
            self.emission_action_current[system] = 0.0
            self.emission_threshold_action[system] = _draw_exponential(
                self._emission_rngs[system],
                self.emission_minimum_threshold,
            )
            resolved[system] = True

            _sum_numeric(diagnostics, partial)
            _sum_numeric(diagnostics, burst)
            elapsed += tau
            remaining = remaining_after_onset

        self.emission_internal_substeps_total += substeps
        self.emission_transport_time_total_s += elapsed
        final_hazard = _emission_hazards(
            self.mpz_state,
            T_K=T_K,
            opening_stress_Pa=opening_stress_Pa,
            drive_factors=drive_factors,
            tau_signed_Pa=tau_signed_Pa,
            rate_function=self._emission_rate_per_site,
        )
        geometry = self.mpz_state.source_geometry()
        total_activations_by_system = (
            trigger_activations_by_system + burst_activations_by_system
        )
        emitted_lines = float(np.sum(line_content_by_system))
        result: dict[str, Any] = {
            **diagnostics,
            "dN_emit": emitted_lines,
            "source_activations": float(np.sum(total_activations_by_system)),
            "activations_by_system": total_activations_by_system,
            "trigger_activations_by_system": trigger_activations_by_system,
            "post_onset_burst_activations_by_system": (
                burst_activations_by_system
            ),
            "line_content_by_system": line_content_by_system,
            "stochastic_emission_events": onset_events,
            "stochastic_emission_event_count": len(onset_events),
            "stochastic_emission_event_count_total": int(
                self.emission_event_count_total
            ),
            "stochastic_emission_threshold_action": (
                self.emission_threshold_action.copy()
            ),
            "stochastic_emission_action_current": (
                self.emission_action_current.copy()
            ),
            "stochastic_emission_event_index": self.emission_event_index.copy(),
            "stochastic_emission_stream_seeds": self.emission_stream_seeds.copy(),
            "aggregate_hazard_initial_by_system_s": np.asarray(
                (first_hazard or final_hazard)["aggregate_hazard_s"], dtype=float
            ),
            "aggregate_hazard_final_by_system_s": np.asarray(
                final_hazard["aggregate_hazard_s"], dtype=float
            ),
            "rate_final_by_system_s": np.asarray(
                final_hazard["rate_per_site_s"], dtype=float
            ),
            "persistent_site_geometry": geometry,
            "persistent_site_multiplicity_per_system": float(
                geometry["multiplicity_per_system"]
            ),
            "stochastic_emission_schema": EMISSION_HAZARD_SCHEMA,
            "stochastic_emission_transport_schema": EMISSION_TRANSPORT_SCHEMA,
            "stochastic_onset_exact": True,
            "post_onset_burst_mean_field": True,
            "individual_post_onset_activations_explicitly_sampled": False,
            "microscopic_exact_SSA_replaced_due_to_stiffness": True,
            "microscopic_predecessor_schema": MICROSCOPIC_EMISSION_HAZARD_SCHEMA,
            "emission_internal_substeps": substeps,
            "emission_transport_time_s": elapsed,
            "finite_source_inventory_active": False,
            "source_depletion_active": False,
            "source_refresh_active": False,
            "explicit_recovery_active": False,
            "transport_integrator": PF_UPDATE_MAP,
            "pf_reference_commit": PF_REFERENCE_COMMIT,
        }
        self.mpz_state.last_emission = copy.deepcopy(result)
        return result

    @classmethod
    def audit_payload(cls) -> dict[str, Any]:
        payload = super().audit_payload()
        payload["model_id"] = MODEL_ID
        payload["stochastic_emission"] = {
            "schema": EMISSION_HAZARD_SCHEMA,
            "transport_schema": EMISSION_TRANSPORT_SCHEMA,
            "case_seed": int(cls._emission_case_seed_default),
            "stream_seeds": [
                int(seed)
                for seed in getattr(cls, "_audit_stream_seeds", [])
            ],
            "independent_signed_channel_streams": True,
            "onset_distribution": "exponential_unit_mean_integrated_action",
            "stochastic_onset_exact": True,
            "post_onset_burst_mean_field": True,
            "post_onset_burst_law": (
                "existing_backstress_limited_persistent_site_activation_solver"
            ),
            "individual_post_onset_activations_explicitly_sampled": False,
            "reason_for_hybridization": (
                "stiff high-rate episodes contain many microscopic activations "
                "before mechanical blocking"
            ),
            "transport_interleaved_with_onset_localization": True,
            "trial_rollback_includes_emission_rng_and_clock_state": True,
            "persistent_sites": True,
            "finite_source_inventory": False,
            "source_refresh": False,
            "explicit_recovery": False,
        }
        return payload


__all__ = [
    "MODEL_ID",
    "EMISSION_HAZARD_SCHEMA",
    "EMISSION_TRANSPORT_SCHEMA",
    "PersistentSiteStochasticOnsetBurstMovingTipFrontEngineV1005181",
]
