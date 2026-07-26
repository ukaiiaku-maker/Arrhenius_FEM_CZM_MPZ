"""Single-layer local K-ramp integration for exact stochastic emission.

v10.0.5.18.3 correctly removed the deterministic global-FEM emission-state
predictor, but its outer load-ramp cell also limited the *expected emission
action* to ``EMISSION_MAX_ACTION_SUBSTEP``.  The inner exact stochastic
emission transport already performs that action control, localizes every
threshold crossing, deposits one activation packet, and recomputes backstress
and source geometry after each event.  Applying the same action limiter again
at the outer ramp level forced roughly twenty outer cells per stochastic event
and could exhaust ``max_internal_steps`` before the first accepted FEM step.

This point repair retains the outer K(t)-variation and cleavage/moving-tip
limits, while leaving emission-event action control exclusively in the exact
inner stochastic transport.  No constitutive parameter, event law, packet
normalization, source closure, or moving-tip coupling is changed.
"""
from __future__ import annotations

import math
from typing import Any

import numpy as np

from .persistent_site_load_ramp_stochastic_emission_v1005183 import (
    PROBE_FALLBACK_SCHEMA,
    PersistentSiteLoadRampStochasticEmissionFrontEngineV1005183,
    _maximum_log_change,
    two_channel_absolute_opening_drives_v1005183,
)
from .persistent_site_stochastic_emission_v100518 import (
    EMISSION_HAZARD_SCHEMA,
    EMISSION_TRANSPORT_SCHEMA,
    _emission_hazards,
)

MODEL_ID = "FEM_CZM_exact_stochastic_emission_single_layer_load_ramp_v10_0_5_18_3_1"
LOAD_RAMP_SCHEMA = (
    "v10.0.5.18.3.1_linear_K_ramp_single_layer_exact_stochastic_kinetics"
)
OUTER_EMISSION_LIMITER_SCHEMA = (
    "v10.0.5.18.3.1_outer_K_variation_inner_exact_event_action_control"
)


class PersistentSiteSingleLayerLoadRampStochasticEmissionFrontEngineV10051831(
    PersistentSiteLoadRampStochasticEmissionFrontEngineV1005183
):
    """Exact signed emission with one, nonduplicated local adaptivity layer."""

    outer_emission_action_limiter_active = False

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

        # Cleavage action, fractional MPZ translation, and the remaining
        # cleavage clock still constrain the outer moving-tip cell.
        h = self._substep_limit_stochastic(remaining, progress0)
        refinements = 0
        diagnostics = {
            "cleavage_action": 0.0,
            "emission_action": 0.0,
            "cleavage_log_change": 0.0,
            "emission_log_change": 0.0,
            "outer_emission_action_limiter_active": 0.0,
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
                "outer_emission_action_limiter_active": 0.0,
            }

            # The outer cell controls cleavage/moving-tip action.  Emission
            # action is diagnostic here because the inner exact stochastic
            # transport already limits integrated action and localizes each
            # threshold crossing.  Keeping it out of action_ratio removes the
            # duplicated twenty-cells-per-event refinement.
            action_ratio = cleavage_action / max(
                float(self.max_action_substep), 1.0e-12
            )

            # Retain K-ramp resolution where either active rate changes steeply.
            # The action floor avoids refining mathematically large log changes
            # when the corresponding integrated hazard is negligible.
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

    @classmethod
    def audit_payload(cls) -> dict[str, Any]:
        payload = super().audit_payload()
        payload["model_id"] = MODEL_ID
        emission = dict(payload.get("stochastic_emission", {}) or {})
        emission.update(
            {
                "load_ramp_schema": LOAD_RAMP_SCHEMA,
                "outer_emission_action_limiter_active": False,
                "outer_emission_limiter_schema": OUTER_EMISSION_LIMITER_SCHEMA,
                "inner_exact_emission_action_control": True,
                "global_FEM_state_change_predictor": False,
                "local_linear_K_ramp": True,
                "accepted_update": "exact_event_localized_one_activation_packets",
                "post_onset_mean_field_burst": False,
            }
        )
        payload["stochastic_emission"] = emission
        return payload


__all__ = [
    "EMISSION_HAZARD_SCHEMA",
    "EMISSION_TRANSPORT_SCHEMA",
    "LOAD_RAMP_SCHEMA",
    "MODEL_ID",
    "OUTER_EMISSION_LIMITER_SCHEMA",
    "PROBE_FALLBACK_SCHEMA",
    "PersistentSiteSingleLayerLoadRampStochasticEmissionFrontEngineV10051831",
    "two_channel_absolute_opening_drives_v1005183",
]
