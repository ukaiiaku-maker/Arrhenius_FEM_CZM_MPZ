"""Action-weighted refinement of the v10.0.5.18.2 emission predictor."""
from __future__ import annotations

from typing import Any

import numpy as np

from .persistent_site_state_coupled_stochastic_emission_v1005182 import (
    ADAPTIVE_SCHEMA,
    PersistentSiteStateCoupledStochasticEmissionFrontEngineV1005182,
)

MODEL_ID = "FEM_CZM_action_weighted_state_coupled_stochastic_emission_v10_0_5_18_2"
ACTION_WEIGHTING_SCHEMA = (
    "v10.0.5.18.2_log_hazard_change_weighted_by_expected_action"
)


class PersistentSiteActionWeightedStochasticEmissionFrontEngineV1005182(
    PersistentSiteStateCoupledStochasticEmissionFrontEngineV1005182
):
    """Ignore large relative rate changes when their integrated action is negligible."""

    def _state_change_prediction(
        self,
        *,
        K_emit: float,
        T_K: float,
        dt_s: float,
        drive_factors: np.ndarray,
        tau_signed_Pa: np.ndarray,
    ) -> dict[str, Any]:
        prediction = super()._state_change_prediction(
            K_emit=K_emit,
            T_K=T_K,
            dt_s=dt_s,
            drive_factors=drive_factors,
            tau_signed_Pa=tau_signed_Pa,
        )
        dt = max(float(dt_s), 0.0)
        start = np.asarray(
            prediction["aggregate_hazard_start_by_system_s"], dtype=float
        )
        endpoint = np.asarray(
            prediction["aggregate_hazard_endpoint_by_system_s"], dtype=float
        )
        feedback = np.asarray(
            prediction["aggregate_hazard_feedback_by_system_s"], dtype=float
        )
        load_action = float(np.max(0.5 * (start + endpoint) * dt))
        feedback_action = float(np.max(0.5 * (endpoint + feedback) * dt))
        reference = max(float(self.emission_adaptive_reference_target), 1.0e-12)
        load_weight = min(load_action / reference, 1.0)
        feedback_weight = min(feedback_action / reference, 1.0)
        raw_load = float(prediction["load_log_hazard_change"])
        raw_feedback = float(prediction["feedback_log_hazard_change"])
        weighted_load = raw_load * load_weight
        weighted_feedback = raw_feedback * feedback_weight
        weighted_log_change = max(weighted_load, weighted_feedback)

        ratios = dict(prediction["component_ratios"])
        ratios["log_hazard"] = (
            weighted_log_change / self.emission_adaptive_max_log_hazard_change
        )
        controlling = max(ratios, key=ratios.get)
        normalized = self.emission_adaptive_reference_target * max(
            max(ratios.values()), 0.0
        )
        prediction.update(
            {
                "adaptive_increment": float(normalized),
                "controlling_component": controlling,
                "component_ratios": ratios,
                "load_expected_action": load_action,
                "feedback_expected_action": feedback_action,
                "load_log_hazard_action_weight": load_weight,
                "feedback_log_hazard_action_weight": feedback_weight,
                "load_log_hazard_change_raw": raw_load,
                "feedback_log_hazard_change_raw": raw_feedback,
                "load_log_hazard_change": weighted_load,
                "feedback_log_hazard_change": weighted_feedback,
                "action_weighted_log_hazard_change": weighted_log_change,
                "action_weighting_schema": ACTION_WEIGHTING_SCHEMA,
            }
        )
        return prediction

    @classmethod
    def audit_payload(cls) -> dict[str, Any]:
        payload = super().audit_payload()
        payload["model_id"] = MODEL_ID
        emission = dict(payload.get("stochastic_emission", {}) or {})
        emission.update(
            {
                "adaptive_schema": ADAPTIVE_SCHEMA,
                "action_weighting_schema": ACTION_WEIGHTING_SCHEMA,
                "relative_hazard_change_weighted_by_expected_action": True,
                "negligible_hazard_does_not_force_load_subdivision": True,
                "backstress_and_geometry_limits_remain_unweighted": True,
            }
        )
        payload["stochastic_emission"] = emission
        return payload


__all__ = [
    "ACTION_WEIGHTING_SCHEMA",
    "MODEL_ID",
    "PersistentSiteActionWeightedStochasticEmissionFrontEngineV1005182",
]
