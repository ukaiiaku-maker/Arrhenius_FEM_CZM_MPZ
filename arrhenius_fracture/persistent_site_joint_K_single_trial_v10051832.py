"""Single-trial joint-K exact stochastic emission for v10.0.5.18.3.2.

The PF v10.2.22 transport map is retained as the authoritative discrete local
operator.  Each event horizon performs one transactional endpoint transport
trial; if a stochastic threshold is crossed, the threshold is localized and
the transport map is recomputed only to that event time before depositing one
activation packet.  Logarithmic hazard change is recorded diagnostically but
is not used to reject or recursively compose the validated PF transport map.
"""
from __future__ import annotations

import math
from typing import Any

from .persistent_site_joint_K_ramp_emission_v10051832 import (
    EMISSION_HAZARD_SCHEMA,
    EMISSION_TRANSPORT_SCHEMA,
    EVENT_DRIVEN_SCHEMA,
    LOAD_RAMP_SCHEMA,
    OUTER_EMISSION_LIMITER_SCHEMA,
    PROBE_FALLBACK_SCHEMA,
    PersistentSiteJointKRampEventDrivenFrontEngineV10051832,
    two_channel_absolute_opening_drives_v1005183,
)

MODEL_ID = "FEM_CZM_joint_K_single_trial_exact_emission_v10_0_5_18_3_2"
SINGLE_TRIAL_SCHEMA = (
    "v10.0.5.18.3.2_one_endpoint_trial_one_event_commit_PF_operator"
)


class PersistentSiteJointKSingleTrialFrontEngineV10051832(
    PersistentSiteJointKRampEventDrivenFrontEngineV10051832
):
    """Joint-K event localization without log-rate rejection of the PF map."""

    state_dependent_log_rate_rejection_active = False
    PF_transport_operator_preserved = True

    @classmethod
    def configure_stochastic(cls, **kwargs) -> None:
        super().configure_stochastic(**kwargs)
        cls._inner_max_log_hazard_change_default = math.inf

    def reset(self):
        super().reset()
        self.inner_max_log_hazard_change = math.inf
        self.state_dependent_log_rate_rejection_active = False

    @classmethod
    def audit_payload(cls) -> dict[str, Any]:
        payload = super().audit_payload()
        payload["model_id"] = MODEL_ID
        emission = dict(payload.get("stochastic_emission", {}) or {})
        emission.update(
            {
                "single_trial_event_horizon_active": True,
                "single_trial_event_horizon_schema": SINGLE_TRIAL_SCHEMA,
                "state_dependent_log_rate_rejection_active": False,
                "log_rate_change_is_diagnostic_only": True,
                "PF_transport_operator_preserved": True,
                "PF_transport_recursive_composition_error_control": False,
                "one_endpoint_transport_trial_per_event_horizon": True,
                "one_committed_transport_update_to_localized_event": True,
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
    "SINGLE_TRIAL_SCHEMA",
    "PersistentSiteJointKSingleTrialFrontEngineV10051832",
    "two_channel_absolute_opening_drives_v1005183",
]
