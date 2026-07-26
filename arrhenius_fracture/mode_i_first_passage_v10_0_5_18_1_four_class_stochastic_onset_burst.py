"""v10.0.5.18.1 four-class FEM/CZM with stochastic emission onset.

This point release keeps the v10.0.5.18 parameter transfer, stochastic
cleavage threshold, threshold-correlated crack distance, moving-MPZ advection,
and seed domains.  It replaces the stiff exact microscopic emission SSA with
exact stochastic onset clocks followed by the existing backstress-limited
mean-field burst relaxation.
"""
from __future__ import annotations

import copy
import numpy as np

from . import mode_i_first_passage_v10_0_5_18_four_class_stochastic_emission as _base
from .persistent_site_stochastic_emission_v100518 import derive_stream_seed
from .persistent_site_stochastic_onset_burst_v1005181 import (
    EMISSION_HAZARD_SCHEMA,
    EMISSION_TRANSPORT_SCHEMA,
    PersistentSiteStochasticOnsetBurstMovingTipFrontEngineV1005181,
)

POINT_RELEASE = "10.0.5.18.1"
MODEL_ID = "FEM_CZM_four_class_stochastic_cleavage_emission_onset_burst_v10_0_5_18_1"
PRODUCTION_MANIFEST = "persistent_site_production_manifest_v10_0_5_18_1.json"
SELECTION_MANIFEST = "persistent_site_parameter_selection_v10_0_5_18_1.json"
TRANSFER_MANIFEST = "four_class_parameter_transfer_v10_0_5_18_1.json"
SEED_MANIFEST = "stochastic_seed_manifest_v10_0_5_18_1.json"
DEFAULT_PARAMETER_SOURCE_ROOT = _base.DEFAULT_PARAMETER_SOURCE_ROOT


class ProductionStochasticOnsetBurstFrontEngineV1005181(
    PersistentSiteStochasticOnsetBurstMovingTipFrontEngineV1005181
):
    """Production diagnostics for exact-onset, mean-field-burst emission."""

    def step_drives(self, K_cleave, K_emit, T, dt, metadata=None):
        out = super().step_drives(
            K_cleave,
            K_emit,
            T,
            dt,
            metadata=metadata,
        )
        last = dict(self.mpz_state.last_emission or {})
        aggregate = np.asarray(
            last.get(
                "aggregate_hazard_final_by_system_s",
                last.get("aggregate_hazard_initial_by_system_s", np.zeros(2)),
            ),
            dtype=float,
        ).reshape(-1)
        out.update(
            {
                "lambda_e": float(np.sum(aggregate)),
                "lambda_e_semantics": (
                    "instantaneous_sum_signed_channel_emission_onset_hazards"
                ),
                "stochastic_emission_schema": EMISSION_HAZARD_SCHEMA,
                "stochastic_emission_transport_schema": EMISSION_TRANSPORT_SCHEMA,
                "emission_event_packet_semantics": (
                    "one_exact_trigger_activation_then_backstress_limited_"
                    "mean_field_burst"
                ),
                "stochastic_emission_onset_exact": True,
                "post_onset_burst_mean_field": True,
                "individual_post_onset_activations_explicitly_sampled": False,
            }
        )
        return out

    @classmethod
    def audit_payload(cls):
        payload = copy.deepcopy(super().audit_payload())
        case_seed = int(cls._emission_case_seed_default)
        payload["model_id"] = MODEL_ID
        payload["stochastic_emission"] = {
            "schema": EMISSION_HAZARD_SCHEMA,
            "transport_schema": EMISSION_TRANSPORT_SCHEMA,
            "case_seed": case_seed,
            "stream_seeds": [
                derive_stream_seed(case_seed, "signed_emission", channel)
                for channel in range(2)
            ],
            "independent_signed_channel_streams": True,
            "onset_distribution": "exponential_unit_mean_integrated_action",
            "stochastic_onset_exact": True,
            "trigger_packet": "one_activation_times_mechanical_line_conversion",
            "post_onset_burst_mean_field": True,
            "post_onset_burst_law": (
                "existing_backstress_limited_persistent_site_activation_solver"
            ),
            "individual_post_onset_activations_explicitly_sampled": False,
            "stiffness_resolution_changes_mean_constitutive_law": False,
            "transport_interleaved_with_onset_localization": True,
            "trial_rejection_restores_rng_and_clock_state": True,
            "dynamic_tip_radius": True,
            "dynamic_front_width": True,
            "persistent_sites": True,
            "finite_source_inventory": False,
            "source_refresh": False,
            "explicit_recovery": False,
        }
        return payload


def main(argv=None):
    saved = {
        "point_release": _base.POINT_RELEASE,
        "model_id": _base.MODEL_ID,
        "production_manifest": _base.PRODUCTION_MANIFEST,
        "selection_manifest": _base.SELECTION_MANIFEST,
        "transfer_manifest": _base.TRANSFER_MANIFEST,
        "seed_manifest": _base.SEED_MANIFEST,
        "hazard_schema": _base.EMISSION_HAZARD_SCHEMA,
        "transport_schema": _base.EMISSION_TRANSPORT_SCHEMA,
        "engine": _base.ProductionStochasticEmissionFrontEngineV100518,
    }
    _base.POINT_RELEASE = POINT_RELEASE
    _base.MODEL_ID = MODEL_ID
    _base.PRODUCTION_MANIFEST = PRODUCTION_MANIFEST
    _base.SELECTION_MANIFEST = SELECTION_MANIFEST
    _base.TRANSFER_MANIFEST = TRANSFER_MANIFEST
    _base.SEED_MANIFEST = SEED_MANIFEST
    _base.EMISSION_HAZARD_SCHEMA = EMISSION_HAZARD_SCHEMA
    _base.EMISSION_TRANSPORT_SCHEMA = EMISSION_TRANSPORT_SCHEMA
    _base.ProductionStochasticEmissionFrontEngineV100518 = (
        ProductionStochasticOnsetBurstFrontEngineV1005181
    )
    try:
        return _base.main(argv)
    finally:
        _base.POINT_RELEASE = saved["point_release"]
        _base.MODEL_ID = saved["model_id"]
        _base.PRODUCTION_MANIFEST = saved["production_manifest"]
        _base.SELECTION_MANIFEST = saved["selection_manifest"]
        _base.TRANSFER_MANIFEST = saved["transfer_manifest"]
        _base.SEED_MANIFEST = saved["seed_manifest"]
        _base.EMISSION_HAZARD_SCHEMA = saved["hazard_schema"]
        _base.EMISSION_TRANSPORT_SCHEMA = saved["transport_schema"]
        _base.ProductionStochasticEmissionFrontEngineV100518 = saved["engine"]


if __name__ == "__main__":
    main()


__all__ = [
    "DEFAULT_PARAMETER_SOURCE_ROOT",
    "MODEL_ID",
    "POINT_RELEASE",
    "PRODUCTION_MANIFEST",
    "ProductionStochasticOnsetBurstFrontEngineV1005181",
    "SEED_MANIFEST",
    "SELECTION_MANIFEST",
    "TRANSFER_MANIFEST",
    "main",
]
