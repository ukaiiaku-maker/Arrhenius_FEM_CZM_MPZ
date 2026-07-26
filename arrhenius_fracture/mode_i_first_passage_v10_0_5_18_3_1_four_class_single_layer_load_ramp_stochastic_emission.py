"""v10.0.5.18.3.1 four-class single-layer stochastic load-ramp FEM/CZM."""
from __future__ import annotations

import json
from pathlib import Path
import sys

from . import crack_backend as _crack_backend
from . import mixed_mode_first_passage_v8 as _mm
from . import mode_i_first_passage_v10_0_5_18_four_class_stochastic_emission as _v18
from .numerical_resilience_v1005183 import (
    CZM_RETRY_SCHEMA,
    PROBE_RETRY_SCHEMA,
    RetryingAdaptiveCZMBackendV1005183,
    make_robust_process_zone_traction_probe,
)
from .persistent_site_load_ramp_stochastic_emission_v10051831 import (
    LOAD_RAMP_SCHEMA,
    OUTER_EMISSION_LIMITER_SCHEMA,
    PROBE_FALLBACK_SCHEMA,
    PersistentSiteSingleLayerLoadRampStochasticEmissionFrontEngineV10051831,
    two_channel_absolute_opening_drives_v1005183,
)

POINT_RELEASE = "10.0.5.18.3.1"
MODEL_ID = "FEM_CZM_four_class_exact_single_layer_load_ramp_v10_0_5_18_3_1"
PRODUCTION_MANIFEST = "persistent_site_production_manifest_v10_0_5_18_3_1.json"
SELECTION_MANIFEST = "persistent_site_parameter_selection_v10_0_5_18_3_1.json"
TRANSFER_MANIFEST = "four_class_parameter_transfer_v10_0_5_18_3_1.json"
SEED_MANIFEST = "stochastic_seed_manifest_v10_0_5_18_3_1.json"


class ProductionSingleLayerLoadRampStochasticEmissionFrontEngineV10051831(
    PersistentSiteSingleLayerLoadRampStochasticEmissionFrontEngineV10051831
):
    """Production adapter retaining the instantaneous aggregate hazard audit."""

    def step_drives(self, K_cleave, K_emit, T, dt, metadata=None):
        out = super().step_drives(
            K_cleave,
            K_emit,
            T,
            dt,
            metadata=metadata,
        )
        last = dict(self.mpz_state.last_emission or {})
        aggregate = last.get(
            "aggregate_hazard_final_by_system_s",
            last.get("aggregate_hazard_initial_by_system_s", [0.0, 0.0]),
        )
        out["lambda_e"] = float(sum(float(value) for value in aggregate))
        out["lambda_e_semantics"] = (
            "instantaneous_sum_exact_stochastic_signed_channel_aggregate_hazards"
        )
        return out


def _update_json(path: Path, update) -> None:
    if not path.is_file():
        return
    payload = json.loads(path.read_text())
    update(payload)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n")


def _rewrite_point_release_metadata(
    original_rewrite,
    out: Path | None,
    audit: dict,
    solver_args: list[str],
) -> None:
    original_rewrite(out, audit, solver_args)
    if out is None:
        return

    def production(payload: dict) -> None:
        payload["schema"] = "persistent_site_production_manifest_v10_0_5_18_3_1"
        payload["model"] = MODEL_ID
        payload["point_release"] = POINT_RELEASE
        physics = dict(payload.get("physics_contract", {}) or {})
        physics.update(
            {
                "parameter_transfer_exact": True,
                "stochastic_emission": True,
                "stochastic_emission_accepted_update": (
                    "exact_event_localized_one_activation_packets"
                ),
                "stochastic_emission_post_onset_mean_field_burst": False,
                "global_FEM_deterministic_emission_state_predictor": False,
                "local_linear_K_ramp_kinetics": True,
                "local_load_ramp_schema": LOAD_RAMP_SCHEMA,
                "outer_emission_action_limiter_active": False,
                "outer_emission_limiter_schema": OUTER_EMISSION_LIMITER_SCHEMA,
                "inner_exact_emission_action_control": True,
                "emission_scalar_directional_factor_applied_to_opening_K": False,
                "emission_signed_channel_factors_applied_once": True,
                "emission_backstress_is_physical_limiter": True,
                "dynamic_tip_radius": True,
                "dynamic_front_width": True,
                "dynamic_source_area": True,
                "dynamic_site_multiplicity": True,
                "robust_tensor_probe_retry": True,
                "tensor_probe_retry_schema": PROBE_RETRY_SCHEMA,
                "bounded_last_reliable_tensor_probe_fallback": True,
                "tensor_probe_fallback_schema": PROBE_FALLBACK_SCHEMA,
                "adaptive_CZM_atomic_partition_retry": True,
                "adaptive_CZM_partition_retry_schema": CZM_RETRY_SCHEMA,
                "physical_crack_event_length_preserved_across_topology_subsegments": True,
            }
        )
        payload["physics_contract"] = physics
        payload["front_engine"] = (
            ProductionSingleLayerLoadRampStochasticEmissionFrontEngineV10051831.audit_payload()
        )

    def selection(payload: dict) -> None:
        payload["schema"] = MODEL_ID
        payload["point_release"] = POINT_RELEASE
        policy = dict(payload.get("policy", {}) or {})
        policy.update(
            {
                "exact_stochastic_emission_events": True,
                "post_onset_mean_field_burst": False,
                "global_deterministic_state_predictor": False,
                "local_linear_K_ramp_kinetics": True,
                "load_ramp_schema": LOAD_RAMP_SCHEMA,
                "outer_emission_action_limiter_active": False,
                "inner_exact_emission_action_control": True,
                "two_channel_absolute_opening_drive": True,
                "tensor_probe_retry_schema": PROBE_RETRY_SCHEMA,
                "CZM_partition_retry_schema": CZM_RETRY_SCHEMA,
            }
        )
        payload["policy"] = policy

    def transfer(payload: dict) -> None:
        payload["schema"] = MODEL_ID
        payload["model"] = MODEL_ID
        payload["point_release"] = POINT_RELEASE
        payload["exact_stochastic_emission_events"] = True
        payload["post_onset_mean_field_burst"] = False
        payload["global_deterministic_state_predictor"] = False
        payload["local_linear_K_ramp_kinetics"] = True
        payload["load_ramp_schema"] = LOAD_RAMP_SCHEMA
        payload["outer_emission_action_limiter_active"] = False
        payload["inner_exact_emission_action_control"] = True
        payload["scalar_emission_factor_applied_to_opening_K"] = False
        payload["signed_channel_factors_applied_once"] = True
        payload["tensor_probe_retry_schema"] = PROBE_RETRY_SCHEMA
        payload["CZM_partition_retry_schema"] = CZM_RETRY_SCHEMA

    def seeds(payload: dict) -> None:
        payload["schema"] = (
            "v10.0.5.18.3.1_domain_separated_stochastic_seed_manifest"
        )
        payload["point_release"] = POINT_RELEASE

    _update_json(out / PRODUCTION_MANIFEST, production)
    _update_json(out / SELECTION_MANIFEST, selection)
    _update_json(out / TRANSFER_MANIFEST, transfer)
    _update_json(out / SEED_MANIFEST, seeds)


def main(argv: list[str] | None = None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    saved = {
        "point_release": _v18.POINT_RELEASE,
        "model_id": _v18.MODEL_ID,
        "production_manifest": _v18.PRODUCTION_MANIFEST,
        "selection_manifest": _v18.SELECTION_MANIFEST,
        "transfer_manifest": _v18.TRANSFER_MANIFEST,
        "seed_manifest": _v18.SEED_MANIFEST,
        "engine": _v18.ProductionStochasticEmissionFrontEngineV100518,
        "rewrite": _v18._rewrite_output_metadata,
        "mm_drives": _mm.CalibratedTipEngineMixin._mm_drives,
        "traction_probe": _mm.process_zone_traction_probe,
        "adaptive_czm": _crack_backend.AdaptiveCZMBackend,
    }

    _v18.POINT_RELEASE = POINT_RELEASE
    _v18.MODEL_ID = MODEL_ID
    _v18.PRODUCTION_MANIFEST = PRODUCTION_MANIFEST
    _v18.SELECTION_MANIFEST = SELECTION_MANIFEST
    _v18.TRANSFER_MANIFEST = TRANSFER_MANIFEST
    _v18.SEED_MANIFEST = SEED_MANIFEST
    _v18.ProductionStochasticEmissionFrontEngineV100518 = (
        ProductionSingleLayerLoadRampStochasticEmissionFrontEngineV10051831
    )
    _mm.CalibratedTipEngineMixin._mm_drives = (
        two_channel_absolute_opening_drives_v1005183
    )
    _mm.process_zone_traction_probe = make_robust_process_zone_traction_probe(
        saved["traction_probe"]
    )
    _crack_backend.AdaptiveCZMBackend = RetryingAdaptiveCZMBackendV1005183
    original_rewrite = saved["rewrite"]
    _v18._rewrite_output_metadata = lambda out, audit, solver_args: (
        _rewrite_point_release_metadata(
            original_rewrite,
            out,
            audit,
            solver_args,
        )
    )

    try:
        return _v18.main(user_args)
    finally:
        _v18.POINT_RELEASE = saved["point_release"]
        _v18.MODEL_ID = saved["model_id"]
        _v18.PRODUCTION_MANIFEST = saved["production_manifest"]
        _v18.SELECTION_MANIFEST = saved["selection_manifest"]
        _v18.TRANSFER_MANIFEST = saved["transfer_manifest"]
        _v18.SEED_MANIFEST = saved["seed_manifest"]
        _v18.ProductionStochasticEmissionFrontEngineV100518 = saved["engine"]
        _v18._rewrite_output_metadata = saved["rewrite"]
        _mm.CalibratedTipEngineMixin._mm_drives = saved["mm_drives"]
        _mm.process_zone_traction_probe = saved["traction_probe"]
        _crack_backend.AdaptiveCZMBackend = saved["adaptive_czm"]


if __name__ == "__main__":
    main()


__all__ = [
    "MODEL_ID",
    "POINT_RELEASE",
    "PRODUCTION_MANIFEST",
    "ProductionSingleLayerLoadRampStochasticEmissionFrontEngineV10051831",
    "SEED_MANIFEST",
    "SELECTION_MANIFEST",
    "TRANSFER_MANIFEST",
    "main",
]
