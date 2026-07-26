"""v10.0.5.18.2 four-class state-coupled exact stochastic emission."""
from __future__ import annotations

import json
from pathlib import Path
import sys

from . import mixed_mode_first_passage_v8 as _mm
from . import mode_i_first_passage_v10_0_5_18_four_class_stochastic_emission as _v18
from .persistent_site_state_coupled_stochastic_emission_v1005182 import (
    ADAPTIVE_SCHEMA,
    PersistentSiteStateCoupledStochasticEmissionFrontEngineV1005182,
    two_channel_absolute_opening_drives,
)

POINT_RELEASE = "10.0.5.18.2"
MODEL_ID = "FEM_CZM_four_class_state_coupled_exact_stochastic_emission_v10_0_5_18_2"
PRODUCTION_MANIFEST = "persistent_site_production_manifest_v10_0_5_18_2.json"
SELECTION_MANIFEST = "persistent_site_parameter_selection_v10_0_5_18_2.json"
TRANSFER_MANIFEST = "four_class_parameter_transfer_v10_0_5_18_2.json"
SEED_MANIFEST = "stochastic_seed_manifest_v10_0_5_18_2.json"


class ProductionStateCoupledStochasticEmissionFrontEngineV1005182(
    PersistentSiteStateCoupledStochasticEmissionFrontEngineV1005182
):
    """Production adapter preserving the instantaneous lambda_e diagnostic."""

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
        payload["schema"] = "persistent_site_production_manifest_v10_0_5_18_2"
        payload["model"] = MODEL_ID
        payload["point_release"] = POINT_RELEASE
        physics = dict(payload.get("physics_contract", {}) or {})
        physics.update(
            {
                "parameter_transfer_exact": True,
                "four_class_parameter_transfer_only": False,
                "stochastic_emission": True,
                "stochastic_emission_accepted_update": (
                    "exact_event_localized_one_activation_packets"
                ),
                "stochastic_emission_post_onset_mean_field_burst": False,
                "emission_outer_load_state_adaptivity": True,
                "emission_adaptive_schema": ADAPTIVE_SCHEMA,
                "emission_scalar_directional_factor_applied_to_opening_K": False,
                "emission_signed_channel_factors_applied_once": True,
                "emission_backstress_is_physical_limiter": True,
                "emission_count_cap_applied_to_physics": False,
                "dynamic_tip_radius": True,
                "dynamic_front_width": True,
                "dynamic_source_area": True,
                "dynamic_site_multiplicity": True,
            }
        )
        payload["physics_contract"] = physics
        payload["front_engine"] = (
            ProductionStateCoupledStochasticEmissionFrontEngineV1005182.audit_payload()
        )

    def selection(payload: dict) -> None:
        payload["schema"] = MODEL_ID
        payload["point_release"] = POINT_RELEASE
        policy = dict(payload.get("policy", {}) or {})
        policy.update(
            {
                "exact_stochastic_emission_events": True,
                "post_onset_mean_field_burst": False,
                "emission_outer_load_state_adaptivity": True,
                "emission_adaptive_schema": ADAPTIVE_SCHEMA,
                "two_channel_absolute_opening_drive": True,
            }
        )
        payload["policy"] = policy

    def transfer(payload: dict) -> None:
        payload["schema"] = MODEL_ID
        payload["model"] = MODEL_ID
        payload["point_release"] = POINT_RELEASE
        payload["exact_stochastic_emission_events"] = True
        payload["post_onset_mean_field_burst"] = False
        payload["emission_outer_load_state_adaptivity"] = True
        payload["emission_adaptive_schema"] = ADAPTIVE_SCHEMA
        payload["scalar_emission_factor_applied_to_opening_K"] = False
        payload["signed_channel_factors_applied_once"] = True

    def seeds(payload: dict) -> None:
        payload["schema"] = (
            "v10.0.5.18.2_domain_separated_stochastic_seed_manifest"
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
    }

    _v18.POINT_RELEASE = POINT_RELEASE
    _v18.MODEL_ID = MODEL_ID
    _v18.PRODUCTION_MANIFEST = PRODUCTION_MANIFEST
    _v18.SELECTION_MANIFEST = SELECTION_MANIFEST
    _v18.TRANSFER_MANIFEST = TRANSFER_MANIFEST
    _v18.SEED_MANIFEST = SEED_MANIFEST
    _v18.ProductionStochasticEmissionFrontEngineV100518 = (
        ProductionStateCoupledStochasticEmissionFrontEngineV1005182
    )
    _mm.CalibratedTipEngineMixin._mm_drives = two_channel_absolute_opening_drives
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


if __name__ == "__main__":
    main()


__all__ = [
    "MODEL_ID",
    "POINT_RELEASE",
    "PRODUCTION_MANIFEST",
    "ProductionStateCoupledStochasticEmissionFrontEngineV1005182",
    "SEED_MANIFEST",
    "SELECTION_MANIFEST",
    "TRANSFER_MANIFEST",
    "main",
]
