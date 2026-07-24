"""v10.0.5.15: PF update-map and continuous moving-tip parity release."""
from __future__ import annotations

import json
from pathlib import Path
import sys

from . import mode_i_first_passage_v10_0_5_14_1_persistent_site_family as _base
from .adaptive_czm_tip_support_v1005144 import (
    MODEL_ID as TIP_SUPPORT_MODEL,
    installed_tip_support_repair_v1005144,
    tip_support_audit_v1005144,
)
from .persistent_site_diagnostics_v1005144 import (
    DIAGNOSTIC_MODEL,
    installed_persistent_diagnostics_v1005144,
)
from .persistent_site_moving_tip_v100515 import (
    COUPLING_SCHEME,
    MODEL_ID as ENGINE_MODEL,
    PersistentSitePFMovingTipFrontEngineV100515,
)
from .persistent_site_pf_update_v100515 import (
    PF_REFERENCE_COMMIT,
    PF_UPDATE_MAP,
)
from .mode_i_first_passage_v10_0_5_14_4_persistent_site_family import (
    _final_line_content_accounting,
)

POINT_RELEASE = "10.0.5.15"
MODEL_ID = (
    "FEM_CZM_full_2D_PF_v10_2_22_exact_update_map_"
    "continuous_moving_tip_v10_0_5_15"
)
PRODUCTION_MANIFEST = "persistent_site_production_manifest_v10_0_5_15.json"
SELECTION_MANIFEST = "persistent_site_parameter_selection_v10_0_5_15.json"


def _out_path(argv: list[str]) -> Path | None:
    if "--out" not in argv:
        return None
    index = argv.index("--out")
    if index + 1 >= len(argv):
        return None
    return Path(argv[index + 1]).expanduser().resolve()


def _rewrite_release_metadata(out: Path | None) -> None:
    if out is None or not out.exists():
        return
    inherited_manifest = out / _base.PRODUCTION_MANIFEST
    manifest = out / PRODUCTION_MANIFEST
    if inherited_manifest.is_file():
        payload = json.loads(inherited_manifest.read_text())
        payload["schema"] = "persistent_site_production_manifest_v10_0_5_15"
        payload["model"] = MODEL_ID
        payload["point_release"] = POINT_RELEASE
        physics = dict(payload.get("physics_contract", {}))
        physics.update(
            {
                "PF_reference_commit": PF_REFERENCE_COMMIT,
                "PF_local_update_map": PF_UPDATE_MAP,
                "transport_operator_order": (
                    "persistent_emission_then_exact_exchange_then_zero_recovery_"
                    "then_population_weighted_scalar_advection"
                ),
                "transport_forest_population": "unsigned_retained_only",
                "transport_velocity": "population_weighted_scalar_per_active_MPZ",
                "coupled_crack_tip_scheme": COUPLING_SCHEME,
                "plastic_half_step_before_cleavage_rate": True,
                "fractional_MPZ_translation_per_cleavage_progress": True,
                "plastic_half_step_after_fractional_translation": True,
                "microstructure_advance_precedes_cohesive_checkpoint": True,
                "checkpoint_MPZ_advance_repeated": False,
                "nonmutating_coupled_clock_prediction": True,
                "cohesive_geometry_checkpoint_distance": "front_config.da",
                "continuous_microadvance_distance": "front_config.da*dB",
                "geometry_veto_policy": "restore_renewal_origin_and_fail_closed",
                "adaptive_czm_tip_support_model": TIP_SUPPORT_MODEL,
                "adaptive_czm_authoritative_tip_pair": True,
                "adaptive_czm_quality_veto_retained": True,
                "persistent_state_diagnostic_model": DIAGNOSTIC_MODEL,
                "bulk_plasticity_mode": "tip_only",
                "constitutive_parameters_changed_from_PF": False,
            }
        )
        payload["physics_contract"] = physics
        payload["front_engine"] = PersistentSitePFMovingTipFrontEngineV100515.audit_payload()
        payload["adaptive_czm_tip_support_audit"] = tip_support_audit_v1005144()
        final_accounting = _final_line_content_accounting(out)
        if final_accounting is not None:
            final_accounting["schema"] = (
                "final_persistent_line_content_accounting_v10_0_5_15"
            )
            payload["final_line_content_accounting"] = final_accounting
        manifest.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        if inherited_manifest != manifest:
            inherited_manifest.unlink()

    old_selection = out / "persistent_site_parameter_selection_v10_0_5_14_1.json"
    if old_selection.is_file():
        selection = json.loads(old_selection.read_text())
        selection["schema"] = MODEL_ID
        selection["point_release"] = POINT_RELEASE
        policy = dict(selection.get("policy", {}))
        policy.update(
            {
                "PF_reference_commit": PF_REFERENCE_COMMIT,
                "PF_local_update_map": PF_UPDATE_MAP,
                "kinetic_tip_cell_active": True,
                "kinetic_coupling_scheme": COUPLING_SCHEME,
                "continuous_fractional_MPZ_translation": True,
                "microstructure_advance_precedes_cohesive_checkpoint": True,
                "checkpoint_MPZ_advance_repeated": False,
                "candidate_source_or_shielding_closure_applied": True,
                "barrier_only_transfer": False,
            }
        )
        selection["policy"] = policy
        (out / SELECTION_MANIFEST).write_text(
            json.dumps(selection, indent=2, sort_keys=True) + "\n"
        )
        old_selection.unlink()


def main(argv: list[str] | None = None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    out = _out_path(user_args)
    saved = {
        "point_release": _base.POINT_RELEASE,
        "model_id": _base.MODEL_ID,
        "manifest": _base.PRODUCTION_MANIFEST,
        "engine_class": _base.PersistentSiteMovingProcessZoneFrontEngineV100514,
    }
    _base.POINT_RELEASE = POINT_RELEASE
    _base.MODEL_ID = MODEL_ID
    _base.PRODUCTION_MANIFEST = PRODUCTION_MANIFEST
    _base.PersistentSiteMovingProcessZoneFrontEngineV100514 = (
        PersistentSitePFMovingTipFrontEngineV100515
    )
    try:
        with (
            installed_tip_support_repair_v1005144(),
            installed_persistent_diagnostics_v1005144(),
        ):
            return _base.main(user_args)
    finally:
        _rewrite_release_metadata(out)
        _base.POINT_RELEASE = saved["point_release"]
        _base.MODEL_ID = saved["model_id"]
        _base.PRODUCTION_MANIFEST = saved["manifest"]
        _base.PersistentSiteMovingProcessZoneFrontEngineV100514 = saved[
            "engine_class"
        ]


if __name__ == "__main__":
    main()


__all__ = ["POINT_RELEASE", "MODEL_ID", "PRODUCTION_MANIFEST", "main"]
