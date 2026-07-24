"""v10.0.5.14.5: asymptotic-preserving low-temperature production release."""
from __future__ import annotations

import json
from pathlib import Path
import sys

from . import mode_i_first_passage_v10_0_5_14_1_persistent_site_family as _base
from . import persistent_site_front_engine_v100514 as _engine
from .adaptive_czm_tip_support_v1005144 import (
    MODEL_ID as TIP_SUPPORT_MODEL,
    installed_tip_support_repair_v1005144,
    tip_support_audit_v1005144,
)
from .mode_i_first_passage_v10_0_5_14_4_persistent_site_family import (
    _final_line_content_accounting,
)
from .persistent_site_diagnostics_v1005144 import (
    DIAGNOSTIC_MODEL,
    installed_persistent_diagnostics_v1005144,
)
from .persistent_site_transport_v1005145 import (
    ASYMPTOTIC_MODEL,
    TRANSPORT_INTEGRATOR,
    installed_asymptotic_transport_v1005145,
)

POINT_RELEASE = "10.0.5.14.5"
MODEL_ID = (
    "FEM_CZM_full_2D_PF_v10_2_22_persistent_site_kernel_family_"
    "asymptotic_lowT_transport_supported_tip_v10_0_5_14_5"
)
PRODUCTION_MANIFEST = "persistent_site_production_manifest_v10_0_5_14_5.json"
SELECTION_MANIFEST = "persistent_site_parameter_selection_v10_0_5_14_5.json"


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
    manifest = out / PRODUCTION_MANIFEST
    if manifest.is_file():
        payload = json.loads(manifest.read_text())
        payload["schema"] = "persistent_site_production_manifest_v10_0_5_14_5"
        payload["model"] = MODEL_ID
        payload["point_release"] = POINT_RELEASE
        physics = dict(payload.get("physics_contract", {}))
        physics.update(
            {
                "transport_integrator": TRANSPORT_INTEGRATOR,
                "low_temperature_asymptotic_model": ASYMPTOTIC_MODEL,
                "asymptotic_selection": (
                    "measured_mobile_fast_courant_and_fast_to_release_separation"
                ),
                "asymptotic_mobile_competing_risks": (
                    "encounter_storage_vs_forward_transport_and_escape"
                ),
                "asymptotic_retained_dynamics": (
                    "exact_substochastic_retained_generator_exponential"
                ),
                "state_dependent_rate_closure": "damped_fixed_point",
                "nonasymptotic_fallback": (
                    "adaptive_physical_backward_euler_tail_control_v10_0_5_14_4"
                ),
                "recursive_resolution_of_separated_fast_mobile_layer": False,
                "transport_equations_changed": False,
                "transport_constitutive_rates_changed": False,
                "transport_time_integrator_changed": True,
                "adaptive_czm_tip_support_model": TIP_SUPPORT_MODEL,
                "adaptive_czm_authoritative_tip_pair": True,
                "adaptive_czm_quality_veto_retained": True,
                "persistent_state_diagnostic_model": DIAGNOSTIC_MODEL,
                "legacy_N_em_is_cumulative_emission": False,
            }
        )
        payload["physics_contract"] = physics
        payload["adaptive_czm_tip_support_audit"] = tip_support_audit_v1005144()
        final_accounting = _final_line_content_accounting(out)
        if final_accounting is not None:
            final_accounting["schema"] = (
                "final_persistent_line_content_accounting_v10_0_5_14_5"
            )
            payload["final_line_content_accounting"] = final_accounting
        manifest.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

    old_selection = out / "persistent_site_parameter_selection_v10_0_5_14_1.json"
    new_selection = out / SELECTION_MANIFEST
    if old_selection.is_file():
        payload = json.loads(old_selection.read_text())
        payload["schema"] = MODEL_ID
        payload["point_release"] = POINT_RELEASE
        policy = dict(payload.get("policy", {}))
        policy["transport_integrator"] = TRANSPORT_INTEGRATOR
        policy["transport_cfl_limited"] = False
        policy["low_temperature_asymptotic_model"] = ASYMPTOTIC_MODEL
        policy["adaptive_czm_tip_support_model"] = TIP_SUPPORT_MODEL
        policy["persistent_state_diagnostic_model"] = DIAGNOSTIC_MODEL
        policy["candidate_source_or_shielding_closure_applied"] = True
        policy["barrier_only_transfer"] = False
        payload["policy"] = policy
        new_selection.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        old_selection.unlink()


def main(argv: list[str] | None = None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    out = _out_path(user_args)
    saved = {
        "point_release": _base.POINT_RELEASE,
        "model_id": _base.MODEL_ID,
        "manifest": _base.PRODUCTION_MANIFEST,
        "engine_model": _engine.MODEL_ID,
    }
    _base.POINT_RELEASE = POINT_RELEASE
    _base.MODEL_ID = MODEL_ID
    _base.PRODUCTION_MANIFEST = PRODUCTION_MANIFEST
    _engine.MODEL_ID = "FEM_CZM_persistent_site_front_engine_v10_0_5_14_5"
    try:
        with (
            installed_asymptotic_transport_v1005145(),
            installed_tip_support_repair_v1005144(),
            installed_persistent_diagnostics_v1005144(),
        ):
            return _base.main(user_args)
    finally:
        _rewrite_release_metadata(out)
        _base.POINT_RELEASE = saved["point_release"]
        _base.MODEL_ID = saved["model_id"]
        _base.PRODUCTION_MANIFEST = saved["manifest"]
        _engine.MODEL_ID = saved["engine_model"]


if __name__ == "__main__":
    main()


__all__ = ["POINT_RELEASE", "MODEL_ID", "PRODUCTION_MANIFEST", "main"]
