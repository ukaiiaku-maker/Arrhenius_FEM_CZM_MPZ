"""v10.0.5.14.4: long-growth persistent-site family production correction."""
from __future__ import annotations

import json
from pathlib import Path
import sys
from typing import Any

import numpy as np

from . import mode_i_first_passage_v10_0_5_14_1_persistent_site_family as _base
from . import persistent_site_front_engine_v100514 as _engine
from .adaptive_czm_tip_support_v1005144 import (
    MODEL_ID as TIP_SUPPORT_MODEL,
    installed_tip_support_repair_v1005144,
    tip_support_audit_v1005144,
)
from .persistent_site_diagnostics_v1005144 import (
    DIAGNOSTIC_MODEL,
    installed_persistent_diagnostics_v1005144,
)
from .persistent_site_transport_v1005144 import (
    TRANSPORT_INTEGRATOR,
    installed_split_transport_v1005144,
)

POINT_RELEASE = "10.0.5.14.4"
MODEL_ID = (
    "FEM_CZM_full_2D_PF_v10_2_22_persistent_site_kernel_family_"
    "physical_Lstable_transport_supported_tip_v10_0_5_14_4"
)
PRODUCTION_MANIFEST = "persistent_site_production_manifest_v10_0_5_14_4.json"
SELECTION_MANIFEST = "persistent_site_parameter_selection_v10_0_5_14_4.json"


def _out_path(argv: list[str]) -> Path | None:
    if "--out" not in argv:
        return None
    index = argv.index("--out")
    if index + 1 >= len(argv):
        return None
    return Path(argv[index + 1]).expanduser().resolve()


def _sum_nested(value: Any) -> float:
    try:
        return float(np.sum(np.asarray(value, dtype=float)))
    except Exception:
        return 0.0


def _final_line_content_accounting(out: Path) -> dict[str, Any] | None:
    files = sorted(out.glob("mpz_state_snapshots_*K.json"))
    if not files:
        return None
    records: list[dict[str, Any]] = []
    for path in files:
        try:
            payload = json.loads(path.read_text())
        except Exception:
            continue
        for front in payload.get("final_fronts", []):
            state = dict(front.get("state", {}) or {})
            active_mobile = _sum_nested(state.get("mobile_positive", [])) + _sum_nested(
                state.get("mobile_negative", [])
            )
            active_retained = _sum_nested(
                state.get("retained_positive", [])
            ) + _sum_nested(state.get("retained_negative", []))
            wake_mobile = _sum_nested(
                state.get("wake_mobile_positive", [])
            ) + _sum_nested(state.get("wake_mobile_negative", []))
            wake_retained = _sum_nested(
                state.get("wake_retained_positive", [])
            ) + _sum_nested(state.get("wake_retained_negative", []))
            emitted = float(state.get("emitted_total", 0.0) or 0.0)
            escaped = float(state.get("escaped_total", 0.0) or 0.0)
            recovered = float(state.get("recovered_total", 0.0) or 0.0)
            discarded = float(state.get("wake_discarded_total", 0.0) or 0.0)
            accounted = (
                active_mobile
                + active_retained
                + wake_mobile
                + wake_retained
                + escaped
                + recovered
                + discarded
            )
            signed_balance = emitted - accounted
            scale = max(abs(emitted), abs(accounted), 1.0e-300)
            records.append(
                {
                    "temperature_K": payload.get("temperature_K"),
                    "front_id": front.get("front_id"),
                    "active_mobile_total": active_mobile,
                    "active_retained_total": active_retained,
                    "wake_mobile_total": wake_mobile,
                    "wake_retained_total": wake_retained,
                    "emitted_total": emitted,
                    "escaped_total": escaped,
                    "recovered_total": recovered,
                    "wake_discarded_total": discarded,
                    "accounted_line_content_total": accounted,
                    "line_content_balance_signed": signed_balance,
                    "line_content_balance_error": abs(signed_balance),
                    "line_content_balance_relative_error": abs(signed_balance) / scale,
                    "advance_total_m": float(state.get("advance_total_m", 0.0) or 0.0),
                }
            )
    if not records:
        return None
    return {
        "schema": "final_persistent_line_content_accounting_v10_0_5_14_4",
        "legacy_N_em_semantics": "instantaneous_active_retained_line_content",
        "authoritative_cumulative_emission_field": "emitted_total",
        "fronts": records,
        "maximum_relative_balance_error": max(
            row["line_content_balance_relative_error"] for row in records
        ),
    }


def _rewrite_release_metadata(out: Path | None) -> None:
    if out is None or not out.exists():
        return
    manifest = out / PRODUCTION_MANIFEST
    if manifest.is_file():
        payload = json.loads(manifest.read_text())
        payload["schema"] = "persistent_site_production_manifest_v10_0_5_14_4"
        payload["model"] = MODEL_ID
        payload["point_release"] = POINT_RELEASE
        physics = dict(payload.get("physics_contract", {}))
        physics.update(
            {
                "transport_integrator": TRANSPORT_INTEGRATOR,
                "transport_linear_state": "physical_mobile_and_retained_only",
                "transport_linear_solution": "sparse_L_stable_backward_euler",
                "nonlinear_transport_error_control": "step_doubling",
                "negligible_active_tail_acceptance": True,
                "physical_state_only_conservation_audit": True,
                "diagnostic_accumulators_excluded_from_transport_solve": True,
                "explicit_transport_CFL_microstepping": False,
                "transport_equations_changed": False,
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
    _engine.MODEL_ID = "FEM_CZM_persistent_site_front_engine_v10_0_5_14_4"
    try:
        with (
            installed_split_transport_v1005144(),
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
