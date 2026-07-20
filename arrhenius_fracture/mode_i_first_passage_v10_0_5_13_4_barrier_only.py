"""v10.0.5.13.4 tip-only policy propagation repair.

v10.0.5.13.3 correctly selected ``tip_only`` in the campaign command, but the
original v10.0.5.13 entry still compared that command against its module-local
``bulk_same_pt_km`` policy before delegating to the established v9.11 solver.
This point release propagates the tip-only policy through the complete wrapper
chain for the duration of one run and restores every imported policy afterward.
No FEM, CZM, MPZ, barrier, source, shielding, transport, or crack-advance law is
changed.
"""
from __future__ import annotations

from datetime import datetime, timezone
import json
from pathlib import Path
import sys
from typing import Any, Mapping

from . import barrier_only_response_registry_v100513 as _registry
from . import mode_i_first_passage_v10_0_5_13_3_barrier_only as _base
from .mode_i_first_passage_v10_0 import _option_value

POINT_RELEASE = "10.0.5.13.4"
MODEL_ID = "FEM_CZM_full_2D_barrier_only_tip_source_MPZ_v10_0_5_13_4"
PRODUCTION_MANIFEST = _base.PRODUCTION_MANIFEST

# Wrapper chain: 13.4 -> 13.3 -> 13.2 -> 13.1 -> 13.0.
_PRESERVED_ENTRY = _base._base._base
_CORE_ENTRY = _base._base._base._base
_ORIGINAL_COMPACT_AUDIT = _CORE_ENTRY._compact_audit

TIP_ONLY_POLICY: dict[str, Any] = {
    **dict(_registry.TWO_D_STATE_POLICY),
    "policy_id": "preserve_existing_tip_only_moving_mpz_v1005134",
    "bulk_plasticity_mode": "tip_only",
    "state_configuration_source": "existing_tip_only_2d_solver_and_explicit_cli",
    "state_evolution_source": "existing_moving_crack_tip_MPZ",
    "continuum_bulk_role": "elastic_fem_only",
    "uniform_bulk_mobile_retained_state_active": False,
}


def _tip_only_compact_audit(row: Mapping[str, Any]) -> dict[str, Any]:
    payload = dict(_ORIGINAL_COMPACT_AUDIT(row))
    payload.update(
        {
            "two_d_state_policy": dict(TIP_ONLY_POLICY),
            "bulk_state_evolves_in_fem": False,
            "uniform_bulk_mobile_retained_state_active": False,
            "moving_crack_tip_mpz_active": True,
            "shielding_derived_from_evolving_state": True,
            "shielding_state_location": "moving_crack_tip_MPZ",
        }
    )
    return payload


def _update_manifest(out: Path, completed: bool) -> None:
    path = out / PRODUCTION_MANIFEST
    if not path.is_file():
        return
    try:
        payload = json.loads(path.read_text())
    except Exception:
        return
    barrier_option = payload.get("barrier_option")
    if isinstance(barrier_option, dict):
        barrier_option["two_d_state_policy"] = dict(TIP_ONLY_POLICY)
    payload.update(
        {
            "schema": "barrier_only_tip_source_manifest_v10_0_5_13_4",
            "model": MODEL_ID,
            "point_release": POINT_RELEASE,
            "two_d_state_policy": dict(TIP_ONLY_POLICY),
            "bulk_plasticity_mode": "tip_only",
            "tip_only_policy_propagation_repair": {
                "active": True,
                "core_entry_policy_overridden_for_run": True,
                "preserved_state_entry_policy_overridden_for_run": True,
                "registry_audit_policy_overridden_for_run": True,
                "all_policy_objects_restored_after_run": True,
                "constitutive_physics_changed": False,
                "recorded_utc": datetime.now(timezone.utc).isoformat(),
            },
            "plasticity_scope": {
                "bulk_plasticity_mode": "tip_only",
                "continuum_bulk_role": "elastic_fem_only",
                "moving_crack_tip_mpz_active": True,
                "uniform_bulk_mobile_retained_state_active": False,
                "candidate_source_or_shielding_closure_applied": False,
                "recorded_utc": datetime.now(timezone.utc).isoformat(),
            },
        }
    )
    if completed:
        payload["status"] = "complete"
        payload["run_completed_without_exception"] = True
    path.write_text(json.dumps(payload, indent=2, default=str))


def main(argv: list[str] | None = None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    out_value = _option_value(user_args, "--out")
    if out_value is None:
        raise SystemExit("v10.0.5.13.4 requires --out")
    mode = _option_value(user_args, "--bulk-plasticity-mode")
    if mode != "tip_only":
        raise SystemExit(
            "v10.0.5.13.4 requires --bulk-plasticity-mode tip_only; "
            f"received {mode!r}"
        )
    out = Path(out_value).resolve()

    saved_core_policy = _CORE_ENTRY.TWO_D_STATE_POLICY
    saved_preserved_policy = _PRESERVED_ENTRY.TWO_D_STATE_POLICY
    saved_registry_policy = _registry.TWO_D_STATE_POLICY
    saved_compact_audit = _CORE_ENTRY._compact_audit
    active_policy = dict(TIP_ONLY_POLICY)
    _CORE_ENTRY.TWO_D_STATE_POLICY = active_policy
    _PRESERVED_ENTRY.TWO_D_STATE_POLICY = active_policy
    _registry.TWO_D_STATE_POLICY = active_policy
    _CORE_ENTRY._compact_audit = _tip_only_compact_audit
    try:
        result = _base.main(user_args)
        _update_manifest(out, completed=True)
        return result
    except BaseException:
        _update_manifest(out, completed=False)
        raise
    finally:
        _CORE_ENTRY.TWO_D_STATE_POLICY = saved_core_policy
        _PRESERVED_ENTRY.TWO_D_STATE_POLICY = saved_preserved_policy
        _registry.TWO_D_STATE_POLICY = saved_registry_policy
        _CORE_ENTRY._compact_audit = saved_compact_audit


if __name__ == "__main__":
    main()


__all__ = [
    "POINT_RELEASE",
    "MODEL_ID",
    "PRODUCTION_MANIFEST",
    "TIP_ONLY_POLICY",
    "_tip_only_compact_audit",
    "main",
]
