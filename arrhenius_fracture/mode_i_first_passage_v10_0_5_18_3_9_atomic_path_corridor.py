"""v10.0.5.18.3.9 atomic target-aware path-corridor FEM/CZM entry."""
from __future__ import annotations

import json
from pathlib import Path
import sys
from typing import Any

from . import (
    mode_i_first_passage_v10_0_5_18_3_2_four_class_joint_K_ramp_stochastic_emission
    as _ramp
)
from . import (
    mode_i_first_passage_v10_0_5_18_3_4_four_class_path_aware_growth_envelope
    as _base
)
from .atomic_path_corridor_czm_v10051839 import (
    AtomicPathCorridorCZMBackendV10051839,
    MODEL_ID as CORRIDOR_MODEL,
    SCHEMA as CORRIDOR_SCHEMA,
    audit_payload,
    reset_audit,
)


POINT_RELEASE = "10.0.5.18.3.9"
MODEL_ID = "FEM_CZM_atomic_target_aware_path_corridor_v10_0_5_18_3_9"
PRODUCTION_MANIFEST = "persistent_site_production_manifest_v10_0_5_18_3_9.json"
SELECTION_MANIFEST = "persistent_site_parameter_selection_v10_0_5_18_3_9.json"
TRANSFER_MANIFEST = "four_class_parameter_transfer_v10_0_5_18_3_9.json"
SEED_MANIFEST = "stochastic_seed_manifest_v10_0_5_18_3_9.json"
CORRIDOR_AUDIT = "atomic_path_corridor_audit_v10_0_5_18_3_9.json"

_OLD_FILES = {
    _base.PRODUCTION_MANIFEST: PRODUCTION_MANIFEST,
    _base.SELECTION_MANIFEST: SELECTION_MANIFEST,
    _base.TRANSFER_MANIFEST: TRANSFER_MANIFEST,
    _base.SEED_MANIFEST: SEED_MANIFEST,
}


def _option_value(argv: list[str], name: str) -> str | None:
    for index, token in enumerate(argv):
        if token == name and index + 1 < len(argv):
            return argv[index + 1]
        if token.startswith(name + "="):
            return token.split("=", 1)[1]
    return None


def _fields() -> dict[str, Any]:
    return {
        "point_release": POINT_RELEASE,
        "model": MODEL_ID,
        "atomic_path_corridor_active": True,
        "atomic_path_corridor_schema": CORRIDOR_SCHEMA,
        "atomic_path_corridor_model": CORRIDOR_MODEL,
        "complete_event_remeshed_before_cohesive_commit": True,
        "pslg_constrained_exact_event_path": True,
        "exact_stochastic_event_endpoint_preserved": True,
        "exact_selected_crack_direction_preserved": True,
        "atomic_physical_event_transaction": True,
        "numerical_path_subsegments_are_independent_hazard_events": False,
        "equal_length_physical_partition_retry": False,
        "triangle_quality_floor_relaxed": False,
        "child_area_ratio_floor_relaxed": False,
        "tip_h_over_da_role": "audit_warning_only_not_event_veto",
        "hazard_or_event_length_changed": False,
        "constitutive_physics_changed": False,
        "bulk_plasticity_mode_required": "tip_only_elastic_FEM_bulk",
        "external_constrained_triangulator": "triangle_20250106",
    }


def _rewrite_payload(old_name: str, payload: dict[str, Any], audit: dict[str, Any]):
    fields = _fields()
    payload["point_release"] = POINT_RELEASE
    payload["model"] = MODEL_ID
    if old_name == _base.PRODUCTION_MANIFEST:
        payload["schema"] = "persistent_site_production_manifest_v10_0_5_18_3_9"
        physics = dict(payload.get("physics_contract", {}) or {})
        physics.update(fields)
        payload["physics_contract"] = physics
        payload["atomic_path_corridor_audit"] = {
            "event_count": int(audit.get("event_count", 0)),
            "success_count": int(audit.get("success_count", 0)),
            "failure_count": int(audit.get("failure_count", 0)),
        }
    elif old_name == _base.SELECTION_MANIFEST:
        payload["schema"] = MODEL_ID
        policy = dict(payload.get("policy", {}) or {})
        policy.update(fields)
        payload["policy"] = policy
    else:
        payload.update(fields)
    return payload


def _rewrite_outputs(out: Path | None, error: BaseException | None) -> None:
    if out is None:
        return
    out.mkdir(parents=True, exist_ok=True)
    audit = audit_payload()
    audit.update(
        {
            **_fields(),
            "run_completed_without_exception": error is None,
            "runtime_error_type": None if error is None else type(error).__name__,
            "runtime_error": None if error is None else str(error),
        }
    )
    (out / CORRIDOR_AUDIT).write_text(
        json.dumps(audit, indent=2, sort_keys=True, default=str) + "\n"
    )

    corridor_path = out / "compact_corridor_mesh_v91852.json"
    if corridor_path.is_file():
        payload = dict(json.loads(corridor_path.read_text()) or {})
        payload.update(_fields())
        corridor_path.write_text(
            json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n"
        )

    for old_name, new_name in _OLD_FILES.items():
        old_path = out / old_name
        if not old_path.is_file():
            continue
        payload = _rewrite_payload(old_name, json.loads(old_path.read_text()), audit)
        new_path = out / new_name
        new_path.write_text(
            json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n"
        )
        if new_path != old_path:
            old_path.unlink()


def main(argv: list[str] | None = None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    out_raw = _option_value(user_args, "--out")
    out = None if out_raw is None else Path(out_raw).expanduser().resolve()

    reset_audit()
    saved_backend = _ramp.RetryingAdaptiveCZMBackendV1005183
    _ramp.RetryingAdaptiveCZMBackendV1005183 = AtomicPathCorridorCZMBackendV10051839

    error: BaseException | None = None
    try:
        return _base.main(user_args)
    except BaseException as exc:
        error = exc
        raise
    finally:
        _ramp.RetryingAdaptiveCZMBackendV1005183 = saved_backend
        _rewrite_outputs(out, error)


if __name__ == "__main__":
    main()


__all__ = [
    "CORRIDOR_AUDIT",
    "MODEL_ID",
    "POINT_RELEASE",
    "PRODUCTION_MANIFEST",
    "SEED_MANIFEST",
    "SELECTION_MANIFEST",
    "TRANSFER_MANIFEST",
    "main",
]
