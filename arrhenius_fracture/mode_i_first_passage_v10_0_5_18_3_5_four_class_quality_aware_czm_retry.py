"""v10.0.5.18.3.5 production entry with quality-aware CZM retries."""
from __future__ import annotations

import json
from pathlib import Path
import sys
from typing import Any

from . import mode_i_first_passage_v9_18_5_6 as _v91856
from . import (
    mode_i_first_passage_v10_0_5_18_3_4_four_class_path_aware_growth_envelope
    as _base
)
from .quality_aware_czm_retry_v10051835 import (
    MODEL_ID as QUALITY_MODEL,
    SCHEMA as QUALITY_SCHEMA,
    audit_payload,
    configure_legacy_quality_wrapper,
    quality_aware_strict_advance_v10051835,
    reset_audit,
)
from .quality_aware_czm_ray_handoff_v10051835 import (
    install as install_ray_handoff,
)


POINT_RELEASE = "10.0.5.18.3.5"
MODEL_ID = (
    "FEM_CZM_four_class_joint_K_single_trial_path_aware_"
    "quality_aware_CZM_retry_v10_0_5_18_3_5"
)
PRODUCTION_MANIFEST = "persistent_site_production_manifest_v10_0_5_18_3_5.json"
SELECTION_MANIFEST = "persistent_site_parameter_selection_v10_0_5_18_3_5.json"
TRANSFER_MANIFEST = "four_class_parameter_transfer_v10_0_5_18_3_5.json"
SEED_MANIFEST = "stochastic_seed_manifest_v10_0_5_18_3_5.json"
QUALITY_AUDIT = "quality_aware_czm_retry_v10_0_5_18_3_5.json"

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


def _out_path(argv: list[str]) -> Path | None:
    raw = _option_value(argv, "--out")
    return None if raw is None else Path(raw).expanduser().resolve()


def _quality_fields() -> dict[str, Any]:
    return {
        "quality_aware_CZM_retry_active": True,
        "quality_aware_CZM_retry_schema": QUALITY_SCHEMA,
        "quality_aware_CZM_retry_model": QUALITY_MODEL,
        "quality_gate_evaluated_inside_retry_transaction": True,
        "shape_regular_local_tip_patch_bisection": True,
        "tip_fan_ray_handoff_after_safe_split": True,
        "adjacent_ray_target_cavity_refinement": True,
        "target_may_leave_tip_one_ring_during_refinement": True,
        "exact_stochastic_event_endpoint_preserved": True,
        "exact_selected_crack_direction_preserved": True,
        "collinear_partition_retry_after_patch_refinement": True,
        "immediate_parent_child_area_ratio_enforced": True,
        "cumulative_child_area_ratio_role": "audit_only",
        "triangle_quality_floor_relaxed": False,
        "child_area_ratio_floor_relaxed": False,
        "tip_h_over_da_role": "audit_warning_only",
        "constitutive_physics_changed": False,
    }


def _rewrite_payload(old_name: str, payload: dict[str, Any]) -> dict[str, Any]:
    payload["point_release"] = POINT_RELEASE
    payload["model"] = MODEL_ID
    fields = _quality_fields()

    if old_name == _base.PRODUCTION_MANIFEST:
        payload["schema"] = (
            "persistent_site_production_manifest_v10_0_5_18_3_5"
        )
        physics = dict(payload.get("physics_contract", {}) or {})
        physics.update(fields)
        payload["physics_contract"] = physics
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

    fields = _quality_fields()
    corridor_path = out / "compact_corridor_mesh_v91852.json"
    if corridor_path.is_file():
        corridor = dict(json.loads(corridor_path.read_text()) or {})
        corridor.update(fields)
        corridor_path.write_text(
            json.dumps(corridor, indent=2, sort_keys=True, default=str) + "\n"
        )

    for old_name, new_name in _OLD_FILES.items():
        old_path = out / old_name
        if not old_path.is_file():
            continue
        payload = _rewrite_payload(old_name, json.loads(old_path.read_text()))
        new_path = out / new_name
        new_path.write_text(
            json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n"
        )
        if new_path != old_path:
            old_path.unlink()

    quality = audit_payload()
    quality.update(fields)
    quality.update(
        {
            "point_release": POINT_RELEASE,
            "model": MODEL_ID,
            "run_completed_without_exception": error is None,
            "runtime_error_type": None if error is None else type(error).__name__,
            "runtime_error": None if error is None else str(error),
        }
    )
    (out / QUALITY_AUDIT).write_text(
        json.dumps(quality, indent=2, sort_keys=True, default=str) + "\n"
    )


def main(argv: list[str] | None = None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    out = _out_path(user_args)
    saved_wrapper = _v91856._strict_quality_advance_v91856
    install_ray_handoff()
    configure_legacy_quality_wrapper(saved_wrapper)
    reset_audit()
    _v91856._strict_quality_advance_v91856 = (
        quality_aware_strict_advance_v10051835
    )
    error: BaseException | None = None
    try:
        return _base.main(user_args)
    except BaseException as exc:
        error = exc
        raise
    finally:
        _v91856._strict_quality_advance_v91856 = saved_wrapper
        _rewrite_outputs(out, error)


if __name__ == "__main__":
    main()


__all__ = [
    "MODEL_ID",
    "POINT_RELEASE",
    "PRODUCTION_MANIFEST",
    "QUALITY_AUDIT",
    "SEED_MANIFEST",
    "SELECTION_MANIFEST",
    "TRANSFER_MANIFEST",
    "main",
]
