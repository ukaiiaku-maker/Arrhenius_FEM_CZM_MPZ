"""v10.0.5.18.3.7 exact-ray event-resolved local cavity entry."""
from __future__ import annotations

import json
import math
import os
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
from . import (
    mode_i_first_passage_v10_0_5_18_3_6_committed_tip_resolution_audit
    as _entry36
)
from . import committed_tip_resolution_audit_v10051836 as _audit36
from . import mode_i_first_passage_v9_18_5_6 as _v91856
from .event_resolved_czm_retry_v10051837 import (
    EventResolvedAuditedAdaptiveCZMBackendV10051837,
    MODEL_ID as RETRY_MODEL,
    SCHEMA as RETRY_SCHEMA,
    audit_payload as retry_audit_payload,
    configure_legacy_quality_wrapper,
    event_resolved_strict_advance_v10051837,
    reset_audit as reset_retry_audit,
)
from .quality_aware_czm_ray_handoff_v10051835 import install as install_ray_handoff


POINT_RELEASE = "10.0.5.18.3.7"
MODEL_ID = (
    "FEM_CZM_four_class_path_aware_event_resolved_exact_ray_"
    "local_cavity_v10_0_5_18_3_7"
)
PRODUCTION_MANIFEST = "persistent_site_production_manifest_v10_0_5_18_3_7.json"
SELECTION_MANIFEST = "persistent_site_parameter_selection_v10_0_5_18_3_7.json"
TRANSFER_MANIFEST = "four_class_parameter_transfer_v10_0_5_18_3_7.json"
SEED_MANIFEST = "stochastic_seed_manifest_v10_0_5_18_3_7.json"
RESOLUTION_AUDIT = "committed_tip_resolution_audit_v10_0_5_18_3_7.json"
RETRY_AUDIT = "event_resolved_czm_retry_v10_0_5_18_3_7.json"

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


def _float_option(argv: list[str], name: str, default: float) -> float:
    raw = _option_value(argv, name)
    value = float(default if raw is None else raw)
    if not math.isfinite(value):
        raise SystemExit(f"v10.0.5.18.3.7 requires finite {name}")
    return value


def _float_env(name: str, default: float) -> float:
    try:
        value = float(os.environ.get(name, default))
    except (TypeError, ValueError):
        value = float(default)
    return value if math.isfinite(value) else float(default)


def _out_path(argv: list[str]) -> Path | None:
    raw = _option_value(argv, "--out")
    return None if raw is None else Path(raw).expanduser().resolve()


def _fields() -> dict[str, Any]:
    return {
        "event_resolved_local_cavity_retry_active": True,
        "event_resolved_local_cavity_retry_schema": RETRY_SCHEMA,
        "event_resolved_local_cavity_retry_model": RETRY_MODEL,
        "committed_tip_resolution_audit_active": True,
        "exact_ray_crossing_decomposition": True,
        "equal_length_partition_retry": False,
        "exact_stochastic_event_endpoint_preserved": True,
        "exact_selected_crack_direction_preserved": True,
        "atomic_physical_event_transaction": True,
        "tip_radial_refinement_active": True,
        "adjacent_target_cavity_refinement_active": True,
        "triangle_quality_floor_relaxed": False,
        "child_area_ratio_floor_relaxed": False,
        "tip_h_over_da_role": "refinement_trigger_and_audit_only_not_event_veto",
        "hazard_or_event_length_changed": False,
        "constitutive_physics_changed": False,
    }


def _rewrite_payload(
    old_name: str,
    payload: dict[str, Any],
    resolution: dict[str, Any],
    retry: dict[str, Any],
) -> dict[str, Any]:
    payload["point_release"] = POINT_RELEASE
    payload["model"] = MODEL_ID
    fields = _fields()
    if old_name == _base.PRODUCTION_MANIFEST:
        payload["schema"] = "persistent_site_production_manifest_v10_0_5_18_3_7"
        physics = dict(payload.get("physics_contract", {}) or {})
        physics.update(fields)
        payload["physics_contract"] = physics
        payload["committed_tip_resolution_audit"] = {
            "event_count": resolution.get("event_count", 0),
            "joint_margin_class_counts": resolution.get("joint_margin_class_counts", {}),
            "tip_h_fine_contract": resolution.get("tip_h_fine_contract", {}),
            "stochastic_event_enrichment": resolution.get("stochastic_event_enrichment", {}),
        }
        payload["event_resolved_retry_audit"] = {
            "event_count": retry.get("event_count", 0),
            "success_count": retry.get("success_count", 0),
            "failure_count": retry.get("failure_count", 0),
            "equal_length_partition_retry": False,
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

    resolution = _audit36.audit_payload()
    resolution["stochastic_event_enrichment"] = _entry36._enrich_stochastic_events(
        out, resolution
    )
    resolution.update(
        {
            "point_release": POINT_RELEASE,
            "model": MODEL_ID,
            "run_completed_without_exception": error is None,
            "runtime_error_type": None if error is None else type(error).__name__,
            "runtime_error": None if error is None else str(error),
            **_fields(),
        }
    )
    (out / RESOLUTION_AUDIT).write_text(
        json.dumps(resolution, indent=2, sort_keys=True, default=str) + "\n"
    )

    retry = retry_audit_payload()
    retry.update(
        {
            "point_release": POINT_RELEASE,
            "model": MODEL_ID,
            "run_completed_without_exception": error is None,
            "runtime_error_type": None if error is None else type(error).__name__,
            "runtime_error": None if error is None else str(error),
            **_fields(),
        }
    )
    (out / RETRY_AUDIT).write_text(
        json.dumps(retry, indent=2, sort_keys=True, default=str) + "\n"
    )

    corridor_path = out / "compact_corridor_mesh_v91852.json"
    if corridor_path.is_file():
        corridor = dict(json.loads(corridor_path.read_text()) or {})
        corridor.update(_fields())
        corridor["tip_h_fine_contract"] = resolution.get("tip_h_fine_contract", {})
        corridor_path.write_text(
            json.dumps(corridor, indent=2, sort_keys=True, default=str) + "\n"
        )

    for old_name, new_name in _OLD_FILES.items():
        old_path = out / old_name
        if not old_path.is_file():
            continue
        payload = _rewrite_payload(
            old_name, json.loads(old_path.read_text()), resolution, retry
        )
        new_path = out / new_name
        new_path.write_text(
            json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n"
        )
        if new_path != old_path:
            old_path.unlink()


def main(argv: list[str] | None = None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    out = _out_path(user_args)
    tip_h_fine = _float_option(user_args, "--tip-h-fine", 0.0)
    tip_ratio = _float_option(user_args, "--tip-ratio", 1.15)
    if tip_h_fine <= 0.0:
        raise SystemExit("v10.0.5.18.3.7 requires positive --tip-h-fine")

    _audit36.configure(
        tip_h_fine_m=tip_h_fine,
        tip_ratio=tip_ratio,
        triangle_quality_floor=_float_env(
            "ARRHENIUS_MIN_ACCEPTED_TRIANGLE_QUALITY", 0.035
        ),
        child_area_ratio_floor=_float_env(
            "ARRHENIUS_MIN_ACCEPTED_CHILD_AREA_RATIO", 0.08
        ),
        event_minimum_factor=_float_env("CLEAVAGE_EVENT_MIN_FACTOR", 0.5),
        event_maximum_factor=_float_env("CLEAVAGE_EVENT_MAX_FACTOR", 4.0),
    )
    _audit36.reset_audit()
    reset_retry_audit()

    saved_backend = _ramp.RetryingAdaptiveCZMBackendV1005183
    saved_wrapper = _v91856._strict_quality_advance_v91856
    _ramp.RetryingAdaptiveCZMBackendV1005183 = (
        EventResolvedAuditedAdaptiveCZMBackendV10051837
    )
    configure_legacy_quality_wrapper(saved_wrapper)
    install_ray_handoff()
    _v91856._strict_quality_advance_v91856 = (
        event_resolved_strict_advance_v10051837
    )

    error: BaseException | None = None
    try:
        return _base.main(user_args)
    except BaseException as exc:
        error = exc
        raise
    finally:
        _v91856._strict_quality_advance_v91856 = saved_wrapper
        _ramp.RetryingAdaptiveCZMBackendV1005183 = saved_backend
        _rewrite_outputs(out, error)


if __name__ == "__main__":
    main()


__all__ = [
    "MODEL_ID",
    "POINT_RELEASE",
    "PRODUCTION_MANIFEST",
    "RESOLUTION_AUDIT",
    "RETRY_AUDIT",
    "SEED_MANIFEST",
    "SELECTION_MANIFEST",
    "TRANSFER_MANIFEST",
    "main",
]
