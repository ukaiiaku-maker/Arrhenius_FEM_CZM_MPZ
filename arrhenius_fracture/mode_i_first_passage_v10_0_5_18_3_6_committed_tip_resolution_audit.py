"""v10.0.5.18.3.6 telemetry-only committed-tip resolution audit entry."""
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
from .committed_tip_resolution_audit_v10051836 import (
    AuditedRetryingAdaptiveCZMBackendV10051836,
    MODEL_ID as AUDIT_MODEL,
    SCHEMA as AUDIT_SCHEMA,
    audit_payload,
    configure,
    event_length_diagnostics,
    reset_audit,
)


POINT_RELEASE = "10.0.5.18.3.6"
MODEL_ID = (
    "FEM_CZM_four_class_path_aware_committed_tip_resolution_audit_"
    "v10_0_5_18_3_6"
)
PRODUCTION_MANIFEST = "persistent_site_production_manifest_v10_0_5_18_3_6.json"
SELECTION_MANIFEST = "persistent_site_parameter_selection_v10_0_5_18_3_6.json"
TRANSFER_MANIFEST = "four_class_parameter_transfer_v10_0_5_18_3_6.json"
SEED_MANIFEST = "stochastic_seed_manifest_v10_0_5_18_3_6.json"
AUDIT_FILE = "committed_tip_resolution_audit_v10_0_5_18_3_6.json"

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
        raise SystemExit(f"v10.0.5.18.3.6 requires finite {name}")
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


def _load_geometry_events(path: Path) -> tuple[list[dict[str, Any]], str]:
    payload = json.loads(path.read_text())
    if isinstance(payload, list):
        return [dict(row) for row in payload], "list"
    if isinstance(payload, dict):
        for key in ("events", "geometry_events", "records"):
            if isinstance(payload.get(key), list):
                return [dict(row) for row in payload[key]], key
    return [], "unknown"


def _write_geometry_events(path: Path, rows: list[dict[str, Any]], container: str) -> None:
    if container == "list":
        payload: Any = rows
    elif container in {"events", "geometry_events", "records"}:
        original = json.loads(path.read_text())
        original[container] = rows
        payload = original
    else:
        payload = rows
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n")


def _enrich_stochastic_events(out: Path, audit: dict[str, Any]) -> dict[str, Any]:
    path = out / "stochastic_geometry_events_v10_0_5_16.json"
    if not path.is_file():
        return {
            "available": False,
            "reason": "stochastic_geometry_events_file_missing",
        }
    try:
        rows, container = _load_geometry_events(path)
    except Exception as exc:
        return {
            "available": False,
            "reason": f"{type(exc).__name__}: {exc}",
        }

    clip_counts = {"lower": 0, "none": 0, "upper": 0, "unknown": 0}
    for row in rows:
        threshold = row.get("threshold_action")
        nominal = row.get("requested_fixed_length_m")
        if threshold is None or nominal is None:
            row["v10051836_event_length_diagnostics_available"] = False
            clip_counts["unknown"] += 1
            continue
        diag = event_length_diagnostics(float(threshold), float(nominal))
        row.update(diag)
        row["v10051836_event_length_diagnostics_available"] = True
        state = str(diag["event_length_clip_state"])
        clip_counts[state] = clip_counts.get(state, 0) + 1
    _write_geometry_events(path, rows, container)

    telemetry = list(audit.get("events", []))
    for index, record in enumerate(telemetry):
        if index < len(rows):
            geometry = rows[index]
            record["stochastic_geometry_record_index"] = int(index)
            record["stochastic_event"] = {
                key: geometry.get(key)
                for key in (
                    "threshold_action",
                    "event_length_factor",
                    "requested_fixed_length_m",
                    "requested_stochastic_length_m",
                    "event_length_clip_state",
                    "unclipped_event_length_factor",
                    "realized_event_length_factor",
                    "inserted",
                    "moved_m",
                    "reason",
                )
            }
    audit["events"] = telemetry
    return {
        "available": True,
        "record_count": int(len(rows)),
        "telemetry_record_count": int(len(telemetry)),
        "record_counts_match": bool(len(rows) == len(telemetry)),
        "event_length_clip_counts": clip_counts,
        "event_length_distribution": (
            "continuous threshold-correlated map with atoms created by lower/upper clipping"
        ),
        "event_length_depends_on_mesh_resolution": False,
    }


def _point_release_fields() -> dict[str, Any]:
    return {
        "committed_tip_resolution_audit_active": True,
        "committed_tip_resolution_audit_schema": AUDIT_SCHEMA,
        "committed_tip_resolution_audit_model": AUDIT_MODEL,
        "telemetry_only": True,
        "geometry_changed": False,
        "hazard_or_event_length_changed": False,
        "quality_floors_changed": False,
        "partition_retry_changed": False,
        "constitutive_physics_changed": False,
        "committed_tip_and_candidate_endpoint_metrics_separated": True,
        "joint_triangle_quality_and_child_area_margins_recorded": True,
        "event_length_clip_state_recorded": True,
        "tip_h_fine_contract_recorded": True,
    }


def _rewrite_payload(old_name: str, payload: dict[str, Any], audit: dict[str, Any]) -> dict[str, Any]:
    payload["point_release"] = POINT_RELEASE
    payload["model"] = MODEL_ID
    fields = _point_release_fields()
    if old_name == _base.PRODUCTION_MANIFEST:
        payload["schema"] = "persistent_site_production_manifest_v10_0_5_18_3_6"
        physics = dict(payload.get("physics_contract", {}) or {})
        physics.update(fields)
        payload["physics_contract"] = physics
        payload["committed_tip_resolution_audit"] = {
            "schema": AUDIT_SCHEMA,
            "event_count": audit.get("event_count", 0),
            "joint_margin_class_counts": audit.get(
                "joint_margin_class_counts", {}
            ),
            "tip_h_fine_contract": audit.get("tip_h_fine_contract", {}),
            "stochastic_event_enrichment": audit.get(
                "stochastic_event_enrichment", {}
            ),
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
    audit["stochastic_event_enrichment"] = _enrich_stochastic_events(out, audit)
    audit.update(
        {
            "point_release": POINT_RELEASE,
            "model": MODEL_ID,
            "run_completed_without_exception": error is None,
            "runtime_error_type": None if error is None else type(error).__name__,
            "runtime_error": None if error is None else str(error),
        }
    )
    (out / AUDIT_FILE).write_text(
        json.dumps(audit, indent=2, sort_keys=True, default=str) + "\n"
    )

    corridor_path = out / "compact_corridor_mesh_v91852.json"
    if corridor_path.is_file():
        corridor = dict(json.loads(corridor_path.read_text()) or {})
        corridor.update(_point_release_fields())
        corridor["tip_h_fine_contract"] = audit.get("tip_h_fine_contract", {})
        corridor_path.write_text(
            json.dumps(corridor, indent=2, sort_keys=True, default=str) + "\n"
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
    out = _out_path(user_args)
    tip_h_fine = _float_option(user_args, "--tip-h-fine", 0.0)
    tip_ratio = _float_option(user_args, "--tip-ratio", 1.15)
    if tip_h_fine <= 0.0:
        raise SystemExit("v10.0.5.18.3.6 requires positive --tip-h-fine")
    configure(
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
    reset_audit()

    saved_backend = _ramp.RetryingAdaptiveCZMBackendV1005183
    _ramp.RetryingAdaptiveCZMBackendV1005183 = (
        AuditedRetryingAdaptiveCZMBackendV10051836
    )
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
    "AUDIT_FILE",
    "MODEL_ID",
    "POINT_RELEASE",
    "PRODUCTION_MANIFEST",
    "SEED_MANIFEST",
    "SELECTION_MANIFEST",
    "TRANSFER_MANIFEST",
    "main",
]
