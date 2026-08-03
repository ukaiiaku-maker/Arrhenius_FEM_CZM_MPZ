"""Audit-only metadata normalization for the theta-zero full-field parity run.

Several inherited compatibility summaries retain historical labels such as a
45-degree orientation, tip-only bulk state, or a one-second Taylor renewal
window even when the live v10.0.5.18.4.0 overlay used theta=0, evolved the
full FEM bulk, and installed the PF v10.4.1 one-nanosecond closure. This module
normalizes those labels after execution and records every replaced value.

No physical array, FEM state, cohesive state, hazard stream, event geometry,
or constitutive parameter is modified.
"""
from __future__ import annotations

from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any

MODEL_ID = "theta0_pf_full_field_metadata_normalization_v10_0_5_18_4_0"
AUDIT_FILE = "theta0_pf_full_field_metadata_normalization_v10_0_5_18_4_0.json"
SEMANTIC_BULK_MODE = "full_field"
SOLVER_BULK_MODE = "bulk_same_pt_km"
BULK_MODEL = "v10.4.1_bulk_peierls_taylor_detailed_balance_exact_port"


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _load_json(path: Path) -> Any:
    try:
        return json.loads(path.read_text())
    except Exception:
        return None


def _write_json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True, default=str) + "\n")


def _replace(
    mapping: dict[str, Any],
    key: str,
    value: Any,
    changes: list[dict[str, Any]],
    location: str,
) -> None:
    previous = mapping.get(key)
    if previous == value:
        return
    mapping[key] = value
    changes.append(
        {
            "location": location,
            "field": key,
            "previous": previous,
            "authoritative": value,
        }
    )


def _normalize_mapping(
    mapping: dict[str, Any],
    changes: list[dict[str, Any]],
    location: str,
) -> None:
    _replace(mapping, "crystal_theta_deg", 0.0, changes, location)
    _replace(mapping, "material_class", "peak", changes, location)
    _replace(mapping, "target_class", "peak", changes, location)
    _replace(mapping, "bulk_state_evolves_in_fem", True, changes, location)
    _replace(mapping, "bulk_plasticity_mode", SOLVER_BULK_MODE, changes, location)
    _replace(
        mapping,
        "bulk_plasticity_semantic_mode",
        SEMANTIC_BULK_MODE,
        changes,
        location,
    )
    _replace(mapping, "bulk_model", BULK_MODEL, changes, location)
    _replace(mapping, "taylor_renewal_time_s", 1.0e-9, changes, location)
    _replace(mapping, "metadata_normalization_model", MODEL_ID, changes, location)


def normalize_theta0_full_field_outputs(out: Path | str) -> dict[str, Any]:
    root = Path(out).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    changes: list[dict[str, Any]] = []
    files_updated: list[str] = []

    summary_path = root / "summary.json"
    summary = _load_json(summary_path)
    if isinstance(summary, list):
        for index, row in enumerate(summary):
            if isinstance(row, dict):
                _normalize_mapping(row, changes, f"summary.json[{index}]")
        _write_json(summary_path, summary)
        files_updated.append(str(summary_path))
    elif isinstance(summary, dict):
        _normalize_mapping(summary, changes, "summary.json")
        _write_json(summary_path, summary)
        files_updated.append(str(summary_path))

    manifest_names = (
        "persistent_site_production_manifest_v10_0_5_18_3_9.json",
        "persistent_site_production_manifest_v10_0_5_16.json",
        "persistent_site_production_manifest_v10_0_5_18.json",
    )
    for name in manifest_names:
        path = root / name
        payload = _load_json(path)
        if not isinstance(payload, dict):
            continue
        _normalize_mapping(payload, changes, name)
        physics = payload.get("physics_contract")
        if not isinstance(physics, dict):
            physics = {}
            payload["physics_contract"] = physics
        _normalize_mapping(physics, changes, f"{name}.physics_contract")
        policy = payload.get("two_d_state_policy")
        if isinstance(policy, dict):
            _normalize_mapping(policy, changes, f"{name}.two_d_state_policy")
        payload["metadata_normalization_audit_file"] = AUDIT_FILE
        _write_json(path, payload)
        files_updated.append(str(path))

    for path in sorted(root.glob("persistent_site_parameter_selection*.json")):
        payload = _load_json(path)
        if not isinstance(payload, dict):
            continue
        policy = payload.get("policy")
        if not isinstance(policy, dict):
            policy = {}
            payload["policy"] = policy
        _normalize_mapping(policy, changes, f"{path.name}.policy")
        payload["metadata_normalization_audit_file"] = AUDIT_FILE
        _write_json(path, payload)
        files_updated.append(str(path))

    run_args_path = root / "run_args.json"
    run_args = _load_json(run_args_path)
    if isinstance(run_args, dict):
        _replace(run_args, "crystal_theta_deg", 0.0, changes, "run_args.json")
        _replace(
            run_args,
            "bulk_plasticity_mode",
            SOLVER_BULK_MODE,
            changes,
            "run_args.json",
        )
        _replace(
            run_args,
            "bulk_plasticity_semantic_mode",
            SEMANTIC_BULK_MODE,
            changes,
            "run_args.json",
        )
        _replace(run_args, "bulk_model", BULK_MODEL, changes, "run_args.json")
        _write_json(run_args_path, run_args)
        files_updated.append(str(run_args_path))

    audit = {
        "schema": MODEL_ID,
        "created_or_updated_utc": _utc_now(),
        "output_root": str(root),
        "metadata_only": True,
        "physics_or_mesh_modified": False,
        "FEM_state_modified": False,
        "cohesive_state_modified": False,
        "hazard_or_event_geometry_modified": False,
        "authoritative_contract": {
            "crystal_theta_deg": 0.0,
            "material_class": "peak",
            "bulk_plasticity_semantic_mode": SEMANTIC_BULK_MODE,
            "inherited_solver_bulk_mode_token": SOLVER_BULK_MODE,
            "bulk_state_evolves_in_fem": True,
            "bulk_model": BULK_MODEL,
            "taylor_renewal_time_s": 1.0e-9,
        },
        "files_updated": sorted(set(files_updated)),
        "changes": changes,
    }
    _write_json(root / AUDIT_FILE, audit)
    return audit


__all__ = [
    "AUDIT_FILE",
    "BULK_MODEL",
    "MODEL_ID",
    "SEMANTIC_BULK_MODE",
    "SOLVER_BULK_MODE",
    "normalize_theta0_full_field_outputs",
]
