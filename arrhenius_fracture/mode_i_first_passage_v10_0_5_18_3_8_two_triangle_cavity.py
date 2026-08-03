"""v10.0.5.18.3.8 exact-endpoint two-triangle cavity entry."""
from __future__ import annotations

import json
from pathlib import Path
import sys
from typing import Any

from . import event_resolved_czm_retry_v10051837 as _retry
from . import mode_i_first_passage_v10_0_5_18_3_7_event_resolved_edge_connected as _edge37
from . import mode_i_first_passage_v10_0_5_18_3_7_event_resolved_local_cavity as _local37
from .quality_aware_czm_two_triangle_cavity_v10051838 import (
    MODEL_ID as CAVITY_MODEL,
    install as install_two_triangle_cavity,
)

POINT_RELEASE = "10.0.5.18.3.8"
MODEL_ID = "FEM_CZM_event_resolved_exact_endpoint_two_triangle_cavity_v10_0_5_18_3_8"
PRODUCTION_MANIFEST = "persistent_site_production_manifest_v10_0_5_18_3_8.json"
SELECTION_MANIFEST = "persistent_site_parameter_selection_v10_0_5_18_3_8.json"
TRANSFER_MANIFEST = "four_class_parameter_transfer_v10_0_5_18_3_8.json"
SEED_MANIFEST = "stochastic_seed_manifest_v10_0_5_18_3_8.json"
RESOLUTION_AUDIT = "committed_tip_resolution_audit_v10_0_5_18_3_8.json"
RETRY_AUDIT = "event_resolved_czm_retry_v10_0_5_18_3_8.json"

_RENAMES = {
    _local37.PRODUCTION_MANIFEST: PRODUCTION_MANIFEST,
    _local37.SELECTION_MANIFEST: SELECTION_MANIFEST,
    _local37.TRANSFER_MANIFEST: TRANSFER_MANIFEST,
    _local37.SEED_MANIFEST: SEED_MANIFEST,
    _local37.RESOLUTION_AUDIT: RESOLUTION_AUDIT,
    _local37.RETRY_AUDIT: RETRY_AUDIT,
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
        "two_triangle_exact_endpoint_cavity_active": True,
        "two_triangle_exact_endpoint_cavity_model": CAVITY_MODEL,
        "two_triangle_cavity_boundary_preserved": True,
        "old_near_diagonal_removed": True,
        "exact_endpoint_preinserted_as_bulk_node": True,
        "exact_endpoint_zero_motion_backend_commit": True,
        "two_triangle_cavity_requires_convex_four_vertex_union": True,
        "two_triangle_cavity_final_quality_gate": True,
        "two_triangle_cavity_parent_transfer_policy": "centroid_classified_old_parent_elastic_bulk_only",
        "continuum_bulk_history_transfer_required": False,
        "direct_endpoint_refinement_requires_tip_fan_shared_edge": True,
        "one_node_target_cavity_routes_to_exact_ray_crossing": True,
        "equal_length_partition_retry": False,
        "triangle_quality_floor_relaxed": False,
        "child_area_ratio_floor_relaxed": False,
        "hazard_or_event_length_changed": False,
        "constitutive_physics_changed": False,
    }


def _annotate(payload: dict[str, Any], error: BaseException | None) -> dict[str, Any]:
    fields = _fields()
    payload.update(fields)
    payload["run_completed_without_exception"] = error is None
    payload["runtime_error_type"] = None if error is None else type(error).__name__
    payload["runtime_error"] = None if error is None else str(error)
    if isinstance(payload.get("physics_contract"), dict):
        payload["physics_contract"].update(fields)
    if isinstance(payload.get("policy"), dict):
        payload["policy"].update(fields)
    schema = payload.get("schema")
    if isinstance(schema, str):
        payload["schema"] = schema.replace("v10_0_5_18_3_7", "v10_0_5_18_3_8").replace(
            "v10.0.5.18.3.7", "v10.0.5.18.3.8"
        )
    return payload


def _rewrite_outputs(out: Path | None, error: BaseException | None) -> None:
    if out is None or not out.is_dir():
        return
    for old_name, new_name in _RENAMES.items():
        old_path = out / old_name
        if not old_path.is_file():
            continue
        payload = _annotate(dict(json.loads(old_path.read_text()) or {}), error)
        if new_name == RETRY_AUDIT:
            payload["two_triangle_cavity_after_single_edge_exhaustion"] = True
        new_path = out / new_name
        new_path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n")
        if old_path != new_path:
            old_path.unlink()
    corridor = out / "compact_corridor_mesh_v91852.json"
    if corridor.is_file():
        payload = _annotate(dict(json.loads(corridor.read_text()) or {}), error)
        corridor.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n")


def main(argv: list[str] | None = None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    out_raw = _option_value(user_args, "--out")
    out = None if out_raw is None else Path(out_raw).expanduser().resolve()
    saved_install = _local37.install_ray_handoff
    saved_refine = _retry._local._refine_once
    _local37.install_ray_handoff = install_two_triangle_cavity
    error: BaseException | None = None
    try:
        return _edge37.main(user_args)
    except BaseException as exc:
        error = exc
        raise
    finally:
        _local37.install_ray_handoff = saved_install
        _retry._local._refine_once = saved_refine
        _rewrite_outputs(out, error)


if __name__ == "__main__":
    main()


__all__ = [
    "MODEL_ID", "POINT_RELEASE", "PRODUCTION_MANIFEST", "RESOLUTION_AUDIT",
    "RETRY_AUDIT", "SEED_MANIFEST", "SELECTION_MANIFEST", "TRANSFER_MANIFEST", "main"
]
