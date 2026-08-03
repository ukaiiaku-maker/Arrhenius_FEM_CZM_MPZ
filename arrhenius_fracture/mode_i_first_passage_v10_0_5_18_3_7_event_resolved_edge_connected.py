"""Edge-connected endpoint classification overlay for v10.0.5.18.3.7.

A target-containing triangle that shares only one node with the current tip fan
is not treated as a local endpoint cavity.  The exact physical segment is first
marched to its intervening mesh-ray crossing.  Direct endpoint refinement is
reserved for the current tip triangle or an adjacent triangle sharing a full
edge (two fan nodes).

The overlay also records committed-tip telemetry at the strict transaction
layer.  This avoids losing the audit when the global production wrapper replaces
the backend subclass ``advance`` method.

No hazard, event-length, direction, quality floor, material, MPZ, source,
shielding, or constitutive law is changed.
"""
from __future__ import annotations

import copy
import json
from pathlib import Path
import sys
from typing import Any

import numpy as np

from . import committed_tip_resolution_audit_v10051836 as _audit36
from . import event_resolved_czm_retry_v10051837 as _retry
from . import mode_i_first_passage_v10_0_5_18_3_7_event_resolved_local_cavity as _base
from . import quality_aware_czm_patch_v10051835 as _patch
from . import quality_aware_czm_ray_handoff_v10051835 as _ray


MODEL_ID = (
    "FEM_CZM_four_class_path_aware_event_resolved_exact_ray_"
    "edge_connected_local_cavity_v10_0_5_18_3_7"
)
CONTROLLER_REVISION = "edge_connected_local_endpoint_v2"


def edge_connected_local_endpoint_reachable(
    self: Any,
    mesh: Any,
    p0,
    p1,
    front_id: int,
) -> bool:
    """Return true only for a tip triangle or a shared-edge target cavity."""
    if _patch._target_tip_triangle(self, mesh, p0, p1, front_id) is not None:
        return True
    adjacent, shared = _ray._adjacent_target_triangle(
        self, mesh, p0, p1, front_id
    )
    return adjacent is not None and int(shared) >= 2


def edge_connected_audited_advance_v10051837(self: Any, *args, **kwargs):
    """Run the v10.0.5.18.3.7 transaction and retain pre-event telemetry."""
    original = getattr(edge_connected_audited_advance_v10051837, "_original", None)
    if original is None:
        raise RuntimeError("edge-connected strict wrapper has no bound raw advance")
    _retry.event_resolved_strict_advance_v10051837._original = original

    nested = int(kwargs.get("_subdepth", 0) or 0) > 0
    if nested:
        return _retry.event_resolved_strict_advance_v10051837(
            self, *args, **kwargs
        )

    mesh = kwargs["mesh"]
    p0 = np.asarray(kwargs["p0"], dtype=float)
    p1 = np.asarray(kwargs["p1"], dtype=float)
    front_id = int(kwargs.get("front_id", 0))
    requested = float(np.linalg.norm(p1 - p0))
    audit_count_before = len(_audit36._AUDIT["events"])
    retry_count_before = len(_retry._AUDIT["events"])

    committed = _audit36._tip_one_ring(self, mesh, p0, front_id)
    record: dict[str, Any] = {
        "physical_event_index": int(audit_count_before),
        "front_id": int(front_id),
        "committed_tip_point_m": p0.tolist(),
        "requested_endpoint_m": p1.tolist(),
        "requested_event_length_m": float(requested),
        "requested_direction": np.asarray(
            kwargs.get("direction", []), dtype=float
        ).tolist(),
        "committed_tip_one_ring": committed,
        "committed_tip_h_over_requested_da": (
            None
            if not committed.get("available") or requested <= 0.0
            else float(committed["h_mean_m"]) / requested
        ),
        "target_parent": _audit36._target_metrics(self, mesh, p1),
        "ray_corridor": _audit36._ray_corridor(self, mesh, p0, p1),
    }

    result = _retry.event_resolved_strict_advance_v10051837(
        self, *args, **kwargs
    )
    endpoint = p1
    if bool(getattr(result, "inserted", False)) and front_id in getattr(
        self, "tip_nodes", {}
    ):
        endpoint = np.asarray(self.tip_nodes[front_id][2], dtype=float)
    record["backend_result"] = _audit36._candidate_quality(
        self, mesh, result, endpoint, front_id
    )
    record["event_resolved_retry"] = (
        copy.deepcopy(_retry._AUDIT["events"][-1])
        if len(_retry._AUDIT["events"]) > retry_count_before
        else None
    )
    record["strict_wrapper_telemetry_active"] = True
    record["direct_endpoint_refinement_requires_tip_fan_shared_edge"] = True

    # The subclass method is normally bypassed by the global production wrapper.
    # Guard against duplicate records if a future call path invokes both layers.
    if len(_audit36._AUDIT["events"]) == audit_count_before:
        _audit36._AUDIT["events"].append(record)
    return result


def _option_value(argv: list[str], name: str) -> str | None:
    for index, token in enumerate(argv):
        if token == name and index + 1 < len(argv):
            return argv[index + 1]
        if token.startswith(name + "="):
            return token.split("=", 1)[1]
    return None


def _annotate_payload(payload: dict[str, Any]) -> dict[str, Any]:
    fields = {
        "event_resolved_controller_revision": CONTROLLER_REVISION,
        "direct_endpoint_refinement_requires_tip_fan_shared_edge": True,
        "minimum_tip_fan_shared_node_count_for_direct_endpoint": 2,
        "one_node_target_cavity_routes_to_exact_ray_crossing": True,
        "committed_tip_telemetry_recorded_at_strict_wrapper": True,
        "equal_length_partition_retry": False,
        "hazard_or_event_length_changed": False,
        "constitutive_physics_changed": False,
    }
    payload.update(fields)
    if isinstance(payload.get("physics_contract"), dict):
        payload["physics_contract"].update(fields)
    if isinstance(payload.get("policy"), dict):
        payload["policy"].update(fields)
    return payload


def _annotate_outputs(out: Path | None, error: BaseException | None) -> None:
    if out is None or not out.is_dir():
        return
    names = (
        _base.RETRY_AUDIT,
        _base.RESOLUTION_AUDIT,
        _base.PRODUCTION_MANIFEST,
        _base.SELECTION_MANIFEST,
        _base.TRANSFER_MANIFEST,
        _base.SEED_MANIFEST,
        "compact_corridor_mesh_v91852.json",
    )
    for name in names:
        path = out / name
        if not path.is_file():
            continue
        payload = _annotate_payload(dict(json.loads(path.read_text()) or {}))
        payload["model"] = MODEL_ID
        payload["edge_connected_overlay_completed_without_exception"] = (
            error is None
        )
        path.write_text(
            json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n"
        )


def main(argv: list[str] | None = None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    out_raw = _option_value(user_args, "--out")
    out = None if out_raw is None else Path(out_raw).expanduser().resolve()

    saved_reachable = _retry._local_endpoint_reachable
    saved_wrapper = _base.event_resolved_strict_advance_v10051837
    _retry._local_endpoint_reachable = edge_connected_local_endpoint_reachable
    _base.event_resolved_strict_advance_v10051837 = (
        edge_connected_audited_advance_v10051837
    )
    error: BaseException | None = None
    try:
        return _base.main(user_args)
    except BaseException as exc:
        error = exc
        raise
    finally:
        _base.event_resolved_strict_advance_v10051837 = saved_wrapper
        _retry._local_endpoint_reachable = saved_reachable
        _annotate_outputs(out, error)


if __name__ == "__main__":
    main()


__all__ = [
    "CONTROLLER_REVISION",
    "MODEL_ID",
    "edge_connected_audited_advance_v10051837",
    "edge_connected_local_endpoint_reachable",
    "main",
]
