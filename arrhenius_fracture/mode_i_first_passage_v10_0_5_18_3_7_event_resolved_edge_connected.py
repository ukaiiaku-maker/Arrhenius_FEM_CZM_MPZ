"""Edge-connected endpoint classification overlay for v10.0.5.18.3.7.

A target-containing triangle that shares only one node with the current tip fan
is not treated as a local endpoint cavity.  The exact physical segment is first
marched to its intervening mesh-ray crossing.  Direct endpoint refinement is
reserved for the current tip triangle or an adjacent triangle sharing a full
edge (two fan nodes).

No hazard, event-length, direction, quality floor, material, MPZ, source,
shielding, or constitutive law is changed.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys
from typing import Any

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

    saved = _retry._local_endpoint_reachable
    _retry._local_endpoint_reachable = edge_connected_local_endpoint_reachable
    error: BaseException | None = None
    try:
        return _base.main(user_args)
    except BaseException as exc:
        error = exc
        raise
    finally:
        _retry._local_endpoint_reachable = saved
        _annotate_outputs(out, error)


if __name__ == "__main__":
    main()


__all__ = [
    "CONTROLLER_REVISION",
    "MODEL_ID",
    "edge_connected_local_endpoint_reachable",
    "main",
]
