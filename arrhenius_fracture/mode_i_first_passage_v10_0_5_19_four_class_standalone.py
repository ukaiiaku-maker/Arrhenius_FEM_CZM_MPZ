"""Standalone frozen four-class FEM/CZM with adaptive quality partitioning."""
from __future__ import annotations

import json
from pathlib import Path
import sys

from . import mode_i_first_passage_v10_0_5_18_four_class_standalone as _base
from .adaptive_czm_quality_partition_v100519 import (
    MODEL_ID as ADAPTIVE_PARTITION_MODEL,
    installed_adaptive_quality_partition_v100519,
)

POINT_RELEASE = "10.0.5.19-frozen-four-class-adaptive-quality-partition"
MODEL_ID = "FEM_CZM_full_2D_frozen_four_class_stochastic_adaptive_quality_partition_v10_0_5_19"
PRODUCTION_MANIFEST = _base.PRODUCTION_MANIFEST
SELECTION_MANIFEST = _base.SELECTION_MANIFEST
TRANSFER_MANIFEST = _base.TRANSFER_MANIFEST


def _out_path(argv: list[str]) -> Path | None:
    if "--out" not in argv:
        return None
    index = argv.index("--out")
    if index + 1 >= len(argv):
        return None
    return Path(argv[index + 1]).expanduser().resolve()


def _rewrite_metadata(out: Path | None) -> None:
    if out is None:
        return
    manifest = out / PRODUCTION_MANIFEST
    if manifest.is_file():
        payload = json.loads(manifest.read_text())
        payload["schema"] = "persistent_site_production_manifest_v10_0_5_19_frozen_four_class"
        payload["model"] = MODEL_ID
        payload["point_release"] = POINT_RELEASE
        physics = dict(payload.get("physics_contract", {}) or {})
        physics.update(
            {
                "adaptive_czm_quality_partition_model": ADAPTIVE_PARTITION_MODEL,
                "adaptive_czm_recursive_segment_bisection_active": True,
                "adaptive_czm_uniform_subdivision_retained_as_fast_path": True,
                "adaptive_czm_exact_event_endpoint_preserved": True,
                "adaptive_czm_exact_event_length_preserved": True,
                "adaptive_czm_physical_event_atomic": True,
                "recursive_geometry_veto_counted_once_per_physical_event": True,
                "triangle_quality_floor_relaxed": False,
                "child_area_ratio_floor_relaxed": False,
                "constitutive_physics_changed_by_geometry_repair": False,
                "stochastic_cleavage_law_changed_by_geometry_repair": False,
                "event_length_changed_by_geometry_repair": False,
            }
        )
        payload["physics_contract"] = physics
        manifest.write_text(
            json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n"
        )

    selection = out / SELECTION_MANIFEST
    if selection.is_file():
        payload = json.loads(selection.read_text())
        payload["schema"] = MODEL_ID
        payload["point_release"] = POINT_RELEASE
        policy = dict(payload.get("policy", {}) or {})
        policy.update(
            {
                "adaptive_czm_quality_partition_model": ADAPTIVE_PARTITION_MODEL,
                "quality_floors_unchanged": True,
                "physical_event_endpoint_preserved": True,
                "physical_event_length_preserved": True,
                "physical_event_atomic": True,
            }
        )
        payload["policy"] = policy
        selection.write_text(
            json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n"
        )


def main(argv: list[str] | None = None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    out = _out_path(user_args)
    try:
        with installed_adaptive_quality_partition_v100519():
            return _base.main(user_args)
    finally:
        _rewrite_metadata(out)


if __name__ == "__main__":
    main()


__all__ = [
    "MODEL_ID",
    "POINT_RELEASE",
    "PRODUCTION_MANIFEST",
    "SELECTION_MANIFEST",
    "TRANSFER_MANIFEST",
    "main",
]
