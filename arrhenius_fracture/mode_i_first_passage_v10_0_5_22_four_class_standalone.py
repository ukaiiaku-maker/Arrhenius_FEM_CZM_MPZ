"""Standalone four-class FEM/CZM with balanced exact-endpoint prerefinement."""
from __future__ import annotations

import json
from pathlib import Path
import sys

from . import mode_i_first_passage_v10_0_5_20_four_class_standalone as _base20
from . import mode_i_first_passage_v10_0_5_21_four_class_standalone as _base21
from .adaptive_czm_forward_support_prerefine_v100521 import (
    MODEL_ID as FORWARD_SUPPORT_MODEL,
    installed_forward_support_prerefine_v100521,
)
from .adaptive_czm_balanced_support_prerefine_v100522 import (
    MODEL_ID as BALANCED_SUPPORT_MODEL,
    installed_balanced_support_prerefine_v100522,
)

POINT_RELEASE = "10.0.5.22-frozen-four-class-balanced-support-prerefine"
MODEL_ID = "FEM_CZM_full_2D_frozen_four_class_balanced_support_prerefine_v10_0_5_22"
PRODUCTION_MANIFEST = _base20.PRODUCTION_MANIFEST
SELECTION_MANIFEST = _base20.SELECTION_MANIFEST
TRANSFER_MANIFEST = _base20.TRANSFER_MANIFEST


def _out_path(argv: list[str]) -> Path | None:
    if "--out" not in argv:
        return None
    index = argv.index("--out")
    if index + 1 >= len(argv):
        return None
    return Path(argv[index + 1]).expanduser().resolve()


def _rewrite_metadata(out: Path | None) -> None:
    _base21._rewrite_metadata(out)
    if out is None:
        return
    manifest = out / PRODUCTION_MANIFEST
    if manifest.is_file():
        payload = json.loads(manifest.read_text())
        payload["schema"] = "persistent_site_production_manifest_v10_0_5_22_frozen_four_class"
        payload["model"] = MODEL_ID
        payload["point_release"] = POINT_RELEASE
        physics = dict(payload.get("physics_contract", {}) or {})
        physics.update({
            "adaptive_czm_forward_support_prerefine_model": FORWARD_SUPPORT_MODEL,
            "adaptive_czm_balanced_support_prerefine_model": BALANCED_SUPPORT_MODEL,
            "exact_crack_endpoint_preserved": True,
            "exact_crack_direction_preserved": True,
            "exact_stochastic_event_length_preserved": True,
            "intact_off_ray_support_node_allowed": True,
            "triangle_quality_floor_relaxed": False,
            "child_area_ratio_floor_relaxed": False,
            "constitutive_physics_changed_by_geometry_repair": False,
            "hazard_law_changed_by_geometry_repair": False,
        })
        payload["physics_contract"] = physics
        manifest.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n")


def main(argv: list[str] | None = None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    out = _out_path(user_args)
    try:
        with installed_forward_support_prerefine_v100521():
            with installed_balanced_support_prerefine_v100522():
                return _base20.main(user_args)
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
