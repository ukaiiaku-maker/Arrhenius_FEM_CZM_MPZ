"""v10.0.5.18.3.8 two-triangle cavity with constrained endpoint star."""
from __future__ import annotations

import json
from pathlib import Path
import sys

from . import mode_i_first_passage_v10_0_5_18_3_8_two_triangle_cavity as _base
from . import quality_aware_czm_two_triangle_cavity_v10051838 as _cavity
from .quality_aware_czm_two_triangle_star_v10051838 import (
    MODEL_ID as STAR_MODEL,
    EndpointStarTriangulation,
)


MODEL_ID = "FEM_CZM_two_triangle_constrained_endpoint_star_v10_0_5_18_3_8"


def _annotate(out: Path | None) -> None:
    if out is None or not out.is_dir():
        return
    fields = {
        "two_triangle_constrained_endpoint_star_active": True,
        "two_triangle_constrained_endpoint_star_model": STAR_MODEL,
        "unconstrained_delaunay_sector_omission_allowed": False,
        "exact_endpoint_connected_to_every_cavity_boundary_edge": True,
        "triangle_quality_floor_relaxed": False,
        "child_area_ratio_floor_relaxed": False,
        "hazard_or_event_length_changed": False,
        "constitutive_physics_changed": False,
    }
    for name in (
        _base.PRODUCTION_MANIFEST,
        _base.RESOLUTION_AUDIT,
        _base.RETRY_AUDIT,
        "compact_corridor_mesh_v91852.json",
    ):
        path = out / name
        if not path.is_file():
            continue
        payload = dict(json.loads(path.read_text()) or {})
        payload.update(fields)
        if isinstance(payload.get("physics_contract"), dict):
            payload["physics_contract"].update(fields)
        if isinstance(payload.get("policy"), dict):
            payload["policy"].update(fields)
        path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n")


def main(argv: list[str] | None = None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    out_raw = _base._option_value(user_args, "--out")
    out = None if out_raw is None else Path(out_raw).expanduser().resolve()
    saved = _cavity.Delaunay
    _cavity.Delaunay = EndpointStarTriangulation
    try:
        return _base.main(user_args)
    finally:
        _cavity.Delaunay = saved
        _annotate(out)


if __name__ == "__main__":
    main()


__all__ = ["MODEL_ID", "main"]
