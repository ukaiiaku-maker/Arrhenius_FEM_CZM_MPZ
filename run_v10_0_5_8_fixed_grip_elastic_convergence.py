#!/usr/bin/env python3
"""Run the v10.0.5.8 fixed-grip elastic FEM convergence audit.

The audit runner installs an audit-only mesh constructor that keeps the physical
refined region fixed across ``tip_h`` values and inserts every perturbed crack-tip
position as an explicit node. Production fracture meshes are not changed.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import arrhenius_fracture.fixed_grip_elastic_audit_v10058 as _audit
from arrhenius_fracture.fixed_grip_audit_mesh_v10058 import (
    clear_fixed_grip_audit_mesh,
    configure_fixed_grip_audit_mesh,
    make_fixed_radius_audit_mesh,
)


def _production_args(argv: list[str]) -> list[str]:
    """Supply increments that are resolvable on the default 10 um coarse mesh."""
    out = list(argv)
    if "--crack-increment-um" not in out:
        out.extend(["--crack-increment-um", "40 20 10"])
    return out


def main(argv=None) -> int:
    arg_list = _production_args(list(argv) if argv is not None else sys.argv[1:])
    args = _audit.build_parser().parse_args(arg_list)
    increments = _audit._values(args.crack_increment_um, 1.0e-6)
    contours = _audit._values(args.contour_outer_um, 1.0e-6)
    crack_center = float(args.crack_mm) * 1.0e-3
    width = float(args.width_mm) * 1.0e-3
    height = float(args.height_mm) * 1.0e-3
    nearest = min(crack_center, width - crack_center, 0.5 * height)
    refinement_radius = min(1.10 * max(contours), 0.80 * nearest)
    spec = configure_fixed_grip_audit_mesh(
        crack_center_m=crack_center,
        crack_increments_m=increments,
        refinement_radius_m=refinement_radius,
    )

    original_make_mesh = _audit.make_tri_mesh
    _audit.make_tri_mesh = make_fixed_radius_audit_mesh
    try:
        status = int(_audit.main(arg_list) or 0)
    finally:
        _audit.make_tri_mesh = original_make_mesh
        clear_fixed_grip_audit_mesh()

    summary_path = Path(args.out).resolve() / _audit.SUMMARY_JSON
    if summary_path.exists():
        payload = json.loads(summary_path.read_text())
        payload["audit_mesh"] = {
            "policy": "fixed_physical_refinement_radius_with_explicit_perturbed_tip_nodes",
            "refinement_radius_m": float(spec.refinement_radius_m),
            "crack_tip_nodes_m": [
                float(spec.crack_center_m + offset) for offset in spec.crack_offsets_m
            ],
            "production_mesh_changed": False,
        }
        summary_path.write_text(json.dumps(payload, indent=2, default=str))
        print(
            "AUDIT MESH: fixed refinement radius "
            f"{1.0e6 * spec.refinement_radius_m:.6g} um; "
            f"{len(spec.crack_offsets_m)} explicit crack-tip nodes"
        )
    return status


if __name__ == "__main__":
    raise SystemExit(main())
