"""v10.0.5.11 production J probe with same-mesh fixed-grip energy release."""
from __future__ import annotations

import inspect
from typing import Any

from . import mode_i_first_passage_v10_0_5_10_refinement_probe as _v100510
from .production_j_same_mesh_energy_v100511 import PROBE_JSON

POINT_RELEASE = "10.0.5.11"
MODEL_ID = "FEM_CZM_production_J_same_mesh_fixed_grip_energy_probe_v10_0_5_11"

_BASE_IMPORT = (
    "from arrhenius_fracture.production_j_refinement_support_v100510 "
    "import record_production_j_refinement_probe_v100510"
)
_NEW_IMPORT = (
    "from arrhenius_fracture.production_j_same_mesh_energy_v100511 "
    "import record_production_j_same_mesh_probe_v100511"
)
_BASE_CALL = "record_production_j_refinement_probe_v100510("
_NEW_CALL = "record_production_j_same_mesh_probe_v100511("
_CALL_TAIL = """                        crystal_theta_deg=float(getattr(args, 'crystal_theta_deg', 0.0) or 0.0),
                    )
"""
_CALL_TAIL_NEW = """                        crystal_theta_deg=float(getattr(args, 'crystal_theta_deg', 0.0) or 0.0),
                        boundary_data=bnd,
                        total_grip_opening_m=Uapp,
                    )
"""


def _replace_unique(source: str, old: str, new: str, name: str) -> str:
    count = source.count(old)
    if count != 1:
        raise RuntimeError(f"v10.0.5.11 expected one {name}; found {count}")
    return source.replace(old, new)


def patch_run_2d_source_v100511(source: str) -> str:
    patched = _v100510.patch_run_2d_source_v100510(source)
    patched = _replace_unique(patched, _BASE_IMPORT, _NEW_IMPORT, "recorder import")
    patched = _replace_unique(patched, _BASE_CALL, _NEW_CALL, "recorder call")
    patched = _replace_unique(patched, _CALL_TAIL, _CALL_TAIL_NEW, "same-mesh arguments")
    return patched


def validate_source_transform_v100511() -> dict[str, Any]:
    source = inspect.getsource(
        __import__("arrhenius_fracture.sharp_front", fromlist=["run_2d"]).run_2d
    )
    patched = patch_run_2d_source_v100511(source)
    compile(patched, "<v10.0.5.11-same-mesh-probe>", "exec")
    checks = {
        "same_mesh_recorder": "record_production_j_same_mesh_probe_v100511(" in patched,
        "boundary_data_supplied": "boundary_data=bnd" in patched,
        "fixed_grip_opening_supplied": "total_grip_opening_m=Uapp" in patched,
        "v10_0_5_10_contour_recorder_composition_preserved": "contours_v10059" in patched,
        "v10_0_5_9_production_path_preserved": "straight_progressive_cluster_no_exclusion" in patched,
        "full_audited_v10055_stack": "cohesive_elements" in patched,
    }
    failed = [name for name, passed in checks.items() if not passed]
    if failed:
        raise RuntimeError("v10.0.5.11 source-transform preflight failed: " + ", ".join(failed))
    return {
        "point_release": POINT_RELEASE,
        "model": MODEL_ID,
        "source_transform_preflight_passed": True,
        "constitutive_physics_changed": False,
        **checks,
    }


def main(argv: list[str] | None = None):
    saved_patch = _v100510.patch_run_2d_source_v100510
    saved_validate = _v100510.validate_source_transform_v100510
    saved_probe_json = _v100510.PROBE_JSON
    _v100510.patch_run_2d_source_v100510 = patch_run_2d_source_v100511
    _v100510.validate_source_transform_v100510 = validate_source_transform_v100511
    _v100510.PROBE_JSON = PROBE_JSON
    try:
        return _v100510.main(argv)
    finally:
        _v100510.patch_run_2d_source_v100510 = saved_patch
        _v100510.validate_source_transform_v100510 = saved_validate
        _v100510.PROBE_JSON = saved_probe_json


if __name__ == "__main__":
    main()


__all__ = [
    "POINT_RELEASE",
    "MODEL_ID",
    "patch_run_2d_source_v100511",
    "validate_source_transform_v100511",
    "main",
]
