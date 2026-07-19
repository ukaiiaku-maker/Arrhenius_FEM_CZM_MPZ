"""One-step v10.0.5.x production J probe with fixed physical mesh support."""
from __future__ import annotations

import inspect
import json
import os
from pathlib import Path
import sys
from typing import Any

from . import mesh as _mesh
from . import mode_i_first_passage_v10_0_5_5_stochastic_vhcf_audited as _v10055
from . import plasticity as _plasticity
from .mode_i_first_passage_v10_0_5_9_production_j_probe import (
    _elastic_update_plasticity,
    _ensure_v911_probe_contract,
    _option_value,
    patch_run_2d_source_v10059,
)
from .physical_refinement_mesh_v100510 import (
    clear_physical_refinement_v100510,
    configure_physical_refinement_v100510,
    make_physical_refinement_mesh_v100510,
)
from .production_j_refinement_support_v100510 import PROBE_JSON

POINT_RELEASE = "10.0.5.10"
MODEL_ID = "FEM_CZM_production_J_fixed_physical_refinement_probe_v10_0_5_10"

_BASE_IMPORT = (
    "from arrhenius_fracture.production_j_parity_v10059 "
    "import record_production_j_probe_v10059"
)
_NEW_IMPORT = (
    "from arrhenius_fracture.production_j_refinement_support_v100510 "
    "import record_production_j_refinement_probe_v100510"
)


def _replace_unique(source: str, old: str, new: str, name: str) -> str:
    count = source.count(old)
    if count != 1:
        raise RuntimeError(f"v10.0.5.10 expected one {name}; found {count}")
    return source.replace(old, new)


def patch_run_2d_source_v100510(source: str) -> str:
    patched = patch_run_2d_source_v10059(source)
    patched = _replace_unique(patched, _BASE_IMPORT, _NEW_IMPORT, "recorder import")
    patched = _replace_unique(
        patched,
        "record_production_j_probe_v10059(",
        "record_production_j_refinement_probe_v100510(",
        "recorder call",
    )
    return patched


def validate_source_transform_v100510() -> dict[str, Any]:
    source = inspect.getsource(
        __import__("arrhenius_fracture.sharp_front", fromlist=["run_2d"]).run_2d
    )
    patched = patch_run_2d_source_v100510(source)
    compile(patched, "<v10.0.5.10-refinement-probe>", "exec")
    checks = {
        "refinement_recorder": "record_production_j_refinement_probe_v100510(" in patched,
        "v10_0_5_9_production_path_preserved": "straight_progressive_cluster_no_exclusion" in patched,
        "full_audited_v10055_stack": "cohesive_elements" in patched,
    }
    failed = [name for name, passed in checks.items() if not passed]
    if failed:
        raise RuntimeError("v10.0.5.10 source-transform preflight failed: " + ", ".join(failed))
    return {
        "point_release": POINT_RELEASE,
        "model": MODEL_ID,
        "source_transform_preflight_passed": True,
        "constitutive_physics_changed": False,
        **checks,
    }


def main(argv: list[str] | None = None):
    args = _ensure_v911_probe_contract(list(sys.argv[1:] if argv is None else argv))
    out_value = _option_value(args, "--out")
    if out_value is None:
        raise SystemExit("v10.0.5.10 refinement probe requires --out")
    out = Path(out_value).resolve()
    out.mkdir(parents=True, exist_ok=True)
    probe_path = out / PROBE_JSON
    if probe_path.exists():
        probe_path.unlink()

    radius_text = os.environ.get("ARRHENIUS_V100510_REFINEMENT_RADIUS_M", "")
    if not radius_text:
        raise SystemExit("v10.0.5.10 requires ARRHENIUS_V100510_REFINEMENT_RADIUS_M")
    radius_m = float(radius_text)
    configure_physical_refinement_v100510(radius_m)
    validate_source_transform_v100510()

    os.environ["ARRHENIUS_V10059_PROBE_PATH"] = str(probe_path)
    os.environ.setdefault("ARRHENIUS_V10059_CONTOURS_UM", "100 140 180 240 300")
    os.environ.setdefault("ARRHENIUS_EVENT_STATISTICS", "mean_field")
    os.environ.setdefault("ARRHENIUS_STOCHASTIC_EMISSION", "0")
    os.environ.setdefault("ARRHENIUS_VHCF_FEM_CACHE", "0")

    saved_patch = _v10055.patch_run_2d_source_v10055_audited
    saved_update = _plasticity.update_plasticity
    saved_mesh = _mesh.make_tri_mesh
    _v10055.patch_run_2d_source_v10055_audited = patch_run_2d_source_v100510
    _plasticity.update_plasticity = _elastic_update_plasticity
    _mesh.make_tri_mesh = make_physical_refinement_mesh_v100510
    try:
        result = _v10055.main(args)
    finally:
        _v10055.patch_run_2d_source_v10055_audited = saved_patch
        _plasticity.update_plasticity = saved_update
        _mesh.make_tri_mesh = saved_mesh
        clear_physical_refinement_v100510()
        os.environ.pop("ARRHENIUS_V10059_PROBE_PATH", None)

    if not probe_path.exists():
        raise RuntimeError("production path completed without writing the v10.0.5.10 probe")
    payload = json.loads(probe_path.read_text())
    payload["source_transform"] = validate_source_transform_v100510()
    payload["fixed_physical_refinement"] = {
        "radius_m": radius_m,
        "radius_um": radius_m * 1.0e6,
        "audit_only": True,
        "production_constitutive_physics_changed": False,
    }
    payload["base_run_returned"] = True
    probe_path.write_text(json.dumps(payload, indent=2, default=str))
    print(probe_path)
    return result


if __name__ == "__main__":
    main()


__all__ = [
    "POINT_RELEASE",
    "MODEL_ID",
    "patch_run_2d_source_v100510",
    "validate_source_transform_v100510",
    "main",
]
