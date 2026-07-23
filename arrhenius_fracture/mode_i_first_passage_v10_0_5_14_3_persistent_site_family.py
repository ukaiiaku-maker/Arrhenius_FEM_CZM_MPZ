"""v10.0.5.14.3: kernel-family parity with exact frozen-generator transport."""
from __future__ import annotations

import json
from pathlib import Path
import sys

from . import mode_i_first_passage_v10_0_5_14_1_persistent_site_family as _base
from . import persistent_site_front_engine_v100514 as _engine
from .persistent_site_transport_v1005143 import (
    TRANSPORT_INTEGRATOR,
    installed_exponential_transport_v1005143,
)

POINT_RELEASE = "10.0.5.14.3"
MODEL_ID = (
    "FEM_CZM_full_2D_PF_v10_2_22_persistent_site_kernel_family_"
    "exponential_transport_v10_0_5_14_3"
)
PRODUCTION_MANIFEST = "persistent_site_production_manifest_v10_0_5_14_3.json"
SELECTION_MANIFEST = "persistent_site_parameter_selection_v10_0_5_14_3.json"


def _out_path(argv: list[str]) -> Path | None:
    if "--out" not in argv:
        return None
    index = argv.index("--out")
    if index + 1 >= len(argv):
        return None
    return Path(argv[index + 1]).expanduser().resolve()


def _rewrite_release_metadata(out: Path | None) -> None:
    if out is None or not out.exists():
        return
    manifest = out / PRODUCTION_MANIFEST
    if manifest.is_file():
        payload = json.loads(manifest.read_text())
        payload["schema"] = "persistent_site_production_manifest_v10_0_5_14_3"
        payload["model"] = MODEL_ID
        payload["point_release"] = POINT_RELEASE
        physics = dict(payload.get("physics_contract", {}))
        physics.update(
            {
                "transport_integrator": TRANSPORT_INTEGRATOR,
                "frozen_linear_transport_solution": "matrix_exponential",
                "nonlinear_transport_error_control": "step_doubling",
                "explicit_transport_CFL_microstepping": False,
                "backward_euler_stiff_tail_refinement": False,
                "transport_equations_changed": False,
                "transport_time_integrator_changed": True,
            }
        )
        payload["physics_contract"] = physics
        manifest.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

    old_selection = out / "persistent_site_parameter_selection_v10_0_5_14_1.json"
    new_selection = out / SELECTION_MANIFEST
    if old_selection.is_file():
        payload = json.loads(old_selection.read_text())
        payload["schema"] = MODEL_ID
        payload["point_release"] = POINT_RELEASE
        policy = dict(payload.get("policy", {}))
        policy["transport_integrator"] = TRANSPORT_INTEGRATOR
        policy["transport_cfl_limited"] = False
        policy["frozen_linear_transport_solution"] = "matrix_exponential"
        payload["policy"] = policy
        new_selection.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        old_selection.unlink()


def main(argv: list[str] | None = None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    out = _out_path(user_args)
    saved = {
        "point_release": _base.POINT_RELEASE,
        "model_id": _base.MODEL_ID,
        "manifest": _base.PRODUCTION_MANIFEST,
        "engine_model": _engine.MODEL_ID,
    }
    _base.POINT_RELEASE = POINT_RELEASE
    _base.MODEL_ID = MODEL_ID
    _base.PRODUCTION_MANIFEST = PRODUCTION_MANIFEST
    _engine.MODEL_ID = "FEM_CZM_persistent_site_front_engine_v10_0_5_14_3"
    try:
        with installed_exponential_transport_v1005143():
            return _base.main(user_args)
    finally:
        _rewrite_release_metadata(out)
        _base.POINT_RELEASE = saved["point_release"]
        _base.MODEL_ID = saved["model_id"]
        _base.PRODUCTION_MANIFEST = saved["manifest"]
        _engine.MODEL_ID = saved["engine_model"]


if __name__ == "__main__":
    main()


__all__ = ["POINT_RELEASE", "MODEL_ID", "PRODUCTION_MANIFEST", "main"]
