"""Phase-C v10.0.5.12.3 production entry.

This point release preserves the fixed physical-refinement audit metadata across
CZM topology rebuilds.  The underlying mesh coordinates, connectivity, FEM,
kinetics, and crack-advance algorithms are unchanged; only dynamic audit
attributes are copied onto each rebuilt ``TriMesh``.
"""
from __future__ import annotations

from datetime import datetime, timezone
import json
import math
from pathlib import Path
import sys
from typing import Any

from . import crack_backend as _crack_backend
from . import mode_i_first_passage_v10_0_5_12_phase_c as _base
from .mode_i_first_passage_v10_0 import _option_value

POINT_RELEASE = "10.0.5.12.3"
MODEL_ID = "FEM_CZM_Phase_C_four_option_monotonic_v10_0_5_12_3"
PRODUCTION_MANIFEST = _base.PRODUCTION_MANIFEST


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _required_float(argv: list[str], name: str) -> float:
    value = _option_value(argv, name)
    if value is None:
        raise SystemExit(f"v10.0.5.12.3 requires {name}")
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise SystemExit(f"{name} must be finite and positive")
    return number


def _tip_centers(args: tuple[Any, ...], kwargs: dict[str, Any]) -> Any:
    if "tip_centers" in kwargs:
        return kwargs["tip_centers"]
    return args[2] if len(args) >= 3 else None


def _annotate_mesh(mesh: Any, radius_m: float, centers: Any) -> Any:
    mesh.production_refinement_radius_m = float(radius_m)
    mesh.production_refinement_policy = "fixed_physical_radius_same_radial_ring_law"
    if centers is None:
        prior = getattr(mesh, "production_refinement_centers_m", None)
        mesh.production_refinement_centers_m = prior
    else:
        try:
            mesh.production_refinement_centers_m = centers.tolist()
        except AttributeError:
            mesh.production_refinement_centers_m = centers
    mesh.production_refinement_metadata_propagated_v1005123 = True
    return mesh


def _update_manifest(out: Path, *, completed: bool) -> None:
    path = out / PRODUCTION_MANIFEST
    if not path.is_file():
        return
    try:
        payload = json.loads(path.read_text())
    except Exception:
        return
    payload.update(
        {
            "schema": "phase_c_production_manifest_v10_0_5_12_3",
            "model": MODEL_ID,
            "point_release": POINT_RELEASE,
            "metadata_propagation_fix": {
                "active": True,
                "scope": "dynamic_audit_attributes_only",
                "physics_changed": False,
                "topology_rebuild_metadata_preserved": True,
                "recorded_utc": _utc_now(),
            },
        }
    )
    if completed:
        payload["status"] = "complete"
        payload["run_completed_without_exception"] = True
    path.write_text(json.dumps(payload, indent=2, default=str))


def main(argv: list[str] | None = None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    radius_m = _required_float(user_args, "--tip-refinement-radius-um") * 1.0e-6
    out_value = _option_value(user_args, "--out")
    if out_value is None:
        raise SystemExit("v10.0.5.12.3 requires --out")
    out = Path(out_value).resolve()

    original_rebuild = _crack_backend.rebuild_tri_mesh

    def rebuild_with_refinement_metadata(*args, **kwargs):
        mesh = original_rebuild(*args, **kwargs)
        return _annotate_mesh(mesh, radius_m, _tip_centers(args, kwargs))

    _crack_backend.rebuild_tri_mesh = rebuild_with_refinement_metadata
    try:
        result = _base.main(user_args)
        _update_manifest(out, completed=True)
        return result
    except BaseException:
        _update_manifest(out, completed=False)
        raise
    finally:
        _crack_backend.rebuild_tri_mesh = original_rebuild


if __name__ == "__main__":
    main()


__all__ = ["POINT_RELEASE", "MODEL_ID", "PRODUCTION_MANIFEST", "main"]
