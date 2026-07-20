"""v10.0.5.13.2 barrier-only point release.

The inherited v9.18.5.3 corridor selector still treated sampled ``h_tip/da`` as a
fatal startup condition even though the production v9.18.5.6 crack-advance gate
already treats the same metric as an audit warning.  This point release makes the
startup and runtime policies consistent:

* finite positive element areas and the initial triangle-quality floor remain
  mandatory;
* sampled ``h_tip/da`` is recorded and warned on, but is not a startup veto;
* no mesh, mechanics, constitutive, barrier, MPZ, shielding, or crack-advance law
  is changed.
"""
from __future__ import annotations

from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import sys
from typing import Any

import numpy as np

from . import mode_i_first_passage_v9_18_5 as _v9185
from . import mode_i_first_passage_v9_18_5_2 as _v91852
from . import mode_i_first_passage_v9_18_5_3 as _v91853
from . import mode_i_first_passage_v10_0_5_13_1_barrier_only as _base
from .mode_i_first_passage_v10_0 import _option_value

POINT_RELEASE = "10.0.5.13.2"
MODEL_ID = "FEM_CZM_full_2D_barrier_only_startup_resolution_warning_v10_0_5_13_2"
PRODUCTION_MANIFEST = _base.PRODUCTION_MANIFEST


def _quality_selected_corridor_mesh_v1005132(
    geom: Any,
    mesh_cfg: Any,
    seed: int | None = None,
    tip_center: Any | None = None,
):
    """Select a valid-quality corridor while auditing, not vetoing, h_tip/da."""

    original = _quality_selected_corridor_mesh_v1005132._original
    graded = float(getattr(mesh_cfg, "tip_h_fine", 0.0) or 0.0) > 0.0
    enabled = os.environ.get("ARRHENIUS_PREFINED_MODE_I_CORRIDOR", "1") not in {
        "0", "false", "False", "no", "NO"
    }
    if tip_center is not None or not graded or not enabled:
        return original(geom, mesh_cfg, seed=seed, tip_center=tip_center)

    target_um = max(
        _v91853._float_env("ARRHENIUS_COMMITTED_TARGET_EXTENSION_UM", 0.0), 0.0
    )
    guard_um = max(_v91853._float_env("ARRHENIUS_CORRIDOR_GUARD_UM", 10.0), 0.0)
    max_gap_um = max(
        _v91853._float_env("ARRHENIUS_CORRIDOR_MAX_CENTER_GAP_UM", 35.0), 1.0
    )
    da_um = max(_v91853._float_env("ARRHENIUS_PHYSICAL_DA_UM", 5.0), 1.0e-6)
    da_m = da_um * 1.0e-6
    qfloor = _v91853._float_env(
        "ARRHENIUS_MIN_INITIAL_TRIANGLE_QUALITY",
        _v91853._float_env("ARRHENIUS_MIN_ACCEPTED_TRIANGLE_QUALITY", 0.035),
    )
    warning_ratio = _v91853._float_env("ARRHENIUS_MAX_TIP_H_OVER_DA", 0.75)

    start = float(geom.a0)
    stop = min(float(geom.Lx), start + (target_um + guard_um) * 1.0e-6)
    length_m = max(stop - start, 0.0)
    counts = _v91853._candidate_counts(length_m * 1.0e6, max_gap_um)

    candidates: list[dict[str, Any]] = []
    accepted: list[
        tuple[tuple[float, float, float], Any, dict[str, Any], np.ndarray]
    ] = []
    for count in counts:
        centers = _v91853._centers_for_count(geom, length_m, count)
        try:
            raw = original(geom, mesh_cfg, seed=seed, tip_center=centers)
            compact, audit = _v91853._compact_without_quality_abort(raw, centers)
            resolution = _v91853._corridor_resolution(compact, start, stop, da_m)
            qmin = float(audit["minimum_initial_triangle_quality"])
            hratio = float(resolution["maximum_sampled_hbar_tip_over_da"])
            resolution_warning = bool(
                math.isfinite(hratio) and hratio > warning_ratio
            )
            ok = bool(qmin >= qfloor)
            row = {
                "center_count": int(len(centers)),
                "center_gap_um": float(
                    length_m * 1.0e6 / max(len(centers) - 1, 1)
                ),
                "node_count": int(compact.nn),
                "triangle_count": int(compact.ne),
                "minimum_initial_triangle_quality": qmin,
                "maximum_sampled_hbar_tip_over_da": hratio,
                "requested_maximum_tip_h_over_da": warning_ratio,
                "tip_h_over_da_enforced_as_veto": False,
                "resolution_warning": resolution_warning,
                "accepted": ok,
                "error": None,
            }
            candidates.append(row)
            if ok:
                # Retain the existing preference for quality margin, then choose
                # the better-resolved and smaller valid mesh.
                score = (qmin - qfloor, -hratio, -float(compact.nn))
                accepted.append((score, compact, {**audit, **resolution, **row}, centers))
        except Exception as exc:
            candidates.append(
                {
                    "center_count": int(count),
                    "accepted": False,
                    "error": f"{type(exc).__name__}: {exc}",
                }
            )

    if not accepted:
        _v91852._STARTUP_AUDIT.clear()
        _v91852._STARTUP_AUDIT.update(
            {
                "schema": "quality_selected_corridor_v1005132_v1",
                "candidate_corridors": candidates,
                "minimum_initial_triangle_quality_required": qfloor,
                "requested_maximum_tip_h_over_da": warning_ratio,
                "tip_h_over_da_role": "audit_warning_only",
                "tip_h_over_da_enforced_as_veto": False,
                "constitutive_physics_changed": False,
            }
        )
        raise RuntimeError(
            "v10.0.5.13.2 found no corridor satisfying initial triangle quality "
            f">= {qfloor:.6g}"
        )

    _, selected, selected_audit, centers = max(accepted, key=lambda item: item[0])
    hratio = float(selected_audit["maximum_sampled_hbar_tip_over_da"])
    payload = {
        **selected_audit,
        "schema": "quality_selected_corridor_v1005132_v1",
        "candidate_corridors": candidates,
        "selected_center_count": int(len(centers)),
        "selected_corridor_centers_m": centers.tolist(),
        "minimum_initial_triangle_quality_required": qfloor,
        "requested_maximum_tip_h_over_da": warning_ratio,
        "tip_h_over_da_role": "audit_warning_only",
        "tip_h_over_da_enforced_as_veto": False,
        "startup_resolution_warning": bool(
            math.isfinite(hratio) and hratio > warning_ratio
        ),
        "physical_da_um": da_um,
        "corridor_target_extension_um": target_um,
        "corridor_guard_um": guard_um,
        "constitutive_physics_changed": False,
    }
    _v91852._STARTUP_AUDIT.clear()
    _v91852._STARTUP_AUDIT.update(payload)
    _v9185._RUNTIME["corridor_centers"] = centers.tolist()
    _v9185._RUNTIME["mesh"] = selected
    return selected


def _update_manifest(out: Path, completed: bool) -> None:
    path = out / PRODUCTION_MANIFEST
    if not path.is_file():
        return
    try:
        payload = json.loads(path.read_text())
    except Exception:
        return

    startup_path = out / "compact_corridor_mesh_v91852.json"
    startup = {}
    if startup_path.is_file():
        try:
            startup = json.loads(startup_path.read_text())
        except Exception:
            startup = {}

    payload.update(
        {
            "schema": "barrier_only_production_manifest_v10_0_5_13_2",
            "model": MODEL_ID,
            "point_release": POINT_RELEASE,
            "startup_resolution_policy_repair": {
                "active": True,
                "initial_triangle_quality_veto_retained": True,
                "finite_positive_area_veto_retained": True,
                "tip_h_over_da_role": "audit_warning_only",
                "tip_h_over_da_enforced_as_startup_veto": False,
                "runtime_quality_policy_matched": True,
                "mesh_or_physics_changed": False,
                "startup_audit_schema": startup.get("schema"),
                "startup_resolution_warning": startup.get(
                    "startup_resolution_warning"
                ),
                "recorded_utc": datetime.now(timezone.utc).isoformat(),
            },
        }
    )
    if completed:
        payload["status"] = "complete"
        payload["run_completed_without_exception"] = True
    path.write_text(json.dumps(payload, indent=2, default=str))


def main(argv: list[str] | None = None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    out_value = _option_value(user_args, "--out")
    if out_value is None:
        raise SystemExit("v10.0.5.13.2 requires --out")
    out = Path(out_value).resolve()

    saved = _v91853._quality_selected_corridor_mesh
    _v91853._quality_selected_corridor_mesh = _quality_selected_corridor_mesh_v1005132
    try:
        result = _base.main(user_args)
        _update_manifest(out, completed=True)
        return result
    except BaseException:
        _update_manifest(out, completed=False)
        raise
    finally:
        _v91853._quality_selected_corridor_mesh = saved


if __name__ == "__main__":
    main()


__all__ = [
    "POINT_RELEASE",
    "MODEL_ID",
    "PRODUCTION_MANIFEST",
    "_quality_selected_corridor_mesh_v1005132",
    "main",
]
