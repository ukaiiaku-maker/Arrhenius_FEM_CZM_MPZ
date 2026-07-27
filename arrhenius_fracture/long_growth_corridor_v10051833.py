"""Target-aware long-growth corridor mesh policy for v10.0.5.18.3.3.

This module changes no constitutive, hazard, event-length, shielding, or cohesive
law.  It corrects the numerical mesh contract used by long-growth campaigns:

* the requested committed crack extension is propagated to the corridor builder;
* the static refinement corridor spans that extension plus a guard distance;
* corridor acceptance is based on triangle quality and process-zone resolution
  ``h_tip / L_pz`` rather than the stochastic event-length ratio ``h_tip / da``;
* ``h_tip / da`` remains recorded as an audit warning only, matching v9.18.5.6.
"""
from __future__ import annotations

import math
import os
from typing import Any

import numpy as np

from . import mode_i_first_passage_v9_18_5 as _v9185
from . import mode_i_first_passage_v9_18_5_2 as _v91852
from . import mode_i_first_passage_v9_18_5_3 as _v91853

MODEL_ID = "target_aware_process_zone_corridor_v10_0_5_18_3_3"
CORRIDOR_SCHEMA = "v10.0.5.18.3.3_target_aware_process_zone_corridor"


def _float_env(name: str, default: float) -> float:
    try:
        value = float(os.environ.get(name, default))
    except (TypeError, ValueError):
        value = float(default)
    return value if math.isfinite(value) else float(default)


def _required_positive_env(name: str) -> float:
    value = _float_env(name, float("nan"))
    if not math.isfinite(value) or value <= 0.0:
        raise RuntimeError(
            f"v10.0.5.18.3.3 requires positive {name}; observed {os.environ.get(name)!r}"
        )
    return value


def target_aware_long_growth_corridor_mesh(
    geom,
    mesh_cfg,
    seed=None,
    tip_center=None,
):
    """Build a quality-selected corridor spanning the requested growth target."""
    original = target_aware_long_growth_corridor_mesh._original
    graded = float(getattr(mesh_cfg, "tip_h_fine", 0.0) or 0.0) > 0.0
    enabled = os.environ.get("ARRHENIUS_PREFINED_MODE_I_CORRIDOR", "1") not in {
        "0", "false", "False", "no", "NO"
    }
    if tip_center is not None or not graded or not enabled:
        return original(geom, mesh_cfg, seed=seed, tip_center=tip_center)

    target_um = _required_positive_env("ARRHENIUS_COMMITTED_TARGET_EXTENSION_UM")
    guard_um = max(_float_env("ARRHENIUS_CORRIDOR_GUARD_UM", 10.0), 0.0)
    max_gap_um = max(_float_env("ARRHENIUS_CORRIDOR_MAX_CENTER_GAP_UM", 100.0), 1.0)
    da_um = _required_positive_env("ARRHENIUS_PHYSICAL_DA_UM")
    process_zone_um = _required_positive_env("ARRHENIUS_CORRIDOR_PROCESS_ZONE_UM")
    max_h_over_lpz = max(
        _float_env("ARRHENIUS_MAX_CORRIDOR_H_OVER_LPZ", 0.25),
        1.0e-6,
    )
    audit_h_over_da = max(
        _float_env("ARRHENIUS_MAX_TIP_H_OVER_DA", 0.75),
        1.0e-6,
    )
    qfloor = _float_env(
        "ARRHENIUS_MIN_INITIAL_TRIANGLE_QUALITY",
        _float_env("ARRHENIUS_MIN_ACCEPTED_TRIANGLE_QUALITY", 0.035),
    )

    da_m = da_um * 1.0e-6
    lpz_m = process_zone_um * 1.0e-6
    start = float(geom.a0)
    stop = min(float(geom.Lx), start + (target_um + guard_um) * 1.0e-6)
    requested_stop = start + (target_um + guard_um) * 1.0e-6
    if stop + 1.0e-15 < requested_stop:
        raise RuntimeError(
            "v10.0.5.18.3.3 corridor target exceeds specimen ligament: "
            f"requested_stop={requested_stop:.9e} domain_stop={float(geom.Lx):.9e}"
        )
    length_m = max(stop - start, 0.0)
    counts = _v91853._candidate_counts(length_m * 1.0e6, max_gap_um)

    candidates: list[dict[str, Any]] = []
    accepted: list[tuple[tuple[float, float, float], Any, dict[str, Any], np.ndarray]] = []
    for count in counts:
        centers = _v91853._centers_for_count(geom, length_m, count)
        try:
            raw = original(geom, mesh_cfg, seed=seed, tip_center=centers)
            compact, audit = _v91853._compact_without_quality_abort(raw, centers)
            resolution = _v91853._corridor_resolution(compact, start, stop, da_m)
            qmin = float(audit["minimum_initial_triangle_quality"])
            hmax = float(resolution["maximum_sampled_hbar_tip_m"])
            h_over_da = hmax / max(da_m, 1.0e-300)
            h_over_lpz = hmax / max(lpz_m, 1.0e-300)
            center_gap_um = float(
                length_m * 1.0e6 / max(len(centers) - 1, 1)
            )
            gap_ok = bool(center_gap_um <= max_gap_um + 1.0e-12)
            quality_ok = bool(qmin >= qfloor)
            process_zone_ok = bool(h_over_lpz <= max_h_over_lpz)
            ok = bool(gap_ok and quality_ok and process_zone_ok)
            row = {
                "center_count": int(len(centers)),
                "center_gap_um": center_gap_um,
                "maximum_center_gap_um_required": max_gap_um,
                "center_gap_requirement_passed": gap_ok,
                "node_count": int(compact.nn),
                "triangle_count": int(compact.ne),
                "minimum_initial_triangle_quality": qmin,
                "triangle_quality_required": qfloor,
                "maximum_sampled_hbar_tip_m": hmax,
                "maximum_sampled_hbar_tip_over_da": h_over_da,
                "maximum_sampled_hbar_tip_over_L_pz": h_over_lpz,
                "maximum_hbar_tip_over_L_pz_required": max_h_over_lpz,
                "tip_h_over_da_audit_threshold": audit_h_over_da,
                "tip_h_over_da_enforced_as_veto": False,
                "tip_h_over_da_resolution_warning": bool(h_over_da > audit_h_over_da),
                "accepted": ok,
                "error": None,
            }
            candidates.append(row)
            if ok:
                score = (
                    -float(compact.nn),
                    qmin - qfloor,
                    max_h_over_lpz - h_over_lpz,
                )
                accepted.append((score, compact, {**audit, **resolution, **row}, centers))
        except Exception as exc:
            candidates.append({
                "center_count": int(count),
                "accepted": False,
                "error": f"{type(exc).__name__}: {exc}",
            })

    if not accepted:
        _v91852._STARTUP_AUDIT.clear()
        _v91852._STARTUP_AUDIT.update({
            "schema": CORRIDOR_SCHEMA,
            "candidate_corridors": candidates,
            "corridor_target_extension_um": target_um,
            "corridor_guard_um": guard_um,
            "corridor_max_center_gap_um": max_gap_um,
            "process_zone_length_um": process_zone_um,
            "maximum_hbar_tip_over_L_pz_required": max_h_over_lpz,
            "minimum_initial_triangle_quality_required": qfloor,
            "tip_h_over_da_role": "audit_warning_only",
            "constitutive_physics_changed": False,
        })
        raise RuntimeError(
            "v10.0.5.18.3.3 found no corridor satisfying maximum center gap "
            f"<= {max_gap_um:.6g} um, initial triangle quality >= {qfloor:.6g}, "
            f"and h_tip/L_pz <= {max_h_over_lpz:.6g} over "
            f"{target_um + guard_um:.6g} um"
        )

    _, selected, selected_audit, centers = max(accepted, key=lambda item: item[0])
    payload = {
        **selected_audit,
        "schema": CORRIDOR_SCHEMA,
        "model": MODEL_ID,
        "candidate_corridors": candidates,
        "selected_center_count": int(len(centers)),
        "selected_corridor_centers_m": centers.tolist(),
        "corridor_target_extension_um": target_um,
        "corridor_guard_um": guard_um,
        "corridor_max_center_gap_um": max_gap_um,
        "process_zone_length_um": process_zone_um,
        "maximum_hbar_tip_over_L_pz_required": max_h_over_lpz,
        "minimum_initial_triangle_quality_required": qfloor,
        "physical_da_um": da_um,
        "tip_h_over_da_role": "audit_warning_only",
        "target_propagated_before_mesh_construction": True,
        "full_requested_corridor_covered": True,
        "constitutive_physics_changed": False,
    }
    _v91852._STARTUP_AUDIT.clear()
    _v91852._STARTUP_AUDIT.update(payload)
    _v9185._RUNTIME["corridor_centers"] = centers.tolist()
    _v9185._RUNTIME["mesh"] = selected
    return selected


__all__ = [
    "CORRIDOR_SCHEMA",
    "MODEL_ID",
    "target_aware_long_growth_corridor_mesh",
]
