"""v9.18.5.3 quality-selected Mode-I corridor mesh.

v9.18.5.2 correctly rejected a fixed-spacing multi-cloud corridor whose worst
triangle quality was 0.030357, below the production floor 0.035.  This revision
constructs several deterministic candidate corridors, compacts each mesh, and
selects the best candidate that simultaneously satisfies the initial triangle
quality floor and the configured tip-resolution ratio h_tip / da.

No barrier, hazard, cohesive, MPZ, wake, shielding, or material law changes.
"""
from __future__ import annotations

import math
import os
from typing import Any

import numpy as np

from . import mesh as _mesh
from . import mode_i_first_passage_v9_18_5 as _v9185
from . import mode_i_first_passage_v9_18_5_2 as _v91852


def _float_env(name: str, default: float) -> float:
    try:
        value = float(os.environ.get(name, default))
    except (TypeError, ValueError):
        value = float(default)
    return value if math.isfinite(value) else float(default)


def _candidate_counts(length_um: float, max_gap_um: float) -> list[int]:
    base = max(2, int(math.ceil(max(length_um, 0.0) / max(max_gap_um, 1.0))) + 1)
    lo = max(2, base - 2)
    hi = base + 3
    return list(range(lo, hi + 1))


def _centers_for_count(geom: Any, length_m: float, count: int) -> np.ndarray:
    start = float(geom.a0)
    stop = min(float(geom.Lx), start + max(float(length_m), 0.0))
    if stop <= start + 1.0e-15:
        return np.array([[start, 0.0]], dtype=float)
    xs = np.linspace(start, stop, max(int(count), 2))
    return np.column_stack([xs, np.zeros_like(xs)])


def _compact_without_quality_abort(raw: Any, centers: np.ndarray):
    key = "ARRHENIUS_MIN_INITIAL_TRIANGLE_QUALITY"
    old = os.environ.get(key)
    os.environ[key] = "0"
    try:
        return _v91852._compact_mesh(raw, centers)
    finally:
        if old is None:
            os.environ.pop(key, None)
        else:
            os.environ[key] = old


def _corridor_resolution(mesh: Any, start: float, stop: float, da_m: float) -> dict[str, Any]:
    if stop <= start + 1.0e-15:
        samples = np.array([start], dtype=float)
    else:
        n = max(int(math.ceil((stop - start) / max(da_m, 1.0e-12))) + 1, 2)
        samples = np.linspace(start, stop, n)
    h = np.array([
        _mesh._estimate_hbar_tip(mesh.nodes, mesh.elems, float(x), 0.0)
        for x in samples
    ], dtype=float)
    return {
        "sample_x_m": samples.tolist(),
        "sample_hbar_tip_m": h.tolist(),
        "maximum_sampled_hbar_tip_m": float(np.max(h)),
        "maximum_sampled_hbar_tip_over_da": float(np.max(h) / max(da_m, 1.0e-300)),
    }


def _quality_selected_corridor_mesh(geom, mesh_cfg, seed=None, tip_center=None):
    original = _quality_selected_corridor_mesh._original
    graded = float(getattr(mesh_cfg, "tip_h_fine", 0.0) or 0.0) > 0.0
    enabled = os.environ.get("ARRHENIUS_PREFINED_MODE_I_CORRIDOR", "1") not in {
        "0", "false", "False", "no", "NO"
    }
    if tip_center is not None or not graded or not enabled:
        return original(geom, mesh_cfg, seed=seed, tip_center=tip_center)

    target_um = max(_float_env("ARRHENIUS_COMMITTED_TARGET_EXTENSION_UM", 0.0), 0.0)
    guard_um = max(_float_env("ARRHENIUS_CORRIDOR_GUARD_UM", 10.0), 0.0)
    max_gap_um = max(_float_env("ARRHENIUS_CORRIDOR_MAX_CENTER_GAP_UM", 35.0), 1.0)
    da_um = max(_float_env("ARRHENIUS_PHYSICAL_DA_UM", 5.0), 1.0e-6)
    da_m = da_um * 1.0e-6
    qfloor = _float_env(
        "ARRHENIUS_MIN_INITIAL_TRIANGLE_QUALITY",
        _float_env("ARRHENIUS_MIN_ACCEPTED_TRIANGLE_QUALITY", 0.035),
    )
    max_h_ratio = _float_env("ARRHENIUS_MAX_TIP_H_OVER_DA", 0.75)

    start = float(geom.a0)
    stop = min(float(geom.Lx), start + (target_um + guard_um) * 1.0e-6)
    length_m = max(stop - start, 0.0)
    counts = _candidate_counts(length_m * 1.0e6, max_gap_um)

    candidates: list[dict[str, Any]] = []
    accepted: list[tuple[tuple[float, float, float], Any, dict[str, Any], np.ndarray]] = []
    for count in counts:
        centers = _centers_for_count(geom, length_m, count)
        try:
            raw = original(geom, mesh_cfg, seed=seed, tip_center=centers)
            compact, audit = _compact_without_quality_abort(raw, centers)
            resolution = _corridor_resolution(compact, start, stop, da_m)
            qmin = float(audit["minimum_initial_triangle_quality"])
            hratio = float(resolution["maximum_sampled_hbar_tip_over_da"])
            ok = bool(qmin >= qfloor and hratio <= max_h_ratio)
            row = {
                "center_count": int(len(centers)),
                "center_gap_um": float(length_m * 1.0e6 / max(len(centers) - 1, 1)),
                "node_count": int(compact.nn),
                "triangle_count": int(compact.ne),
                "minimum_initial_triangle_quality": qmin,
                "maximum_sampled_hbar_tip_over_da": hratio,
                "accepted": ok,
                "error": None,
            }
            candidates.append(row)
            if ok:
                # Prefer quality margin first, then lower h/da, then lower node count.
                score = (qmin - qfloor, max_h_ratio - hratio, -float(compact.nn))
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
            "schema": "quality_selected_corridor_v91853_v1",
            "candidate_corridors": candidates,
            "minimum_initial_triangle_quality_required": qfloor,
            "maximum_tip_h_over_da_required": max_h_ratio,
            "constitutive_physics_changed": False,
        })
        raise RuntimeError(
            "v9.18.5.3 found no corridor satisfying both initial triangle quality "
            f">= {qfloor:.6g} and sampled h_tip/da <= {max_h_ratio:.6g}"
        )

    _, selected, selected_audit, centers = max(accepted, key=lambda item: item[0])
    payload = {
        **selected_audit,
        "schema": "quality_selected_corridor_v91853_v1",
        "candidate_corridors": candidates,
        "selected_center_count": int(len(centers)),
        "selected_corridor_centers_m": centers.tolist(),
        "minimum_initial_triangle_quality_required": qfloor,
        "maximum_tip_h_over_da_required": max_h_ratio,
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


def main(argv=None):
    original = _v91852._compact_corridor_mesh
    _v91852._compact_corridor_mesh = _quality_selected_corridor_mesh
    try:
        return _v91852.main(argv)
    finally:
        _v91852._compact_corridor_mesh = original


if __name__ == "__main__":
    main()
