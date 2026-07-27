"""Target-aware long-growth corridor mesh policy for v10.0.5.18.3.3.

This module changes no constitutive, hazard, event-length, shielding, or cohesive
law.  It corrects the numerical mesh contract used by long-growth campaigns:

* the requested committed crack extension is propagated to the corridor builder;
* the static refinement corridor spans that extension plus a guard distance;
* corridor acceptance is based on triangle quality and process-zone resolution
  ``h_tip / L_pz`` rather than the stochastic event-length ratio ``h_tip / da``;
* ``h_tip / da`` remains recorded as an audit warning only, matching v9.18.5.6;
* the deterministic center-count search starts at the minimum count required by
  the maximum-gap contract and continues until the first admissible mesh is found;
* when the production 330 um physical-refinement provider is active, one swept
  capsule is generated around the complete path instead of superposing a full
  radial disk around every closely spaced corridor center.
"""
from __future__ import annotations

import math
import os
from typing import Any

import numpy as np
from scipy.spatial import Delaunay

from . import mesh as _mesh
from . import mode_i_first_passage_v9_18_5 as _v9185
from . import mode_i_first_passage_v9_18_5_2 as _v91852
from . import mode_i_first_passage_v9_18_5_3 as _v91853
from . import physical_refinement_mesh_v100510 as _physical

MODEL_ID = "target_aware_process_zone_corridor_v10_0_5_18_3_3"
CORRIDOR_SCHEMA = "v10.0.5.18.3.3_target_aware_process_zone_corridor"
SWEPT_PHYSICAL_POLICY = "fixed_physical_radius_swept_capsule_same_size_law"


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


def _production_candidate_counts(length_um: float, max_gap_um: float) -> list[int]:
    """Return an ascending bounded search satisfying the center-gap contract."""
    minimum = max(
        2,
        int(math.ceil(max(float(length_um), 0.0) / max(float(max_gap_um), 1.0))) + 1,
    )
    extra = int(max(_float_env("ARRHENIUS_CORRIDOR_EXTRA_CENTER_COUNTS", 18.0), 0.0))
    return list(range(minimum, minimum + extra + 1))


def _is_production_physical_constructor(function: Any) -> bool:
    return (
        getattr(function, "__module__", "") == _physical.__name__
        and getattr(function, "__name__", "") == "make_physical_refinement_mesh_v100510"
    )


def _distance_to_horizontal_segment(
    x: float,
    y: float,
    start: float,
    stop: float,
) -> float:
    if x < start:
        dx = start - x
    elif x > stop:
        dx = x - stop
    else:
        dx = 0.0
    return math.hypot(dx, y)


def _swept_capsule_nodes(
    geom: Any,
    mesh_cfg: Any,
    start: float,
    stop: float,
    radius: float,
) -> np.ndarray:
    """Generate a graded point cloud from distance to the requested path segment.

    The local spacing law is exactly the production radial law,
    ``h(d)=min(h_fine + slope*d, h_far)``, with ``d`` measured from the complete
    path segment rather than from several overlapping point centers.
    """
    Lx = float(geom.Lx)
    Ly = float(geom.Ly)
    h_fine = float(getattr(mesh_cfg, "tip_h_fine", 0.0) or 0.0)
    if h_fine <= 0.0:
        raise ValueError("swept physical corridor requires a graded tip mesh")
    ratio = float(getattr(mesh_cfg, "tip_ratio", 1.15) or 1.15)
    slope = max(ratio - 1.0, 0.02)
    h_far = max(Lx, Ly) / 40.0

    def size(distance: float) -> float:
        return min(h_fine + slope * max(float(distance), 0.0), h_far)

    offsets = [0.0]
    distance = h_fine
    while distance < radius:
        offsets.extend((-distance, distance))
        distance += size(distance)
    offsets = sorted(set(offsets))

    rows: list[np.ndarray] = []
    for row_index, y in enumerate(offsets):
        cap = math.sqrt(max(radius * radius - y * y, 0.0))
        left = max(0.0, start - cap)
        right = min(Lx, stop + cap)
        if right <= left:
            continue

        points: list[tuple[float, float]] = []
        x = left
        if row_index % 2:
            x += 0.5 * size(abs(y))
        while x <= right + 1.0e-15:
            points.append((x, y))
            local_distance = _distance_to_horizontal_segment(x, y, start, stop)
            x += size(local_distance)
        if points:
            rows.append(np.asarray(points, dtype=float))

    nbx = max(2, int(Lx / h_far))
    nby = max(2, int(Ly / h_far))
    gx, gy = np.meshgrid(
        np.linspace(0.0, Lx, nbx + 1),
        np.linspace(-0.5 * Ly, 0.5 * Ly, nby + 1),
    )
    background = np.column_stack((gx.ravel(), gy.ravel()))
    background_distance = np.array(
        [
            _distance_to_horizontal_segment(float(x), float(y), start, stop)
            for x, y in background
        ],
        dtype=float,
    )
    background = background[background_distance > 1.05 * radius]

    mandatory = np.asarray(
        [
            [start, 0.0],
            [stop, 0.0],
            [0.0, -0.5 * Ly],
            [0.0, 0.5 * Ly],
            [Lx, -0.5 * Ly],
            [Lx, 0.5 * Ly],
        ],
        dtype=float,
    )
    nodes = np.vstack([*rows, background, mandatory])
    key = np.round(nodes / max(1.0e-12, 0.05 * h_fine)).astype(np.int64)
    _, keep = np.unique(key, axis=0, return_index=True)
    return nodes[np.sort(keep)]


def _make_swept_physical_corridor_mesh(
    geom: Any,
    mesh_cfg: Any,
    centers: np.ndarray,
    start: float,
    stop: float,
    seed: int | None,
):
    spec = _physical._ACTIVE_SPEC
    if spec is None:
        raise RuntimeError("production physical refinement mesh was not configured")
    radius = float(spec.radius_m)
    spec.validate(
        geom,
        np.asarray([[start, 0.0], [stop, 0.0]], dtype=float),
    )
    if seed is not None:
        np.random.seed(seed)

    nodes = _swept_capsule_nodes(geom, mesh_cfg, start, stop, radius)
    triangulation = Delaunay(nodes)
    elems = triangulation.simplices
    centroids = nodes[elems].mean(axis=1)
    inside = (
        (centroids[:, 0] >= 0.0)
        & (centroids[:, 0] <= float(geom.Lx))
        & (centroids[:, 1] >= -0.5 * float(geom.Ly))
        & (centroids[:, 1] <= 0.5 * float(geom.Ly))
    )
    mesh = _mesh.rebuild_tri_mesh(
        nodes,
        elems[inside],
        tip_centers=np.asarray(centers, dtype=float),
        validate=True,
    )
    mesh.production_refinement_radius_m = radius
    mesh.production_refinement_policy = SWEPT_PHYSICAL_POLICY
    mesh.production_refinement_centers_m = np.asarray(centers, dtype=float).tolist()
    mesh.production_refinement_h_fine_m = float(
        getattr(mesh_cfg, "tip_h_fine", 0.0) or 0.0
    )
    mesh.production_refinement_tip_ratio = float(
        getattr(mesh_cfg, "tip_ratio", 1.15) or 1.15
    )
    mesh.production_refinement_corridor_start_m = float(start)
    mesh.production_refinement_corridor_stop_m = float(stop)
    mesh.production_refinement_swept_capsule = True
    return mesh


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
    counts = _production_candidate_counts(length_m * 1.0e6, max_gap_um)
    physical_provider = _is_production_physical_constructor(original)

    candidates: list[dict[str, Any]] = []
    selected = None
    selected_audit: dict[str, Any] | None = None
    selected_centers: np.ndarray | None = None

    for count in counts:
        centers = _v91853._centers_for_count(geom, length_m, count)
        try:
            if physical_provider:
                raw = _make_swept_physical_corridor_mesh(
                    geom,
                    mesh_cfg,
                    centers,
                    start,
                    stop,
                    seed,
                )
            else:
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
                "mesh_provider": (
                    SWEPT_PHYSICAL_POLICY if physical_provider else "inherited_constructor"
                ),
                "production_refinement_radius_m": getattr(
                    compact, "production_refinement_radius_m", None
                ),
                "accepted": ok,
                "error": None,
            }
            candidates.append(row)
            if ok:
                selected = compact
                selected_audit = {**audit, **resolution, **row}
                selected_centers = centers
                break
        except Exception as exc:
            candidates.append({
                "center_count": int(count),
                "mesh_provider": (
                    SWEPT_PHYSICAL_POLICY if physical_provider else "inherited_constructor"
                ),
                "accepted": False,
                "error": f"{type(exc).__name__}: {exc}",
            })

    constructor_module = getattr(original, "__module__", None)
    constructor_qualname = getattr(original, "__qualname__", None)
    if selected is None or selected_audit is None or selected_centers is None:
        _v91852._STARTUP_AUDIT.clear()
        _v91852._STARTUP_AUDIT.update({
            "schema": CORRIDOR_SCHEMA,
            "candidate_corridors": candidates,
            "candidate_center_counts": counts,
            "candidate_search_stopped_after_first_admissible": True,
            "corridor_target_extension_um": target_um,
            "corridor_guard_um": guard_um,
            "corridor_max_center_gap_um": max_gap_um,
            "process_zone_length_um": process_zone_um,
            "maximum_hbar_tip_over_L_pz_required": max_h_over_lpz,
            "minimum_initial_triangle_quality_required": qfloor,
            "tip_h_over_da_role": "audit_warning_only",
            "mesh_constructor_module": constructor_module,
            "mesh_constructor_qualname": constructor_qualname,
            "production_physical_provider_detected": physical_provider,
            "swept_physical_refinement_active": physical_provider,
            "constitutive_physics_changed": False,
        })
        raise RuntimeError(
            "v10.0.5.18.3.3 found no corridor satisfying maximum center gap "
            f"<= {max_gap_um:.6g} um, initial triangle quality >= {qfloor:.6g}, "
            f"and h_tip/L_pz <= {max_h_over_lpz:.6g} over "
            f"{target_um + guard_um:.6g} um after center counts "
            f"{counts[0]} through {counts[-1]}"
        )

    payload = {
        **selected_audit,
        "schema": CORRIDOR_SCHEMA,
        "model": MODEL_ID,
        "candidate_corridors": candidates,
        "candidate_center_counts": counts,
        "candidate_search_stopped_after_first_admissible": True,
        "selected_center_count": int(len(selected_centers)),
        "selected_corridor_centers_m": selected_centers.tolist(),
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
        "mesh_constructor_module": constructor_module,
        "mesh_constructor_qualname": constructor_qualname,
        "production_physical_provider_detected": physical_provider,
        "swept_physical_refinement_active": physical_provider,
        "swept_physical_refinement_policy": (
            SWEPT_PHYSICAL_POLICY if physical_provider else None
        ),
        "constitutive_physics_changed": False,
    }
    _v91852._STARTUP_AUDIT.clear()
    _v91852._STARTUP_AUDIT.update(payload)
    _v9185._RUNTIME["corridor_centers"] = selected_centers.tolist()
    _v9185._RUNTIME["mesh"] = selected
    return selected


__all__ = [
    "CORRIDOR_SCHEMA",
    "MODEL_ID",
    "SWEPT_PHYSICAL_POLICY",
    "target_aware_long_growth_corridor_mesh",
]
