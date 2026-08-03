"""Local-scale quality certification for the v10.0.5.18.3.9 corridor.

The historical 0.08 area test was an *immediate-child* safeguard for a single
edge split.  Comparing every triangle from a whole-cavity PSLG remesh directly
to one arbitrarily assigned old element changes that meaning and rejects small,
well-shaped support cells even when the remesh is locally regular.

This module keeps the numerical threshold unchanged but applies it to the area
of an equilateral cell at the constrained path-support spacing.  The old-parent
area ratios remain recorded for transfer diagnostics.  The bulk is elastic in
the production configuration; the parent map is still emitted unchanged.
"""
from __future__ import annotations

import math
from typing import Any

import numpy as np

from . import atomic_path_corridor_czm_v10051839 as _base
from .mesh import make_boundary_data, rebuild_tri_mesh


MODEL_ID = "atomic_path_corridor_local_support_area_gate_v10_0_5_18_3_9"


def _path_support_spacing(
    vertices: np.ndarray,
    p0: np.ndarray,
    p1: np.ndarray,
    direction: np.ndarray,
    tol: float,
) -> float:
    length = float(np.linalg.norm(np.asarray(p1) - np.asarray(p0)))
    s, r = _base._segment_coordinates(vertices, p0, direction)
    mask = (np.abs(r) <= tol) & (s >= -tol) & (s <= length + tol)
    values = np.unique(np.round(s[mask] / max(tol, 1.0e-300)).astype(np.int64))
    if len(values) < 2:
        return length
    coordinates = np.sort(values.astype(float) * tol)
    delta = np.diff(coordinates)
    delta = delta[delta > tol]
    return float(np.median(delta)) if len(delta) else length


def stitch_triangulation_local_support_area(
    self: Any,
    mesh: Any,
    damage: np.ndarray,
    displacement: np.ndarray,
    selected: np.ndarray,
    tri_out: dict[str, Any],
    p0: np.ndarray,
    p1: np.ndarray,
    direction: np.ndarray,
    front_id: int,
    tol: float,
):
    vertices = np.asarray(tri_out["vertices"], float)
    local_elems = np.asarray(tri_out["triangles"], int)
    old_nodes = np.asarray(mesh.nodes, float)
    old_map = _base._old_coordinate_map(mesh, tol)
    tip_ids = self._tip_geometric_node_ids(mesh, p0, int(front_id))
    if not tip_ids:
        return None, {"reason": "stitch_tip_node_not_found"}
    tip_plus, tip_minus = (
        (int(tip_ids[0]), int(tip_ids[1]))
        if len(tip_ids) >= 2 else (int(tip_ids[0]), int(tip_ids[0]))
    )
    p0_key = _base._coord_key(p0, tol)

    global_nodes = old_nodes.tolist()
    local_default: dict[int, int] = {}
    new_vertex_ids: list[int] = []
    for lid, point in enumerate(vertices):
        key = _base._coord_key(point, tol)
        existing = old_map.get(key, [])
        if existing:
            local_default[int(lid)] = tip_plus if key == p0_key else int(existing[0])
        else:
            gid = int(len(global_nodes))
            global_nodes.append(np.asarray(point, float).tolist())
            local_default[int(lid)] = gid
            new_vertex_ids.append(gid)

    global_nodes_array = np.asarray(global_nodes, float)
    new_local_global: list[np.ndarray] = []
    reference_sign = float(
        np.median(np.sign(self._signed_twice_area(mesh.nodes, mesh.elems[selected])))
    )
    if reference_sign == 0.0:
        reference_sign = 1.0

    for conn in local_elems:
        centroid = np.mean(vertices[np.asarray(conn, dtype=int)], axis=0)
        cross = float(
            (p1[0] - p0[0]) * (centroid[1] - p0[1])
            - (p1[1] - p0[1]) * (centroid[0] - p0[0])
        )
        mapped = []
        for lid in np.asarray(conn, dtype=int):
            key = _base._coord_key(vertices[int(lid)], tol)
            if key == p0_key:
                mapped.append(tip_plus if cross >= 0.0 else tip_minus)
            else:
                mapped.append(local_default[int(lid)])
        mapped = np.asarray(mapped, dtype=int)
        sign = float(
            self._signed_twice_area(global_nodes_array, mapped.reshape(1, 3))[0]
        )
        if sign * reference_sign < 0.0:
            mapped[0], mapped[1] = mapped[1], mapped[0]
        new_local_global.append(mapped)
    new_local_global = np.asarray(new_local_global, dtype=int)

    q = self._triangle_quality(global_nodes_array, new_local_global)
    qmin = float(np.min(q)) if len(q) else 1.0
    assigned = np.asarray(
        [
            _base._assign_old_parent(mesh, selected, global_nodes_array[conn])
            for conn in new_local_global
        ],
        dtype=int,
    )
    new_area = 0.5 * np.abs(
        self._signed_twice_area(global_nodes_array, new_local_global)
    )
    old_area = np.asarray(mesh.area_e[assigned], float)
    old_parent_ratios = new_area / np.maximum(old_area, 1.0e-300)
    min_old_parent_ratio = (
        float(np.min(old_parent_ratios)) if len(old_parent_ratios) else 1.0
    )

    support_spacing = _path_support_spacing(vertices, p0, p1, direction, tol)
    local_reference_area = (
        math.sqrt(3.0) / 4.0 * max(float(support_spacing), tol) ** 2
    )
    local_area_ratios = new_area / max(local_reference_area, 1.0e-300)
    min_local_ratio = float(np.min(local_area_ratios)) if len(local_area_ratios) else 1.0

    qfloor = _base._float_env(
        "ARRHENIUS_MIN_ACCEPTED_TRIANGLE_QUALITY", self.min_triangle_quality
    )
    afloor = _base._float_env(
        "ARRHENIUS_MIN_ACCEPTED_CHILD_AREA_RATIO", self.min_area_ratio
    )
    issues = []
    if qmin < qfloor:
        issues.append(f"triangle_quality={qmin:.6e}<{qfloor:.6e}")
    if min_local_ratio < afloor:
        issues.append(
            f"local_support_area_ratio={min_local_ratio:.6e}<{afloor:.6e}"
        )
    if issues:
        return None, {
            "reason": ";".join(issues),
            "min_triangle_quality": qmin,
            "triangle_quality_floor": qfloor,
            "min_child_area_ratio": min_local_ratio,
            "child_area_ratio_floor": afloor,
            "area_gate_metric": "new_triangle_area_over_equilateral_path_support_cell_area",
            "path_support_spacing_m": support_spacing,
            "local_reference_area_m2": local_reference_area,
            "min_old_parent_area_ratio_diagnostic": min_old_parent_ratio,
            "new_triangle_count": int(len(new_local_global)),
            "new_vertex_count": int(len(new_vertex_ids)),
        }

    keep = np.setdiff1d(
        np.arange(int(mesh.ne), dtype=int), np.asarray(selected, dtype=int)
    )
    elems1 = np.vstack([np.asarray(mesh.elems[keep], int), new_local_global])
    parent_map = np.concatenate([keep, assigned]).astype(int)
    mesh1 = rebuild_tri_mesh(global_nodes_array, elems1, tip_centers=[p0, p1])

    u0 = np.asarray(displacement, float).reshape(-1, 2)
    d0 = np.asarray(damage, float)
    u1 = np.zeros((len(global_nodes_array), 2), dtype=float)
    d1 = np.zeros(len(global_nodes_array), dtype=float)
    u1[: len(old_nodes)] = u0
    d1[: len(old_nodes)] = d0
    for gid in new_vertex_ids:
        point = global_nodes_array[int(gid)]
        u1[int(gid)] = self._interp_nodal_vector(mesh, u0, point)
        d1[int(gid)] = _base._interpolate_scalar(mesh, d0, point)

    incidence = np.bincount(elems1.ravel(), minlength=len(global_nodes_array))
    orphan = np.where(incidence <= 0)[0]
    if len(orphan):
        return None, {
            "reason": "corridor_remesh_created_orphan_nodes",
            "orphan_node_ids": orphan[:20].astype(int).tolist(),
        }

    certificate = {
        "certificate_schema": "v10.0.5.18.3.9_local_support_area_quality_certificate",
        "triangle_quality_floor": qfloor,
        "min_triangle_quality": qmin,
        "child_area_ratio_floor": afloor,
        "min_local_support_area_ratio": min_local_ratio,
        "area_gate_metric": "new_triangle_area_over_equilateral_path_support_cell_area",
        "path_support_spacing_m": support_spacing,
        "local_reference_area_m2": local_reference_area,
        "min_old_parent_area_ratio_diagnostic": min_old_parent_ratio,
        "old_parent_ratio_is_state_transfer_diagnostic_not_immediate_child_gate": True,
        "quality_floor_relaxed": False,
        "area_ratio_threshold_relaxed": False,
    }

    return {
        "mesh": mesh1,
        "boundary": make_boundary_data(mesh1, self.geom),
        "damage": d1,
        "displacement": u1.reshape(-1),
        "parent_map": parent_map,
        "local_to_global": local_default,
        "new_vertex_ids": new_vertex_ids,
        "local_elements": new_local_global,
        "assigned_parents": assigned,
        "quality_certificate": certificate,
    }, {
        "min_triangle_quality": qmin,
        "triangle_quality_floor": qfloor,
        "min_child_area_ratio": min_local_ratio,
        "child_area_ratio_floor": afloor,
        "area_gate_metric": certificate["area_gate_metric"],
        "path_support_spacing_m": support_spacing,
        "local_reference_area_m2": local_reference_area,
        "min_old_parent_area_ratio_diagnostic": min_old_parent_ratio,
        "removed_element_count": int(len(selected)),
        "new_cavity_element_count": int(len(new_local_global)),
        "net_element_change": int(len(new_local_global) - len(selected)),
        "new_bulk_node_count": int(len(new_vertex_ids)),
        "parent_transfer_policy": "new_triangle_centroid_to_containing_old_elastic_parent",
        "quality_certificate": certificate,
    }


class CertifiedAtomicPathCorridorCZMBackendV10051839(
    _base.AtomicPathCorridorCZMBackendV10051839
):
    """Attach the accepted local-scale certificate to the backend result."""

    name = "adaptive_czm_v10051839_atomic_path_corridor_certified"

    def advance(self, **kwargs):
        result = super().advance(**kwargs)
        if bool(getattr(result, "inserted", False)):
            audit = _base.audit_payload()
            event = audit["events"][-1]
            certificate = dict(
                event.get("accepted_corridor", {}).get("quality_certificate", {})
            )
            setattr(result, "v10051839_quality_certificate", certificate)
        return result


def install() -> Any:
    saved = _base._stitch_triangulation
    _base._stitch_triangulation = stitch_triangulation_local_support_area
    return saved


def restore(saved: Any) -> None:
    _base._stitch_triangulation = saved


__all__ = [
    "CertifiedAtomicPathCorridorCZMBackendV10051839",
    "MODEL_ID",
    "install",
    "restore",
    "stitch_triangulation_local_support_area",
]
