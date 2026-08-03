"""Exact-endpoint convex two-triangle cavity retriangulation for v10.0.5.18.3.8.

This fallback is used only after the v10.0.5.18.3.7 one-edge midpoint
operation family is exhausted for an edge-connected target.  The target parent
and the triangle across the endpoint-nearest edge form a four-vertex cavity.
The old internal diagonal is removed, the exact hazard-selected endpoint is
inserted as a bulk node, and the convex cavity is Delaunay-retriangulated.

The cavity boundary is unchanged.  The exact endpoint, direction, and event
length are unchanged.  All final cavity triangles must satisfy the retained
quality and immediate area floors.  No equal-length physical partition retry is
introduced.

The current production contract has elastic FEM bulk (tip-only MPZ state), so
no bulk constitutive Gauss-point history is active.  A deterministic
centroid-classified parent map is nevertheless emitted for the existing state
transfer interface and audited explicitly.
"""
from __future__ import annotations

from typing import Any

import numpy as np
from scipy.spatial import ConvexHull, Delaunay, QhullError

from .mesh import make_boundary_data, rebuild_tri_mesh
from . import quality_aware_czm_patch_v10051835 as _base
from . import quality_aware_czm_ray_handoff_v10051835 as _handoff
from . import quality_aware_czm_retry_v10051835 as _retry


MODEL_ID = "exact_endpoint_convex_two_triangle_cavity_v10_0_5_18_3_8"
_ORIGINAL_REFINE_ONCE = _handoff.refine_once


def _edge_key(i: int, j: int) -> tuple[int, int]:
    return (int(min(i, j)), int(max(i, j)))


def _triangle_contains(nodes: np.ndarray, conn: np.ndarray, point: np.ndarray) -> tuple[bool, float]:
    weights = _base._barycentric(
        np.asarray(nodes[np.asarray(conn, dtype=int)], dtype=float),
        np.asarray(point, dtype=float),
    )
    if weights is None:
        return False, -np.inf
    margin = float(np.min(weights))
    return bool(margin >= -1.0e-10), margin


def _target_parent_and_nearest_edge(mesh: Any, target: np.ndarray):
    parent = _handoff._target_triangle_any(mesh, np.asarray(target, dtype=float))
    if parent is None:
        return None
    conn = np.asarray(mesh.elems[int(parent)], dtype=int)
    tri = np.asarray(mesh.nodes[conn], dtype=float)
    weights = _base._barycentric(tri, np.asarray(target, dtype=float))
    if weights is None or float(np.min(weights)) < -1.0e-10:
        return None
    opposite_local = int(np.argmin(weights))
    edge_nodes = [int(conn[i]) for i in range(3) if i != opposite_local]
    return int(parent), conn, np.asarray(weights, dtype=float), _edge_key(*edge_nodes)


def _adjacent_across_edge(mesh: Any, parent: int, edge: tuple[int, int]) -> int | None:
    elems = np.asarray(mesh.elems, dtype=int)
    i, j = edge
    parents = np.where(
        np.any(elems == int(i), axis=1) & np.any(elems == int(j), axis=1)
    )[0]
    other = [int(e) for e in parents if int(e) != int(parent)]
    return other[0] if len(other) == 1 else None


def _boundary_edges(elems: np.ndarray) -> set[tuple[int, int]]:
    counts: dict[tuple[int, int], int] = {}
    for conn in np.asarray(elems, dtype=int):
        for edge in (
            _edge_key(conn[0], conn[1]),
            _edge_key(conn[1], conn[2]),
            _edge_key(conn[2], conn[0]),
        ):
            counts[edge] = counts.get(edge, 0) + 1
    return {edge for edge, count in counts.items() if count == 1}


def _oriented(self: Any, nodes: np.ndarray, conn: np.ndarray, reference_sign: float) -> np.ndarray:
    out = np.asarray(conn, dtype=int).copy()
    sign = float(self._signed_twice_area(nodes, out.reshape(1, 3))[0])
    if reference_sign * sign < 0.0:
        out[0], out[1] = out[1], out[0]
    return out


def _assign_old_parent(
    nodes: np.ndarray,
    new_conn: np.ndarray,
    old_elems: np.ndarray,
    old_ids: tuple[int, int],
) -> int:
    centroid = np.mean(nodes[np.asarray(new_conn, dtype=int)], axis=0)
    ranked: list[tuple[int, float, int]] = []
    for old_id in old_ids:
        inside, margin = _triangle_contains(nodes, old_elems[int(old_id)], centroid)
        ranked.append((int(inside), float(margin), -int(old_id)))
    best = max(range(len(old_ids)), key=lambda k: ranked[k])
    return int(old_ids[best])


def _two_triangle_cavity_candidate(
    self: Any,
    state: _base.State,
    p0: np.ndarray,
    target: np.ndarray,
    front_id: int,
    level: int,
):
    mesh0 = state.mesh
    nodes0 = np.asarray(mesh0.nodes, dtype=float)
    elems0 = np.asarray(mesh0.elems, dtype=int)
    target = np.asarray(target, dtype=float)
    p0 = np.asarray(p0, dtype=float)

    located = _target_parent_and_nearest_edge(mesh0, target)
    if located is None:
        return None, {
            "level": int(level),
            "accepted": False,
            "reason": "two_triangle_cavity_target_parent_not_found",
        }
    parent, parent_conn, target_weights, near_edge = located
    adjacent = _adjacent_across_edge(mesh0, parent, near_edge)
    if adjacent is None:
        return None, {
            "level": int(level),
            "accepted": False,
            "reason": "endpoint_nearest_edge_has_no_unique_adjacent_triangle",
            "target_parent_element": int(parent),
            "endpoint_nearest_edge": list(near_edge),
        }

    old_ids = (int(parent), int(adjacent))
    old_local = np.asarray(elems0[list(old_ids)], dtype=int)
    cavity_vertices = np.unique(old_local.ravel())
    if len(cavity_vertices) != 4:
        return None, {
            "level": int(level),
            "accepted": False,
            "reason": "two_triangle_cavity_requires_four_unique_boundary_vertices",
            "target_parent_element": int(parent),
            "adjacent_element": int(adjacent),
            "cavity_vertex_count": int(len(cavity_vertices)),
        }

    boundary_before = _boundary_edges(old_local)
    local_xy = np.asarray(nodes0[cavity_vertices], dtype=float)
    try:
        hull = ConvexHull(local_xy)
    except QhullError as exc:
        return None, {
            "level": int(level),
            "accepted": False,
            "reason": f"two_triangle_cavity_convex_hull_error:{exc}",
        }
    if len(hull.vertices) != 4:
        return None, {
            "level": int(level),
            "accepted": False,
            "reason": "two_triangle_cavity_is_not_strictly_convex",
            "convex_hull_vertex_count": int(len(hull.vertices)),
        }

    local_points = np.vstack([local_xy, target[None, :]])
    target_local = 4
    try:
        tri = Delaunay(local_points)
    except QhullError as exc:
        return None, {
            "level": int(level),
            "accepted": False,
            "reason": f"two_triangle_cavity_delaunay_error:{exc}",
        }

    target_id = int(len(nodes0))
    local_to_global = np.concatenate(
        [np.asarray(cavity_vertices, dtype=int), np.asarray([target_id], dtype=int)]
    )
    local_simplices = np.asarray(tri.simplices, dtype=int)
    if len(local_simplices) < 3:
        return None, {
            "level": int(level),
            "accepted": False,
            "reason": "two_triangle_cavity_has_too_few_retriangulated_elements",
            "new_local_element_count": int(len(local_simplices)),
        }
    if not np.all(np.any(local_simplices == target_local, axis=1)):
        return None, {
            "level": int(level),
            "accepted": False,
            "reason": "exact_endpoint_not_connected_to_every_cavity_sector",
        }

    nodes1 = np.vstack([nodes0, target[None, :]])
    reference_sign = float(self._signed_twice_area(nodes0, old_local[[0]])[0])
    local_global = np.asarray(
        [
            _oriented(
                self,
                nodes1,
                local_to_global[np.asarray(conn, dtype=int)],
                reference_sign,
            )
            for conn in local_simplices
        ],
        dtype=int,
    )

    boundary_after = _boundary_edges(local_global)
    if boundary_after != boundary_before:
        return None, {
            "level": int(level),
            "accepted": False,
            "reason": "two_triangle_cavity_boundary_not_preserved",
            "boundary_before": [list(edge) for edge in sorted(boundary_before)],
            "boundary_after": [list(edge) for edge in sorted(boundary_after)],
        }

    tip_ids = set(map(int, self._tip_geometric_node_ids(mesh0, p0, int(front_id))))
    target_neighbors: set[int] = set()
    for conn in local_global:
        if target_id in conn:
            target_neighbors.update(int(v) for v in conn if int(v) != target_id)
    connected_tip_ids = sorted(tip_ids.intersection(target_neighbors))
    if not connected_tip_ids:
        return None, {
            "level": int(level),
            "accepted": False,
            "reason": "retriangulated_exact_endpoint_not_connected_to_tip_copy",
            "target_neighbor_ids": sorted(target_neighbors),
            "tip_geometric_node_ids": sorted(tip_ids),
        }

    qfloor = _base.float_env(
        "ARRHENIUS_MIN_ACCEPTED_TRIANGLE_QUALITY",
        float(self.min_triangle_quality),
    )
    afloor = _base.float_env(
        "ARRHENIUS_MIN_ACCEPTED_CHILD_AREA_RATIO",
        float(self.min_area_ratio),
    )
    quality = self._triangle_quality(nodes1, local_global)
    qmin = float(np.min(quality)) if quality.size else 1.0

    assigned = np.asarray(
        [
            _assign_old_parent(nodes1, conn, elems0, old_ids)
            for conn in local_global
        ],
        dtype=int,
    )
    new_area_twice = np.abs(self._signed_twice_area(nodes1, local_global))
    old_area_twice = np.abs(self._signed_twice_area(nodes0, elems0[assigned]))
    area_ratios = new_area_twice / np.maximum(old_area_twice, 1.0e-300)
    amin = float(np.min(area_ratios)) if area_ratios.size else 1.0
    issues = _base.floor_issues(qmin, amin, qfloor, afloor)
    if issues:
        return None, {
            "level": int(level),
            "accepted": False,
            "reason": ";".join(issues),
            "refinement_kind": "exact_endpoint_convex_two_triangle_delaunay",
            "target_parent_element": int(parent),
            "adjacent_element": int(adjacent),
            "endpoint_nearest_edge": list(near_edge),
            "target_barycentric_weights_before": target_weights.tolist(),
            "min_triangle_quality": qmin,
            "triangle_quality_floor": qfloor,
            "min_immediate_child_area_ratio": amin,
            "child_area_ratio_floor": afloor,
            "new_triangle_parent_ids": assigned.tolist(),
        }

    ordered_local = sorted(
        range(len(local_global)),
        key=lambda idx: (
            int(assigned[idx]),
            tuple(sorted(map(int, local_global[idx]))),
        ),
    )
    local_global = local_global[np.asarray(ordered_local, dtype=int)]
    assigned = assigned[np.asarray(ordered_local, dtype=int)]
    area_ratios = area_ratios[np.asarray(ordered_local, dtype=int)]

    elems1 = elems0.copy()
    replacement_ids = sorted(old_ids)
    elems1[replacement_ids[0]] = local_global[0]
    elems1[replacement_ids[1]] = local_global[1]
    appended = local_global[2:]
    if len(appended):
        elems1 = np.vstack([elems1, appended])
    mesh1 = rebuild_tri_mesh(
        nodes1,
        elems1,
        tip_centers=[p0, target],
    )

    u0 = np.asarray(state.displacement, dtype=float).reshape(-1, 2)
    u_target = np.asarray(target_weights, dtype=float) @ u0[parent_conn]
    u1 = np.vstack([u0, u_target[None, :]]).reshape(-1)
    d0 = np.asarray(state.damage, dtype=float)
    d_target = float(np.asarray(target_weights, dtype=float) @ d0[parent_conn])
    d1 = np.concatenate([d0, np.asarray([d_target], dtype=float)])

    elem_parent = np.arange(int(mesh0.ne), dtype=int)
    elem_parent[replacement_ids[0]] = int(assigned[0])
    elem_parent[replacement_ids[1]] = int(assigned[1])
    if len(appended):
        elem_parent = np.concatenate([elem_parent, assigned[2:]])

    next_state = _base.State(
        mesh=mesh1,
        boundary=make_boundary_data(mesh1, self.geom),
        damage=d1,
        displacement=u1,
        parent_to_root=state.parent_to_root[elem_parent],
    )
    record = {
        "level": int(level),
        "accepted": True,
        "refinement_kind": "exact_endpoint_convex_two_triangle_delaunay",
        "target_parent_element": int(parent),
        "adjacent_element": int(adjacent),
        "removed_internal_edge": list(near_edge),
        "cavity_boundary_edges": [list(edge) for edge in sorted(boundary_before)],
        "cavity_vertex_ids": list(map(int, cavity_vertices)),
        "exact_endpoint_node_id": int(target_id),
        "exact_endpoint_m": target.tolist(),
        "exact_endpoint_connected_tip_ids": connected_tip_ids,
        "target_barycentric_weights_before": target_weights.tolist(),
        "min_triangle_quality": qmin,
        "triangle_quality_floor": qfloor,
        "min_immediate_child_area_ratio": amin,
        "child_area_ratio_floor": afloor,
        "immediate_parent_area_ratio_enforced": True,
        "n_removed_elements": 2,
        "n_new_cavity_elements": int(len(local_global)),
        "n_new_nodes": 1,
        "n_net_new_elements": int(len(local_global) - 2),
        "new_triangle_parent_ids": assigned.tolist(),
        "new_triangle_parent_area_ratios": area_ratios.tolist(),
        "parent_transfer_policy": "new_triangle_centroid_classified_to_old_elastic_bulk_parent",
        "bulk_constitutive_state_required": False,
        "exact_endpoint_preinserted_as_bulk_node": True,
        "backend_zero_motion_neighbor_commit_expected": True,
    }
    return next_state, record


def refine_once(
    self: Any,
    state: _base.State,
    p0: np.ndarray,
    target: np.ndarray,
    front_id: int,
    level: int,
):
    refined, record = _ORIGINAL_REFINE_ONCE(
        self,
        state,
        np.asarray(p0, dtype=float),
        np.asarray(target, dtype=float),
        int(front_id),
        int(level),
    )
    if refined is not None:
        return refined, record

    cavity = record.get("target_cavity_failure")
    cavity_reason = None if not isinstance(cavity, dict) else cavity.get("reason")
    if (
        str(record.get("reason")) == "no_target_tip_triangle"
        and str(cavity_reason) == "no_quality_safe_adjacent_target_cavity_bisection"
    ):
        retriangulated, cavity_record = _two_triangle_cavity_candidate(
            self,
            state,
            np.asarray(p0, dtype=float),
            np.asarray(target, dtype=float),
            int(front_id),
            int(level),
        )
        if retriangulated is not None:
            return retriangulated, cavity_record
        merged = dict(record)
        merged["two_triangle_cavity_failure"] = cavity_record
        return None, merged
    return refined, record


def install() -> None:
    _handoff.install()
    _retry._refine_once = refine_once


__all__ = [
    "MODEL_ID",
    "_two_triangle_cavity_candidate",
    "install",
    "refine_once",
]
