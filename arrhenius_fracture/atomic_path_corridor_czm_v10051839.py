"""Atomic target-aware constrained path-corridor remeshing for v10.0.5.18.3.9.

A completed Arrhenius cleavage renewal is one physical event with one exact
start point, endpoint, direction, and total length.  This backend does not
march through whatever edges happen to exist in the inherited mesh.  Instead,
it removes a bounded connected set of forward bulk triangles, reconstructs
that cavity with a Planar Straight Line Graph (PSLG), and embeds the complete
exact physical segment as a constrained path before any cohesive topology is
committed.

The constrained path may contain numerical support vertices, but those vertices
are topology support for one atomic physical event; they are not independent
hazard events and do not alter the selected endpoint or total length.  The
whole remesh plus every cohesive subedge is committed or rolled back together.

All inherited cavity boundary vertices are retained.  Every final triangle is
checked against the unchanged production quality and mapped-parent area floors.
The production bulk is elastic, but an explicit new-element -> old-element map
is still emitted for the existing state-transfer contract.
"""
from __future__ import annotations

from collections import defaultdict, deque
import copy
import math
import os
from typing import Any

import numpy as np
import triangle as _triangle

from . import crack_backend as _cb
from .mesh import make_boundary_data, rebuild_tri_mesh


MODEL_ID = "atomic_target_aware_path_corridor_czm_v10_0_5_18_3_9"
SCHEMA = "v10.0.5.18.3.9_atomic_pslg_path_corridor"

_AUDIT: dict[str, Any] = {
    "schema": SCHEMA,
    "model_id": MODEL_ID,
    "events": [],
}


def reset_audit() -> None:
    _AUDIT["events"] = []


def _float_env(name: str, default: float) -> float:
    try:
        value = float(os.environ.get(name, default))
    except (TypeError, ValueError):
        value = float(default)
    return value if math.isfinite(value) else float(default)


def _int_env(name: str, default: int) -> int:
    try:
        value = int(os.environ.get(name, default))
    except (TypeError, ValueError):
        value = int(default)
    return max(value, 1)


def _coord_key(point: np.ndarray, tol: float) -> tuple[int, int]:
    return tuple(np.round(np.asarray(point, float) / tol).astype(np.int64))


def _edge_key(a: tuple[int, int], b: tuple[int, int]):
    return tuple(sorted((a, b)))


def _segment_coordinates(points: np.ndarray, p0: np.ndarray, direction: np.ndarray):
    rel = np.asarray(points, float) - np.asarray(p0, float)[None, :]
    normal = np.array([-direction[1], direction[0]], dtype=float)
    return rel @ direction, rel @ normal


def _barycentric(triangle: np.ndarray, point: np.ndarray) -> np.ndarray | None:
    tri = np.asarray(triangle, float)
    p = np.asarray(point, float)
    matrix = np.array(
        [
            [tri[0, 0], tri[1, 0], tri[2, 0]],
            [tri[0, 1], tri[1, 1], tri[2, 1]],
            [1.0, 1.0, 1.0],
        ],
        dtype=float,
    )
    try:
        return np.linalg.solve(matrix, np.array([p[0], p[1], 1.0]))
    except np.linalg.LinAlgError:
        return None


def _point_in_triangle(triangle: np.ndarray, point: np.ndarray, tol: float = 1.0e-10):
    weights = _barycentric(triangle, point)
    if weights is None:
        return False, -np.inf, None
    margin = float(np.min(weights))
    return bool(margin >= -tol and float(np.max(weights)) <= 1.0 + tol), margin, weights


def _unique_edge_lengths(mesh: Any, element_ids: np.ndarray) -> np.ndarray:
    edges: set[tuple[int, int]] = set()
    for conn in np.asarray(mesh.elems[np.asarray(element_ids, dtype=int)], dtype=int):
        for i, j in ((conn[0], conn[1]), (conn[1], conn[2]), (conn[2], conn[0])):
            edges.add((int(min(i, j)), int(max(i, j))))
    if not edges:
        return np.zeros(0, dtype=float)
    pairs = np.asarray(sorted(edges), dtype=int)
    return np.linalg.norm(mesh.nodes[pairs[:, 1]] - mesh.nodes[pairs[:, 0]], axis=1)


def _element_edge_adjacency(elems: np.ndarray) -> list[set[int]]:
    owners: dict[tuple[int, int], list[int]] = defaultdict(list)
    for e, conn in enumerate(np.asarray(elems, dtype=int)):
        for i, j in ((conn[0], conn[1]), (conn[1], conn[2]), (conn[2], conn[0])):
            owners[(int(min(i, j)), int(max(i, j)))].append(int(e))
    adjacency = [set() for _ in range(len(elems))]
    for rows in owners.values():
        if len(rows) == 2:
            a, b = rows
            adjacency[a].add(b)
            adjacency[b].add(a)
    return adjacency


def _expand_elements(seed: set[int], adjacency: list[set[int]], rings: int) -> set[int]:
    selected = set(map(int, seed))
    frontier = set(selected)
    for _ in range(max(int(rings), 0)):
        nxt: set[int] = set()
        for e in frontier:
            nxt.update(adjacency[int(e)])
        nxt.difference_update(selected)
        if not nxt:
            break
        selected.update(nxt)
        frontier = nxt
    return selected


def _locate_parent(mesh: Any, point: np.ndarray, candidates: np.ndarray | None = None):
    ids = np.arange(int(mesh.ne), dtype=int) if candidates is None else np.asarray(candidates, dtype=int)
    best = None
    for e in ids:
        inside, margin, weights = _point_in_triangle(mesh.nodes[mesh.elems[int(e)]], point)
        if inside and (best is None or margin > best[0]):
            best = (float(margin), int(e), np.asarray(weights, float))
    return best


def _selected_corridor_elements(
    self: Any,
    mesh: Any,
    p0: np.ndarray,
    p1: np.ndarray,
    direction: np.ndarray,
    front_id: int,
    width: float,
    back: float,
    forward_pad: float,
    rings: int,
) -> tuple[np.ndarray | None, dict[str, Any]]:
    nodes = np.asarray(mesh.nodes, float)
    elems = np.asarray(mesh.elems, int)
    X = nodes[elems]
    s, r = _segment_coordinates(X.reshape(-1, 2), p0, direction)
    s = s.reshape(len(elems), 3)
    r = r.reshape(len(elems), 3)
    length = float(np.linalg.norm(p1 - p0))

    mask = (
        (s.max(axis=1) >= -float(back))
        & (s.min(axis=1) <= length + float(forward_pad))
        & (r.max(axis=1) >= -float(width))
        & (r.min(axis=1) <= float(width))
    )
    seed = set(np.where(mask)[0].astype(int).tolist())

    tip_ids = self._tip_geometric_node_ids(mesh, p0, int(front_id))
    if not tip_ids:
        return None, {"reason": "corridor_tip_node_not_found"}
    incident = self._incident_elements(elems, set(map(int, tip_ids)))
    cent = X.mean(axis=1)
    cent_s, _ = _segment_coordinates(cent, p0, direction)
    seed.update(int(e) for e in incident if float(cent_s[int(e)]) >= -0.35 * float(back))

    target = _locate_parent(mesh, p1)
    if target is None:
        return None, {"reason": "corridor_target_parent_not_found"}
    seed.add(int(target[1]))

    selected = _expand_elements(seed, _element_edge_adjacency(elems), rings)
    # Never consume a wake element whose centroid is materially behind the
    # current physical tip.  This keeps the existing cohesive crack boundary
    # outside the remeshed cavity.
    selected = {e for e in selected if float(cent_s[int(e)]) >= -0.75 * float(back)}
    selected.update(int(e) for e in incident if float(cent_s[int(e)]) >= -0.35 * float(back))
    selected.add(int(target[1]))
    if not selected:
        return None, {"reason": "empty_path_corridor"}

    ids = np.asarray(sorted(selected), dtype=int)
    return ids, {
        "selected_element_count": int(len(ids)),
        "target_parent_element": int(target[1]),
        "tip_incident_selected_count": int(sum(int(e) in selected for e in incident)),
        "corridor_width_m": float(width),
        "corridor_back_m": float(back),
        "corridor_forward_pad_m": float(forward_pad),
        "corridor_adjacency_rings": int(rings),
    }


def _geometric_boundary(
    mesh: Any,
    selected: np.ndarray,
    tol: float,
) -> tuple[dict[str, Any] | None, dict[str, Any]]:
    nodes = np.asarray(mesh.nodes, float)
    elems = np.asarray(mesh.elems[np.asarray(selected, dtype=int)], int)
    coord_ids: dict[tuple[int, int], list[int]] = defaultdict(list)
    edge_rows: dict[tuple[tuple[int, int], tuple[int, int]], list[tuple[int, int]]] = defaultdict(list)

    for conn in elems:
        for nid in conn:
            key = _coord_key(nodes[int(nid)], tol)
            coord_ids[key].append(int(nid))
        for i, j in ((conn[0], conn[1]), (conn[1], conn[2]), (conn[2], conn[0])):
            ki = _coord_key(nodes[int(i)], tol)
            kj = _coord_key(nodes[int(j)], tol)
            if ki == kj:
                continue
            edge_rows[_edge_key(ki, kj)].append((int(i), int(j)))

    boundary_keys = [edge for edge, rows in edge_rows.items() if len(rows) == 1]
    if len(boundary_keys) < 3:
        return None, {"reason": "corridor_boundary_has_too_few_edges"}

    graph: dict[tuple[int, int], set[tuple[int, int]]] = defaultdict(set)
    for ka, kb in boundary_keys:
        graph[ka].add(kb)
        graph[kb].add(ka)
    degrees = {key: len(value) for key, value in graph.items()}
    if any(value != 2 for value in degrees.values()):
        return None, {
            "reason": "corridor_boundary_not_closed_degree_two",
            "boundary_degree_histogram": {
                str(k): int(sum(v == k for v in degrees.values()))
                for k in sorted(set(degrees.values()))
            },
        }

    components: list[set[tuple[int, int]]] = []
    remaining = set(graph)
    while remaining:
        root = next(iter(remaining))
        comp = {root}
        queue = [root]
        while queue:
            key = queue.pop()
            for other in graph[key]:
                if other not in comp:
                    comp.add(other)
                    queue.append(other)
        components.append(comp)
        remaining.difference_update(comp)
    if len(components) != 1:
        return None, {
            "reason": "corridor_boundary_has_multiple_loops",
            "boundary_component_sizes": [int(len(c)) for c in components],
        }

    unique_keys = sorted(coord_ids)
    key_to_local = {key: index for index, key in enumerate(unique_keys)}
    vertices = np.asarray(
        [np.mean(nodes[np.asarray(coord_ids[key], dtype=int)], axis=0) for key in unique_keys],
        dtype=float,
    )
    segments = np.asarray(
        [[key_to_local[a], key_to_local[b]] for a, b in boundary_keys],
        dtype=int,
    )
    boundary_node_choices = {
        key: sorted(set(map(int, ids))) for key, ids in coord_ids.items()
    }
    boundary_edge_original_ids = {
        edge: tuple(map(int, rows[0])) for edge, rows in edge_rows.items() if len(rows) == 1
    }
    return {
        "vertices": vertices,
        "segments": segments,
        "keys": unique_keys,
        "key_to_local": key_to_local,
        "coord_ids": boundary_node_choices,
        "boundary_keys": boundary_keys,
        "boundary_edge_original_ids": boundary_edge_original_ids,
    }, {
        "boundary_vertex_count": int(len(vertices)),
        "boundary_edge_count": int(len(segments)),
        "boundary_component_count": 1,
    }


def _append_unique_vertex(vertices: list[np.ndarray], point: np.ndarray, tol: float) -> int:
    p = np.asarray(point, float)
    for index, value in enumerate(vertices):
        if float(np.linalg.norm(np.asarray(value) - p)) <= tol:
            return int(index)
    vertices.append(p.copy())
    return int(len(vertices) - 1)


def _input_pslg(
    mesh: Any,
    selected: np.ndarray,
    boundary: dict[str, Any],
    p0: np.ndarray,
    p1: np.ndarray,
    direction: np.ndarray,
    h_local: float,
    tol: float,
    path_subsegments: int,
):
    nodes = np.asarray(mesh.nodes, float)
    vertices = [np.asarray(value, float).copy() for value in boundary["vertices"]]
    segments = [tuple(map(int, value)) for value in np.asarray(boundary["segments"], int)]
    markers = [1 for _ in segments]

    # Preserve every inherited cavity coordinate as a PSLG vertex, including
    # interior nodes.  This prevents orphaned old bulk nodes without compacting
    # or renumbering the global cohesive network.
    selected_nodes = sorted(set(map(int, np.asarray(mesh.elems[selected], int).ravel())))
    for nid in selected_nodes:
        _append_unique_vertex(vertices, nodes[int(nid)], tol)

    nseg = max(int(path_subsegments), 1)
    path_ids: list[int] = []
    for i in range(nseg + 1):
        point = np.asarray(p0, float) + (float(i) / float(nseg)) * (np.asarray(p1) - np.asarray(p0))
        path_ids.append(_append_unique_vertex(vertices, point, tol))
    for a, b in zip(path_ids[:-1], path_ids[1:]):
        if int(a) == int(b):
            return None, {"reason": "collapsed_input_path_subsegment"}
        segments.append((int(a), int(b)))
        markers.append(2)

    return {
        "vertices": np.asarray(vertices, dtype=float),
        "segments": np.asarray(segments, dtype=int),
        "segment_markers": np.asarray(markers, dtype=int).reshape(-1, 1),
        "path_input_ids": path_ids,
        "selected_old_node_ids": selected_nodes,
    }, {
        "input_vertex_count": int(len(vertices)),
        "input_segment_count": int(len(segments)),
        "input_path_subsegment_count": int(nseg),
        "input_path_support_spacing_m": float(np.linalg.norm(p1 - p0) / nseg),
        "input_path_support_spacing_over_h_local": float(np.linalg.norm(p1 - p0) / max(nseg * h_local, 1.0e-300)),
    }


def _nearest_vertex_id(vertices: np.ndarray, point: np.ndarray, tol: float) -> int | None:
    dist = np.linalg.norm(np.asarray(vertices, float) - np.asarray(point, float)[None, :], axis=1)
    index = int(np.argmin(dist))
    return index if float(dist[index]) <= tol else None


def _path_chain_from_output(
    vertices: np.ndarray,
    edges: np.ndarray,
    p0: np.ndarray,
    p1: np.ndarray,
    direction: np.ndarray,
    tol: float,
):
    start = _nearest_vertex_id(vertices, p0, tol)
    end = _nearest_vertex_id(vertices, p1, tol)
    if start is None or end is None:
        return None, {"reason": "exact_path_endpoint_missing_after_triangulation"}

    length = float(np.linalg.norm(p1 - p0))
    s, r = _segment_coordinates(vertices, p0, direction)
    on_line = (
        (np.abs(r) <= tol)
        & (s >= -tol)
        & (s <= length + tol)
    )
    graph: dict[int, set[int]] = defaultdict(set)
    for a, b in np.asarray(edges, dtype=int):
        a = int(a); b = int(b)
        if not (bool(on_line[a]) and bool(on_line[b])):
            continue
        if abs(float(s[a] - s[b])) <= tol:
            continue
        graph[a].add(b)
        graph[b].add(a)

    queue = deque([int(start)])
    parent = {int(start): None}
    while queue:
        current = queue.popleft()
        if current == int(end):
            break
        neighbors = sorted(
            graph.get(current, set()),
            key=lambda idx: (float(s[idx]), int(idx)),
        )
        for nxt in neighbors:
            if nxt in parent:
                continue
            # The exact path is monotone; reject backward graph branches.
            if float(s[nxt]) < float(s[current]) - tol:
                continue
            parent[nxt] = current
            queue.append(nxt)
    if int(end) not in parent:
        return None, {
            "reason": "constrained_exact_path_not_recovered_as_mesh_edges",
            "on_line_vertex_count": int(np.sum(on_line)),
        }

    chain = []
    current: int | None = int(end)
    while current is not None:
        chain.append(int(current))
        current = parent[current]
    chain.reverse()
    if len(chain) < 2:
        return None, {"reason": "recovered_path_chain_has_no_edge"}
    projections = np.asarray([s[idx] for idx in chain], dtype=float)
    if np.any(np.diff(projections) <= tol):
        return None, {"reason": "recovered_path_chain_is_not_strictly_monotone"}
    return chain, {
        "output_path_vertex_count": int(len(chain)),
        "output_path_subsegment_count": int(len(chain) - 1),
        "output_path_subsegment_lengths_m": [
            float(np.linalg.norm(vertices[b] - vertices[a]))
            for a, b in zip(chain[:-1], chain[1:])
        ],
    }


def _old_coordinate_map(mesh: Any, tol: float) -> dict[tuple[int, int], list[int]]:
    out: dict[tuple[int, int], list[int]] = defaultdict(list)
    for nid, point in enumerate(np.asarray(mesh.nodes, float)):
        out[_coord_key(point, tol)].append(int(nid))
    return out


def _interpolate_scalar(mesh: Any, values: np.ndarray, point: np.ndarray) -> float:
    p = np.asarray(point, float)
    vals = np.asarray(values, float)
    X = mesh.nodes[mesh.elems]
    eps = max(1.0e-14, 1.0e-9 * max(float(mesh.hbar_tip), float(mesh.hbar), 1.0e-12))
    mask = (
        (X[:, :, 0].min(axis=1) - eps <= p[0])
        & (X[:, :, 0].max(axis=1) + eps >= p[0])
        & (X[:, :, 1].min(axis=1) - eps <= p[1])
        & (X[:, :, 1].max(axis=1) + eps >= p[1])
    )
    for e in np.where(mask)[0]:
        inside, _, weights = _point_in_triangle(X[int(e)], p, tol=1.0e-8)
        if inside and weights is not None:
            return float(np.asarray(weights) @ vals[mesh.elems[int(e)]])
    nearest = int(np.argmin(np.linalg.norm(mesh.nodes - p[None, :], axis=1)))
    return float(vals[nearest])


def _assign_old_parent(mesh: Any, selected: np.ndarray, triangle_xy: np.ndarray) -> int:
    centroid = np.mean(np.asarray(triangle_xy, float), axis=0)
    best = None
    for e in np.asarray(selected, dtype=int):
        inside, margin, _ = _point_in_triangle(mesh.nodes[mesh.elems[int(e)]], centroid)
        if inside and (best is None or float(margin) > best[0]):
            best = (float(margin), int(e))
    if best is not None:
        return int(best[1])
    old_centroids = mesh.nodes[mesh.elems[selected]].mean(axis=1)
    return int(selected[int(np.argmin(np.linalg.norm(old_centroids - centroid[None, :], axis=1)))])


def _stitch_triangulation(
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
    old_map = _old_coordinate_map(mesh, tol)
    tip_ids = self._tip_geometric_node_ids(mesh, p0, int(front_id))
    if not tip_ids:
        return None, {"reason": "stitch_tip_node_not_found"}
    tip_plus, tip_minus = (
        (int(tip_ids[0]), int(tip_ids[1]))
        if len(tip_ids) >= 2 else (int(tip_ids[0]), int(tip_ids[0]))
    )
    p0_key = _coord_key(p0, tol)

    global_nodes = old_nodes.tolist()
    local_default: dict[int, int] = {}
    new_vertex_ids: list[int] = []
    for lid, point in enumerate(vertices):
        key = _coord_key(point, tol)
        existing = old_map.get(key, [])
        if existing:
            if key == p0_key:
                local_default[int(lid)] = tip_plus
            else:
                local_default[int(lid)] = int(existing[0])
        else:
            gid = int(len(global_nodes))
            global_nodes.append(np.asarray(point, float).tolist())
            local_default[int(lid)] = gid
            new_vertex_ids.append(gid)

    global_nodes_array = np.asarray(global_nodes, float)
    new_local_global: list[np.ndarray] = []
    reference_sign = float(np.median(np.sign(self._signed_twice_area(mesh.nodes, mesh.elems[selected]))))
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
            key = _coord_key(vertices[int(lid)], tol)
            if key == p0_key:
                mapped.append(tip_plus if cross >= 0.0 else tip_minus)
            else:
                mapped.append(local_default[int(lid)])
        mapped = np.asarray(mapped, dtype=int)
        sign = float(self._signed_twice_area(global_nodes_array, mapped.reshape(1, 3))[0])
        if sign * reference_sign < 0.0:
            mapped[0], mapped[1] = mapped[1], mapped[0]
        new_local_global.append(mapped)
    new_local_global = np.asarray(new_local_global, dtype=int)

    q = self._triangle_quality(global_nodes_array, new_local_global)
    qmin = float(np.min(q)) if len(q) else 1.0
    assigned = np.asarray(
        [
            _assign_old_parent(mesh, selected, global_nodes_array[conn])
            for conn in new_local_global
        ],
        dtype=int,
    )
    new_area = 0.5 * np.abs(self._signed_twice_area(global_nodes_array, new_local_global))
    old_area = np.asarray(mesh.area_e[assigned], float)
    ratios = new_area / np.maximum(old_area, 1.0e-300)
    amin = float(np.min(ratios)) if len(ratios) else 1.0
    qfloor = _float_env("ARRHENIUS_MIN_ACCEPTED_TRIANGLE_QUALITY", self.min_triangle_quality)
    afloor = _float_env("ARRHENIUS_MIN_ACCEPTED_CHILD_AREA_RATIO", self.min_area_ratio)
    issues = []
    if qmin < qfloor:
        issues.append(f"triangle_quality={qmin:.6e}<{qfloor:.6e}")
    if amin < afloor:
        issues.append(f"child_area_ratio={amin:.6e}<{afloor:.6e}")
    if issues:
        return None, {
            "reason": ";".join(issues),
            "min_triangle_quality": qmin,
            "triangle_quality_floor": qfloor,
            "min_child_area_ratio": amin,
            "child_area_ratio_floor": afloor,
            "new_triangle_count": int(len(new_local_global)),
            "new_vertex_count": int(len(new_vertex_ids)),
        }

    keep = np.setdiff1d(np.arange(int(mesh.ne), dtype=int), np.asarray(selected, dtype=int))
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
        d1[int(gid)] = _interpolate_scalar(mesh, d0, point)

    incidence = np.bincount(elems1.ravel(), minlength=len(global_nodes_array))
    orphan = np.where(incidence <= 0)[0]
    if len(orphan):
        return None, {
            "reason": "corridor_remesh_created_orphan_nodes",
            "orphan_node_ids": orphan[:20].astype(int).tolist(),
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
    }, {
        "min_triangle_quality": qmin,
        "triangle_quality_floor": qfloor,
        "min_child_area_ratio": amin,
        "child_area_ratio_floor": afloor,
        "removed_element_count": int(len(selected)),
        "new_cavity_element_count": int(len(new_local_global)),
        "net_element_change": int(len(new_local_global) - len(selected)),
        "new_bulk_node_count": int(len(new_vertex_ids)),
        "parent_transfer_policy": "new_triangle_centroid_to_containing_old_elastic_parent",
    }


def _global_path_chain(
    stitched: dict[str, Any],
    tri_out: dict[str, Any],
    local_chain: list[int],
    p0: np.ndarray,
    direction: np.ndarray,
    front_id: int,
    self: Any,
    tol: float,
):
    vertices = np.asarray(tri_out["vertices"], float)
    old_map = _old_coordinate_map(stitched["mesh"], tol)
    tip_ids = self._tip_geometric_node_ids(stitched["mesh"], p0, int(front_id))
    if not tip_ids:
        return None
    path: list[tuple[int, np.ndarray]] = [(int(tip_ids[0]), np.asarray(p0, float).copy())]
    for lid in local_chain[1:]:
        point = np.asarray(vertices[int(lid)], float)
        ids = old_map.get(_coord_key(point, tol), [])
        if not ids:
            return None
        path.append((int(ids[0]), point.copy()))
    return path


def _corridor_attempt(
    self: Any,
    mesh: Any,
    damage: np.ndarray,
    displacement: np.ndarray,
    p0: np.ndarray,
    p1: np.ndarray,
    direction: np.ndarray,
    front_id: int,
    width_factor: float,
    rings: int,
    quality_angle_deg: float,
    path_subsegments: int,
):
    tip_ids = self._tip_geometric_node_ids(mesh, p0, int(front_id))
    incident = self._incident_elements(mesh.elems, set(map(int, tip_ids))) if tip_ids else np.zeros(0, int)
    lengths = _unique_edge_lengths(mesh, incident)
    h_local = float(np.median(lengths)) if len(lengths) else float(max(mesh.hbar_tip, 1.0e-12))
    length = float(np.linalg.norm(p1 - p0))
    width = float(max(width_factor * h_local, 0.35 * length, 1.0e-12))
    back = float(max(0.35 * h_local, 1.0e-12))
    forward_pad = float(max(0.35 * h_local, 0.10 * length))
    tol = max(1.0e-12, 1.0e-8 * max(h_local, length, 1.0e-12))

    selected, selection_record = _selected_corridor_elements(
        self, mesh, p0, p1, direction, front_id,
        width, back, forward_pad, rings,
    )
    if selected is None:
        return None, {**selection_record, "stage": "select_corridor"}

    boundary, boundary_record = _geometric_boundary(mesh, selected, tol)
    if boundary is None:
        return None, {**selection_record, **boundary_record, "stage": "extract_boundary"}

    pslg, pslg_record = _input_pslg(
        mesh, selected, boundary, p0, p1, direction, h_local, tol, path_subsegments
    )
    if pslg is None:
        return None, {**selection_record, **boundary_record, **pslg_record, "stage": "build_pslg"}

    options = f"pYq{float(quality_angle_deg):.6g}neCQ"
    try:
        tri_out = _triangle.triangulate(
            {
                "vertices": pslg["vertices"],
                "segments": pslg["segments"],
                "segment_markers": pslg["segment_markers"],
            },
            options,
        )
    except Exception as exc:
        return None, {
            **selection_record,
            **boundary_record,
            **pslg_record,
            "stage": "triangle_pslg",
            "reason": f"triangle_exception:{type(exc).__name__}:{exc}",
            "triangle_options": options,
        }
    if "triangles" not in tri_out or "vertices" not in tri_out:
        return None, {
            **selection_record,
            **boundary_record,
            **pslg_record,
            "stage": "triangle_pslg",
            "reason": "triangle_returned_no_mesh",
            "triangle_output_keys": sorted(tri_out),
            "triangle_options": options,
        }
    edges = tri_out.get("edges")
    if edges is None:
        local_elems = np.asarray(tri_out["triangles"], int)
        edge_set = set()
        for conn in local_elems:
            for a, b in ((conn[0], conn[1]), (conn[1], conn[2]), (conn[2], conn[0])):
                edge_set.add((int(min(a, b)), int(max(a, b))))
        edges = np.asarray(sorted(edge_set), dtype=int)

    chain, chain_record = _path_chain_from_output(
        np.asarray(tri_out["vertices"], float), np.asarray(edges, int),
        p0, p1, direction, tol,
    )
    if chain is None:
        return None, {
            **selection_record,
            **boundary_record,
            **pslg_record,
            **chain_record,
            "stage": "recover_exact_path",
            "triangle_options": options,
        }

    stitched, stitch_record = _stitch_triangulation(
        self, mesh, damage, displacement, selected, tri_out,
        p0, p1, direction, front_id, tol,
    )
    if stitched is None:
        return None, {
            **selection_record,
            **boundary_record,
            **pslg_record,
            **chain_record,
            **stitch_record,
            "stage": "stitch_and_gate",
            "triangle_options": options,
        }

    path = _global_path_chain(stitched, tri_out, chain, p0, direction, front_id, self, tol)
    if path is None:
        return None, {
            **selection_record,
            **boundary_record,
            **pslg_record,
            **chain_record,
            **stitch_record,
            "stage": "map_exact_path",
            "reason": "failed_to_map_local_path_to_global_nodes",
        }

    record = {
        **selection_record,
        **boundary_record,
        **pslg_record,
        **chain_record,
        **stitch_record,
        "stage": "accepted",
        "accepted": True,
        "triangle_options": options,
        "quality_angle_deg": float(quality_angle_deg),
        "width_factor": float(width_factor),
        "h_local_m": h_local,
        "exact_endpoint_m": np.asarray(p1, float).tolist(),
        "exact_direction": np.asarray(direction, float).tolist(),
        "exact_requested_length_m": length,
        "numerical_path_subsegments_are_not_physical_events": True,
    }
    return {**stitched, "path": path}, record


class AtomicPathCorridorCZMBackendV10051839(_cb.AdaptiveCZMBackend):
    """Complete-event constrained-cavity backend for one physical renewal."""

    name = "adaptive_czm_v10051839_atomic_path_corridor"

    def _build_corridor(
        self,
        mesh: Any,
        damage: np.ndarray,
        displacement: np.ndarray,
        p0: np.ndarray,
        p1: np.ndarray,
        direction: np.ndarray,
        front_id: int,
    ):
        tip_ids = self._tip_geometric_node_ids(mesh, p0, int(front_id))
        incident = self._incident_elements(mesh.elems, set(map(int, tip_ids))) if tip_ids else np.zeros(0, int)
        lengths = _unique_edge_lengths(mesh, incident)
        h_local = float(np.median(lengths)) if len(lengths) else float(max(mesh.hbar_tip, 1.0e-12))
        length = float(np.linalg.norm(p1 - p0))
        base_n = max(1, int(math.ceil(length / max(0.75 * h_local, 1.0e-12))))
        max_path = _int_env("ARRHENIUS_CORRIDOR_MAX_PATH_SUBSEGMENTS", 16)
        n_values = []
        for value in (1, base_n, base_n + 1, 2 * base_n):
            value = min(max(int(value), 1), max_path)
            if value not in n_values:
                n_values.append(value)

        attempts = []
        width_values = (0.75, 1.0, 1.5, 2.25, 3.25)
        ring_values = (0, 1, 2)
        angle_values = (1.0, 2.0, 5.0)
        for nseg in n_values:
            for width in width_values:
                for rings in ring_values:
                    for angle in angle_values:
                        state, record = _corridor_attempt(
                            self, mesh, damage, displacement,
                            p0, p1, direction, front_id,
                            width, rings, angle, nseg,
                        )
                        record.update(
                            attempt_index=int(len(attempts)),
                            accepted=state is not None,
                        )
                        attempts.append(copy.deepcopy(record))
                        if state is not None:
                            return state, record, attempts
        return None, {
            "reason": "no_feasible_atomic_path_corridor",
            "attempt_count": int(len(attempts)),
        }, attempts

    @staticmethod
    def _extend_damage(old_mesh: Any, new_mesh: Any, damage: np.ndarray) -> np.ndarray:
        d = np.asarray(damage, float)
        if int(new_mesh.nn) <= len(d):
            return d.copy()
        values = list(d)
        for point in np.asarray(new_mesh.nodes[len(d):], float):
            nearest = int(np.argmin(np.linalg.norm(old_mesh.nodes - point[None, :], axis=1)))
            values.append(float(d[nearest]))
        return np.asarray(values, float)

    def advance(
        self,
        *, mesh: Any, boundary: Any, damage: np.ndarray, displacement: np.ndarray,
        p0: np.ndarray, p1: np.ndarray, direction: np.ndarray, front_id: int,
        **kwargs,
    ):
        p0 = np.asarray(p0, float)
        p1 = np.asarray(p1, float)
        direction = np.asarray(direction, float)
        norm = float(np.linalg.norm(direction))
        requested = float(np.linalg.norm(p1 - p0))
        event_index = int(len(_AUDIT["events"]))
        base_record: dict[str, Any] = {
            "physical_event_index": event_index,
            "front_id": int(front_id),
            "p0_m": p0.tolist(),
            "p1_m": p1.tolist(),
            "requested_length_m": requested,
            "exact_endpoint_preserved": True,
            "exact_direction_preserved": True,
            "atomic_transaction": True,
            "equal_length_physical_partition_retry": False,
            "numerical_subsegments_are_topology_only": True,
        }
        if norm <= 0.0 or requested <= 0.0:
            record = {**base_record, "success": False, "failure": "zero_direction_or_length"}
            _AUDIT["events"].append(record)
            return _cb.CrackAdvanceResult(mesh, boundary, damage, displacement, 0.0, False,
                                          reason="v10051839_zero_direction_or_length")
        direction = direction / norm
        alignment = float(np.dot((p1 - p0) / requested, direction))
        if alignment < 1.0 - 1.0e-8:
            record = {**base_record, "success": False, "failure": "endpoint_direction_mismatch",
                      "endpoint_direction_cosine": alignment}
            _AUDIT["events"].append(record)
            return _cb.CrackAdvanceResult(mesh, boundary, damage, displacement, 0.0, False,
                                          reason="v10051839_endpoint_direction_mismatch")

        tip_ids = self._tip_geometric_node_ids(mesh, p0, int(front_id))
        if not tip_ids:
            record = {**base_record, "success": False, "failure": "tip_node_not_found"}
            _AUDIT["events"].append(record)
            return _cb.CrackAdvanceResult(mesh, boundary, damage, displacement, 0.0, False,
                                          reason="v10051839_tip_node_not_found")
        p0_mesh = np.asarray(mesh.nodes[int(tip_ids[0])], float)
        point_tol = max(1.0e-12, 1.0e-8 * max(requested, mesh.hbar_tip, 1.0e-12))
        if float(np.linalg.norm(p0_mesh - p0)) > point_tol:
            record = {**base_record, "success": False, "failure": "physical_tip_mesh_tip_mismatch",
                      "tip_mesh_error_m": float(np.linalg.norm(p0_mesh - p0))}
            _AUDIT["events"].append(record)
            return _cb.CrackAdvanceResult(mesh, boundary, damage, displacement, 0.0, False,
                                          reason="v10051839_physical_tip_mesh_tip_mismatch")

        snap = self._transaction_snapshot()
        corridor, accepted_record, attempts = self._build_corridor(
            mesh, damage, displacement, p0_mesh, p1, direction, int(front_id)
        )
        if corridor is None:
            self._transaction_rollback(snap)
            record = {
                **base_record,
                "success": False,
                "failure": str(accepted_record.get("reason")),
                "corridor_attempt_count": int(len(attempts)),
                "corridor_attempts": attempts,
            }
            _AUDIT["events"].append(copy.deepcopy(record))
            return _cb.CrackAdvanceResult(
                mesh, boundary, damage, displacement, 0.0, False,
                reason="v10051839_no_feasible_atomic_path_corridor",
            )

        current_mesh = corridor["mesh"]
        current_boundary = corridor["boundary"]
        current_damage = np.asarray(corridor["damage"], float)
        current_u = np.asarray(corridor["displacement"], float)
        path = corridor["path"]
        subsegments = []
        moved_total = 0.0
        try:
            current_point = np.asarray(path[0][1], float)
            for index, (_, next_point) in enumerate(path[1:]):
                next_point = np.asarray(next_point, float)
                before_mesh = current_mesh
                new_mesh, new_u, ce = self._split_segment_topology(
                    current_mesh,
                    current_u,
                    current_point,
                    next_point,
                    int(front_id),
                )
                current_damage = self._extend_damage(before_mesh, new_mesh, current_damage)
                current_mesh = new_mesh
                current_u = np.asarray(new_u, float)
                current_boundary = make_boundary_data(current_mesh, self.geom)
                length_i = float(np.linalg.norm(next_point - current_point))
                moved_total += length_i
                row = {
                    "front_id": int(front_id),
                    "event_index": int(ce.event_index),
                    "physical_event_index": event_index,
                    "topology_subsegment_index": int(index),
                    "topology_subsegment_count": int(len(path) - 1),
                    "x0": float(current_point[0]),
                    "y0": float(current_point[1]),
                    "x1": float(next_point[0]),
                    "y1": float(next_point[1]),
                    "length_m": length_i,
                    "angle_error_deg": 0.0,
                    "damage": float(ce.damage),
                    "reason": "v10051839_atomic_path_corridor",
                    "numerical_subsegment_is_not_physical_event": True,
                }
                self.advance_log.append(row)
                subsegments.append(copy.deepcopy(row))
                current_point = next_point
        except Exception as exc:
            self._transaction_rollback(snap)
            record = {
                **base_record,
                "success": False,
                "failure": f"cohesive_path_commit_error:{type(exc).__name__}:{exc}",
                "accepted_corridor": accepted_record,
                "corridor_attempt_count": int(len(attempts)),
                "topology_subsegments_before_failure": subsegments,
            }
            _AUDIT["events"].append(copy.deepcopy(record))
            return _cb.CrackAdvanceResult(
                mesh, boundary, damage, displacement, 0.0, False,
                reason="v10051839_atomic_path_commit_error",
            )

        endpoint = np.asarray(path[-1][1], float)
        endpoint_error = float(np.linalg.norm(endpoint - p1))
        length_error = abs(float(moved_total) - requested)
        length_tol = max(1.0e-12, 1.0e-8 * requested)
        if endpoint_error > length_tol or length_error > length_tol:
            self._transaction_rollback(snap)
            record = {
                **base_record,
                "success": False,
                "failure": "atomic_path_endpoint_or_length_mismatch",
                "endpoint_error_m": endpoint_error,
                "length_error_m": length_error,
                "accepted_corridor": accepted_record,
                "topology_subsegments": subsegments,
            }
            _AUDIT["events"].append(copy.deepcopy(record))
            return _cb.CrackAdvanceResult(
                mesh, boundary, damage, displacement, 0.0, False,
                reason="v10051839_atomic_path_mismatch",
            )

        record = {
            **base_record,
            "success": True,
            "failure": None,
            "moved_length_m": float(moved_total),
            "length_error_m": float(length_error),
            "endpoint_error_m": float(endpoint_error),
            "accepted_corridor": accepted_record,
            "corridor_attempt_count": int(len(attempts)),
            "rejected_corridor_attempts": attempts[:-1],
            "topology_subsegment_count": int(len(subsegments)),
            "topology_subsegments": subsegments,
        }
        _AUDIT["events"].append(copy.deepcopy(record))
        return _cb.CrackAdvanceResult(
            current_mesh,
            current_boundary,
            current_damage,
            current_u,
            float(moved_total),
            True,
            angle_error_deg=0.0,
            selected_edge_length=float(moved_total),
            reason="v10051839_atomic_target_aware_path_corridor",
            elem_parent_map=np.asarray(corridor["parent_map"], dtype=int),
        )


def audit_payload() -> dict[str, Any]:
    events = copy.deepcopy(_AUDIT["events"])
    success_count = sum(bool(row.get("success")) for row in events)
    return {
        "schema": SCHEMA,
        "model_id": MODEL_ID,
        "atomic_target_aware_corridor_remesh": True,
        "pslg_constrained_exact_path": True,
        "complete_event_remeshed_before_cohesive_commit": True,
        "exact_physical_endpoint_preserved": True,
        "exact_selected_direction_preserved": True,
        "equal_length_physical_partition_retry": False,
        "numerical_path_subsegments_are_independent_hazard_events": False,
        "triangle_quality_floor_relaxed": False,
        "child_area_ratio_floor_relaxed": False,
        "bulk_constitutive_physics_changed": False,
        "event_count": int(len(events)),
        "success_count": int(success_count),
        "failure_count": int(len(events) - success_count),
        "events": events,
    }


__all__ = [
    "AtomicPathCorridorCZMBackendV10051839",
    "MODEL_ID",
    "SCHEMA",
    "audit_payload",
    "reset_audit",
]
