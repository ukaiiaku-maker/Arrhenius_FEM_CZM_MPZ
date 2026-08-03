"""Endpoint-centered shape-regular patch helpers for v10.0.5.18.3.5."""
from __future__ import annotations

from dataclasses import dataclass
import math
import os
from typing import Any

import numpy as np

from . import mode_i_first_passage_v9_18_5 as _v9185
from .mesh import make_boundary_data


@dataclass
class State:
    mesh: Any
    boundary: Any
    damage: np.ndarray
    displacement: np.ndarray
    parent_to_root: np.ndarray


def float_env(name: str, default: float) -> float:
    try:
        value = float(os.environ.get(name, default))
    except (TypeError, ValueError):
        value = float(default)
    return value if math.isfinite(value) else float(default)


def floor_issues(qmin: float, amin: float, qfloor: float, afloor: float) -> list[str]:
    issues: list[str] = []
    if qmin < qfloor:
        issues.append(f"triangle_quality={qmin:.6e}<{qfloor:.6e}")
    if amin < afloor:
        issues.append(f"child_area_ratio={amin:.6e}<{afloor:.6e}")
    return issues


def compose_parent_map(state: State, result: Any) -> np.ndarray:
    parent = getattr(result, "elem_parent_map", None)
    if parent is None:
        if int(result.mesh.ne) != len(state.parent_to_root):
            raise RuntimeError("CZM changed element count without an element-parent map")
        return state.parent_to_root.copy()
    return state.parent_to_root[np.asarray(parent, dtype=int)]


def _extend_damage(old_mesh: Any, new_mesh: Any, damage: np.ndarray, edge_i: int, edge_j: int) -> np.ndarray:
    old = np.asarray(damage, dtype=float)
    if int(new_mesh.nn) <= len(old):
        return old.copy()
    values = list(old)
    a = np.asarray(old_mesh.nodes[edge_i], dtype=float)
    b = np.asarray(old_mesh.nodes[edge_j], dtype=float)
    edge = b - a
    den = float(edge @ edge)
    for index in range(len(old), int(new_mesh.nn)):
        point = np.asarray(new_mesh.nodes[index], dtype=float)
        xi = 0.5 if den <= 1.0e-300 else float(np.clip(((point - a) @ edge) / den, 0.0, 1.0))
        values.append((1.0 - xi) * old[edge_i] + xi * old[edge_j])
    return np.asarray(values, dtype=float)


def _barycentric(triangle: np.ndarray, point: np.ndarray):
    matrix = np.array([
        [triangle[0, 0], triangle[1, 0], triangle[2, 0]],
        [triangle[0, 1], triangle[1, 1], triangle[2, 1]],
        [1.0, 1.0, 1.0],
    ])
    try:
        return np.linalg.solve(matrix, [point[0], point[1], 1.0])
    except np.linalg.LinAlgError:
        return None


def _target_tip_triangle(self: Any, mesh: Any, p0: np.ndarray, target: np.ndarray, front_id: int) -> int | None:
    tip_ids = self._tip_geometric_node_ids(mesh, p0, front_id)
    if not tip_ids:
        return None
    incident = self._incident_elements(mesh.elems, set(map(int, tip_ids)))
    choices: list[tuple[float, int]] = []
    for elem_id in np.asarray(incident, dtype=int):
        triangle = np.asarray(mesh.nodes[mesh.elems[elem_id]], dtype=float)
        weights = _barycentric(triangle, target)
        if weights is not None and float(np.min(weights)) >= -1.0e-10:
            choices.append((float(np.min(weights)), int(elem_id)))
    if not choices:
        return None
    choices.sort(reverse=True)
    return choices[0][1]


def _tip_node_in_triangle(mesh: Any, elem_id: int, p0: np.ndarray) -> int | None:
    conn = np.asarray(mesh.elems[int(elem_id)], dtype=int)
    distance = np.linalg.norm(np.asarray(mesh.nodes[conn], dtype=float) - np.asarray(p0)[None, :], axis=1)
    local = int(np.argmin(distance))
    scale = max(float(getattr(mesh, "hbar_tip", 0.0)), float(getattr(mesh, "hbar", 0.0)), 1.0e-12)
    tol = max(1.0e-12, 1.0e-6 * scale)
    if float(distance[local]) > tol:
        return None
    return int(conn[local])


def _predicted_exact_target_metrics(self: Any, mesh: Any, elem_id: int, target: np.ndarray) -> dict[str, Any] | None:
    conn = np.asarray(mesh.elems[int(elem_id)], dtype=int)
    triangle = np.asarray(mesh.nodes[conn], dtype=float)
    weights = _barycentric(triangle, np.asarray(target, dtype=float))
    if weights is None or float(np.min(weights)) < -1.0e-10:
        return None
    nodes = np.vstack([np.asarray(mesh.nodes, dtype=float), target[None, :]])
    target_id = int(len(mesh.nodes))
    children = np.array([
        [int(conn[0]), int(conn[1]), target_id],
        [int(conn[1]), int(conn[2]), target_id],
        [int(conn[2]), int(conn[0]), target_id],
    ], dtype=int)
    quality = self._triangle_quality(nodes, children)
    return {
        "target_parent_element": int(elem_id),
        "target_barycentric_weights": np.asarray(weights, dtype=float).tolist(),
        "predicted_min_child_area_ratio": float(np.min(weights)),
        "predicted_min_triangle_quality": float(np.min(quality)) if quality.size else 1.0,
    }


def _midpoint_split_candidate(self: Any, state: State, edge_i: int, edge_j: int, qfloor: float, afloor: float):
    midpoint = 0.5 * (np.asarray(state.mesh.nodes[edge_i], dtype=float) + np.asarray(state.mesh.nodes[edge_j], dtype=float))
    try:
        mesh1, u1, reason, meta, parent = self._insert_point_on_edge(
            state.mesh, state.displacement, midpoint, edge_i, edge_j
        )
    except Exception as exc:
        return None, {"edge_i": int(edge_i), "edge_j": int(edge_j), "reason": f"{type(exc).__name__}:{exc}"}
    if mesh1 is None or parent is None:
        return None, {"edge_i": int(edge_i), "edge_j": int(edge_j), "reason": str(reason)}
    parent = np.asarray(parent, dtype=int)
    affected = _v9185._affected_elements(state.mesh, mesh1)
    quality = self._triangle_quality(mesh1.nodes, mesh1.elems[affected])
    qmin = float(np.min(quality)) if quality.size else 1.0
    ratios = [
        float(mesh1.area_e[e]) / max(float(state.mesh.area_e[parent[e]]), 1.0e-300)
        for e in np.asarray(affected, dtype=int)
        if e < len(parent) and 0 <= parent[e] < int(state.mesh.ne)
    ]
    amin = min(ratios) if ratios else 1.0
    issues = floor_issues(qmin, amin, qfloor, afloor)
    valid = np.all(np.isfinite(mesh1.area_e)) and np.all(np.asarray(mesh1.area_e) > 0.0) and not issues
    if not valid:
        return None, {
            "edge_i": int(edge_i), "edge_j": int(edge_j),
            "reason": ";".join(issues) or "invalid_refined_mesh",
            "min_triangle_quality": qmin,
            "min_immediate_child_area_ratio": amin,
        }
    next_state = State(
        mesh=mesh1,
        boundary=make_boundary_data(mesh1, self.geom),
        damage=_extend_damage(state.mesh, mesh1, state.damage, edge_i, edge_j),
        displacement=np.asarray(u1, dtype=float),
        parent_to_root=state.parent_to_root[parent],
    )
    record = {
        "edge_i": int(edge_i), "edge_j": int(edge_j),
        "midpoint_m": np.asarray(midpoint, dtype=float).tolist(),
        "min_triangle_quality": qmin,
        "triangle_quality_floor": qfloor,
        "min_immediate_child_area_ratio": amin,
        "child_area_ratio_floor": afloor,
        "immediate_parent_area_ratio_enforced": True,
        "n_new_nodes": int(mesh1.nn - state.mesh.nn),
        "n_new_elements": int(mesh1.ne - state.mesh.ne),
        **dict(meta or {}),
    }
    return next_state, record


def refine_once(self: Any, state: State, p0: np.ndarray, target: np.ndarray, front_id: int, level: int):
    elem_id = _target_tip_triangle(self, state.mesh, p0, target, front_id)
    if elem_id is None:
        return None, {"level": level, "accepted": False, "reason": "no_target_tip_triangle"}
    tip_node = _tip_node_in_triangle(state.mesh, elem_id, p0)
    if tip_node is None:
        return None, {
            "level": level, "accepted": False,
            "reason": "target_triangle_has_no_tip_vertex",
            "target_parent_element": int(elem_id),
        }
    qfloor = float_env("ARRHENIUS_MIN_ACCEPTED_TRIANGLE_QUALITY", float(self.min_triangle_quality))
    afloor = float_env("ARRHENIUS_MIN_ACCEPTED_CHILD_AREA_RATIO", float(self.min_area_ratio))
    conn = [int(value) for value in state.mesh.elems[int(elem_id)]]
    radial_edges = [(int(tip_node), int(node)) for node in conn if int(node) != int(tip_node)]
    radial_edges.sort(key=lambda edge: (edge[1], edge[0]))
    candidates = []
    errors = []
    for edge_i, edge_j in radial_edges:
        candidate_state, split = _midpoint_split_candidate(self, state, edge_i, edge_j, qfloor, afloor)
        if candidate_state is None:
            errors.append(split)
            continue
        candidate_elem = _target_tip_triangle(self, candidate_state.mesh, p0, target, front_id)
        if candidate_elem is None:
            errors.append({**split, "reason": "midpoint_split_moved_target_outside_tip_patch"})
            continue
        predicted = _predicted_exact_target_metrics(self, candidate_state.mesh, candidate_elem, np.asarray(target, dtype=float))
        if predicted is None:
            errors.append({**split, "reason": "cannot_predict_exact_target_insertion"})
            continue
        pq = float(predicted["predicted_min_triangle_quality"])
        pa = float(predicted["predicted_min_child_area_ratio"])
        qmargin = pq / max(qfloor, 1.0e-300)
        amargin = pa / max(afloor, 1.0e-300)
        limiting_margin = min(qmargin, amargin)
        exact_pass_predicted = bool(pq >= qfloor and pa >= afloor)
        score = (
            int(exact_pass_predicted), float(limiting_margin),
            float(min(qmargin, 10.0)), float(min(amargin, 10.0)), -int(edge_j),
        )
        candidates.append((score, candidate_state, {
            **split, **predicted,
            "predicted_exact_target_pass": exact_pass_predicted,
            "predicted_quality_margin": qmargin,
            "predicted_area_margin": amargin,
            "predicted_limiting_margin": limiting_margin,
        }))
    if not candidates:
        return None, {
            "level": level, "accepted": False,
            "reason": "no_quality_safe_tip_radial_midpoint_bisection",
            "target_parent_element": int(elem_id), "errors": errors,
        }
    candidates.sort(key=lambda item: item[0], reverse=True)
    _, selected_state, selected = candidates[0]
    record = {
        "level": int(level), "accepted": True,
        "refinement_kind": "endpoint_centered_tip_radial_midpoint",
        "target_parent_element_before": int(elem_id),
        "tip_node": int(tip_node),
        "candidate_count": int(len(candidates)),
        **selected,
    }
    return selected_state, record


__all__ = ["State", "compose_parent_map", "floor_issues", "float_env", "refine_once"]
