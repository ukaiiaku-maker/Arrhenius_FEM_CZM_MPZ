"""Shape-regular local patch helpers for v10.0.5.18.3.5."""
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


def floor_issues(
    qmin: float,
    amin: float,
    qfloor: float,
    afloor: float,
) -> list[str]:
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
            raise RuntimeError(
                "CZM changed element count without an element-parent map"
            )
        return state.parent_to_root.copy()
    return state.parent_to_root[np.asarray(parent, dtype=int)]


def _extend_damage(
    old_mesh: Any,
    new_mesh: Any,
    damage: np.ndarray,
    edge_i: int,
    edge_j: int,
) -> np.ndarray:
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
        xi = 0.5 if den <= 1.0e-300 else float(
            np.clip(((point - a) @ edge) / den, 0.0, 1.0)
        )
        values.append((1.0 - xi) * old[edge_i] + xi * old[edge_j])
    return np.asarray(values, dtype=float)


def _barycentric(triangle: np.ndarray, point: np.ndarray):
    matrix = np.array(
        [
            [triangle[0, 0], triangle[1, 0], triangle[2, 0]],
            [triangle[0, 1], triangle[1, 1], triangle[2, 1]],
            [1.0, 1.0, 1.0],
        ]
    )
    try:
        return np.linalg.solve(matrix, [point[0], point[1], 1.0])
    except np.linalg.LinAlgError:
        return None


def _target_tip_triangle(
    self: Any,
    mesh: Any,
    p0: np.ndarray,
    target: np.ndarray,
    front_id: int,
) -> int | None:
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


def _candidate_edges(mesh: Any, elem_id: int, p0: np.ndarray):
    conn = list(map(int, mesh.elems[elem_id]))
    scale = max(float(mesh.hbar_tip), float(mesh.hbar), 1.0e-12)
    non_tip = [
        node
        for node in conn
        if np.linalg.norm(mesh.nodes[node] - p0) > max(1.0e-12, 1.0e-7 * scale)
    ]
    edges: list[tuple[int, int]] = []
    if len(non_tip) >= 2:
        edges.append((non_tip[0], non_tip[1]))
    all_edges = [(conn[0], conn[1]), (conn[1], conn[2]), (conn[2], conn[0])]
    all_edges.sort(
        key=lambda edge: np.linalg.norm(
            mesh.nodes[edge[1]] - mesh.nodes[edge[0]]
        ),
        reverse=True,
    )
    for edge in all_edges:
        if tuple(sorted(edge)) not in {tuple(sorted(item)) for item in edges}:
            edges.append(edge)
    return edges


def refine_once(
    self: Any,
    state: State,
    p0: np.ndarray,
    target: np.ndarray,
    front_id: int,
    level: int,
):
    elem_id = _target_tip_triangle(self, state.mesh, p0, target, front_id)
    if elem_id is None:
        return None, {
            "level": level,
            "accepted": False,
            "reason": "no_target_tip_triangle",
        }

    qfloor = float_env(
        "ARRHENIUS_MIN_ACCEPTED_TRIANGLE_QUALITY",
        float(self.min_triangle_quality),
    )
    afloor = float_env(
        "ARRHENIUS_MIN_ACCEPTED_CHILD_AREA_RATIO",
        float(self.min_area_ratio),
    )
    accepted = []
    errors = []
    for edge_i, edge_j in _candidate_edges(state.mesh, elem_id, p0):
        midpoint = 0.5 * (
            state.mesh.nodes[edge_i] + state.mesh.nodes[edge_j]
        )
        try:
            mesh1, u1, reason, meta, parent = self._insert_point_on_edge(
                state.mesh,
                state.displacement,
                midpoint,
                edge_i,
                edge_j,
            )
        except Exception as exc:
            errors.append(f"{edge_i},{edge_j}:{type(exc).__name__}:{exc}")
            continue
        if mesh1 is None or parent is None:
            errors.append(f"{edge_i},{edge_j}:{reason}")
            continue

        parent = np.asarray(parent, dtype=int)
        affected = _v9185._affected_elements(state.mesh, mesh1)
        quality = self._triangle_quality(mesh1.nodes, mesh1.elems[affected])
        qmin = float(np.min(quality)) if quality.size else 1.0
        ratios = [
            float(mesh1.area_e[e])
            / max(float(state.mesh.area_e[parent[e]]), 1.0e-300)
            for e in np.asarray(affected, dtype=int)
            if e < len(parent) and 0 <= parent[e] < int(state.mesh.ne)
        ]
        amin = min(ratios) if ratios else 1.0
        issues = floor_issues(qmin, amin, qfloor, afloor)
        valid = (
            np.all(np.isfinite(mesh1.area_e))
            and np.all(np.asarray(mesh1.area_e) > 0.0)
            and not issues
        )
        if not valid:
            errors.append(
                f"{edge_i},{edge_j}:q={qmin:.6e},area={amin:.6e}:"
                + ";".join(issues)
            )
            continue

        next_state = State(
            mesh=mesh1,
            boundary=make_boundary_data(mesh1, self.geom),
            damage=_extend_damage(
                state.mesh, mesh1, state.damage, edge_i, edge_j
            ),
            displacement=np.asarray(u1, dtype=float),
            parent_to_root=state.parent_to_root[parent],
        )
        record = {
            "level": level,
            "accepted": True,
            "refined_parent_element": elem_id,
            "split_edge_i": edge_i,
            "split_edge_j": edge_j,
            "midpoint_m": np.asarray(midpoint).tolist(),
            "min_triangle_quality": qmin,
            "triangle_quality_floor": qfloor,
            "min_immediate_child_area_ratio": amin,
            "child_area_ratio_floor": afloor,
            "immediate_parent_area_ratio_enforced": True,
            "n_new_nodes": int(mesh1.nn - state.mesh.nn),
            "n_new_elements": int(mesh1.ne - state.mesh.ne),
            **dict(meta or {}),
        }
        accepted.append((qmin, amin, next_state, record))

    if not accepted:
        return None, {
            "level": level,
            "accepted": False,
            "reason": "no_quality_safe_midpoint_bisection",
            "errors": errors,
        }
    accepted.sort(key=lambda item: (item[0], item[1]), reverse=True)
    return accepted[0][2], accepted[0][3]


__all__ = [
    "State",
    "compose_parent_map",
    "floor_issues",
    "float_env",
    "refine_once",
]
