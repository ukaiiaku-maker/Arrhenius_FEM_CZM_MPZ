"""Local ray-handoff and target-cavity refinement for v10.0.5.18.3.5.

A shape-regular midpoint split can legitimately move a nearby exact endpoint
from the tip-incident child into the adjacent child. That is a valid handoff
to the existing exact-ray edge marcher, but the adjacent target triangle can
still need additional shape-regular refinement before exact insertion. This
module changes only those two narrowly diagnosed local cases.
"""
from __future__ import annotations

from typing import Any

import numpy as np

from . import quality_aware_czm_patch_v10051835 as _base
from . import quality_aware_czm_retry_v10051835 as _retry


MODEL_ID = "quality_aware_tip_fan_ray_handoff_v10_0_5_18_3_5"
_ORIGINAL_REFINE_ONCE = _base.refine_once


def _target_triangle_any(mesh: Any, target: np.ndarray) -> int | None:
    choices: list[tuple[float, int]] = []
    for elem_id, conn in enumerate(np.asarray(mesh.elems, dtype=int)):
        triangle = np.asarray(mesh.nodes[conn], dtype=float)
        weights = _base._barycentric(triangle, target)
        if weights is not None and float(np.min(weights)) >= -1.0e-10:
            choices.append((float(np.min(weights)), int(elem_id)))
    if not choices:
        return None
    choices.sort(reverse=True)
    return choices[0][1]


def _adjacent_target_triangle(
    self: Any,
    mesh: Any,
    p0: np.ndarray,
    target: np.ndarray,
    front_id: int,
):
    elem_id = _target_triangle_any(mesh, target)
    if elem_id is None:
        return None, 0

    tip_ids = self._tip_geometric_node_ids(mesh, p0, front_id)
    if not tip_ids:
        return None, 0
    incident = self._incident_elements(
        mesh.elems, set(map(int, tip_ids))
    )
    if len(incident) == 0:
        return None, 0

    fan_nodes = set(
        map(
            int,
            np.asarray(mesh.elems, dtype=int)[
                np.asarray(incident, dtype=int)
            ].ravel(),
        )
    )
    target_nodes = set(
        map(int, np.asarray(mesh.elems[int(elem_id)], dtype=int))
    )
    shared = len(fan_nodes.intersection(target_nodes))
    if shared == 0:
        return None, 0
    return int(elem_id), int(shared)


def _ray_handoff_candidate(
    self: Any,
    state: _base.State,
    p0: np.ndarray,
    target: np.ndarray,
    front_id: int,
    level: int,
):
    elem_id = _base._target_tip_triangle(
        self, state.mesh, p0, target, front_id
    )
    if elem_id is None:
        return None, {
            "level": int(level),
            "accepted": False,
            "reason": "no_target_tip_triangle",
        }

    tip_node = _base._tip_node_in_triangle(
        state.mesh, elem_id, p0
    )
    if tip_node is None:
        return None, {
            "level": int(level),
            "accepted": False,
            "reason": "target_triangle_has_no_tip_vertex",
            "target_parent_element": int(elem_id),
        }

    qfloor = _base.float_env(
        "ARRHENIUS_MIN_ACCEPTED_TRIANGLE_QUALITY",
        float(self.min_triangle_quality),
    )
    afloor = _base.float_env(
        "ARRHENIUS_MIN_ACCEPTED_CHILD_AREA_RATIO",
        float(self.min_area_ratio),
    )

    ray = np.asarray(target, dtype=float) - np.asarray(p0, dtype=float)
    ray_norm = float(np.linalg.norm(ray))
    if ray_norm <= 0.0:
        return None, {
            "level": int(level),
            "accepted": False,
            "reason": "zero_length_ray_handoff",
            "target_parent_element": int(elem_id),
        }
    ray_unit = ray / ray_norm

    conn = [int(value) for value in state.mesh.elems[int(elem_id)]]
    radial_edges = [
        (int(tip_node), int(node))
        for node in conn
        if int(node) != int(tip_node)
    ]
    radial_edges.sort(key=lambda edge: (edge[1], edge[0]))

    candidates = []
    errors = []
    for edge_i, edge_j in radial_edges:
        candidate_state, split = _base._midpoint_split_candidate(
            self,
            state,
            edge_i,
            edge_j,
            qfloor,
            afloor,
        )
        if candidate_state is None:
            errors.append(split)
            continue

        candidate_elem = _base._target_tip_triangle(
            self,
            candidate_state.mesh,
            p0,
            target,
            front_id,
        )
        if candidate_elem is not None:
            continue

        adjacent_elem, shared = _adjacent_target_triangle(
            self,
            candidate_state.mesh,
            p0,
            target,
            front_id,
        )
        if adjacent_elem is None:
            errors.append(
                {
                    **split,
                    "reason": "split_target_not_adjacent_to_tip_fan",
                }
            )
            continue

        midpoint = np.asarray(split["midpoint_m"], dtype=float)
        vector = midpoint - np.asarray(p0, dtype=float)
        vector_norm = float(np.linalg.norm(vector))
        forward = float(vector @ ray_unit)
        if vector_norm <= 0.0 or forward <= 0.0:
            errors.append(
                {
                    **split,
                    "reason": "ray_handoff_midpoint_not_forward",
                    "ray_forward_projection_m": forward,
                }
            )
            continue

        perpendicular = float(
            np.linalg.norm(vector - forward * ray_unit)
        )
        alignment = float(
            np.clip(forward / vector_norm, -1.0, 1.0)
        )
        qmargin = float(split["min_triangle_quality"]) / max(
            qfloor, 1.0e-300
        )
        amargin = float(split["min_immediate_child_area_ratio"]) / max(
            afloor, 1.0e-300
        )
        score = (
            alignment,
            -perpendicular,
            min(qmargin, 10.0),
            min(amargin, 10.0),
            -int(edge_j),
        )
        candidates.append(
            (
                score,
                candidate_state,
                {
                    **split,
                    "level": int(level),
                    "accepted": True,
                    "refinement_kind": (
                        "endpoint_centered_tip_radial_ray_handoff"
                    ),
                    "target_parent_element_before": int(elem_id),
                    "target_parent_element_after": int(adjacent_elem),
                    "target_fan_shared_node_count": int(shared),
                    "tip_node": int(tip_node),
                    "target_remains_in_tip_triangle": False,
                    "exact_ray_march_required_after_refinement": True,
                    "predicted_exact_target_pass": False,
                    "predicted_min_child_area_ratio": None,
                    "predicted_min_triangle_quality": None,
                    "predicted_quality_margin": None,
                    "predicted_area_margin": None,
                    "predicted_limiting_margin": None,
                    "ray_alignment_cosine": alignment,
                    "ray_forward_projection_m": forward,
                    "ray_perpendicular_distance_m": perpendicular,
                },
            )
        )

    if not candidates:
        return None, {
            "level": int(level),
            "accepted": False,
            "reason": "no_quality_safe_tip_radial_ray_handoff",
            "target_parent_element": int(elem_id),
            "errors": errors,
        }

    candidates.sort(key=lambda item: item[0], reverse=True)
    _, selected_state, selected = candidates[0]
    selected["candidate_count"] = int(len(candidates))
    return selected_state, selected


def _target_cavity_candidate(
    self: Any,
    state: _base.State,
    p0: np.ndarray,
    target: np.ndarray,
    front_id: int,
    level: int,
):
    elem_id, shared = _adjacent_target_triangle(
        self, state.mesh, p0, target, front_id
    )
    if elem_id is None:
        return None, {
            "level": int(level),
            "accepted": False,
            "reason": "no_adjacent_ray_target_triangle",
        }

    qfloor = _base.float_env(
        "ARRHENIUS_MIN_ACCEPTED_TRIANGLE_QUALITY",
        float(self.min_triangle_quality),
    )
    afloor = _base.float_env(
        "ARRHENIUS_MIN_ACCEPTED_CHILD_AREA_RATIO",
        float(self.min_area_ratio),
    )

    conn = [int(value) for value in state.mesh.elems[int(elem_id)]]
    edges = [
        tuple(sorted((conn[0], conn[1]))),
        tuple(sorted((conn[1], conn[2]))),
        tuple(sorted((conn[2], conn[0]))),
    ]
    edges = sorted(set(edges))

    candidates = []
    errors = []
    for edge_i, edge_j in edges:
        candidate_state, split = _base._midpoint_split_candidate(
            self,
            state,
            edge_i,
            edge_j,
            qfloor,
            afloor,
        )
        if candidate_state is None:
            errors.append(split)
            continue

        candidate_elem, candidate_shared = _adjacent_target_triangle(
            self,
            candidate_state.mesh,
            p0,
            target,
            front_id,
        )
        if candidate_elem is None:
            errors.append(
                {
                    **split,
                    "reason": "target_cavity_lost_tip_fan_adjacency",
                }
            )
            continue

        predicted = _base._predicted_exact_target_metrics(
            self,
            candidate_state.mesh,
            candidate_elem,
            np.asarray(target, dtype=float),
        )
        if predicted is None:
            errors.append(
                {
                    **split,
                    "reason": "cannot_predict_target_cavity_insertion",
                }
            )
            continue

        pq = float(predicted["predicted_min_triangle_quality"])
        pa = float(predicted["predicted_min_child_area_ratio"])
        qmargin = pq / max(qfloor, 1.0e-300)
        amargin = pa / max(afloor, 1.0e-300)
        limiting = min(qmargin, amargin)
        exact_pass = bool(pq >= qfloor and pa >= afloor)
        # For non-passing candidates, preserve an edge connection (two shared
        # nodes) between the target triangle and the current tip fan before
        # optimizing the one-step quality margin. The archived step-71 case
        # otherwise greedily chose a slightly better margin that reduced the
        # connection to one vertex and produced a deterministic dead end.
        score = (
            int(exact_pass),
            int(candidate_shared),
            float(limiting),
            float(min(qmargin, 10.0)),
            float(min(amargin, 10.0)),
            -int(edge_j),
            -int(edge_i),
        )
        candidates.append(
            (
                score,
                candidate_state,
                {
                    **split,
                    **predicted,
                    "level": int(level),
                    "accepted": True,
                    "refinement_kind": (
                        "endpoint_centered_adjacent_target_cavity_midpoint"
                    ),
                    "target_parent_element_before": int(elem_id),
                    "target_fan_shared_node_count_before": int(shared),
                    "target_fan_shared_node_count_after": int(
                        candidate_shared
                    ),
                    "target_remains_in_tip_triangle": False,
                    "exact_ray_march_required_after_refinement": True,
                    "predicted_exact_target_pass": exact_pass,
                    "predicted_quality_margin": qmargin,
                    "predicted_area_margin": amargin,
                    "predicted_limiting_margin": limiting,
                },
            )
        )

    if not candidates:
        return None, {
            "level": int(level),
            "accepted": False,
            "reason": "no_quality_safe_adjacent_target_cavity_bisection",
            "target_parent_element": int(elem_id),
            "target_fan_shared_node_count": int(shared),
            "errors": errors,
        }

    candidates.sort(key=lambda item: item[0], reverse=True)
    _, selected_state, selected = candidates[0]
    selected["candidate_count"] = int(len(candidates))
    return selected_state, selected


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
        p0,
        target,
        front_id,
        level,
    )
    if refined is not None:
        return refined, record

    reason = str(record.get("reason"))
    if reason == "no_quality_safe_tip_radial_midpoint_bisection":
        handoff_state, handoff = _ray_handoff_candidate(
            self,
            state,
            np.asarray(p0, dtype=float),
            np.asarray(target, dtype=float),
            int(front_id),
            int(level),
        )
        if handoff_state is not None:
            return handoff_state, handoff
        merged = dict(record)
        merged["ray_handoff_failure"] = handoff
        return None, merged

    if reason == "no_target_tip_triangle":
        cavity_state, cavity = _target_cavity_candidate(
            self,
            state,
            np.asarray(p0, dtype=float),
            np.asarray(target, dtype=float),
            int(front_id),
            int(level),
        )
        if cavity_state is not None:
            return cavity_state, cavity
        merged = dict(record)
        merged["target_cavity_failure"] = cavity
        return None, merged

    return refined, record


def install() -> None:
    """Install the local ray-handoff/cavity policy into retry transactions."""
    _retry._refine_once = refine_once


__all__ = ["MODEL_ID", "install", "refine_once"]
