"""Atomic centroid-Steiner cavity refinement for v10.0.5.18.3.8.

When every conforming one-edge midpoint split of an edge-connected target
triangle violates the retained quality floor, split only that target triangle at
its centroid.  The three final children each inherit one third of the parent
area, no neighbouring triangle is touched, and the operation is accepted only
when all final children satisfy the unchanged production quality and area
floors.  Repeated accepted centroid splits can increase the barycentric
clearance of an exact endpoint that lies close to a parent edge.
"""
from __future__ import annotations

from typing import Any

import numpy as np

from .mesh import make_boundary_data, rebuild_tri_mesh
from . import quality_aware_czm_patch_v10051835 as _base
from . import quality_aware_czm_ray_handoff_v10051835 as _handoff
from . import quality_aware_czm_retry_v10051835 as _retry


MODEL_ID = "atomic_centroid_steiner_target_cavity_v10_0_5_18_3_8"
_ORIGINAL_REFINE_ONCE = _handoff.refine_once


def _centroid_split_candidate(
    self: Any,
    state: _base.State,
    p0: np.ndarray,
    target: np.ndarray,
    front_id: int,
    level: int,
):
    """Split the edge-connected target parent into three centroid children."""
    elem_id, shared = _handoff._adjacent_target_triangle(
        self,
        state.mesh,
        np.asarray(p0, dtype=float),
        np.asarray(target, dtype=float),
        int(front_id),
    )
    if elem_id is None:
        return None, {
            "level": int(level),
            "accepted": False,
            "reason": "no_edge_connected_target_parent_for_centroid_split",
        }
    if int(shared) < 2:
        return None, {
            "level": int(level),
            "accepted": False,
            "reason": "centroid_split_requires_shared_edge_connectivity",
            "target_parent_element": int(elem_id),
            "target_fan_shared_node_count": int(shared),
        }

    qfloor = _base.float_env(
        "ARRHENIUS_MIN_ACCEPTED_TRIANGLE_QUALITY",
        float(self.min_triangle_quality),
    )
    afloor = _base.float_env(
        "ARRHENIUS_MIN_ACCEPTED_CHILD_AREA_RATIO",
        float(self.min_area_ratio),
    )

    mesh0 = state.mesh
    nodes0 = np.asarray(mesh0.nodes, dtype=float)
    elems0 = np.asarray(mesh0.elems, dtype=int)
    conn = np.asarray(elems0[int(elem_id)], dtype=int)
    triangle = np.asarray(nodes0[conn], dtype=float)
    centroid = np.mean(triangle, axis=0)
    centroid_id = int(len(nodes0))

    nodes1 = np.vstack([nodes0, centroid[None, :]])
    a, b, c = [int(value) for value in conn]
    children = np.asarray(
        [
            [a, b, centroid_id],
            [b, c, centroid_id],
            [c, a, centroid_id],
        ],
        dtype=int,
    )
    old_sign = float(self._signed_twice_area(nodes0, elems0[[int(elem_id)]])[0])
    for child in children:
        sign = float(self._signed_twice_area(nodes1, child.reshape(1, 3))[0])
        if old_sign * sign < 0.0:
            child[0], child[1] = child[1], child[0]

    elems1 = elems0.copy()
    elems1[int(elem_id)] = children[0]
    elems1 = np.vstack([elems1, children[1], children[2]])
    mesh1 = rebuild_tri_mesh(
        nodes1,
        elems1,
        tip_centers=[np.asarray(p0, dtype=float), np.asarray(target, dtype=float)],
    )

    child_ids = np.asarray(
        [int(elem_id), int(mesh0.ne), int(mesh0.ne) + 1],
        dtype=int,
    )
    quality = self._triangle_quality(mesh1.nodes, mesh1.elems[child_ids])
    qmin = float(np.min(quality)) if quality.size else 1.0
    parent_area = max(float(mesh0.area_e[int(elem_id)]), 1.0e-300)
    area_ratios = np.asarray(mesh1.area_e[child_ids], dtype=float) / parent_area
    amin = float(np.min(area_ratios)) if area_ratios.size else 1.0
    issues = _base.floor_issues(qmin, amin, qfloor, afloor)
    if issues:
        return None, {
            "level": int(level),
            "accepted": False,
            "reason": ";".join(issues),
            "refinement_kind": "atomic_target_parent_centroid_steiner",
            "target_parent_element": int(elem_id),
            "target_fan_shared_node_count_before": int(shared),
            "centroid_m": centroid.tolist(),
            "min_triangle_quality": qmin,
            "triangle_quality_floor": qfloor,
            "min_immediate_child_area_ratio": amin,
            "child_area_ratio_floor": afloor,
        }

    displacement0 = np.asarray(state.displacement, dtype=float).reshape(-1, 2)
    displacement1 = np.vstack(
        [displacement0, np.mean(displacement0[conn], axis=0)[None, :]]
    ).reshape(-1)
    damage0 = np.asarray(state.damage, dtype=float)
    damage1 = np.concatenate(
        [damage0, np.asarray([float(np.mean(damage0[conn]))], dtype=float)]
    )

    parent_map = np.concatenate(
        [
            np.arange(int(mesh0.ne), dtype=int),
            np.asarray([int(elem_id), int(elem_id)], dtype=int),
        ]
    )
    next_state = _base.State(
        mesh=mesh1,
        boundary=make_boundary_data(mesh1, self.geom),
        damage=damage1,
        displacement=displacement1,
        parent_to_root=state.parent_to_root[parent_map],
    )

    candidate_elem, candidate_shared = _handoff._adjacent_target_triangle(
        self,
        mesh1,
        np.asarray(p0, dtype=float),
        np.asarray(target, dtype=float),
        int(front_id),
    )
    if candidate_elem is None or int(candidate_shared) < 2:
        return None, {
            "level": int(level),
            "accepted": False,
            "reason": "centroid_split_lost_shared_edge_tip_fan_connectivity",
            "refinement_kind": "atomic_target_parent_centroid_steiner",
            "target_parent_element": int(elem_id),
            "target_fan_shared_node_count_before": int(shared),
            "target_fan_shared_node_count_after": int(candidate_shared),
            "centroid_m": centroid.tolist(),
            "min_triangle_quality": qmin,
            "min_immediate_child_area_ratio": amin,
        }

    predicted = _base._predicted_exact_target_metrics(
        self,
        mesh1,
        int(candidate_elem),
        np.asarray(target, dtype=float),
    )
    if predicted is None:
        return None, {
            "level": int(level),
            "accepted": False,
            "reason": "cannot_predict_endpoint_after_centroid_split",
            "refinement_kind": "atomic_target_parent_centroid_steiner",
            "target_parent_element": int(elem_id),
            "centroid_m": centroid.tolist(),
        }

    pq = float(predicted["predicted_min_triangle_quality"])
    pa = float(predicted["predicted_min_child_area_ratio"])
    qmargin = pq / max(qfloor, 1.0e-300)
    amargin = pa / max(afloor, 1.0e-300)
    exact_pass = bool(pq >= qfloor and pa >= afloor)
    record = {
        "level": int(level),
        "accepted": True,
        "refinement_kind": "atomic_target_parent_centroid_steiner",
        "target_parent_element_before": int(elem_id),
        "target_parent_element_after": int(candidate_elem),
        "target_fan_shared_node_count_before": int(shared),
        "target_fan_shared_node_count_after": int(candidate_shared),
        "centroid_node_id": int(centroid_id),
        "centroid_m": centroid.tolist(),
        "min_triangle_quality": qmin,
        "triangle_quality_floor": qfloor,
        "min_immediate_child_area_ratio": amin,
        "child_area_ratio_floor": afloor,
        "immediate_parent_area_ratio_enforced": True,
        "n_new_nodes": 1,
        "n_new_elements": 2,
        "predicted_exact_target_pass": exact_pass,
        "predicted_quality_margin": qmargin,
        "predicted_area_margin": amargin,
        "predicted_limiting_margin": min(qmargin, amargin),
        **predicted,
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
        centroid_state, centroid_record = _centroid_split_candidate(
            self,
            state,
            np.asarray(p0, dtype=float),
            np.asarray(target, dtype=float),
            int(front_id),
            int(level),
        )
        if centroid_state is not None:
            return centroid_state, centroid_record
        merged = dict(record)
        merged["centroid_cavity_failure"] = centroid_record
        return None, merged
    return refined, record


def install() -> None:
    _handoff.install()
    _retry._refine_once = refine_once


__all__ = [
    "MODEL_ID",
    "_centroid_split_candidate",
    "install",
    "refine_once",
]
