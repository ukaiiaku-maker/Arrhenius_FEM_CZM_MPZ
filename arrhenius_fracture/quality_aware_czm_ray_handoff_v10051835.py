"""Ray-handoff extension for the v10.0.5.18.3.5 quality-aware tip patch.

A shape-regular midpoint split can legitimately move a nearby exact endpoint
from the tip-incident child into the adjacent child.  That is not a failed
refinement: it is a valid handoff to the existing exact-ray edge marcher.
This module changes only that narrowly diagnosed case.
"""
from __future__ import annotations

from typing import Any

import numpy as np

from . import quality_aware_czm_patch_v10051835 as _base
from . import quality_aware_czm_retry_v10051835 as _retry


MODEL_ID = "quality_aware_tip_fan_ray_handoff_v10_0_5_18_3_5"
_ORIGINAL_REFINE_ONCE = _base.refine_once


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
            # The base policy handles all candidates that remain in the
            # tip-incident triangle and can predict the exact insertion.
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

    if str(record.get("reason")) != (
        "no_quality_safe_tip_radial_midpoint_bisection"
    ):
        return refined, record

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


def install() -> None:
    """Install the narrow ray-handoff policy into the retry transaction."""
    _retry._refine_once = refine_once


__all__ = ["MODEL_ID", "install", "refine_once"]
