"""Balanced intact-support prerefinement for exact adaptive-CZM endpoints.

The v10.0.5.21 collinear support construction can be infeasible when an exact
crack endpoint lies close to a side of the active tip triangle and the ray
leaves that triangle before a quality-safe forward support point exists.

This additive fallback inserts one intact off-ray support point inside the same
tip triangle, then inserts the exact requested crack endpoint into a resulting
tip-adjacent child.  Both topology subtransactions independently satisfy the
unchanged immediate-parent child-area and triangle-quality floors.  The crack
endpoint, direction, stochastic event length, constitutive state, and hazard
law remain unchanged.
"""
from __future__ import annotations

from contextlib import contextmanager
import math
from typing import Any, Iterator

import numpy as np

from . import crack_backend as _cb
from . import adaptive_czm_quality_subdivision_v100518 as _v100518
from . import adaptive_czm_forward_support_prerefine_v100521 as _v100521


MODEL_ID = "adaptive_CZM_balanced_support_prerefine_v10_0_5_22"
_RAW_INSERT = _cb.AdaptiveCZMBackend._insert_target_in_incident_triangle


def _local_quality(backend: Any, triangle: np.ndarray) -> float:
    nodes = np.asarray(triangle, float).reshape(3, 2)
    elems = np.array([[0, 1, 2]], dtype=int)
    return float(backend._triangle_quality(nodes, elems)[0])


def _candidate_supports(
    backend: Any,
    triangle: np.ndarray,
    target: np.ndarray,
    tip_local_index: int,
    area_floor: float,
    quality_floor: float,
):
    """Return ranked barycentric support candidates with a feasible second split."""
    guard = min(max(float(area_floor) + 0.015, 1.15 * float(area_floor)), 0.18)
    second_guard = min(max(float(area_floor) + 0.005, 1.05 * float(area_floor)), 0.15)
    values = np.linspace(guard, 1.0 - 2.0 * guard, 25)
    candidates = []

    for s0 in values:
        for s1 in values:
            s2 = 1.0 - float(s0) - float(s1)
            support_weights = np.array([s0, s1, s2], dtype=float)
            if float(np.min(support_weights)) < guard - 1.0e-12:
                continue
            support = support_weights @ triangle
            first_children = [
                np.array([triangle[0], triangle[1], support]),
                np.array([triangle[1], triangle[2], support]),
                np.array([triangle[2], triangle[0], support]),
            ]
            first_q = min(_local_quality(backend, child) for child in first_children)
            if first_q < quality_floor:
                continue

            for child_index, child in enumerate(first_children):
                child_vertices = ((0, 1), (1, 2), (2, 0))[child_index]
                if int(tip_local_index) not in child_vertices:
                    continue
                weights = _v100521._barycentric(child, target)
                if weights is None:
                    continue
                minimum = float(np.min(weights))
                maximum = float(np.max(weights))
                if minimum < second_guard or maximum > 1.0 - second_guard:
                    continue
                second_children = [
                    np.array([child[0], child[1], target]),
                    np.array([child[1], child[2], target]),
                    np.array([child[2], child[0], target]),
                ]
                second_q = min(_local_quality(backend, kid) for kid in second_children)
                if second_q < quality_floor:
                    continue
                score = (
                    min(float(np.min(support_weights)), minimum),
                    min(first_q, second_q),
                    -float(np.linalg.norm(support_weights - 1.0 / 3.0)),
                )
                candidates.append((score, support.copy(), support_weights.copy()))
                break

    candidates.sort(key=lambda item: item[0], reverse=True)
    return candidates


def _balanced_support_insert(
    self: Any,
    original: Any,
    mesh: Any,
    displacement: np.ndarray,
    p0: np.ndarray,
    target: np.ndarray,
    front_id: int,
):
    direct = original(self, mesh, displacement, p0, target, int(front_id))
    dmesh, _, _, _, _, dmap = direct
    if dmesh is None or dmap is None:
        return direct
    direct_q, direct_a = _v100521._immediate_metrics(self, mesh, dmesh, dmap)
    qfloor = float(self.min_triangle_quality)
    afloor = float(self.min_area_ratio)
    if direct_q >= qfloor and direct_a >= afloor:
        return direct

    p0 = np.asarray(p0, float).reshape(2)
    target = np.asarray(target, float).reshape(2)
    selected = _v100521._containing_tip_triangle(
        self, mesh, p0, target, int(front_id)
    )
    if selected is None:
        return direct
    _, parent, conn, _, tip_local_index = selected
    triangle = np.asarray(mesh.nodes[np.asarray(conn, dtype=int)], float)

    candidates = _candidate_supports(
        self,
        triangle,
        target,
        int(tip_local_index),
        afloor,
        qfloor,
    )
    for _, support, support_weights in candidates[:64]:
        support_result = _RAW_INSERT(
            self,
            mesh,
            displacement,
            p0,
            support,
            int(front_id),
        )
        smesh, su, sq, sreason, smeta, smap = support_result
        if smesh is None or smap is None or str(sreason) != "ok":
            continue
        support_q, support_a = _v100521._immediate_metrics(
            self, mesh, smesh, smap
        )
        if support_q < qfloor or support_a < afloor:
            continue

        endpoint_result = _RAW_INSERT(
            self,
            smesh,
            su,
            p0,
            target,
            int(front_id),
        )
        emesh, eu, endpoint, ereason, emeta, emap = endpoint_result
        if emesh is None or emap is None or str(ereason) != "ok":
            continue
        if endpoint is None or not np.allclose(
            np.asarray(endpoint, float), target, rtol=0.0, atol=1.0e-13
        ):
            continue
        endpoint_q, endpoint_a = _v100521._immediate_metrics(
            self, smesh, emesh, emap
        )
        if endpoint_q < qfloor or endpoint_a < afloor:
            continue

        smap_array = np.asarray(smap, dtype=int)
        emap_array = np.asarray(emap, dtype=int)
        if np.any(emap_array < 0) or np.any(emap_array >= len(smap_array)):
            continue
        combined = smap_array[emap_array]
        minimum_q = float(min(support_q, endpoint_q))
        minimum_a = float(min(support_a, endpoint_a))
        meta = dict(smeta or {})
        meta.update(emeta or {})
        meta.update(
            {
                "v100522_balanced_support_prerefine": True,
                "v100522_original_parent_element": int(parent),
                "v100522_support_x_m": float(support[0]),
                "v100522_support_y_m": float(support[1]),
                "v100522_support_barycentric_weights": support_weights.tolist(),
                "v100522_min_immediate_triangle_quality": minimum_q,
                "v100522_min_immediate_child_area_ratio": minimum_a,
                "min_triangle_quality": minimum_q,
                "min_area_ratio": minimum_a,
                "topology_subtransaction_count": 2,
                "quality_prerefine_subtransaction_count": 2,
                "exact_target_preserved": True,
                "support_node_is_intact_off_ray": True,
            }
        )
        return emesh, eu, target.copy(), "ok_balanced_support_prerefine", meta, combined

    return direct


@contextmanager
def installed_balanced_support_prerefine_v100522() -> Iterator[None]:
    original_insert = _cb.AdaptiveCZMBackend._insert_target_in_incident_triangle
    original_assess = _v100518._assess

    def repaired_insert(self, mesh, displacement, p0, target, front_id):
        return _balanced_support_insert(
            self,
            original_insert,
            mesh,
            displacement,
            p0,
            target,
            int(front_id),
        )

    def repaired_assess(self, old_mesh, result, kwargs, log_start):
        issues, row, warning = original_assess(
            self, old_mesh, result, kwargs, log_start
        )
        logs = list(getattr(self, "advance_log", []) or [])[int(log_start) :]
        repair_rows = [
            item
            for item in logs
            if isinstance(item, dict)
            and bool(item.get("v100522_balanced_support_prerefine"))
        ]
        if not repair_rows:
            return issues, row, warning
        repair = repair_rows[-1]
        amin = float(repair["v100522_min_immediate_child_area_ratio"])
        afloor = float(row["child_area_ratio_floor"])
        revised = [
            issue for issue in issues if not str(issue).startswith("child_area_ratio=")
        ]
        if amin < afloor:
            revised.append(f"child_area_ratio={amin:.6e}<{afloor:.6e}")
        row.update(
            {
                "min_child_area_ratio": amin,
                "child_area_ratio_reference": (
                    "minimum_immediate_balanced_support_prerefine_subtransaction"
                ),
                "multi_generation_atomic_event": True,
                "topology_subtransaction_count": int(
                    repair.get("quality_prerefine_subtransaction_count", 2)
                ),
                "v100522_balanced_support_prerefine": True,
                "issues": revised,
                "accepted": not bool(revised),
            }
        )
        return revised, row, warning

    _cb.AdaptiveCZMBackend._insert_target_in_incident_triangle = repaired_insert
    _v100518._assess = repaired_assess
    try:
        yield
    finally:
        _v100518._assess = original_assess
        _cb.AdaptiveCZMBackend._insert_target_in_incident_triangle = original_insert


__all__ = [
    "MODEL_ID",
    "installed_balanced_support_prerefine_v100522",
]
