"""Telemetry-only committed-tip and endpoint audit for v10.0.5.18.3.6.

This module does not modify crack geometry, retry ordering, quality floors,
hazard thresholds, event lengths, constitutive state, or transaction rollback.
It observes one existing physical ``advance`` call before and after delegating to
v10.0.5.18.3's retry backend.
"""
from __future__ import annotations

import copy
import math
from typing import Any

import numpy as np

from . import mode_i_first_passage_v9_18_5 as _v9185
from .numerical_resilience_v1005183 import RetryingAdaptiveCZMBackendV1005183


SCHEMA = "v10.0.5.18.3.6_committed_tip_resolution_and_joint_margin_audit"
MODEL_ID = "committed_tip_resolution_audit_v10_0_5_18_3_6"

_CONFIG: dict[str, float] = {
    "tip_h_fine_m": math.nan,
    "tip_ratio": math.nan,
    "triangle_quality_floor": 0.035,
    "child_area_ratio_floor": 0.08,
    "event_minimum_factor": 0.5,
    "event_maximum_factor": 4.0,
}
_AUDIT: dict[str, Any] = {
    "schema": SCHEMA,
    "model_id": MODEL_ID,
    "events": [],
}


def configure(
    *,
    tip_h_fine_m: float,
    tip_ratio: float,
    triangle_quality_floor: float,
    child_area_ratio_floor: float,
    event_minimum_factor: float,
    event_maximum_factor: float,
) -> None:
    _CONFIG.update(
        {
            "tip_h_fine_m": float(tip_h_fine_m),
            "tip_ratio": float(tip_ratio),
            "triangle_quality_floor": float(triangle_quality_floor),
            "child_area_ratio_floor": float(child_area_ratio_floor),
            "event_minimum_factor": float(event_minimum_factor),
            "event_maximum_factor": float(event_maximum_factor),
        }
    )


def reset_audit() -> None:
    _AUDIT["events"] = []


def _barycentric(triangle: np.ndarray, point: np.ndarray) -> np.ndarray | None:
    tri = np.asarray(triangle, dtype=float)
    p = np.asarray(point, dtype=float)
    matrix = np.column_stack((tri[0] - tri[2], tri[1] - tri[2]))
    det = float(np.linalg.det(matrix))
    scale = max(float(np.max(np.linalg.norm(tri - tri.mean(axis=0), axis=1))), 1.0e-30)
    if abs(det) <= 1.0e-14 * scale * scale:
        return None
    first = np.linalg.solve(matrix, p - tri[2])
    return np.asarray([first[0], first[1], 1.0 - first.sum()], dtype=float)


def _target_element(mesh: Any, point: np.ndarray) -> tuple[int | None, np.ndarray | None]:
    choices: list[tuple[float, int, np.ndarray]] = []
    for elem_id, conn in enumerate(np.asarray(mesh.elems, dtype=int)):
        weights = _barycentric(np.asarray(mesh.nodes, float)[conn], point)
        if weights is None:
            continue
        minimum = float(np.min(weights))
        if minimum >= -1.0e-10:
            choices.append((minimum, int(elem_id), weights))
    if not choices:
        return None, None
    choices.sort(key=lambda row: (row[0], -row[1]), reverse=True)
    _, elem_id, weights = choices[0]
    return elem_id, weights


def _point_segment_distance(point: np.ndarray, a: np.ndarray, b: np.ndarray) -> float:
    ab = np.asarray(b, float) - np.asarray(a, float)
    denom = float(ab @ ab)
    if denom <= 1.0e-300:
        return float(np.linalg.norm(np.asarray(point, float) - np.asarray(a, float)))
    t = float(np.clip(((np.asarray(point, float) - a) @ ab) / denom, 0.0, 1.0))
    return float(np.linalg.norm(np.asarray(point, float) - (a + t * ab)))


def _tip_one_ring(self: Any, mesh: Any, point: np.ndarray, front_id: int) -> dict[str, Any]:
    point = np.asarray(point, dtype=float)
    tip_ids = list(map(int, self._tip_geometric_node_ids(mesh, point, int(front_id))))
    if not tip_ids:
        return {"available": False, "reason": "no_geometric_tip_nodes"}
    incident = np.asarray(
        self._incident_elements(np.asarray(mesh.elems, int), set(tip_ids)),
        dtype=int,
    )
    if incident.size == 0:
        return {"available": False, "reason": "no_incident_tip_elements"}

    nodes = np.asarray(mesh.nodes, dtype=float)
    elems = np.asarray(mesh.elems, dtype=int)
    scale = max(float(getattr(mesh, "hbar", 0.0) or 0.0), 1.0e-12)
    quant = max(1.0e-12, 1.0e-7 * scale)
    edges: dict[tuple[tuple[int, int], tuple[int, int]], float] = {}
    for tri in elems[incident]:
        for ia, ib in ((0, 1), (1, 2), (2, 0)):
            a = nodes[int(tri[ia])]
            b = nodes[int(tri[ib])]
            length = float(np.linalg.norm(b - a))
            if length <= quant * 1.0e-3:
                continue
            ka = tuple(np.round(a / quant).astype(np.int64))
            kb = tuple(np.round(b / quant).astype(np.int64))
            edges[tuple(sorted((ka, kb)))] = length
    lengths = np.asarray(list(edges.values()), dtype=float)
    quality = np.asarray(
        self._triangle_quality(nodes, elems[incident]),
        dtype=float,
    )
    if lengths.size == 0:
        return {"available": False, "reason": "no_nonzero_unique_edges"}
    return {
        "available": True,
        "point_m": point.tolist(),
        "geometric_node_ids": tip_ids,
        "geometric_node_count": int(len(tip_ids)),
        "incident_element_ids": incident.tolist(),
        "incident_element_count": int(incident.size),
        "unique_nonzero_edge_count": int(lengths.size),
        "h_mean_m": float(np.mean(lengths)),
        "h_median_m": float(np.median(lengths)),
        "h_p90_m": float(np.quantile(lengths, 0.90)),
        "h_max_m": float(np.max(lengths)),
        "triangle_quality_min": float(np.min(quality)) if quality.size else None,
        "triangle_quality_max": float(np.max(quality)) if quality.size else None,
        "stored_mesh_hbar_tip_m": float(getattr(mesh, "hbar_tip", math.nan)),
        "global_mesh_hbar_m": float(getattr(mesh, "hbar", math.nan)),
    }


def _orientation(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    return float((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]))


def _segments_intersect(a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray, tol: float) -> bool:
    o1 = _orientation(a, b, c)
    o2 = _orientation(a, b, d)
    o3 = _orientation(c, d, a)
    o4 = _orientation(c, d, b)
    if ((o1 > tol and o2 < -tol) or (o1 < -tol and o2 > tol)) and (
        (o3 > tol and o4 < -tol) or (o3 < -tol and o4 > tol)
    ):
        return True
    return min(
        _point_segment_distance(c, a, b),
        _point_segment_distance(d, a, b),
        _point_segment_distance(a, c, d),
        _point_segment_distance(b, c, d),
    ) <= tol


def _segment_intersects_triangle(p0: np.ndarray, p1: np.ndarray, triangle: np.ndarray, tol: float) -> bool:
    for point in (p0, p1):
        weights = _barycentric(triangle, point)
        if weights is not None and float(np.min(weights)) >= -1.0e-10:
            return True
    for ia, ib in ((0, 1), (1, 2), (2, 0)):
        if _segments_intersect(p0, p1, triangle[ia], triangle[ib], tol):
            return True
    return False


def _ray_corridor(self: Any, mesh: Any, p0: np.ndarray, p1: np.ndarray) -> dict[str, Any]:
    nodes = np.asarray(mesh.nodes, dtype=float)
    elems = np.asarray(mesh.elems, dtype=int)
    tol = max(1.0e-14, 1.0e-9 * max(float(getattr(mesh, "hbar", 0.0) or 0.0), 1.0e-12))
    selected = [
        int(elem_id)
        for elem_id, conn in enumerate(elems)
        if _segment_intersects_triangle(p0, p1, nodes[conn], tol)
    ]
    if not selected:
        return {"available": False, "reason": "no_intersected_elements"}
    selected_array = np.asarray(selected, dtype=int)
    edge_ids = np.vstack(
        [
            elems[selected_array][:, [0, 1]],
            elems[selected_array][:, [1, 2]],
            elems[selected_array][:, [2, 0]],
        ]
    )
    edge_ids = np.unique(np.sort(edge_ids, axis=1), axis=0)
    lengths = np.linalg.norm(nodes[edge_ids[:, 1]] - nodes[edge_ids[:, 0]], axis=1)
    quality = np.asarray(self._triangle_quality(nodes, elems[selected_array]), dtype=float)
    return {
        "available": True,
        "element_ids": selected,
        "element_count": int(len(selected)),
        "h_mean_m": float(np.mean(lengths)),
        "h_max_m": float(np.max(lengths)),
        "triangle_quality_min": float(np.min(quality)),
    }


def _target_metrics(self: Any, mesh: Any, target: np.ndarray) -> dict[str, Any]:
    elem_id, weights = _target_element(mesh, target)
    if elem_id is None or weights is None:
        return {"available": False, "reason": "target_not_in_bulk_triangle"}
    conn = np.asarray(mesh.elems[int(elem_id)], dtype=int)
    triangle = np.asarray(mesh.nodes, dtype=float)[conn]
    quality = np.asarray(
        self._triangle_quality(np.asarray(mesh.nodes, float), conn.reshape(1, 3)),
        dtype=float,
    )
    vertex_distances = np.linalg.norm(triangle - np.asarray(target, float)[None, :], axis=1)
    edge_distances = np.asarray(
        [
            _point_segment_distance(target, triangle[ia], triangle[ib])
            for ia, ib in ((0, 1), (1, 2), (2, 0))
        ],
        dtype=float,
    )
    return {
        "available": True,
        "element_id": int(elem_id),
        "node_ids": conn.tolist(),
        "parent_triangle_quality": float(quality[0]),
        "barycentric_coordinates": weights.tolist(),
        "minimum_barycentric_coordinate": float(np.min(weights)),
        "distance_to_nearest_vertex_m": float(np.min(vertex_distances)),
        "distance_to_nearest_edge_m": float(np.min(edge_distances)),
        "parent_area_m2": float(mesh.area_e[int(elem_id)]),
    }


def _candidate_quality(self: Any, old_mesh: Any, result: Any, endpoint: np.ndarray, front_id: int) -> dict[str, Any]:
    if not bool(getattr(result, "inserted", False)):
        return {
            "inserted": False,
            "reason": str(getattr(result, "reason", "unknown")),
        }
    new_mesh = result.mesh
    affected = np.asarray(_v9185._affected_elements(old_mesh, new_mesh), dtype=int)
    quality = np.asarray(
        self._triangle_quality(new_mesh.nodes, new_mesh.elems[affected]),
        dtype=float,
    )
    ratios: list[float] = []
    parent_map = getattr(result, "elem_parent_map", None)
    if parent_map is not None:
        parent_map = np.asarray(parent_map, dtype=int)
        for elem_id in affected:
            if elem_id < len(parent_map) and 0 <= parent_map[elem_id] < int(old_mesh.ne):
                ratios.append(
                    float(new_mesh.area_e[elem_id])
                    / max(float(old_mesh.area_e[parent_map[elem_id]]), 1.0e-300)
                )
    qmin = float(np.min(quality)) if quality.size else 1.0
    amin = float(min(ratios)) if ratios else 1.0
    qfloor = float(_CONFIG["triangle_quality_floor"])
    afloor = float(_CONFIG["child_area_ratio_floor"])
    q_margin = qmin / max(qfloor, 1.0e-300)
    area_margin = amin / max(afloor, 1.0e-300)
    if qmin >= qfloor and amin >= afloor:
        classification = "joint_pass"
    elif qmin < qfloor and amin >= afloor:
        classification = "quality_only_failure"
    elif qmin >= qfloor and amin < afloor:
        classification = "area_only_failure"
    else:
        classification = "joint_failure"
    return {
        "inserted": True,
        "reason": str(getattr(result, "reason", "ok")),
        "moved_m": float(getattr(result, "moved", 0.0)),
        "affected_element_ids": affected.tolist(),
        "affected_element_count": int(affected.size),
        "min_triangle_quality": qmin,
        "min_child_area_ratio": amin,
        "triangle_quality_floor": qfloor,
        "child_area_ratio_floor": afloor,
        "quality_margin": q_margin,
        "area_margin": area_margin,
        "limiting_normalized_margin": min(q_margin, area_margin),
        "joint_margin_classification": classification,
        "candidate_endpoint_one_ring": _tip_one_ring(
            self, new_mesh, endpoint, front_id
        ),
    }


def clipped_exponential_mean(minimum_factor: float, maximum_factor: float) -> float:
    a = max(float(minimum_factor), 0.0)
    b = max(float(maximum_factor), a)
    return max(a + math.exp(-a) - math.exp(-b), 1.0e-300)


def event_length_diagnostics(threshold_action: float, nominal_length_m: float) -> dict[str, Any]:
    minimum = float(_CONFIG["event_minimum_factor"])
    maximum = float(_CONFIG["event_maximum_factor"])
    mean = clipped_exponential_mean(minimum, maximum)
    threshold = float(threshold_action)
    clipped = min(max(threshold, minimum), maximum)
    if threshold <= minimum:
        clip_state = "lower"
    elif threshold >= maximum:
        clip_state = "upper"
    else:
        clip_state = "none"
    return {
        "hazard_threshold_action": threshold,
        "event_length_minimum_factor": minimum,
        "event_length_maximum_factor": maximum,
        "clipped_exponential_mean": mean,
        "unclipped_event_length_factor": threshold / mean,
        "realized_event_length_factor": clipped / mean,
        "event_length_clip_state": clip_state,
        "lower_clip_active": clip_state == "lower",
        "upper_clip_active": clip_state == "upper",
        "nominal_event_length_m": float(nominal_length_m),
        "realized_event_length_m": float(nominal_length_m) * clipped / mean,
    }


def tip_h_fine_contract() -> dict[str, Any]:
    return {
        "tip_h_fine_m": float(_CONFIG["tip_h_fine_m"]),
        "tip_ratio": float(_CONFIG["tip_ratio"]),
        "tip_h_fine_semantics": (
            "first radial-ring increment and target radial/angular point-cloud "
            "spacing at the refinement center"
        ),
        "graded_size_law": "h(r)=min(tip_h_fine+(tip_ratio-1)*r,h_far)",
        "initial_hbar_tip_semantics": (
            "mean edge length of elements in the nearest max(4,2_percent) "
            "centroid patch; not required to equal tip_h_fine"
        ),
        "tip_h_fine_equals_one_ring_mean_edge": False,
    }


class AuditedRetryingAdaptiveCZMBackendV10051836(
    RetryingAdaptiveCZMBackendV1005183
):
    """Existing retry backend with observation-only transaction telemetry."""

    name = "adaptive_czm_v1005183_retry_v10051836_audited"

    def advance(self, **kwargs):
        # Internal ray-march or partition subsegments are implementation details;
        # only the outer physical event receives a committed-tip audit record.
        nested = int(kwargs.get("_subdepth", 0) or 0) > 0 or bool(
            kwargs.get("_partition_retry_active", False)
        )
        if nested:
            return super().advance(**kwargs)

        mesh = kwargs["mesh"]
        p0 = np.asarray(kwargs["p0"], dtype=float)
        p1 = np.asarray(kwargs["p1"], dtype=float)
        front_id = int(kwargs.get("front_id", 0))
        requested = float(np.linalg.norm(p1 - p0))
        committed = _tip_one_ring(self, mesh, p0, front_id)
        record: dict[str, Any] = {
            "physical_event_index": int(len(_AUDIT["events"])),
            "front_id": front_id,
            "committed_tip_point_m": p0.tolist(),
            "requested_endpoint_m": p1.tolist(),
            "requested_event_length_m": requested,
            "requested_direction": np.asarray(kwargs.get("direction", []), float).tolist(),
            "committed_tip_one_ring": committed,
            "committed_tip_h_over_requested_da": (
                None
                if not committed.get("available") or requested <= 0.0
                else float(committed["h_mean_m"]) / requested
            ),
            "target_parent": _target_metrics(self, mesh, p1),
            "ray_corridor": _ray_corridor(self, mesh, p0, p1),
        }
        result = super().advance(**kwargs)
        endpoint = p1
        if bool(getattr(result, "inserted", False)) and front_id in getattr(self, "tip_nodes", {}):
            endpoint = np.asarray(self.tip_nodes[front_id][2], dtype=float)
        record["backend_result"] = _candidate_quality(
            self, mesh, result, endpoint, front_id
        )
        _AUDIT["events"].append(record)
        return result


def audit_payload() -> dict[str, Any]:
    events = copy.deepcopy(list(_AUDIT["events"]))
    classes: dict[str, int] = {}
    for row in events:
        classification = str(
            row.get("backend_result", {}).get(
                "joint_margin_classification", "not_inserted"
            )
        )
        classes[classification] = classes.get(classification, 0) + 1
    return {
        "schema": SCHEMA,
        "model_id": MODEL_ID,
        "telemetry_only": True,
        "geometry_changed": False,
        "hazard_or_event_length_changed": False,
        "quality_floors_changed": False,
        "partition_policy_changed": False,
        "tip_h_fine_contract": tip_h_fine_contract(),
        "configuration": copy.deepcopy(_CONFIG),
        "event_count": int(len(events)),
        "joint_margin_class_counts": classes,
        "events": events,
    }


__all__ = [
    "AuditedRetryingAdaptiveCZMBackendV10051836",
    "MODEL_ID",
    "SCHEMA",
    "audit_payload",
    "configure",
    "event_length_diagnostics",
    "reset_audit",
    "tip_h_fine_contract",
]
