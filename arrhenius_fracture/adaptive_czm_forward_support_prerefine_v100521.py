"""Quality-safe forward support pre-refinement for exact CZM endpoints.

A target lying slightly too close to one side of a tip-adjacent triangle can
produce a well-shaped but undersized child in the direct one-to-three target
insertion.  Uniform crack-event subdivision cannot necessarily change that
final endpoint geometry.

This repair inserts one intact support node farther along the same exact ray,
then inserts the requested crack endpoint on the new tip-to-support edge.  The
support split and endpoint edge split are checked independently against the
unchanged immediate-parent child-area and triangle-quality floors.  The crack
endpoint, crack direction, event length, constitutive state, and hazard law are
unchanged.
"""
from __future__ import annotations

from contextlib import contextmanager
import math
from typing import Any, Iterator

import numpy as np

from . import crack_backend as _cb
from . import adaptive_czm_quality_subdivision_v100518 as _v100518


MODEL_ID = "adaptive_CZM_forward_support_prerefine_v10_0_5_21"


def _barycentric(triangle: np.ndarray, point: np.ndarray) -> np.ndarray | None:
    tri = np.asarray(triangle, float).reshape(3, 2)
    p = np.asarray(point, float).reshape(2)
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


def _changed_elements(old_mesh: Any, new_mesh: Any) -> np.ndarray:
    changed: list[int] = []
    common = min(int(old_mesh.ne), int(new_mesh.ne))
    for e in range(common):
        if not np.array_equal(new_mesh.elems[e], old_mesh.elems[e]):
            changed.append(e)
            continue
        if not math.isclose(
            float(new_mesh.area_e[e]),
            float(old_mesh.area_e[e]),
            rel_tol=1.0e-12,
            abs_tol=1.0e-30,
        ):
            changed.append(e)
    changed.extend(range(common, int(new_mesh.ne)))
    return np.asarray(sorted(set(changed)), dtype=int)


def _immediate_metrics(
    backend: Any,
    old_mesh: Any,
    new_mesh: Any,
    parent_map: Any,
) -> tuple[float, float]:
    changed = _changed_elements(old_mesh, new_mesh)
    if changed.size == 0:
        return 1.0, 1.0
    quality = backend._triangle_quality(new_mesh.nodes, new_mesh.elems[changed])
    qmin = float(np.min(quality)) if quality.size else 1.0
    pm = np.asarray(parent_map, dtype=int)
    if len(pm) != int(new_mesh.ne) or np.any(pm < 0) or np.any(pm >= int(old_mesh.ne)):
        return qmin, 0.0
    ratios = [
        float(new_mesh.area_e[e])
        / max(float(old_mesh.area_e[int(pm[e])]), 1.0e-300)
        for e in changed
    ]
    return qmin, float(min(ratios)) if ratios else 1.0


def _containing_tip_triangle(
    backend: Any,
    mesh: Any,
    p0: np.ndarray,
    target: np.ndarray,
    front_id: int,
):
    tip_ids = [int(i) for i in backend._tip_geometric_node_ids(mesh, p0, front_id)]
    if not tip_ids:
        return None
    incident = backend._incident_elements(mesh.elems, set(tip_ids))
    choices = []
    for e in incident:
        conn = np.asarray(mesh.elems[int(e)], dtype=int)
        w = _barycentric(mesh.nodes[conn], target)
        if w is None or np.min(w) < -1.0e-10 or np.max(w) > 1.0 + 1.0e-10:
            continue
        local_tip = [k for k, nid in enumerate(conn) if int(nid) in tip_ids]
        if not local_tip:
            continue
        ktip = min(
            local_tip,
            key=lambda k: float(np.linalg.norm(mesh.nodes[int(conn[k])] - p0)),
        )
        choices.append((float(np.min(w)), int(e), conn, w, int(ktip)))
    if not choices:
        return None
    choices.sort(key=lambda item: item[0], reverse=True)
    return choices[0]


def _edge_tip_id(mesh: Any, p0: np.ndarray, support_id: int) -> int | None:
    tol = max(1.0e-12, 1.0e-6 * max(float(mesh.hbar_tip), float(mesh.hbar), 1.0e-12))
    tip_ids = np.where(np.linalg.norm(mesh.nodes - p0[None, :], axis=1) <= tol)[0]
    for tip_id in tip_ids.astype(int):
        parents = np.where(
            np.any(mesh.elems == int(tip_id), axis=1)
            & np.any(mesh.elems == int(support_id), axis=1)
        )[0]
        if len(parents):
            return int(tip_id)
    return None


def _forward_support_insert(
    self: Any,
    original: Any,
    mesh: Any,
    displacement: np.ndarray,
    p0: np.ndarray,
    target: np.ndarray,
    front_id: int,
):
    direct = original(self, mesh, displacement, p0, target, front_id)
    dmesh, _, _, _, _, dmap = direct
    if dmesh is None or dmap is None:
        return direct
    direct_q, direct_a = _immediate_metrics(self, mesh, dmesh, dmap)
    qfloor = float(self.min_triangle_quality)
    afloor = float(self.min_area_ratio)
    if direct_q >= qfloor and direct_a >= afloor:
        return direct

    p0 = np.asarray(p0, float).reshape(2)
    target = np.asarray(target, float).reshape(2)
    ray = target - p0
    length = float(np.linalg.norm(ray))
    if not math.isfinite(length) or length <= 1.0e-14:
        return direct

    selected = _containing_tip_triangle(self, mesh, p0, target, int(front_id))
    if selected is None:
        return direct
    _, _, conn, weights, ktip = selected
    wtip = float(weights[ktip])
    other = [float(weights[k]) for k in range(3) if k != ktip]
    if min(other) <= 0.0 or wtip >= 1.0 - 1.0e-12:
        return direct

    guard = min(max(1.10 * afloor, afloor + 0.005), 0.20)
    alpha_min = max(
        1.0 / max(1.0 - guard, 1.0e-12),
        *(guard / max(value, 1.0e-300) for value in other),
        1.0 + 1.0e-8,
    )
    alpha_max = (1.0 - guard) / max(1.0 - wtip, 1.0e-300)
    if not math.isfinite(alpha_min) or not math.isfinite(alpha_max):
        return direct
    if alpha_min >= alpha_max * (1.0 - 1.0e-10):
        return direct

    direction = ray / length
    lo = alpha_min * (1.0 + 1.0e-6)
    hi = alpha_max * (1.0 - 1.0e-6)
    candidates = np.unique(
        np.concatenate(
            [
                np.linspace(lo, hi, 16),
                np.asarray([0.5 * (lo + hi), math.sqrt(lo * hi)]),
            ]
        )
    )

    for alpha in candidates:
        alpha = float(alpha)
        if alpha <= 1.0 or alpha >= alpha_max:
            continue
        support = p0 + alpha * ray
        support_result = original(
            self, mesh, displacement, p0, support, int(front_id)
        )
        smesh, su, sq, sreason, smeta, smap = support_result
        if smesh is None or smap is None or str(sreason) != "ok":
            continue
        support_q, support_a = _immediate_metrics(self, mesh, smesh, smap)
        if support_q < qfloor or support_a < afloor:
            continue

        new_nodes = np.arange(int(mesh.nn), int(smesh.nn), dtype=int)
        if new_nodes.size == 0:
            continue
        support_id = int(
            new_nodes[
                np.argmin(np.linalg.norm(smesh.nodes[new_nodes] - support[None, :], axis=1))
            ]
        )
        if float(np.linalg.norm(smesh.nodes[support_id] - support)) > 1.0e-9:
            continue
        tip_id = _edge_tip_id(smesh, p0, support_id)
        if tip_id is None:
            continue

        emesh, eu, ereason, emeta, emap = self._insert_point_on_edge(
            smesh,
            su,
            target,
            int(tip_id),
            int(support_id),
        )
        if emesh is None or emap is None or str(ereason) != "ok":
            continue
        edge_q, edge_a = _immediate_metrics(self, smesh, emesh, emap)
        if edge_q < qfloor or edge_a < afloor:
            continue

        smap_array = np.asarray(smap, dtype=int)
        emap_array = np.asarray(emap, dtype=int)
        if np.any(emap_array < 0) or np.any(emap_array >= len(smap_array)):
            continue
        combined = smap_array[emap_array]
        minimum_q = float(min(support_q, edge_q))
        minimum_a = float(min(support_a, edge_a))
        meta = dict(smeta or {})
        meta.update(emeta or {})
        meta.update(
            {
                "v100521_forward_support_prerefine": True,
                "v100521_support_x_m": float(support[0]),
                "v100521_support_y_m": float(support[1]),
                "v100521_support_distance_m": float(alpha * length),
                "v100521_target_fraction_of_support_edge": float(1.0 / alpha),
                "v100521_min_immediate_triangle_quality": minimum_q,
                "v100521_min_immediate_child_area_ratio": minimum_a,
                "min_triangle_quality": minimum_q,
                "min_area_ratio": minimum_a,
                "topology_subtransaction_count": 2,
                "quality_prerefine_subtransaction_count": 2,
                "exact_target_preserved": True,
                "support_node_is_intact_ahead_of_crack": True,
            }
        )
        return emesh, eu, target.copy(), "ok_forward_support_prerefine", meta, combined

    return direct


@contextmanager
def installed_forward_support_prerefine_v100521() -> Iterator[None]:
    original_insert = _cb.AdaptiveCZMBackend._insert_target_in_incident_triangle
    original_assess = _v100518._assess

    def repaired_insert(self, mesh, displacement, p0, target, front_id):
        return _forward_support_insert(
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
            and bool(item.get("v100521_forward_support_prerefine"))
        ]
        if not repair_rows:
            return issues, row, warning
        repair = repair_rows[-1]
        amin = float(repair["v100521_min_immediate_child_area_ratio"])
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
                    "minimum_immediate_forward_support_prerefine_subtransaction"
                ),
                "multi_generation_atomic_event": True,
                "topology_subtransaction_count": int(
                    repair.get("quality_prerefine_subtransaction_count", 2)
                ),
                "v100521_forward_support_prerefine": True,
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
    "installed_forward_support_prerefine_v100521",
]
