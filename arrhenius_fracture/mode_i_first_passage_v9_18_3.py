"""v9.18.3 exact-ray recovery for target-on-edge geometry stalls.

The constitutive/event stack is unchanged from v9.18.2/v9.18.1.  This wrapper
patches only the adaptive CZM geometry backend while the solver runs.

Failure reproduced in the 700 K ceramic 100 um gate:

    local_hrefine_error:degenerate elements after topology update: [2126]

The exact requested target lay numerically on an existing triangle edge.  The
legacy incident-triangle fallback accepted a barycentric coordinate near zero
as an interior point, split the parent into three children, and necessarily
created a zero-area child.  v9.18.3 classifies target locations before topology
construction:

* vertex target -> reuse the existing geometric vertex;
* edge target   -> split the edge exactly and inherit parent histories;
* interior      -> perform the original three-child triangle split after an
  explicit pre-rebuild area/quality check.

Repeated identical geometry vetoes also fail fast instead of consuming tens of
thousands of solver steps.  No fracture event is consumed on a failed geometry
transaction.
"""
from __future__ import annotations

import json
import math
import os
from pathlib import Path
import sys
from typing import Any

import numpy as np

from . import crack_backend as _cb
from . import mode_i_first_passage_v9_18_1 as _v9181
from .mesh import rebuild_tri_mesh


def _edge_aware_insert_target_in_incident_triangle(
    self,
    mesh,
    displacement: np.ndarray,
    p0: np.ndarray,
    target: np.ndarray,
    front_id: int,
):
    """Insert an exact-ray target without treating edge points as interior.

    Return signature matches ``AdaptiveCZMBackend._insert_target_in_incident_triangle``.
    """
    p0 = np.asarray(p0, float)
    target = np.asarray(target, float)
    tip_ids = self._tip_geometric_node_ids(mesh, p0, int(front_id))
    if not tip_ids:
        return None, None, None, "tip_node_not_found", {}, None
    incident = self._incident_elements(mesh.elems, set(int(i) for i in tip_ids))
    if len(incident) == 0:
        return None, None, None, "no_tip_incident_elements", {}, None

    choices: list[tuple[float, int, np.ndarray, np.ndarray]] = []
    rhs = np.array([target[0], target[1], 1.0], dtype=float)
    for e in incident:
        conn = np.asarray(mesh.elems[int(e)], dtype=int)
        tri = np.asarray(mesh.nodes[conn], float)
        M = np.array([
            [tri[0, 0], tri[1, 0], tri[2, 0]],
            [tri[0, 1], tri[1, 1], tri[2, 1]],
            [1.0, 1.0, 1.0],
        ], dtype=float)
        try:
            w = np.linalg.solve(M, rhs)
        except np.linalg.LinAlgError:
            continue
        if np.min(w) >= -1.0e-10 and np.max(w) <= 1.0 + 1.0e-10:
            choices.append((float(np.min(w)), int(e), conn.copy(), w.copy()))
    if not choices:
        return None, None, None, "target_not_in_tip_incident_triangle", {}, None

    # Prefer the most interior representation if the target lies on a shared edge.
    choices.sort(key=lambda z: z[0], reverse=True)
    min_w, parent, conn, w = choices[0]
    edge_tol = float(os.environ.get("ARRHENIUS_GEOMETRY_BARYCENTRIC_EDGE_TOL", "1e-8"))
    edge_tol = max(edge_tol, 1.0e-12)
    small = np.where(w <= edge_tol)[0]
    identity_map = np.arange(mesh.ne, dtype=int)

    # Exact/near vertex: use the existing node.  This avoids an endpoint edge
    # split, which would produce a zero-length child triangle.
    if len(small) >= 2:
        local = int(np.argmax(w))
        vid = int(conn[local])
        q = np.asarray(mesh.nodes[vid], float).copy()
        meta = {
            "refined_parent_element": int(parent),
            "n_new_bulk_elements": 0,
            "min_triangle_quality": float("nan"),
            "min_area_ratio": 1.0,
            "node_move_m": 0.0,
            "n_moved_geometric_copies": 0,
            "n_incident_elements": int(len(incident)),
            "target_location_case": "existing_vertex",
            "target_vertex_id": int(vid),
            "target_barycentric_min": float(min_w),
            "target_snap_distance_m": float(np.linalg.norm(q - target)),
        }
        return mesh, np.asarray(displacement, float).copy(), q, "ok", meta, identity_map

    # Exact/near edge: split the opposite edge at the requested target.  This is
    # the required recovery for the element-2126 failure.
    if len(small) == 1:
        opposite = int(small[0])
        edge_local = [k for k in range(3) if k != opposite]
        edge_i = int(conn[edge_local[0]])
        edge_j = int(conn[edge_local[1]])
        a = np.asarray(mesh.nodes[edge_i], float)
        b = np.asarray(mesh.nodes[edge_j], float)
        edge = b - a
        den = float(edge @ edge)
        xi = 0.5 if den <= 1.0e-300 else float(np.clip(((target - a) @ edge) / den, 0.0, 1.0))
        vertex_tol = max(edge_tol, 1.0e-8)
        if xi <= vertex_tol or xi >= 1.0 - vertex_tol:
            vid = edge_i if xi <= vertex_tol else edge_j
            q = np.asarray(mesh.nodes[int(vid)], float).copy()
            meta = {
                "refined_parent_element": int(parent),
                "n_new_bulk_elements": 0,
                "min_triangle_quality": float("nan"),
                "min_area_ratio": 1.0,
                "node_move_m": 0.0,
                "n_moved_geometric_copies": 0,
                "n_incident_elements": int(len(incident)),
                "target_location_case": "edge_endpoint_existing_vertex",
                "target_vertex_id": int(vid),
                "target_edge_xi": float(xi),
                "target_barycentric_min": float(min_w),
                "target_snap_distance_m": float(np.linalg.norm(q - target)),
            }
            return mesh, np.asarray(displacement, float).copy(), q, "ok", meta, identity_map

        try:
            refined, u1, reason, meta, parent_map = self._insert_point_on_edge(
                mesh, displacement, target, edge_i, edge_j
            )
        except Exception as exc:
            return None, None, None, f"target_edge_split_error:{exc}", {}, None
        if refined is None:
            return None, None, None, f"target_edge_split:{reason}", meta, parent_map
        meta = dict(meta)
        meta.update({
            "refined_parent_element": int(parent),
            "target_location_case": "existing_edge_exact_split",
            "target_edge_i": int(edge_i),
            "target_edge_j": int(edge_j),
            "target_edge_xi": float(xi),
            "target_barycentric_min": float(min_w),
        })
        return refined, u1, target.copy(), "ok", meta, parent_map

    # Strictly interior target: retain the original three-child refinement, but
    # validate child areas and quality before calling rebuild_tri_mesh.
    nodes0 = np.asarray(mesh.nodes, float)
    elems0 = np.asarray(mesh.elems, int)
    u0 = np.asarray(displacement, float).reshape(-1, 2)
    new_id = int(len(nodes0))
    nodes1 = np.vstack([nodes0, target[None, :]])
    u_target = w @ u0[conn]
    u1 = np.vstack([u0, u_target[None, :]])

    a_id, b_id, c_id = [int(x) for x in conn]
    children = np.array([
        [a_id, b_id, new_id],
        [b_id, c_id, new_id],
        [c_id, a_id, new_id],
    ], dtype=int)
    old_sign = float(self._signed_twice_area(nodes0, elems0[[parent]])[0])
    for k in range(3):
        sgn = float(self._signed_twice_area(nodes1, children[[k]])[0])
        if old_sign * sgn < 0.0:
            children[k, 0], children[k, 1] = children[k, 1], children[k, 0]

    child_area = np.abs(self._signed_twice_area(nodes1, children))
    old_area = max(abs(old_sign), 1.0e-300)
    min_ratio = float(np.min(child_area) / old_area)
    child_quality = self._triangle_quality(nodes1, children)
    min_quality = float(np.min(child_quality))
    if (not np.all(np.isfinite(child_area)) or not np.all(np.isfinite(child_quality))
            or min_ratio <= 1.0e-10 or min_quality <= 1.0e-10):
        return None, None, None, (
            f"interior_split_precheck_failed:min_area_ratio={min_ratio:.6e};"
            f"min_quality={min_quality:.6e}"
        ), {
            "target_location_case": "interior_precheck_failed",
            "target_barycentric_min": float(min_w),
            "min_area_ratio": min_ratio,
            "min_triangle_quality": min_quality,
        }, None

    elems1 = elems0.copy()
    elems1[parent] = children[0]
    elems1 = np.vstack([elems1, children[1], children[2]])
    try:
        refined = rebuild_tri_mesh(nodes1, elems1, tip_centers=[target])
    except Exception as exc:
        return None, None, None, f"interior_rebuild_error:{exc}", {
            "target_location_case": "interior_rebuild_failed",
            "target_barycentric_min": float(min_w),
            "min_area_ratio": min_ratio,
            "min_triangle_quality": min_quality,
        }, None

    parent_map = np.concatenate([
        np.arange(mesh.ne, dtype=int),
        np.array([parent, parent], dtype=int),
    ])
    meta = {
        "refined_parent_element": int(parent),
        "n_new_bulk_elements": 2,
        "min_triangle_quality": min_quality,
        "min_area_ratio": min_ratio,
        "node_move_m": 0.0,
        "n_moved_geometric_copies": 0,
        "n_incident_elements": int(len(incident)),
        "target_location_case": "strict_triangle_interior",
        "target_barycentric_min": float(min_w),
    }
    return refined, u1.reshape(-1), target.copy(), "ok", meta, parent_map


def _advance_with_identical_veto_guard(self, *args, **kwargs):
    """Abort repeated deterministic geometry vetoes instead of burning steps."""
    original = _advance_with_identical_veto_guard._original
    result = original(self, *args, **kwargs)
    if bool(getattr(result, "inserted", False)):
        self._v9183_last_veto_signature = None
        self._v9183_identical_veto_count = 0
        return result

    reason = str(getattr(result, "reason", "unknown"))
    front_id = int(kwargs.get("front_id", -1) or -1)
    p0 = np.asarray(kwargs.get("p0", [math.nan, math.nan]), float)
    p1 = np.asarray(kwargs.get("p1", [math.nan, math.nan]), float)
    scale = max(float(getattr(kwargs.get("mesh", None), "hbar_tip", 1.0e-12)), 1.0e-12)
    signature = (
        front_id,
        tuple(np.round(p0 / scale, 8)),
        tuple(np.round(p1 / scale, 8)),
        reason,
    )
    if signature == getattr(self, "_v9183_last_veto_signature", None):
        count = int(getattr(self, "_v9183_identical_veto_count", 0)) + 1
    else:
        count = 1
    self._v9183_last_veto_signature = signature
    self._v9183_identical_veto_count = count

    limit = max(int(os.environ.get("ARRHENIUS_MAX_IDENTICAL_GEOMETRY_VETOES", "12")), 1)
    if count >= limit:
        raise RuntimeError(
            "v9.18.3 repeated identical geometry veto; physical event remains "
            f"unconsumed: front={front_id} count={count}/{limit} reason={reason} "
            f"p0={p0.tolist()} p1={p1.tolist()}"
        )
    return result


def _option_value(args: list[str], name: str) -> str | None:
    for i, token in enumerate(args):
        if token == name and i + 1 < len(args):
            return args[i + 1]
        if token.startswith(name + "="):
            return token.split("=", 1)[1]
    return None


def main(argv=None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    original_insert = _cb.AdaptiveCZMBackend._insert_target_in_incident_triangle
    original_advance = _cb.AdaptiveCZMBackend.advance
    _advance_with_identical_veto_guard._original = original_advance
    _cb.AdaptiveCZMBackend._insert_target_in_incident_triangle = (
        _edge_aware_insert_target_in_incident_triangle
    )
    _cb.AdaptiveCZMBackend.advance = _advance_with_identical_veto_guard
    try:
        results = _v9181.main(user_args)
    finally:
        _cb.AdaptiveCZMBackend._insert_target_in_incident_triangle = original_insert
        _cb.AdaptiveCZMBackend.advance = original_advance

    out_value = _option_value(user_args, "--out")
    if out_value is not None:
        out = Path(out_value)
        out.mkdir(parents=True, exist_ok=True)
        payload = {
            "schema": "edge_aware_geometry_recovery_v9183_v1",
            "target_classification": [
                "existing_vertex",
                "existing_edge_exact_split",
                "strict_triangle_interior",
            ],
            "edge_target_never_split_as_triangle_interior": True,
            "identical_geometry_veto_fail_fast_enabled": True,
            "max_identical_geometry_vetoes": int(
                os.environ.get("ARRHENIUS_MAX_IDENTICAL_GEOMETRY_VETOES", "12")
            ),
            "constitutive_physics_changed": False,
        }
        (out / "geometry_recovery_v9183.json").write_text(
            json.dumps(payload, indent=2, default=str)
        )
    return results


if __name__ == "__main__":
    main()
