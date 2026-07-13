"""Replaceable crack-geometry backends for the Arrhenius FEM solver.

The physics above this layer is intentionally unchanged: FrontEngine, fatigue
cycle integration, directional hazard clocks, plasticity, process-zone fields,
and branch bookkeeping all remain authoritative.  This module only changes how
a completed crack-opening renewal is represented mechanically.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Optional

import numpy as np

from .cohesive import CohesiveElement, CohesiveNetwork
from .mesh import TriMesh, make_boundary_data, rebuild_tri_mesh


@dataclass
class CrackAdvanceResult:
    mesh: TriMesh
    boundary: object
    damage: np.ndarray
    displacement: np.ndarray
    moved: float
    inserted: bool
    angle_error_deg: float = 0.0
    selected_edge_length: float = 0.0
    reason: str = "ok"
    # Optional map from each NEW bulk element index to its parent OLD element
    # index.  Identity/topology-only backends leave this as None; local
    # h-refinement returns a map so Gauss-point history fields can be inherited
    # conservatively without global remeshing or interpolation.
    elem_parent_map: Optional[np.ndarray] = None


class CrackBackend:
    name = "base"

    def __init__(self) -> None:
        self.cohesive_network: Optional[CohesiveNetwork] = None

    def advance(self, **kwargs) -> CrackAdvanceResult:  # pragma: no cover - interface
        raise NotImplementedError

    def write_diagnostics(self, out_dir: str) -> None:
        return


class SharpWakeBackend(CrackBackend):
    """Compatibility backend: preserve the original stiffness-kill geometry."""

    name = "sharp_wake"

    def advance(
        self,
        *, mesh: TriMesh, boundary, damage: np.ndarray, displacement: np.ndarray,
        p0: np.ndarray, p1: np.ndarray, kill_r: float, **kwargs,
    ) -> CrackAdvanceResult:
        p0 = np.asarray(p0, float); p1 = np.asarray(p1, float)
        seg = p1 - p0
        L = float(np.linalg.norm(seg))
        if L <= 0.0:
            return CrackAdvanceResult(mesh, boundary, damage, displacement, 0.0, False,
                                      reason="zero_length")
        c = mesh.nodes[mesh.elems].mean(axis=1)
        L2 = float(seg @ seg) + 1e-30
        tt = np.clip(((c[:, 0] - p0[0]) * seg[0] + (c[:, 1] - p0[1]) * seg[1]) / L2,
                     0.0, 1.0)
        proj = p0[None, :] + tt[:, None] * seg[None, :]
        dist2 = np.sum((c - proj) ** 2, axis=1)
        erad = np.sqrt(np.maximum(mesh.area_e, 1e-30))
        rad = np.maximum(float(kill_r), 0.7 * erad)
        dnew = damage.copy()
        dnew[mesh.elems[dist2 <= rad ** 2]] = 1.0
        return CrackAdvanceResult(mesh, boundary, dnew, displacement, L, True,
                                  selected_edge_length=L)


class EdgeSplitCZMBackend(CrackBackend):
    """Topology-preserving discrete CZM backend.

    A completed Arrhenius renewal selects a physical direction.  The backend
    chooses the best existing edge in that direction, duplicates the edge nodes
    by classifying the node stars on the two crack sides, and inserts a
    zero-thickness cohesive interface.  Bulk element count and integration-point
    histories are unchanged, which is the key migration property.

    This backend deliberately records the angular steering error.  A future
    LocalPatchCZMBackend can replace only ``_select_edge``/topology construction
    while the hazard/fatigue/plasticity architecture remains unchanged.
    """

    name = "edge_split_czm"

    def __init__(
        self,
        *, geom,
        penalty_normal_Pa_per_m: float = 1.0e18,
        penalty_tangent_Pa_per_m: float = 1.0e18,
        max_angle_error_deg: float = 35.0,
        event_damage: float = 1.0,
    ) -> None:
        super().__init__()
        self.geom = geom
        self.max_angle_error_deg = float(max_angle_error_deg)
        self.event_damage = float(np.clip(event_damage, 0.0, 1.0))
        self.cohesive_network = CohesiveNetwork(
            penalty_normal_Pa_per_m=float(penalty_normal_Pa_per_m),
            penalty_tangent_Pa_per_m=float(penalty_tangent_Pa_per_m),
        )
        # front_id -> (plus_tip_node, minus_tip_node, geometric xy)
        self.tip_nodes: Dict[int, tuple[int, int, np.ndarray]] = {}
        self.event_counter = 0
        self.advance_log = []

    @staticmethod
    def _geometric_groups(nodes: np.ndarray, tol: float = 1e-14):
        scale = max(tol, 1e-15)
        keys = np.round(nodes / scale).astype(np.int64)
        groups: Dict[tuple[int, int], list[int]] = {}
        for i, k in enumerate(keys):
            groups.setdefault((int(k[0]), int(k[1])), []).append(i)
        return groups

    @staticmethod
    def _incident_elements(elems: np.ndarray, node_ids: set[int]) -> np.ndarray:
        mask = np.zeros(len(elems), dtype=bool)
        for nid in node_ids:
            mask |= np.any(elems == int(nid), axis=1)
        return np.where(mask)[0]

    def _tip_geometric_node_ids(self, mesh: TriMesh, p0: np.ndarray, front_id: int) -> list[int]:
        if front_id in self.tip_nodes:
            a, b, _ = self.tip_nodes[front_id]
            return [int(a), int(b)]
        dist = np.linalg.norm(mesh.nodes - p0[None, :], axis=1)
        dmin = float(np.min(dist))
        tol = max(1e-12, 1e-6 * max(mesh.hbar_tip, mesh.hbar, 1e-12))
        return np.where(dist <= dmin + tol)[0].astype(int).tolist()

    def _select_edge(self, mesh: TriMesh, p0: np.ndarray, direction: np.ndarray,
                     front_id: int) -> tuple[Optional[np.ndarray], float, str]:
        direction = np.asarray(direction, float)
        nd = float(np.linalg.norm(direction))
        if nd <= 0.0:
            return None, 180.0, "zero_direction"
        direction /= nd
        tip_ids = self._tip_geometric_node_ids(mesh, p0, front_id)
        tip_set = set(tip_ids)
        incident = self._incident_elements(mesh.elems, tip_set)
        neigh = set()
        for e in incident:
            for nid in mesh.elems[e]:
                if int(nid) not in tip_set:
                    neigh.add(int(nid))
        if not neigh:
            return None, 180.0, "no_tip_neighbors"

        # Deduplicate coincident neighbor nodes by geometry.
        candidates = []
        seen = set()
        for nid in neigh:
            q = mesh.nodes[nid]
            key = tuple(np.round(q / max(1e-14, 1e-8 * max(mesh.hbar_tip, mesh.hbar))).astype(np.int64))
            if key in seen:
                continue
            seen.add(key)
            v = q - p0
            L = float(np.linalg.norm(v))
            if L <= 1e-14:
                continue
            t = v / L
            dot = float(direction @ t)
            if dot <= 1e-8:
                continue
            ang = float(np.degrees(np.arccos(np.clip(dot, -1.0, 1.0))))
            # Favor direction first, then a local edge length close to tip scale.
            candidates.append((ang, abs(L - max(mesh.hbar_tip, 1e-30)), L, q.copy()))
        if not candidates:
            return None, 180.0, "no_forward_edge"
        candidates.sort(key=lambda x: (x[0], x[1]))
        ang, _, _, q = candidates[0]
        if ang > self.max_angle_error_deg:
            return None, ang, "angle_error_exceeds_limit"
        return q, ang, "ok"

    @staticmethod
    def _coincident_ids(nodes: np.ndarray, point: np.ndarray, tol: float) -> list[int]:
        d = np.linalg.norm(nodes - point[None, :], axis=1)
        return np.where(d <= tol)[0].astype(int).tolist()

    @staticmethod
    def _classify_element_sides(mesh: TriMesh, p0: np.ndarray, p1: np.ndarray,
                                element_ids: np.ndarray) -> np.ndarray:
        seg = p1 - p0
        c = mesh.nodes[mesh.elems[element_ids]].mean(axis=1)
        rel = c - p0[None, :]
        cross = seg[0] * rel[:, 1] - seg[1] * rel[:, 0]
        # Exact zeros are assigned by a tiny centroid-normal projection fallback.
        return cross

    def _split_segment_topology(self, mesh: TriMesh, displacement: np.ndarray,
                                p0: np.ndarray, p1: np.ndarray, front_id: int):
        nodes = mesh.nodes.copy()
        elems = mesh.elems.copy()
        u = displacement.reshape(-1, 2).copy()
        tol = max(1e-12, 1e-6 * max(mesh.hbar_tip, mesh.hbar, 1e-12))

        # Endpoint node populations. p0 may already be duplicated; p1 should be
        # intact before the extension but the routine tolerates prior duplicates.
        ids0 = self._coincident_ids(nodes, p0, tol)
        ids1 = self._coincident_ids(nodes, p1, tol)
        if not ids0 or not ids1:
            raise RuntimeError("CZM edge endpoints are not represented by mesh nodes")

        # Create/identify a plus/minus pair at each endpoint.
        def ensure_pair(ids: list[int], point: np.ndarray):
            nonlocal nodes, u
            if len(ids) >= 2:
                return int(ids[0]), int(ids[1])
            base = int(ids[0])
            dup = len(nodes)
            nodes = np.vstack([nodes, point[None, :]])
            u = np.vstack([u, u[base][None, :]])
            return base, dup

        p0_plus, p0_minus = ensure_pair(ids0, p0)
        p1_plus, p1_minus = ensure_pair(ids1, p1)

        # Reassign node stars according to side of the new crack tangent.  This
        # changes connectivity only; element count and material GP histories stay
        # exactly aligned with their original element indices.
        for point, plus_id, minus_id in ((p0, p0_plus, p0_minus),
                                         (p1, p1_plus, p1_minus)):
            coinc = set(self._coincident_ids(nodes, point, tol)) | {plus_id, minus_id}
            incident = self._incident_elements(elems, coinc)
            side = self._classify_element_sides(
                rebuild_tri_mesh(nodes, elems, tip_centers=[p1], validate=False),
                p0, p1, incident)
            for loc, e in enumerate(incident):
                # Replace any coincident endpoint id by the selected side copy.
                target = plus_id if side[loc] >= 0.0 else minus_id
                for a in range(3):
                    if int(elems[e, a]) in coinc:
                        elems[e, a] = target

        new_mesh = rebuild_tri_mesh(nodes, elems, tip_centers=[p1])
        new_u = u.reshape(-1)
        tangent = p1 - p0
        L = float(np.linalg.norm(tangent))
        tangent /= max(L, 1e-300)
        normal = np.array([-tangent[1], tangent[0]], dtype=float)
        elem = CohesiveElement(
            plus_nodes=(p0_plus, p1_plus),
            minus_nodes=(p0_minus, p1_minus),
            normal=normal,
            tangent=tangent,
            length=L,
            damage=self.event_damage,
            front_id=int(front_id),
            event_index=int(self.event_counter),
            barrier_kind="exp_floor",
        )
        self.cohesive_network.add(elem)
        self.tip_nodes[int(front_id)] = (p1_plus, p1_minus, p1.copy())
        self.event_counter += 1
        return new_mesh, new_u, elem

    def advance(
        self,
        *, mesh: TriMesh, boundary, damage: np.ndarray, displacement: np.ndarray,
        p0: np.ndarray, p1: np.ndarray, direction: np.ndarray, front_id: int,
        **kwargs,
    ) -> CrackAdvanceResult:
        p0 = np.asarray(p0, float)
        direction = np.asarray(direction, float)
        tip_ids = self._tip_geometric_node_ids(mesh, p0, int(front_id))
        if not tip_ids:
            return CrackAdvanceResult(mesh, boundary, damage, displacement, 0.0, False,
                                      angle_error_deg=180.0, reason="tip_node_not_found")
        # The front geometry is continuous/physical, whereas an edge-split backend
        # must anchor the topological cut to an actual mesh vertex. Snap only the
        # topology operation; subsequent front state follows the actual inserted edge.
        p0_mesh = mesh.nodes[int(tip_ids[0])].copy()
        q, angle_error, reason = self._select_edge(mesh, p0_mesh, direction, int(front_id))
        if q is None:
            return CrackAdvanceResult(mesh, boundary, damage, displacement, 0.0, False,
                                      angle_error_deg=float(angle_error), reason=reason)
        try:
            new_mesh, new_u, ce = self._split_segment_topology(
                mesh, displacement, p0_mesh, q, int(front_id))
        except Exception as exc:
            return CrackAdvanceResult(mesh, boundary, damage, displacement, 0.0, False,
                                      angle_error_deg=float(angle_error), reason=f"topology_error:{exc}")

        # Extend nodal damage for duplicated nodes. New CZM crack surfaces are
        # represented topologically, not by bulk stiffness killing.
        dnew = np.asarray(damage, float)
        if new_mesh.nn > len(dnew):
            # duplicate endpoint nodes inherit the nearest original damage value
            extra = []
            for x in new_mesh.nodes[len(dnew):]:
                j = int(np.argmin(np.linalg.norm(mesh.nodes - x[None, :], axis=1)))
                extra.append(float(dnew[j]))
            dnew = np.concatenate([dnew, np.asarray(extra, float)])
        new_bnd = make_boundary_data(new_mesh, self.geom)
        moved = float(np.linalg.norm(q - p0_mesh))
        self.advance_log.append({
            "front_id": int(front_id), "event_index": int(ce.event_index),
            "x0": float(p0_mesh[0]), "y0": float(p0_mesh[1]),
            "x1": float(q[0]), "y1": float(q[1]),
            "length_m": moved, "angle_error_deg": float(angle_error),
            "damage": float(ce.damage), "reason": "ok",
        })
        return CrackAdvanceResult(new_mesh, new_bnd, dnew, new_u, moved, True,
                                  angle_error_deg=float(angle_error),
                                  selected_edge_length=moved, reason="ok")

    def register_branch_front(self, parent_id: int, child_id: int, xy: np.ndarray) -> None:
        if int(parent_id) in self.tip_nodes:
            a, b, _ = self.tip_nodes[int(parent_id)]
            self.tip_nodes[int(child_id)] = (int(a), int(b), np.asarray(xy, float).copy())

    def write_diagnostics(self, out_dir: str) -> None:
        import csv
        import json
        from pathlib import Path
        out = Path(out_dir)
        out.mkdir(parents=True, exist_ok=True)
        rows = self.cohesive_network.to_rows()
        np.savetxt(out / "cohesive_elements.csv", rows, delimiter=",",
                   header=("index,front_id,event_index,plus0,plus1,minus0,minus1,"
                           "length_m,damage,clock,tx,ty,nx,ny"), comments="")
        with (out / "czm_advance_log.json").open("w") as f:
            json.dump(self.advance_log, f, indent=2)


class AdaptiveCZMBackend(EdgeSplitCZMBackend):
    """Angle-faithful CZM backend using local tip-node r-adaptation.

    Before the topology split, one existing geometric neighbor of the active
    tip is relocated onto the *exact* hazard-selected ray at the requested
    physical advance length.  Only the one-ring star of that vertex changes
    geometry; element connectivity, element count, and element ordering remain
    unchanged.  Consequently all element/Gauss-point history arrays retain
    their indices.  The move is accepted only when no incident triangle flips
    orientation and local element quality remains above configured limits.

    This is deliberately an r-adaptive local patch update rather than global
    remeshing.  It removes mesh-edge steering while preserving the mature
    history-rich fatigue/plasticity state without projection.
    """

    name = "adaptive_czm"

    def __init__(
        self,
        *,
        geom,
        penalty_normal_Pa_per_m: float = 1.0e18,
        penalty_tangent_Pa_per_m: float = 1.0e18,
        max_angle_error_deg: float = 35.0,
        event_damage: float = 1.0,
        min_area_ratio: float = 0.08,
        min_triangle_quality: float = 0.035,
        max_node_move_factor: float = 1.75,
        max_hrefine_subsegments: int = 512,
    ) -> None:
        super().__init__(
            geom=geom,
            penalty_normal_Pa_per_m=penalty_normal_Pa_per_m,
            penalty_tangent_Pa_per_m=penalty_tangent_Pa_per_m,
            max_angle_error_deg=max_angle_error_deg,
            event_damage=event_damage,
        )
        self.min_area_ratio = float(max(min_area_ratio, 1e-6))
        self.min_triangle_quality = float(max(min_triangle_quality, 1e-8))
        self.max_node_move_factor = float(max(max_node_move_factor, 0.0))
        self.max_hrefine_subsegments = int(max(max_hrefine_subsegments, 1))

    @staticmethod
    def _signed_twice_area(nodes: np.ndarray, elems: np.ndarray) -> np.ndarray:
        X = nodes[elems]
        return ((X[:, 1, 0] - X[:, 0, 0]) * (X[:, 2, 1] - X[:, 0, 1])
                - (X[:, 1, 1] - X[:, 0, 1]) * (X[:, 2, 0] - X[:, 0, 0]))

    @staticmethod
    def _triangle_quality(nodes: np.ndarray, elems: np.ndarray) -> np.ndarray:
        """Return 4*sqrt(3)*A/sum(l_i^2), in [0,1] for valid triangles."""
        X = nodes[elems]
        l2 = (np.sum((X[:, 1] - X[:, 0]) ** 2, axis=1)
              + np.sum((X[:, 2] - X[:, 1]) ** 2, axis=1)
              + np.sum((X[:, 0] - X[:, 2]) ** 2, axis=1))
        a2 = np.abs(AdaptiveCZMBackend._signed_twice_area(nodes, elems))
        return 2.0 * np.sqrt(3.0) * a2 / np.maximum(l2, 1e-300)

    @staticmethod
    def _interp_nodal_vector(mesh: TriMesh, values: np.ndarray, point: np.ndarray) -> np.ndarray:
        """Piecewise-linear interpolation at ``point`` with nearest-node fallback."""
        p = np.asarray(point, float)
        vals = np.asarray(values, float)
        X = mesh.nodes[mesh.elems]
        # Bounding-box prefilter keeps this cheap even for long-growth meshes.
        eps = max(1e-14, 1e-9 * max(mesh.hbar_tip, mesh.hbar, 1e-12))
        mask = ((X[:, :, 0].min(axis=1) - eps <= p[0])
                & (X[:, :, 0].max(axis=1) + eps >= p[0])
                & (X[:, :, 1].min(axis=1) - eps <= p[1])
                & (X[:, :, 1].max(axis=1) + eps >= p[1]))
        for e in np.where(mask)[0]:
            tri = X[e]
            M = np.array([[tri[0, 0], tri[1, 0], tri[2, 0]],
                          [tri[0, 1], tri[1, 1], tri[2, 1]],
                          [1.0, 1.0, 1.0]])
            rhs = np.array([p[0], p[1], 1.0])
            try:
                w = np.linalg.solve(M, rhs)
            except np.linalg.LinAlgError:
                continue
            if np.min(w) >= -1e-8 and np.max(w) <= 1.0 + 1e-8:
                return w @ vals[mesh.elems[e]]
        j = int(np.argmin(np.linalg.norm(mesh.nodes - p[None, :], axis=1)))
        return vals[j].copy()

    def _candidate_neighbor_points(self, mesh: TriMesh, p0: np.ndarray, front_id: int):
        tip_ids = self._tip_geometric_node_ids(mesh, p0, front_id)
        tip_set = set(int(i) for i in tip_ids)
        incident = self._incident_elements(mesh.elems, tip_set)
        neigh = set()
        for e in incident:
            for nid in mesh.elems[e]:
                if int(nid) not in tip_set:
                    neigh.add(int(nid))
        tol = max(1e-12, 1e-6 * max(mesh.hbar_tip, mesh.hbar, 1e-12))
        out = []
        seen = set()
        for nid in neigh:
            q = mesh.nodes[nid]
            key = tuple(np.round(q / tol).astype(np.int64))
            if key in seen:
                continue
            seen.add(key)
            ids = self._coincident_ids(mesh.nodes, q, tol)
            out.append((q.copy(), ids))
        return out

    def _ray_exit_edge(
        self,
        mesh: TriMesh,
        p0: np.ndarray,
        direction: np.ndarray,
        front_id: int,
        max_distance: float,
    ):
        """Return the first mesh edge crossed by an exact forward ray.

        The fast path searches the geometric one-ring around the current tip.
        After repeated node duplication/refinement, however, the topological
        one-ring can contain split tip copies whose triangle fan no longer
        satisfies the old ``one tip node + two opposite nodes`` assumption.
        A geometry-based edge scan is therefore used, with a global fallback
        limited to the requested physical segment length.
        """
        p0 = np.asarray(p0, float)
        d = np.asarray(direction, float)
        d /= max(float(np.linalg.norm(d)), 1e-300)
        max_distance = float(max(max_distance, 0.0))
        geom_tol = max(1e-12, 1e-7 * max(mesh.hbar_tip, mesh.hbar, 1e-12))
        t_eps = max(1e-12, 1e-8 * max(max_distance, mesh.hbar_tip, 1e-12))

        tip_ids = self._tip_geometric_node_ids(mesh, p0, int(front_id))
        incident = self._incident_elements(mesh.elems, set(int(i) for i in tip_ids)) \
            if tip_ids else np.zeros(0, dtype=int)

        def unique_edges(element_ids):
            seen = set()
            out = []
            for e in np.asarray(element_ids, dtype=int):
                conn = [int(x) for x in mesh.elems[int(e)]]
                for ia, ib in ((0, 1), (1, 2), (2, 0)):
                    i, j = conn[ia], conn[ib]
                    a = mesh.nodes[i]; b = mesh.nodes[j]
                    # Deduplicate coincident topological edges geometrically.
                    ka = tuple(np.round(a / geom_tol).astype(np.int64))
                    kb = tuple(np.round(b / geom_tol).astype(np.int64))
                    key = tuple(sorted((ka, kb)))
                    if key in seen:
                        continue
                    seen.add(key)
                    out.append((i, j, int(e)))
            return out

        def intersections(edges):
            candidates = []
            for i, j, e in edges:
                a = mesh.nodes[int(i)]; b = mesh.nodes[int(j)]
                edge = b - a
                M = np.array([[d[0], -edge[0]], [d[1], -edge[1]]], dtype=float)
                rhs = a - p0
                try:
                    t, xi = np.linalg.solve(M, rhs)
                except np.linalg.LinAlgError:
                    # Collinear edge: use the nearest forward endpoint as a
                    # vertex hit rather than declaring that the ray has no exit.
                    rel_a = a - p0; rel_b = b - p0
                    cross_a = abs(d[0] * rel_a[1] - d[1] * rel_a[0])
                    cross_b = abs(d[0] * rel_b[1] - d[1] * rel_b[0])
                    if cross_a <= geom_tol and cross_b <= geom_tol:
                        ta = float(rel_a @ d); tb = float(rel_b @ d)
                        forward = [(ta, 0.0, a), (tb, 1.0, b)]
                        forward = [z for z in forward if z[0] > t_eps and z[0] <= max_distance + geom_tol]
                        if forward:
                            tt, xxi, qq = min(forward, key=lambda z: z[0])
                            candidates.append((tt, int(i), int(j), qq.copy(), int(e), xxi))
                    continue
                if t > t_eps and t <= max_distance + geom_tol \
                        and xi >= -1e-9 and xi <= 1.0 + 1e-9:
                    q = p0 + float(t) * d
                    candidates.append((float(t), int(i), int(j), q, int(e), float(xi)))
            if not candidates:
                return None
            candidates.sort(key=lambda z: z[0])
            return candidates[0]

        # Prefer the local fan.
        hit = intersections(unique_edges(incident)) if len(incident) else None
        if hit is not None:
            return hit

        # Robust fallback: scan only edges whose bounding boxes overlap the
        # requested ray-segment bounding box.  This handles split/duplicated tip
        # fans without making crack propagation depend on topological valence.
        qend = p0 + max_distance * d
        lo = np.minimum(p0, qend) - geom_tol
        hi = np.maximum(p0, qend) + geom_tol
        X = mesh.nodes[mesh.elems]
        mask = ((X[:, :, 0].max(axis=1) >= lo[0])
                & (X[:, :, 0].min(axis=1) <= hi[0])
                & (X[:, :, 1].max(axis=1) >= lo[1])
                & (X[:, :, 1].min(axis=1) <= hi[1]))
        return intersections(unique_edges(np.where(mask)[0]))

    def _insert_point_on_edge(
        self,
        mesh: TriMesh,
        displacement: np.ndarray,
        q: np.ndarray,
        edge_i: int,
        edge_j: int,
    ):
        """Split all triangles sharing one topological edge at point ``q``."""
        nodes0 = np.asarray(mesh.nodes, float)
        elems0 = np.asarray(mesh.elems, int)
        u0 = np.asarray(displacement, float).reshape(-1, 2)
        i = int(edge_i); j = int(edge_j)
        parents = np.where(np.any(elems0 == i, axis=1) & np.any(elems0 == j, axis=1))[0]
        if len(parents) == 0:
            return None, None, "edge_has_no_parent_elements", {}, None

        new_id = int(len(nodes0))
        nodes1 = np.vstack([nodes0, np.asarray(q, float)[None, :]])
        edge = nodes0[j] - nodes0[i]
        den = float(edge @ edge)
        xi = 0.5 if den <= 1e-300 else float(np.clip(((q - nodes0[i]) @ edge) / den, 0.0, 1.0))
        uq = (1.0 - xi) * u0[i] + xi * u0[j]
        u1 = np.vstack([u0, uq[None, :]])

        elems1 = elems0.copy()
        appended = []
        appended_parent = []
        for parent in parents:
            conn = [int(x) for x in elems0[int(parent)]]
            k = [nid for nid in conn if nid not in (i, j)]
            if len(k) != 1:
                return None, None, "invalid_edge_parent_triangle", {}, None
            k = int(k[0])
            old_sign = float(self._signed_twice_area(nodes0, elems0[[int(parent)]])[0])
            kids = [np.array([i, new_id, k], dtype=int),
                    np.array([new_id, j, k], dtype=int)]
            for kid in kids:
                sgn = float(self._signed_twice_area(nodes1, kid.reshape(1, 3))[0])
                if old_sign * sgn < 0.0:
                    kid[0], kid[1] = kid[1], kid[0]
            elems1[int(parent)] = kids[0]
            appended.append(kids[1])
            appended_parent.append(int(parent))

        elems1 = np.vstack([elems1, np.asarray(appended, dtype=int)])
        refined = rebuild_tri_mesh(nodes1, elems1, tip_centers=[q])
        parent_map = np.concatenate([np.arange(mesh.ne, dtype=int),
                                     np.asarray(appended_parent, dtype=int)])
        new_ids = [int(p) for p in parents] + list(range(mesh.ne, refined.ne))
        qmin = float(np.min(self._triangle_quality(refined.nodes, refined.elems[new_ids])))
        meta = {
            "n_new_bulk_elements": int(len(appended)),
            "min_triangle_quality": qmin,
            "min_area_ratio": 1.0,
            "node_move_m": 0.0,
            "n_moved_geometric_copies": 0,
            "n_incident_elements": int(len(parents)),
            "split_edge_i": i,
            "split_edge_j": j,
        }
        return refined, u1.reshape(-1), "ok", meta, parent_map

    def _insert_target_in_incident_triangle(
        self,
        mesh: TriMesh,
        displacement: np.ndarray,
        p0: np.ndarray,
        target: np.ndarray,
        front_id: int,
    ):
        """Insert ``target`` into a tip-adjacent triangle and split it locally.

        This is the quality-safe fallback when repeated r-adaptation has
        exhausted the one-ring geometry.  The requested target lies on the
        exact hazard-selected ray.  If it is inside a triangle incident to one
        of the coincident tip-node copies, that parent triangle is replaced by
        three children.  Only two new bulk elements are appended.

        Returns a new mesh, interpolated displacement, target point, status,
        metadata, and ``parent_map`` where ``parent_map[e_new]`` identifies the
        old element whose integration-point history should be inherited.
        """
        p0 = np.asarray(p0, float)
        target = np.asarray(target, float)
        tip_ids = self._tip_geometric_node_ids(mesh, p0, int(front_id))
        if not tip_ids:
            return None, None, None, "tip_node_not_found", {}, None
        incident = self._incident_elements(mesh.elems, set(int(i) for i in tip_ids))
        if len(incident) == 0:
            return None, None, None, "no_tip_incident_elements", {}, None

        # Locate target by barycentric coordinates.  Prefer an element where the
        # target is comfortably interior; this avoids nearly degenerate child
        # triangles when the point lies numerically on an existing edge.
        choices = []
        for e in incident:
            conn = mesh.elems[int(e)]
            tri = mesh.nodes[conn]
            M = np.array([[tri[0, 0], tri[1, 0], tri[2, 0]],
                          [tri[0, 1], tri[1, 1], tri[2, 1]],
                          [1.0, 1.0, 1.0]], dtype=float)
            rhs = np.array([target[0], target[1], 1.0], dtype=float)
            try:
                w = np.linalg.solve(M, rhs)
            except np.linalg.LinAlgError:
                continue
            if np.min(w) >= -1e-10 and np.max(w) <= 1.0 + 1e-10:
                choices.append((float(np.min(w)), int(e), conn.copy(), w))
        if not choices:
            return None, None, None, "target_not_in_tip_incident_triangle", {}, None
        choices.sort(key=lambda z: z[0], reverse=True)
        _, parent, conn, w = choices[0]

        nodes0 = np.asarray(mesh.nodes, float)
        elems0 = np.asarray(mesh.elems, int)
        u0 = np.asarray(displacement, float).reshape(-1, 2)
        new_id = int(len(nodes0))
        nodes1 = np.vstack([nodes0, target[None, :]])
        u_target = w @ u0[conn]
        u1 = np.vstack([u0, u_target[None, :]])

        a, b, c = [int(x) for x in conn]
        children = np.array([[a, b, new_id],
                             [b, c, new_id],
                             [c, a, new_id]], dtype=int)
        old_sign = float(self._signed_twice_area(nodes0, elems0[[parent]])[0])
        for k in range(3):
            sgn = float(self._signed_twice_area(nodes1, children[[k]])[0])
            if old_sign * sgn < 0.0:
                children[k, 0], children[k, 1] = children[k, 1], children[k, 0]

        elems1 = elems0.copy()
        elems1[parent] = children[0]
        elems1 = np.vstack([elems1, children[1], children[2]])
        refined = rebuild_tri_mesh(nodes1, elems1, tip_centers=[target])

        parent_map = np.concatenate([np.arange(mesh.ne, dtype=int),
                                     np.array([parent, parent], dtype=int)])
        qnew = self._triangle_quality(refined.nodes, refined.elems[[parent, mesh.ne, mesh.ne + 1]])
        meta = {
            "refined_parent_element": int(parent),
            "n_new_bulk_elements": 2,
            "min_triangle_quality": float(np.min(qnew)),
            "min_area_ratio": 1.0,
            "node_move_m": 0.0,
            "n_moved_geometric_copies": 0,
            "n_incident_elements": int(len(incident)),
        }
        return refined, u1.reshape(-1), target.copy(), "ok", meta, parent_map

    def _steer_neighbor_to_target(
        self,
        mesh: TriMesh,
        displacement: np.ndarray,
        p0: np.ndarray,
        target: np.ndarray,
        direction: np.ndarray,
        front_id: int,
    ):
        """Relocate the safest tip-neighbor geometric vertex to ``target``."""
        nodes0 = mesh.nodes
        elems = mesh.elems
        u0 = np.asarray(displacement, float).reshape(-1, 2)
        req_L = float(np.linalg.norm(target - p0))
        if req_L <= 1e-14:
            return None, None, None, "zero_requested_length", {}

        # Target must remain inside the physical specimen.
        if not (0.0 <= target[0] <= float(self.geom.Lx)
                and -0.5 * float(self.geom.Ly) <= target[1] <= 0.5 * float(self.geom.Ly)):
            return None, None, None, "target_outside_domain", {}

        q_candidates = self._candidate_neighbor_points(mesh, p0, front_id)
        if not q_candidates:
            return None, None, None, "no_tip_neighbors", {}

        candidates = []
        for qold, ids in q_candidates:
            v = qold - p0
            Lold = float(np.linalg.norm(v))
            if Lold <= 1e-14:
                continue
            forward = float(np.dot(v / Lold, direction))
            if forward <= 1e-8:
                continue
            move = float(np.linalg.norm(target - qold))
            if self.max_node_move_factor > 0.0 and move > self.max_node_move_factor * req_L:
                continue

            moved_set = set(int(i) for i in ids)
            incident = self._incident_elements(elems, moved_set)
            if len(incident) == 0:
                continue
            trial_nodes = nodes0.copy()
            trial_nodes[np.asarray(ids, dtype=int)] = target[None, :]
            old_s = self._signed_twice_area(nodes0, elems[incident])
            new_s = self._signed_twice_area(trial_nodes, elems[incident])
            # Preserve orientation and reject near-collapse.
            if np.any(old_s * new_s <= 0.0):
                continue
            ratio = np.abs(new_s) / np.maximum(np.abs(old_s), 1e-300)
            if float(np.min(ratio)) < self.min_area_ratio:
                continue
            qnew = self._triangle_quality(trial_nodes, elems[incident])
            qmin = float(np.min(qnew))
            if qmin < self.min_triangle_quality:
                continue
            # Prefer less mesh motion, then better worst-element quality.
            candidates.append((move / req_L, -qmin, qold, ids, incident, trial_nodes, ratio, qnew))

        if not candidates:
            return None, None, None, "no_quality_safe_exact_ray_move", {}
        candidates.sort(key=lambda z: (z[0], z[1]))
        _, _, qold, ids, incident, nodes1, ratio, qnew = candidates[0]

        u1 = u0.copy()
        u_target = self._interp_nodal_vector(mesh, u0, target)
        u1[np.asarray(ids, dtype=int)] = u_target[None, :]
        meta = {
            "steered_from_x": float(qold[0]),
            "steered_from_y": float(qold[1]),
            "node_move_m": float(np.linalg.norm(target - qold)),
            "min_area_ratio": float(np.min(ratio)),
            "min_triangle_quality": float(np.min(qnew)),
            "n_moved_geometric_copies": int(len(ids)),
            "n_incident_elements": int(len(incident)),
        }
        steered = rebuild_tri_mesh(nodes1, elems, tip_centers=[target])
        return steered, u1.reshape(-1), target.copy(), "ok", meta

    def _transaction_snapshot(self):
        """Capture mutable backend state for one physical crack event.

        Multi-triangle exact-ray marching may add several cohesive subsegments
        before the final requested physical endpoint is reached.  The complete
        physical Arrhenius renewal must be atomic: if any subsegment fails, the
        mesh returned to the solver is the original mesh and the cohesive
        network/backend bookkeeping must be rolled back to the same state.
        """
        tips = {
            int(fid): (int(a), int(b), np.asarray(xy, float).copy())
            for fid, (a, b, xy) in self.tip_nodes.items()
        }
        return {
            "n_cohesive": len(self.cohesive_network.elements),
            "tip_nodes": tips,
            "event_counter": int(self.event_counter),
            "n_log": len(self.advance_log),
        }

    def _transaction_rollback(self, snap) -> None:
        del self.cohesive_network.elements[int(snap["n_cohesive"]):]
        self.tip_nodes = {
            int(fid): (int(a), int(b), np.asarray(xy, float).copy())
            for fid, (a, b, xy) in snap["tip_nodes"].items()
        }
        self.event_counter = int(snap["event_counter"])
        del self.advance_log[int(snap["n_log"]):]

    def advance(
        self,
        *, mesh: TriMesh, boundary, damage: np.ndarray, displacement: np.ndarray,
        p0: np.ndarray, p1: np.ndarray, direction: np.ndarray, front_id: int,
        **kwargs,
    ) -> CrackAdvanceResult:
        depth = int(kwargs.get("_subdepth", 0) or 0)
        if depth > 0:
            return self._advance_impl(
                mesh=mesh, boundary=boundary, damage=damage, displacement=displacement,
                p0=p0, p1=p1, direction=direction, front_id=front_id, **kwargs)

        snap = self._transaction_snapshot()
        result = self._advance_impl(
            mesh=mesh, boundary=boundary, damage=damage, displacement=displacement,
            p0=p0, p1=p1, direction=direction, front_id=front_id, **kwargs)
        if not result.inserted:
            self._transaction_rollback(snap)
            # Return the authoritative pre-event state.  A failed physical
            # renewal must not expose a partially refined mesh/displacement.
            return CrackAdvanceResult(
                mesh, boundary, damage, displacement, 0.0, False,
                angle_error_deg=result.angle_error_deg, reason=result.reason,
                elem_parent_map=None)
        return result

    def _advance_impl(
        self,
        *, mesh: TriMesh, boundary, damage: np.ndarray, displacement: np.ndarray,
        p0: np.ndarray, p1: np.ndarray, direction: np.ndarray, front_id: int,
        **kwargs,
    ) -> CrackAdvanceResult:
        p0 = np.asarray(p0, float)
        p1 = np.asarray(p1, float)
        direction = np.asarray(direction, float)
        nd = float(np.linalg.norm(direction))
        if nd <= 0.0:
            return CrackAdvanceResult(mesh, boundary, damage, displacement, 0.0, False,
                                      angle_error_deg=180.0, reason="zero_direction")
        direction = direction / nd

        tip_ids = self._tip_geometric_node_ids(mesh, p0, int(front_id))
        if not tip_ids:
            return CrackAdvanceResult(mesh, boundary, damage, displacement, 0.0, False,
                                      angle_error_deg=180.0, reason="tip_node_not_found")
        p0_mesh = mesh.nodes[int(tip_ids[0])].copy()
        req_L = float(np.linalg.norm(p1 - p0))
        target = p0_mesh + req_L * direction

        parent_map = None
        geometry_update = "local_r_adapt_exact_ray"
        try:
            steered_mesh, steered_u, q, reason, meta = self._steer_neighbor_to_target(
                mesh, displacement, p0_mesh, target, direction, int(front_id))
        except Exception as exc:
            return CrackAdvanceResult(mesh, boundary, damage, displacement, 0.0, False,
                                      angle_error_deg=180.0, reason=f"steering_error:{exc}")

        # Repeated vertex steering eventually exhausts one-ring mesh quality.
        # When that happens, insert the exact target point into the adjacent
        # triangle and locally h-refine it instead of consuming a physical
        # fracture event without geometric advance.
        if steered_mesh is None:
            try:
                steered_mesh, steered_u, q, href_reason, href_meta, parent_map = \
                    self._insert_target_in_incident_triangle(
                        mesh, displacement, p0_mesh, target, int(front_id))
            except Exception as exc:
                return CrackAdvanceResult(mesh, boundary, damage, displacement, 0.0, False,
                                          angle_error_deg=180.0,
                                          reason=f"local_hrefine_error:{exc}")
            if steered_mesh is None:
                # The requested physical increment can extend beyond the current
                # one-ring after previous refinement.  March exactly to the first
                # crossed one-ring edge, split that edge locally, insert a CZM
                # subsegment, then recurse on the remaining collinear distance.
                depth = int(kwargs.get("_subdepth", 0) or 0)
                if depth >= self.max_hrefine_subsegments:
                    return CrackAdvanceResult(mesh, boundary, damage, displacement, 0.0, False,
                                              angle_error_deg=180.0,
                                              reason=(f"local_hrefine_recursion_limit:"
                                                      f"{depth}/{self.max_hrefine_subsegments}"))
                hit = self._ray_exit_edge(mesh, p0_mesh, direction, int(front_id), req_L)
                if hit is None:
                    return CrackAdvanceResult(mesh, boundary, damage, displacement, 0.0, False,
                                              angle_error_deg=180.0,
                                              reason=f"{reason};hrefine:{href_reason};no_ray_exit")
                t_hit, edge_i, edge_j, q_hit, _, xi_hit = hit
                if t_hit >= req_L - 1e-10:
                    return CrackAdvanceResult(mesh, boundary, damage, displacement, 0.0, False,
                                              angle_error_deg=180.0,
                                              reason=f"{reason};hrefine:{href_reason};exit_at_target")

                # When the exact ray crosses an existing mesh vertex, do not
                # split an edge at its endpoint (which would create a zero-length
                # child).  Use the existing vertex directly and let the cohesive
                # topology splitter duplicate it as needed.
                vertex_tol = 1e-8
                if xi_hit <= vertex_tol or xi_hit >= 1.0 - vertex_tol:
                    vid = int(edge_i if xi_hit <= vertex_tol else edge_j)
                    q_hit = mesh.nodes[vid].copy()
                    edge_mesh = mesh
                    edge_u = displacement
                    edge_reason = "ok"
                    edge_map = np.arange(mesh.ne, dtype=int)
                    edge_meta = {
                        "n_new_bulk_elements": 0,
                        "min_triangle_quality": float("nan"),
                        "min_area_ratio": 1.0,
                        "node_move_m": 0.0,
                        "n_moved_geometric_copies": 0,
                        "n_incident_elements": 0,
                        "ray_hit_existing_vertex": vid,
                    }
                else:
                    edge_mesh, edge_u, edge_reason, edge_meta, edge_map = self._insert_point_on_edge(
                        mesh, displacement, q_hit, edge_i, edge_j)
                    if edge_mesh is None:
                        return CrackAdvanceResult(mesh, boundary, damage, displacement, 0.0, False,
                                                  angle_error_deg=180.0, reason=edge_reason)
                try:
                    first_mesh, first_u, ce_first = self._split_segment_topology(
                        edge_mesh, edge_u, p0_mesh, q_hit, int(front_id))
                except Exception as exc:
                    return CrackAdvanceResult(mesh, boundary, damage, displacement, 0.0, False,
                                              angle_error_deg=0.0,
                                              reason=f"topology_error_edge_march:{exc}")

                dfirst = np.asarray(damage, float)
                if first_mesh.nn > len(dfirst):
                    extra = []
                    for xx in first_mesh.nodes[len(dfirst):]:
                        jn = int(np.argmin(np.linalg.norm(mesh.nodes - xx[None, :], axis=1)))
                        extra.append(float(dfirst[jn]))
                    dfirst = np.concatenate([dfirst, np.asarray(extra, float)])
                bfirst = make_boundary_data(first_mesh, self.geom)
                moved_first = float(np.linalg.norm(q_hit - p0_mesh))
                row = {
                    "front_id": int(front_id), "event_index": int(ce_first.event_index),
                    "x0": float(p0_mesh[0]), "y0": float(p0_mesh[1]),
                    "x1": float(q_hit[0]), "y1": float(q_hit[1]),
                    "length_m": moved_first, "angle_error_deg": 0.0,
                    "damage": float(ce_first.damage), "reason": "ok",
                    "geometry_update": "local_edge_march_exact_ray",
                }
                row.update(edge_meta)
                self.advance_log.append(row)

                rem = req_L - moved_first
                final_target = p0_mesh + req_L * direction
                rr2 = self.advance(
                    mesh=first_mesh, boundary=bfirst, damage=dfirst, displacement=first_u,
                    p0=q_hit, p1=final_target, direction=direction, front_id=int(front_id),
                    _subdepth=depth + 1, **{k: v for k, v in kwargs.items() if k != "_subdepth"})
                if not rr2.inserted:
                    return CrackAdvanceResult(
                        rr2.mesh, rr2.boundary, rr2.damage, rr2.displacement, moved_first, False,
                        angle_error_deg=rr2.angle_error_deg, reason=rr2.reason,
                        elem_parent_map=edge_map)
                if rr2.elem_parent_map is None:
                    composed = edge_map
                else:
                    composed = edge_map[np.asarray(rr2.elem_parent_map, dtype=int)]
                return CrackAdvanceResult(
                    rr2.mesh, rr2.boundary, rr2.damage, rr2.displacement,
                    moved_first + float(rr2.moved), True, angle_error_deg=0.0,
                    selected_edge_length=moved_first + float(rr2.moved), reason="ok",
                    elem_parent_map=composed)
            meta = href_meta
            geometry_update = "local_h_refine_exact_ray"

        try:
            new_mesh, new_u, ce = self._split_segment_topology(
                steered_mesh, steered_u, p0_mesh, q, int(front_id))
        except Exception as exc:
            return CrackAdvanceResult(mesh, boundary, damage, displacement, 0.0, False,
                                      angle_error_deg=0.0, reason=f"topology_error:{exc}")

        dnew = np.asarray(damage, float)
        if new_mesh.nn > len(dnew):
            extra = []
            for x in new_mesh.nodes[len(dnew):]:
                j = int(np.argmin(np.linalg.norm(mesh.nodes - x[None, :], axis=1)))
                extra.append(float(dnew[j]))
            dnew = np.concatenate([dnew, np.asarray(extra, float)])
        new_bnd = make_boundary_data(new_mesh, self.geom)
        moved = float(np.linalg.norm(q - p0_mesh))
        row = {
            "front_id": int(front_id), "event_index": int(ce.event_index),
            "x0": float(p0_mesh[0]), "y0": float(p0_mesh[1]),
            "x1": float(q[0]), "y1": float(q[1]),
            "length_m": moved, "angle_error_deg": 0.0,
            "damage": float(ce.damage), "reason": "ok",
            "geometry_update": geometry_update,
        }
        row.update(meta)
        self.advance_log.append(row)
        return CrackAdvanceResult(new_mesh, new_bnd, dnew, new_u, moved, True,
                                  angle_error_deg=0.0,
                                  selected_edge_length=moved, reason="ok",
                                  elem_parent_map=parent_map)


def build_crack_backend(args, geom) -> CrackBackend:
    kind = str(getattr(args, "crack_backend", "sharp_wake") or "sharp_wake").lower()
    if kind in ("sharp", "sharp_wake", "legacy"):
        return SharpWakeBackend()
    common = dict(
        geom=geom,
        penalty_normal_Pa_per_m=float(getattr(args, "czm_penalty_normal", 1.0e18)),
        penalty_tangent_Pa_per_m=float(getattr(args, "czm_penalty_tangent", 1.0e18)),
        max_angle_error_deg=float(getattr(args, "czm_max_angle_error_deg", 35.0)),
        event_damage=float(getattr(args, "czm_event_damage", 1.0)),
    )
    if kind in ("edge_czm", "edge_split_czm"):
        return EdgeSplitCZMBackend(**common)
    if kind == "adaptive_czm":
        return AdaptiveCZMBackend(
            **common,
            min_area_ratio=float(getattr(args, "czm_min_area_ratio", 0.08)),
            min_triangle_quality=float(getattr(args, "czm_min_triangle_quality", 0.035)),
            max_node_move_factor=float(getattr(args, "czm_max_node_move_factor", 1.75)),
            max_hrefine_subsegments=int(getattr(args, "czm_max_hrefine_subsegments", 512)),
        )
    raise ValueError(f"unknown crack backend: {kind}")
