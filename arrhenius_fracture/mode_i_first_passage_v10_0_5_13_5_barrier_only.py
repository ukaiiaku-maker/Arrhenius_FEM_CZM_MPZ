"""v10.0.5.13.5 long-corridor startup repair.

The 100 um production corridor exposed a geometric defect in the inherited
multi-cloud physical-refinement mesh.  Its rounded-coordinate deduplication can
leave points closer than the nominal tolerance when they fall on opposite sides
of a rounding-bin boundary.  Delaunay then creates sliver triangles even though
the underlying radial-ring placement is well shaped.

This point release changes only point deduplication:

* merge points by true Euclidean distance using a deterministic union-find;
* preserve exact explicit corridor centers and specimen-boundary points;
* retain the same 330 um radial-ring law, Delaunay triangulation, quality floor,
  FEM/CZM mechanics, barriers, MPZ state, hazards, and crack-growth law.
"""
from __future__ import annotations

import math
import sys
from typing import Any

import numpy as np
from scipy.spatial import Delaunay, cKDTree

from . import mode_i_first_passage_v9_18_5_3 as _v91853
from . import mode_i_first_passage_v10_0_5_13_barrier_only as _core
from . import mode_i_first_passage_v10_0_5_13_4_barrier_only as _base
from . import physical_refinement_mesh_v100510 as _phys

POINT_RELEASE = "10.0.5.13.5"
MODEL_ID = "FEM_CZM_full_2D_barrier_only_tip_source_long_corridor_v10_0_5_13_5"
PRODUCTION_MANIFEST = _base.PRODUCTION_MANIFEST


def expanded_candidate_counts_v1005135(length_um: float, max_gap_um: float) -> list[int]:
    """Search all deterministic center counts from two through the old upper bound."""
    base = max(
        2,
        int(math.ceil(max(float(length_um), 0.0) / max(float(max_gap_um), 1.0))) + 1,
    )
    return list(range(2, base + 4))


def _union_find_keep_indices_v1005135(
    points: np.ndarray,
    tolerance_m: float,
    protected_start: int,
    geom: Any,
) -> np.ndarray:
    """Return one deterministic representative for each radius-connected cluster."""
    pts = np.asarray(points, dtype=float)
    n = int(len(pts))
    if n == 0:
        return np.zeros(0, dtype=int)
    tol = max(float(tolerance_m), 1.0e-15)
    parent = np.arange(n, dtype=int)
    rank = np.zeros(n, dtype=np.int8)

    def find(i: int) -> int:
        j = int(i)
        while parent[j] != j:
            parent[j] = parent[parent[j]]
            j = int(parent[j])
        return j

    def union(i: int, j: int) -> None:
        ri, rj = find(i), find(j)
        if ri == rj:
            return
        if rank[ri] < rank[rj]:
            ri, rj = rj, ri
        parent[rj] = ri
        if rank[ri] == rank[rj]:
            rank[ri] += 1

    pairs = cKDTree(pts).query_pairs(r=tol, output_type="ndarray")
    for i, j in np.asarray(pairs, dtype=int).reshape(-1, 2):
        union(int(i), int(j))

    groups: dict[int, list[int]] = {}
    for i in range(n):
        groups.setdefault(find(i), []).append(i)

    boundary_tol = max(1.0e-14, 1.0e-8 * tol)

    def is_boundary(i: int) -> bool:
        x, y = pts[int(i)]
        return bool(
            abs(x) <= boundary_tol
            or abs(x - float(geom.Lx)) <= boundary_tol
            or abs(y + 0.5 * float(geom.Ly)) <= boundary_tol
            or abs(y - 0.5 * float(geom.Ly)) <= boundary_tol
        )

    keep: list[int] = []
    for members in groups.values():
        protected = [i for i in members if i >= int(protected_start)]
        if protected:
            chosen = min(protected)
        else:
            boundary = [i for i in members if is_boundary(i)]
            chosen = min(boundary) if boundary else min(members)
        keep.append(int(chosen))
    return np.asarray(sorted(keep), dtype=int)


def make_physical_refinement_mesh_v1005135(
    geom, mesh_cfg, seed=None, tip_center=None
):
    """Build the existing physical-refinement mesh with robust radius deduplication."""
    if _phys._ACTIVE_SPEC is None:
        raise RuntimeError("v10.0.5.13.5 physical refinement mesh was not configured")
    if seed is not None:
        np.random.seed(seed)

    h_fine = float(getattr(mesh_cfg, "tip_h_fine", 0.0) or 0.0)
    if h_fine <= 0.0:
        raise ValueError("v10.0.5.13.5 requires a graded tip mesh")
    ratio = float(getattr(mesh_cfg, "tip_ratio", 1.15) or 1.15)
    slope = max(ratio - 1.0, 0.02)
    h_far = max(float(geom.Lx), float(geom.Ly)) / 40.0
    centers = _phys._centers(geom, tip_center)
    spec = _phys._ACTIVE_SPEC.validate(geom, centers)

    clouds = [
        _phys._radial_ring_nodes(
            float(geom.Lx),
            float(geom.Ly),
            float(center[0]),
            float(center[1]),
            h_fine,
            slope,
            h_far,
            float(spec.radius_m),
        )
        for center in centers
    ]
    raw_nodes = np.vstack([*clouds, centers])
    protected_start = len(raw_nodes) - len(centers)
    tolerance_m = max(1.0e-12, 0.05 * h_fine)
    keep = _union_find_keep_indices_v1005135(
        raw_nodes, tolerance_m, protected_start, geom
    )
    nodes = raw_nodes[keep]

    triangulation = Delaunay(nodes)
    elems = triangulation.simplices
    centroids = nodes[elems].mean(axis=1)
    inside = (
        (centroids[:, 0] >= 0.0)
        & (centroids[:, 0] <= float(geom.Lx))
        & (centroids[:, 1] >= -0.5 * float(geom.Ly))
        & (centroids[:, 1] <= 0.5 * float(geom.Ly))
    )
    mesh = _phys.rebuild_tri_mesh(nodes, elems[inside], tip_centers=centers)

    for center in centers:
        distance = np.linalg.norm(mesh.nodes - center[None, :], axis=1)
        if float(np.min(distance)) > max(1.0e-12, 1.0e-8 * h_fine):
            raise RuntimeError(
                "v10.0.5.13.5 robust deduplication lost an explicit tip center: "
                f"{center.tolist()}"
            )

    mesh.production_refinement_radius_m = float(spec.radius_m)
    mesh.production_refinement_policy = (
        "fixed_physical_radius_same_radial_ring_law_radius_dedup_v1005135"
    )
    mesh.production_refinement_centers_m = centers.tolist()
    mesh.production_refinement_h_fine_m = h_fine
    mesh.production_refinement_tip_ratio = ratio
    mesh.production_refinement_dedup_tolerance_m = tolerance_m
    mesh.production_refinement_input_node_count = int(len(raw_nodes))
    mesh.production_refinement_output_node_count = int(len(nodes))
    mesh.production_refinement_removed_near_duplicates = int(len(raw_nodes) - len(nodes))
    return mesh


def main(argv: list[str] | None = None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    saved_counts = _v91853._candidate_counts
    saved_mesh = _core.make_physical_refinement_mesh_v100510
    _v91853._candidate_counts = expanded_candidate_counts_v1005135
    _core.make_physical_refinement_mesh_v100510 = make_physical_refinement_mesh_v1005135
    try:
        return _base.main(user_args)
    finally:
        _v91853._candidate_counts = saved_counts
        _core.make_physical_refinement_mesh_v100510 = saved_mesh


if __name__ == "__main__":
    main()


__all__ = [
    "POINT_RELEASE",
    "MODEL_ID",
    "PRODUCTION_MANIFEST",
    "expanded_candidate_counts_v1005135",
    "_union_find_keep_indices_v1005135",
    "make_physical_refinement_mesh_v1005135",
    "main",
]
