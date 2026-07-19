"""Audit-only mesh construction for the v10.0.5.8 fixed-grip benchmark.

The production graded mesh intentionally shrinks its refined disk with ``tip_h``.
That is efficient for crack-following simulations, but it confounds an h-convergence
study because both element size and the physical refinement extent change.  This
module holds the refinement radius fixed and inserts every perturbed crack-tip
position as an explicit mesh node.  It is installed only by the v10.0.5.8 audit
runner and does not change production fracture calculations.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np
from scipy.spatial import Delaunay

from .mesh import _radial_ring_nodes, rebuild_tri_mesh


@dataclass(frozen=True)
class FixedGripAuditMeshSpec:
    crack_center_m: float
    crack_offsets_m: tuple[float, ...]
    refinement_radius_m: float

    def validate(self, width_m: float, height_m: float) -> "FixedGripAuditMeshSpec":
        if not (0.0 < self.crack_center_m < float(width_m)):
            raise ValueError("audit crack center must lie inside the specimen")
        if self.refinement_radius_m <= 0.0:
            raise ValueError("audit refinement radius must be positive")
        for offset in self.crack_offsets_m:
            tip = self.crack_center_m + float(offset)
            if not (0.0 < tip < float(width_m)):
                raise ValueError("an audit crack perturbation leaves the specimen")
        if self.refinement_radius_m >= min(
            self.crack_center_m,
            float(width_m) - self.crack_center_m,
            0.5 * float(height_m),
        ):
            raise ValueError("audit refinement radius must close inside the specimen")
        return self


_ACTIVE_SPEC: FixedGripAuditMeshSpec | None = None


def configure_fixed_grip_audit_mesh(
    *,
    crack_center_m: float,
    crack_increments_m: Iterable[float],
    refinement_radius_m: float,
) -> FixedGripAuditMeshSpec:
    increments = sorted(set(abs(float(value)) for value in crack_increments_m))
    if not increments or any(value <= 0.0 for value in increments):
        raise ValueError("positive crack increments are required")
    offsets = tuple(sorted({0.0, *increments, *(-value for value in increments)}))
    spec = FixedGripAuditMeshSpec(
        crack_center_m=float(crack_center_m),
        crack_offsets_m=offsets,
        refinement_radius_m=float(refinement_radius_m),
    )
    global _ACTIVE_SPEC
    _ACTIVE_SPEC = spec
    return spec


def clear_fixed_grip_audit_mesh() -> None:
    global _ACTIVE_SPEC
    _ACTIVE_SPEC = None


def _configured_centers(geom) -> np.ndarray:
    if _ACTIVE_SPEC is None:
        raise RuntimeError("fixed-grip audit mesh was not configured")
    _ACTIVE_SPEC.validate(geom.Lx, geom.Ly)
    return np.asarray(
        [[_ACTIVE_SPEC.crack_center_m + offset, 0.0] for offset in _ACTIVE_SPEC.crack_offsets_m],
        dtype=float,
    )


def make_fixed_radius_audit_mesh(geom, mesh_cfg, seed=None, tip_center=None):
    """Build one common mesh for every crack perturbation at a given ``tip_h``.

    All requested crack tips are explicit nodes.  The physical refined radius is
    independent of ``tip_h``, so successive meshes change resolution without
    shrinking the J-domain support.
    """
    del tip_center  # the configured perturbation inventory is authoritative
    if seed is not None:
        np.random.seed(seed)
    h_fine = float(getattr(mesh_cfg, "tip_h_fine", 0.0) or 0.0)
    if h_fine <= 0.0:
        raise ValueError("the fixed-grip convergence audit requires a graded tip mesh")
    ratio = float(getattr(mesh_cfg, "tip_ratio", 1.15) or 1.15)
    slope = max(ratio - 1.0, 0.02)
    h_far = max(float(geom.Lx), float(geom.Ly)) / 40.0
    centers = _configured_centers(geom)
    radius = float(_ACTIVE_SPEC.refinement_radius_m)

    clouds = [
        _radial_ring_nodes(
            float(geom.Lx),
            float(geom.Ly),
            float(center[0]),
            float(center[1]),
            h_fine,
            slope,
            h_far,
            radius,
        )
        for center in centers
    ]
    nodes = np.vstack([*clouds, centers])
    key = np.round(nodes / max(1.0e-12, 0.05 * h_fine)).astype(np.int64)
    _, keep = np.unique(key, axis=0, return_index=True)
    nodes = nodes[np.sort(keep)]

    triangulation = Delaunay(nodes)
    elems = triangulation.simplices
    centroids = nodes[elems].mean(axis=1)
    inside = (
        (centroids[:, 0] >= 0.0)
        & (centroids[:, 0] <= float(geom.Lx))
        & (centroids[:, 1] >= -0.5 * float(geom.Ly))
        & (centroids[:, 1] <= 0.5 * float(geom.Ly))
    )
    mesh = rebuild_tri_mesh(nodes, elems[inside], tip_centers=centers)

    # Fail closed if duplicate filtering or triangulation lost an exact tip node.
    for center in centers:
        distance = np.linalg.norm(mesh.nodes - center[None, :], axis=1)
        if float(np.min(distance)) > max(1.0e-12, 1.0e-8 * h_fine):
            raise RuntimeError(f"audit mesh is missing explicit crack-tip node {center.tolist()}")
    return mesh


__all__ = [
    "FixedGripAuditMeshSpec",
    "configure_fixed_grip_audit_mesh",
    "clear_fixed_grip_audit_mesh",
    "make_fixed_radius_audit_mesh",
]
