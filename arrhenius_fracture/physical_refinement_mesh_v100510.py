"""Audit-only fixed physical refinement mesh for v10.0.5.10.

The production graded mesh normally sets its refined disk from ``40*h_tip``.
That couples J-domain support to the requested tip spacing. This module reuses the
same production radial-ring generator, coarsening law, Delaunay triangulation and
boundary clipping, but replaces only the disk extent with one explicit physical
radius. It is installed only by the v10.0.5.10 parity probe.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy.spatial import Delaunay

from .mesh import _radial_ring_nodes, rebuild_tri_mesh


@dataclass(frozen=True)
class PhysicalRefinementSpecV100510:
    radius_m: float

    def validate(self, geom: Any, centers: np.ndarray) -> "PhysicalRefinementSpecV100510":
        radius = float(self.radius_m)
        if not np.isfinite(radius) or radius <= 0.0:
            raise ValueError("v10.0.5.10 physical refinement radius must be finite and positive")
        for center in np.asarray(centers, dtype=float):
            clearance = min(
                float(center[0]),
                float(geom.Lx) - float(center[0]),
                0.5 * float(geom.Ly) - abs(float(center[1])),
            )
            if radius >= clearance:
                raise ValueError(
                    "v10.0.5.10 physical refinement disk must close inside the specimen: "
                    f"radius={radius:.16g}, center={center.tolist()}, clearance={clearance:.16g}"
                )
        return self


_ACTIVE_SPEC: PhysicalRefinementSpecV100510 | None = None


def configure_physical_refinement_v100510(radius_m: float) -> PhysicalRefinementSpecV100510:
    spec = PhysicalRefinementSpecV100510(float(radius_m))
    if not np.isfinite(spec.radius_m) or spec.radius_m <= 0.0:
        raise ValueError("v10.0.5.10 physical refinement radius must be finite and positive")
    global _ACTIVE_SPEC
    _ACTIVE_SPEC = spec
    return spec


def clear_physical_refinement_v100510() -> None:
    global _ACTIVE_SPEC
    _ACTIVE_SPEC = None


def _centers(geom: Any, tip_center: Any) -> np.ndarray:
    if tip_center is None:
        return np.asarray([[float(geom.a0), 0.0]], dtype=float)
    value = np.asarray(tip_center, dtype=float)
    return value.reshape(1, 2) if value.ndim == 1 else value[:, :2]


def make_physical_refinement_mesh_v100510(geom, mesh_cfg, seed=None, tip_center=None):
    """Build the production radial mesh with one fixed physical support radius."""
    if _ACTIVE_SPEC is None:
        raise RuntimeError("v10.0.5.10 physical refinement mesh was not configured")
    if seed is not None:
        np.random.seed(seed)

    h_fine = float(getattr(mesh_cfg, "tip_h_fine", 0.0) or 0.0)
    if h_fine <= 0.0:
        raise ValueError("v10.0.5.10 physical refinement audit requires a graded tip mesh")
    ratio = float(getattr(mesh_cfg, "tip_ratio", 1.15) or 1.15)
    slope = max(ratio - 1.0, 0.02)
    h_far = max(float(geom.Lx), float(geom.Ly)) / 40.0
    centers = _centers(geom, tip_center)
    spec = _ACTIVE_SPEC.validate(geom, centers)

    clouds = [
        _radial_ring_nodes(
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

    for center in centers:
        distance = np.linalg.norm(mesh.nodes - center[None, :], axis=1)
        if float(np.min(distance)) > max(1.0e-12, 1.0e-8 * h_fine):
            raise RuntimeError(
                "v10.0.5.10 physical refinement mesh lost an explicit tip center: "
                f"{center.tolist()}"
            )

    # Dynamic attributes are audit metadata only; TriMesh is intentionally not
    # changed, so production modules outside this probe retain their exact API.
    mesh.production_refinement_radius_m = float(spec.radius_m)
    mesh.production_refinement_policy = "fixed_physical_radius_same_radial_ring_law"
    mesh.production_refinement_centers_m = centers.tolist()
    mesh.production_refinement_h_fine_m = h_fine
    mesh.production_refinement_tip_ratio = ratio
    return mesh


__all__ = [
    "PhysicalRefinementSpecV100510",
    "configure_physical_refinement_v100510",
    "clear_physical_refinement_v100510",
    "make_physical_refinement_mesh_v100510",
]
