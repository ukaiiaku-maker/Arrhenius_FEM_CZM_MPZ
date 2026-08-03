"""Constrained exact-endpoint star triangulation for the v10.0.5.18.3.8 cavity.

The unconstrained five-point Delaunay triangulation may retain a hull-to-hull
internal diagonal, producing one cavity triangle that does not contain the
exact endpoint.  That topology is valid as a generic Delaunay mesh but is not
the desired crack-tip cavity: the exact endpoint must be the interior hub of
the four preserved cavity boundary edges.

This overlay replaces only the local triangulator used by the existing
v10.0.5.18.3.8 two-triangle cavity candidate.  The existing implementation
continues to enforce convexity, cavity-boundary preservation, direct connection
to a geometric tip copy, triangle quality, child-area ratio, state transfer,
and atomic rollback.
"""
from __future__ import annotations

import numpy as np
from scipy.spatial import ConvexHull, QhullError

from . import quality_aware_czm_two_triangle_cavity_v10051838 as _cavity


MODEL_ID = "constrained_exact_endpoint_star_v10_0_5_18_3_8"


class EndpointStarTriangulation:
    """Return the four-triangle fan from a convex quadrilateral to point 4."""

    def __init__(self, points):
        xy = np.asarray(points, dtype=float)
        if xy.shape != (5, 2):
            raise QhullError(
                "endpoint-star triangulation requires four cavity vertices "
                "plus one exact endpoint"
            )

        hull = ConvexHull(xy[:4])
        order = np.asarray(hull.vertices, dtype=int)
        if len(order) != 4:
            raise QhullError("cavity boundary is not a strict convex quadrilateral")

        # The target must be inside the unchanged cavity boundary.  The original
        # cavity implementation already locates it inside the target parent; the
        # hull test here makes the constrained triangulator independently safe.
        signed = hull.equations[:, :2] @ xy[4] + hull.equations[:, 2]
        scale = max(1.0, float(np.max(np.linalg.norm(xy[:4], axis=1))))
        if float(np.max(signed)) > 1.0e-10 * scale:
            raise QhullError("exact endpoint lies outside the convex cavity")

        self.simplices = np.asarray(
            [
                [int(order[i]), int(order[(i + 1) % 4]), 4]
                for i in range(4)
            ],
            dtype=int,
        )


def install() -> None:
    """Install the endpoint-star triangulator and the existing cavity fallback."""
    _cavity.Delaunay = EndpointStarTriangulation
    _cavity.install()


__all__ = ["EndpointStarTriangulation", "MODEL_ID", "install"]
