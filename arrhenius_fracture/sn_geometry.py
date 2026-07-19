"""Geometry utilities for S-N crack-initiation simulations.

The 2-D S-N driver starts from a *geometrically blunt* free-surface feature,
not from an initial sharp precrack.  The helper below generates a triangular
mesh of a rectangular plate with a smooth half-elliptical edge notch removed
from the left free surface.

For an ellipse x=a cos(theta), y=b sin(theta), theta in [-pi/2, pi/2], the
notch root is (a, 0) and its local radius of curvature is rho=b^2/a.  Choosing
b^2/a much larger than the local discretization length gives a genuinely blunt stress concentrator.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from scipy.spatial import Delaunay

from .mesh import (
    TriMesh,
    BoundaryData,
    _estimate_hbar,
    _estimate_hbar_tip,
    _precompute_element_data,
)


@dataclass
class BluntNotchGeometry:
    Lx: float = 2.0e-3
    Ly: float = 4.0e-3
    depth_a: float = 0.15e-3
    half_height_b: float = 0.30e-3

    @property
    def root_xy(self):
        return (float(self.depth_a), 0.0)

    @property
    def root_radius(self):
        return float(self.half_height_b**2 / max(self.depth_a, 1e-30))


def _inside_void(points: np.ndarray, geom: BluntNotchGeometry, margin: float = 0.0) -> np.ndarray:
    x = points[:, 0]
    y = points[:, 1]
    a = max(geom.depth_a + margin, 1e-30)
    b = max(geom.half_height_b + margin, 1e-30)
    return (x >= -1e-15) & ((x / a) ** 2 + (y / b) ** 2 < 1.0 - 1e-10)


def _triangle_crosses_void(nodes: np.ndarray, elems: np.ndarray, geom: BluntNotchGeometry) -> np.ndarray:
    """Conservative triangle filter for the non-convex edge-notch domain."""
    pts = nodes[elems]  # ne,3,2
    samples = [pts.mean(axis=1)]
    # sample each edge at 1/4, 1/2 and 3/4 to catch chords spanning the void
    for i, j in ((0, 1), (1, 2), (2, 0)):
        pi, pj = pts[:, i, :], pts[:, j, :]
        for t in (0.25, 0.5, 0.75):
            samples.append((1.0 - t) * pi + t * pj)
    bad = np.zeros(len(elems), dtype=bool)
    for s in samples:
        bad |= _inside_void(s, geom)
    return bad


def make_blunt_edge_notch_mesh(
    geom: BluntNotchGeometry,
    nx: int = 60,
    ny: int = 120,
    jitter: float = 0.10,
    root_h_fine: float = 20e-6,
    seed: int = 42,
):
    """Create a Delaunay triangular mesh with a smooth half-elliptical edge notch.

    Returns
    -------
    mesh, boundary_data, root_xy
    """
    rng = np.random.default_rng(seed)
    xv = np.linspace(0.0, geom.Lx, nx + 1)
    yv = np.linspace(-geom.Ly / 2.0, geom.Ly / 2.0, ny + 1)
    hx = xv[1] - xv[0]
    hy = yv[1] - yv[0]

    pts = []
    for j, y in enumerate(yv):
        for i, x in enumerate(xv):
            xx, yy = x, y
            if 0 < i < nx and 0 < j < ny:
                xx += jitter * hx * rng.uniform(-1.0, 1.0)
                yy += jitter * hy * rng.uniform(-1.0, 1.0)
            pts.append((xx, yy))
    pts = np.asarray(pts, float)
    pts = pts[~_inside_void(pts, geom)]

    # Explicitly resolve the free notch arc.
    arc_len_scale = max(root_h_fine, min(hx, hy) * 0.35)
    n_arc = max(40, int(np.pi * max(geom.depth_a, geom.half_height_b) / arc_len_scale))
    theta = np.linspace(-0.5 * np.pi, 0.5 * np.pi, n_arc)
    arc = np.column_stack([
        geom.depth_a * np.cos(theta),
        geom.half_height_b * np.sin(theta),
    ])

    # Local isotropic rings around the root to resolve initiation.
    root = np.asarray(geom.root_xy)
    rings = []
    r = max(root_h_fine, 0.25 * min(hx, hy))
    rmax = min(0.35e-3, 0.25 * geom.Ly)
    while r <= rmax:
        nr = max(12, int(2 * np.pi * r / max(root_h_fine, 1e-30)))
        th = np.linspace(0, 2 * np.pi, nr, endpoint=False)
        q = root[None, :] + r * np.column_stack([np.cos(th), np.sin(th)])
        inside_domain = (
            (q[:, 0] >= 0.0) & (q[:, 0] <= geom.Lx) &
            (q[:, 1] >= -geom.Ly / 2.0) & (q[:, 1] <= geom.Ly / 2.0)
        )
        q = q[inside_domain]
        q = q[~_inside_void(q, geom)]
        rings.append(q)
        r *= 1.45

    pts = np.vstack([pts, arc] + rings)
    # De-duplicate at a tolerance tied to the fine mesh.
    tol = max(root_h_fine * 0.05, 1e-12)
    key = np.round(pts / tol).astype(np.int64)
    _, keep = np.unique(key, axis=0, return_index=True)
    nodes = pts[np.sort(keep)]

    tri = Delaunay(nodes)
    elems = tri.simplices.copy()
    cent = nodes[elems].mean(axis=1)
    in_box = (
        (cent[:, 0] >= 0.0) & (cent[:, 0] <= geom.Lx) &
        (cent[:, 1] >= -geom.Ly / 2.0) & (cent[:, 1] <= geom.Ly / 2.0)
    )
    elems = elems[in_box]
    bad = _triangle_crosses_void(nodes, elems, geom)
    elems = elems[~bad]

    # Remove unused nodes and remap connectivity.
    used = np.unique(elems.ravel())
    remap = -np.ones(len(nodes), dtype=int)
    remap[used] = np.arange(len(used))
    nodes = nodes[used]
    elems = remap[elems]

    nn = len(nodes)
    ne = len(elems)
    area_e, dNdx_e, B_e = _precompute_element_data(nodes, elems)
    hbar = _estimate_hbar(nodes, elems)
    hbar_tip = _estimate_hbar_tip(nodes, elems, geom.depth_a, 0.0)
    mesh = TriMesh(
        nodes=nodes,
        elems=elems,
        nn=nn,
        ne=ne,
        ndof=2 * nn,
        hbar=hbar,
        area_e=area_e,
        dNdx_e=dNdx_e,
        B_e=B_e,
        hbar_tip=hbar_tip,
    )

    x, y = nodes[:, 0], nodes[:, 1]
    tol_b = max(0.75 * hbar, 1.5 * root_h_fine)
    top = np.where(np.abs(y - geom.Ly / 2.0) < tol_b)[0]
    bot = np.where(np.abs(y + geom.Ly / 2.0) < tol_b)[0]
    left_bot = int(np.argmin((x - 0.0) ** 2 + (y + geom.Ly / 2.0) ** 2))
    right_bot = int(np.argmin((x - geom.Lx) ** 2 + (y + geom.Ly / 2.0) ** 2))
    bnd = BoundaryData(
        top_nodes=top,
        bot_nodes=bot,
        left_bot=left_bot,
        right_bot=right_bot,
        notch_nodes=np.array([], dtype=int),  # no initial crack/damage field
    )
    return mesh, bnd, geom.root_xy

# -----------------------------------------------------------------------------
# Evolving blunt-feature geometry helpers used by the intact-FEM initiation workflow
# -----------------------------------------------------------------------------


def identify_feature_surface_nodes(mesh: TriMesh, geom: BluntNotchGeometry, n_sample: int = 400) -> np.ndarray:
    """Identify the explicit analytical ellipse-arc nodes, ordered by y.

    The mesh generator inserts the arc coordinates exactly.  Selecting by the
    ellipse residual is more reliable than nearest-neighbour sampling because a
    locally refined background ring can lie closer to a sample point than the
    corresponding arc node.
    """
    x = mesh.nodes[:, 0]
    y = mesh.nodes[:, 1]
    a = max(float(geom.depth_a), 1e-30)
    b = max(float(geom.half_height_b), 1e-30)
    q = (x/a)**2 + (y/b)**2
    mask = (x >= -1e-12) & (x <= 1.001*a) & (np.abs(q - 1.0) < 1e-7)
    idx = np.where(mask)[0]
    if len(idx) < 5:
        # Conservative fallback for aggressively rounded coordinate I/O.
        mask = (x >= -1e-12) & (x <= 1.01*a) & (np.abs(q - 1.0) < 5e-4)
        idx = np.where(mask)[0]
    return idx[np.argsort(mesh.nodes[idx, 1])]


def feature_tangent_normal(mesh: TriMesh, feature_nodes: np.ndarray):
    """Return unit tangent and outward-to-void unit normal on the feature arc.

    The feature node list must be ordered along the free surface.  The normal is
    oriented toward the removed elliptical void; near the root it therefore
    points approximately in the -x direction.
    """
    idx = np.asarray(feature_nodes, dtype=int)
    q = mesh.nodes[idx]
    npt = len(idx)
    t = np.zeros((npt, 2), dtype=float)
    if npt < 2:
        return t, t.copy()
    t[0] = q[1] - q[0]
    t[-1] = q[-1] - q[-2]
    if npt > 2:
        t[1:-1] = q[2:] - q[:-2]
    L = np.linalg.norm(t, axis=1)
    t /= np.maximum(L[:, None], 1e-30)
    n = np.column_stack([-t[:, 1], t[:, 0]])
    # Choose the sign that points toward the void center (0,0).
    toward_void = -q
    flip = np.sum(n * toward_void, axis=1) < 0.0
    n[flip] *= -1.0
    return t, n


def rebuild_mesh_geometry(mesh: TriMesh, root_xy=None) -> TriMesh:
    """Recompute element geometry after an ALE node-coordinate update.

    Connectivity and node numbering are preserved, so material state arrays
    remain attached to the same elements and no field interpolation is needed.
    """
    area_e, dNdx_e, B_e = _precompute_element_data(mesh.nodes, mesh.elems)
    mesh.area_e = area_e
    mesh.dNdx_e = dNdx_e
    mesh.B_e = B_e
    mesh.hbar = _estimate_hbar(mesh.nodes, mesh.elems)
    if root_xy is None:
        root_xy = mesh.nodes[np.argmax(mesh.nodes[:, 0])]
    mesh.hbar_tip = _estimate_hbar_tip(
        mesh.nodes, mesh.elems, float(root_xy[0]), float(root_xy[1])
    )
    return mesh


def local_root_xy(mesh: TriMesh, feature_nodes: np.ndarray) -> tuple[float, float]:
    """Current geometric root, taken as the maximum-x feature node."""
    idx = np.asarray(feature_nodes, dtype=int)
    j = int(idx[np.argmax(mesh.nodes[idx, 0])])
    return float(mesh.nodes[j, 0]), float(mesh.nodes[j, 1])


def local_root_radius(mesh: TriMesh, feature_nodes: np.ndarray) -> float:
    """Estimate local curvature radius at the root by a three-point circle fit."""
    idx = np.asarray(feature_nodes, dtype=int)
    q = mesh.nodes[idx]
    ir = int(np.argmax(q[:, 0]))
    if ir <= 0 or ir >= len(q) - 1:
        return float('nan')
    p1, p2, p3 = q[ir - 1], q[ir], q[ir + 1]
    a = np.linalg.norm(p2 - p1)
    b = np.linalg.norm(p3 - p2)
    c = np.linalg.norm(p3 - p1)
    v21 = p2 - p1
    v31 = p3 - p1
    # 2-D scalar cross product: twice the signed triangle area.
    area2 = abs(v21[0] * v31[1] - v21[1] * v31[0])
    if area2 <= 1e-30:
        return float('inf')
    # Circumradius R = abc/(4A) = abc/(2*area2), where area2=2A.
    return float(a * b * c / (2.0 * area2))


def apply_local_ale_surface_update(
    mesh: TriMesh,
    feature_nodes: np.ndarray,
    feature_normal_displacement: np.ndarray,
    decay_length: float,
    fixed_nodes: np.ndarray | None = None,
    max_move: float | None = None,
    min_area_fraction: float = 0.15,
):
    """Move the free feature surface and smoothly propagate motion into the mesh.

    Parameters
    ----------
    feature_normal_displacement
        Signed scalar normal displacement at each ordered feature node.
        Positive displacement is toward the void; negative is into the solid.
    decay_length
        Exponential ALE smoothing length used to propagate boundary motion into
        nearby interior nodes.

    Notes
    -----
    This is an updated-mesh kinematic step, not a remeshing/topology change.
    A backtracking line search prevents inverted or severely collapsed elements.
    """
    from scipy.spatial import cKDTree

    idx = np.asarray(feature_nodes, dtype=int)
    dh = np.asarray(feature_normal_displacement, dtype=float)
    if len(idx) == 0 or np.max(np.abs(dh)) <= 0.0:
        return 0.0
    _, normal = feature_tangent_normal(mesh, idx)
    disp_feature = normal * dh[:, None]

    if max_move is not None and np.isfinite(max_move) and max_move > 0:
        m = float(np.max(np.linalg.norm(disp_feature, axis=1)))
        if m > max_move:
            disp_feature *= max_move / max(m, 1e-30)

    tree = cKDTree(mesh.nodes[idx])
    dist, nearest = tree.query(mesh.nodes, k=1)
    w = np.exp(-dist / max(float(decay_length), 1e-30))
    disp = disp_feature[nearest] * w[:, None]
    disp[idx] = disp_feature
    if fixed_nodes is not None:
        disp[np.asarray(fixed_nodes, dtype=int)] = 0.0

    old_nodes = mesh.nodes.copy()
    old_area = mesh.area_e.copy()
    scale = 1.0
    accepted = False
    for _ in range(12):
        trial = old_nodes + scale * disp
        try:
            area_trial, _, _ = _precompute_element_data(trial, mesh.elems)
            if np.all(np.isfinite(area_trial)) and np.min(area_trial) > min_area_fraction * np.min(old_area):
                mesh.nodes[:] = trial
                accepted = True
                break
        except Exception:
            pass
        scale *= 0.5
    if not accepted:
        mesh.nodes[:] = old_nodes
        return 0.0
    return float(scale)
