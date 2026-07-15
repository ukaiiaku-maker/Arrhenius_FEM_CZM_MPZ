from __future__ import annotations

import numpy as np

from arrhenius_fracture.config import make_emergent_config
from arrhenius_fracture.event_remesh_czm_v914 import EventRemeshCZMBackend
from arrhenius_fracture.mesh import make_tri_mesh


def test_forward_patch_refinement_has_exact_conservative_parent_map():
    cfg = make_emergent_config()
    cfg.mesh.nx = 12
    cfg.mesh.ny = 24
    cfg.mesh.tip_h_fine = 0.0
    mesh = make_tri_mesh(cfg.geometry, cfg.mesh, seed=42)
    backend = EventRemeshCZMBackend(
        geom=cfg.geometry,
        target_h_m=max(mesh.hbar / 4.0, 1.0e-7),
        patch_radius_m=max(8.0 * mesh.hbar, 2.0e-5),
        max_edge_splits_per_event=12,
        target_edge_factor=1.25,
        min_remesh_triangle_quality=1.0e-4,
        penalty_normal_Pa_per_m=1.0e18,
        penalty_tangent_Pa_per_m=1.0e18,
    )
    u = np.zeros(mesh.ndof)
    d = np.zeros(mesh.nn)
    tip = np.array([cfg.geometry.a0, 0.0])
    new_mesh, new_bnd, new_d, new_u, parent, audit = backend._refine_forward_patch(
        mesh, u, d, tip, np.array([1.0, 0.0])
    )
    assert new_mesh.ne >= mesh.ne
    assert parent.shape == (new_mesh.ne,)
    assert np.min(parent) >= 0
    assert np.max(parent) < mesh.ne
    assert new_u.shape == (new_mesh.ndof,)
    assert new_d.shape == (new_mesh.nn,)
    assert new_bnd.top_nodes.size > 0
    assert audit["n_edge_splits"] > 0

    rho = np.linspace(2.0e12, 8.0e12, mesh.ne)
    ep = np.vstack([
        np.linspace(0.0, 1.0e-4, mesh.ne),
        np.linspace(1.0e-4, 0.0, mesh.ne),
        np.linspace(-2.0e-5, 2.0e-5, mesh.ne),
    ])
    old_rho = float(np.sum(rho * mesh.area_e))
    new_rho = float(np.sum(rho[parent] * new_mesh.area_e))
    old_ep = np.sum(ep * mesh.area_e[None, :], axis=1)
    new_ep = np.sum(ep[:, parent] * new_mesh.area_e[None, :], axis=1)
    assert np.isclose(new_rho, old_rho, rtol=1.0e-12, atol=1.0e-20)
    assert np.allclose(new_ep, old_ep, rtol=1.0e-12, atol=1.0e-20)
    assert audit["max_parent_relative_area_conservation_error"] < 1.0e-12
    assert audit["relative_total_area_error"] < 1.0e-12
