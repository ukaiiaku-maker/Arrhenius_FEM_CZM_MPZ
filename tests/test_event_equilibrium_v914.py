from __future__ import annotations

import numpy as np

from arrhenius_fracture.config import make_emergent_config
from arrhenius_fracture.event_equilibrium_v914 import EventEquilibriumContext
from arrhenius_fracture.event_remesh_czm_v914 import EventRemeshCZMBackend
from arrhenius_fracture.fem import assemble_mechanics, plane_strain_D, solve_dirichlet
from arrhenius_fracture.mesh import make_boundary_data, make_tri_mesh


def test_post_event_equilibrium_preserves_load_and_conservative_state():
    cfg = make_emergent_config()
    cfg.mesh.nx = 12
    cfg.mesh.ny = 24
    cfg.mesh.tip_h_fine = 0.0
    mesh = make_tri_mesh(cfg.geometry, cfg.mesh, seed=42)
    bnd = make_boundary_data(mesh, cfg.geometry)
    D = plane_strain_D(cfg.material)
    ep = np.zeros((3, mesh.ne))
    rho = np.full(mesh.ne, 5.0e12)
    d = np.zeros(mesh.nn)
    u0 = np.zeros(mesh.ndof)
    K0, R0, *_ = assemble_mechanics(mesh, u0, ep, rho, d, D, cfg.material)
    ueq, _ = solve_dirichlet(K0, R0, u0, bnd, 1.0e-5, -1.0e-5)

    ctx = EventEquilibriumContext(
        original_assemble=assemble_mechanics,
        solve_callback=solve_dirichlet,
    )
    ctx.record_mechanics(mesh, ueq, ep, rho, d, D, cfg.material, 1.0e-6, None)

    def fake_j(mesh_, u_, sigma_, psi_, damage_, tip_, direction_, mat_, ell_,
               cfg=None, crack_segments=None, exclude_radius=0.0):
        return 100.0, 2.0e7, {"n_active_elements": int(mesh_.ne)}

    ctx.record_j_call(fake_j, cfg.material, 1.0e-5, None, 0.0, 90.0, 1.8e7)
    backend = EventRemeshCZMBackend(
        geom=cfg.geometry,
        target_h_m=max(mesh.hbar / 4.0, 1.0e-7),
        patch_radius_m=max(8.0 * mesh.hbar, 2.0e-5),
        max_edge_splits_per_event=8,
        target_edge_factor=1.25,
        min_remesh_triangle_quality=1.0e-4,
        penalty_normal_Pa_per_m=1.0e18,
        penalty_tangent_Pa_per_m=1.0e18,
    )
    tip = np.array([cfg.geometry.a0, 0.0])
    new_mesh, new_bnd, new_d, new_u, parent, _ = backend._refine_forward_patch(
        mesh, ueq, d, tip, np.array([1.0, 0.0])
    )
    u_after, record = ctx.equilibrate(
        pre_mesh=mesh,
        pre_boundary=bnd,
        pre_displacement=ueq,
        new_mesh=new_mesh,
        new_boundary=new_bnd,
        new_damage=new_d,
        new_displacement=new_u,
        parent_map=parent,
        cohesive_network=None,
        new_tip=tip,
        direction=np.array([1.0, 0.0]),
        crack_segments=[],
        event_index=0,
        front_id=0,
    )
    assert u_after.shape == (new_mesh.ndof,)
    assert record["physical_time_increment_s"] == 0.0
    assert record["hazard_action_increment"] == 0.0
    assert record["max_relative_boundary_displacement_drift"] < 1.0e-13
    assert record["relative_rho_area_integral_error"] < 1.0e-12
    assert record["max_relative_ep_area_integral_error"] < 1.0e-12
    assert record["relative_total_mesh_area_error"] < 1.0e-12
    assert record["J_after_event_status"] == "ok"
    assert record["KJ_after_event_equilibrium_Pa_sqrt_m"] == 2.0e7
    assert record["J_after_event_active_elements"] == new_mesh.ne
