#!/usr/bin/env python3
"""Generate a common elastic projected-crack K/F calibration for DBTT/1000 K.

This is a read-only mechanics reference, not a fracture simulation.  It uses the
same specimen dimensions, theta=0 cubic W elasticity, mesh scale and symmetric
fixed-grip boundary conditions as the authoritative FEM/CZM parity run, but it
removes all tip MPZ, stochastic kinetics, bulk plasticity and evolved crack-tip
morphology.  The crack is represented as the production initial notch plus a
straight thin post-notch wake along y=0.

The output K/F(a_projected) is intended ONLY as a common experimental-style
mapping for comparing PF and FEM/CZM remote loads:

    K_app(a) = F_actual(a) * [K/F]_common_reference(a)

It is deliberately independent of either solver's live J-to-K mapping.
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import numpy as np

from arrhenius_fracture.config import make_emergent_config
from arrhenius_fracture.mesh import make_tri_mesh, make_boundary_data
from arrhenius_fracture.fem import assemble_mechanics, solve_dirichlet
from arrhenius_fracture.j_integral import compute_J_integral
from arrhenius_fracture.crystal import cubic_plane_strain_D


C11 = 523.0e9
C12 = 203.0e9
C44 = 160.0e9
THETA_DEG = 0.0
A0_M = 0.5e-3
LX_M = 2.0e-3
LY_M = 4.0e-3
TIP_H_M = 1.0e-6
TIP_RATIO = 1.20
R_J_ELL_M = 30.0e-6  # production cluster outer radius = 8*ell = 240 um
RHO0_M2 = 5.0e12

# Dense enough to interpolate the actual PF/FEM event histories without using a
# handbook geometry correction, but still a tiny mechanics-only campaign.
CRACK_LENGTHS_M = np.asarray([
    0.500, 0.525, 0.550, 0.600, 0.650, 0.700, 0.800, 0.900,
    1.000, 1.100, 1.200, 1.300, 1.400, 1.500, 1.525,
], dtype=float) * 1e-3


def straight_reference_damage(mesh, geom, a_m: float) -> np.ndarray:
    """Initial notch + straight thin post-notch wake, no evolved tip geometry."""
    x = mesh.nodes[:, 0]
    y = mesh.nodes[:, 1]
    d = np.zeros(mesh.nn, dtype=float)
    d[(x <= geom.a0) & (np.abs(y) <= geom.notch_half_thickness)] = 1.0
    if a_m <= geom.a0 + 1.0e-15:
        return d

    c = mesh.nodes[mesh.elems].mean(axis=1)
    p0 = np.asarray([geom.a0, 0.0])
    p1 = np.asarray([a_m, 0.0])
    seg = p1-p0
    L2 = float(seg@seg) + 1.0e-30
    tt = np.clip(((c[:,0]-p0[0])*seg[0] + (c[:,1]-p0[1])*seg[1])/L2, 0.0, 1.0)
    proj = p0[None,:] + tt[:,None]*seg[None,:]
    dist2 = np.sum((c-proj)**2, axis=1)
    erad = np.sqrt(np.maximum(mesh.area_e, 1.0e-30))
    kill_r = max(mesh.hbar_tip, 0.5e-6)
    rad = np.maximum(kill_r, 0.7*erad)
    mask = (dist2 <= rad**2) & (c[:,0] >= geom.a0-2.0*kill_r) & (c[:,0] <= a_m+2.0*kill_r)
    d[mesh.elems[mask].ravel()] = 1.0
    return d


def solve_reference(a_m: float, Uapp_m: float) -> dict:
    cfg = make_emergent_config()
    cfg.geometry.Lx = LX_M
    cfg.geometry.Ly = LY_M
    cfg.geometry.a0 = A0_M
    cfg.mesh.nx = 36
    cfg.mesh.ny = 72
    cfg.mesh.tip_h_fine = TIP_H_M
    cfg.mesh.tip_ratio = TIP_RATIO
    mat = cfg.material

    mesh = make_tri_mesh(cfg.geometry, cfg.mesh, seed=42, tip_center=(a_m, 0.0))
    bnd = make_boundary_data(mesh, cfg.geometry)
    d = straight_reference_damage(mesh, cfg.geometry, a_m)
    D = cubic_plane_strain_D(C11, C12, C44, THETA_DEG)

    u = np.zeros(mesh.ndof)
    ep_gp = np.zeros((3, mesh.ne))
    rho_gp = np.full(mesh.ne, RHO0_M2)
    Kmat, Rint, sigma_gp, seq_gp, s1_gp, psi_gp = assemble_mechanics(
        mesh, u, ep_gp, rho_gp, d, D, mat, cohesive_network=None
    )
    u, Ftop = solve_dirichlet(Kmat, Rint, u, bnd, 0.5*Uapp_m, -0.5*Uapp_m)
    Kmat, Rint, sigma_gp, seq_gp, s1_gp, psi_gp = assemble_mechanics(
        mesh, u, ep_gp, rho_gp, d, D, mat, cohesive_network=None
    )
    J, KJ, Jinfo = compute_J_integral(
        mesh, u, sigma_gp, psi_gp, d,
        np.asarray([a_m, 0.0]), np.asarray([1.0, 0.0]), mat,
        ell=max(R_J_ELL_M, 3.0*mesh.hbar_tip),
        crack_segments=[(np.asarray([0.0,0.0]), np.asarray([a_m,0.0]))],
    )
    return {
        "projected_crack_length_m": float(a_m),
        "projected_extension_m": float(a_m-A0_M),
        "Uapp_m": float(Uapp_m),
        "Ftop_N": float(Ftop),
        "Fabs_N": abs(float(Ftop)),
        "J_J_m2": float(J),
        "K_reference_Pa_sqrtm": float(KJ),
        "K_over_F_Pa_sqrtm_per_N": float(KJ)/max(abs(float(Ftop)),1e-300),
        "mesh_hbar_tip_m": float(mesh.hbar_tip),
        "mesh_nn": int(mesh.nn),
        "mesh_ne": int(mesh.ne),
    }


def main() -> None:
    out = Path("dbtt1000_common_reference")
    out.mkdir(exist_ok=True)
    rows = []
    linearity = []
    for a in CRACK_LENGTHS_M:
        r1 = solve_reference(float(a), 1.0e-6)
        r2 = solve_reference(float(a), 0.5e-6)
        ratio1 = r1["K_over_F_Pa_sqrtm_per_N"]
        ratio2 = r2["K_over_F_Pa_sqrtm_per_N"]
        rel = abs(ratio1-ratio2)/max(abs(ratio1),abs(ratio2),1e-300)
        r1["linearity_rel_difference"] = rel
        rows.append(r1)
        linearity.append(rel)
        print(f"a={a*1e6:8.3f} um  K/F={ratio1:.9e}  linearity={rel:.3e}")

    csv_path = out/"dbtt1000_common_reference_K_over_F.csv"
    with csv_path.open("w", newline="") as fp:
        w = csv.DictWriter(fp, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    manifest = {
        "schema": "dbtt1000_common_projected_crack_reference_v1",
        "diagnostic_only": True,
        "physics_state": "linear_elastic_no_MPZ_no_plasticity_no_hazard",
        "crack_representation": "production_initial_notch_plus_straight_thin_post_notch_wake",
        "geometry_m": {"Lx":LX_M,"Ly":LY_M,"a0":A0_M},
        "elasticity": {"kind":"cubic_plane_strain","C11_Pa":C11,"C12_Pa":C12,"C44_Pa":C44,"theta_deg":THETA_DEG},
        "mesh": {"nx":36,"ny":72,"tip_h_fine_m":TIP_H_M,"tip_ratio":TIP_RATIO,"seed":42},
        "J_ell_m": R_J_ELL_M,
        "maximum_linearity_rel_difference": max(linearity),
        "formula": "K_app(a)=F_actual(a)*interpolated[K_reference/F_reference](a_projected)",
        "caveat": "Common projected-crack reference, not a claim of ASTM specimen validity; it intentionally removes either solver's evolved microscopic tip/wake state.",
    }
    (out/"manifest.json").write_text(json.dumps(manifest, indent=2))
    if max(linearity) > 2e-3:
        raise RuntimeError(f"reference K/F failed linearity check: max rel difference {max(linearity):.3e}")


if __name__ == "__main__":
    main()
