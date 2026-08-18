#!/usr/bin/env python3
"""Frozen-geometry DBTT/1000 K adaptive-CZM compliance audit.

No hazard, RNG, MPZ transport, plasticity, or fracture evolution is run.  The
accepted 179-event geometry ledger is replayed mechanically, with all inserted
cohesive segments fully failed (damage=1).  At selected extensions a small
linear-elastic opening solve measures global compliance and stored energy.

The production manifest retained event lengths and projected dx but not the
sign of each small transverse increment.  We therefore replay three bounded
geometries: alternating signs, all-positive signs, and straight projected
increments.  This tests whether the missing transverse signs can affect the
compliance conclusion.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path
import numpy as np

from arrhenius_fracture.config import make_emergent_config
from arrhenius_fracture.mesh import make_tri_mesh, make_boundary_data
from arrhenius_fracture.fem import assemble_mechanics, solve_dirichlet
from arrhenius_fracture.crystal import cubic_plane_strain_D
from arrhenius_fracture.crack_backend import AdaptiveCZMBackend

ROOT = Path(__file__).resolve().parents[1]
GEOM = ROOT / "reference_inputs/dbtt1000_frozen_audit/event_geometry_projected.json"
OUT = ROOT / "dbtt1000_frozen_czm_compliance_audit"
TARGETS_UM = (0.0, 100.0, 300.0, 500.0, 750.0, 1000.0)
UAPP = 1.0e-7
C11, C12, C44 = 523e9, 203e9, 160e9
RHO0 = 5e12


def vectors(data, mode):
    L = np.asarray(data["event_lengths_m"], float)
    dx = np.asarray(data["projected_dx_m"], float)
    dyabs = np.sqrt(np.maximum(L*L - dx*dx, 0.0))
    if mode == "alternating":
        signs = np.where(np.arange(len(L)) % 2 == 0, 1.0, -1.0)
        dy = signs * dyabs
    elif mode == "one_sided":
        dy = dyabs
    elif mode == "straight_projected":
        dy = np.zeros_like(dx)
    else:
        raise ValueError(mode)
    return np.column_stack([dx, dy])


def solve_elastic(mesh, bnd, damage, network, D, mat):
    ep = np.zeros((3, mesh.ne))
    rho = np.full(mesh.ne, RHO0)
    u = np.zeros(mesh.ndof)
    K, R, sig, seq, s1, psi = assemble_mechanics(
        mesh, u, ep, rho, damage, D, mat, cohesive_network=network
    )
    u, F = solve_dirichlet(K, R, u, bnd, 0.5*UAPP, -0.5*UAPP)
    K, R, sig, seq, s1, psi = assemble_mechanics(
        mesh, u, ep, rho, damage, D, mat, cohesive_network=network
    )
    C = UAPP / abs(float(F))
    energy = 0.5 * float(u @ (K @ u))
    return abs(float(F)), C, energy


def initial_state():
    cfg = make_emergent_config()
    cfg.geometry.Lx = 2e-3
    cfg.geometry.Ly = 4e-3
    cfg.geometry.a0 = 0.5e-3
    cfg.mesh.nx = 36
    cfg.mesh.ny = 72
    cfg.mesh.tip_h_fine = 1e-6
    cfg.mesh.tip_ratio = 1.20
    mesh = make_tri_mesh(cfg.geometry, cfg.mesh, seed=42)
    bnd = make_boundary_data(mesh, cfg.geometry)
    x, y = mesh.nodes[:,0], mesh.nodes[:,1]
    damage = np.zeros(mesh.nn)
    damage[(x <= cfg.geometry.a0) & (np.abs(y) <= cfg.geometry.notch_half_thickness)] = 1.0
    u = np.zeros(mesh.ndof)
    backend = AdaptiveCZMBackend(
        geom=cfg.geometry,
        penalty_normal_Pa_per_m=1e18,
        penalty_tangent_Pa_per_m=1e18,
        max_angle_error_deg=35.0,
        event_damage=1.0,
    )
    return cfg, mesh, bnd, damage, u, backend


def replay(mode, data):
    cfg, mesh, bnd, damage, u, backend = initial_state()
    mat = cfg.material
    D = cubic_plane_strain_D(C11, C12, C44, 0.0)
    vecs = vectors(data, mode)
    p = np.array([cfg.geometry.a0, 0.0], float)
    projected = 0.0
    rows = []
    next_target = 0

    def record(label):
        F, C, W = solve_elastic(mesh, bnd, damage, backend.cohesive_network, D, mat)
        rows.append({
            "path_mode": mode,
            "checkpoint_um": label,
            "actual_projected_extension_um": projected*1e6,
            "tip_x_m": float(p[0]),
            "tip_y_m": float(p[1]),
            "nodes": int(mesh.nn),
            "elements": int(mesh.ne),
            "cohesive_segments": int(len(backend.cohesive_network.elements)),
            "F_per_thickness_N_per_m": F,
            "Uapp_m": UAPP,
            "compliance_m2_per_N": C,
            "elastic_energy_J_per_m": W,
        })

    record(0.0)
    next_target = 1
    for i, v in enumerate(vecs):
        req = float(np.linalg.norm(v))
        if req <= 0:
            continue
        direction = v / req
        p1 = p + v
        rr = backend.advance(
            mesh=mesh, boundary=bnd, damage=damage, displacement=u,
            p0=p, p1=p1, direction=direction, front_id=0,
            kill_r=max(mesh.hbar_tip, 0.5e-6),
        )
        if not rr.inserted or rr.moved <= 0.0:
            raise RuntimeError(
                f"{mode}: frozen replay failed at event {i+1}: {rr.reason}; "
                f"requested={req:.9e} m"
            )
        mesh, bnd, damage, u = rr.mesh, rr.boundary, rr.damage, rr.displacement
        log = backend.advance_log[-1]
        p = np.array([log["x1"], log["y1"]], float)
        projected += float(v[0])
        while next_target < len(TARGETS_UM) and projected*1e6 >= TARGETS_UM[next_target]-1e-9:
            record(TARGETS_UM[next_target])
            next_target += 1
    if next_target < len(TARGETS_UM):
        record(TARGETS_UM[-1])
    return rows


def main():
    data = json.loads(GEOM.read_text())
    if int(data["n_events"]) != 179:
        raise RuntimeError("frozen event count is not 179")
    OUT.mkdir(exist_ok=True)
    all_rows = []
    for mode in ("alternating", "one_sided", "straight_projected"):
        all_rows.extend(replay(mode, data))
    with (OUT/"frozen_adaptive_czm_compliance.csv").open("w", newline="") as fp:
        w = csv.DictWriter(fp, fieldnames=list(all_rows[0].keys()))
        w.writeheader(); w.writerows(all_rows)
    # Normalize compliance within each path reconstruction for the decisive audit.
    summary = {}
    for mode in sorted({r["path_mode"] for r in all_rows}):
        rr = [r for r in all_rows if r["path_mode"] == mode]
        c0 = rr[0]["compliance_m2_per_N"]
        summary[mode] = [
            {
                "target_um": r["checkpoint_um"],
                "actual_um": r["actual_projected_extension_um"],
                "C_over_C0": r["compliance_m2_per_N"]/c0,
                "F_over_F0": rr[0]["F_per_thickness_N_per_m"] / r["F_per_thickness_N_per_m"],
                "tip_y_um": r["tip_y_m"]*1e6,
                "cohesive_segments": r["cohesive_segments"],
            } for r in rr
        ]
    (OUT/"summary.json").write_text(json.dumps({
        "schema":"dbtt1000_frozen_adaptive_czm_compliance_audit_v1",
        "no_fracture_simulation_rerun": True,
        "physics_disabled":["hazard","RNG","MPZ","plasticity","backstress evolution","shielding evolution"],
        "missing_original_transverse_signs": True,
        "path_reconstructions": list(summary),
        "summary": summary,
    }, indent=2))
    print(json.dumps(summary, indent=2))

if __name__ == "__main__":
    main()
