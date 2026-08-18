#!/usr/bin/env python3
# CI trigger: frozen-mechanics audit; no physics changes.
"""Frozen mechanics audit of adaptive-CZM crack topology.

No Arrhenius hazard, RNG, MPZ state, backstress evolution, plasticity, or fracture
kinetics are executed.  The purpose is only to ask whether a sequence of fully
failed adaptive-CZM interfaces produces the expected loss of global stiffness as
projected crack length grows.

Three elastic representations are compared on the same specimen/control mesh:

1. adaptive_czm_failed: production AdaptiveCZMBackend, each inserted interface
   immediately fully failed (damage=1).
2. adaptive_split_no_network: identical split-node topology, but the cohesive
   network is omitted from assembly entirely.  It MUST match (1) to numerical
   tolerance because fully failed cohesive elements have zero tensile/tangential
   stiffness.
3. sharp_wake: production SharpWakeBackend stiffness-kill representation.

The historical DBTT run is theta=0 and single-front.  This controlled replay uses
straight 5 um physical advances to 1 mm.  It is therefore a topology/implementation
audit, not an exact historical-path replay; the archived run's maximum segment
angle was only ~7 deg, so this is an appropriate first isolation test.
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
from arrhenius_fracture.crack_backend import AdaptiveCZMBackend, SharpWakeBackend

NX = 36
NY = 72
TIP_H_FINE_M = 1.0e-6
TIP_RATIO = 1.20
A0_M = 0.5e-3
DA_M = 5.0e-6
TARGET_EXTENSION_M = 1.0e-3
UAPP_M = 1.0e-6
RHO0_M2 = 5.0e12
C11 = 523.0e9
C12 = 203.0e9
C44 = 160.0e9
THETA_DEG = 0.0
CHECKPOINT_UM = (0, 100, 300, 500, 750, 1000)
OUT = Path("dbtt1000_frozen_czm_topology_audit")


def initial_state():
    cfg = make_emergent_config()
    cfg.geometry.a0 = A0_M
    cfg.mesh.nx = NX
    cfg.mesh.ny = NY
    cfg.mesh.tip_h_fine = TIP_H_FINE_M
    cfg.mesh.tip_ratio = TIP_RATIO
    cfg.mesh.jitter = 0.0
    mesh = make_tri_mesh(cfg.geometry, cfg.mesh, seed=42)
    bnd = make_boundary_data(mesh, cfg.geometry)
    d = np.zeros(mesh.nn)
    x = mesh.nodes[:, 0]
    y = mesh.nodes[:, 1]
    d[(x <= cfg.geometry.a0) & (np.abs(y) <= cfg.geometry.notch_half_thickness)] = 1.0
    u = np.zeros(mesh.ndof)
    return cfg, mesh, bnd, d, u


def frozen_solve(mesh, bnd, d, cohesive_network, cfg):
    mat = cfg.material
    D = cubic_plane_strain_D(C11, C12, C44, THETA_DEG)
    ep = np.zeros((3, mesh.ne))
    rho = np.full(mesh.ne, RHO0_M2)
    u = np.zeros(mesh.ndof)
    Kmat, Rint, sigma, seq, s1, psi = assemble_mechanics(
        mesh, u, ep, rho, d, D, mat, cohesive_network=cohesive_network
    )
    u, Ftop = solve_dirichlet(Kmat, Rint, u, bnd, 0.5 * UAPP_M, -0.5 * UAPP_M)
    Kmat, Rint, sigma, seq, s1, psi = assemble_mechanics(
        mesh, u, ep, rho, d, D, mat, cohesive_network=cohesive_network
    )
    stored = float(np.sum(psi * mesh.area_e))
    return {
        "Ftop_N_per_m": float(Ftop),
        "abs_Ftop_N_per_m": abs(float(Ftop)),
        "Uapp_m": UAPP_M,
        "compliance_m2_per_N": UAPP_M / max(abs(float(Ftop)), 1e-300),
        "stiffness_N_per_m2": abs(float(Ftop)) / UAPP_M,
        "bulk_elastic_energy_J_per_m": stored,
        "node_count": int(mesh.nn),
        "element_count": int(mesh.ne),
    }


def append_row(rows, representation, extension_m, result, extra=None):
    row = {"representation": representation, "extension_um": extension_m * 1e6, **result}
    if extra:
        row.update(extra)
    rows.append(row)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cfg_a, mesh_a, bnd_a, d_a, u_a = initial_state()
    cfg_s, mesh_s, bnd_s, d_s, u_s = initial_state()
    adaptive = AdaptiveCZMBackend(geom=cfg_a.geometry, event_damage=1.0, max_angle_error_deg=35.0)
    sharp = SharpWakeBackend()
    rows = []
    append_row(rows, "adaptive_czm_failed", 0.0, frozen_solve(mesh_a, bnd_a, d_a, adaptive.cohesive_network, cfg_a), {"failed_interfaces": 0})
    append_row(rows, "adaptive_split_no_network", 0.0, frozen_solve(mesh_a, bnd_a, d_a, None, cfg_a), {"failed_interfaces": 0})
    append_row(rows, "sharp_wake", 0.0, frozen_solve(mesh_s, bnd_s, d_s, None, cfg_s), {"failed_interfaces": 0})
    target_set = set(CHECKPOINT_UM[1:])
    p0_a = np.array([A0_M, 0.0]); p0_s = np.array([A0_M, 0.0])
    ext_a = 0.0; ext_s = 0.0
    n_events = int(round(TARGET_EXTENSION_M / DA_M))
    adaptive_fail = None
    for iev in range(1, n_events + 1):
        p1_req_a = np.array([A0_M + iev * DA_M, 0.0])
        rr = adaptive.advance(mesh=mesh_a, boundary=bnd_a, damage=d_a, displacement=u_a,
            p0=p0_a, p1=p1_req_a, direction=np.array([1.0, 0.0]), front_id=0,
            kill_r=max(mesh_a.hbar_tip, 0.5e-6))
        if not rr.inserted:
            adaptive_fail = {"event": iev, "requested_extension_um": iev * DA_M * 1e6,
                             "reason": rr.reason, "angle_error_deg": rr.angle_error_deg}
            break
        mesh_a, bnd_a, d_a, u_a = rr.mesh, rr.boundary, rr.damage, rr.displacement
        p0_a = np.array([adaptive.advance_log[-1]["x1"], adaptive.advance_log[-1]["y1"]])
        ext_a = float(p0_a[0] - A0_M)
        p1_s = np.array([A0_M + iev * DA_M, 0.0])
        rs = sharp.advance(mesh=mesh_s, boundary=bnd_s, damage=d_s, displacement=u_s,
            p0=p0_s, p1=p1_s, direction=np.array([1.0, 0.0]), front_id=0,
            kill_r=max(mesh_s.hbar_tip, 0.5e-6))
        if not rs.inserted:
            raise RuntimeError(f"sharp-wake advance failed at event {iev}: {rs.reason}")
        mesh_s, bnd_s, d_s, u_s = rs.mesh, rs.boundary, rs.damage, rs.displacement
        p0_s = p1_s; ext_s = float(p0_s[0] - A0_M)
        requested_um = int(round(iev * DA_M * 1e6))
        if requested_um in target_set:
            append_row(rows, "adaptive_czm_failed", ext_a,
                       frozen_solve(mesh_a, bnd_a, d_a, adaptive.cohesive_network, cfg_a),
                       {"failed_interfaces": adaptive.cohesive_network.failed_count()})
            append_row(rows, "adaptive_split_no_network", ext_a,
                       frozen_solve(mesh_a, bnd_a, d_a, None, cfg_a),
                       {"failed_interfaces": adaptive.cohesive_network.failed_count()})
            append_row(rows, "sharp_wake", ext_s, frozen_solve(mesh_s, bnd_s, d_s, None, cfg_s),
                       {"failed_interfaces": 0})
    if adaptive_fail is not None:
        (OUT / "adaptive_failure.json").write_text(json.dumps(adaptive_fail, indent=2))
    initial = {}
    for row in rows:
        initial.setdefault(row["representation"], row["compliance_m2_per_N"])
    for row in rows:
        row["normalized_compliance"] = row["compliance_m2_per_N"] / initial[row["representation"]]
        row["normalized_stiffness"] = 1.0 / row["normalized_compliance"]
    with (OUT / "frozen_topology_compliance.csv").open("w", newline="") as fp:
        writer = csv.DictWriter(fp, fieldnames=list(rows[0].keys())); writer.writeheader(); writer.writerows(rows)
    by_key = {(r["representation"], r["extension_um"]): r for r in rows}
    identity = []
    for ext in sorted({r["extension_um"] for r in rows}):
        a = by_key.get(("adaptive_czm_failed", ext)); n = by_key.get(("adaptive_split_no_network", ext))
        if a is None or n is None: continue
        rel = abs(a["abs_Ftop_N_per_m"] - n["abs_Ftop_N_per_m"]) / max(a["abs_Ftop_N_per_m"], n["abs_Ftop_N_per_m"], 1e-300)
        identity.append({"extension_um": ext, "relative_reaction_difference": rel})
    with (OUT / "failed_network_zero_stiffness_identity.csv").open("w", newline="") as fp:
        writer = csv.DictWriter(fp, fieldnames=list(identity[0].keys())); writer.writeheader(); writer.writerows(identity)
    max_identity = max(x["relative_reaction_difference"] for x in identity)
    summary = {
        "schema": "dbtt1000_frozen_czm_topology_audit_v1", "diagnostic_only": True,
        "stochastic_simulation_rerun": False, "hazard_active": False, "MPZ_active": False,
        "plasticity_active": False, "backstress_active": False, "theta_deg": THETA_DEG,
        "straight_control_path": True, "physical_step_um": DA_M * 1e6,
        "target_extension_um": TARGET_EXTENSION_M * 1e6,
        "matched_mesh_controls": {"nx": NX, "ny": NY, "tip_h_fine_m": TIP_H_FINE_M, "tip_ratio": TIP_RATIO},
        "adaptive_failure": adaptive_fail,
        "max_failed_network_vs_no_network_relative_reaction_difference": max_identity,
        "failed_network_zero_stiffness_identity_pass": bool(max_identity < 1e-10), "rows": rows,
        "interpretation_contract": "Controlled topology replay, not exact historical path replay."
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2))
    try:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8.5, 5.8))
        for name in ("adaptive_czm_failed", "adaptive_split_no_network", "sharp_wake"):
            rr = [r for r in rows if r["representation"] == name]
            ax.plot([r["extension_um"] for r in rr], [r["normalized_compliance"] for r in rr], marker="o", label=name)
        ax.set_xlabel("Projected crack extension (um)"); ax.set_ylabel("Compliance / initial-notch compliance")
        ax.set_title("Frozen elastic crack-topology audit"); ax.grid(alpha=0.25); ax.legend(); fig.tight_layout()
        fig.savefig(OUT / "frozen_topology_compliance.png", dpi=180); plt.close(fig)
    except Exception as exc:
        (OUT / "plot_error.txt").write_text(repr(exc))
    print(json.dumps({"adaptive_failure": adaptive_fail, "identity_max_rel_diff": max_identity,
        "final_rows": [r for r in rows if abs(r["extension_um"] - max(x["extension_um"] for x in rows)) < 1e-9]}, indent=2))

if __name__ == "__main__":
    main()
