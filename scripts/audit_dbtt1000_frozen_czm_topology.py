#!/usr/bin/env python3
"""Frozen elastic audit of the DBTT/1000 K adaptive-CZM crack topology.

No hazard, RNG, MPZ evolution, backstress, plasticity, or fracture kinetics are
executed.  The startup mesh and crack-topology algorithms are the production
v10.0.5.18.3.9 path-aware corridor and atomic CZM routines.

The historical manifest retains all 179 accepted physical event lengths but not
the transverse sign of each endpoint.  This audit therefore uses the exact event
length sequence on a straight theta=0 path.  It is a topology control, not an
exact historical-path replay.
"""
from __future__ import annotations

import csv, json, os
from pathlib import Path
import numpy as np

from arrhenius_fracture.config import make_emergent_config
from arrhenius_fracture.mesh import make_boundary_data
from arrhenius_fracture.fem import assemble_mechanics, solve_dirichlet
from arrhenius_fracture.crystal import cubic_plane_strain_D
from arrhenius_fracture.crack_backend import SharpWakeBackend
from arrhenius_fracture.physical_refinement_mesh_v100510 import (
    configure_physical_refinement_v100510,
    clear_physical_refinement_v100510,
    make_physical_refinement_mesh_v100510,
)
from arrhenius_fracture import path_aware_growth_envelope_v10051834 as envelope
from arrhenius_fracture.atomic_path_corridor_adaptive_quality_v10051839 import (
    AdaptiveQualityAtomicPathCorridorCZMBackendV10051839,
)
from arrhenius_fracture.atomic_path_corridor_local_scale_v10051839 import (
    install as install_local_scale_gate,
    restore as restore_local_scale_gate,
)
from arrhenius_fracture.atomic_path_corridor_czm_v10051839 import reset_audit

ROOT = Path(__file__).resolve().parents[1]
EVENT_LEDGER = ROOT / "audit_inputs/dbtt1000_tip_only_event_lengths_m.json"
OUT = Path("dbtt1000_frozen_czm_topology_audit")

NX, NY = 36, 72
TIP_H_FINE_M, TIP_RATIO = 1.0e-6, 1.20
A0_M = 0.5e-3
UAPP_M = 1.0e-6
RHO0_M2 = 5.0e12
C11, C12, C44 = 523.0e9, 203.0e9, 160.0e9
THETA_DEG = 0.0
REFINEMENT_RADIUS_M = 330.0e-6
LPZ_UM = 50.0
TARGET_UM = 1000.0
GUARD_UM = 10.0
CHECKPOINTS_UM = (0.0, 100.0, 300.0, 500.0, 750.0, 1000.0)


def build_production_corridor(cfg):
    env = {
        "ARRHENIUS_PREFINED_MODE_I_CORRIDOR": "1",
        "ARRHENIUS_COMMITTED_TARGET_EXTENSION_UM": str(TARGET_UM),
        "ARRHENIUS_CORRIDOR_GUARD_UM": str(GUARD_UM),
        "ARRHENIUS_PHYSICAL_DA_UM": "5.0",
        "ARRHENIUS_CORRIDOR_PROCESS_ZONE_UM": str(LPZ_UM),
        "ARRHENIUS_CRYSTAL_THETA_DEG": str(THETA_DEG),
        "ARRHENIUS_MIN_GLOBAL_FORWARD": "0.05",
        "ARRHENIUS_MAX_CORRIDOR_H_OVER_LPZ": "0.25",
        "ARRHENIUS_ENVELOPE_SUPPORT_H_OVER_LPZ": "0.245",
        "ARRHENIUS_ENVELOPE_SUPPORT_REFINEMENT_FACTOR": "0.9",
        "ARRHENIUS_ENVELOPE_SUPPORT_CANDIDATES": "4",
        "ARRHENIUS_MIN_INITIAL_TRIANGLE_QUALITY": "0.035",
        "ARRHENIUS_MIN_ACCEPTED_TRIANGLE_QUALITY": "0.035",
        "ARRHENIUS_MIN_ACCEPTED_CHILD_AREA_RATIO": "0.08",
    }
    saved_env = {k: os.environ.get(k) for k in env}
    for k, v in env.items(): os.environ[k] = v
    configure_physical_refinement_v100510(REFINEMENT_RADIUS_M)
    old_detector = envelope._is_production_physical_constructor
    had_original = hasattr(envelope.path_aware_long_growth_envelope_mesh, "_original")
    old_original = getattr(envelope.path_aware_long_growth_envelope_mesh, "_original", None)
    envelope._is_production_physical_constructor = lambda _fn: True
    envelope.path_aware_long_growth_envelope_mesh._original = make_physical_refinement_mesh_v100510
    try:
        mesh = envelope.path_aware_long_growth_envelope_mesh(
            cfg.geometry, cfg.mesh, seed=42, tip_center=None
        )
    finally:
        envelope._is_production_physical_constructor = old_detector
        if had_original:
            envelope.path_aware_long_growth_envelope_mesh._original = old_original
        else:
            delattr(envelope.path_aware_long_growth_envelope_mesh, "_original")
        clear_physical_refinement_v100510()
        for k, v in saved_env.items():
            if v is None: os.environ.pop(k, None)
            else: os.environ[k] = v
    return mesh


def initial_state():
    cfg = make_emergent_config()
    cfg.geometry.a0 = A0_M
    cfg.mesh.nx, cfg.mesh.ny = NX, NY
    cfg.mesh.tip_h_fine, cfg.mesh.tip_ratio = TIP_H_FINE_M, TIP_RATIO
    cfg.mesh.jitter = 0.0
    mesh = build_production_corridor(cfg)
    bnd = make_boundary_data(mesh, cfg.geometry)
    d = np.zeros(mesh.nn)
    x, y = mesh.nodes[:, 0], mesh.nodes[:, 1]
    d[(x <= A0_M) & (np.abs(y) <= cfg.geometry.notch_half_thickness)] = 1.0
    return cfg, mesh, bnd, d, np.zeros(mesh.ndof)


def frozen_solve(mesh, bnd, d, network, cfg):
    mat = cfg.material
    D = cubic_plane_strain_D(C11, C12, C44, THETA_DEG)
    ep = np.zeros((3, mesh.ne)); rho = np.full(mesh.ne, RHO0_M2)
    u = np.zeros(mesh.ndof)
    K, R, *_ = assemble_mechanics(mesh, u, ep, rho, d, D, mat, cohesive_network=network)
    u, F = solve_dirichlet(K, R, u, bnd, 0.5*UAPP_M, -0.5*UAPP_M)
    K, R, sigma, seq, s1, psi = assemble_mechanics(
        mesh, u, ep, rho, d, D, mat, cohesive_network=network
    )
    return {
        "abs_Ftop_N_per_m": abs(float(F)),
        "compliance_m2_per_N": UAPP_M/max(abs(float(F)), 1e-300),
        "stiffness_N_per_m2": abs(float(F))/UAPP_M,
        "bulk_elastic_energy_J_per_m": float(np.sum(psi*mesh.area_e)),
        "node_count": int(mesh.nn), "element_count": int(mesh.ne),
    }


def add(rows, name, ext, result, failed=0):
    rows.append({"representation": name, "extension_um": 1e6*ext,
                 "failed_interfaces": int(failed), **result})


def crossed_checkpoint(previous_um, current_um, pending):
    return [x for x in pending if previous_um < x <= current_um + 1e-9]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    ledger = json.loads(EVENT_LEDGER.read_text())
    lengths = np.asarray(ledger["event_lengths_m"], float)
    if len(lengths) != 179 or np.any(lengths <= 0):
        raise RuntimeError("invalid retained DBTT event-length ledger")

    cfg_a, mesh_a, bnd_a, d_a, u_a = initial_state()
    cfg_s, mesh_s, bnd_s, d_s, u_s = initial_state()
    startup = {
        "nodes": int(mesh_a.nn), "triangles": int(mesh_a.ne),
        "hbar_tip_m": float(mesh_a.hbar_tip),
        "support_spacing_m": float(getattr(mesh_a, "production_refinement_support_spacing_m", np.nan)),
    }
    # Fail closed against the retained production startup manifest.
    if startup["nodes"] != 20321 or startup["triangles"] != 40502:
        raise RuntimeError(f"startup corridor does not reproduce production accepted candidate: {startup}")

    reset_audit()
    saved_stitch = install_local_scale_gate()
    adaptive = AdaptiveQualityAtomicPathCorridorCZMBackendV10051839(
        geom=cfg_a.geometry, event_damage=1.0, max_angle_error_deg=35.0
    )
    sharp = SharpWakeBackend()
    rows = []
    for name, network in (("adaptive_czm_failed", adaptive.cohesive_network),
                          ("adaptive_split_no_network", None)):
        add(rows, name, 0.0, frozen_solve(mesh_a,bnd_a,d_a,network,cfg_a), 0)
    add(rows, "sharp_wake", 0.0, frozen_solve(mesh_s,bnd_s,d_s,None,cfg_s), 0)

    p_a = np.array([A0_M,0.0]); p_s=p_a.copy(); ext_a=ext_s=0.0
    pending=list(CHECKPOINTS_UM[1:]); adaptive_fail=None
    try:
        for iev, length in enumerate(lengths, start=1):
            prev_um=ext_a*1e6
            p1=p_a+np.array([float(length),0.0])
            rr=adaptive.advance(mesh=mesh_a,boundary=bnd_a,damage=d_a,displacement=u_a,
                p0=p_a,p1=p1,direction=np.array([1.0,0.0]),front_id=0)
            if not rr.inserted:
                adaptive_fail={"event":iev,"developed_before_um":prev_um,
                    "requested_length_um":float(length*1e6),"reason":rr.reason}
                break
            mesh_a,bnd_a,d_a,u_a=rr.mesh,rr.boundary,rr.damage,rr.displacement
            p_a=p1; ext_a=float(p_a[0]-A0_M)

            # Sharp-wake control uses the same physical straight segment length.
            p1s=p_s+np.array([float(length),0.0])
            rs=sharp.advance(mesh=mesh_s,boundary=bnd_s,damage=d_s,displacement=u_s,
                p0=p_s,p1=p1s,direction=np.array([1.0,0.0]),front_id=0,
                kill_r=max(mesh_s.hbar_tip,0.5e-6))
            if not rs.inserted: raise RuntimeError(f"sharp wake failed event {iev}: {rs.reason}")
            mesh_s,bnd_s,d_s,u_s=rs.mesh,rs.boundary,rs.damage,rs.displacement
            p_s=p1s; ext_s=float(p_s[0]-A0_M)

            crossed=crossed_checkpoint(prev_um,ext_a*1e6,pending)
            for cp in crossed:
                add(rows,"adaptive_czm_failed",ext_a,
                    frozen_solve(mesh_a,bnd_a,d_a,adaptive.cohesive_network,cfg_a),
                    adaptive.cohesive_network.failed_count())
                add(rows,"adaptive_split_no_network",ext_a,
                    frozen_solve(mesh_a,bnd_a,d_a,None,cfg_a),
                    adaptive.cohesive_network.failed_count())
                add(rows,"sharp_wake",ext_s,frozen_solve(mesh_s,bnd_s,d_s,None,cfg_s),0)
                pending.remove(cp)
        # Always record final attained state even if it did not exactly cross 1000.
        if adaptive_fail is None and (not rows or abs(rows[-3]["extension_um"]-ext_a*1e6)>1e-6):
            add(rows,"adaptive_czm_failed",ext_a,frozen_solve(mesh_a,bnd_a,d_a,adaptive.cohesive_network,cfg_a),adaptive.cohesive_network.failed_count())
            add(rows,"adaptive_split_no_network",ext_a,frozen_solve(mesh_a,bnd_a,d_a,None,cfg_a),adaptive.cohesive_network.failed_count())
            add(rows,"sharp_wake",ext_s,frozen_solve(mesh_s,bnd_s,d_s,None,cfg_s),0)
    finally:
        restore_local_scale_gate(saved_stitch)

    initial={}
    for r in rows: initial.setdefault(r["representation"],r["compliance_m2_per_N"])
    for r in rows:
        r["normalized_compliance"]=r["compliance_m2_per_N"]/initial[r["representation"]]
        r["normalized_stiffness"]=1.0/r["normalized_compliance"]
    with (OUT/"frozen_topology_compliance.csv").open("w",newline="") as fp:
        w=csv.DictWriter(fp,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

    # Failed cohesive interfaces must be mechanically identical to omitting the network.
    A=[r for r in rows if r["representation"]=="adaptive_czm_failed"]
    N=[r for r in rows if r["representation"]=="adaptive_split_no_network"]
    identity=[]
    for a,n in zip(A,N):
        rel=abs(a["abs_Ftop_N_per_m"]-n["abs_Ftop_N_per_m"])/max(a["abs_Ftop_N_per_m"],n["abs_Ftop_N_per_m"],1e-300)
        identity.append({"extension_um":a["extension_um"],"relative_reaction_difference":rel})
    with (OUT/"failed_network_zero_stiffness_identity.csv").open("w",newline="") as fp:
        w=csv.DictWriter(fp,fieldnames=list(identity[0])); w.writeheader(); w.writerows(identity)
    max_identity=max(x["relative_reaction_difference"] for x in identity)

    summary={
        "schema":"dbtt1000_frozen_czm_topology_audit_v2",
        "diagnostic_only":True,"stochastic_simulation_rerun":False,
        "hazard_active":False,"MPZ_active":False,"plasticity_active":False,"backstress_active":False,
        "startup_corridor":startup,"event_length_source":ledger["source_run"],
        "exact_event_lengths_used":True,"straight_path_control":True,
        "historical_transverse_signs_available":False,
        "adaptive_failure":adaptive_fail,
        "failed_network_zero_stiffness_identity_pass":bool(max_identity<1e-10),
        "max_failed_network_vs_no_network_relative_reaction_difference":max_identity,
        "rows":rows,
    }
    (OUT/"summary.json").write_text(json.dumps(summary,indent=2))
    try:
        import matplotlib.pyplot as plt
        fig,ax=plt.subplots(figsize=(8.5,5.8))
        for name in ("adaptive_czm_failed","adaptive_split_no_network","sharp_wake"):
            rr=[r for r in rows if r["representation"]==name]
            ax.plot([r["extension_um"] for r in rr],[r["normalized_compliance"] for r in rr],marker="o",label=name)
        ax.set_xlabel("Developed straight-path extension (um)"); ax.set_ylabel("Compliance / initial compliance")
        ax.set_title("Frozen elastic DBTT crack-topology audit"); ax.grid(alpha=.25); ax.legend(); fig.tight_layout()
        fig.savefig(OUT/"frozen_topology_compliance.png",dpi=180); plt.close(fig)
    except Exception as exc: (OUT/"plot_error.txt").write_text(repr(exc))
    print(json.dumps({"startup":startup,"adaptive_failure":adaptive_fail,
        "identity_max_rel_diff":max_identity,"last_adaptive":A[-1],
        "last_sharp":[r for r in rows if r["representation"]=="sharp_wake"][-1]},indent=2))

if __name__=="__main__": main()
