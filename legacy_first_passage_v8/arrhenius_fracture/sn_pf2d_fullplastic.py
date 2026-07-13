"""2-D S-N initiation with fully Arrhenius plasticity and evolving geometry.

This driver starts from a blunt free-surface feature and d(x)=0 everywhere.
The cyclic plastic constitutive response is the same rate-space mechanism family
used by the fatigue framework:

    surface emission -> Peierls glide -> Taylor junction depinning.

Each step uses a scaled EXP-floor free-energy barrier derived from the surface
nucleation family.  Peierls and Taylor are sequential residence-time bottlenecks;
emission is in series with that mobility branch.  The Taylor density dependence
enters through phi_T(rho)=min[1/(2 b sqrt(rho)), phi_max].

There is no quasi-static Taylor stress, no additive Peierls stress, no athermal
Taylor floor, and no hard yield gate.  The completed Arrhenius event chain drives
actual FEM plastic eigenstrain.  rho evolves from the accepted Arrhenius plastic
strain, the retained eigenstrain generates residual stress on re-equilibration,
and surface-localized plastic strain drives ALE evolution of the blunt feature.
The resulting intrusion/extrusion-like relief can sharpen geometrically before
phase-field crack formation.

A representative physical cycle is integrated explicitly, then block-accelerated.
The accepted block size is limited by plastic strain, rho change, surface motion,
and crack-clock increments so the evolving state is re-solved frequently enough.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Iterable

import numpy as np
from scipy.special import gammainc
from scipy.sparse.linalg import spsolve

from .config import KB, EV_TO_J, make_emergent_config
from .fem import (
    plane_strain_D, assemble_mechanics, assemble_pf_matrices,
    project_gp_to_nodes, stress_state,
)
from .materials import FractureModel
from .sn_arrhenius_chain import build_chain_from_namespace, ArrheniusPlasticChain
from .phase_field import update_phase_field, at2_surface_energy
from .sn_geometry import (
    BluntNotchGeometry, make_blunt_edge_notch_mesh,
    identify_feature_surface_nodes, feature_tangent_normal,
    rebuild_mesh_geometry, local_root_xy, local_root_radius,
    apply_local_ale_surface_update,
)
from .sn_v1 import make_barriers, KBEV


def _write_csv(path: Path, rows):
    if not rows:
        return
    keys = sorted(set().union(*(r.keys() for r in rows)))
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader(); w.writerows(rows)


def _solve_symmetric_tension(K, Rint, u, bnd, mesh, Uy_top, Uy_bot):
    """Incremental symmetric-displacement solve with one x anchor."""
    ndof = mesh.ndof
    prescribed = np.zeros(ndof, dtype=bool)
    up = np.zeros(ndof)
    prescribed[2*bnd.top_nodes + 1] = True
    up[2*bnd.top_nodes + 1] = Uy_top
    prescribed[2*bnd.bot_nodes + 1] = True
    up[2*bnd.bot_nodes + 1] = Uy_bot
    x, y = mesh.nodes[:, 0], mesh.nodes[:, 1]
    anchor = int(np.argmin((x - x.max())**2 + y**2))
    prescribed[2*anchor] = True
    up[2*anchor] = 0.0
    free = ~prescribed
    Kc = K.tocsr()
    du_p = up[prescribed] - u[prescribed]
    rhs = -Rint[free] - Kc[np.ix_(free, prescribed)] @ du_p
    un = u.copy()
    un[free] = u[free] + spsolve(Kc[np.ix_(free, free)], rhs)
    un[prescribed] = up[prescribed]
    Rfull = Rint + Kc @ (un - u)
    Ftop = float(np.sum(Rfull[2*bnd.top_nodes + 1]))
    return un, Ftop


def _affine_stress_control_displacements(mesh, bnd, cfg, Dmat, ep_gp, rho_gp, d,
                                         sigma_max_Pa, sigma_min_Pa, u_guess):
    """Calibrate top displacement for target nominal max/min stresses.

    At fixed internal state the small-strain equilibrium problem is linear and
    reaction force is affine in imposed displacement because plastic strain acts
    as an eigenstrain.  Two solves therefore determine the exact affine map.
    """
    K, Rint, *_ = assemble_mechanics(mesh, u_guess, ep_gp, rho_gp, d, Dmat, cfg.material)
    uz, F0 = _solve_symmetric_tension(K, Rint, u_guess, bnd, mesh, 0.0, 0.0)
    K0, R0, *_ = assemble_mechanics(mesh, uz, ep_gp, rho_gp, d, Dmat, cfg.material)
    Uprobe = max(mesh.hbar_tip, 1e-8) * 1e-3
    up, Fp = _solve_symmetric_tension(K0, R0, uz, bnd, mesh, Uprobe, -Uprobe)
    slope = (Fp - F0) / max(Uprobe, 1e-30)
    if not np.isfinite(slope) or abs(slope) < 1e-30:
        raise RuntimeError("stress-control displacement calibration failed")
    area2d = max(float(mesh.nodes[:,0].max() - mesh.nodes[:,0].min()), 1e-30)
    Fmax = sigma_max_Pa * area2d
    Fmin = sigma_min_Pa * area2d
    Umax = (Fmax - F0) / slope
    Umin = (Fmin - F0) / slope
    return float(Umax), float(Umin), uz, float(F0), float(slope)


def _equivalent_dep_from_tensor(dep):
    """Equivalent plastic increment from in-plane Voigt engineering strain."""
    exx, eyy, gxy = dep[0], dep[1], dep[2]
    # Approximate plane-strain J2 norm with incompressible plastic ezz=-(exx+eyy).
    ezz = -(exx + eyy)
    tensor_sq = exx**2 + eyy**2 + ezz**2 + 0.5 * gxy**2
    return np.sqrt(np.maximum(2.0/3.0 * tensor_sq, 0.0))




def _von_mises_flow_direction(sigma_gp, nu):
    sx, sy, txy = sigma_gp[0], sigma_gp[1], sigma_gp[2]
    szz = nu * (sx + sy)
    mean = (sx + sy + szz) / 3.0
    sxx = sx - mean; syy = sy - mean; szzd = szz - mean
    seq = np.sqrt(np.maximum(1.5*(sxx*sxx + syy*syy + szzd*szzd + 2.0*txy*txy), 0.0))
    norm = np.sqrt(np.maximum(sxx*sxx + syy*syy + szzd*szzd + 2.0*txy*txy, 1e-60))
    return seq, sxx/norm, syy/norm, txy/norm


def _arrhenius_chain_plastic_update(
    ep_gp, rho_gp, sigma_gp, mat, T_K, dt_s,
    chain: ArrheniusPlasticChain,
    k_store: float, k_dyn: float,
    rho_floor: float, rho_cap: float,
    max_dep_phase: float, max_rho_rel_phase: float,
):
    """Fully Arrhenius local plasticity update for one resolved cycle phase.

    No quasi-static flow surface, additive Peierls stress, athermal Taylor
    floor, or alpha*G*b*sqrt(rho) subtraction is used.  The constitutive rate
    comes directly from the serial Arrhenius chain

        emission -> Peierls glide -> Taylor depinning.

    Taylor density dependence enters through the Arrhenius node amplification
    phi_T(rho)=min[1/(2 b sqrt(rho)), phi_max].  The only cap is a local
    elastic-energy relaxation-distance projection, which prevents a numerical
    phase update from overshooting available elastic energy but does not switch
    plasticity off below a yield stress.
    """
    seq, nxx, nyy, nxy = _von_mises_flow_direction(sigma_gp, mat.nu)
    rho = np.clip(np.asarray(rho_gp, float), rho_floor, rho_cap)
    rr = chain.rates(seq, rho, T_K)
    dep_prop = np.asarray(rr["dot_ep"], float) * max(float(dt_s), 0.0)

    sqrt23 = np.sqrt(2.0/3.0)
    dep_relax = 0.999 * seq / np.maximum(3.0*mat.G*sqrt23, 1e-30)
    dep = np.minimum(dep_prop, dep_relax)
    if np.isfinite(max_dep_phase) and max_dep_phase > 0:
        dep = np.minimum(dep, max_dep_phase)
    dep = np.maximum(dep, 0.0)

    dgamma = sqrt23 * dep
    ep_new = ep_gp.copy()
    ep_new[0] += 1.5 * dgamma * nxx
    ep_new[1] += 1.5 * dgamma * nyy
    ep_new[2] += 1.5 * dgamma * nxy

    # Kocks-Mecking state evolution is driven by accepted Arrhenius strain.
    # It is a transition rule for rho, not a stress law or yield criterion.
    drho = float(k_store) * np.sqrt(np.maximum(rho, 1e-30)) / max(mat.b, 1e-30) * dep \
           - float(k_dyn) * rho * dep
    if np.isfinite(max_rho_rel_phase) and max_rho_rel_phase > 0:
        drho = np.clip(drho, -max_rho_rel_phase*rho, max_rho_rel_phase*rho)
    rho_new = np.clip(rho + drho, rho_floor, rho_cap)

    seq_after = np.maximum(seq - 3.0*mat.G*dgamma, 0.0)
    dWp = 0.5*(seq + seq_after)*dep

    bd = chain.barrier_diagnostics(seq, rho, T_K)
    diag = {**rr, **bd, "seq_Pa": seq, "dep_eq": dep,
            "drho": rho_new-rho, "dWp_J_m3": dWp}
    return ep_new, rho_new, dep, dWp, diag


def _representative_plastic_cycle(
    mesh, bnd, cfg, Dmat, d, ep_gp, rho_gp,
    Umax, Umin, T_K, frequency_Hz, n_phase,
    plast_chain: ArrheniusPlasticChain, u_start,
    k_store, k_dyn, rho_floor, rho_cap,
    max_dep_phase, max_rho_rel_phase,
):
    """Integrate one physical cycle through the fully Arrhenius plastic chain."""
    phase = np.linspace(0.0, 2.0*np.pi, int(n_phase), endpoint=False)
    f = 0.5 * (1.0 + np.cos(phase))
    Uhist = Umin + (Umax - Umin) * f
    dt_phase = 1.0 / max(frequency_Hz * len(Uhist), 1e-30)

    ep0 = ep_gp.copy(); rho0 = rho_gp.copy()
    ep = ep_gp.copy(); rho = rho_gp.copy(); u = u_start.copy()
    dep_acc = np.zeros(mesh.ne)
    Wp_acc = np.zeros(mesh.ne)
    max_seq = np.zeros(mesh.ne)
    mu_emit = np.zeros(mesh.ne); mu_P = np.zeros(mesh.ne)
    mu_T = np.zeros(mesh.ne); mu_escape = np.zeros(mesh.ne)
    mu_flow = np.zeros(mesh.ne)
    phi_T_sum = np.zeros(mesh.ne)
    Gemit_sum = np.zeros(mesh.ne); GP_sum = np.zeros(mesh.ne); GT_sum = np.zeros(mesh.ne)

    for U in Uhist:
        K, Rint, *_ = assemble_mechanics(mesh, u, ep, rho, d, Dmat, cfg.material)
        u, _ = _solve_symmetric_tension(K, Rint, u, bnd, mesh, U, -U)
        sig, seq, _, _ = stress_state(mesh, u, ep, d, Dmat, cfg.material)
        max_seq = np.maximum(max_seq, seq)
        ep, rho, dep_phase, dWp_phase, diag = _arrhenius_chain_plastic_update(
            ep, rho, sig, cfg.material, T_K, dt_phase,
            plast_chain, k_store, k_dyn, rho_floor, rho_cap,
            max_dep_phase, max_rho_rel_phase,
        )
        dep_acc += np.asarray(dep_phase, float)
        Wp_acc += np.asarray(dWp_phase, float)
        mu_emit += np.asarray(diag["lambda_emit"], float) * dt_phase
        mu_P += np.asarray(diag["lambda_peierls"], float) * dt_phase
        mu_T += np.asarray(diag["lambda_taylor"], float) * dt_phase
        mu_escape += np.asarray(diag["lambda_escape"], float) * dt_phase
        mu_flow += np.asarray(diag["lambda_flow"], float) * dt_phase
        phi_T_sum += np.asarray(diag["phi_taylor"], float)
        Gemit_sum += np.asarray(diag["G_emit_eV"], float)
        GP_sum += np.asarray(diag["G_peierls_eV"], float)
        GT_sum += np.asarray(diag["G_taylor_eV"], float)

        # Re-equilibrate after the accepted plastic eigenstrain increment.
        K, Rint, *_ = assemble_mechanics(mesh, u, ep, rho, d, Dmat, cfg.material)
        u, _ = _solve_symmetric_tension(K, Rint, u, bnd, mesh, U, -U)

    nph = max(len(Uhist), 1)
    return {
        "dep_tensor_cycle": ep - ep0,
        "dep_eq_cycle": dep_acc,
        "drho_cycle": rho - rho0,
        "Wp_cycle_gp": Wp_acc,
        "u_end": u,
        "max_seq_gp": max_seq,
        "mu_emit_cycle_gp": mu_emit,
        "mu_peierls_cycle_gp": mu_P,
        "mu_taylor_cycle_gp": mu_T,
        "mu_escape_cycle_gp": mu_escape,
        "mu_flow_cycle_gp": mu_flow,
        "phi_taylor_mean_gp": phi_T_sum/nph,
        "G_emit_mean_eV_gp": Gemit_sum/nph,
        "G_peierls_mean_eV_gp": GP_sum/nph,
        "G_taylor_mean_eV_gp": GT_sum/nph,
    }


def _cycle_stress_histories(mesh, bnd, cfg, Dmat, d, ep_gp, rho_gp,
                            Umax, Umin, n_phase, u_start):
    """Elastic stress history over one cycle with current internal state fixed."""
    phase = np.linspace(0.0, 2.0*np.pi, int(n_phase), endpoint=False)
    f = 0.5 * (1.0 + np.cos(phase))
    Uhist = Umin + (Umax - Umin) * f
    u = u_start.copy()
    seq_nodes = []
    s1_nodes = []
    psi_nodes = []
    Ftop_hist = []
    for U in Uhist:
        K, Rint, *_ = assemble_mechanics(mesh, u, ep_gp, rho_gp, d, Dmat, cfg.material)
        u, Ft = _solve_symmetric_tension(K, Rint, u, bnd, mesh, U, -U)
        _, seq, s1, psi = stress_state(mesh, u, ep_gp, d, Dmat, cfg.material)
        seq_nodes.append(project_gp_to_nodes(mesh, seq))
        s1_nodes.append(project_gp_to_nodes(mesh, s1))
        psi_nodes.append(project_gp_to_nodes(mesh, psi))
        Ftop_hist.append(Ft)
    return {
        "seq_node": np.asarray(seq_nodes),
        "s1_node": np.asarray(s1_nodes),
        "psi_node": np.asarray(psi_nodes),
        "Ftop": np.asarray(Ftop_hist),
        "u_end": u,
    }


def _cycle_nucleation_hazard(crack_barrier, s1_hist, T_K, frequency_Hz,
                             state_shift_eV, sigma_back_node, chi,
                             multihit_m, tau_c_s):
    sig = np.maximum(np.asarray(s1_hist) - chi * sigma_back_node[None, :], 0.0)
    G = crack_barrier.deltaG_eV(sig, T_K)
    G = np.maximum(G + state_shift_eV[None, :], 1e-12)
    lam = crack_barrier.rate_prefactor * np.exp(
        np.clip(-G / max(KBEV*T_K, 1e-30), -700.0, 0.0)
    )
    if multihit_m > 1.0 + 1e-12:
        lam = gammainc(multihit_m, np.minimum(lam * max(tau_c_s, 1e-30), 1e12)) / max(tau_c_s, 1e-30)
    return np.mean(lam, axis=0) / max(frequency_Hz, 1e-30)


def _project_plastic_state(mesh, epsp_acc_gp, rho_gp, epsp_shield_scale, epsp_damage_scale):
    epsp_node = np.maximum(project_gp_to_nodes(mesh, epsp_acc_gp), 0.0)
    rho_node = np.maximum(project_gp_to_nodes(mesh, rho_gp), 0.0)
    P = 1.0 - np.exp(-epsp_node / max(epsp_shield_scale, 1e-30))
    Dloc = 1.0 - np.exp(-epsp_node / max(epsp_damage_scale, 1e-30))
    return epsp_node, rho_node, P, Dloc


def _surface_morphology_proposal(mesh, feature_nodes, dep_tensor_block,
                                 morph_band_length, normal_weight, shear_weight):
    """Convert block plastic strain into signed free-surface normal motion.

    The kinematic surrogate is
        dh = L_m [ w_n d eps_nn^p + w_s d gamma_nt^p ].
    The signed shear term naturally produces opposite-signed motion on the two
    sides of a symmetric blunt root and can create paired intrusion/extrusion
    morphology without imposing a pre-crack or a sharpening threshold.
    """
    idx = np.asarray(feature_nodes, dtype=int)
    t, n = feature_tangent_normal(mesh, idx)
    exx_n = project_gp_to_nodes(mesh, dep_tensor_block[0])[idx]
    eyy_n = project_gp_to_nodes(mesh, dep_tensor_block[1])[idx]
    gxy_n = project_gp_to_nodes(mesh, dep_tensor_block[2])[idx]
    E = np.zeros((len(idx), 2, 2), dtype=float)
    E[:,0,0] = exx_n; E[:,1,1] = eyy_n
    E[:,0,1] = 0.5*gxy_n; E[:,1,0] = 0.5*gxy_n
    Enn = np.einsum('ni,nij,nj->n', n, E, n)
    Ent = np.einsum('ni,nij,nj->n', n, E, t)
    gamma_nt = 2.0 * Ent

    # J2 plasticity does not resolve discrete slip bands, so use the local
    # equivalent plastic increment as the slip magnitude and gamma_nt only to
    # select the signed surface step direction.  This is an isotropic surrogate
    # for a persistent-slip-band offset: magnitude from accumulated slip, sign
    # from the resolved surface shear.  The tanh regularization avoids sign
    # chatter where gamma_nt is numerically tiny.
    dep_eq_gp = _equivalent_dep_from_tensor(dep_tensor_block)
    dep_eq_node = project_gp_to_nodes(mesh, dep_eq_gp)[idx]
    sign_reg = np.tanh(gamma_nt / np.maximum(0.05*dep_eq_node, 1e-16))
    signed_slip = dep_eq_node * sign_reg
    dh = morph_band_length * (normal_weight * Enn + shear_weight * signed_slip)
    return dh, Enn, gamma_nt


def _adjacency(mesh):
    adj = [set() for _ in range(mesh.nn)]
    for tri in mesh.elems:
        a,b,c = map(int, tri)
        adj[a].update((b,c)); adj[b].update((a,c)); adj[c].update((a,b))
    return adj


def _connected_crack_extent(mesh, d, root_xy, dthr, root_seed_radius, adj):
    root = np.asarray(root_xy)
    active = np.asarray(d >= dthr, bool)
    dist = np.linalg.norm(mesh.nodes - root[None,:], axis=1)
    seeds = np.where(active & (dist <= root_seed_radius))[0]
    if len(seeds) == 0:
        return 0.0, 0
    seen = set(map(int, seeds)); stack = list(map(int, seeds))
    while stack:
        i = stack.pop()
        for j in adj[i]:
            if active[j] and j not in seen:
                seen.add(j); stack.append(j)
    ii = np.fromiter(seen, dtype=int)
    extent = max(float(np.max(mesh.nodes[ii,0]) - root[0]), 0.0)
    return extent, len(ii)


def _plot_fields(mesh, fields, out_png: Path, root_xy, feature_nodes=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri
    tri = mtri.Triangulation(mesh.nodes[:,0]*1e3, mesh.nodes[:,1]*1e3, mesh.elems)
    n = len(fields)
    fig, axes = plt.subplots(1, n, figsize=(4.3*n, 4.0), constrained_layout=True)
    if n == 1: axes = [axes]
    for ax, (title, vals) in zip(axes, fields):
        pc = ax.tripcolor(tri, vals, shading="gouraud")
        ax.plot(root_xy[0]*1e3, root_xy[1]*1e3, "kx", ms=6)
        if feature_nodes is not None:
            q = mesh.nodes[np.asarray(feature_nodes, dtype=int)]
            ax.plot(q[:,0]*1e3, q[:,1]*1e3, "k-", lw=0.8)
        ax.set_aspect("equal"); ax.set_title(title); ax.set_xlabel("x (mm)"); ax.set_ylabel("y (mm)")
        fig.colorbar(pc, ax=ax)
    fig.savefig(out_png, dpi=220); plt.close(fig)


def run_case_stress(args, case_name: str, sigma_a_MPa: float):
    shield_on = case_name == "shielded"
    cfg = make_emergent_config()
    cfg.geometry.Lx = args.Lx; cfg.geometry.Ly = args.Ly
    geom = BluntNotchGeometry(args.Lx, args.Ly, args.notch_depth_m, args.notch_half_height_m)
    mesh, bnd, root_xy0 = make_blunt_edge_notch_mesh(
        geom, nx=args.nx, ny=args.ny, jitter=args.jitter,
        root_h_fine=args.root_h_fine, seed=args.seed,
    )
    feature_nodes = identify_feature_surface_nodes(mesh, geom)
    root_xy = local_root_xy(mesh, feature_nodes)
    root_radius0 = local_root_radius(mesh, feature_nodes)
    fixed_mesh_nodes = np.unique(np.r_[bnd.top_nodes, bnd.bot_nodes])

    Dmat = plane_strain_D(cfg.material)
    Md, Kd = assemble_pf_matrices(mesh)
    ell = args.ell_m if args.ell_m is not None else args.ell_factor * mesh.hbar_tip
    fracture_model = FractureModel(cfg.fracture_barrier, cfg.material, pf=cfg.phase_field)
    Gc_base = float(fracture_model.Gc_of_T(args.T, ell, method=args.Gc_method))
    plast_chain = build_chain_from_namespace(args, cfg.material.b)
    _, crack = make_barriers(-40.0, args.S_crack_kB, args.emit_energy_scale)

    sigma_max = 2.0 * sigma_a_MPa * 1e6 / max(1.0 - args.R, 1e-30)
    sigma_min = args.R * sigma_max

    u = np.zeros(mesh.ndof)
    ep_gp = np.zeros((3, mesh.ne))
    rho_gp = np.full(mesh.ne, args.rho0)
    epsp_acc_gp = np.zeros(mesh.ne)
    d = np.zeros(mesh.nn)
    Hhist = np.zeros(mesh.nn)
    B_nuc = np.zeros(mesh.nn)
    Gc_eff = np.full(mesh.nn, Gc_base)
    adj = _adjacency(mesh)

    outdir = Path(args.out) / case_name / (f"sigmaA_{sigma_a_MPa:g}MPa".replace(".", "p"))
    outdir.mkdir(parents=True, exist_ok=True)
    with (outdir / "run_args.json").open("w") as f:
        json.dump(vars(args) | {"case":case_name, "sigma_a_MPa":sigma_a_MPa,
                                "root_radius_initial_m":root_radius0, "ell_m":ell},
                  f, indent=2, sort_keys=True)

    cycles = 0.0; cycles_clock = None; cycles_pf = None
    Wp_total = 0.0; rows = []
    last_s1 = np.zeros(mesh.nn); last_residual = np.zeros(mesh.nn)

    for ib in range(args.max_blocks):
        if cycles >= args.cycles_max:
            break

        # Stress-control displacement calibration for the current geometry and eigenstrain state.
        Umax, Umin, u_zero, F0, F_U_slope = _affine_stress_control_displacements(
            mesh, bnd, cfg, Dmat, ep_gp, rho_gp, d, sigma_max, sigma_min, u
        )

        # Actual constitutive representative cycle.  This updates copies of ep and rho;
        # the resulting one-cycle increments are block-extrapolated below.
        cyc = _representative_plastic_cycle(
            mesh, bnd, cfg, Dmat, d, ep_gp, rho_gp, Umax, Umin,
            args.T, args.frequency_Hz, args.plastic_n_phase,
            plast_chain, u_zero,
            args.k_store, args.k_dyn, args.rho_floor, args.rho_cap,
            args.max_dep_phase, args.max_rho_rel_phase,
        )
        dep_tensor_cycle = cyc["dep_tensor_cycle"]
        dep_eq_cycle = np.maximum(cyc["dep_eq_cycle"], 0.0)
        drho_cycle = cyc["drho_cycle"]

        # Pre-block crack hazard with current internal state and true residual-stress-shifted cycle.
        epsp_node, rho_node, P, Dloc = _project_plastic_state(
            mesh, epsp_acc_gp, rho_gp, args.epsp_shield_scale, args.epsp_damage_scale
        )
        chi = args.shield_chi if shield_on else 0.0
        Gsh = args.Gshield_eV if shield_on else 0.0
        sigma_back = args.sigma_back_max_GPa * 1e9 * P
        state_shift = Gsh * P - args.Gstored_eV * Dloc
        hist_pre = _cycle_stress_histories(
            mesh, bnd, cfg, Dmat, d, ep_gp, rho_gp,
            Umax, Umin, args.hazard_n_phase, u_zero,
        )
        mu_nuc_pre = _cycle_nucleation_hazard(
            crack, hist_pre["s1_node"], args.T, args.frequency_Hz,
            state_shift, sigma_back, chi, args.multihit_m, args.multihit_tau_s,
        )

        remaining = args.cycles_max - cycles
        dN = min(args.block_cycles, remaining)
        max_dep_cycle = float(np.max(dep_eq_cycle)) if dep_eq_cycle.size else 0.0
        if max_dep_cycle > 0:
            dN = min(dN, args.target_dep_eq_block / max_dep_cycle)
        rel_rho_cycle = np.max(np.abs(drho_cycle) / np.maximum(rho_gp, args.rho0)) if rho_gp.size else 0.0
        if rel_rho_cycle > 0:
            dN = min(dN, args.target_rho_rel_block / rel_rho_cycle)
        max_mu_n = float(np.max(mu_nuc_pre)) if mu_nuc_pre.size else 0.0
        if max_mu_n > 0:
            dN = min(dN, args.target_dB_nuc / max_mu_n)

        # Geometry-resolution limiter from the signed morphology proposal.
        if args.enable_geometry_evolution and max_dep_cycle > 0:
            dh_cycle, _, _ = _surface_morphology_proposal(
                mesh, feature_nodes, dep_tensor_cycle,
                args.morph_band_length_m, args.morph_normal_weight, args.morph_shear_weight,
            )
            max_dh_cycle = float(np.max(np.abs(dh_cycle))) if dh_cycle.size else 0.0
            if max_dh_cycle > 0:
                dN = min(dN, args.target_surface_move_fraction * mesh.hbar_tip / max_dh_cycle)

        dN = max(min(dN, remaining), args.min_block_cycles)
        dN = min(dN, remaining)
        if dN <= 0:
            break

        # Accept actual plastic eigenstrain and dislocation-density increments.
        dep_tensor_block = dep_tensor_cycle * dN
        dep_eq_block = dep_eq_cycle * dN
        ep_gp = ep_gp + dep_tensor_block
        rho_gp = np.clip(rho_gp + drho_cycle * dN, args.rho_floor, args.rho_cap)
        epsp_acc_gp = epsp_acc_gp + dep_eq_block
        Wp_total += float(np.sum(cyc["Wp_cycle_gp"] * mesh.area_e) * dN)

        # Evolve the free-surface geometry from the signed surface plastic strain.
        mesh_scale = 0.0
        dh_block = np.zeros(len(feature_nodes))
        gamma_nt = np.zeros(len(feature_nodes))
        if args.enable_geometry_evolution:
            dh_block, _, gamma_nt = _surface_morphology_proposal(
                mesh, feature_nodes, dep_tensor_block,
                args.morph_band_length_m, args.morph_normal_weight, args.morph_shear_weight,
            )
            mesh_scale = apply_local_ale_surface_update(
                mesh, feature_nodes, dh_block,
                decay_length=args.morph_decay_length_m,
                fixed_nodes=fixed_mesh_nodes,
                max_move=args.max_surface_move_fraction * mesh.hbar_tip,
                min_area_fraction=args.min_area_fraction,
            )
            if mesh_scale > 0:
                root_xy = local_root_xy(mesh, feature_nodes)
                rebuild_mesh_geometry(mesh, root_xy)
                Md, Kd = assemble_pf_matrices(mesh)
                adj = _adjacency(mesh)

        # Recalibrate stress control after state/geometry evolution.
        Umax2, Umin2, u_zero2, F02, slope2 = _affine_stress_control_displacements(
            mesh, bnd, cfg, Dmat, ep_gp, rho_gp, d, sigma_max, sigma_min, u_zero
        )

        # Zero-load residual stress field: the eigenstrain is now part of FEM equilibrium.
        _, seq_res, s1_res, _ = stress_state(mesh, u_zero2, ep_gp, d, Dmat, cfg.material)
        residual_node = project_gp_to_nodes(mesh, s1_res)
        last_residual = residual_node.copy()

        # Post-block cycle hazard with updated plasticity, hardening, residual stress and geometry.
        epsp_node, rho_node, P, Dloc = _project_plastic_state(
            mesh, epsp_acc_gp, rho_gp, args.epsp_shield_scale, args.epsp_damage_scale
        )
        sigma_back = args.sigma_back_max_GPa * 1e9 * P
        state_shift = Gsh * P - args.Gstored_eV * Dloc
        hist_post = _cycle_stress_histories(
            mesh, bnd, cfg, Dmat, d, ep_gp, rho_gp,
            Umax2, Umin2, args.hazard_n_phase, u_zero2,
        )
        mu_nuc_post = _cycle_nucleation_hazard(
            crack, hist_post["s1_node"], args.T, args.frequency_Hz,
            state_shift, sigma_back, chi, args.multihit_m, args.multihit_tau_s,
        )
        dB_nuc = 0.5 * (mu_nuc_pre + mu_nuc_post) * dN
        B_old = B_nuc.copy(); B_nuc += dB_nuc

        # Update phase-field resistance from actual plastic strain state.
        log_ratio = (args.shield_log_gain * P if shield_on else 0.0) - args.damage_log_drop * Dloc
        Gc_eff = Gc_base * np.exp(log_ratio)

        # Phase-field crack nucleation/extension at the actual evolved geometry.
        p_block = 1.0 - np.exp(-np.clip(dB_nuc, 0.0, 700.0))
        psi_node = np.max(hist_post["psi_node"], axis=0)
        for _ in range(args.n_stagger):
            d, Hhist = update_phase_field(
                d=d, Hhist=Hhist, psi_e_node=psi_node,
                Md=Md, Kd=Kd, notch_nodes=bnd.notch_nodes,
                Gc_eff=Gc_eff, ell=ell, Gamma0=1.0, dt=1.0,
                crack_hazard_probability=p_block,
                use_kinetic_drive=False,
                max_damage_increment=args.max_damage_increment,
                damage_drive_cap=args.damage_drive_cap,
            )

        cycles_old = cycles; cycles += dN
        if cycles_clock is None and float(np.max(B_nuc)) >= 1.0:
            crossing = (B_old < 1.0) & (B_nuc >= 1.0) & (dB_nuc > 0)
            if np.any(crossing):
                frac = (1.0 - B_old[crossing]) / np.maximum(dB_nuc[crossing], 1e-300)
                cycles_clock = cycles_old + float(np.clip(np.min(frac), 0.0, 1.0))*dN
            else:
                cycles_clock = cycles

        root_xy = local_root_xy(mesh, feature_nodes)
        root_radius = local_root_radius(mesh, feature_nodes)
        extent, nconn = _connected_crack_extent(
            mesh, d, root_xy, args.pf_damage_threshold,
            args.root_seed_radius_factor * ell, adj,
        )
        if cycles_pf is None and extent >= args.pf_crack_extent_factor * ell:
            cycles_pf = cycles

        last_s1 = np.max(hist_post["s1_node"], axis=0)
        row = {
            "block": ib, "case": case_name, "sigma_a_MPa": sigma_a_MPa,
            "cycles_total": cycles, "dN": dN,
            "Umax_m": Umax2, "Umin_m": Umin2,
            "F0_residual_N_per_m": F02,
            "dep_eq_cycle_max": max_dep_cycle,
            "dep_eq_block_max": float(np.max(dep_eq_block)),
            "mu_emit_cycle_max": float(np.max(cyc["mu_emit_cycle_gp"])),
            "mu_peierls_cycle_max": float(np.max(cyc["mu_peierls_cycle_gp"])),
            "mu_taylor_cycle_max": float(np.max(cyc["mu_taylor_cycle_gp"])),
            "mu_escape_cycle_max": float(np.max(cyc["mu_escape_cycle_gp"])),
            "mu_flow_cycle_max": float(np.max(cyc["mu_flow_cycle_gp"])),
            "phi_taylor_mean_max": float(np.max(cyc["phi_taylor_mean_gp"])),
            "G_emit_mean_eV_min": float(np.min(cyc["G_emit_mean_eV_gp"])),
            "G_peierls_mean_eV_min": float(np.min(cyc["G_peierls_mean_eV_gp"])),
            "G_taylor_mean_eV_min": float(np.min(cyc["G_taylor_mean_eV_gp"])),
            "epsp_acc_max": float(np.max(epsp_acc_gp)),
            "rho_min_m2": float(np.min(rho_gp)), "rho_max_m2": float(np.max(rho_gp)),
            "rho_mean_m2": float(np.mean(rho_gp)),
            "residual_sigma1_max_Pa": float(np.max(residual_node)),
            "residual_sigma1_min_Pa": float(np.min(residual_node)),
            "sigma1_cycle_max_Pa": float(np.max(last_s1)),
            "P_max": float(np.max(P)), "Dloc_max": float(np.max(Dloc)),
            "mu_nuc_pre_max": float(np.max(mu_nuc_pre)),
            "mu_nuc_post_max": float(np.max(mu_nuc_post)),
            "dB_nuc_max": float(np.max(dB_nuc)), "B_nuc_max": float(np.max(B_nuc)),
            "Gc_min_J_m2": float(np.min(Gc_eff)), "Gc_max_J_m2": float(np.max(Gc_eff)),
            "d_max": float(np.max(d)), "connected_crack_extent_m": extent,
            "connected_crack_nodes": nconn,
            "root_x_m": root_xy[0], "root_y_m": root_xy[1],
            "root_radius_m": root_radius,
            "root_radius_over_initial": root_radius / max(root_radius0, 1e-30) if np.isfinite(root_radius) else np.nan,
            "precursor_sharpness_ell_over_rho": ell / max(root_radius, 1e-30) if np.isfinite(root_radius) else np.nan,
            "surface_move_max_m": float(np.max(np.abs(dh_block))) if dh_block.size else 0.0,
            "surface_shear_gamma_nt_max": float(np.max(np.abs(gamma_nt))) if gamma_nt.size else 0.0,
            "ale_accept_scale": mesh_scale,
            "plastic_work_J_per_m": Wp_total,
            "cycles_to_nucleation_clock": cycles_clock if cycles_clock is not None else np.nan,
            "cycles_to_pf_crack": cycles_pf if cycles_pf is not None else np.nan,
            "at2_surface_energy_J_per_m": float(at2_surface_energy(mesh, d, ell, Gc_eff)),
        }
        rows.append(row)

        if args.print_every and ib % args.print_every == 0:
            print(
                f"FULLPLASTIC {case_name} sigma_a={sigma_a_MPa:g} MPa block={ib} "
                f"N={cycles:.3e} dN={dN:.3g} dep={row['dep_eq_block_max']:.2e} "
                f"rho={row['rho_max_m2']:.2e} res={row['residual_sigma1_max_Pa']*1e-6:.1f}MPa "
                f"r/r0={row['root_radius_over_initial']:.3g} Bnu={row['B_nuc_max']:.3g} d={row['d_max']:.3g}"
            )

        if ib == 0 or (args.snapshot_every > 0 and ib % args.snapshot_every == 0):
            _plot_fields(
                mesh,
                [("phase field d", d), ("accumulated eps_p", epsp_node),
                 ("rho (m^-2)", rho_node), ("residual sigma1 (MPa)", residual_node*1e-6)],
                outdir / f"fields_block_{ib:05d}.png", root_xy, feature_nodes,
            )

        u = hist_post["u_end"].copy()
        if args.stop_after_pf_crack and cycles_pf is not None:
            break

    _write_csv(outdir / "sn_pf2d_fullplastic_history.csv", rows)
    epsp_node, rho_node, P, Dloc = _project_plastic_state(
        mesh, epsp_acc_gp, rho_gp, args.epsp_shield_scale, args.epsp_damage_scale
    )
    _plot_fields(
        mesh,
        [("phase field d", d), ("accumulated eps_p", epsp_node),
         ("rho (m^-2)", rho_node), ("residual sigma1 (MPa)", last_residual*1e-6)],
        outdir / "fields_final.png", root_xy, feature_nodes,
    )

    summary = {
        "model":"SN_PF2D_fully_Arrhenius_emission_Peierls_Taylor_evolving_surface",
        "case":case_name, "sigma_a_MPa":sigma_a_MPa,
        "T_K":args.T, "R":args.R, "frequency_Hz":args.frequency_Hz,
        "cycles_total":cycles, "cycles_to_nucleation_clock":cycles_clock,
        "cycles_to_pf_crack":cycles_pf,
        "status":"pf_crack" if cycles_pf is not None else ("clock_first_passage" if cycles_clock is not None else "right_censored"),
        "root_radius_initial_m":root_radius0,
        "root_radius_final_m":local_root_radius(mesh, feature_nodes),
        "root_x_final_m":root_xy[0], "root_y_final_m":root_xy[1],
        "ell_m":ell, "Gc_base_J_m2":Gc_base,
        "chi_back":chi, "Gshield_eV":Gsh,
        "emit_energy_scale":args.emit_energy_scale,
        "emit_entropy_scale":args.emit_entropy_scale,
        "peierls_energy_scale":args.peierls_energy_scale,
        "peierls_entropy_scale":args.peierls_entropy_scale,
        "taylor_energy_scale":args.taylor_energy_scale,
        "taylor_entropy_scale":args.taylor_entropy_scale,
        "epsp_acc_final_max":float(np.max(epsp_acc_gp)),
        "rho_final_max_m2":float(np.max(rho_gp)),
        "rho_final_mean_m2":float(np.mean(rho_gp)),
        "B_nuc_final_max":float(np.max(B_nuc)),
        "P_final_max":float(np.max(P)), "Dloc_final_max":float(np.max(Dloc)),
        "d_final_max":float(np.max(d)),
        "residual_sigma1_final_max_Pa":float(np.max(last_residual)),
        "plastic_work_final_J_per_m":Wp_total,
        "history_csv":str(outdir / "sn_pf2d_fullplastic_history.csv"),
    }
    with (outdir / "summary.json").open("w") as f:
        json.dump(summary, f, indent=2, sort_keys=True)
    return summary


def run_sweep(args):
    Path(args.out).mkdir(parents=True, exist_ok=True)
    summaries = []
    for case in ("no_shield", "shielded"):
        for s in args.sigma_a_MPa:
            print(f"=== FULL-PLASTIC 2D case={case} sigma_a={s:g} MPa ===")
            summaries.append(run_case_stress(args, case, float(s)))
    _write_csv(Path(args.out) / "sn_pf2d_fullplastic_summary.csv", summaries)
    return summaries


def build_parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", default="runs/sn_pf2d_fullplastic_two_case")
    p.add_argument("--T", type=float, default=300.0)
    p.add_argument("--sigma-a-MPa", nargs="+", type=float, default=[500,600,700], dest="sigma_a_MPa")
    p.add_argument("--R", type=float, default=0.1)
    p.add_argument("--frequency-Hz", type=float, default=1000.0, dest="frequency_Hz")
    p.add_argument("--cycles-max", type=float, default=1e9, dest="cycles_max")
    p.add_argument("--block-cycles", type=float, default=1e7, dest="block_cycles")
    p.add_argument("--min-block-cycles", type=float, default=1e-6, dest="min_block_cycles")
    p.add_argument("--max-blocks", type=int, default=3000, dest="max_blocks")

    p.add_argument("--target-dep-eq-block", type=float, default=2e-4, dest="target_dep_eq_block")
    p.add_argument("--target-rho-rel-block", type=float, default=0.05, dest="target_rho_rel_block")
    p.add_argument("--target-dB-nuc", type=float, default=0.05, dest="target_dB_nuc")
    p.add_argument("--plastic-n-phase", type=int, default=12, dest="plastic_n_phase")
    p.add_argument("--hazard-n-phase", type=int, default=16, dest="hazard_n_phase")

    p.add_argument("--Lx", type=float, default=2e-3)
    p.add_argument("--Ly", type=float, default=4e-3)
    p.add_argument("--notch-depth-m", type=float, default=0.15e-3, dest="notch_depth_m")
    p.add_argument("--notch-half-height-m", type=float, default=0.30e-3, dest="notch_half_height_m")
    p.add_argument("--nx", type=int, default=36); p.add_argument("--ny", type=int, default=72)
    p.add_argument("--jitter", type=float, default=0.08)
    p.add_argument("--root-h-fine", type=float, default=30e-6, dest="root_h_fine")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--ell-m", type=float, default=None, dest="ell_m")
    p.add_argument("--ell-factor", type=float, default=3.0, dest="ell_factor")
    p.add_argument("--Gc-method", choices=["lambertw","hazard"], default="lambertw", dest="Gc_method")

    p.add_argument("--rho0", type=float, default=1e12)
    p.add_argument("--rho-floor", type=float, default=1e8, dest="rho_floor")
    p.add_argument("--rho-cap", type=float, default=1e17, dest="rho_cap")
    p.add_argument("--k-store", type=float, default=np.sqrt(2.0), dest="k_store")
    p.add_argument("--k-dyn", type=float, default=1.0, dest="k_dyn")

    # Fully Arrhenius emission -> Peierls -> Taylor plastic-event chain.
    # Defaults match the selected case-64-M1 fatigue scaling.
    p.add_argument("--exp-system", default="W[100]", dest="exp_system")
    p.add_argument("--exp-G00-eV", type=float, default=None, dest="exp_G00_eV")
    p.add_argument("--exp-gT-eV-per-K", type=float, default=None, dest="exp_gT_eV_per_K")
    p.add_argument("--exp-sigc0-GPa", type=float, default=None, dest="exp_sigc0_GPa")
    p.add_argument("--exp-sT-MPa-per-K", type=float, default=None, dest="exp_sT_MPa_per_K")
    p.add_argument("--exp-Tref-K", type=float, default=None, dest="exp_Tref_K")
    p.add_argument("--exp-floor-frac", type=float, default=None, dest="exp_floor_frac")
    p.add_argument("--exp-a", type=float, default=None, dest="exp_a")
    p.add_argument("--exp-n", type=float, default=None, dest="exp_n")
    p.add_argument("--emit-energy-scale", type=float, default=0.75, dest="emit_energy_scale")
    p.add_argument("--emit-entropy-scale", type=float, default=0.75, dest="emit_entropy_scale")
    p.add_argument("--emit-stress-scale", type=float, default=1.0, dest="emit_stress_scale")
    p.add_argument("--peierls-energy-scale", type=float, default=0.00375, dest="peierls_energy_scale")
    p.add_argument("--peierls-entropy-scale", type=float, default=0.00375, dest="peierls_entropy_scale")
    p.add_argument("--peierls-stress-scale", type=float, default=1.0, dest="peierls_stress_scale")
    p.add_argument("--taylor-energy-scale", type=float, default=0.015, dest="taylor_energy_scale")
    p.add_argument("--taylor-entropy-scale", type=float, default=0.015, dest="taylor_entropy_scale")
    p.add_argument("--taylor-stress-scale", type=float, default=1.0, dest="taylor_stress_scale")
    p.add_argument("--nu0-emit-pz", type=float, default=1e11, dest="nu0_emit_pz")
    p.add_argument("--nu0-peierls", type=float, default=1e11, dest="nu0_peierls")
    p.add_argument("--nu0-taylor", type=float, default=1e11, dest="nu0_taylor")
    p.add_argument("--plastic-event-strain", type=float, default=1e-5, dest="plastic_event_strain")
    p.add_argument("--phi-taylor-max", type=float, default=20.0, dest="phi_taylor_max")
    p.add_argument("--max-dep-phase", type=float, default=2e-5, dest="max_dep_phase")
    p.add_argument("--max-rho-rel-phase", type=float, default=0.02, dest="max_rho_rel_phase")

    p.add_argument("--epsp-shield-scale", type=float, default=5e-3, dest="epsp_shield_scale")
    p.add_argument("--epsp-damage-scale", type=float, default=2e-2, dest="epsp_damage_scale")
    p.add_argument("--S-crack-kB", type=float, default=0.0, dest="S_crack_kB")
    p.add_argument("--sigma-back-max-GPa", type=float, default=1.0, dest="sigma_back_max_GPa")
    p.add_argument("--shield-chi", type=float, default=0.6, dest="shield_chi")
    p.add_argument("--Gshield-eV", type=float, default=0.35, dest="Gshield_eV")
    p.add_argument("--Gstored-eV", type=float, default=0.25, dest="Gstored_eV")
    p.add_argument("--shield-log-gain", type=float, default=0.8, dest="shield_log_gain")
    p.add_argument("--damage-log-drop", type=float, default=1.5, dest="damage_log_drop")
    p.add_argument("--multihit-m", type=float, default=3.0, dest="multihit_m")
    p.add_argument("--multihit-tau-s", type=float, default=1e-6, dest="multihit_tau_s")

    p.add_argument("--enable-geometry-evolution", action="store_true", default=True, dest="enable_geometry_evolution")
    p.add_argument("--disable-geometry-evolution", action="store_false", dest="enable_geometry_evolution")
    p.add_argument("--morph-band-length-m", type=float, default=100e-6, dest="morph_band_length_m")
    p.add_argument("--morph-normal-weight", type=float, default=0.25, dest="morph_normal_weight")
    p.add_argument("--morph-shear-weight", type=float, default=1.0, dest="morph_shear_weight")
    p.add_argument("--morph-decay-length-m", type=float, default=150e-6, dest="morph_decay_length_m")
    p.add_argument("--target-surface-move-fraction", type=float, default=0.10, dest="target_surface_move_fraction")
    p.add_argument("--max-surface-move-fraction", type=float, default=0.20, dest="max_surface_move_fraction")
    p.add_argument("--min-area-fraction", type=float, default=0.15, dest="min_area_fraction")

    p.add_argument("--n-stagger", type=int, default=2, dest="n_stagger")
    p.add_argument("--max-damage-increment", type=float, default=0.05, dest="max_damage_increment")
    p.add_argument("--damage-drive-cap", type=float, default=20.0, dest="damage_drive_cap")
    p.add_argument("--pf-damage-threshold", type=float, default=0.5, dest="pf_damage_threshold")
    p.add_argument("--root-seed-radius-factor", type=float, default=3.0, dest="root_seed_radius_factor")
    p.add_argument("--pf-crack-extent-factor", type=float, default=3.0, dest="pf_crack_extent_factor")
    p.add_argument("--stop-after-pf-crack", action="store_true", default=True, dest="stop_after_pf_crack")
    p.add_argument("--snapshot-every", type=int, default=50, dest="snapshot_every")
    p.add_argument("--print-every", type=int, default=1, dest="print_every")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    run_sweep(args)


if __name__ == "__main__":
    main()
