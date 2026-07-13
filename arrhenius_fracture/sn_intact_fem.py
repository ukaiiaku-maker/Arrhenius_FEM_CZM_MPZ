"""Shared intact-FEM helpers for blunt-feature S--N and stateful-PD models.

This module restores the public helper layer expected by
``sn_pd2d_stateful.py`` and ``stateful_peridynamics.py``.  The algorithms are
factored from the existing, tested intact-FEM implementation in
``sn_pf2d_fullplastic.py``; no phase-field degradation is used here.  The
module intentionally keeps the prior anisotropic/spatial drivers separate from
crack-tip MPZ constitutive changes.
"""
from __future__ import annotations

from typing import Dict, Tuple

import numpy as np
from scipy.sparse.linalg import spsolve

from .config import ElasticProperties
from .fem import (
    assemble_mechanics,
    plane_strain_D,
    project_gp_to_nodes,
    stress_state,
)
from .sn_arrhenius_chain import ArrheniusPlasticChain
from .sn_geometry import feature_tangent_normal

__all__ = [
    "affine_stress_control_displacements",
    "cycle_stress_histories",
    "lumped_nodal_area",
    "plane_strain_D",
    "project_gp_to_nodes",
    "project_plastic_state",
    "representative_plastic_cycle",
    "stress_state_intact",
    "surface_morphology_proposal",
]


def _intact_damage(mesh) -> np.ndarray:
    """Return the identically intact nodal degradation field."""
    return np.zeros(mesh.nn, dtype=float)


def _rho_placeholder(mesh) -> np.ndarray:
    """Return a finite GP array required by the common FEM assembly API.

    ``assemble_mechanics`` currently accepts ``rho_gp`` for interface
    compatibility but does not use it in the elastic assembly.
    """
    return np.zeros(mesh.ne, dtype=float)


def lumped_nodal_area(mesh) -> np.ndarray:
    """Area associated with each node by conservative P1 mass lumping.

    Each linear triangular element contributes one third of its physical area
    to each of its three vertices.  Consequently ``sum(A_node)`` equals the
    total mesh area to floating-point precision, which is required when a
    finite density of candidate sites is mapped onto the local PD patch.
    """
    area = np.zeros(mesh.nn, dtype=float)
    conn = np.asarray(mesh.elems, dtype=int)
    contrib = np.asarray(mesh.area_e, dtype=float) / 3.0
    np.add.at(area, conn[:, 0], contrib)
    np.add.at(area, conn[:, 1], contrib)
    np.add.at(area, conn[:, 2], contrib)
    if np.any(~np.isfinite(area)) or np.any(area < 0.0):
        raise ValueError("invalid lumped nodal area")
    return area


def _solve_symmetric_tension(K, Rint, u, bnd, mesh, Uy_top, Uy_bot):
    """Incremental symmetric-displacement solve with one horizontal anchor."""
    prescribed = np.zeros(mesh.ndof, dtype=bool)
    target = np.zeros(mesh.ndof, dtype=float)

    prescribed[2 * np.asarray(bnd.top_nodes, dtype=int) + 1] = True
    target[2 * np.asarray(bnd.top_nodes, dtype=int) + 1] = float(Uy_top)
    prescribed[2 * np.asarray(bnd.bot_nodes, dtype=int) + 1] = True
    target[2 * np.asarray(bnd.bot_nodes, dtype=int) + 1] = float(Uy_bot)

    x, y = mesh.nodes[:, 0], mesh.nodes[:, 1]
    anchor = int(np.argmin((x - np.max(x)) ** 2 + y**2))
    prescribed[2 * anchor] = True
    target[2 * anchor] = 0.0

    free = ~prescribed
    Kc = K.tocsr()
    du_prescribed = target[prescribed] - u[prescribed]
    rhs = -Rint[free] - Kc[np.ix_(free, prescribed)] @ du_prescribed

    u_new = np.asarray(u, dtype=float).copy()
    u_new[free] = u[free] + spsolve(Kc[np.ix_(free, free)], rhs)
    u_new[prescribed] = target[prescribed]
    Rfull = Rint + Kc @ (u_new - u)
    Ftop = float(np.sum(Rfull[2 * np.asarray(bnd.top_nodes, dtype=int) + 1]))
    return u_new, Ftop


def _assemble_intact(mesh, u, ep_gp, Dmat, mat):
    return assemble_mechanics(
        mesh,
        u,
        ep_gp,
        _rho_placeholder(mesh),
        _intact_damage(mesh),
        Dmat,
        mat,
    )


def affine_stress_control_displacements(
    mesh,
    bnd,
    mat: ElasticProperties,
    Dmat: np.ndarray,
    ep_gp: np.ndarray,
    sigma_max_Pa: float,
    sigma_min_Pa: float,
    u_guess: np.ndarray,
) -> Tuple[float, float, np.ndarray, float, float]:
    """Calibrate symmetric boundary displacement to target nominal stresses.

    At fixed plastic eigenstrain the small-strain problem is affine in imposed
    displacement.  A zero-displacement equilibrium and one probe solve define
    the reaction/displacement map exactly.
    """
    K, Rint, *_ = _assemble_intact(mesh, u_guess, ep_gp, Dmat, mat)
    u_zero, F0 = _solve_symmetric_tension(K, Rint, u_guess, bnd, mesh, 0.0, 0.0)

    K0, R0, *_ = _assemble_intact(mesh, u_zero, ep_gp, Dmat, mat)
    Uprobe = max(float(mesh.hbar_tip), 1.0e-8) * 1.0e-3
    _, Fprobe = _solve_symmetric_tension(
        K0, R0, u_zero, bnd, mesh, Uprobe, -Uprobe
    )
    slope = (Fprobe - F0) / max(Uprobe, 1.0e-30)
    if not np.isfinite(slope) or abs(slope) < 1.0e-30:
        raise RuntimeError("stress-control displacement calibration failed")

    # In the 2-D plane-strain model, reaction has units N/m.  Division by the
    # specimen width normal to the loading direction produces nominal stress.
    width_2d = max(
        float(np.max(mesh.nodes[:, 0]) - np.min(mesh.nodes[:, 0])), 1.0e-30
    )
    Fmax = float(sigma_max_Pa) * width_2d
    Fmin = float(sigma_min_Pa) * width_2d
    Umax = (Fmax - F0) / slope
    Umin = (Fmin - F0) / slope
    return float(Umax), float(Umin), u_zero, float(F0), float(slope)


def stress_state_intact(mesh, u, ep_gp, Dmat, mat):
    """Stress-state wrapper for an intact body (no phase-field degradation)."""
    return stress_state(mesh, u, ep_gp, _intact_damage(mesh), Dmat, mat)


def _von_mises_flow_direction(sigma_gp: np.ndarray, nu: float):
    sx, sy, txy = sigma_gp[0], sigma_gp[1], sigma_gp[2]
    szz = nu * (sx + sy)
    mean = (sx + sy + szz) / 3.0
    sxx = sx - mean
    syy = sy - mean
    szzd = szz - mean
    seq = np.sqrt(
        np.maximum(
            1.5 * (sxx * sxx + syy * syy + szzd * szzd + 2.0 * txy * txy),
            0.0,
        )
    )
    norm = np.sqrt(
        np.maximum(sxx * sxx + syy * syy + szzd * szzd + 2.0 * txy * txy, 1e-60)
    )
    return seq, sxx / norm, syy / norm, txy / norm


def _arrhenius_chain_plastic_update(
    ep_gp,
    rho_gp,
    sigma_gp,
    mat,
    T_K,
    dt_s,
    chain: ArrheniusPlasticChain,
    k_store,
    k_dyn,
    rho_floor,
    rho_cap,
    max_dep_phase,
    max_rho_rel_phase,
):
    """Advance one resolved phase using the serial Arrhenius plastic chain."""
    seq, nxx, nyy, nxy = _von_mises_flow_direction(sigma_gp, mat.nu)
    rho = np.clip(np.asarray(rho_gp, dtype=float), rho_floor, rho_cap)
    rates = chain.rates(seq, rho, T_K)
    dep_proposed = np.asarray(rates["dot_ep"], dtype=float) * max(float(dt_s), 0.0)

    sqrt23 = np.sqrt(2.0 / 3.0)
    dep_relax = 0.999 * seq / np.maximum(3.0 * mat.G * sqrt23, 1e-30)
    dep = np.minimum(dep_proposed, dep_relax)
    if np.isfinite(max_dep_phase) and max_dep_phase > 0.0:
        dep = np.minimum(dep, max_dep_phase)
    dep = np.maximum(dep, 0.0)

    dgamma = sqrt23 * dep
    ep_new = np.asarray(ep_gp, dtype=float).copy()
    ep_new[0] += 1.5 * dgamma * nxx
    ep_new[1] += 1.5 * dgamma * nyy
    ep_new[2] += 1.5 * dgamma * nxy

    drho = (
        float(k_store)
        * np.sqrt(np.maximum(rho, 1e-30))
        / max(mat.b, 1e-30)
        * dep
        - float(k_dyn) * rho * dep
    )
    if np.isfinite(max_rho_rel_phase) and max_rho_rel_phase > 0.0:
        drho = np.clip(
            drho, -max_rho_rel_phase * rho, max_rho_rel_phase * rho
        )
    rho_new = np.clip(rho + drho, rho_floor, rho_cap)

    seq_after = np.maximum(seq - 3.0 * mat.G * dgamma, 0.0)
    dWp = 0.5 * (seq + seq_after) * dep
    barrier = chain.barrier_diagnostics(seq, rho, T_K)
    diag = {
        **rates,
        **barrier,
        "seq_Pa": seq,
        "dep_eq": dep,
        "drho": rho_new - rho,
        "dWp_J_m3": dWp,
    }
    return ep_new, rho_new, dep, dWp, diag


def representative_plastic_cycle(
    mesh,
    bnd,
    mat: ElasticProperties,
    Dmat,
    ep_gp,
    rho_gp,
    Umax,
    Umin,
    T_K,
    frequency_Hz,
    n_phase,
    plast_chain: ArrheniusPlasticChain,
    u_start,
    k_store,
    k_dyn,
    rho_floor,
    rho_cap,
    max_dep_phase,
    max_rho_rel_phase,
) -> Dict[str, np.ndarray]:
    """Resolve one physical cycle and return accepted state increments."""
    n_phase = max(int(n_phase), 2)
    phase = np.linspace(0.0, 2.0 * np.pi, n_phase, endpoint=False)
    load_fraction = 0.5 * (1.0 + np.cos(phase))
    Uhist = Umin + (Umax - Umin) * load_fraction
    dt_phase = 1.0 / max(float(frequency_Hz) * len(Uhist), 1e-30)

    ep0 = np.asarray(ep_gp, dtype=float).copy()
    rho0 = np.asarray(rho_gp, dtype=float).copy()
    ep = ep0.copy()
    rho = rho0.copy()
    u = np.asarray(u_start, dtype=float).copy()

    dep_acc = np.zeros(mesh.ne)
    Wp_acc = np.zeros(mesh.ne)
    max_seq = np.zeros(mesh.ne)
    mu_emit = np.zeros(mesh.ne)
    mu_peierls = np.zeros(mesh.ne)
    mu_taylor = np.zeros(mesh.ne)
    mu_escape = np.zeros(mesh.ne)
    mu_flow = np.zeros(mesh.ne)
    phi_sum = np.zeros(mesh.ne)
    Ge_sum = np.zeros(mesh.ne)
    Gp_sum = np.zeros(mesh.ne)
    Gt_sum = np.zeros(mesh.ne)

    for U in Uhist:
        K, Rint, *_ = _assemble_intact(mesh, u, ep, Dmat, mat)
        u, _ = _solve_symmetric_tension(K, Rint, u, bnd, mesh, U, -U)
        sigma, seq, _, _ = stress_state_intact(mesh, u, ep, Dmat, mat)
        max_seq = np.maximum(max_seq, seq)
        ep, rho, dep_phase, dWp_phase, diag = _arrhenius_chain_plastic_update(
            ep,
            rho,
            sigma,
            mat,
            T_K,
            dt_phase,
            plast_chain,
            k_store,
            k_dyn,
            rho_floor,
            rho_cap,
            max_dep_phase,
            max_rho_rel_phase,
        )
        dep_acc += dep_phase
        Wp_acc += dWp_phase
        mu_emit += np.asarray(diag["lambda_emit"]) * dt_phase
        mu_peierls += np.asarray(diag["lambda_peierls"]) * dt_phase
        mu_taylor += np.asarray(diag["lambda_taylor"]) * dt_phase
        mu_escape += np.asarray(diag["lambda_escape"]) * dt_phase
        mu_flow += np.asarray(diag["lambda_flow"]) * dt_phase
        phi_sum += np.asarray(diag["phi_taylor"])
        Ge_sum += np.asarray(diag["G_emit_eV"])
        Gp_sum += np.asarray(diag["G_peierls_eV"])
        Gt_sum += np.asarray(diag["G_taylor_eV"])

        # Re-equilibrate after accepting the plastic eigenstrain increment.
        K, Rint, *_ = _assemble_intact(mesh, u, ep, Dmat, mat)
        u, _ = _solve_symmetric_tension(K, Rint, u, bnd, mesh, U, -U)

    nph = float(len(Uhist))
    return {
        "dep_tensor_cycle": ep - ep0,
        "dep_eq_cycle": dep_acc,
        "drho_cycle": rho - rho0,
        "Wp_cycle_gp": Wp_acc,
        "u_end": u,
        "max_seq_gp": max_seq,
        "mu_emit_cycle_gp": mu_emit,
        "mu_peierls_cycle_gp": mu_peierls,
        "mu_taylor_cycle_gp": mu_taylor,
        "mu_escape_cycle_gp": mu_escape,
        "mu_flow_cycle_gp": mu_flow,
        "phi_taylor_mean_gp": phi_sum / nph,
        "G_emit_mean_eV_gp": Ge_sum / nph,
        "G_peierls_mean_eV_gp": Gp_sum / nph,
        "G_taylor_mean_eV_gp": Gt_sum / nph,
    }


def cycle_stress_histories(
    mesh,
    bnd,
    mat: ElasticProperties,
    Dmat,
    ep_gp,
    Umax,
    Umin,
    n_phase,
    u_start,
):
    """Elastic stress tensor history over a cycle at fixed internal state."""
    n_phase = max(int(n_phase), 2)
    phase = np.linspace(0.0, 2.0 * np.pi, n_phase, endpoint=False)
    load_fraction = 0.5 * (1.0 + np.cos(phase))
    Uhist = Umin + (Umax - Umin) * load_fraction
    u = np.asarray(u_start, dtype=float).copy()

    sigma_nodes = []
    seq_nodes = []
    s1_nodes = []
    psi_nodes = []
    Ftop_hist = []
    u_at_max = None
    imax = int(np.argmax(Uhist))

    for i, U in enumerate(Uhist):
        K, Rint, *_ = _assemble_intact(mesh, u, ep_gp, Dmat, mat)
        u, Ftop = _solve_symmetric_tension(K, Rint, u, bnd, mesh, U, -U)
        sigma, seq, s1, psi = stress_state_intact(mesh, u, ep_gp, Dmat, mat)
        sigma_nodes.append(project_gp_to_nodes(mesh, sigma))
        seq_nodes.append(project_gp_to_nodes(mesh, seq))
        s1_nodes.append(project_gp_to_nodes(mesh, s1))
        psi_nodes.append(project_gp_to_nodes(mesh, psi))
        Ftop_hist.append(Ftop)
        if i == imax:
            u_at_max = u.copy()

    if u_at_max is None:
        u_at_max = u.copy()
    return {
        "sigma_node": np.asarray(sigma_nodes),
        "seq_node": np.asarray(seq_nodes),
        "s1_node": np.asarray(s1_nodes),
        "psi_node": np.asarray(psi_nodes),
        "Ftop": np.asarray(Ftop_hist),
        "u_max": u_at_max,
        "u_end": u,
    }


def project_plastic_state(
    mesh,
    epsp_acc_gp,
    rho_gp,
    epsp_shield_scale,
    epsp_damage_scale,
):
    """Project accumulated plastic state and bounded closure coordinates."""
    epsp_node = np.maximum(project_gp_to_nodes(mesh, epsp_acc_gp), 0.0)
    rho_node = np.maximum(project_gp_to_nodes(mesh, rho_gp), 0.0)
    P = 1.0 - np.exp(-epsp_node / max(float(epsp_shield_scale), 1e-30))
    Dloc = 1.0 - np.exp(-epsp_node / max(float(epsp_damage_scale), 1e-30))
    return epsp_node, rho_node, P, Dloc


def _equivalent_dep_from_tensor(dep):
    exx, eyy, gxy = dep[0], dep[1], dep[2]
    ezz = -(exx + eyy)
    tensor_sq = exx**2 + eyy**2 + ezz**2 + 0.5 * gxy**2
    return np.sqrt(np.maximum((2.0 / 3.0) * tensor_sq, 0.0))


def surface_morphology_proposal(
    mesh,
    feature_nodes,
    dep_tensor_block,
    morph_band_length,
    normal_weight,
    shear_weight,
):
    """Convert accepted plastic strain into signed free-surface motion."""
    idx = np.asarray(feature_nodes, dtype=int)
    tangent, normal = feature_tangent_normal(mesh, idx)
    exx = project_gp_to_nodes(mesh, dep_tensor_block[0])[idx]
    eyy = project_gp_to_nodes(mesh, dep_tensor_block[1])[idx]
    gxy = project_gp_to_nodes(mesh, dep_tensor_block[2])[idx]

    strain = np.zeros((len(idx), 2, 2), dtype=float)
    strain[:, 0, 0] = exx
    strain[:, 1, 1] = eyy
    strain[:, 0, 1] = 0.5 * gxy
    strain[:, 1, 0] = 0.5 * gxy
    Enn = np.einsum("ni,nij,nj->n", normal, strain, normal)
    Ent = np.einsum("ni,nij,nj->n", normal, strain, tangent)
    gamma_nt = 2.0 * Ent

    dep_eq_gp = _equivalent_dep_from_tensor(dep_tensor_block)
    dep_eq_node = project_gp_to_nodes(mesh, dep_eq_gp)[idx]
    sign_regularized = np.tanh(
        gamma_nt / np.maximum(0.05 * dep_eq_node, 1e-16)
    )
    signed_slip = dep_eq_node * sign_regularized
    dh = float(morph_band_length) * (
        float(normal_weight) * Enn + float(shear_weight) * signed_slip
    )
    return dh, Enn, gamma_nt
