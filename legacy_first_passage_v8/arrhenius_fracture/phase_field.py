"""
Phase-field fracture with AT2 functional and Model-A kinetics.

Key fixes from original:
1. Clean Model-A (Allen-Cahn) kinetics: Gamma * (dF/dd)
   No hybrid AT2 + hazard overlay. Temperature enters ONLY through Gc(T).
2. Proper irreversibility via history field H_hist
3. Crack-front detection for diagnostics
4. AT2 surface energy functional for crack length measurement
"""

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
from typing import Tuple, Optional

from .mesh import TriMesh
from .config import TipMemoryConfig


def update_phase_field(
    d: np.ndarray,               # current damage field
    Hhist: np.ndarray,           # history field (max elastic energy density)
    psi_e_node: np.ndarray,      # current elastic energy density at nodes
    Md: sparse.csr_matrix,       # PF mass matrix
    Kd: sparse.csr_matrix,       # PF stiffness matrix
    notch_nodes: np.ndarray,     # nodes in initial notch (d=1)
    Gc_eff,                      # fracture energy: scalar or (nn,) array
    ell: float,                  # PF length scale
    Gamma0: float,               # constant kinetic mobility (1/s)
    dt: float,                   # time step
    cohesive_gate: Optional[np.ndarray] = None,  # optional cohesive traction gate
    drive_multiplier: Optional[np.ndarray] = None,  # optional M_tip^2-like stress/energy amplification
    crack_hazard_probability: Optional[np.ndarray] = None,  # optional Arrhenius crack-growth probability
    crack_fired_nodes: Optional[np.ndarray] = None,  # bool mask: completed first-passage events -> d driven to 1 (source mode)
    use_kinetic_drive: bool = False,
    max_damage_increment: Optional[float] = None,
    damage_drive_cap: float = 20.0,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Update damage field using AT2 phase-field with Model-A kinetics.

    Gc_eff can be a scalar (uniform) or a per-node array (spatially varying).
    When spatially varying, regions with higher Gc (e.g., from plastic dissipation)
    resist fracture more strongly, producing emergent toughening.

    Returns
    -------
    d_new : updated damage field
    Hhist : updated history field
    """
    nn = len(d)

    # Update history field (ensures irreversibility of driving force)
    Hhist = np.maximum(Hhist, psi_e_node)

    # Handle scalar or vector Gc
    Gc = np.broadcast_to(np.asarray(Gc_eff, dtype=float), (nn,)).copy()
    Gc = np.maximum(Gc, 1e-30)

    # Build and solve AT2 system with spatially varying Gc
    Hc = Gc / (2 * ell)  # (nn,) critical energy density

    # Effective driving force (with optional cohesive gate)
    drive = np.maximum(Hhist / Hc - 1, 0)
    if drive_multiplier is not None:
        dm = np.broadcast_to(np.asarray(drive_multiplier, dtype=float), (nn,))
        dm = np.clip(dm, 0.0, 1e3)
        drive = drive * dm

    if cohesive_gate is not None:
        drive = drive * cohesive_gate

    drive = np.clip(drive, 0.0, damage_drive_cap)

    # Variational AT2 solve with spatially varying Gc:
    # [diag(Gc/ell) @ Md + diag(Gc*ell) @ Kd + 2*diag(Hhist) @ Md] d = 2*Md*Hhist
    A_pf = (sparse.diags(Gc / ell) @ Md +
            sparse.diags(Gc * ell) @ Kd +
            2 * sparse.diags(Hhist) @ Md)

    rhs = 2 * Md @ Hhist

    # Free nodes (not in notch)
    free = np.ones(nn, dtype=bool)
    free[notch_nodes] = False

    d_pf = d.copy()
    d_pf[free] = spsolve(A_pf[np.ix_(free, free)], rhs[free])
    d_pf[notch_nodes] = 1.0

    # Model-A semi-implicit kinetic update.  In the quasi-static AT2 mode
    # use_kinetic_drive=False, the update is simply the irreversible
    # variational AT2 solution d_pf.  The explicit drive term is retained as
    # an optional dynamic/hazard-like accelerator, but it must be used with a
    # physically meaningful dt*Gamma0; otherwise d jumps directly to 1.
    if use_kinetic_drive:
        dt_eff = dt * Gamma0
        Gamma_eff = drive * dt_eff
        d_trial = (d_pf + Gamma_eff) / (1 + Gamma_eff)
    else:
        d_trial = d_pf

    # Optional Arrhenius/Eyring crack-growth hazard.  The variational AT2
    # solution remains the thermodynamic target, but crack advance toward that
    # target is accepted only with the local event probability supplied by the
    # shielded process-zone crack-growth law.  This replaces an instantaneous
    # deterministic AT2 jump with a kinetic crack-growth clock without adding a
    # post-hoc toughness gate.
    if crack_hazard_probability is not None:
        P = np.broadcast_to(np.asarray(crack_hazard_probability, dtype=float), (nn,))
        P = np.clip(P, 0.0, 1.0)
        d_trial = d + P * (d_trial - d)
        # SOURCE mode: a completed first-passage event IS the sub-grid bond
        # rupture -- the material point is broken regardless of whether the
        # ell-SMEARED energy H/Hc has reached the variational AT2 threshold.
        # In 'gate' mode (legacy) P only modulates advance toward the
        # variational target d_pf, so a fired front with d_pf ~ d stalls: the
        # de-smeared sigma_tip lives only in the hazard while the AT2 target
        # still sees the smoothed field.  Fired nodes are NOT driven to d=1 by
        # fiat here -- that injects unpaid AT2 surface energy (an unrelaxed
        # d-spike costs ~Gc*ell in gradient energy) and breaks the global
        # audit.  Instead 'source' mode collapses Gc_local at fired nodes
        # (cleavage = persistent loss of resistance, applied in main.py), and
        # THIS variational solve breaks them only if the elastic field can
        # pay -- which makes the Griffith condition emergent.

    # Irreversibility + bounds
    d_new = np.clip(d_trial, d, 1.0)

    # Optional damage-rate cap per stagger iteration.  This is a numerical
    # continuation device: it does not change the variational target, but it
    # prevents losing the entire propagation history in one load step.
    if max_damage_increment is not None and np.isfinite(max_damage_increment):
        if max_damage_increment > 0:
            d_capped = np.minimum(d_new, d + max_damage_increment)
            # Fired first-passage nodes bypass the continuation cap: the
            # event is discrete (the sub-grid bond is broken this step).
            if crack_fired_nodes is not None:
                fired = np.asarray(crack_fired_nodes, dtype=bool)
                d_new = np.where(fired, d_new, d_capped)
            else:
                d_new = d_capped

    d_new = np.clip(d_new, d, 1.0)
    d_new[notch_nodes] = 1.0

    return d_new, Hhist


def at2_surface_energy(
    mesh: TriMesh, d: np.ndarray, ell: float, Gc
) -> float:
    """
    Compute AT2 crack surface energy functional.

    E = ∫ Gc(x) * [d²/(2ℓ) + (ℓ/2)|∇d|²] dV  (per unit thickness)

    Gc can be scalar or per-node array.
    """
    Gc_arr = np.broadcast_to(np.asarray(Gc, dtype=float), (mesh.nn,))

    E = 0.0
    for e in range(mesh.ne):
        conn = mesh.elems[e]
        de = d[conn]
        A = mesh.area_e[e]
        Gc_e = np.mean(Gc_arr[conn])

        d1, d2, d3 = de
        int_d2 = A / 6 * (d1**2 + d2**2 + d3**2 +
                          d1*d2 + d2*d3 + d3*d1)

        dNdx = mesh.dNdx_e[e]
        grad_d = dNdx @ de
        grad2 = grad_d[0]**2 + grad_d[1]**2
        int_g2 = A * grad2

        E += Gc_e * ((1 / (2 * ell)) * int_d2 + (ell / 2) * int_g2)

    return E


def crack_front_mask(
    mesh: TriMesh, d: np.ndarray, ell: float
) -> np.ndarray:
    """
    Identify crack-front nodes: high |∇d| and intermediate d.

    Returns a smooth weight field (0 = far from front, 1 = at front).
    """
    nn = mesh.nn

    # Compute nodal |∇d| by area-weighted projection
    gd_acc = np.zeros(nn)
    w_acc = np.zeros(nn)

    for e in range(mesh.ne):
        conn = mesh.elems[e]
        de = d[conn]
        A = mesh.area_e[e]
        dNdx = mesh.dNdx_e[e]
        grad_d = dNdx @ de
        gdmag = np.sqrt(grad_d[0]**2 + grad_d[1]**2)

        for a in range(3):
            gd_acc[conn[a]] += gdmag * A / 3
            w_acc[conn[a]] += A / 3

    gd_node = gd_acc / np.maximum(w_acc, 1e-30)

    # Front indicator: |∇d| * (1-d)
    front_ind = gd_node * np.maximum(1 - d, 0)
    fmax = np.max(front_ind)

    if fmax < 1e-12:
        return np.zeros(nn)

    # Smooth mask around front seeds
    seeds = np.where(front_ind > 0.3 * fmax)[0]
    if len(seeds) == 0:
        return np.zeros(nn)

    Xs = mesh.nodes[seeds]
    dist = np.full(nn, np.inf)

    # Chunk-wise distance computation
    chunk = 2000
    for i in range(0, nn, chunk):
        idx = slice(i, min(i + chunk, nn))
        Xi = mesh.nodes[idx]
        dx = Xi[:, 0:1] - Xs[:, 0].T
        dy = Xi[:, 1:2] - Xs[:, 1].T
        dist[idx] = np.sqrt(np.min(dx**2 + dy**2, axis=1))

    R_front = 3 * ell
    mask = _smoothstep(1 - dist / R_front)
    return mask



def memory_energy_density(
    rtip_state: np.ndarray,
    shield_state: np.ndarray,
    front_w: np.ndarray,
    ell: float,
    rtip_ref: float,
    Gc_node,
    cfg: TipMemoryConfig,
) -> np.ndarray:
    """Local memory storage density [J/m^3].

    Minimal thermodynamic closure for the crack-tip memory state.  The energy
    scale is Gc/ell, and the dimensionless storage is quadratic in the blunting
    radius state and shielding/backstress state.  This supplies conjugate forces
    A_r = -dPsi/dr and A_z = -dPsi/dz for diagnostics and energy auditing.
    """
    Gc = np.broadcast_to(np.asarray(Gc_node if Gc_node is not None else 1.0, dtype=float), rtip_state.shape)
    e_scale = np.maximum(Gc, 1e-30) / max(ell, 1e-30)
    kr = max(float(getattr(cfg, 'memory_energy_r_coeff', 0.0)), 0.0)
    kz = max(float(getattr(cfg, 'memory_energy_z_coeff', 0.0)), 0.0)
    x = (np.asarray(rtip_state, dtype=float) - rtip_ref) / max(ell, 1e-30)
    z = np.asarray(shield_state, dtype=float)
    fw = np.clip(np.asarray(front_w, dtype=float), 0.0, 1.0)
    return fw * e_scale * (0.5 * kr * x*x + 0.5 * kz * z*z)

def update_tip_memory(
    d: np.ndarray, d_old: np.ndarray,
    dot_ep_node: np.ndarray, dt: float,
    mesh: TriMesh, ell: float,
    rtip_state: np.ndarray,      # (nn,) local tip radius
    shield_state: np.ndarray,    # (nn,) local shielding
    cfg: TipMemoryConfig,
    rtip_ref: float,
    dwp_node: Optional[np.ndarray] = None,       # plastic work increment [J/m^3]
    Gc_node=None,                              # scalar or (nn,) local fracture energy [J/m^2]
    emit_prob_node: Optional[np.ndarray] = None, # optional crack-tip emission probability [0,1]
    crack_advance: float = 0.0,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, dict]:
    """
    Update node-local crack-tip memory variables.

    This is the reduced Stage-1 memory closure.  Plastic flow and plastic
    work blunt the tip and build shielding; crack/damage advance sharpens the
    tip and can add wake/branch shielding.  Optional dislocation-emission
    probability contributes to the same blunting/shielding state.

    The state modifies only M_tip, i.e. the local fracture stress/energy
    amplification used by the phase-field or cohesive gate.
    """
    nn = len(d)

    mode = getattr(cfg, 'mode', 'stage1')
    enabled = bool(getattr(cfg, 'enabled', True)) and mode != 'off'
    if not enabled:
        info = {
            'front_weight_mean': 0.0, 'dwp_norm_front_mean': 0.0,
            'emit_prob_front_mean': 0.0, 'state_gain': 0.0,
            'dEmem_node': np.zeros(nn), 'dDmem_node': np.zeros(nn),
            'Emem_node': np.zeros(nn), 'memory_A_r_mean': 0.0, 'memory_A_z_mean': 0.0,
        }
        return (rtip_state, shield_state, np.ones(nn), np.ones(nn), info)

    # Global ablation gain.  The explicit mode setting overrides state_gain
    # so parameter sweeps can use a clean off/weak/base comparison.
    if mode == 'weak_stage1':
        gain = 0.25
    else:
        gain = float(getattr(cfg, 'state_gain', 1.0))
    gain = max(gain, 0.0)

    # Front weighting
    front_w = crack_front_mask(mesh, d, ell)
    front_w = np.maximum(front_w, ((d > 0.02) & (d < 0.98)).astype(float))
    front_w = np.clip(front_w, 0, 1)

    # Memory storage before the update.  The conjugate is defined with respect
    # to this quadratic storage plus the way M_tip modifies elastic/fracture
    # driving.  The latter is represented implicitly by the accepted update;
    # the explicit storage/dissipation audit prevents memory from acting as a
    # free, non-dissipative toughness multiplier.
    rtip_old = rtip_state.copy()
    shield_old = shield_state.copy()
    if bool(getattr(cfg, 'use_memory_energetics', True)):
        Emem_old_node = memory_energy_density(rtip_old, shield_old, front_w, ell, rtip_ref,
                                              Gc_node if Gc_node is not None else 1.0, cfg)
    else:
        Emem_old_node = np.zeros(nn)

    # Kinematic increments
    dEp = dt * np.maximum(dot_ep_node, 0.0)
    dDamage = np.maximum(d - d_old, 0.0)

    # Normalized plastic work increment: dWp*ell/Gc.  This is dimensionless
    # and ties memory growth to the ratio of plastic work to cleavage work.
    if dwp_node is None:
        dWp_norm = np.zeros(nn)
    else:
        dWp = np.maximum(np.asarray(dwp_node, dtype=float), 0.0)
        Gc = np.broadcast_to(np.asarray(Gc_node if Gc_node is not None else 1.0, dtype=float), (nn,))
        dWp_norm = dWp * ell / np.maximum(Gc, 1e-30)
        dWp_norm = np.clip(dWp_norm, 0.0, 10.0)

    if emit_prob_node is None:
        P_emit = np.zeros(nn)
    else:
        P_emit = np.clip(np.asarray(emit_prob_node, dtype=float), 0.0, 1.0)

    # Optional wake-memory decay under crack advance.  This prevents old,
    # fully damaged material points from dominating global mean diagnostics.
    Lwake = max(getattr(cfg, 'wake_length_factor', 8.0) * ell, 1e-30)
    if np.isfinite(Lwake) and crack_advance > 0:
        decay = np.exp(-max(crack_advance, 0.0) / Lwake)
        shield_state = shield_state * decay

    # Wake relaxation OUTSIDE the active front: nodes that were once "tip"
    # (sharpened rtip, accumulated shield) but are no longer at the front
    # relax back toward the reference state at fractional rate
    # wake_relax*(1-front_w) per step.  This removes the M_tip ratchet web
    # (old front rings retaining amplification forever and re-firing).
    wr = max(float(getattr(cfg, 'wake_relax', 0.0)), 0.0)
    if wr > 0.0:
        relax = np.clip(wr * (1.0 - np.clip(front_w, 0.0, 1.0)), 0.0, 1.0)
        rtip_state = rtip_state + relax * (rtip_ref - rtip_state)
        shield_state = shield_state * (1.0 - relax)

    # Tip radius update
    rtip_min = max(cfg.rtip_min_factor * ell, 1e-12)
    rtip_max = max(cfg.rtip_max_factor * ell, 1.01 * rtip_min)

    blunt_drive = (
        cfg.blunt_per_plastic_strain * dEp +
        getattr(cfg, 'blunt_per_work', 0.0) * dWp_norm +
        getattr(cfg, 'blunt_per_emission', 0.0) * P_emit
    )
    sharp_drive = cfg.sharpen_per_damage * dDamage

    dr_tip = gain * ell * front_w * (blunt_drive - sharp_drive)
    rtip_state = np.clip(rtip_state + dr_tip, rtip_min, rtip_max)

    # Shielding update
    dz = gain * front_w * (
        cfg.shield_from_plastic * dEp +
        getattr(cfg, 'shield_from_work', 0.0) * dWp_norm +
        cfg.shield_from_damage * dDamage +
        getattr(cfg, 'shield_from_emission', 0.0) * P_emit
    )
    shield_state = np.clip(shield_state + dz, 0.0, cfg.shield_max)

    # Memory storage and dissipation increments.  D_m is rate-independent here
    # (linear in |dr| and |dz|) so it behaves like a toughening/process-zone
    # dissipation rather than a purely elastic storage term.
    if bool(getattr(cfg, 'use_memory_energetics', True)):
        Gc_arr = np.broadcast_to(np.asarray(Gc_node if Gc_node is not None else 1.0, dtype=float), (nn,))
        e_scale = np.maximum(Gc_arr, 1e-30) / max(ell, 1e-30)
        Emem_node = memory_energy_density(rtip_state, shield_state, front_w, ell, rtip_ref, Gc_arr, cfg)
        dEmem_node = Emem_node - Emem_old_node
        Rr = max(float(getattr(cfg, 'memory_dissipation_r_coeff', 0.0)), 0.0)
        Rz = max(float(getattr(cfg, 'memory_dissipation_z_coeff', 0.0)), 0.0)
        dDmem_node = front_w * e_scale * (Rr * np.abs(rtip_state - rtip_old) / max(ell, 1e-30)
                                          + Rz * np.abs(shield_state - shield_old))
        kr = max(float(getattr(cfg, 'memory_energy_r_coeff', 0.0)), 0.0)
        kz = max(float(getattr(cfg, 'memory_energy_z_coeff', 0.0)), 0.0)
        A_r = -front_w * e_scale * kr * ((rtip_state - rtip_ref) / max(ell, 1e-30)) / max(ell, 1e-30)
        A_z = -front_w * e_scale * kz * shield_state
    else:
        Emem_node = np.zeros(nn); dEmem_node = np.zeros(nn); dDmem_node = np.zeros(nn)
        A_r = np.zeros(nn); A_z = np.zeros(nn)

    # Amplification factor
    rtip_amp = np.sqrt(max(rtip_ref, 1e-12) / np.maximum(rtip_state, rtip_min))
    rtip_amp = np.clip(rtip_amp, cfg.amp_min, cfg.amp_max)

    # Shield factor attenuates only excess amplification, not the baseline M=1.
    shield_min = 1 - cfg.shield_max
    shield_factor = np.maximum(1 - shield_state, shield_min)

    # Small diagnostics dictionary for the caller.
    fw_sum = np.sum(front_w)
    if fw_sum > 0:
        dwp_front = float(np.sum(front_w * dWp_norm) / fw_sum)
        emit_front = float(np.sum(front_w * P_emit) / fw_sum)
    else:
        dwp_front = 0.0
        emit_front = 0.0
    info = {
        'front_weight_mean': float(np.mean(front_w)),
        'dwp_norm_front_mean': dwp_front,
        'emit_prob_front_mean': emit_front,
        'state_gain': gain,
        'Emem_node': Emem_node,
        'dEmem_node': dEmem_node,
        'dDmem_node': dDmem_node,
        'memory_A_r_mean': float(np.sum(front_w * A_r) / max(np.sum(front_w), 1e-30)),
        'memory_A_z_mean': float(np.sum(front_w * A_z) / max(np.sum(front_w), 1e-30)),
    }

    return rtip_state, shield_state, rtip_amp, shield_factor, info


def tip_emission_probability(
    sigma_pos: np.ndarray,
    T: float,
    dt: float,
    plast_model,
    cohesive_cfg,
) -> np.ndarray:
    """
    Arrhenius probability for a local dislocation-emission/blunting event.

    This deliberately reuses the plasticity barrier shape, scaled by the
    cohesive/tip-emission factors.  Thus the emission memory is tied mainly
    to H0 and v*, rather than introducing an independent fitted barrier.
    """
    sigma = np.maximum(np.asarray(sigma_pos, dtype=float), 0.0)
    if T <= 0 or not getattr(cohesive_cfg, 'use_emission', False):
        return np.zeros_like(sigma)

    H_emit = getattr(cohesive_cfg, 'emit_H_factor', 0.55) * plast_model.H(sigma)
    v_emit = getattr(cohesive_cfg, 'emit_v_factor', 1.0) * plast_model.v(sigma, T)
    G_emit = np.maximum(H_emit - sigma * v_emit, 0.0)

    eta0 = max(getattr(cohesive_cfg, 'emit_eta0', 1e13), 1e-300)
    log_rate = np.log(eta0) - G_emit / (1.380649e-23 * T)
    log_rate = np.clip(log_rate, -745.0, 80.0)
    rate = np.exp(log_rate)
    P = 1.0 - np.exp(-np.maximum(dt, 0.0) * rate)
    return np.clip(P, 0.0, 1.0)

def compute_fracture_amplification(
    sigma1_node: np.ndarray,
    d: np.ndarray,
    mesh: TriMesh,
    ell: float,
    rtip_amp: np.ndarray,
    shield_factor: np.ndarray,
    M_max: float = 4.0,
    lambda_tip: float = 5.0,
    kappa_tip_max: float = 4.0,
) -> np.ndarray:
    """
    Compute local stress amplification for fracture driving force.

    M_fract = 1 + (M_base - 1) * rtip_amp * shield_factor

    where M_base comes from the gradient-based crack-tip indicator.
    """
    nn = mesh.nn

    # Gradient-based tip indicator
    gd_node = _nodal_grad_d(mesh, d)
    tip_ind = gd_node * np.maximum(1 - d, 0)
    k_tip = np.minimum(1 + lambda_tip * tip_ind, kappa_tip_max)

    # Base amplification
    M_base = np.minimum(k_tip, M_max)

    # Apply memory
    M_fract = 1 + (M_base - 1) * rtip_amp * shield_factor
    M_fract = np.clip(M_fract, 1.0, M_max)

    return M_fract


def cohesive_gate(
    sigma_pos: np.ndarray, sigma_coh: float,
    width: float = 0.08, floor: float = 0.0
) -> np.ndarray:
    """
    Smooth logistic gate for cohesive traction criterion.

    Returns ~0 when sigma < sigma_coh, ~1 when sigma > sigma_coh.
    """
    ratio = sigma_pos / max(sigma_coh, 1e-30)
    x = (ratio - 1.0) / max(width, 1e-12)
    x = np.clip(x, -80, 80)
    gate_raw = 1.0 / (1 + np.exp(-x))

    # Normalize
    g0 = 1.0 / (1 + np.exp(1.0 / max(width, 1e-12)))
    gate_raw = np.clip((gate_raw - g0) / max(1 - g0, 1e-12), 0, 1)

    return floor + (1 - floor) * gate_raw


def _nodal_grad_d(mesh: TriMesh, d: np.ndarray) -> np.ndarray:
    """Compute |∇d| projected to nodes."""
    nn = mesh.nn
    acc = np.zeros(nn)
    wacc = np.zeros(nn)

    for e in range(mesh.ne):
        conn = mesh.elems[e]
        de = d[conn]
        A = mesh.area_e[e]
        dNdx = mesh.dNdx_e[e]
        grad_d = dNdx @ de
        gdmag = np.sqrt(grad_d[0]**2 + grad_d[1]**2)

        for a in range(3):
            acc[conn[a]] += gdmag * A / 3
            wacc[conn[a]] += A / 3

    return acc / np.maximum(wacc, 1e-30)


def _smoothstep(x: np.ndarray) -> np.ndarray:
    """Hermite smoothstep: 0 for x<=0, 1 for x>=1."""
    x = np.clip(x, 0, 1)
    return x**2 * (3 - 2 * x)

