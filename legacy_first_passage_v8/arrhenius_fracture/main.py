"""
Main simulation driver for Arrhenius phase-field fracture.

This is the entry point that:
1. Builds the mesh (once)
2. Loops over temperatures
3. At each T: runs the staggered PF+mechanics+plasticity solve
4. Computes proper J-integral toughness
5. Saves diagnostics and results

Usage:
    python -m arrhenius_fracture.main                    # default DBTT preset
    python -m arrhenius_fracture.main --preset ceramic   # ceramic preset
    python -m arrhenius_fracture.main --preset cohesive  # cohesive DBTT
"""

import numpy as np
import os
import time
import argparse
import json
from typing import Dict
import copy

from .config import (
    SimulationConfig, make_dbtt_config, make_ceramic_config,
    make_cohesive_dbtt_config, make_emergent_config, KB, EV_TO_J,
)
from .mesh import make_tri_mesh, make_boundary_data
from .materials import PlasticityModel, FractureModel
from .fem import (
    plane_strain_D, assemble_mechanics, solve_dirichlet,
    project_gp_to_nodes, assemble_pf_matrices,
    boundary_reaction_forces, elastic_energy_densities,
    stress_state)
from .plasticity import update_plasticity, plastic_flow_diagnostics, arrhenius_taylor_flow_stress, calibrate_peierls_floor, peierls_flow_stress
from .phase_field import (
    update_phase_field, at2_surface_energy, crack_front_mask,
    update_tip_memory, compute_fracture_amplification, cohesive_gate,
    tip_emission_probability,
)
from .at1 import update_phase_field_at1, at1_critical_stress, at1_surface_energy
from .j_integral import compute_J_integral, find_crack_tip, compute_crack_advance
from .diagnostics import (
    StepDiagnostics, SimulationHistory,
    save_history, plot_diagnostics, plot_toughness_vs_T,
    save_step_table, save_summary_json, plot_field_snapshots,
    save_results_summary_table,
)


def edge_crack_force_K(F_per_thickness, geometry, material):
    """Force-based LEFM sanity-check K for a single-edge crack in tension.

    Parameters
    ----------
    F_per_thickness : float
        Reaction force per out-of-plane thickness [N/m].
    geometry : GeometryConfig
        Uses Lx as specimen width W and a0 as the edge-crack length a.
    material : ElasticProperties
        Included for signature consistency; not used directly because K is
        computed from force and geometry.

    Returns
    -------
    float
        |K_I| estimate [Pa*sqrt(m)]. This is a diagnostic only, not the
        selected toughness metric.
    """
    W = max(float(getattr(geometry, 'Lx', 0.0)), 1e-30)
    a = max(float(getattr(geometry, 'a0', 0.0)), 1e-30)
    alpha = np.clip(a / W, 1e-6, 0.95)

    # Tada/Paris/Irwin polynomial for a single-edge crack in a finite-width
    # plate under uniform tension, valid as a robust diagnostic for alpha<~0.6.
    Y = (1.12 - 0.231*alpha + 10.55*alpha**2
         - 21.72*alpha**3 + 30.39*alpha**4)

    sigma_nom = abs(float(F_per_thickness)) / W  # Pa for unit thickness
    return float(Y * sigma_nom * np.sqrt(np.pi * a))


def _node_lumped_area(mesh):
    """Lumped nodal area weights for integrating node fields per unit thickness."""
    w = np.zeros(mesh.nn)
    for e, conn in enumerate(mesh.elems):
        w[conn] += mesh.area_e[e] / 3.0
    return w

def _toughening_weight(mesh, d, ell, front_only=True, wake_weight=0.25):
    """Smooth process-zone weight for plastic-work-driven fracture resistance.

    The weight localizes retained toughening to the crack front/process zone so
    remote plastic work does not directly raise Gc everywhere.  A small wake
    term can retain shielding behind the front, but the fully intact far-field
    and fully damaged crack wake are suppressed.
    """
    d = np.asarray(d, dtype=float)
    front = crack_front_mask(mesh, d, ell)
    band = ((d > 0.02) & (d < 0.98)).astype(float)
    if front_only:
        w = np.maximum(front, band)
    else:
        wake = np.clip(d * (1.0 - d), 0.0, 0.25) / 0.25
        w = np.maximum(np.maximum(front, band), max(float(wake_weight), 0.0) * wake)
    return np.clip(w, 0.0, 1.0)



def _write_progress(outdir, status):
    """Atomically write a progress.json snapshot and append a progress.log line.

    Written so an external shell can poll a long run without attaching to the
    process:  ``cat <outdir>/progress.json``  or  ``tail -f <outdir>/progress.log``.
    The JSON is written to a temp file and renamed so a concurrent reader never
    sees a half-written file.  Failures here must never interrupt the solve.
    """
    try:
        os.makedirs(outdir, exist_ok=True)
        jpath = os.path.join(outdir, 'progress.json')
        tmp = jpath + '.tmp'
        with open(tmp, 'w') as f:
            json.dump(status, f, indent=2, default=str)
        os.replace(tmp, jpath)
        with open(os.path.join(outdir, 'progress.log'), 'a') as f:
            f.write(json.dumps(status, default=str) + '\n')
    except Exception:
        pass


def _strict_crack_front_drive_weight(mesh, d, ell, cfg):
    """Return a narrow front-localization weight for degraded crack drive.

    This is deliberately more restrictive than the retained-toughening weight.
    Stored defect energy should lower the crack-front resistance when the crack
    cuts through the process zone; it should not become a broad bulk damage
    source throughout the full process-zone cloud.
    """
    pz = getattr(cfg, 'process_zone', None)
    d = np.asarray(d, dtype=float)
    nn = d.size

    # Compute nodal |grad d| by area-weighted projection.
    gd_acc = np.zeros(nn)
    w_acc = np.zeros(nn)
    for e in range(mesh.ne):
        conn = mesh.elems[e]
        de = d[conn]
        A = mesh.area_e[e]
        grad_d = mesh.dNdx_e[e] @ de
        gdmag = float(np.sqrt(grad_d[0]**2 + grad_d[1]**2))
        gd_acc[conn] += gdmag * A / 3.0
        w_acc[conn] += A / 3.0
    gd_node = gd_acc / np.maximum(w_acc, 1e-30)

    # Use the phase-field gradient to seed the crack front, but suppress the
    # fully damaged wake.  Allow intact material immediately ahead of the front
    # to receive drive through the distance mask below.
    dmax = float(getattr(pz, 'crack_front_dmax', 0.92)) if pz is not None else 0.92
    fthr = float(getattr(pz, 'crack_front_grad_threshold', 0.50)) if pz is not None else 0.50
    front_ind = gd_node * np.maximum(1.0 - d, 0.0) * (d < dmax)
    fmax = float(np.max(front_ind)) if front_ind.size else 0.0
    if fmax < 1e-14:
        return np.zeros(nn)
    seeds = np.where(front_ind > np.clip(fthr, 1e-6, 0.999) * fmax)[0]
    if len(seeds) == 0:
        return np.zeros(nn)

    radius_factor = float(getattr(pz, 'crack_front_radius_factor', 1.25)) if pz is not None else 1.25
    R_front = max(radius_factor * float(ell), 1e-30)
    Xs = mesh.nodes[seeds]
    dist = np.full(nn, np.inf)
    chunk = 2000
    for i in range(0, nn, chunk):
        idx = slice(i, min(i + chunk, nn))
        Xi = mesh.nodes[idx]
        dx = Xi[:, 0:1] - Xs[:, 0][None, :]
        dy = Xi[:, 1:2] - Xs[:, 1][None, :]
        dist[idx] = np.sqrt(np.min(dx*dx + dy*dy, axis=1))

    wdist = np.clip(1.0 - dist / R_front, 0.0, 1.0)
    wdist = wdist*wdist*(3.0 - 2.0*wdist)  # smoothstep

    # Suppress the fully damaged wake.  Do not require d>0; the crack target
    # must be allowed to act just ahead of the front in intact material.
    wake_cut = float(getattr(pz, 'crack_front_wake_dmax', 0.85)) if pz is not None else 0.85
    wake_width = max(float(getattr(pz, 'crack_front_wake_width', 0.10)) if pz is not None else 0.10, 1e-6)
    x = np.clip((d - wake_cut) / wake_width, 0.0, 1.0)
    wake_supp = 1.0 - x*x*(3.0 - 2.0*x)

    return np.clip(wdist * wake_supp, 0.0, 1.0)

def _compute_toughening_energy(q_node, weight, mesh, ell, Gc_eff, cfg):
    """Stored energy of the retained process-zone toughening state [J/m]."""
    q = np.maximum(np.asarray(q_node, dtype=float), 0.0)  # J/m^2
    w = np.clip(np.asarray(weight, dtype=float), 0.0, 1.0)
    coeff = max(float(getattr(cfg.phase_field, 'toughening_storage_coeff', 0.0)), 0.0)
    if coeff <= 0:
        return 0.0
    # energy density scale: q^2/Gc has J/m^2, divide by ell -> J/m^3.
    edens = w * 0.5 * coeff * q*q / max(Gc_eff, 1e-30) / max(ell, 1e-30)
    node_area = getattr(mesh, '_node_area_cache', None)
    if node_area is None:
        node_area = _node_lumped_area(mesh)
        try:
            setattr(mesh, '_node_area_cache', node_area)
        except Exception:
            pass
    return float(np.sum(edens * node_area))


def _update_toughening_state(q_node, dwp_step_node, d, mesh, ell, Gc_eff, cfg, emit_prob_node=None):
    """Update retained plasticity-driven local fracture resistance q_Gc.

    q_node has units J/m^2 and enters Gc_local = Gc_eff + q_node.  It is driven
    by accepted plastic work only near the crack/front process zone.  The update
    is intentionally separate from gross Wp so the energy audit can report a
    partition between immediate plastic dissipation and retained toughening.
    """
    mode = getattr(cfg.phase_field, 'wp_gc_coupling_mode', 'state')
    q_old = np.asarray(q_node, dtype=float).copy()
    q_new = q_old.copy()
    nn = q_old.size
    info = {
        'dq_node': np.zeros(nn),
        'weight': np.zeros(nn),
        'Etough': 0.0,
        'dEtough': 0.0,
        'dDtough': 0.0,
        'qmax': 0.0,
    }
    if mode in ('off', 'direct') or float(getattr(cfg.phase_field, 'plastic_work_to_Gc_efficiency', 0.0)) <= 0.0:
        return q_new, info

    cap = float(getattr(cfg.phase_field, 'Gc_local_cap_factor', np.inf)) * float(Gc_eff)
    q_cap = max(cap - float(Gc_eff), 0.0) if np.isfinite(cap) else np.inf
    info['qmax'] = q_cap if np.isfinite(q_cap) else 0.0

    w = _toughening_weight(
        mesh, d, ell,
        front_only=bool(getattr(cfg.phase_field, 'toughening_front_only', True)),
        wake_weight=float(getattr(cfg.phase_field, 'toughening_wake_weight', 0.25)),
    )
    dwp = np.maximum(np.asarray(dwp_step_node, dtype=float), 0.0)  # J/m^3
    eta = max(float(getattr(cfg.phase_field, 'plastic_work_to_Gc_efficiency', 0.0)), 0.0)

    # q_Gc can be driven by accepted plastic work, by crack-tip emission
    # events, or by a mixture of both.  The emission branch is the preferred
    # physical closure: dislocation nucleation at the tip builds a retained
    # process-zone state.  Plastic-work drive remains available as an ablation.
    pz = getattr(cfg, 'process_zone', None)
    driver = str(getattr(pz, 'qgc_driver', 'plastic_work')).lower() if pz is not None else 'plastic_work'
    dq_plastic = eta * w * dwp * max(ell, 1e-30)  # J/m^2
    if emit_prob_node is None:
        P_emit = np.zeros_like(q_old)
    else:
        P_emit = np.clip(np.asarray(emit_prob_node, dtype=float), 0.0, 1.0)
    qfac = max(float(getattr(pz, 'qgc_from_emission_factor', 0.0)), 0.0) if pz is not None else 0.0
    dq_emit = qfac * w * float(Gc_eff) * P_emit
    if driver == 'emission':
        dq_drive = dq_emit
    elif driver == 'mixed':
        dq_drive = dq_plastic + dq_emit
    else:
        dq_drive = dq_plastic

    # Optional relaxation of retained shielding/toughening, useful for later
    # memory-kernel studies but zero by default.
    relax = np.clip(float(getattr(cfg.phase_field, 'toughening_relax_per_step', 0.0)), 0.0, 1.0)
    if relax > 0:
        q_new = (1.0 - relax) * q_new
    q_new = q_new + dq_drive
    if np.isfinite(q_cap):
        q_new = np.minimum(q_new, q_cap)
    q_new = np.maximum(q_new, 0.0)
    dq = q_new - q_old

    Eold = _compute_toughening_energy(q_old, w, mesh, ell, Gc_eff, cfg)
    Enew = _compute_toughening_energy(q_new, w, mesh, ell, Gc_eff, cfg)
    cdiss = max(float(getattr(cfg.phase_field, 'toughening_dissipation_coeff', 0.0)), 0.0)
    node_area = getattr(mesh, '_node_area_cache', None)
    if node_area is None:
        node_area = _node_lumped_area(mesh)
        try:
            setattr(mesh, '_node_area_cache', node_area)
        except Exception:
            pass
    dDtough = float(np.sum(w * cdiss * np.abs(dq) / max(ell, 1e-30) * node_area))
    info.update({'dq_node': dq, 'weight': w, 'Etough': Enew, 'dEtough': Enew - Eold, 'dDtough': dDtough})
    return q_new, info


def _crack_frontier_mask(mesh, d, notch_nodes, dthr=0.8, layers=2):
    """Nodes within `layers` graph-layers of the broken component connected to
    the notch.  Cleavage (source mode) may only fire here: a first-passage
    event EXTENDS the connected crack; it does not nucleate remote damage.
    Without this constraint, persistent fired memory percolates through the
    diffuse process-zone cloud (speckle failure, load-bearing debris)."""
    adj = getattr(mesh, '_node_adj', None)
    if adj is None:
        adj = [set() for _ in range(mesh.nn)]
        for conn in mesh.elems:
            a, b, c = int(conn[0]), int(conn[1]), int(conn[2])
            adj[a].update((b, c)); adj[b].update((a, c)); adj[c].update((a, b))
        adj = [np.fromiter(s, dtype=np.int64) for s in adj]
        mesh._node_adj = adj
    d = np.asarray(d, dtype=float)
    broken = d >= float(dthr)
    comp = np.zeros(mesh.nn, dtype=bool)
    seeds = [int(n) for n in np.asarray(notch_nodes, dtype=int) if broken[int(n)]]
    if not seeds:
        seeds = [int(n) for n in np.asarray(notch_nodes, dtype=int)]
    stack = list(seeds)
    comp[stack] = True
    while stack:
        n = stack.pop()
        for m in adj[n]:
            if broken[m] and not comp[m]:
                comp[m] = True
                stack.append(int(m))
    out = comp.copy()
    for _ in range(max(int(layers), 0)):
        grow = np.zeros(mesh.nn, dtype=bool)
        idx = np.where(out)[0]
        for n in idx:
            grow[adj[n]] = True
        out |= grow
    return out


def _local_Gc_from_state(Gc_eff, wp_node, q_tough_node, ell, cfg,
                         rho_node=None, material=None, plast_model=None, T=None,
                         return_terms=False, d=None, mesh=None):
    """Return local net fracture resistance according to process-zone state.

    New process-zone form separates three effects:
      Gc_net = Gc0 + q_blunt - G_stored_release(rho)   # G_shield is drive-side only

    q_blunt is the retained blunting/shielding state; G_shield is a
    crack-extension backstress contribution; G_stored_release is the
    embrittling stored defect energy released when the crack cuts through the
    process zone.  With rho_node=None this reduces to the legacy q_tough form.
    """
    mode = getattr(cfg.phase_field, 'wp_gc_coupling_mode', 'state')
    cap = float(getattr(cfg.phase_field, 'Gc_local_cap_factor', np.inf)) * float(Gc_eff)
    base = np.asarray(q_tough_node, dtype=float) if np.ndim(q_tough_node) else np.asarray([float(q_tough_node)])
    nn = base.size

    if mode == 'off' or float(getattr(cfg.phase_field, 'plastic_work_to_Gc_efficiency', 0.0)) <= 0.0:
        Gc_net = float(Gc_eff) * np.ones(nn)
        terms = {'G_shield': np.zeros(nn), 'G_stored_release': np.zeros(nn), 'e_stored': np.zeros(nn),
                 'tau_back_crack': np.zeros(nn), 'Gc_floor': 0.0, 'q_blunt': np.zeros(nn)}
    elif mode == 'direct':
        eta = float(getattr(cfg.phase_field, 'plastic_work_to_Gc_efficiency', 0.0))
        Gc_net = float(Gc_eff) + eta * np.asarray(wp_node, dtype=float) * ell
        terms = {'G_shield': np.zeros_like(Gc_net), 'G_stored_release': np.zeros_like(Gc_net),
                 'e_stored': np.zeros_like(Gc_net), 'tau_back_crack': np.zeros_like(Gc_net),
                 'Gc_floor': 0.0, 'q_blunt': Gc_net - float(Gc_eff)}
    else:
        q = np.asarray(q_tough_node, dtype=float)
        if rho_node is not None and material is not None:
            terms = _process_zone_crack_terms(rho_node, q, ell, Gc_eff, material, cfg,
                                              plast_model=plast_model, T=T, d=d, mesh=mesh)
            # Shielding/backstress is a crack-driving-force reduction, not an additional
            # local fracture-resistance term. Do not double-count G_shield here.
            Gc_net = float(Gc_eff) + q - terms['G_stored_release']
        else:
            terms = {'G_shield': np.zeros_like(q), 'G_stored_release': np.zeros_like(q),
                     'e_stored': np.zeros_like(q), 'tau_back': np.zeros_like(q),
                     'Gc_floor': 0.0, 'q_blunt': q}
            Gc_net = float(Gc_eff) + q

    if np.isfinite(cap):
        Gc_net = np.minimum(Gc_net, cap)
    floor = float(terms.get('Gc_floor', 0.0)) if isinstance(terms, dict) else 0.0
    if floor > 0:
        Gc_net = np.maximum(Gc_net, floor)
    terms = dict(terms)
    terms['Gc_net'] = Gc_net
    # naming consistency for diagnostics
    if 'tau_back' in terms and 'tau_back_crack' not in terms:
        terms['tau_back_crack'] = terms['tau_back']
    if return_terms:
        return Gc_net, terms
    return Gc_net


def _node_to_gp_mean(mesh, node_vals):
    """Element-wise mean of a nodal field."""
    v = np.asarray(node_vals, dtype=float)
    return np.mean(v[mesh.elems], axis=1)

def _process_zone_dislocation_backstress(rho_node, T, plast_model, material, cfg):
    """Stress-like process-zone backstress opposing dislocation emission.

    Preferred model is the Arrhenius-Taylor inversion used for plastic flow,
    evaluated at the local process-zone density and a reference rate.  The
    fallback is alpha*G*b*sqrt(rho).  The returned quantity is a positive stress
    [Pa] and should be subtracted from the crack-tip resolved/emission stress.
    """
    pz = getattr(cfg, 'process_zone', None)
    rho = np.maximum(np.asarray(rho_node, dtype=float), 1e6)
    alpha = max(float(getattr(pz, 'backstress_alpha', 0.35)), 0.0) if pz is not None else 0.35
    tau_sqrt = alpha * material.G * material.b * np.sqrt(rho)

    tau_arr = np.zeros_like(tau_sqrt)
    if pz is not None and str(getattr(pz, 'backstress_model', 'arrhenius_taylor')).lower() in ('arrhenius_taylor', 'max'):
        try:
            eps_ref = max(float(getattr(pz, 'backstress_rate_ref', 1e-4)), 1e-30)
            phi_max = float(getattr(cfg.dislocations, 'phi_plastic_max', 20.0))
            tau_arr = arrhenius_taylor_flow_stress(
                rho, T, eps_ref, plast_model, material.b, phi_plastic_max=phi_max
            )
            tau_arr = np.maximum(np.asarray(tau_arr, dtype=float), 0.0)
        except Exception:
            tau_arr = np.zeros_like(tau_sqrt)

    model = str(getattr(pz, 'backstress_model', 'arrhenius_taylor')).lower() if pz is not None else 'sqrt_taylor'
    if model == 'sqrt_taylor':
        tau = tau_sqrt
    elif model == 'max':
        tau = np.maximum(tau_sqrt, tau_arr)
    else:
        # If the Arrhenius inversion produces zeros because the barrier is not
        # resolvable at the local state, retain the Taylor sqrt floor.
        tau = np.maximum(tau_arr, 0.25 * tau_sqrt)
    scale = max(float(getattr(pz, 'backstress_scale', 1.0)), 0.0) if pz is not None else 1.0
    return scale * tau




def _process_zone_source_availability(rho_node, cfg):
    """Finite source/pile-up availability for crack-tip emission/storage.

    This is not a numerical gate on fracture.  It is a coarse-grained process-zone
    state law: as the local dislocation density approaches the source/pile-up
    storage capacity, additional independent crack-tip emission events become
    less available because the zone already carries Taylor/Peierls backstress and
    nearby depinning events become correlated.  The Hill form is a smooth
    surrogate for the high-density multi-hit depinning hazard discussed in the
    model notes.
    """
    pz = getattr(cfg, 'process_zone', None)
    rho = np.maximum(np.asarray(rho_node, dtype=float), 1e6)
    if pz is None or not bool(getattr(pz, 'source_availability_enabled', True)):
        return np.ones_like(rho), np.full_like(rho, np.inf)

    disl = getattr(cfg, 'dislocations', None)
    rho_cap = float(getattr(disl, 'rho_cap', 1e16)) if disl is not None else 1e16
    rho_sat = float(getattr(pz, 'source_rho_sat', 0.0))
    if not np.isfinite(rho_sat) or rho_sat <= 0:
        frac = float(getattr(pz, 'source_rho_sat_fraction', 0.20))
        rho_sat = max(frac * rho_cap, 1e8)
    rho_sat_arr = np.full_like(rho, rho_sat)

    m = max(float(getattr(pz, 'source_availability_power', 2.0)), 0.1)
    floor = np.clip(float(getattr(pz, 'source_availability_floor', 0.0)), 0.0, 1.0)
    x = np.maximum(rho / np.maximum(rho_sat_arr, 1e6), 0.0)
    avail = floor + (1.0 - floor) / (1.0 + np.power(x, m))
    return np.clip(avail, floor, 1.0), rho_sat_arr


def _process_zone_multihit_state(rho_node, material, cfg):
    """Return correlated-depinning multi-hit count for each node.

    The obstacle count is geometric: an emitted crack-tip dislocation segment
    must clear obstacles over a finite slip/depinning path L_path.  The obstacle
    spacing is taken as L_rho = 1/sqrt(rho_pz).  This avoids introducing a
    phenomenological rho_sat while preserving temperature dependence through the
    elementary Arrhenius barrier used below.
    """
    pz = getattr(cfg, 'process_zone', None)
    rho = np.maximum(np.asarray(rho_node, dtype=float), 1e6)
    spacing = 1.0 / np.sqrt(rho)
    if pz is None or not bool(getattr(pz, 'multihit_enabled', True)):
        n = np.ones_like(rho, dtype=float)
        return n, spacing, np.full_like(rho, np.nan)

    L_nm = float(getattr(pz, 'multihit_path_length_nm', 50.0))
    if np.isfinite(L_nm) and L_nm > 0:
        L_path = L_nm * 1e-9
    else:
        L_path = max(float(getattr(pz, 'multihit_path_length_b', 200.0)), 1.0) * material.b
    power = max(float(getattr(pz, 'multihit_density_power', 1.0)), 0.1)
    nmin = max(int(getattr(pz, 'multihit_min_hits', 1)), 1)
    nmax = max(int(getattr(pz, 'multihit_max_hits', 12)), nmin)
    obstacle_count = np.power(np.maximum(L_path / np.maximum(spacing, 1e-30), 0.0), power)
    n = np.ceil(np.maximum(float(nmin), obstacle_count))
    n = np.clip(n, nmin, nmax)
    return n.astype(float), spacing, np.full_like(rho, L_path)


def _process_zone_multihit_log_rate(log_eta0, G_eV, T, rho_node, material, cfg, target):
    """Apply correlated multi-hit Arrhenius/Taylor depinning to a hazard rate.

    For a single elementary event, log(lambda)=log(eta0)-G/kT.  If n correlated
    hits are required, the thermodynamic activation cost is n*G.  This is the
    low-cost implementation of the multi-hit Arrhenius-Taylor idea: high rho
    increases the number of correlated obstacles rather than imposing an
    arbitrary density saturation or non-thermodynamic gate.
    """
    pz = getattr(cfg, 'process_zone', None)
    G = np.maximum(np.asarray(G_eV, dtype=float), 0.0)
    nn = G.size
    out = {
        'n_hits': np.ones(nn),
        'spacing_m': np.full(nn, np.nan),
        'path_m': np.full(nn, np.nan),
        'log_suppression': np.zeros(nn),
    }
    base = np.asarray(log_eta0 - G / max(KB * T, 1e-30), dtype=float)
    if pz is None or not bool(getattr(pz, 'multihit_enabled', True)):
        return base, out
    apply_to = str(getattr(pz, 'multihit_apply_to', 'both')).lower()
    if apply_to == 'off' or (apply_to != 'both' and apply_to != str(target).lower()):
        return base, out

    n_hits, spacing, path = _process_zone_multihit_state(rho_node, material, cfg)
    kT = max(KB * T, 1e-30)
    log_rate = log_eta0 - n_hits * G / kT
    out.update({
        'n_hits': n_hits,
        'spacing_m': spacing,
        'path_m': path,
        'log_suppression': -np.maximum(n_hits - 1.0, 0.0) * G / kT,
    })
    return np.asarray(log_rate, dtype=float), out


def _process_zone_crack_terms(rho_node, q_tough_node, ell, Gc_eff, material, cfg,
                              plast_model=None, T=None, d=None, mesh=None):
    """Return crack-shielding and stored-energy release terms [J/m^2].

    The process zone affects crack extension in two opposing ways:
      1. shielding/backstress/blunting increases apparent resistance,
      2. stored defect energy decreases net resistance when the crack cuts
         through the dislocation-rich process zone.

    This function returns local nodal arrays so that phase-field crack growth
    sees Gc_net = Gc0 + q_tough - G_stored_release; G_shield enters only G_eff.
    """
    pz = getattr(cfg, 'process_zone', None)
    rho = np.maximum(np.asarray(rho_node, dtype=float), 1e6)
    nn = rho.size
    q = np.asarray(q_tough_node, dtype=float) if q_tough_node is not None else np.zeros(nn)

    Eprime = material.E / max(1.0 - material.nu**2, 1e-12)
    tau_back = _process_zone_dislocation_backstress(
        rho, float(T) if T is not None else 300.0,
        plast_model if plast_model is not None else None,
        material, cfg
    ) if plast_model is not None and T is not None else (
        max(float(getattr(pz, 'backstress_alpha', 0.35)), 0.0) * material.G * material.b * np.sqrt(rho)
        if pz is not None else np.zeros(nn)
    )

    if pz is not None and bool(getattr(pz, 'crack_shielding_enabled', True)):
        cshield = max(float(getattr(pz, 'crack_shielding_coeff', 1.0)), 0.0)
        G_shield = cshield * (tau_back**2 / max(Eprime, 1e-30)) * max(float(ell), 1e-30)
    else:
        G_shield = np.zeros(nn)

    if pz is not None and bool(getattr(pz, 'stored_energy_enabled', True)):
        cstore = max(float(getattr(pz, 'stored_energy_coeff', 0.5)), 0.0)
        chi = max(float(getattr(pz, 'stored_energy_release_efficiency', 0.25)), 0.0)
        # Dislocation line energy density ~ C*G*b^2*rho*log(R/r0).  Use the
        # local spacing as R and b as the core scale, bounded to avoid singulars.
        logarg = np.maximum(1.0 / np.maximum(material.b * np.sqrt(rho), 1e-30), 1.0)
        line_log = np.clip(np.log(logarg), 1.0, 20.0)
        e_stored = cstore * material.G * material.b**2 * rho * line_log  # J/m^3
        G_available = e_stored * max(float(ell), 1e-30)
        # Release is bounded by available stored defect energy through chi.
        # The old cap G_release <= factor*Gc0 is retained only as an explicit
        # optional ablation; the default cap factor is now infinity.
        G_release = chi * G_available
        capf = float(getattr(pz, 'stored_energy_release_cap_factor', np.inf))
        if np.isfinite(capf) and capf > 0:
            G_release = np.minimum(G_release, capf * float(Gc_eff))
    else:
        e_stored = np.zeros(nn)
        G_release = np.zeros(nn)

    Gc_floor = max(float(getattr(pz, 'Gc_net_floor_factor', 0.05)), 0.0) * float(Gc_eff) if pz is not None else 0.0

    # Front-localize the stored-energy embrittlement.  The bulk release is kept
    # for diagnostics, but the value that actually lowers Gc_net (and feeds the
    # crack drive) is restricted to the connected, advancing crack front via the
    # same strict front weight used for the drive.  This prevents a broad
    # high-rho cloud from collapsing the toughness everywhere and producing a
    # diffuse damage field instead of a crack.
    G_release_bulk = np.maximum(G_release, 0.0)
    front_masked = bool(getattr(pz, 'crack_stored_release_front_masked', True)) if pz is not None else False
    if front_masked and (d is not None) and (mesh is not None):
        fw = np.clip(_strict_crack_front_drive_weight(mesh, d, ell, cfg), 0.0, 1.0)
        G_release_eff = fw * G_release_bulk
    else:
        G_release_eff = G_release_bulk

    return {
        'tau_back': tau_back,
        'G_shield': np.maximum(G_shield, 0.0),
        'e_stored': np.maximum(e_stored, 0.0),
        'G_stored_available': np.maximum((e_stored * max(float(ell), 1e-30)), 0.0),
        'G_stored_release': np.maximum(G_release_eff, 0.0),
        'G_stored_release_bulk': G_release_bulk,
        'Gc_floor': Gc_floor,
        'q_blunt': q,
    }




def _stash_stress_diags(hist, step, mesh, s1_node=None, crack_hazard_info=None, B_crack_node=None):
    """Store elastic/hazard diagnostic fields for the snapshot plots."""
    try:
        nn = mesh.nn
        if s1_node is not None:
            hist.sig1_fields[step] = np.asarray(s1_node, dtype=float).copy()
        if crack_hazard_info is not None:
            st = np.asarray(crack_hazard_info.get('sigma_tip_crack', np.zeros(nn)), dtype=float)
            fw = np.asarray(crack_hazard_info.get('front_weight', np.zeros(nn)), dtype=float)
            hist.sigtip_fields[step] = st.copy()
            hist.fw_fields[step] = fw.copy()
        if B_crack_node is not None:
            hist.B_fields[step] = np.asarray(B_crack_node, dtype=float).copy()
    except Exception:
        pass

def _compute_process_zone_crack_hazard(Hhist_node, Gc_net_node, pz_crack_info, d, mesh, ell, T, dt, cfg,
                                       sigma_tip_node=None, B_crack_node=None, B_emit_node=None):
    """Crack-growth hazard based on either G_eff/Gc or resolved crack-tip stress.

    Two modes are supported:
      * ratio (legacy): phenomenological barrier collapse with R=G_eff/Gc_net.
      * resolved_stress: paper-aligned hazard using the fracture barrier
        G*_f(sigma_tip,T) evaluated on the resolved/amplified near-tip stress.

    If ``crack_first_passage`` is enabled, this routine carries an accumulated
    action field B_crack = integral(lambda dt) and accepts crack advance when
    B/B_target approaches unity.  This avoids re-rolling independent crack
    probabilities each load step and is closer to a renewal/first-passage law.
    """
    pz = getattr(cfg, 'process_zone', None)
    H = np.maximum(np.asarray(Hhist_node, dtype=float), 0.0)
    nn = H.size
    z = np.zeros(nn)
    out = {
        'P_crack': z.copy(), 'G_app': z.copy(), 'G_eff': z.copy(),
        'barrier_eV': np.full(nn, np.inf), 'hazard': z.copy(),
        'hazard_raw': z.copy(), 'front_weight': z.copy(),
        'B_crack': z.copy(), 'sigma_tip_crack': z.copy(),
    }
    if pz is None or not bool(getattr(pz, 'enabled', True)) or not bool(getattr(pz, 'crack_hazard_enabled', False)):
        return out
    if T <= 0:
        return out

    Gc_arg = np.asarray(pz_crack_info.get('Gc_net', Gc_net_node), dtype=float)
    Gc = np.maximum(np.broadcast_to(Gc_arg, (nn,)), 1e-30)
    G_app = 2.0 * max(ell, 1e-30) * H
    G_shield = np.asarray(pz_crack_info.get('G_shield', z), dtype=float)
    G_release = np.asarray(pz_crack_info.get('G_stored_release', z), dtype=float)
    G_eff = np.maximum(G_app - G_shield + G_release, 0.0)

    # --- Crack-tip dislocation EMISSION competition (Rice-Thomson) ---------
    # An emission first-passage hazard, driven by the SAME stress-concentrated
    # tip stress as cleavage but using the native nanopillar (exp_floor) barrier,
    # injects emitted-dislocation density at the front.  That density shields the
    # tip via the standard elastic backstress (G_shield_em = cshield*tau_back^2/
    # E'*ell, tau_back = alpha*G*b*sqrt(rho_emit)), lowering G_eff and hence the
    # cleavage drive.  Emission wins -> tip shields -> ductile; cleavage wins ->
    # brittle.  The crossover temperature is the predicted DBTT.
    B_emit_new = np.zeros(nn)
    rho_emit = np.zeros(nn)
    G_shield_em = np.zeros(nn)
    emission_model = getattr(cfg, 'emission_model', None)
    if (bool(getattr(pz, 'crack_emission_enabled', False)) and emission_model is not None
            and sigma_tip_node is not None):
        # de-smear the tip stress exactly as the cleavage drive does, including
        # FRONT LOCALIZATION (the sub-grid K-field exists only at the front)
        sig_em = np.maximum(np.broadcast_to(np.asarray(sigma_tip_node, dtype=float), (nn,)), 0.0)
        sig_em = sig_em * max(float(getattr(pz, 'crack_resolved_stress_scale', 1.0)), 0.0)
        fw_em = np.clip(_strict_crack_front_drive_weight(mesh, d, ell, cfg), 0.0, 1.0)
        r_pz_e = float(getattr(pz, 'crack_process_zone_r_pz_m', 0.0) or 0.0)
        if r_pz_e > 0.0:
            chi_e = float(np.clip(np.sqrt(max(ell, 1e-30) / r_pz_e), 1.0,
                                  max(float(getattr(pz, 'crack_desmear_max', 1e3)), 1.0)))
            sig_em = sig_em * (1.0 + (chi_e - 1.0) * fw_em)
            scap = float(getattr(pz, 'crack_sigma_cap_Pa', 0.0) or 0.0)
            if scap > 0.0:
                sig_em = np.minimum(sig_em, scap)
        G_emit = emission_model.G_barrier(sig_em, T)
        eta0_em = max(float(getattr(pz, 'crack_emission_eta0', 1e12)), 1e-300)
        lr_em = np.clip(np.log(eta0_em) - G_emit / max(KB * T, 1e-30), -745.0, 80.0)
        Hinc_em = np.clip(np.exp(lr_em) * fw_em * max(dt, 0.0), 0.0, 80.0)
        Bprev_em = np.zeros(nn) if B_emit_node is None else np.broadcast_to(
            np.asarray(B_emit_node, dtype=float), (nn,)).copy()
        reset_d = float(getattr(pz, 'crack_B_reset_damage', 0.95))
        Bprev_em = np.where(np.asarray(d, dtype=float) >= reset_d, 0.0, Bprev_em)
        tau_rel_em = float(getattr(pz, 'crack_B_relax_time_s', 0.0) or 0.0)
        if tau_rel_em > 0.0 and dt > 0.0:
            Bprev_em = Bprev_em * float(np.exp(-min(dt / tau_rel_em, 80.0)))
        Bcap = max(float(getattr(pz, 'crack_B_cap', 5.0)), 1.0)
        B_emit_new = np.clip(Bprev_em + Hinc_em, 0.0, Bcap)
        Btar_em = max(float(getattr(pz, 'crack_emission_B_target', 1.0)), 1e-12)
        rho_max = max(float(getattr(pz, 'crack_emission_rho_max', 1e15)), 0.0)
        # saturating emitted density, localized to the advancing front
        rho_emit = rho_max * np.clip(B_emit_new / Btar_em, 0.0, 1.0) * fw_em
        # elastic shielding from the emitted dislocations (same path as bulk PZ)
        alpha_em = max(float(getattr(pz, 'crack_emission_backstress_alpha', 0.35)), 0.0)
        cshield = max(float(getattr(pz, 'crack_shielding_coeff', 1.0)), 0.0)
        Eprime = cfg.material.E / max(1.0 - cfg.material.nu**2, 1e-12)
        tau_back_em = alpha_em * cfg.material.G * cfg.material.b * np.sqrt(np.maximum(rho_emit, 0.0))
        G_shield_em = cshield * (tau_back_em**2 / max(Eprime, 1e-30)) * max(float(ell), 1e-30)
        # emission shielding reduces the crack driving force (blunting)
        G_eff = np.maximum(G_eff - G_shield_em, 0.0)

    H_eff_drive = G_eff / (2.0 * max(ell, 1e-30))

    fw = _strict_crack_front_drive_weight(mesh, d, ell, cfg)
    drive_scale = max(float(getattr(pz, 'crack_drive_scale', 1.0)), 0.0)
    R = drive_scale * G_eff / Gc

    eta0 = max(float(getattr(pz, 'crack_eta0', 1e12)), 1e-300)
    model = str(getattr(pz, 'crack_hazard_model', 'arrhenius')).lower()
    drive_mode = str(getattr(pz, 'crack_hazard_drive', 'ratio')).lower()

    if drive_mode == 'resolved_stress':
        sig_arg = z if sigma_tip_node is None else np.asarray(sigma_tip_node, dtype=float)
        sigma_tip = np.maximum(np.broadcast_to(sig_arg, (nn,)), 0.0)
        sigma_tip = sigma_tip * max(float(getattr(pz, 'crack_resolved_stress_scale', 1.0)), 0.0)
        # De-smear the AT2-regularized tip stress to the physical process-zone
        # scale.  The FEM resolves the near-tip field smoothed over ell, so the
        # peak resolved stress is ~ K/sqrt(2*pi*ell); the stress that actually
        # drives the Arrhenius fracture instability is the field at the physical
        # process-zone radius r_pz, ~ K/sqrt(2*pi*r_pz).  The intensification is
        #     chi = sqrt(ell / r_pz),
        # which RECOVERS the un-regularized near-tip stress and is mesh-OBJECTIVE
        # (chi * sigma_res = K/sqrt(2*pi*r_pz) is independent of ell).  r_pz sets
        # the stress at which the tip barrier collapses, i.e. the emergent brittle
        # toughness.  Disabled (chi=1) when r_pz<=0.
        r_pz = float(getattr(pz, 'crack_process_zone_r_pz_m', 0.0) or 0.0)
        if r_pz > 0.0:
            chi = np.sqrt(max(ell, 1e-30) / r_pz)
            chi = float(np.clip(chi, 1.0, max(float(getattr(pz, 'crack_desmear_max', 1e3)), 1.0)))
            # FRONT-LOCALIZED de-smearing.  The intensification recovers the
            # sub-grid singular K-field, which exists ONLY at the crack front; a
            # generic bulk point has no hidden singularity to recover.  Applying
            # chi globally turns moderate bulk stresses (~50-80 MPa) into fake
            # GPa tip stresses once chi ~ 50, firing the hazard over the whole
            # stress lobe (the wide damage tongue / remote Gc blobs).  Blend with
            # the strict front weight so chi_node -> chi at the front and -> 1
            # in the bulk:
            fw_chi = np.clip(_strict_crack_front_drive_weight(mesh, d, ell, cfg), 0.0, 1.0)
            chi_node = 1.0 + (chi - 1.0) * fw_chi
            sigma_tip = sigma_tip * chi_node
            # Physical ceiling: the near-tip stress cannot exceed the cohesive /
            # spinodal scale (the material fails first).  Optional; off when <=0.
            s_cap = float(getattr(pz, 'crack_sigma_cap_Pa', 0.0) or 0.0)
            if s_cap > 0.0:
                sigma_tip = np.minimum(sigma_tip, s_cap)
        # First-passage crack barrier from the fracture model itself.  This is
        # the FEM-resolved version of lambda = Gamma0 exp[-G*_f(sigma_tip,T)/kT].
        barrier = cfg.fracture_barrier.G_barrier(sigma_tip, T, cfg.material.b)
        # Optional Eyring bias retains detailed-balance-like low-drive behavior
        # without replacing the physical stress barrier.
        log_rate = np.log(eta0) - barrier / max(KB*T, 1e-30)
        if model == 'eyring':
            beta = max(float(getattr(pz, 'crack_eyring_volume_factor', 1.0)), 0.0)
            x = np.clip(beta * np.maximum(R, 0.0), 0.0, 80.0)
            log_rate = log_rate + np.log(np.maximum(2.0*np.sinh(x), 1e-300))
    else:
        sigma_tip = z.copy()
        H0 = max(float(getattr(pz, 'crack_H0_eV', 1.0)), 0.0) * EV_TO_J
        pexp = max(float(getattr(pz, 'crack_drive_exponent', 2.0)), 0.1)
        barrier = H0 * np.power(np.maximum(1.0 - R, 0.0), pexp)
        log_rate = np.log(eta0) - barrier / max(KB*T, 1e-30)
        if model == 'eyring':
            beta = max(float(getattr(pz, 'crack_eyring_volume_factor', 1.0)), 0.0)
            x = np.clip(beta * R, 0.0, 80.0)
            log_rate = log_rate + np.log(np.maximum(2.0*np.sinh(x), 1e-300))

    log_rate = np.clip(log_rate, -745.0, 80.0)
    # --- Correlated (multi-hit) cleavage renewal -------------------------
    # m>1: front advance requires m cooperative events within one renewal
    # window tau_c.  lambda_eff = gammainc(m, lambda*tau_c)/tau_c.  For
    # m=1, gammainc(1,x)=1-exp(-x) -> lambda in the x<<1 limit (exact
    # Poisson recovery); for any m the rate is bounded by 1/tau_c (renewal
    # bound).  Sub-threshold flank rates are suppressed ~ (lam*tau)^m/m!.
    m_hits = max(float(getattr(pz, 'crack_multihit_m', 1.0)), 1.0)
    if m_hits > 1.0 + 1e-12:
        from scipy.special import gammainc
        tau_c = float(getattr(pz, 'crack_multihit_tau_s', 1e-9))
        # tau_c <= 0: ADAPTIVE renewal window = the load-step time dt.  Firing
        # then requires lambda*dt ~ m (nearly the same threshold as the m=1
        # first-passage clock, so Kc shifts only slightly) while sub-critical
        # flank rates are still suppressed combinatorially ~ (lambda*dt)^m/m!.
        # Use this when a constant-entropy shelf caps the saturated tip rate
        # (lambda_max = eta0*exp(-G*_min/kT)) below 1/tau_c for any fixed
        # atomic-scale tau_c, which would otherwise make m hits uncompletable.
        if tau_c <= 0.0:
            tau_c = max(float(dt), 1e-30)
        tau_c = max(tau_c, 1e-30)
        lam = np.exp(log_rate)
        x = np.clip(lam * tau_c, 0.0, 1e12)
        rate_eff = gammainc(m_hits, x) / tau_c
        log_rate = np.log(np.maximum(rate_eff, 1e-300))
        log_rate = np.clip(log_rate, -745.0, 80.0)
    rate = np.exp(log_rate) * np.clip(fw, 0.0, 1.0)
    Hinc_raw = np.clip(rate * max(dt, 0.0), 0.0, 80.0)

    if bool(getattr(pz, 'crack_first_passage', False)):
        Bprev = np.zeros(nn) if B_crack_node is None else np.asarray(B_crack_node, dtype=float)
        Bprev = np.broadcast_to(Bprev, (nn,)).copy()
        # Do not keep integrating action in fully damaged wake; that action has
        # already been consumed by crack advance.
        reset_d = float(getattr(pz, 'crack_B_reset_damage', 0.95))
        Bprev = np.where(np.asarray(d, dtype=float) >= reset_d, 0.0, Bprev)
        # Sub-critical annealing: dB/dt = lambda - B/tau_relax.  Action that
        # never reaches threshold decays instead of ratcheting, so flank
        # nodes sitting at moderate stress for many steps do not eventually
        # fire by pure accumulation.
        tau_rel = float(getattr(pz, 'crack_B_relax_time_s', 0.0) or 0.0)
        if tau_rel > 0.0 and dt > 0.0:
            Bprev = Bprev * float(np.exp(-min(dt / tau_rel, 80.0)))
        Bcap = max(float(getattr(pz, 'crack_B_cap', 5.0)), 1.0)
        Bnew = np.clip(Bprev + Hinc_raw, 0.0, Bcap)
        Btarget = max(float(getattr(pz, 'crack_B_target', 1.0)), 1e-12)
        # Smooth first-passage acceptance: B/Btarget approaches deterministic
        # crack advance at unity while still resolving subcritical accumulation.
        Hhaz_raw = Bnew.copy()
        P = np.clip(Bnew / Btarget, 0.0, 1.0)
    else:
        Bnew = np.zeros(nn)
        Hhaz_raw = Hinc_raw
        Pcap = float(getattr(pz, 'crack_probability_cap', 1.0))
        if np.isfinite(Pcap) and 0 < Pcap < 1.0:
            Hcap = -np.log(max(1.0 - Pcap, 1e-300))
            Hhaz = np.minimum(Hhaz_raw, Hcap)
        else:
            Hhaz = Hhaz_raw
        P = 1.0 - np.exp(-Hhaz)

    if bool(getattr(pz, 'crack_first_passage', False)):
        Hhaz = np.minimum(Hhaz_raw, max(float(getattr(pz, 'crack_B_cap', 5.0)), 1.0))

    out.update({
        'P_crack': np.clip(P, 0.0, 1.0),
        'G_app': G_app, 'G_eff': G_eff, 'Gc_net': Gc,
        'R_crack': R, 'H_eff_drive': H_eff_drive,
        'barrier_eV': barrier / EV_TO_J,
        'hazard': Hhaz, 'hazard_raw': Hinc_raw,
        'front_weight': fw,
        'front_mask_frac': np.mean(fw > 0.05),
        'H_eff_masked': H_eff_drive * np.clip(fw, 0.0, 1.0),
        'B_crack': Bnew,
        'sigma_tip_crack': sigma_tip,
        'B_emit': B_emit_new,
        'rho_emit': rho_emit,
        'G_shield_emit': G_shield_em,
    })
    return out

def _compute_process_zone_emission(sig_pos_node, rho_node, shield_state, T, dt,
                                   plast_model, material, cfg, d=None, mesh=None, ell=None):
    """Arrhenius crack-tip dislocation-emission probability.

    The effective tip stress is shielded by a Taylor/process-zone back stress
    and by the retained shielding memory.  The same plasticity barrier family is
    used for the emission barrier, avoiding an independent arbitrary crack-tip
    activation law.
    """
    pz = getattr(cfg, 'process_zone', None)
    sigma_tip = np.maximum(np.asarray(sig_pos_node, dtype=float), 0.0)
    nn = sigma_tip.size
    out = {
        'P_emit': np.zeros(nn),
        'sigma_tip_raw': sigma_tip,
        'sigma_tip_eff': np.zeros(nn),
        'sigma_back': np.zeros(nn),
        'emission_rate': np.zeros(nn),
        'emission_hazard': np.zeros(nn),
        'front_weight': np.ones(nn),
    }
    if pz is None or not bool(getattr(pz, 'enabled', True)) or not bool(getattr(pz, 'emission_enabled', True)):
        return out
    if T <= 0:
        return out

    rho = np.maximum(np.asarray(rho_node, dtype=float), 1e6)
    # Process-zone backstress for dislocation emission.  This is stress-like and
    # should follow the process-zone Taylor/Peierls resistance, not the crack
    # fracture-energy bookkeeping.  Prefer Arrhenius-Taylor inversion; fall back
    # to alpha*G*b*sqrt(rho).
    sigma_back_disl = _process_zone_dislocation_backstress(rho, T, plast_model, material, cfg)
    z = np.clip(np.asarray(shield_state, dtype=float), 0.0, 1.0)
    sigma_back_mem = max(float(getattr(pz, 'memory_backstress_factor', 0.0)), 0.0) * z * sigma_tip
    sigma_back = sigma_back_disl + sigma_back_mem
    min_frac = max(float(getattr(pz, 'min_effective_stress_frac', 0.0)), 0.0)
    sigma_eff = np.maximum(sigma_tip - sigma_back, min_frac * sigma_tip)

    # Localize emission to the crack/process-zone front unless disabled.
    if mesh is not None and d is not None and ell is not None:
        fw = _toughening_weight(mesh, d, ell, front_only=True, wake_weight=0.0)
    else:
        fw = np.ones(nn)

    # Smooth source/pile-up availability suppresses the independent emission
    # clock before any emitted dislocations are committed to rho or plastic work.
    # This is the low-cost alternative to a full correlated multi-hit depinning
    # hazard at high process-zone density.
    source_avail, rho_sat_local = _process_zone_source_availability(rho, cfg)

    stress_scale = max(float(getattr(pz, 'emission_H_scale', 1.0)), 0.0)
    sigma_barrier = stress_scale * sigma_eff
    try:
        G_emit = plast_model.G_barrier(sigma_barrier, T)
    except Exception:
        # Fall back to rational H-sigma*v if a custom model lacks G_barrier.
        try:
            G_emit = np.maximum(plast_model.H(sigma_barrier) - sigma_barrier * plast_model.v(sigma_barrier, T), 0.0)
        except Exception:
            G_emit = np.full(nn, np.inf)
    eta0 = max(float(getattr(pz, 'emission_eta0', 1e13)), 1e-300)
    log_rate, mh_emit = _process_zone_multihit_log_rate(
        np.log(eta0), np.asarray(G_emit, dtype=float), T, rho, material, cfg, target='emission'
    )
    log_rate = np.clip(log_rate, -745.0, 80.0)
    # Legacy source_availability is disabled by default.  If explicitly enabled,
    # it remains an ablation multiplier; the preferred high-rho physics is the
    # multi-hit Arrhenius/Taylor barrier above.
    rate = np.exp(log_rate) * np.clip(fw, 0.0, 1.0) * np.clip(source_avail, 0.0, 1.0)
    H_raw = np.clip(rate * max(dt, 0.0), 0.0, 80.0)
    Pcap = float(getattr(pz, 'emission_probability_cap', 1.0))
    if np.isfinite(Pcap) and Pcap > 0 and Pcap < 1.0:
        # Convert probability cap to a hazard cap.  This keeps the event-clock
        # increment itself small instead of letting H=80 and only clipping the
        # final probability.  A future version should subcycle the full PZ
        # kinetics, but this prevents saturated emission from masquerading as
        # resolved Arrhenius competition.
        Hcap = -np.log(max(1.0 - Pcap, 1e-300))
        H = np.minimum(H_raw, Hcap)
    else:
        H = H_raw
    P = 1.0 - np.exp(-H)
    out.update({'P_emit': np.clip(P, 0.0, 1.0), 'sigma_tip_eff': sigma_eff,
                'sigma_back': sigma_back, 'sigma_back_disl': sigma_back_disl,
                'sigma_back_mem': sigma_back_mem, 'emission_rate': rate, 'emission_hazard': H,
                'emission_hazard_raw': H_raw, 'front_weight': fw,
                'source_availability': source_avail, 'rho_source_sat': rho_sat_local,
                'multihit_n_emit': mh_emit.get('n_hits', np.ones(nn)),
                'multihit_spacing_emit': mh_emit.get('spacing_m', np.full(nn, np.nan)),
                'multihit_path_emit': mh_emit.get('path_m', np.full(nn, np.nan)),
                'multihit_log_suppression_emit': mh_emit.get('log_suppression', np.zeros(nn))})
    return out



def _compute_process_zone_mobility_partition(P_emit_node, sigma_tip_eff_node, sigma_back_node,
                                             rho_node, T, dt, plast_model, material, cfg):
    """Split emitted crack-tip dislocations into mobile/escaped/stored fractions.

    This is the missing physical separation between dislocation nucleation at
    the crack tip and subsequent mobility through the process zone.  Emission is
    controlled by the tip nucleation barrier; mobility/escape is controlled by a
    second Arrhenius hazard based on the residual driving stress after the
    process-zone Taylor/Peierls backstress.  The immobile/retained fraction is
    stored as rho_pz/pile-up and contributes to shielding and stored energy.
    """
    pz = getattr(cfg, 'process_zone', None)
    P_emit = np.clip(np.asarray(P_emit_node, dtype=float), 0.0, 1.0)
    nn = P_emit.size
    z = np.zeros(nn)
    out = {
        'P_mobile': z.copy(), 'P_escape': z.copy(), 'P_store': P_emit.copy(),
        'P_mobility': z.copy(), 'mobility_hazard': z.copy(), 'mobility_hazard_raw': z.copy(),
        'sigma_mobility_eff': z.copy(), 'storage_fraction': np.ones(nn),
    }
    if pz is None or not bool(getattr(pz, 'enabled', True)) or not bool(getattr(pz, 'mobility_enabled', True)):
        return out
    if T <= 0 or P_emit.size == 0:
        return out

    sigma_eff_emit = np.maximum(np.asarray(sigma_tip_eff_node, dtype=float), 0.0)
    sigma_back = np.maximum(np.asarray(sigma_back_node, dtype=float), 0.0)
    # Mobility should be harder when the process zone already carries a strong
    # Taylor/Peierls backstress.  This lets immobile emitted dislocations pile up
    # and shield the tip before becoming plastic strain/escape.
    fac_back = max(float(getattr(pz, 'mobility_backstress_factor', 1.0)), 0.0)
    sigma_mob = np.maximum(sigma_eff_emit - fac_back * sigma_back, 0.0)

    stress_scale = max(float(getattr(pz, 'mobility_H_scale', 1.0)), 0.0)
    try:
        G_mob = plast_model.G_barrier(stress_scale * sigma_mob, T)
    except Exception:
        try:
            G_mob = np.maximum(plast_model.H(stress_scale * sigma_mob) - stress_scale * sigma_mob * plast_model.v(stress_scale * sigma_mob, T), 0.0)
        except Exception:
            G_mob = np.full(nn, np.inf)

    eta0 = max(float(getattr(pz, 'mobility_eta0', 1e12)), 1e-300)
    log_rate, mh_mob = _process_zone_multihit_log_rate(
        np.log(eta0), np.asarray(G_mob, dtype=float), T, rho_node, material, cfg, target='mobility'
    )
    log_rate = np.clip(log_rate, -745.0, 80.0)
    rate = np.exp(log_rate)
    H_raw = np.clip(rate * max(dt, 0.0), 0.0, 80.0)
    Pcap = float(getattr(pz, 'mobility_probability_cap', 0.10))
    if np.isfinite(Pcap) and Pcap > 0 and Pcap < 1.0:
        Hcap = -np.log(max(1.0 - Pcap, 1e-300))
        H = np.minimum(H_raw, Hcap)
    else:
        H = H_raw
    P_mob = np.clip(1.0 - np.exp(-H), 0.0, 1.0)

    escape_frac = np.clip(float(getattr(pz, 'mobility_escape_fraction', 0.5)), 0.0, 1.0)
    storage_min = np.clip(float(getattr(pz, 'storage_min_fraction', 0.02)), 0.0, 1.0)
    storage_max = np.clip(float(getattr(pz, 'storage_max_fraction', 0.98)), storage_min, 1.0)
    # Base storage is immobile fraction.  A backstress boost makes high-PZ-density
    # regions retain more of the emitted population as pile-up, but finite
    # source/pile-up capacity must also suppress *additional* storage as rho
    # approaches the process-zone saturation density.
    rho = np.maximum(np.asarray(rho_node, dtype=float), 1e6)
    rho_ref = max(float(getattr(cfg.dislocations, 'rho_cap', 1e16)), 1e6)
    boost = max(float(getattr(pz, 'storage_backstress_boost', 1.0)), 0.0) * np.clip(np.sqrt(rho / rho_ref), 0.0, 1.0)
    f_store_raw = np.clip((1.0 - P_mob) + boost * P_mob, storage_min, storage_max)
    storage_capacity, rho_sat_local = _process_zone_source_availability(rho, cfg)
    f_store = np.clip(f_store_raw * np.clip(storage_capacity, 0.0, 1.0), 0.0, storage_max)
    P_store = np.clip(P_emit * f_store, 0.0, 1.0)
    # Events that are emitted but cannot be stored because the local pile-up is
    # saturated are treated as mobile/escaped candidates rather than stored rho.
    P_mobile_total = np.clip(P_emit * (1.0 - f_store), 0.0, 1.0)
    P_escape = np.clip(P_mobile_total * escape_frac, 0.0, 1.0)
    P_mobile = np.clip(P_mobile_total * (1.0 - escape_frac), 0.0, 1.0)

    out.update({
        'P_mobile': P_mobile, 'P_escape': P_escape, 'P_store': P_store,
        'P_mobility': P_mob, 'mobility_hazard': H, 'mobility_hazard_raw': H_raw,
        'sigma_mobility_eff': sigma_mob, 'storage_fraction': f_store,
        'storage_fraction_raw': f_store_raw, 'storage_capacity': storage_capacity,
        'rho_source_sat': rho_sat_local,
        'multihit_n_mobility': mh_mob.get('n_hits', np.ones(nn)),
        'multihit_spacing_mobility': mh_mob.get('spacing_m', np.full(nn, np.nan)),
        'multihit_path_mobility': mh_mob.get('path_m', np.full(nn, np.nan)),
        'multihit_log_suppression_mobility': mh_mob.get('log_suppression', np.zeros(nn)),
    })
    return out

def _apply_process_zone_rho_terms(rho_gp, dot_ep_gp, P_emit_node, d, mesh, ell, T, dt, material, cfg, P_storage_node=None):
    """Add emission-driven storage and Arrhenius recovery to rho.

    This is the physical stabilizer that should replace pure cap tuning.  Tip
    emission locally adds process-zone dislocations; recovery/escape removes
    density with an Arrhenius rate and a plastic-flow-assisted term.
    """
    pz = getattr(cfg, 'process_zone', None)
    if pz is None or not bool(getattr(pz, 'enabled', True)) or bool(getattr(cfg.dislocations, 'freeze_rho', False)):
        return np.asarray(rho_gp, dtype=float), {
            'drho_emit_gp': np.zeros_like(rho_gp), 'drho_rec_gp': np.zeros_like(rho_gp),
            'rho_recovery_rate_gp': np.zeros_like(rho_gp), 'emission_gp': np.zeros_like(rho_gp),
        }
    rho = np.asarray(rho_gp, dtype=float).copy()
    P_emit_gp = _node_to_gp_mean(mesh, P_emit_node)
    P_storage_gp = _node_to_gp_mean(mesh, P_storage_node) if P_storage_node is not None else P_emit_gp
    if bool(getattr(pz, 'rho_emission_front_only', True)):
        fw_gp = _node_to_gp_mean(mesh, _toughening_weight(mesh, d, ell, front_only=True, wake_weight=0.0))
    else:
        fw_gp = np.ones_like(rho)
    # Storage is driven by the retained/immobile part of the emitted population,
    # not by the total nucleation probability.  Mobile/escaped dislocations do
    # not immediately become stored rho_pz.  Finite source/pile-up availability
    # additionally suppresses storage as the local process zone approaches its
    # saturation density.
    source_cap_node, rho_sat_node = _process_zone_source_availability(project_gp_to_nodes(mesh, rho), cfg)
    source_cap_gp = _node_to_gp_mean(mesh, source_cap_node)
    rho_sat_gp = _node_to_gp_mean(mesh, rho_sat_node)
    drho_emit = max(float(getattr(pz, 'rho_increment_per_event', 0.0)), 0.0) * P_storage_gp * fw_gp * source_cap_gp

    # Static Arrhenius recovery/escape.
    rho_safe = np.maximum(rho, 1e6)
    rec_rate = np.zeros_like(rho)
    if bool(getattr(pz, 'recovery_enabled', True)) and T > 0:
        n = max(float(getattr(pz, 'recovery_rho_power', 1.0)), 0.0)
        model = str(getattr(pz, 'recovery_model', 'arrhenius')).lower()
        if model == 'climb_diffusion':
            # Physically motivated recovery scale from the lattice-diffusion
            # parameters already carried by DislocationConfig.  The D/b^2 scale
            # is an atomic-jump frequency; kprime remains a dimensionless
            # coarse-graining factor for climb/annihilation geometry.
            disl = getattr(cfg, 'dislocations', None)
            b = max(float(getattr(material, 'b', 2.5e-10)), 1e-30)
            Dl0a = float(getattr(disl, 'Dl0a', 0.0)) if disl is not None else 0.0
            Dl0b = float(getattr(disl, 'Dl0b', 0.0)) if disl is not None else 0.0
            Ea = float(getattr(disl, 'Ea_eV', 0.0)) * EV_TO_J if disl is not None else 0.0
            Eb = float(getattr(disl, 'Eb_eV', 0.0)) * EV_TO_J if disl is not None else 0.0
            Dl = Dl0a*np.exp(-Ea/max(KB*T,1e-30)) + Dl0b*np.exp(-Eb/max(KB*T,1e-30))
            kprime = max(float(getattr(disl, 'kprime', 1.0)), 0.0) if disl is not None else 1.0
            eta_eff = kprime * max(Dl, 0.0) / (b*b)
            rec_rate = eta_eff * np.power(rho_safe, n)
        else:
            eta = max(float(getattr(pz, 'recovery_eta0', 0.0)), 0.0)
            Q = max(float(getattr(pz, 'recovery_Q_eV', 0.0)), 0.0) * EV_TO_J
            rec_rate = eta * np.exp(-Q / max(KB*T, 1e-30)) * np.power(rho_safe, n)
        dyn = max(float(getattr(pz, 'dynamic_recovery_coeff', 0.0)), 0.0) * rho_safe * np.maximum(dot_ep_gp, 0.0)
        rec_rate = rec_rate + dyn
        if bool(getattr(pz, 'emission_recovery_front_only', True)):
            rec_rate = rec_rate * np.maximum(fw_gp, 1e-3)
    drho_rec = rec_rate * max(dt, 0.0)

    rho_new = rho + drho_emit - drho_rec
    rel_cap = float(getattr(cfg.dislocations, 'max_rho_relative_increment', np.inf))
    if np.isfinite(rel_cap) and rel_cap > 0:
        upper = rho * (1.0 + rel_cap)
        lower = rho / (1.0 + rel_cap)
        rho_new = np.minimum(np.maximum(rho_new, lower), upper)
    rho_upper = float(getattr(cfg.dislocations, 'rho_cap', np.inf))
    if bool(getattr(pz, 'rho_source_saturation_cap_enabled', False)):
        rho_upper_arr = np.minimum(rho_upper, np.maximum(rho_sat_gp, 1e6))
        rho_new = np.minimum(rho_new, rho_upper_arr)
    rho_new = np.clip(rho_new, 1e6, rho_upper)
    return rho_new, {'drho_emit_gp': drho_emit, 'drho_rec_gp': drho_rec,
                     'rho_recovery_rate_gp': rec_rate, 'emission_gp': P_emit_gp,
                     'storage_gp': P_storage_gp, 'source_capacity_gp': source_cap_gp,
                     'rho_source_sat_gp': rho_sat_gp}


def _merge_plastic_infos(infos, ne, total_dt):
    """Merge plastic-info dictionaries from adaptive substeps."""
    z = np.zeros(ne)
    if not infos:
        return {
            'dWp_requested_gp': z.copy(), 'dWp_accepted_gp': z.copy(),
            'dep_eq_requested_gp': z.copy(), 'dep_eq_accepted_gp': z.copy(),
            'thermo_scale_gp': np.ones(ne), 'thermo_admissible_gp': z.copy(),
            'thermo_hazard_gp': z.copy(), 'thermo_mode': 'off',
            'thermo_substeps': 0, 'thermo_dt_min': 0.0, 'thermo_retry_count': 0,
        }
    out = {}
    for name in ['dWp_requested_gp', 'dWp_accepted_gp', 'dep_eq_requested_gp', 'dep_eq_accepted_gp']:
        out[name] = np.sum([np.asarray(info.get(name, z), dtype=float) for info in infos], axis=0)
    out['dep_eq_uncapped_gp'] = np.sum([np.asarray(info.get('dep_eq_uncapped_gp', info.get('dep_eq_accepted_gp', z)), dtype=float) for info in infos], axis=0)
    out['dep_eq_limited_gp'] = np.maximum.reduce([np.asarray(info.get('dep_eq_limited_gp', z), dtype=float) for info in infos])
    out['dep_eq_cap'] = np.minimum.reduce([np.asarray(info.get('dep_eq_cap', np.full(ne, np.nan)), dtype=float) for info in infos])
    out['thermo_scale_gp'] = np.minimum.reduce([np.asarray(info.get('thermo_scale_gp', np.ones(ne)), dtype=float) for info in infos])
    out['thermo_admissible_gp'] = np.maximum.reduce([np.asarray(info.get('thermo_admissible_gp', z), dtype=float) for info in infos])
    out['thermo_hazard_gp'] = np.maximum.reduce([np.asarray(info.get('thermo_hazard_gp', z), dtype=float) for info in infos])
    out['thermo_mode'] = infos[-1].get('thermo_mode', 'off')
    out['thermo_substeps'] = int(sum(int(info.get('_substeps', 1)) for info in infos))
    out['thermo_dt_min'] = float(min(float(info.get('_dt_used', total_dt)) for info in infos))
    out['thermo_retry_count'] = int(sum(int(info.get('_retries', 0)) for info in infos))
    return out


def _update_plasticity_maybe_adaptive(ep_gp, rho_gp, sigma_gp, material, T, dt,
                                      plast_model, disl_cfg):
    """Plastic update with optional adaptive kinetic substepping.

    This is a rollback/subcycling device for stiff Arrhenius kinetics.  It does
    not cap the final state directly; it resolves the hazard/Onsager clock until
    each accepted substep is small enough according to configured thermodynamic
    increment criteria.
    """
    ne = sigma_gp.shape[1]
    if not bool(getattr(disl_cfg, 'thermo_adaptive_substepping', False)):
        ep, rho, dot, info = update_plasticity(
            ep_gp, rho_gp, sigma_gp, material, T, dt, plast_model, disl_cfg, return_info=True
        )
        info['_substeps'] = 1
        info['_dt_used'] = dt
        info['_retries'] = 0
        return ep, rho, dot, info

    max_sub = max(int(getattr(disl_cfg, 'thermo_max_substeps', 64)), 1)
    dep_lim = max(float(getattr(disl_cfg, 'thermo_max_dep_increment', np.inf)), 0.0)
    H_lim = max(float(getattr(disl_cfg, 'thermo_max_hazard_increment', np.inf)), 0.0)

    ep = ep_gp.copy(); rho = rho_gp.copy()
    remain = float(max(dt, 0.0))
    sub_dt = remain
    infos = []
    dep_total = np.zeros(ne)
    retries = 0
    used = 0
    min_dt = remain if remain > 0 else 0.0

    # Tiny floor prevents infinite halving.  The final accept protects progress
    # while recording the residual through thermo_retry_count / thermo_dt_min.
    dt_floor = max(remain * 1e-8, 1e-18)
    while remain > dt_floor and used < max_sub:
        trial_dt = min(sub_dt, remain)
        ep_try, rho_try, dot_try, info_try = update_plasticity(
            ep.copy(), rho.copy(), sigma_gp, material, T, trial_dt, plast_model, disl_cfg, return_info=True
        )
        dep_max = float(np.nanmax(np.asarray(info_try.get('dep_eq_accepted_gp', np.zeros(ne)))))
        H_max = float(np.nanmax(np.asarray(info_try.get('thermo_hazard_gp', np.zeros(ne)))))
        too_big = False
        if np.isfinite(dep_lim) and dep_lim > 0 and dep_max > dep_lim:
            too_big = True
        if np.isfinite(H_lim) and H_lim > 0 and H_max > H_lim:
            too_big = True
        if too_big and trial_dt > dt_floor and used + retries < 50 * max_sub:
            sub_dt = max(0.5 * trial_dt, dt_floor)
            retries += 1
            continue

        info_try['_substeps'] = 1
        info_try['_dt_used'] = trial_dt
        info_try['_retries'] = retries
        infos.append(info_try)
        ep, rho = ep_try, rho_try
        dep_total += np.asarray(info_try.get('dep_eq_accepted_gp', np.zeros(ne)), dtype=float)
        remain -= trial_dt
        used += 1
        min_dt = min(min_dt, trial_dt)
        if remain <= dt_floor:
            break
        sub_dt = min(max(1.25 * trial_dt, dt_floor), remain)

    info = _merge_plastic_infos(infos, ne, max(dt, 1e-30))
    info['thermo_substeps'] = used
    info['thermo_dt_min'] = min_dt
    info['thermo_retry_count'] = retries
    dot_avg = dep_total / max(dt, 1e-30)
    return ep, rho, dot_avg, info


def _invalid_state_reason(cfg, diag, Gc_eff):
    """Return a string reason if a diagnostic run has left the useful regime.

    The invalid-state checks are deliberately delayed for the first few
    increments.  At step 1 the trapezoidal external-work accumulator can
    still be zero and, when Wp->Gc coupling is disabled, Gc_local is exactly
    equal to the intrinsic Gc.  Those are not invalid physics states.
    """
    if not bool(getattr(cfg, 'stop_on_invalid', False)):
        return None

    min_step = int(getattr(cfg, 'invalid_min_step', 3))
    if getattr(diag, 'step', 0) < min_step:
        return None

    if diag.rho_max > float(getattr(cfg, 'invalid_rho_max', 5e17)):
        return f"rho_max={diag.rho_max:.3e} exceeded invalid_rho_max"

    # Wp/Wext is meaningful only after positive external work has accumulated.
    Wext_min = float(getattr(cfg, 'invalid_Wext_min', 1e-12))
    if diag.Wext > Wext_min:
        wp_pct = 100.0 * diag.Wp / diag.Wext
        if wp_pct > float(getattr(cfg, 'invalid_wp_wext_pct', 2000.0)):
            return f"Wp/Wext={wp_pct:.2g}% exceeded invalid_wp_wext_pct"

    # Only use the local-Gc cap as an invalid criterion when plastic-work
    # toughening is actually enabled.  In ablation runs with Wp->Gc disabled,
    # Gc_local == Gc_eff by construction and must not stop the run.
    wp_gc_eff = float(getattr(cfg.phase_field, 'plastic_work_to_Gc_efficiency', 0.0))
    if (wp_gc_eff > 0.0 and getattr(diag, 'Gc_local_max', 0.0) > 0
            and cfg.fracture_mode == 'emergent'):
        cap = float(getattr(cfg.phase_field, 'Gc_local_cap_factor', np.inf)) * Gc_eff
        frac = float(getattr(cfg, 'invalid_Gc_factor', 1.001))
        if np.isfinite(cap) and diag.Gc_local_max > frac * cap:
            return f"Gc_local_max={diag.Gc_local_max:.3g} exceeded local-Gc cap"

    sigma_limit = float(getattr(cfg, 'invalid_sigma_eq_GPa', np.inf)) * 1e9
    if np.isfinite(sigma_limit) and getattr(diag, 'sigma_eq_max', 0.0) > sigma_limit:
        return f"sigma_eq_max={diag.sigma_eq_max/1e9:.3g} GPa exceeded invalid_sigma_eq_GPa"

    dep_limit = float(getattr(cfg, 'invalid_dep_eq_increment', np.inf))
    if np.isfinite(dep_limit) and getattr(diag, 'dep_eq_accepted_max', 0.0) > dep_limit:
        return f"dep_eq_accepted_max={diag.dep_eq_accepted_max:.3g} exceeded invalid_dep_eq_increment"

    K_limit = float(getattr(cfg, 'invalid_K_MPa', 100.0)) * 1e6
    if max(abs(diag.KJ_domain), abs(diag.KJ_global)) > K_limit:
        return "KJ exceeded invalid_K_MPa"

    if diag.d_frac > float(getattr(cfg, 'invalid_d_frac', 0.85)):
        return f"d_frac={diag.d_frac:.3f} exceeded invalid_d_frac"

    return None




def configure_sharp_overlay(cfg, tip_h_fine=1.2e-6, tip_ratio=1.15, ell_m=3.0e-6,
                            nx=100, ny=200, micro_Gc_amp=0.0, micro_Gc_corr_m=3.0e-4,
                            micro_Gc_seed=0):
    """Configure the FEM-driven AT2 run with the SHARP-TIP OVERLAY so the crack
    stays localized (no AT2 smearing): adaptive ~1 um tip mesh, ell sized to the
    tip process zone (NOT ell_factor*global_hbar, which is ~300 um and smears),
    and the source-mode first-passage frontier gating. Optionally a heterogeneous
    microstructure Gc(x) for emergent meander. Returns cfg."""
    cfg.mesh.nx, cfg.mesh.ny = nx, ny
    cfg.mesh.tip_h_fine = tip_h_fine        # adaptive radial refinement at the tip
    cfg.mesh.tip_ratio = tip_ratio
    cfg.mesh.ell_absolute_m = ell_m         # resolve the process zone
    cfg.process_zone.crack_hazard_enabled = True   # sharp-tip source-mode overlay
    cfg.process_zone.crack_first_passage = True
    cfg.micro_Gc_amp = micro_Gc_amp
    cfg.micro_Gc_corr_m = micro_Gc_corr_m
    cfg.micro_Gc_seed = micro_Gc_seed
    return cfg


def _build_micro_Gc(mesh, cfg):
    """Smooth random microstructural Gc(x) FACTOR field on the mesh nodes (mean 1).
    Multiplies the local fracture energy so the crack path is selected by the real
    elastic field against a heterogeneous toughness -- meander/branch emerge from
    energy minimization (update_phase_field), they are NOT imposed. Controlled by
    cfg.micro_Gc_amp (0 = homogeneous/off), cfg.micro_Gc_corr_m, cfg.micro_Gc_seed.
    Built from a few random Fourier modes so it is smooth on any unstructured mesh."""
    amp = float(getattr(cfg, 'micro_Gc_amp', 0.0) or 0.0)
    if amp <= 0.0:
        return np.ones(mesh.nn)
    corr = float(getattr(cfg, 'micro_Gc_corr_m', 0.0) or (0.15 * cfg.geometry.Lx))
    seed = int(getattr(cfg, 'micro_Gc_seed', 0))
    rng = np.random.default_rng(seed)
    x, y = mesh.nodes[:, 0], mesh.nodes[:, 1]
    field = np.zeros(mesh.nn)
    nmodes = 24
    for _ in range(nmodes):
        kdir = rng.normal(size=2); kdir /= (np.linalg.norm(kdir) + 1e-12)
        k = kdir / corr
        field += np.cos(k[0] * x + k[1] * y + rng.uniform(0, 2 * np.pi))
    field = (field - field.mean()) / (field.std() + 1e-12)
    return np.clip(1.0 + amp * field, 0.3, 2.5)


def run_simulation(cfg: SimulationConfig = None) -> Dict[float, SimulationHistory]:
    """
    Run the full Arrhenius fracture simulation across temperatures.

    Parameters
    ----------
    cfg : simulation configuration (default: DBTT preset)

    Returns
    -------
    results : dict mapping temperature -> SimulationHistory
    """
    if cfg is None:
        cfg = make_dbtt_config()

    print("=" * 65)
    print("  ARRHENIUS PHASE-FIELD FRACTURE SIMULATION")
    print(f"  Fracture mode: {cfg.fracture_mode}")
    print(f"  Temperatures: {cfg.T_list}")
    print("=" * 65)

    # ====================== Build mesh (once) ======================
    print("\nBuilding mesh...")
    mesh = make_tri_mesh(cfg.geometry, cfg.mesh, seed=42)
    bnd = make_boundary_data(mesh, cfg.geometry)
    print(f"  Mesh: {mesh.nn} nodes, {mesh.ne} elements, hbar = {mesh.hbar:.4e} m")
    _MICRO_GC = _build_micro_Gc(mesh, cfg)
    if float(getattr(cfg, 'micro_Gc_amp', 0.0) or 0.0) > 0.0:
        print(f"  Heterogeneous Gc(x): amp={cfg.micro_Gc_amp}, "
              f"range [{_MICRO_GC.min():.2f}, {_MICRO_GC.max():.2f}] x Gc")

    if getattr(cfg.mesh, 'ell_absolute_m', None) is not None:
        ell = float(cfg.mesh.ell_absolute_m)
        print(f"  Phase-field length scale: ell = {ell:.4e} m "
              f"(fixed physical ell; ell/hbar = {ell/mesh.hbar:.2f})")
    else:
        ell = cfg.mesh.ell_factor * mesh.hbar
        print(f"  Phase-field length scale: ell = {ell:.4e} m "
              f"(ell/hbar = {cfg.mesh.ell_factor:.1f})")

    # Lumped node areas for integrating node-local memory variables.
    node_area = _node_lumped_area(mesh)

    # Elasticity
    D = plane_strain_D(cfg.material)

    # Phase-field matrices
    Md, Kd = assemble_pf_matrices(mesh)

    # Initial damage (notch)
    d0 = np.zeros(mesh.nn)
    d0[bnd.notch_nodes] = 1.0

    # Build plasticity model from config
    plast_model = PlasticityModel(cfg.plasticity_barrier, cfg.material)
    # Make the dislocation config visible to the flow-stress inversion (used by
    # the optional correlated multi-hit Taylor renewal).
    plast_model._disl_cfg = cfg.dislocations

    # --- Startup auto-calibration of the additive Peierls stress floor ---
    if bool(getattr(cfg.dislocations, 'peierls_autocalibrate', False)):
        T_cal = float(getattr(cfg.dislocations, 'peierls_cal_T_K', 0.0) or 0.0)
        if T_cal <= 0.0:
            T_cal = max(cfg.T_list) if len(cfg.T_list) else 300.0
        info = calibrate_peierls_floor(
            cfg.dislocations, cfg.plasticity_barrier.eta0, cfg.material.b, T_cal,
            eps_ref=cfg.dislocations.flow_epsdot_ref)
        print("\n>> Auto-calibrated additive Peierls stress floor")
        print(f"   target floor    = {info['sigma_min_MPa']:.3g} MPa at T_cal = {info['T_cal_K']:.0f} K")
        print(f"   activation entropy S = {info['S_kB']:.1f} kB  (athermal point {info['S_athermal_kB']:.1f} kB)")
        print(f"   solved enthalpy H_P  = {info['H_P_eV']:.3f} eV,  v_P = {info['v_P_b3']:.2f} b^3")
        if info.get('S_was_clamped'):
            print(f"   NOTE: requested S was more negative than the athermal point "
                  f"({info['S_athermal_kB']:.1f} kB) and was clamped to {info['S_kB']:.1f} kB "
                  f"to keep a physical (decreasing-with-T) floor with H_P > 0.")
        if info['unphysical_negative_H']:
            print("   WARNING: H_P < 0 — entropy is more negative than the athermal point; "
                  "floor INCREASES with T (anomalous).  Use S in [-37,-10] kB for a normal Peierls collapse.")
        # Report the resulting floor across the requested temperatures.
        _Ts = sorted(set(list(cfg.T_list) + [T_cal]))
        _fl = [peierls_flow_stress(T, cfg.dislocations.flow_epsdot_ref, plast_model, cfg.dislocations, cfg.material.b)
               for T in _Ts]
        print("   sigma_Peierls(T): " + "  ".join(f"{T:.0f}K={s/1e6:.1f}MPa" for T, s in zip(_Ts, _fl)))

    # --- Crack-tip dislocation-emission barrier (native nanopillar exp_floor) ---
    # Built separately from the bulk plasticity barrier so emission uses the FULL
    # surface-nucleation barrier (its native physical meaning) while the bulk
    # Taylor branch keeps its down-scaled, reduced-|S| form.
    cfg.emission_model = None
    if bool(getattr(cfg.process_zone, 'crack_emission_enabled', False)):
        try:
            import copy as _copy
            pb_em = _copy.deepcopy(cfg.plasticity_barrier)
            pb_em.exp_energy_scale = float(getattr(cfg.process_zone, 'crack_emission_energy_scale', 1.0))
            pb_em.exp_entropy_scale = float(getattr(cfg.process_zone, 'crack_emission_entropy_scale', 1.0))
            pb_em.exp_stress_scale = float(getattr(cfg.process_zone, 'crack_emission_stress_scale', 1.0))
            cfg.emission_model = PlasticityModel(pb_em, cfg.material)
            G0e = cfg.emission_model.G_barrier(np.zeros(1), max(cfg.T_list))[0] / EV_TO_J
            S_em = -pb_em.exp_entropy_scale * float(getattr(pb_em, 'exp_gT_eV_per_K', 0.0)) * EV_TO_J / KB
            print("\n>> Crack-tip dislocation-emission barrier (Rice-Thomson competition) enabled")
            print(f"   native exp_floor: energy_scale={pb_em.exp_energy_scale}, entropy_scale={pb_em.exp_entropy_scale} (S~{S_em:.1f} kB)")
            print(f"   emission G0(maxT)={G0e:.3f} eV ; rho_emit_max={cfg.process_zone.crack_emission_rho_max:.1e} /m^2")
        except Exception as _e:
            print(f"   WARNING: could not build emission barrier ({_e}); emission disabled")
            cfg.process_zone.crack_emission_enabled = False

    # Build fracture model
    frac_model = FractureModel(
        cfg.fracture_barrier, cfg.material,
        cfg.hazard, cfg.phase_field
    )

    # Cohesive strength (if cohesive branch)
    sigma_coh = 0.0
    if cfg.cohesive.enabled:
        ell_coh = max(cfg.cohesive.length_factor * ell, 1e-12)
        Eprime = cfg.material.Eprime
        sigma_coh = cfg.cohesive.strength_factor * np.sqrt(
            Eprime * max(cfg.cohesive.Gc, 1e-30) / max(ell_coh, 1e-30)
        )
        print(f"  Cohesive: Gc={cfg.cohesive.Gc:.4g} J/m², "
              f"sigma_coh={sigma_coh/1e6:.1f} MPa")

    # Output directory.  In sweeps this is set to a unique parameter-specific
    # directory so individual run histories and snapshots are not overwritten.
    outdir = cfg.output_dir
    if not os.path.isabs(outdir):
        outdir = os.path.join(os.getcwd(), outdir)
    os.makedirs(outdir, exist_ok=True)

    # ====================== Temperature loop ======================
    results = {}

    for iT, T0 in enumerate(cfg.T_list):
        print(f"\n{'='*55}")
        print(f"  T = {T0:.0f} K  ({iT+1}/{len(cfg.T_list)})")
        print(f"{'='*55}")

        t_start = time.time()

        # --- Compute Gc(T) ---
        if cfg.fracture_mode == 'emergent':
            # CONSTANT Gc: toughness emerges from plasticity-fracture competition
            Gc_eff = cfg.phase_field.Gc0_athermal
            Kc_input = np.sqrt(cfg.material.Eprime * Gc_eff)
            print(f"  Gc_eff = {Gc_eff:.4g} J/m² (CONSTANT — apparent Kc emerges from simulation)")
        elif cfg.fracture_mode == 'cohesive_dbtt':
            Gc_eff = cfg.cohesive.Gc
            Kc_input = np.sqrt(cfg.material.Eprime * Gc_eff)
            print(f"  Gc_eff = {Gc_eff:.4g} J/m²")
        else:
            Gc_eff = frac_model.Gc_of_T(T0, ell, method=cfg.toughness_method)
            Kc_input = np.sqrt(cfg.material.Eprime * Gc_eff)
            print(f"  Gc_eff = {Gc_eff:.4g} J/m²")
        print(f"  Kc_input = {Kc_input/1e6:.3f} MPa·√m")

        # --- Initialize state ---
        u = np.zeros(mesh.ndof)
        d = d0.copy()
        Hhist = np.zeros(mesh.nn)
        ep_gp = np.zeros((3, mesh.ne))
        rho_gp = np.ones(mesh.ne) * cfg.dislocations.rho0
        dot_ep_gp = np.zeros(mesh.ne)

        # Tip memory
        rtip_ref = ell / 2
        rtip_state = np.full(mesh.nn, rtip_ref)
        shield_state = np.zeros(mesh.nn)
        rtip_amp = np.ones(mesh.nn)
        shield_factor = np.ones(mesh.nn)

        # Accumulated gross plastic work density per GP and retained
        # process-zone toughening state q_Gc [J/m^2].  q_Gc is the
        # thermodynamic replacement for the old direct Wp -> Gc map.
        wp_gp = np.zeros(mesh.ne)  # cumulative gross plastic work density [J/m³]
        q_tough_node = np.zeros(mesh.nn)  # retained local fracture resistance [J/m²]
        B_crack_node = np.zeros(mesh.nn)  # accumulated crack first-passage action B=int(lambda dt)
        B_emit_node = np.zeros(mesh.nn)   # accumulated crack-tip emission first-passage action
        Wp_tip_cum = 0.0                  # cumulative crack-tip emission dissipation [J/m]
        # History
        hist = SimulationHistory(
            T=T0, Gc_eff=Gc_eff, Kc_input=Kc_input, ell=ell
        )

        # Initial field snapshot (step 0) so every saved run has a baseline.
        if cfg.diagnostics.save_fields:
            hist.d_fields[0] = d.copy()
            hist.u_fields[0] = u.copy()
            hist.rho_fields[0] = project_gp_to_nodes(mesh, rho_gp)
            hist.rtip_fields[0] = rtip_state.copy()
            hist.shield_fields[0] = shield_state.copy()
            hist.wp_fields[0] = project_gp_to_nodes(mesh, wp_gp)
            if cfg.fracture_mode == 'emergent':
                hist.Gc_fields[0] = np.full(mesh.nn, Gc_eff)
            hist.M_fields[0] = np.ones(mesh.nn)

        # Auto-stop trackers
        Fmax = 0.0
        n_quiet = 0

        # Energy accumulators.  Wext_cum is the signed two-boundary
        # external work used for the thermodynamic audit.  Wext_top_cum is the
        # legacy/top-reaction work metric; Wext_abs_cum is a positive-magnitude
        # work diagnostic useful for detecting sign/thickness normalization
        # errors.
        Wext_cum = 0.0
        Wext_top_cum = 0.0
        Wext_abs_cum = 0.0
        Wp_cum = 0.0              # gross plastic work
        Dp_eff_cum = 0.0          # plastic dissipation after retained toughening partition
        Etough_cum = 0.0
        Dtough_cum = 0.0
        Emem_cum = 0.0
        Dmem_cum = 0.0
        Dfrac_cum = 0.0
        U_prev = 0.0
        Uy_top_prev = 0.0
        Uy_bot_prev = 0.0
        F_prev = 0.0
        Ftop_prev = 0.0
        Fbot_prev = 0.0
        Wext_prev_cum = 0.0
        Wext_top_prev_cum = 0.0
        Wext_abs_prev_cum = 0.0
        Wp_prev_cum = 0.0
        Dp_eff_prev_cum = 0.0
        Etough_prev = 0.0
        Dtough_prev = 0.0
        Uel_prev = 0.0
        Uel_drive_prev = 0.0
        # Reference-state surface energy of the initial notch.  This must be
        # measured by the SAME in-loop code path that computes Epf each step
        # (a standalone at2_surface_energy of the clamped notch band does NOT
        # match the relaxed in-loop value).  Captured on the first step below;
        # subtracted from the incremental and cumulative balance so the
        # pre-existing crack is not counted as freshly created by external work.
        Epf0 = None
        Epf_prev = 0.0
        Emem_prev = 0.0
        Dmem_prev = 0.0
        Dfrac_prev = 0.0

        # --- Load stepping ---
        dt = cfg.loading.dt
        n_steps = cfg.loading.n_steps
        dU = cfg.loading.dU_top

        # --- Live progress monitoring state ---
        _prog_on = bool(getattr(cfg.diagnostics, 'progress', True))
        _prog_interval = max(float(getattr(cfg.diagnostics, 'progress_interval_s', 15.0)), 0.0)
        _prog_every = max(int(getattr(cfg.diagnostics, 'progress_every', 1)), 1)
        _prog_t_run0 = time.time()
        _prog = {'t_step0': _prog_t_run0, 'last_beat': _prog_t_run0,
                 'step_times': []}

        def beat(phase):
            """Wall-clock-throttled heartbeat: prints + writes progress.json the
            first time `phase` is reached after `progress_interval_s` has elapsed.
            Reads live loop variables (step, it, Fstep, rho_gp) at call time so a
            single slow step still ticks and localizes the hang to a phase."""
            if not _prog_on:
                return
            now = time.time()
            if now - _prog['last_beat'] < _prog_interval:
                return
            _prog['last_beat'] = now
            try:
                rmax = float(rho_gp.max()) if rho_gp.size else 0.0
            except Exception:
                rmax = 0.0
            elapsed_step = now - _prog['t_step0']
            elapsed_run = now - _prog_t_run0
            try:
                stg = f"{it+1}/{cfg.loading.n_stagger}"
            except Exception:
                stg = "?"
            print(f"    .. T={T0:.0f}K step {step}/{n_steps} stagger {stg} [{phase}] "
                  f"step_elapsed={elapsed_step:5.1f}s run={elapsed_run/60:.1f}m "
                  f"|F|={abs(Fstep):.3g} rho_max={rmax:.2e}", flush=True)
            _write_progress(outdir, {
                'state': 'running', 'phase': phase, 'temperature_K': float(T0),
                'step': int(step), 'n_steps': int(n_steps),
                'stagger': stg,
                'step_elapsed_s': round(elapsed_step, 2),
                'run_elapsed_s': round(elapsed_run, 2),
                'F_abs': float(abs(Fstep)), 'rho_max': rmax,
                'mean_step_time_s': (round(float(np.mean(_prog['step_times'])), 2)
                                     if _prog['step_times'] else None),
                'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
            })

        KJ_license_prev = 0.0
        fired_cum = np.zeros(mesh.nn, dtype=bool)

        for step in range(1, n_steps + 1):
            _prog['t_step0'] = time.time()
            _show_step = _prog_on and (step <= 3 or step % _prog_every == 0 or step == n_steps)

            d_old = d.copy()

            # Symmetric opening
            Uy_top = step * dU / 2
            Uy_bot = -step * dU / 2
            Uapp = Uy_top - Uy_bot

            if _show_step:
                print(f"  [T={T0:.0f}K] step {step}/{n_steps} start  Uapp={Uapp:.3e}", flush=True)
            _consume_n_step = 0
            _consume_dq_step = 0.0

            Fstep = 0.0
            dWp_step_gp = np.zeros(mesh.ne)
            dEtough_step = 0.0
            dDtough_step = 0.0
            dqtough_step_node = np.zeros(mesh.nn)
            tough_weight_node = np.zeros(mesh.nn)
            M_fract = np.ones(mesh.nn)
            tip_emit_prob = np.zeros(mesh.nn)
            pz_emit_info = {
                'P_emit': np.zeros(mesh.nn), 'sigma_tip_eff': np.zeros(mesh.nn),
                'sigma_back': np.zeros(mesh.nn), 'emission_rate': np.zeros(mesh.nn),
                'emission_hazard': np.zeros(mesh.nn), 'front_weight': np.zeros(mesh.nn),
            }
            pz_mobility_info = {
                'P_mobile': np.zeros(mesh.nn), 'P_escape': np.zeros(mesh.nn), 'P_store': np.zeros(mesh.nn),
                'P_mobility': np.zeros(mesh.nn), 'mobility_hazard': np.zeros(mesh.nn),
                'mobility_hazard_raw': np.zeros(mesh.nn), 'sigma_mobility_eff': np.zeros(mesh.nn),
                'storage_fraction': np.ones(mesh.nn),
            }
            pz_rho_info = {
                'drho_emit_gp': np.zeros(mesh.ne), 'drho_rec_gp': np.zeros(mesh.ne),
                'rho_recovery_rate_gp': np.zeros(mesh.ne), 'emission_gp': np.zeros(mesh.ne),
                'storage_gp': np.zeros(mesh.ne),
            }
            pz_crack_info = {
                'G_shield': np.zeros(mesh.nn), 'G_stored_release': np.zeros(mesh.nn),
                'e_stored': np.zeros(mesh.nn), 'tau_back_crack': np.zeros(mesh.nn),
                'Gc_net': np.zeros(mesh.nn),
            }
            crack_hazard_info = {
                'P_crack': np.zeros(mesh.nn), 'G_app': np.zeros(mesh.nn),
                'G_eff': np.zeros(mesh.nn), 'barrier_eV': np.full(mesh.nn, np.inf),
                'hazard': np.zeros(mesh.nn), 'hazard_raw': np.zeros(mesh.nn),
                'B_crack': np.zeros(mesh.nn), 'sigma_tip_crack': np.zeros(mesh.nn),
            }
            drive_multiplier = None
            sigma_eq_gp = np.zeros(mesh.ne)
            sigma1_gp = np.zeros(mesh.ne)
            psi_e_gp = np.zeros(mesh.ne)
            sigma_gp = np.zeros((3, mesh.ne))

            # --- Staggered solve ---
            for it in range(cfg.loading.n_stagger):
                beat('mech')

                # (1) Mechanics
                K, Rint, sigma_gp, sigma_eq_gp, sigma1_gp, psi_e_gp = \
                    assemble_mechanics(
                        mesh, u, ep_gp, rho_gp, d, D, cfg.material,
                        cfg.phase_field.kappa
                    )

                u, Fstep = solve_dirichlet(K, Rint, u, bnd, Uy_top, Uy_bot)

                # POST-SOLVE STRESS: the plastic update, hazard drive and
                # phase-field target must see the equilibrium stress at the
                # NEW displacement.  Feeding the pre-solve field makes the
                # plasticity<->mechanics splitting an unstable ratchet
                # (epsilon_p oscillates against a stress it never relaxes,
                # booking plastic work each cycle -- the high-T runaway).
                sigma_gp, sigma_eq_gp, sigma1_gp, psi_e_gp = stress_state(
                    mesh, u, ep_gp, d, D, cfg.material, cfg.phase_field.kappa)

                beat('plast')
                # (2) Plasticity + dislocation evolution
                ep_gp, rho_gp, dot_ep_gp, plastic_info = _update_plasticity_maybe_adaptive(
                    ep_gp, rho_gp, sigma_gp, cfg.material, T0, dt,
                    plast_model, cfg.dislocations
                )

                # Plastic dissipation (per element).  Thermodynamic update modes
                # return the accepted stress-path work density.  This replaces the
                # old pre-return stress-power estimate sigma_eq*dot_ep*dt, which
                # could report plastic dissipation far larger than the mechanical
                # work available in a quasi-static load step.
                dWp_gp = np.asarray(plastic_info.get('dWp_accepted_gp', sigma_eq_gp * dot_ep_gp * dt), dtype=float)
                wp_gp += dWp_gp  # accumulate local plastic work density
                dWp_step_gp += dWp_gp
                dWp = np.sum(dWp_gp * mesh.area_e)
                if np.isfinite(dWp):
                    Wp_cum += dWp

                # (3) Phase-field update
                beat('phase_field')
                psi_e_node = project_gp_to_nodes(mesh, psi_e_gp)

                # Compute spatially varying Gc for emergent mode.  The default
                # thermodynamic path is a retained process-zone toughening state
                # q_Gc(x), driven by accepted plastic work near the crack/front.
                # The old direct Wp->Gc mapping remains available as
                # --wp-gc-coupling-mode direct for ablation.
                if cfg.fracture_mode == 'emergent':
                    wp_node = project_gp_to_nodes(mesh, wp_gp)
                    dwp_node_for_tough = project_gp_to_nodes(mesh, dWp_gp)
                    q_tough_node, tough_info = _update_toughening_state(
                        q_tough_node, dwp_node_for_tough, d, mesh, ell, Gc_eff, cfg, emit_prob_node=None
                    )
                    dEtough_step += float(tough_info.get('dEtough', 0.0))
                    dDtough_step += float(tough_info.get('dDtough', 0.0))
                    dqtough_step_node = np.maximum(dqtough_step_node, np.abs(np.asarray(tough_info.get('dq_node', np.zeros(mesh.nn)))))
                    tough_weight_node = np.maximum(tough_weight_node, np.asarray(tough_info.get('weight', np.zeros(mesh.nn))))
                    Etough_cum = float(tough_info.get('Etough', Etough_cum))
                    Dtough_cum += float(tough_info.get('dDtough', 0.0))
                    rho_node_for_gc = project_gp_to_nodes(mesh, rho_gp)
                    Gc_local, pz_crack_info = _local_Gc_from_state(
                        Gc_eff, wp_node, q_tough_node, ell, cfg,
                        rho_node=rho_node_for_gc, material=cfg.material,
                        plast_model=plast_model, T=T0, return_terms=True, d=d, mesh=mesh
                    )
                else:
                    Gc_local = Gc_eff

                # Tip amplification and optional cohesive gate.
                # In cohesive_dbtt, M_tip enters the local cohesive-opening
                # stress.  In emergent/arrhenius_gc modes, M_tip^2 multiplies
                # the AT2 damage drive, which is equivalent to applying the
                # local stress amplification to the elastic energy density.
                coh_gate = None
                tip_emit_prob = np.zeros(mesh.nn)
                s1_node = project_gp_to_nodes(mesh, sigma1_gp)
                if cfg.tip_memory.enabled:
                    M_fract = compute_fracture_amplification(
                        s1_node, d, mesh, ell, rtip_amp, shield_factor,
                        cfg.tip_memory.M_max,
                        lambda_tip=cfg.tip_memory.lambda_tip,
                        kappa_tip_max=cfg.tip_memory.kappa_tip_max,
                    )
                else:
                    M_fract = np.ones(mesh.nn)

                sig_pos = M_fract * np.maximum(s1_node, 0.0)

                if cfg.cohesive.enabled and sigma_coh > 0:
                    coh_gate = cohesive_gate(
                        sig_pos, sigma_coh,
                        cfg.cohesive.gate_width, cfg.cohesive.gate_floor
                    )
                    tip_emit_prob = tip_emission_probability(
                        sig_pos, T0, dt, plast_model, cfg.cohesive
                    )
                    drive_multiplier = None  # avoid double-counting with cohesive gate
                elif cfg.tip_memory.enabled and cfg.tip_memory.couple_to_damage_drive:
                    drive_multiplier = M_fract ** cfg.tip_memory.drive_exponent
                else:
                    drive_multiplier = None

                # Physical process-zone emission branch.  This is independent
                # of the cohesive ablation mode and is active in emergent mode.
                # Emission uses the shielded effective tip stress and the same
                # Arrhenius barrier family as plasticity.  It then updates rho,
                # memory, and optionally q_Gc through audited internal states.
                if getattr(cfg, 'process_zone', None) is not None and bool(getattr(cfg.process_zone, 'enabled', True)):
                    rho_node_for_emit = project_gp_to_nodes(mesh, rho_gp)
                    pz_emit_info = _compute_process_zone_emission(
                        sig_pos, rho_node_for_emit, shield_state, T0, dt,
                        plast_model, cfg.material, cfg, d=d, mesh=mesh, ell=ell
                    )
                    pz_emit_prob = np.asarray(pz_emit_info.get('P_emit', np.zeros(mesh.nn)))
                    pz_mobility_info = _compute_process_zone_mobility_partition(
                        pz_emit_prob,
                        np.asarray(pz_emit_info.get('sigma_tip_eff', np.zeros(mesh.nn))),
                        np.asarray(pz_emit_info.get('sigma_back', np.zeros(mesh.nn))),
                        rho_node_for_emit, T0, dt, plast_model, cfg.material, cfg
                    )
                    # Total emitted probability remains useful as a diagnostic and
                    # for cohesive emission competition.  Retained/stored probability
                    # is what builds rho_pz, shielding memory, and q_Gc.
                    tip_emit_prob = np.maximum(tip_emit_prob, pz_emit_prob)
                    pz_store_prob = np.asarray(pz_mobility_info.get('P_store', pz_emit_prob))
                    rho_gp, pz_rho_info = _apply_process_zone_rho_terms(
                        rho_gp, dot_ep_gp, tip_emit_prob, d, mesh, ell, T0, dt, cfg.material, cfg,
                        P_storage_node=pz_store_prob
                    )
                    # If q_Gc is emission-driven, update retained process-zone
                    # toughness before the damage solve so crack growth sees
                    # the new shielded resistance in the current stagger.
                    drv = str(getattr(cfg.process_zone, 'qgc_driver', 'plastic_work')).lower()
                    if cfg.fracture_mode == 'emergent' and drv in ('emission', 'mixed'):
                        q_tough_node, tough_info_emit = _update_toughening_state(
                            q_tough_node, np.zeros(mesh.nn), d, mesh, ell, Gc_eff, cfg,
                            emit_prob_node=pz_store_prob
                        )
                        dEtough_step += float(tough_info_emit.get('dEtough', 0.0))
                        dDtough_step += float(tough_info_emit.get('dDtough', 0.0))
                        dqtough_step_node = np.maximum(dqtough_step_node, np.abs(np.asarray(tough_info_emit.get('dq_node', np.zeros(mesh.nn)))))
                        tough_weight_node = np.maximum(tough_weight_node, np.asarray(tough_info_emit.get('weight', np.zeros(mesh.nn))))
                        Etough_cum = float(tough_info_emit.get('Etough', Etough_cum))
                        Dtough_cum += float(tough_info_emit.get('dDtough', 0.0))
                        wp_node = project_gp_to_nodes(mesh, wp_gp)
                        rho_node_for_gc = project_gp_to_nodes(mesh, rho_gp)
                        Gc_local, pz_crack_info = _local_Gc_from_state(
                            Gc_eff, wp_node, q_tough_node, ell, cfg,
                            rho_node=rho_node_for_gc, material=cfg.material,
                            plast_model=plast_model, T=T0, return_terms=True, d=d, mesh=mesh
                        )

                # Crack-growth Arrhenius/Eyring hazard.  This branch uses
                # the same process-zone state as emission, but projected onto an
                # energy-release-rate-like crack driving force:
                #   G_eff = G_app - G_shield + G_stored_release.
                # The variational AT2 solution remains the target; the hazard
                # controls how much of that target is kinetically accepted.
                Htrial_for_crack = np.maximum(Hhist, psi_e_node)
                crack_hazard_info = _compute_process_zone_crack_hazard(
                    Htrial_for_crack, Gc_local, pz_crack_info, d, mesh, ell, T0, dt, cfg,
                    sigma_tip_node=sig_pos, B_crack_node=B_crack_node, B_emit_node=B_emit_node
                )
                B_crack_node = np.asarray(crack_hazard_info.get('B_crack', B_crack_node), dtype=float)
                B_emit_node = np.asarray(crack_hazard_info.get('B_emit', B_emit_node), dtype=float)
                P_crack = np.asarray(crack_hazard_info.get('P_crack', np.zeros(mesh.nn)))
                if not (getattr(cfg, 'process_zone', None) is not None and bool(getattr(cfg.process_zone, 'crack_hazard_enabled', False))):
                    P_crack_arg = None
                    crack_fired_arg = None
                else:
                    P_crack_arg = P_crack
                    # Completed first-passage events at the connected front.
                    # Used by (a) the EXPERIMENTAL 'source' advance mode and
                    # (b) cleavage consumption of blunting toughening below.
                    Btar_fire = max(float(getattr(cfg.process_zone, 'crack_B_target', 1.0)), 1e-12)
                    B_fire = np.asarray(crack_hazard_info.get('B_crack', np.zeros(mesh.nn)), dtype=float)
                    fw_fire = np.asarray(crack_hazard_info.get('front_weight', np.zeros(mesh.nn)), dtype=float)
                    # NOTE: B can only accumulate where the front weight (and
                    # hence the rate) was significant, so B >= B_target is
                    # already front-localized.  A strict fw floor here is
                    # counterproductive: the strict front weight SUPPRESSES THE
                    # WAKE (d > ~0.85 -> fw ~ 0), and the B-saturated nodes are
                    # exactly the high-d crescent that accumulated B when they
                    # WERE the front.  Keep only a small floor against noise.
                    fw_min = float(getattr(cfg.process_zone, 'crack_fire_front_weight_min', 0.05))
                    fired_front = (B_fire >= Btar_fire) & (fw_fire >= fw_min)

                    # --- CLEAVAGE CONSUMPTION OF BLUNTING TOUGHENING ---------
                    # The q_blunt halo (Gc_local = Gc0 + q) is R-curve
                    # shielding of the BLUNTED tip.  When the cleavage clock
                    # completes, the crack re-initiates sharply ahead of the
                    # blunted tip and cuts through: the halo cannot block a
                    # fired cleavage event.  Without this, the variational
                    # target freezes (Hc = Gc_local/2ell unreachable) and the
                    # front locks with B at cap and P=1 gating toward nothing.
                    # Energy: the consumed q was retained plastic work; it is
                    # re-booked from toughening storage to dissipation, so the
                    # global audit is conserved.
                    cons = float(getattr(cfg.process_zone, 'crack_fire_consume_toughening_frac', 0.0))
                    if cons > 0.0 and fired_front.any() and cfg.fracture_mode == 'emergent':
                        try:
                            w_cons = np.asarray(tough_info.get('weight', np.zeros(mesh.nn)), dtype=float)
                        except NameError:
                            w_cons = _toughening_weight(mesh, d, ell)
                        # CONSUMPTION FOOTPRINT: the q_blunt halo that blocks
                        # the variational target sits in a ring of radius ~ell
                        # AHEAD of the front, not at the fired node itself.  A
                        # fired cleavage event re-initiates ahead over the
                        # process-zone size, so consume q within
                        # R_c = consume_radius_factor * ell of any fired node
                        # (smoothstep falloff), not just pointwise.
                        Rc = max(float(getattr(cfg.process_zone, 'crack_fire_consume_radius_factor', 1.5)) * float(ell), 1e-30)
                        Xs = mesh.nodes[np.where(fired_front)[0]]
                        dist = np.full(mesh.nn, np.inf)
                        chunk = 2000
                        for _i in range(0, mesh.nn, chunk):
                            _idx = slice(_i, min(_i + chunk, mesh.nn))
                            _dx = mesh.nodes[_idx, 0:1] - Xs[:, 0][None, :]
                            _dy = mesh.nodes[_idx, 1:2] - Xs[:, 1][None, :]
                            dist[_idx] = np.sqrt(np.min(_dx*_dx + _dy*_dy, axis=1))
                        w_fp = np.clip(1.0 - dist / Rc, 0.0, 1.0)
                        w_fp = w_fp*w_fp*(3.0 - 2.0*w_fp)
                        Eold_t = _compute_toughening_energy(q_tough_node, w_cons, mesh, ell, Gc_eff, cfg)
                        q_old_cons = q_tough_node
                        q_tough_node = q_tough_node * (1.0 - np.clip(cons, 0.0, 1.0) * w_fp)
                        Enew_t = _compute_toughening_energy(q_tough_node, w_cons, mesh, ell, Gc_eff, cfg)
                        dE_cons = Enew_t - Eold_t  # <= 0
                        dEtough_step += dE_cons
                        Etough_cum = max(Etough_cum + dE_cons, 0.0)
                        # released toughening storage stays dissipated plastic work
                        dDtough_step += max(-dE_cons, 0.0)
                        Dtough_cum += max(-dE_cons, 0.0)
                        # Self-reporting: never fly blind on whether this fired.
                        _n_f = int(np.count_nonzero(fired_front))
                        _dq = float(np.max(q_old_cons - q_tough_node))
                        consume_n_fired = _n_f
                        consume_dq_max = _dq
                        # Rebuild the local Gc with the consumed state so THIS
                        # stagger's variational target unfreezes.
                        Gc_local, pz_crack_info = _local_Gc_from_state(
                            Gc_eff, wp_node, q_tough_node, ell, cfg,
                            rho_node=rho_node_for_gc, material=cfg.material,
                            plast_model=plast_model, T=T0, return_terms=True, d=d, mesh=mesh
                        )
                    else:
                        consume_n_fired = 0
                        consume_dq_max = 0.0
                    _consume_n_step = max(int(locals().get('_consume_n_step', 0) or 0), int(consume_n_fired))
                    _consume_dq_step = max(float(locals().get('_consume_dq_step', 0.0) or 0.0), float(consume_dq_max))

                    # SOURCE advance mode: completed first-passage events break
                    # the material point directly (sub-grid bond rupture).  In
                    # gate mode the crack can only advance where the AT2 target
                    # moves, which requires the SMEARED stress to reach the AT2
                    # regularization strength sigma_c ~ sqrt(E*Gc/ell) (~2.7 GPa
                    # at ell=2.8e-4) -- the de-smeared clock fires at far lower
                    # applied K and then waits forever, so the emergent Kc is
                    # set by ell, not by the cleavage physics.  Source mode
                    # restores the de-smeared Kc: when the clock completes, the
                    # node breaks and AT2 supplies only the regularization.
                    #
                    # GLOBAL GRIFFITH LICENSE: fired nodes may break only when
                    # the APPLIED stress intensity (previous step's J-integral)
                    # satisfies K_J >= frac * K_G, K_G = sqrt(E' * Gc_eff).
                    # Every NODAL energy license tried before is broken near
                    # the regularized notch (Hhist absorbs the embrittled
                    # drive; Gc_net collapses in the dislocation cloud; the
                    # undegraded psi adjacent to d=1 nodes is unbounded).  The
                    # GLOBAL J-integral is the mesh-objective statement that
                    # the far field pays Gc per unit advance; it blocks exactly
                    # the sub-Griffith notch transient that previously
                    # avalanched, and nothing else: the clock fires at
                    # K_clock ~ sigma_fire*sqrt(2*pi*r_pz) > K_G, so in normal
                    # loading the license is transparent and the emergent
                    # Kc = max(K_G, K_clock).
                    crack_fired_arg = None
                    if str(getattr(cfg.process_zone, 'crack_advance_mode', 'gate')).lower() == 'source':
                        crack_fired_arg = fired_front.copy()
                        g_gate = float(getattr(cfg.process_zone, 'crack_fire_griffith_frac', 1.0))
                        if g_gate > 0.0 and crack_fired_arg.any():
                            K_G = float(np.sqrt(max(cfg.material.Eprime * float(Gc_eff), 0.0)))
                            K_now = float(locals().get('KJ_license_prev', 0.0) or 0.0)
                            if K_now < g_gate * K_G:
                                crack_fired_arg[:] = False
                        # PERSISTENT FIRED MEMORY + RESISTANCE COLLAPSE.
                        # Cleaved sub-grid bonds do not heal: once a node has
                        # fired (clock complete AND globally licensed), its
                        # local fracture resistance is permanently collapsed,
                        #   Gc_local <- relief * Gc_local   (relief ~ 0.05),
                        # and the ordinary AT2 solve breaks it variationally,
                        # paying from the actual elastic field.  No d=1 fiat:
                        # the previous fiat break injected ~Gc*ell of unpaid
                        # gradient energy per node (audit deficit -31% in the
                        # burst).  With relief, Hc = relief*Gc/(2*ell) is
                        # reachable at ordinary smeared stress (~20 MPa), so
                        # the de-smeared clock sets WHEN, the license sets the
                        # MINIMUM K, and conservation is automatic.
                        # CURRENT-FRONT GATE: firing requires the node to be
                        # near the PRESENT tip (front weight from current
                        # seeds), not merely to hold stale B.  Wake/flank nodes
                        # keep B >= 1 for ~tau_relax after the tip passes and
                        # sit adjacent to broken material (in the frontier);
                        # without this gate they re-fire and the band widens
                        # ~2 layers/step into a branching swath.  Consumption
                        # keeps its permissive floor (0.05); SOURCE firing is
                        # gated at fw >= source_fw_min (default 0.3).
                        src_fw_min = float(getattr(cfg.process_zone, 'crack_source_fw_min', 0.3))
                        crack_fired_arg &= (fw_fire >= src_fw_min)
                        # CONNECTED-FRONT CONSTRAINT: only the frontier of the
                        # crack component connected to the notch may fire.
                        if crack_fired_arg.any():
                            frontier = _crack_frontier_mask(
                                mesh, d, bnd.notch_nodes,
                                dthr=float(getattr(cfg.process_zone, 'crack_fire_connect_dthr', 0.8)),
                                layers=int(getattr(cfg.process_zone, 'crack_fire_connect_layers', 2)))
                            crack_fired_arg &= frontier
                        fired_cum |= crack_fired_arg
                        crack_fired_arg = fired_cum  # cap bypass for all cleaved nodes
                        if fired_cum.any():
                            relief = float(getattr(cfg.process_zone, 'crack_fired_gc_relief', 0.05))
                            Gc_local = np.where(fired_cum, np.maximum(relief, 0.0) * Gc_local, Gc_local)
                # Use the degraded, shielded crack-driving force directly in
                # the AT2 target when the crack hazard is active.  This is the
                # crack-growth counterpart of the process-zone force balance:
                #   G_eff  = G_app - G_shield + G_stored_release
                #   Gc_net = Gc0 + q_blunt - G_stored_release; G_shield is not double-counted in resistance.
                # Without this, stored-energy embrittlement could increase the
                # crack event probability but still leave the variational target
                # controlled only by the raw elastic history field.
                psi_for_damage_node = psi_e_node
                if (getattr(cfg, 'process_zone', None) is not None
                    and bool(getattr(cfg.process_zone, 'crack_hazard_enabled', False))
                    and bool(getattr(cfg.process_zone, 'crack_use_effective_drive', True))):
                    H_eff_drive = np.asarray(crack_hazard_info.get('H_eff_drive', np.zeros(mesh.nn)), dtype=float)
                    fw_crack = np.asarray(crack_hazard_info.get('front_weight', np.ones(mesh.nn)), dtype=float)
                    mix = np.clip(float(getattr(cfg.process_zone, 'crack_effective_drive_mix', 1.0)), 0.0, 1.0)
                    w_eff = np.clip(mix * fw_crack, 0.0, 1.0)
                    psi_embrittled = (1.0 - w_eff) * psi_e_node + w_eff * H_eff_drive
                    psi_for_damage_node = np.maximum(psi_e_node, psi_embrittled)

                Gc_local = Gc_local * _MICRO_GC      # heterogeneous microstructure (mean 1)
                if str(getattr(cfg, 'phase_field_model', 'at2')).lower() == 'at1':
                    # Clean AT1: linear local dissipation -> real elastic threshold,
                    # sharp crack, NO source-mode caps/hazards/floors. By default
                    # PURE VARIATIONAL (drive_multiplier off) so the baseline crack
                    # is clean; the (A) tip-kinetic coupling is opt-in and must be
                    # tip-localized (a global M_fract>1 over-nucleates the bulk).
                    at1_kinetic = bool(getattr(cfg, 'phase_field_at1_kinetic', False))
                    d, Hhist = update_phase_field_at1(
                        d, Hhist, psi_for_damage_node, Md, Kd, bnd.notch_nodes,
                        Gc_local, ell, cfg.hazard.Gamma0, dt,
                        drive_multiplier=(drive_multiplier if at1_kinetic else None),
                        cohesive_gate=(coh_gate if at1_kinetic else None),
                        use_kinetic_drive=(cfg.phase_field.use_kinetic_damage_drive if at1_kinetic else False),
                        max_damage_increment=cfg.phase_field.max_damage_increment_per_stagger,
                    )
                else:
                    d, Hhist = update_phase_field(
                        d, Hhist, psi_for_damage_node, Md, Kd, bnd.notch_nodes,
                        Gc_local, ell, cfg.hazard.Gamma0, dt,
                        cohesive_gate=coh_gate,
                        drive_multiplier=drive_multiplier,
                        crack_hazard_probability=P_crack_arg,
                        crack_fired_nodes=crack_fired_arg,
                        use_kinetic_drive=cfg.phase_field.use_kinetic_damage_drive,
                        max_damage_increment=cfg.phase_field.max_damage_increment_per_stagger,
                        damage_drive_cap=cfg.phase_field.damage_drive_cap,
                    )

            # Final mechanical re-equilibration after the last accepted
            # plasticity/damage increments.  Diagnostics and energy increments
            # must use the equilibrated state, not the pre-return stress field.
            K, Rint, sigma_gp, sigma_eq_gp, sigma1_gp, psi_e_gp = \
                assemble_mechanics(
                    mesh, u, ep_gp, rho_gp, d, D, cfg.material,
                    cfg.phase_field.kappa
                )
            u, Fstep = solve_dirichlet(K, Rint, u, bnd, Uy_top, Uy_bot)

            # --- Crack advance and tip memory update (once per physical step) ---
            Da_proj, _, _ = compute_crack_advance(
                mesh, d, cfg.geometry.a0
            )
            dot_ep_node = project_gp_to_nodes(mesh, dot_ep_gp)
            dwp_step_node = project_gp_to_nodes(mesh, dWp_step_gp)
            if np.ndim(Gc_local) == 0:
                Gc_node_for_memory = np.full(mesh.nn, float(Gc_local))
            else:
                Gc_node_for_memory = np.asarray(Gc_local, dtype=float)

            rtip_state, shield_state, rtip_amp, shield_factor, tip_mem_info = \
                update_tip_memory(
                    d, d_old, dot_ep_node, dt, mesh, ell,
                    rtip_state, shield_state, cfg.tip_memory, rtip_ref,
                    dwp_node=dwp_step_node,
                    Gc_node=Gc_node_for_memory,
                    emit_prob_node=np.asarray(pz_mobility_info.get('P_store', tip_emit_prob)),
                    crack_advance=Da_proj * float(getattr(cfg.process_zone, 'crack_advance_memory_erasure', 1.0)),
                )

            # Integrate node-local memory storage and dissipation.  This gives
            # the memory/tip-amplification state a conjugate energy and a
            # dissipative cost instead of letting it act as a free multiplier.
            Emem_node = np.asarray(tip_mem_info.get('Emem_node', np.zeros(mesh.nn)), dtype=float)
            dDmem_node = np.asarray(tip_mem_info.get('dDmem_node', np.zeros(mesh.nn)), dtype=float)
            Emem_cum = float(np.sum(Emem_node * node_area))
            Dmem_cum += float(np.sum(np.maximum(dDmem_node, 0.0) * node_area))

            # --- Diagnostics ---
            # Boundary reactions and external work.  The imposed loading moves
            # both top and bottom boundaries, so the thermodynamic work audit
            # must use both signed boundary reactions.  The legacy top-reaction
            # work and a positive-magnitude work are retained as diagnostics.
            Ftop_now, Fbot_now, Fpair_abs_now = boundary_reaction_forces(K, Rint, u, bnd)
            U_now = Uapp
            F_now = Ftop_now
            dUy_top = Uy_top - Uy_top_prev
            dUy_bot = Uy_bot - Uy_bot_prev
            dUapp = U_now - U_prev
            dWext_pair = 0.5 * (Ftop_prev + Ftop_now) * dUy_top + \
                         0.5 * (Fbot_prev + Fbot_now) * dUy_bot
            dWext_top = 0.5 * (F_prev + F_now) * dUapp
            dWext_abs = 0.5 * (abs(Ftop_prev) + abs(Ftop_now)) * abs(dUy_top) + \
                        0.5 * (abs(Fbot_prev) + abs(Fbot_now)) * abs(dUy_bot)
            Wext_cum += dWext_pair
            Wext_top_cum += dWext_top
            Wext_abs_cum += dWext_abs
            U_prev = U_now
            Uy_top_prev = Uy_top
            Uy_bot_prev = Uy_bot
            F_prev = F_now
            Ftop_prev = Ftop_now
            Fbot_prev = Fbot_now

            # Elastic energy diagnostics.  psi_e_gp from assemble_mechanics is
            # the undegraded positive energy used as the phase-field damage
            # drive.  The thermodynamic audit must instead use the stored
            # degraded elastic energy 1/2 eps_e:sigma_degraded.
            psi_store_gp, psi_undeg_gp = elastic_energy_densities(mesh, u, ep_gp, sigma_gp, D)
            Uel = float(np.sum(psi_store_gp * mesh.area_e))
            Uel_drive = float(np.sum(psi_e_gp * mesh.area_e))
            Uel_undegraded = float(np.sum(psi_undeg_gp * mesh.area_e))

            # AT2 surface energy
            if cfg.fracture_mode == 'emergent':
                wp_node_diag = project_gp_to_nodes(mesh, wp_gp)
                rho_node_diag = project_gp_to_nodes(mesh, rho_gp)
                Gc_diag, pz_crack_info_diag = _local_Gc_from_state(
                    Gc_eff, wp_node_diag, q_tough_node, ell, cfg,
                    rho_node=rho_node_diag, material=cfg.material,
                    plast_model=plast_model, T=T0, return_terms=True, d=d, mesh=mesh
                )
                # LEDGER CONSISTENCY: price the AT2 surface energy with the
                # SAME Gc field the variational solve actually used.  In
                # source mode fired nodes carry collapsed resistance
                # (relief*Gc); pricing their profile at the unrelieved Gc
                # books ~1/relief times the energy the field transferred and
                # fabricates a monotone audit deficit.  The macroscopic
                # toughness remains KJ_sel (license/clock); along the cleaved
                # path the BOOKED surface energy is relief*Gc, with the
                # excess release G - relief*Gc remaining in the mechanical
                # ledger (the quasi-static analogue of kinetic dissipation in
                # real cleavage, where G >> 2*gamma_s).
                if (str(getattr(cfg.process_zone, 'crack_advance_mode', 'gate')).lower() == 'source'
                        and bool(np.any(fired_cum))):
                    _relief_d = float(getattr(cfg.process_zone, 'crack_fired_gc_relief', 0.05))
                    Gc_diag = np.where(fired_cum, max(_relief_d, 0.0) * Gc_diag, Gc_diag)
                Epf = at2_surface_energy(mesh, d, ell, Gc_diag)
            else:
                Epf = at2_surface_energy(mesh, d, ell, Gc_eff)

            # Capture the notch reference from the first consistent in-loop
            # measurement, so dEpf at step 1 is ~0 (no crack created yet) and
            # the cumulative balance subtracts a reference on the same scale as
            # the per-step Epf.
            if Epf0 is None:
                Epf0 = Epf
                Epf_prev = Epf

            # Optional irreversible damage dissipation.  This is intentionally
            # separate from the AT2 surface-energy state Epf: it estimates the
            # local damage-production sink associated with d increasing during
            # a hazard/phase-field crack advance, using the agreed leading term
            # D_frac = integral (Gc/(2 ell)) d_mid Delta d dV.
            if bool(getattr(cfg.phase_field, 'include_fracture_dissipation_audit', True)):
                dd_node_frac = np.maximum(d - d_old, 0.0)
                d_mid_frac = 0.5 * (d + d_old)
                if cfg.fracture_mode == 'emergent':
                    Gc_for_Dfrac = np.asarray(Gc_diag, dtype=float)
                else:
                    Gc_for_Dfrac = np.full(mesh.nn, float(Gc_eff))
                dDfrac = float(np.sum(np.maximum(Gc_for_Dfrac, 0.0) / (2.0 * max(ell, 1e-30)) * d_mid_frac * dd_node_frac * node_area))
            else:
                dDfrac = 0.0
            Dfrac_cum += max(dDfrac, 0.0)

            # Incremental thermodynamic energy audit.  Negative residual means
            # the accepted coupled update spent more energy than was supplied
            # by the external work increment plus released stored/free energy.
            dWext = Wext_cum - Wext_prev_cum
            dWext_top_diag = Wext_top_cum - Wext_top_prev_cum
            dWext_abs_diag = Wext_abs_cum - Wext_abs_prev_cum
            dUel = Uel - Uel_prev
            dUel_drive = Uel_drive - Uel_drive_prev
            dEpf = Epf - Epf_prev
            dWp_step_total = Wp_cum - Wp_prev_cum
            dEtough = Etough_cum - Etough_prev
            dDtough = Dtough_cum - Dtough_prev
            # Gross plastic work is partitioned into immediate plastic
            # dissipation plus retained process-zone toughening storage and
            # toughening dissipation.  This prevents the same accepted plastic
            # work from being counted once as Wp and a second time as a free Gc
            # increase.  If the toughening state demands more energy than the
            # plastic work source provides, the residual will expose it.
            if bool(getattr(cfg.phase_field, 'toughening_include_in_energy_audit', True)):
                Dp_eff_cum = max(Wp_cum - Etough_cum - Dtough_cum, 0.0)
            else:
                Dp_eff_cum = Wp_cum
            dDp_eff = Dp_eff_cum - Dp_eff_prev_cum
            dEmem = Emem_cum - Emem_prev
            dDmem = Dmem_cum - Dmem_prev
            dDfrac_inc = Dfrac_cum - Dfrac_prev
            sink_increment = dUel + dEpf + dDfrac_inc + dDp_eff + dEtough + dDtough + dEmem + dDmem
            energy_residual = dWext - sink_increment
            energy_residual_topWext = dWext_top_diag - sink_increment
            energy_residual_absWext = dWext_abs_diag - sink_increment
            energy_cumulative_residual = Wext_cum - Uel - (Epf - Epf0) - Dfrac_cum - Dp_eff_cum - Etough_cum - Dtough_cum - Emem_cum - Dmem_cum
            energy_cumulative_residual_absWext = Wext_abs_cum - Uel - (Epf - Epf0) - Dfrac_cum - Dp_eff_cum - Etough_cum - Dtough_cum - Emem_cum - Dmem_cum
            denom_energy = max(abs(dWext) + abs(dUel) + abs(dEpf) + abs(dDfrac_inc) +
                               abs(dDp_eff) + abs(dEtough) + abs(dDtough) + abs(dEmem) + abs(dDmem), 1e-30)
            energy_residual_rel = energy_residual / denom_energy
            energy_units_ratio = Uel / max(Wext_abs_cum, 1e-30)
            abs_tol = float(getattr(cfg.dislocations, 'thermo_energy_abs_tol', 1e-12))
            rel_tol = float(getattr(cfg.dislocations, 'thermo_energy_rel_tol', 0.05))
            # Cumulative audit scale: a fraction of total work done / energy
            # stored.  Robust to the per-step denominator collapsing to zero
            # once the specimen has failed (Freact -> 0) while damage smearing
            # still nudges Epf.  A step is admissible if EITHER the incremental
            # balance holds OR the running cumulative balance is within rel_tol
            # of the total external work.  The cumulative test is the physically
            # meaningful FE energy-conservation statement.
            audit_scale = max(abs(Wext_abs_cum), abs(Uel) + abs(Epf - Epf0), 1e-30)
            cum_rel = energy_cumulative_residual_absWext / audit_scale
            if bool(getattr(cfg.dislocations, 'thermo_energy_audit', True)):
                inc_ok = (energy_residual >= -abs_tol) or (energy_residual_rel >= -rel_tol)
                cum_ok = (cum_rel >= -rel_tol)
                energy_balance_ok = float(inc_ok or cum_ok)
            else:
                energy_balance_ok = 1.0

            # R-curve/tearing diagnostic: external work minus stored elastic,
            # plastic, and memory costs per projected crack advance.  This is
            # the preferred scalar for soft tearing/mixed modes; K_Ic-like
            # values are only clean for brittle onset.
            numerator_R = Wext_cum - Uel - Dfrac_cum - Dp_eff_cum - Etough_cum - Dtough_cum - Emem_cum - Dmem_cum
            if Da_proj > 1e-12 and numerator_R > 0:
                J_tearing = numerator_R / Da_proj
                KJ_tearing = np.sqrt(max(J_tearing, 0.0) * cfg.material.Eprime)
            else:
                J_tearing = 0.0
                KJ_tearing = 0.0

            # Crack advance was computed before tip-memory update.
            Gamma_total = Epf / max(Gc_eff, 1e-30)
            branch_factor = Gamma_total / max(Da_proj, 1e-12)

            # Domain-integral J
            tip, direction = find_crack_tip(mesh, d, cfg.geometry.a0)
            J_domain, KJ_domain, _ = compute_J_integral(
                mesh, u, sigma_gp, psi_e_gp, d,
                tip, direction, cfg.material, ell, cfg.j_integral
            )

            # Global energy balance J.  This is diagnostic only; it is not
            # used as the default selected toughness because it blows up when
            # Da_projected is tiny or after the crack has crossed the window.
            if Da_proj > 1e-12 and Wext_cum > 0:
                J_global = Wext_cum / Da_proj
                KJ_global = np.sqrt(max(J_global, 0) * cfg.material.Eprime)
            else:
                J_global = 0.0
                KJ_global = 0.0

            # Independent force-based LEFM estimate for elastic calibration.
            K_force = edge_crack_force_K(Fstep, cfg.geometry, cfg.material)

            # Rho node values for diagnostics
            rho_node = project_gp_to_nodes(mesh, rho_gp)

            # Distribution diagnostics for rho. Mean/max alone are misleading
            # when rho evolves in a very small crack-tip/process-zone region.
            rho_p95 = float(np.percentile(rho_gp, 95)) if rho_gp.size else 0.0
            rho_p99 = float(np.percentile(rho_gp, 99)) if rho_gp.size else 0.0
            rho_cap_val = float(getattr(cfg.dislocations, 'rho_cap', np.inf))
            rho_cap_frac = float(np.mean(rho_gp >= 0.999 * rho_cap_val)) if np.isfinite(rho_cap_val) and rho_gp.size else 0.0
            rho_gt_1e14_frac = float(np.mean(rho_gp > 1e14)) if rho_gp.size else 0.0
            rho_gt_1e15_frac = float(np.mean(rho_gp > 1e15)) if rho_gp.size else 0.0
            rho_gt_1e16_frac = float(np.mean(rho_gp > 1e16)) if rho_gp.size else 0.0

            # Plastic flow/yield threshold diagnostics.  This is essential for
            # EXP_floor barriers: it tells whether no plasticity means the local
            # stress is below sigma_y, the inversion is floor-limited, or the
            # whole specimen is yielding.
            flow_diag = plastic_flow_diagnostics(
                rho_gp, sigma_gp, cfg.material, T0, plast_model, cfg.dislocations
            )

            # Record diagnostics
            # --- Crack-tip emission dissipation bookkeeping (Wp_tip) ---------
            # The emission channel IS the sub-grid tip plasticity at this mesh
            # scale; book its dissipation so "plastic work" is visible.  Energy
            # dissipated = shielding removed from the drive per unit new crack
            # surface:  Wp_tip += sum( G_shield_emit * dGamma ), with the AT2
            # crack-surface-density increment dGamma ~ (d^2-d_old^2)/(2*ell)*V_n
            # (local part only; the gradient part adds a comparable amount, so
            # this is a mild underestimate, NOT injected into the energy audit).
            try:
                _Gse = np.asarray(crack_hazard_info.get('G_shield_emit', np.zeros(mesh.nn)), dtype=float)
                if _Gse.any():
                    _Vn = (cfg.geometry.Lx * cfg.geometry.Ly) / mesh.nn
                    _dGam = np.maximum(d**2 - d_old**2, 0.0) / (2.0 * max(ell, 1e-30)) * _Vn
                    Wp_tip_cum += float(np.sum(_Gse * _dGam))
            except Exception:
                pass

            KJ_license_prev = max(float(KJ_domain or 0.0), float(KJ_global or 0.0))
            diag = StepDiagnostics(
                step=step, t=step*dt,
                Uapp=Uapp, Freact=Ftop_now, Ftop=Ftop_now, Fbot=Fbot_now, Fpair_abs=Fpair_abs_now,
                Wext=Wext_cum, Wext_top=Wext_top_cum, Wext_pair=Wext_cum, Wext_abs=Wext_abs_cum,
                Uel=Uel, Uel_drive=Uel_drive, Uel_undegraded=Uel_undegraded,
                Wp=Wp_cum, Wp_tip=Wp_tip_cum, Dp_eff=Dp_eff_cum, Etough=Etough_cum, Dtough=Dtough_cum, Epf_surf=Epf,
                dWext=dWext, dWext_top=dWext_top_diag, dWext_pair=dWext, dWext_abs=dWext_abs_diag,
                dUel=dUel, dUel_drive=dUel_drive, dEpf=dEpf, dWp_step=dWp_step_total,
                dDp_eff=dDp_eff, dEtough=dEtough, dDtough=dDtough,
                Emem=Emem_cum, Dmem=Dmem_cum, Dfrac=Dfrac_cum, dEmem=dEmem, dDmem=dDmem, dDfrac=dDfrac_inc,
                energy_residual=energy_residual, energy_residual_rel=energy_residual_rel,
                energy_residual_absWext=energy_residual_absWext, energy_residual_topWext=energy_residual_topWext,
                energy_cumulative_residual=energy_cumulative_residual,
                energy_cumulative_residual_absWext=energy_cumulative_residual_absWext,
                energy_units_ratio_Uel_over_Wext_abs=energy_units_ratio,
                energy_balance_ok=energy_balance_ok,
                J_tearing=J_tearing, KJ_tearing=KJ_tearing,
                crack_len=Epf / max(Gc_eff, 1e-30),
                Da_projected=Da_proj, Gamma_total=Gamma_total,
                branch_factor=branch_factor,
                J_domain=J_domain, KJ_domain=KJ_domain,
                J_global=J_global, KJ_global=KJ_global,
                K_force=K_force,
                rho_mean=np.mean(rho_gp), rho_p95=rho_p95, rho_p99=rho_p99, rho_max=np.max(rho_gp),
                rho_gt_1e14_frac=rho_gt_1e14_frac, rho_gt_1e15_frac=rho_gt_1e15_frac,
                rho_gt_1e16_frac=rho_gt_1e16_frac, rho_cap_frac=rho_cap_frac,
                dotep_mean=np.mean(dot_ep_gp), dotep_max=np.max(dot_ep_gp),
                dWp_requested=float(np.sum(np.asarray(plastic_info.get('dWp_requested_gp', np.zeros(mesh.ne))) * mesh.area_e)),
                dWp_accepted=float(np.sum(np.asarray(plastic_info.get('dWp_accepted_gp', np.zeros(mesh.ne))) * mesh.area_e)),
                dep_eq_requested_max=float(np.max(np.asarray(plastic_info.get('dep_eq_requested_gp', np.zeros(mesh.ne))))),
                dep_eq_accepted_max=float(np.max(np.asarray(plastic_info.get('dep_eq_accepted_gp', np.zeros(mesh.ne))))),
                dep_eq_uncapped_max=float(np.max(np.asarray(plastic_info.get('dep_eq_uncapped_gp', plastic_info.get('dep_eq_accepted_gp', np.zeros(mesh.ne)))))),
                dep_limited_frac=float(np.mean(np.asarray(plastic_info.get('dep_eq_limited_gp', np.zeros(mesh.ne))) > 0)),
                thermo_scale_min=float(np.min(np.asarray(plastic_info.get('thermo_scale_gp', np.ones(mesh.ne))))),
                thermo_scale_mean=float(np.mean(np.asarray(plastic_info.get('thermo_scale_gp', np.ones(mesh.ne))))),
                thermo_admissible_frac=float(np.mean(np.asarray(plastic_info.get('thermo_admissible_gp', np.zeros(mesh.ne))) > 0)),
                thermo_hazard_max=float(np.max(np.asarray(plastic_info.get('thermo_hazard_gp', np.zeros(mesh.ne))))),
                thermo_substeps=float(plastic_info.get('thermo_substeps', plastic_info.get('_substeps', 1))),
                thermo_dt_min=float(plastic_info.get('thermo_dt_min', plastic_info.get('_dt_used', dt))),
                thermo_retry_count=float(plastic_info.get('thermo_retry_count', plastic_info.get('_retries', 0))),
                memory_energy_increment=float(np.sum(np.asarray(tip_mem_info.get('dEmem_node', np.zeros(mesh.nn))) * node_area)),
                memory_dissipation_increment=float(np.sum(np.asarray(tip_mem_info.get('dDmem_node', np.zeros(mesh.nn))) * node_area)),
                memory_A_r_mean=float(tip_mem_info.get('memory_A_r_mean', 0.0)),
                memory_A_z_mean=float(tip_mem_info.get('memory_A_z_mean', 0.0)),
                d_frac=np.mean(d > cfg.diagnostics.damage_threshold),
                plast_frac=np.mean(dot_ep_gp > cfg.diagnostics.plastic_threshold),
                sigma_eq_mean=flow_diag.get('sigma_eq_mean', 0.0),
                sigma_eq_max=flow_diag.get('sigma_eq_max', 0.0),
                sigma_y_min=flow_diag.get('sigma_y_min', 0.0),
                sigma_y_mean=flow_diag.get('sigma_y_mean', 0.0),
                sigma_y_max=flow_diag.get('sigma_y_max', 0.0),
                sigma_T_min=flow_diag.get('sigma_T_min', 0.0),
                sigma_T_mean=flow_diag.get('sigma_T_mean', 0.0),
                sigma_T_max=flow_diag.get('sigma_T_max', 0.0),
                sigma_Peierls=flow_diag.get('sigma_Peierls', 0.0),
                sigma_eq_over_sigma_y_max=flow_diag.get('sigma_eq_over_sigma_y_max', 0.0),
                yield_frac=flow_diag.get('yield_frac', 0.0),
                flow_dgamma_uncapped_max=flow_diag.get('flow_dgamma_uncapped_max', 0.0),
                flow_dgamma_cap=flow_diag.get('flow_dgamma_cap', 0.0),
                flow_cap_frac=flow_diag.get('flow_cap_frac', 0.0),
                flow_phi_mean=flow_diag.get('flow_phi_mean', 0.0),
                flow_phi_max=flow_diag.get('flow_phi_max', 0.0),
                flow_Gtarget_eV_min=flow_diag.get('flow_Gtarget_eV_min', 0.0),
                flow_Gtarget_eV_mean=flow_diag.get('flow_Gtarget_eV_mean', 0.0),
                flow_Gtarget_eV_max=flow_diag.get('flow_Gtarget_eV_max', 0.0),
                flow_DG0_eV=flow_diag.get('flow_DG0_eV', 0.0),
                flow_DGfloor_eV=flow_diag.get('flow_DGfloor_eV', 0.0),
                flow_vstar_ref_b3=flow_diag.get('flow_vstar_ref_b3', 0.0),
                flow_status_zero_stress_frac=flow_diag.get('flow_status_zero_stress_frac', 0.0),
                flow_status_solved_frac=flow_diag.get('flow_status_solved_frac', 0.0),
                flow_status_floor_limited_frac=flow_diag.get('flow_status_floor_limited_frac', 0.0),
                rtip_mean=np.mean(rtip_state),
                shield_mean=np.mean(shield_state),
                rtip_amp_mean=np.mean(rtip_amp),
                rtip_min=np.min(rtip_state),
                rtip_max=np.max(rtip_state),
                shield_max=np.max(shield_state),
                M_fract_mean=np.mean(M_fract),
                M_fract_max=np.max(M_fract),
                tip_emit_prob_mean=float(np.mean(tip_emit_prob)),
                tip_emit_prob_max=float(np.max(tip_emit_prob)),
                pz_emit_prob_mean=float(np.mean(np.asarray(pz_emit_info.get('P_emit', np.zeros(mesh.nn))))),
                pz_emit_prob_max=float(np.max(np.asarray(pz_emit_info.get('P_emit', np.zeros(mesh.nn))))),
                pz_mobility_prob_mean=float(np.mean(np.asarray(pz_mobility_info.get('P_mobility', np.zeros(mesh.nn))))),
                pz_mobility_prob_max=float(np.max(np.asarray(pz_mobility_info.get('P_mobility', np.zeros(mesh.nn))))),
                pz_mobile_prob_max=float(np.max(np.asarray(pz_mobility_info.get('P_mobile', np.zeros(mesh.nn))))),
                pz_escape_prob_max=float(np.max(np.asarray(pz_mobility_info.get('P_escape', np.zeros(mesh.nn))))),
                pz_store_prob_mean=float(np.mean(np.asarray(pz_mobility_info.get('P_store', np.zeros(mesh.nn))))),
                pz_store_prob_max=float(np.max(np.asarray(pz_mobility_info.get('P_store', np.zeros(mesh.nn))))),
                pz_storage_fraction_mean=float(np.mean(np.asarray(pz_mobility_info.get('storage_fraction', np.ones(mesh.nn))))),
                pz_storage_capacity_mean=float(np.mean(np.asarray(pz_mobility_info.get('storage_capacity', np.ones(mesh.nn))))),
                pz_storage_capacity_min=float(np.min(np.asarray(pz_mobility_info.get('storage_capacity', np.ones(mesh.nn))))),
                pz_source_availability_mean=float(np.mean(np.asarray(pz_emit_info.get('source_availability', np.ones(mesh.nn))))),
                pz_source_availability_min=float(np.min(np.asarray(pz_emit_info.get('source_availability', np.ones(mesh.nn))))),
                pz_rho_source_sat_min=float(np.min(np.asarray(pz_emit_info.get('rho_source_sat', np.full(mesh.nn, np.nan))))),
                pz_rho_source_sat_max=float(np.max(np.asarray(pz_emit_info.get('rho_source_sat', np.full(mesh.nn, np.nan))))),
                pz_multihit_n_emit_mean=float(np.mean(np.asarray(pz_emit_info.get('multihit_n_emit', np.ones(mesh.nn))))),
                pz_multihit_n_emit_max=float(np.max(np.asarray(pz_emit_info.get('multihit_n_emit', np.ones(mesh.nn))))),
                pz_multihit_n_mobility_mean=float(np.mean(np.asarray(pz_mobility_info.get('multihit_n_mobility', np.ones(mesh.nn))))),
                pz_multihit_n_mobility_max=float(np.max(np.asarray(pz_mobility_info.get('multihit_n_mobility', np.ones(mesh.nn))))),
                pz_multihit_spacing_nm_min=float(1e9*np.nanmin(np.asarray(pz_emit_info.get('multihit_spacing_emit', np.full(mesh.nn, np.nan))))),
                pz_multihit_log_suppression_emit_min=float(np.nanmin(np.asarray(pz_emit_info.get('multihit_log_suppression_emit', np.zeros(mesh.nn))))),
                pz_multihit_log_suppression_mobility_min=float(np.nanmin(np.asarray(pz_mobility_info.get('multihit_log_suppression_mobility', np.zeros(mesh.nn))))),
                pz_mobility_hazard_max=float(np.max(np.asarray(pz_mobility_info.get('mobility_hazard', np.zeros(mesh.nn))))),
                pz_mobility_hazard_raw_max=float(np.max(np.asarray(pz_mobility_info.get('mobility_hazard_raw', np.zeros(mesh.nn))))),
                pz_sigma_mobility_eff_max=float(np.max(np.asarray(pz_mobility_info.get('sigma_mobility_eff', np.zeros(mesh.nn))))),
                pz_sigma_tip_eff_max=float(np.max(np.asarray(pz_emit_info.get('sigma_tip_eff', np.zeros(mesh.nn))))),
                pz_sigma_back_max=float(np.max(np.asarray(pz_emit_info.get('sigma_back', np.zeros(mesh.nn))))),
                pz_sigma_back_disl_max=float(np.max(np.asarray(pz_emit_info.get('sigma_back_disl', np.zeros(mesh.nn))))),
                pz_sigma_back_mem_max=float(np.max(np.asarray(pz_emit_info.get('sigma_back_mem', np.zeros(mesh.nn))))),
                pz_sigma_back_crack_max=float(np.max(np.asarray(pz_crack_info.get('tau_back_crack', np.zeros(mesh.nn))))),
                pz_G_shield_max=float(np.max(np.asarray(pz_crack_info.get('G_shield', np.zeros(mesh.nn))))),
                pz_G_stored_release_max=float(np.max(np.asarray(pz_crack_info.get('G_stored_release', np.zeros(mesh.nn))))),
                pz_G_stored_release_p99=float(np.percentile(np.asarray(pz_crack_info.get('G_stored_release', np.zeros(mesh.nn))), 99)),
                pz_e_stored_max=float(np.max(np.asarray(pz_crack_info.get('e_stored', np.zeros(mesh.nn))))),
                pz_Gc_net_min=float(np.min(np.asarray(pz_crack_info.get('Gc_net', np.full(mesh.nn, Gc_eff))))),
                pz_Gc_net_p01=float(np.percentile(np.asarray(pz_crack_info.get('Gc_net', np.full(mesh.nn, Gc_eff))), 1)),
                pz_G_app_max=float(np.max(np.asarray(crack_hazard_info.get('G_app', np.zeros(mesh.nn))))),
                pz_G_eff_max=float(np.max(np.asarray(crack_hazard_info.get('G_eff', np.zeros(mesh.nn))))),
                pz_crack_R_max=float(np.max(np.asarray(crack_hazard_info.get('R_crack', np.zeros(mesh.nn))))),
                pz_crack_R_p99=float(np.percentile(np.asarray(crack_hazard_info.get('R_crack', np.zeros(mesh.nn))), 99)),
                pz_H_eff_drive_max=float(np.max(np.asarray(crack_hazard_info.get('H_eff_drive', np.zeros(mesh.nn))))),
                pz_H_eff_drive_p99=float(np.percentile(np.asarray(crack_hazard_info.get('H_eff_drive', np.zeros(mesh.nn))), 99)),
                pz_front_mask_frac=float(crack_hazard_info.get('front_mask_frac', 0.0)),
                pz_H_eff_masked_max=float(np.max(np.asarray(crack_hazard_info.get('H_eff_masked', np.zeros(mesh.nn))))),
                pz_H_eff_unmasked_max=float(np.max(np.asarray(crack_hazard_info.get('H_eff_drive', np.zeros(mesh.nn))))),
                pz_crack_barrier_min_eV=float(np.nanmin(np.asarray(crack_hazard_info.get('barrier_eV', np.full(mesh.nn, np.inf))))),
                pz_crack_hazard_max=float(np.max(np.asarray(crack_hazard_info.get('hazard', np.zeros(mesh.nn))))),
                pz_crack_prob_mean=float(np.mean(np.asarray(crack_hazard_info.get('P_crack', np.zeros(mesh.nn))))),
                pz_crack_prob_max=float(np.max(np.asarray(crack_hazard_info.get('P_crack', np.zeros(mesh.nn))))),
                pz_crack_hazard_raw_max=float(np.max(np.asarray(crack_hazard_info.get('hazard_raw', np.zeros(mesh.nn))))),
                pz_crack_B_mean=float(np.mean(np.asarray(crack_hazard_info.get('B_crack', np.zeros(mesh.nn))))),
                pz_crack_B_max=float(np.max(np.asarray(crack_hazard_info.get('B_crack', np.zeros(mesh.nn))))),
                pz_emit_B_max=float(np.max(np.asarray(crack_hazard_info.get('B_emit', np.zeros(mesh.nn))))),
                pz_emit_rho_max=float(np.max(np.asarray(crack_hazard_info.get('rho_emit', np.zeros(mesh.nn))))),
                pz_emit_Gshield_max=float(np.max(np.asarray(crack_hazard_info.get('G_shield_emit', np.zeros(mesh.nn))))),
                pz_crack_sigma_tip_max=float(np.max(np.asarray(crack_hazard_info.get('sigma_tip_crack', np.zeros(mesh.nn))))),
                pz_emission_hazard_max=float(np.max(np.asarray(pz_emit_info.get('emission_hazard', np.zeros(mesh.nn))))),
                pz_emission_hazard_raw_max=float(np.max(np.asarray(pz_emit_info.get('emission_hazard_raw', pz_emit_info.get('emission_hazard', np.zeros(mesh.nn)))))),
                pz_drho_emit_max=float(np.max(np.asarray(pz_rho_info.get('drho_emit_gp', np.zeros(mesh.ne))))),
                pz_drho_rec_max=float(np.max(np.asarray(pz_rho_info.get('drho_rec_gp', np.zeros(mesh.ne))))),
                pz_recovery_rate_max=float(np.max(np.asarray(pz_rho_info.get('rho_recovery_rate_gp', np.zeros(mesh.ne))))),
                dwp_norm_front_mean=float(tip_mem_info.get('dwp_norm_front_mean', 0.0)),
                q_tough_mean_front=float(np.sum(tough_weight_node * q_tough_node) / max(np.sum(tough_weight_node), 1e-30)),
                q_tough_mean_all=float(np.mean(q_tough_node)) if q_tough_node.size else 0.0,
                q_tough_p95=float(np.percentile(q_tough_node, 95)) if q_tough_node.size else 0.0,
                q_tough_p99=float(np.percentile(q_tough_node, 99)) if q_tough_node.size else 0.0,
                q_tough_max=float(np.max(q_tough_node)) if q_tough_node.size else 0.0,
                dqtough_max=float(np.max(dqtough_step_node)) if dqtough_step_node.size else 0.0,
                toughening_weight_mean=float(np.mean(tough_weight_node)) if tough_weight_node.size else 0.0,
                toughening_energy_increment=float(dEtough),
                toughening_dissipation_increment=float(dDtough),
            )

            if coh_gate is not None:
                diag.cohesive_gate_max = np.max(coh_gate)
                if sigma_coh > 0:
                    s1_n = project_gp_to_nodes(mesh, sigma1_gp)
                    diag.cohesive_ratio_max = np.max(
                        np.maximum(s1_n, 0) / sigma_coh)

            if cfg.fracture_mode == 'emergent':
                wp_node_snap = project_gp_to_nodes(mesh, wp_gp)
                rho_node_snap = project_gp_to_nodes(mesh, rho_gp)
                Gc_snap, pz_crack_info_snap = _local_Gc_from_state(
                    Gc_eff, wp_node_snap, q_tough_node, ell, cfg,
                    rho_node=rho_node_snap, material=cfg.material,
                    plast_model=plast_model, T=T0, return_terms=True, d=d, mesh=mesh
                )
                diag.Gc_local_mean = float(np.mean(Gc_snap))
                diag.Gc_local_p95 = float(np.percentile(Gc_snap, 95))
                diag.Gc_local_p99 = float(np.percentile(Gc_snap, 99))
                diag.Gc_local_max = float(np.max(Gc_snap))
                # Mean Gc at crack front (d between 0.1 and 0.9)
                front = (d > 0.1) & (d < 0.9)
                if np.any(front):
                    diag.Gc_local_mean_front = float(np.mean(Gc_snap[front]))

            hist.add_step(diag)

            # --- Step-end progress summary (console + progress.json) ---
            _step_dt = time.time() - _prog['t_step0']
            _prog['step_times'].append(_step_dt)
            if _prog_on:
                _wp_ratio = 100.0 * float(diag.Wp) / max(abs(float(diag.Wext_abs)), 1e-30)
                _mean_dt = float(np.mean(_prog['step_times'][-10:]))
                _eta_s = _mean_dt * (n_steps - step)
                if _show_step:
                    print(f"  [T={T0:.0f}K] step {step}/{n_steps} done {_step_dt:5.1f}s "
                          f"(eta {_eta_s/60:4.1f}m) | F={float(diag.Freact):.3g} "
                          f"d_frac={float(diag.d_frac):.3f} rho_max={float(diag.rho_max):.2e} "
                          f"Wp/Wext={_wp_ratio:.0f}% Wtip/Wext={100.0*float(getattr(diag,'Wp_tip',0.0))/max(abs(float(diag.Wext_abs)),1e-30):.1f}% ok={int(float(diag.energy_balance_ok))} "
                          f"substeps={int(getattr(diag,'thermo_substeps',0) or 0)} "
                          f"B_max={float(getattr(diag,'pz_crack_B_max',0.0) or 0.0):.2g}"
                          + (f" CONSUME n={int(_consume_n_step)} dq_max={_consume_dq_step:.2g}" if _consume_n_step > 0 else ""),
                          flush=True)
                _write_progress(outdir, {
                    'state': 'running', 'phase': 'step_end', 'temperature_K': float(T0),
                    'step': int(step), 'n_steps': int(n_steps),
                    'frac_done': round(step / max(n_steps, 1), 4),
                    'step_time_s': round(_step_dt, 2),
                    'mean_step_time_s': round(_mean_dt, 2),
                    'eta_remaining_s': round(_eta_s, 1),
                    'run_elapsed_s': round(time.time() - _prog_t_run0, 1),
                    'Freact': float(diag.Freact), 'd_frac': float(diag.d_frac),
                    'rho_max': float(diag.rho_max),
                    'Wp_over_Wext_pct': round(_wp_ratio, 1),
                    'energy_balance_ok': int(float(diag.energy_balance_ok)),
                    'thermo_substeps': int(getattr(diag, 'thermo_substeps', 0) or 0),
                    'pz_crack_B_max': float(getattr(diag, 'pz_crack_B_max', 0.0) or 0.0),
                    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
                })

            invalid_reason = _invalid_state_reason(cfg, diag, Gc_eff)
            if invalid_reason is not None:
                print(f"    Stop-on-invalid at step {step}: {invalid_reason}", flush=True)
                if _prog_on:
                    _write_progress(outdir, {
                        'state': 'stopped', 'reason': f'invalid: {invalid_reason}',
                        'temperature_K': float(T0), 'step': int(step), 'n_steps': int(n_steps),
                        'run_elapsed_s': round(time.time() - _prog_t_run0, 1),
                        'timestamp': time.strftime('%Y-%m-%d %H:%M:%S')})
                # Keep the final invalid state for inspection, but do not
                # continue accumulating post-failure garbage.
                if cfg.diagnostics.save_fields:
                    hist.d_fields[step] = d.copy()
                    hist.u_fields[step] = u.copy()
                    hist.rho_fields[step] = project_gp_to_nodes(mesh, rho_gp)
                    hist.rtip_fields[step] = rtip_state.copy()
                    hist.shield_fields[step] = shield_state.copy()
                    hist.wp_fields[step] = project_gp_to_nodes(mesh, wp_gp)
                    if cfg.fracture_mode == 'emergent':
                        hist.Gc_fields[step] = Gc_snap.copy()
                    hist.M_fields[step] = M_fract.copy()
                    _stash_stress_diags(hist, step, mesh, s1_node=s1_node,
                                        crack_hazard_info=crack_hazard_info, B_crack_node=B_crack_node)
                break

            # Save field snapshots periodically
            if cfg.diagnostics.save_fields and step % max(cfg.diagnostics.save_every, 1) == 0:
                hist.d_fields[step] = d.copy()
                hist.u_fields[step] = u.copy()
                hist.rho_fields[step] = project_gp_to_nodes(mesh, rho_gp)
                hist.rtip_fields[step] = rtip_state.copy()
                hist.shield_fields[step] = shield_state.copy()
                hist.wp_fields[step] = project_gp_to_nodes(mesh, wp_gp)
                if cfg.fracture_mode == 'emergent':
                    hist.Gc_fields[step] = Gc_snap.copy()
                hist.M_fields[step] = M_fract.copy()
                _stash_stress_diags(hist, step, mesh, s1_node=s1_node,
                                    crack_hazard_info=crack_hazard_info, B_crack_node=B_crack_node)

            # Update previous-step thermodynamic references only after the
            # current coupled increment has been accepted into history.
            Wext_prev_cum = Wext_cum
            Wext_top_prev_cum = Wext_top_cum
            Wext_abs_prev_cum = Wext_abs_cum
            Wp_prev_cum = Wp_cum
            Dp_eff_prev_cum = Dp_eff_cum
            Etough_prev = Etough_cum
            Dtough_prev = Dtough_cum
            Uel_prev = Uel
            Uel_drive_prev = Uel_drive
            Epf_prev = Epf
            Emem_prev = Emem_cum
            Dmem_prev = Dmem_cum
            Dfrac_prev = Dfrac_cum

            # --- Auto-stop ---
            if cfg.auto_stop.enabled:
                Fabs = abs(Fstep)
                Fmax = max(Fmax, Fabs)

                if step >= cfg.auto_stop.min_step and Fmax > 0:
                    if Fabs < cfg.auto_stop.drop_factor * Fmax:
                        n_quiet += 1
                    else:
                        n_quiet = 0

                    if n_quiet >= cfg.auto_stop.n_quiet_required:
                        print(f"    Auto-stop at step {step} "
                              f"(Fmax={Fmax:.3g}, F={Fstep:.3g})", flush=True)
                        if _prog_on:
                            _write_progress(outdir, {
                                'state': 'stopped', 'reason': 'auto_stop_force_drop',
                                'temperature_K': float(T0), 'step': int(step), 'n_steps': int(n_steps),
                                'run_elapsed_s': round(time.time() - _prog_t_run0, 1),
                                'timestamp': time.strftime('%Y-%m-%d %H:%M:%S')})
                        break

        # --- End of step loop ---
        if _prog_on:
            _write_progress(outdir, {
                'state': 'temperature_done', 'temperature_K': float(T0),
                'last_step': int(step), 'n_steps': int(n_steps),
                'run_elapsed_s': round(time.time() - _prog_t_run0, 1),
                'mean_step_time_s': (round(float(np.mean(_prog['step_times'])), 2)
                                     if _prog['step_times'] else None),
                'timestamp': time.strftime('%Y-%m-%d %H:%M:%S')})
        # Always save the final state, even when it is not on the periodic cadence.
        if cfg.diagnostics.save_fields and step not in hist.d_fields:
            hist.d_fields[step] = d.copy()
            hist.u_fields[step] = u.copy()
            hist.rho_fields[step] = project_gp_to_nodes(mesh, rho_gp)
            hist.rtip_fields[step] = rtip_state.copy()
            hist.shield_fields[step] = shield_state.copy()
            hist.wp_fields[step] = project_gp_to_nodes(mesh, wp_gp)
            if cfg.fracture_mode == 'emergent':
                wp_node_snap = project_gp_to_nodes(mesh, wp_gp)
                rho_node_snap_hist = project_gp_to_nodes(mesh, rho_gp)
                hist.Gc_fields[step] = _local_Gc_from_state(
                    Gc_eff, wp_node_snap, q_tough_node, ell, cfg,
                    rho_node=rho_node_snap_hist, material=cfg.material,
                    plast_model=plast_model, T=T0, d=d, mesh=mesh
                )
            hist.M_fields[step] = M_fract.copy()
            _stash_stress_diags(hist, step, mesh, s1_node=s1_node,
                                crack_hazard_info=crack_hazard_info, B_crack_node=B_crack_node)

        hist.n_steps_used = step
        hist.KJ_final = KJ_domain
        hist.J_final = J_domain

        elapsed = time.time() - t_start
        print(f"\n  Results for T = {T0:.0f} K:")
        print(f"    Steps used: {step}")
        print(f"    Gc_eff = {Gc_eff:.4g} J/m²")
        print(f"    Kc_input = {Kc_input/1e6:.3f} MPa·√m")
        print(f"    KJ (domain integral) = {KJ_domain/1e6:.3f} MPa·√m")
        print(f"    KJ (global energy)   = {KJ_global/1e6:.3f} MPa·√m")
        print(f"    Elapsed: {elapsed:.1f} s")

        # Save history, scalar diagnostics, and field snapshots.
        if cfg.save_to_disk:
            tag = f"{T0:04.0f}K"
            T_outdir = os.path.join(outdir, tag)
            os.makedirs(T_outdir, exist_ok=True)
            save_history(hist, os.path.join(T_outdir, f'history_{tag}.npz'),
                        mesh_nodes=mesh.nodes, mesh_elems=mesh.elems)
            save_step_table(hist, os.path.join(T_outdir, f'step_diagnostics_{tag}.csv'))
            save_summary_json(hist, os.path.join(T_outdir, f'summary_{tag}.json'))
            if cfg.diagnostics.save_field_pngs:
                plot_field_snapshots(hist, T_outdir, mesh.nodes, mesh.elems,
                                     max_cols=cfg.diagnostics.max_snapshot_cols)

        # Diagnostic scalar time-series plots
        if cfg.diagnostics.make_plots:
            T_plotdir = os.path.join(outdir, f'plots_{T0:04.0f}K')
            plot_diagnostics(hist, T_plotdir)

        results[T0] = hist

    # --- Summary across temperatures ---
    print("\n" + "=" * 80)
    print("  SUMMARY: Toughness vs Temperature")
    if cfg.fracture_mode == 'emergent':
        print("  Mode: EMERGENT — Gc is constant, apparent Kc from plasticity-fracture competition")
    print("=" * 80)
    gc_hdr = "  {'Gc_max':>8}" if cfg.fracture_mode == 'emergent' else ""
    gc_unit = "  {'(J/m²)':>8}" if cfg.fracture_mode == 'emergent' else ""
    print(f"  {'T (K)':>8}  {'Gc_eff':>8}  {'Kc_in':>8}  "
          f"{'KJ_sel':>8}  {'KJ_fin':>8}  "
          f"{'Wp/Wext':>8}  {'rho_max':>10}  {'steps':>5}"
          + (f"  {'Gc_max':>8}" if cfg.fracture_mode == 'emergent' else "")
          + f"  {'failure_mode':>24}")
    print(f"  {'':>8}  {'(J/m²)':>8}  {'(MPa√m)':>8}  "
          f"{'(MPa√m)':>8}  {'(MPa√m)':>8}  "
          f"{'(%)':>8}  {'(m⁻²)':>10}  {'':>5}"
          + (f"  {'(J/m²)':>8}" if cfg.fracture_mode == 'emergent' else "")
          + f"  {'':>24}")
    print("  " + "-" * (75 + (10 if cfg.fracture_mode == 'emergent' else 0)))

    from .diagnostics import history_summary
    for T in sorted(results.keys()):
        h = results[T]
        summ = history_summary(h)
        KJ_sel = summ.get('KJ_selected_MPa_sqrt_m', 0.0)
        Wext_arr = h.get_array('Wext')
        Wp_arr = h.get_array('Wp')
        Wext_final = Wext_arr[-1] if len(Wext_arr) > 0 else 1e-30
        Wp_final = Wp_arr[-1] if len(Wp_arr) > 0 else 0
        Wp_ratio = 100 * Wp_final / max(Wext_final, 1e-30)
        rho_max_arr = h.get_array('rho_max')
        rho_peak = np.max(rho_max_arr) if len(rho_max_arr) > 0 else 0
        # Emergent mode: show max local Gc
        Gc_max_arr = h.get_array('Gc_local_max') if cfg.fracture_mode == 'emergent' else np.array([0])
        Gc_max = np.max(Gc_max_arr) if len(Gc_max_arr) > 0 else 0
        Gc_str = f"  {Gc_max:8.1f}" if cfg.fracture_mode == 'emergent' else ""
        mode = str(summ.get('failure_mode', 'unknown'))
        print(f"  {T:8.0f}  {h.Gc_eff:8.2f}  {h.Kc_input/1e6:8.3f}  "
              f"{KJ_sel:8.3f}  {h.KJ_final/1e6:8.3f}  "
              f"{Wp_ratio:8.1f}  {rho_peak:10.2e}  {h.n_steps_used:5d}{Gc_str}  {mode:>24}")

    if cfg.save_to_disk:
        save_results_summary_table(results, os.path.join(outdir, 'summary_by_temperature.csv'))

    # Summary plot
    if cfg.diagnostics.make_plots:
        plot_toughness_vs_T(results, outdir)

    return results



def _load_exp_floor_plastic_barrier(cfg, json_path, system_name, energy_scale=0.2,
                                    entropy_scale=1.0, stress_scale=1.0):
    """Load one EXP_floor barrier parameter set from BarrierModel_Export.json."""
    if not os.path.exists(json_path):
        raise FileNotFoundError(
            f"Plastic barrier JSON not found: {json_path}. "
            "Pass --plastic-barrier-json /path/to/BarrierModel_Export.json"
        )
    with open(json_path, 'r') as f:
        data = json.load(f)

    models = data.get('models', [])
    choices = [m.get('system', '') for m in models]
    match = None
    for m in models:
        if str(m.get('system', '')).lower() == str(system_name).lower():
            match = m
            break
    if match is None:
        raise ValueError(
            f"System '{system_name}' not found in {json_path}. Available systems: {choices}"
        )

    par = match.get('parameters', {})
    if str(par.get('model_type', match.get('model_type', ''))).lower() != 'exp_floor':
        raise ValueError(f"Selected barrier system {system_name} is not EXP_floor")

    p = cfg.plasticity_barrier
    p.model_type = 'exp_floor'
    p.exp_system = str(match.get('system', system_name))
    p.exp_Tref_K = float(par.get('Tref_K'))
    p.exp_Tmin_K = float(par.get('Tmin_K', 0.0))
    p.exp_Tmax_K = float(par.get('Tmax_K', 1e99))
    p.exp_G00_eV = float(par.get('G00_eV'))
    p.exp_gT_eV_per_K = float(par.get('gT_eV_per_K'))
    p.exp_sigc0_Pa = float(par.get('sigc0_Pa'))
    p.exp_sT_Pa_per_K = float(par.get('sT_Pa_per_K'))
    p.exp_a = float(par.get('a'))
    p.exp_n = float(par.get('n'))
    p.exp_Gfloor_fraction = float(par.get('Gfloor_fraction', 0.02))
    p.exp_Gfloor_min_eV = float(par.get('Gfloor_min_eV', 1e-4))
    p.exp_Gfloor_max_fraction = float(par.get('Gfloor_max_fraction', 0.95))
    p.exp_energy_scale = float(energy_scale)
    p.exp_entropy_scale = float(entropy_scale)
    p.exp_stress_scale = float(stress_scale)

    implied_S_kB = -p.exp_entropy_scale * p.exp_gT_eV_per_K / 8.617333262e-5
    print("\n>> Loaded EXP_floor plastic barrier")
    print(f"   system        = {p.exp_system}")
    print(f"   energy_scale  = {p.exp_energy_scale:g}")
    print(f"   entropy_scale = {p.exp_entropy_scale:g}  (effective S ≈ {implied_S_kB:.1f} kB)")
    print(f"   stress_scale  = {p.exp_stress_scale:g}")
    print(f"   G0(Tref)      = {p.exp_energy_scale*p.exp_G00_eV:.3g} eV at Tref={p.exp_Tref_K:.1f} K")
    print(f"   sigc(Tref)    = {p.exp_stress_scale*p.exp_sigc0_Pa/1e9:.3g} GPa")



def _available_exp_floor_systems(json_path, exclude_si=True):
    """Return EXP_floor system names from a BarrierModel_Export.json file."""
    with open(json_path, 'r') as f:
        data = json.load(f)
    out = []
    for m in data.get('models', []):
        sysname = str(m.get('system', ''))
        par = m.get('parameters', {})
        mtype = str(par.get('model_type', m.get('model_type', ''))).lower()
        if mtype != 'exp_floor':
            continue
        if exclude_si and sysname.lower().startswith('si'):
            continue
        out.append(sysname)
    return out


def _safe_system_tag(system_name):
    """Filesystem-safe tag for system-specific output directories."""
    return (str(system_name).replace('[','_').replace(']','')
            .replace('/','_').replace(' ','_').replace('.','p'))



def _safe_float_tag(x):
    try:
        val=float(x)
    except Exception:
        return str(x).replace('.', 'p')
    txt=f"{val:.3e}" if (abs(val)>=1000 or (abs(val)>0 and abs(val)<0.01)) else f"{val:g}"
    return txt.replace('-', 'm').replace('+','').replace('.', 'p')

def _plastic_barrier_system_request_is_all(name):
    return str(name).lower() in ('all_non_si','all_nonsi','all_metals','all','all_including_si')

def _systems_from_request(json_path, request):
    r=str(request).lower()
    if r in ('all_non_si','all_nonsi','all_metals'):
        return _available_exp_floor_systems(json_path, exclude_si=True)
    if r in ('all','all_including_si'):
        return _available_exp_floor_systems(json_path, exclude_si=False)
    return [request]

def main():
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description='Arrhenius phase-field fracture simulation'
    )
    parser.add_argument(
        '--preset', choices=['dbtt', 'ceramic', 'cohesive', 'emergent'],
        default='emergent',
        help='Fracture physics preset (default: emergent)'
    )
    parser.add_argument(
        '--temperatures', nargs='+', type=float,
        default=None,
        help='Temperature list in K (e.g., 500 700 900 1100)'
    )
    parser.add_argument(
        '--steps', type=int, default=None,
        help='Number of load steps'
    )
    parser.add_argument('--dU-top', type=float, default=None,
                        help='Override load increment per step [m]. Useful for resolving stiff crack-tip plasticity transients.')
    parser.add_argument('--load-subcycle', type=int, default=None,
                        help='Refine the imposed load path by multiplying n_steps and dividing dU_top by this factor. This is an explicit load-step retry/subcycling surrogate for stiff process-zone cases.')
    parser.add_argument('--nx', type=int, default=None,
                        help='Mesh divisions in x; nodes = nx+1')
    parser.add_argument('--ny', type=int, default=None,
                        help='Mesh divisions in y; nodes = ny+1')
    parser.add_argument('--mesh-jitter', type=float, default=None,
                        help='Interior node jitter fraction, 0 disables jitter')
    parser.add_argument('--ell-factor', type=float, default=None,
                        help='Use ell = ell_factor*hbar unless --ell is supplied')
    parser.add_argument('--ell', type=float, default=None,
                        help='Fixed physical phase-field length ell [m]; overrides --ell-factor')
    parser.add_argument(
        '--no-plots', action='store_true',
        help='Disable diagnostic plots'
    )
    parser.add_argument(
        '--memory-mode', choices=['off', 'weak_stage1', 'stage1'],
        default=None,
        help='Reduced crack-tip memory ablation mode'
    )
    parser.add_argument('--tip-memory-gain', type=float, default=None,
                        help='Override crack-tip memory state_gain')
    parser.add_argument('--tip-M-max', type=float, default=None,
                        help='Override maximum local crack-tip amplification M_max')
    parser.add_argument('--tip-amp-max', type=float, default=None,
                        help='Override r_tip amplification cap amp_max')
    parser.add_argument('--tip-shield-max', type=float, default=None,
                        help='Override maximum shielding z_shield')
    parser.add_argument('--tip-blunt-work', type=float, default=None,
                        help='Override blunting gain from normalized plastic work')
    parser.add_argument('--tip-sharpen-damage', type=float, default=None,
                        help='Override sharpening gain from damage advance')
    parser.add_argument('--tip-drive-exponent', type=float, default=None,
                        help='Override exponent for M_tip coupling to PF damage drive')
    parser.add_argument('--no-tip-drive-coupling', action='store_true',
                        help='Disable M_tip coupling to the phase-field damage drive')
    parser.add_argument('--enable-kinetic-damage-drive', action='store_true',
                        help='Enable extra Model-A kinetic drive in addition to variational AT2')
    parser.add_argument('--pf-damage-cap', type=float, default=None,
                        help='Maximum phase-field damage increment per stagger iteration')
    parser.add_argument('--no-auto-stop', action='store_true',
                        help='Disable force-drop auto-stop so full step history is saved')
    parser.add_argument('--save-every', type=int, default=None,
                        help='Field snapshot cadence in load steps')
    parser.add_argument('--no-progress', action='store_true',
                        help='Disable live progress monitoring (console heartbeat + progress.json)')
    parser.add_argument('--progress-interval', type=float, default=None,
                        help='Wall-clock seconds between intra-step heartbeats (default 15)')
    parser.add_argument('--progress-every', type=int, default=None,
                        help='Console step-summary cadence in steps (default 1; progress.json is written every step regardless)')
    parser.add_argument('--output-dir', type=str, default=None,
                        help='Output directory for this diagnostic run')
    parser.add_argument('--rho-cap', type=float, default=None,
                        help='Maximum dislocation density [m^-2]')
    parser.add_argument('--dot-ep-max', type=float, default=None,
                        help='Maximum plastic strain rate [1/s]')
    parser.add_argument('--plastic-update-mode', choices=['explicit_rate', 'flow_stress'], default=None,
                        help='Plasticity update mode: explicit_rate uses Arrhenius rate over pseudo-time; flow_stress uses rate-dependent yield/return')
    parser.add_argument('--flow-epsdot-ref', type=float, default=None,
                        help='Reference strain rate [1/s] used for Arrhenius flow-stress inversion')
    parser.add_argument('--peierls-autocalibrate', action='store_true',
                        help='At startup, solve the additive Peierls-floor enthalpy so sigma_Peierls(T_cal)=floor-min with a physical entropy (per-material, no tuning)')
    parser.add_argument('--peierls-floor-min-MPa', type=float, default=None,
                        help='Target Peierls stress floor at the hottest operating T (default 1 MPa)')
    parser.add_argument('--peierls-cal-T', type=float, default=None,
                        help='Calibration temperature for the Peierls floor (default: max requested T)')
    parser.add_argument('--peierls-S-kB', type=float, default=None,
                        help='Physical activation entropy for the Peierls floor [units of kB], e.g. -30; ~ -37 is athermal')
    parser.add_argument('--peierls-v0-b3', type=float, default=None,
                        help='Peierls activation volume [units of b^3] (default 5)')
    parser.add_argument('--taylor-multihit', action='store_true',
                        help='Enable correlated multi-hit Taylor renewal (fixes high-density softening; m=1 recovers independent model)')
    parser.add_argument('--taylor-corr-rho-c', type=float, default=None,
                        help='Density where forest spacing ~ correlation length (n_c=1), [1/m^2], default 1e14')
    parser.add_argument('--taylor-renewal-time', type=float, default=None,
                        help='Correlated-segment renewal time t_c [s], default 1e-9')
    parser.add_argument('--taylor-m-max', type=float, default=None,
                        help='Maximum cooperative hit number at high density (default 5)')
    parser.add_argument('--taylor-m-exponent', type=float, default=None,
                        help='Sharpness p of the m(rho) crossover (default 1)')
    parser.add_argument('--thermo-consistency-mode', choices=['off', 'onsager', 'time_cone'], default=None,
                        help='Thermodynamic plasticity coupling: off, Onsager dissipative flow, or hazard/time-cone acceptance')
    parser.add_argument('--thermo-event-strain', type=float, default=None,
                        help='Equivalent plastic strain increment associated with one time-cone hazard event')
    parser.add_argument('--thermo-onsager-max-fraction', type=float, default=None,
                        help='Maximum fraction of local relaxation distance allowed per Onsager update')
    parser.add_argument('--thermo-work-mode', choices=['avg_stress', 'yield_stress'], default=None,
                        help='Plastic work accounting: accepted average stress path or yield/flow stress dissipation')
    parser.add_argument('--thermo-adaptive-substepping', action='store_true',
                        help='Subcycle stiff Arrhenius/Onsager plasticity updates using hazard/strain criteria')
    parser.add_argument('--thermo-max-substeps', type=int, default=None,
                        help='Maximum internal plasticity substeps per stagger')
    parser.add_argument('--thermo-max-dep-increment', type=float, default=None,
                        help='Maximum accepted equivalent plastic strain per internal substep')
    parser.add_argument('--thermo-max-hazard-increment', type=float, default=None,
                        help='Maximum time-cone hazard increment per internal substep')
    parser.add_argument('--no-thermo-energy-audit', action='store_true',
                        help='Disable incremental thermodynamic energy audit classification')
    parser.add_argument('--thermo-energy-rel-tol', type=float, default=None,
                        help='Relative tolerance for incremental energy-balance audit')
    parser.add_argument('--memory-energetics', action='store_true',
                        help='Enable explicit crack-tip memory storage/dissipation audit')
    parser.add_argument('--no-memory-energetics', action='store_true',
                        help='Disable explicit crack-tip memory storage/dissipation audit')
    parser.add_argument('--memory-energy-r-coeff', type=float, default=None,
                        help='Quadratic storage coefficient for tip-radius memory')
    parser.add_argument('--memory-energy-z-coeff', type=float, default=None,
                        help='Quadratic storage coefficient for shielding memory')
    parser.add_argument('--memory-dissipation-r-coeff', type=float, default=None,
                        help='Rate-independent dissipation coefficient for tip-radius memory changes')
    parser.add_argument('--memory-dissipation-z-coeff', type=float, default=None,
                        help='Rate-independent dissipation coefficient for shielding memory changes')
    parser.add_argument('--plastic-H0-eV', type=float, default=None,
                        help='Override plasticity H0 at sigma0 in eV for diagnostic single-case runs')
    parser.add_argument('--plastic-vstar-b3', type=float, default=None,
                        help='Override plasticity v* at sigma0 in units of b^3 for diagnostic single-case runs')
    parser.add_argument('--plastic-barrier-model', choices=['rational_Hv', 'exp_floor'], default=None,
                        help='Plastic barrier model. exp_floor imports full DeltaG(sigma,T) from JSON')
    parser.add_argument('--plastic-barrier-json', type=str, default=None,
                        help='Path to BarrierModel_Export.json for --plastic-barrier-model exp_floor')
    parser.add_argument('--plastic-barrier-system', type=str, default='W[100]',
                        help='System name inside BarrierModel_Export.json, e.g. W[100], Ta[111], Cu; use all_non_si to run each non-Si system')
    parser.add_argument('--list-plastic-barrier-systems', action='store_true',
                        help='Print available EXP_floor systems in the JSON and exit')
    parser.add_argument('--plastic-barrier-scale', type=float, default=0.2,
                        help='Taylor enthalpic barrier scale relative to nanopillar nucleation')
    parser.add_argument('--plastic-exp-entropy-scale', type=float, default=1.0,
                        help='Scale for fitted gT entropy slope; 1 preserves roughly -40 kB')
    parser.add_argument('--plastic-exp-stress-scale', type=float, default=1.0,
                        help='Scale for sigc(T) in EXP_floor plastic model')
    parser.add_argument('--plastic-exp-vmax-b3', type=float, default=None,
                        help='Cap for derived local activation volume from -dG/dsigma, in b^3')
    parser.add_argument('--Gc0', type=float, default=None,
                        help='Override intrinsic/emergent baseline fracture energy Gc0_athermal [J/m^2]')
    parser.add_argument('--Gc-list', nargs='+', type=float, default=None,
                        help='Run a baseline fracture-energy sweep over these Gc0_athermal values [J/m^2]')
    parser.add_argument('--wp-gc-coupling-mode', choices=['off', 'direct', 'state'], default=None,
                        help='Plastic-work toughening: off, legacy direct Gc=Gc0+eta*Wp*ell, or thermodynamic retained q_Gc state')
    parser.add_argument('--wp-gc-efficiency', type=float, default=None,
                        help='Efficiency converting accepted plastic work near the crack into retained Gc state')
    parser.add_argument('--gc-local-cap-factor', type=float, default=None,
                        help='Upper bound Gc_local <= cap_factor*Gc0 for local process-zone toughening')
    parser.add_argument('--toughening-storage-coeff', type=float, default=None,
                        help='Storage coefficient for retained q_Gc toughening state')
    parser.add_argument('--toughening-dissipation-coeff', type=float, default=None,
                        help='Dissipation coefficient for changes in retained q_Gc toughening state')
    parser.add_argument('--toughening-relax-per-step', type=float, default=None,
                        help='Optional fractional relaxation of retained q_Gc per load step')
    parser.add_argument('--no-toughening-front-only', action='store_true',
                        help='Allow wake weighting in addition to crack-front localization for q_Gc toughening')
    parser.add_argument('--no-toughening-energy-audit', action='store_true',
                        help='Do not partition retained q_Gc storage/dissipation in the energy audit')
    parser.add_argument('--process-zone-mode', choices=['off', 'on'], default=None,
                        help='Enable/disable physical process-zone emission/recovery kinetics')
    parser.add_argument('--pz-qgc-driver', choices=['plastic_work', 'emission', 'mixed'], default=None,
                        help='Driver for q_Gc retained toughening state')
    parser.add_argument('--pz-emission-eta0', type=float, default=None,
                        help='Attempt frequency for crack-tip dislocation emission [1/s]')
    parser.add_argument('--pz-emission-H-scale', type=float, default=None,
                        help='Stress/barrier scale factor used in tip-emission barrier evaluation')
    parser.add_argument('--pz-emission-probability-cap', type=float, default=None,
                        help='Maximum accepted tip-emission probability per load step')
    parser.add_argument('--pz-mobility-enabled', dest='pz_mobility_enabled', action='store_true', default=None,
                        help='Enable separated emitted-dislocation mobility/storage partition')
    parser.add_argument('--no-pz-mobility-enabled', dest='pz_mobility_enabled', action='store_false',
                        help='Disable separated mobility/storage partition')
    parser.add_argument('--pz-mobility-eta0', type=float, default=None,
                        help='Attempt frequency for emitted-dislocation mobility/escape [1/s]')
    parser.add_argument('--pz-mobility-H-scale', type=float, default=None,
                        help='Stress/barrier scale factor for emitted-dislocation mobility')
    parser.add_argument('--pz-mobility-probability-cap', type=float, default=None,
                        help='Maximum mobility probability per load step')
    parser.add_argument('--pz-mobility-backstress-factor', type=float, default=None,
                        help='Multiplier on PZ backstress opposing emitted-dislocation mobility')
    parser.add_argument('--pz-mobility-escape-fraction', type=float, default=None,
                        help='Fraction of mobile emitted dislocations that escape without storage')
    parser.add_argument('--pz-storage-min-fraction', type=float, default=None,
                        help='Minimum fraction of emitted dislocations retained as stored PZ density')
    parser.add_argument('--pz-storage-max-fraction', type=float, default=None,
                        help='Maximum fraction of emitted dislocations retained as stored PZ density')
    parser.add_argument('--pz-storage-backstress-boost', type=float, default=None,
                        help='Extra storage/pile-up fraction created by high PZ rho/backstress')
    parser.add_argument('--pz-source-availability-enabled', dest='pz_source_availability_enabled', action='store_true', default=None,
                        help='Enable source/pile-up availability suppression of tip emission and storage')
    parser.add_argument('--no-pz-source-availability-enabled', dest='pz_source_availability_enabled', action='store_false',
                        help='Disable source/pile-up availability suppression')
    parser.add_argument('--pz-source-rho-sat', type=float, default=None,
                        help='Process-zone source/pile-up saturation density [m^-2]')
    parser.add_argument('--pz-source-rho-sat-fraction', type=float, default=None,
                        help='Fallback source saturation as fraction of rho_cap')
    parser.add_argument('--pz-source-availability-power', type=float, default=None,
                        help='Hill exponent for source availability vs rho/rho_sat')
    parser.add_argument('--pz-source-availability-floor', type=float, default=None,
                        help='Residual source availability at high rho')
    parser.add_argument('--pz-source-backstress-scale', type=float, default=None,
                        help='Reserved multiplier for future source-backstress coupling')
    parser.add_argument('--pz-multihit-enabled', dest='pz_multihit_enabled', action='store_true', default=None,
                        help='Enable correlated multi-hit Arrhenius/Taylor depinning for high-rho process-zone events')
    parser.add_argument('--no-pz-multihit-enabled', dest='pz_multihit_enabled', action='store_false',
                        help='Disable correlated multi-hit process-zone depinning')
    parser.add_argument('--pz-multihit-apply-to', choices=['emission','mobility','both','off'], default=None)
    parser.add_argument('--pz-multihit-path-length-nm', type=float, default=None)
    parser.add_argument('--pz-multihit-path-length-b', type=float, default=None)
    parser.add_argument('--pz-multihit-density-power', type=float, default=None)
    parser.add_argument('--pz-multihit-max-hits', type=int, default=None)
    parser.add_argument('--pz-multihit-min-hits', type=int, default=None)
    parser.add_argument('--pz-rho-source-sat-cap', dest='pz_rho_source_sat_cap', action='store_true', default=None,
                        help='Diagnostic: also hard-cap rho at source saturation density')
    parser.add_argument('--no-pz-rho-source-sat-cap', dest='pz_rho_source_sat_cap', action='store_false',
                        help='Do not hard-cap rho at source saturation density')
    parser.add_argument('--pz-rho-increment-per-event', type=float, default=None,
                        help='Local rho increment [m^-2] for unit tip-emission probability')
    parser.add_argument('--pz-qgc-from-emission-factor', type=float, default=None,
                        help='dq_Gc/Gc0 generated by unit tip-emission probability')
    parser.add_argument('--pz-backstress-model', choices=['sqrt_taylor', 'arrhenius_taylor', 'max'], default=None,
                        help='Process-zone backstress model for dislocation emission: sqrt Taylor, Arrhenius-Taylor inversion, or max of both')
    parser.add_argument('--pz-backstress-alpha', type=float, default=None,
                        help='Taylor back-stress coefficient alpha in alpha*G*b*sqrt(rho)')
    parser.add_argument('--pz-backstress-rate-ref', type=float, default=None,
                        help='Reference strain/slip rate for Arrhenius-Taylor process-zone backstress inversion')
    parser.add_argument('--pz-backstress-scale', type=float, default=None,
                        help='Multiplier on process-zone dislocation backstress')
    parser.add_argument('--pz-memory-backstress-factor', type=float, default=None,
                        help='Shielding-memory contribution to emission/crack effective stress')
    parser.add_argument('--pz-crack-shielding-coeff', type=float, default=None,
                        help='Crack-extension shielding coefficient in G_shield ~ c*tau_back^2/Eprime*ell')
    parser.add_argument('--pz-stored-energy-coeff', type=float, default=None,
                        help='Stored defect energy coefficient in e_stored ~ c*G*b^2*rho*log(...)')
    parser.add_argument('--pz-stored-release-efficiency', type=float, default=None,
                        help='Fraction of process-zone stored defect energy released during crack extension')
    parser.add_argument('--pz-stored-release-cap-factor', type=float, default=None,
                        help='Optional legacy cap G_stored_release <= factor*Gc0; default is disabled/infinite')
    parser.add_argument('--pz-Gc-net-floor-factor', type=float, default=None,
                        help='Local net Gc floor as fraction of intrinsic Gc0 after stored-energy degradation')
    parser.add_argument('--no-pz-crack-shielding', action='store_true',
                        help='Disable crack-extension shielding/backstress contribution from process-zone rho')
    parser.add_argument('--no-pz-stored-energy', action='store_true',
                        help='Disable stored-energy embrittlement/degradation term')
    parser.add_argument('--pz-recovery-model', choices=['arrhenius','climb_diffusion'], default=None,
                        help='Recovery model: explicit Arrhenius or climb/diffusion scale from DislocationConfig')
    parser.add_argument('--pz-recovery-eta0', type=float, default=None,
                        help='Arrhenius recovery/escape prefactor [1/s]')
    parser.add_argument('--pz-recovery-Q-eV', type=float, default=None,
                        help='Arrhenius recovery/escape activation energy [eV]')
    parser.add_argument('--pz-dynamic-recovery-coeff', type=float, default=None,
                        help='Plastic-flow-assisted recovery coefficient multiplying rho*dot_ep')
    parser.add_argument('--pz-crack-advance-memory-erasure', type=float, default=None,
                        help='Multiplier for crack-advance erasure/convective loss of tip memory')
    parser.add_argument('--pz-crack-hazard', choices=['off','on'], default=None,
                        help='Enable/disable Arrhenius crack-growth hazard based on shielded G_eff')
    parser.add_argument('--pz-crack-hazard-model', choices=['arrhenius','eyring'], default=None,
                        help='Crack-growth hazard form')
    parser.add_argument('--pz-crack-hazard-drive', choices=['ratio','resolved_stress'], default=None,
                        help='Crack hazard drive: legacy G_eff/Gc barrier or first-principles resolved-stress fracture barrier')
    parser.add_argument('--pz-crack-first-passage', dest='pz_crack_first_passage', action='store_true', default=None,
                        help='Accumulate B=int(lambda dt) and accept crack advance as B/Btarget approaches one')
    parser.add_argument('--no-pz-crack-first-passage', dest='pz_crack_first_passage', action='store_false',
                        help='Use memoryless per-step crack event probability')
    parser.add_argument('--pz-crack-B-target', type=float, default=None,
                        help='Integrated hazard/action threshold for crack first passage')
    parser.add_argument('--pz-crack-B-cap', type=float, default=None,
                        help='Numerical cap for accumulated crack action B')
    parser.add_argument('--pz-crack-multihit-m', type=float, default=None,
                        help='Cooperative hit count m for correlated cleavage renewal (1 = independent Poisson, default). m~3-5 suppresses sub-threshold flank firing combinatorially -> straight cracks')
    parser.add_argument('--pz-crack-multihit-tau', type=float, default=None,
                        help='Renewal window tau_c [s] for multi-hit cleavage (default 1e-9)')
    parser.add_argument('--pz-crack-advance-mode', choices=['gate','source'], default=None,
                        help="Crack advance coupling: 'gate' (default; P gates advance toward the variational AT2 target with Model-A kinetics) or 'source' (EXPERIMENTAL: fired events break nodes directly; avalanche-prone when local Gc collapses from stored release -- prefer gate)")
    parser.add_argument('--pz-crack-fire-fw-min', type=float, default=None,
                        help='Minimum front weight for a fired node to break in source mode (default 0.5)')
    parser.add_argument('--pz-crack-fire-griffith-frac', type=float, default=None,
                        help='Source mode Griffith gate: fired node breaks only if G_eff >= frac*Gc_net (default 1.0; <=0 disables)')
    parser.add_argument('--plastic-taylor-athermal-alpha', type=float, default=None,
                        help='Athermal Taylor floor coefficient alpha in sigma_T >= alpha*G*b*sqrt(rho) (default 0.2; 0 disables)')
    parser.add_argument('--pz-crack-source-fw-min', type=float, default=None,
                        help='Front-weight floor for source-mode firing (default 0.3)')
    parser.add_argument('--pz-crack-fired-gc-relief', type=float, default=None,
                        help='Gc_local multiplier at fired nodes in source mode (default 0.05)')
    parser.add_argument('--pz-crack-consume-radius', type=float, default=None,
                        help='Consumption footprint radius in units of ell around fired nodes (default 1.5)')
    parser.add_argument('--pz-crack-consume-toughening', type=float, default=None,
                        help='Fraction of local blunting toughening q_blunt consumed per stagger where the cleavage clock has fired (B>=B_target on the front). Unfreezes the Gc_local halo stall; energy re-booked as dissipation. Recommended 1.0; default 0 (off)')
    parser.add_argument('--pz-crack-B-relax-time', type=float, default=None,
                        help='Sub-critical action relaxation time tau_relax [s]: dB/dt = lambda - B/tau_relax, so flank action anneals instead of ratcheting. <=0 disables (default)')
    parser.add_argument('--tip-wake-relax', type=float, default=None,
                        help='Per-step fractional relaxation of M_tip memory (rtip->ref, shield->0) outside the active front; kills the M_tip ratchet web. 0 disables (default)')
    parser.add_argument('--no-pz-crack-stored-release-front-mask', action='store_true',
                        help='Diagnostic ablation: subtract G_stored_release from Gc_net over the whole dislocation cloud (legacy broad behavior) instead of only at the connected crack front')
    parser.add_argument('--pz-crack-resolved-stress-scale', type=float, default=None,
                        help='Multiplier on resolved/amplified near-tip stress used in fracture barrier')
    parser.add_argument('--pz-crack-eta0', type=float, default=None,
                        help='Crack-growth hazard attempt frequency [1/s]')
    parser.add_argument('--pz-crack-H0-eV', type=float, default=None,
                        help='Zero-drive crack-growth barrier [eV]')
    parser.add_argument('--pz-crack-r-pz', type=float, default=None,
                        help='Physical process-zone radius r_pz [m] for de-smearing the resolved tip stress (chi=sqrt(ell/r_pz)); sets emergent brittle toughness. <=0 disables')
    parser.add_argument('--pz-crack-sigma-cap-GPa', type=float, default=None,
                        help='Optional cohesive/spinodal ceiling on the de-smeared tip stress [GPa]')
    parser.add_argument('--pz-crack-emission', action='store_true',
                        help='Enable crack-tip dislocation-emission competition (Rice-Thomson): native exp_floor emission vs cleavage, emission shields the tip')
    parser.add_argument('--pz-emission-rho-max', type=float, default=None,
                        help='Saturation emitted-dislocation density at the tip [1/m^2] (default 1e15)')
    parser.add_argument('--pz-emission-energy-scale', type=float, default=None,
                        help='exp_floor energy scale for the emission barrier (1.0 = native nanopillar)')
    parser.add_argument('--pz-emission-entropy-scale', type=float, default=None,
                        help='exp_floor entropy scale for the emission barrier (DBTT knob; <1 reduces |S|)')
    parser.add_argument('--pz-emission-B-target', type=float, default=None,
                        help='Emission first-passage threshold (default 1.0)')
    parser.add_argument('--frac-S0-neg-kB', type=float, default=None,
                        help='Cleavage barrier entropy magnitude |S_c| in kB (default 3; ~28 = athermal-flat shelf; <28 toughness falls with T)')
    parser.add_argument('--frac-sigma0-S-GPa', type=float, default=None,
                        help='Stress scale of cleavage entropy S_f=-S0_neg*(1+sigma/sigma0_S) [GPa]. Set large (e.g. 1e4) for stress-INDEPENDENT entropy so --frac-S0-neg-kB maps directly onto the regime criterion |S| vs k*ln(eta0*t)~28 kB')
    parser.add_argument('--no-frac-monotone-barrier', action='store_true',
                        help='Diagnostic: disable the monotone-in-stress guard on the cleavage barrier (legacy raw form, which RISES again above ~8 GPa and can lock up overshooting tips)')
    parser.add_argument('--pz-crack-drive-exponent', type=float, default=None,
                        help='Exponent for barrier collapse with G_eff/Gc_net')
    parser.add_argument('--pz-crack-probability-cap', type=float, default=None,
                        help='Maximum crack-growth event probability per damage update')
    parser.add_argument('--pz-crack-drive-scale', type=float, default=None,
                        help='Multiplier on G_eff/Gc_net in crack hazard')
    parser.add_argument('--no-pz-crack-effective-drive', action='store_true',
                        help='Diagnostic ablation: compute crack hazard from G_eff but do not feed G_eff/(2ell) into the AT2 damage target')
    parser.add_argument('--pz-crack-effective-drive-mix', type=float, default=None,
                        help='Mixing factor [0,1] for using shielded G_eff/(2ell) as the crack damage target at the front')
    parser.add_argument('--pz-crack-front-radius-factor', type=float, default=None,
                        help='Radius, in units of ell, over which degraded crack drive is allowed around the active crack front')
    parser.add_argument('--pz-crack-front-grad-threshold', type=float, default=None,
                        help='Relative |grad d| threshold for selecting active crack-front seeds')
    parser.add_argument('--pz-crack-front-wake-dmax', type=float, default=None,
                        help='Suppress degraded crack drive in the fully damaged wake above this d value')
    parser.add_argument('--allow-invalid-soft-tearing', action='store_true',
                        help='Disable Wp/Wext fail-fast while still classifying plastic/soft tearing in summaries')
    parser.add_argument('--max-plastic-strain-increment', type=float, default=None,
                        help='Maximum equivalent plastic strain increment per stagger')
    parser.add_argument('--max-rho-relative-increment', type=float, default=None,
                        help='Maximum fractional rho change per stagger')
    parser.add_argument('--disable-plasticity', action='store_true',
                        help='Diagnostic ablation: mechanics/fracture only, no plastic update')
    parser.add_argument('--freeze-rho', action='store_true',
                        help='Diagnostic ablation: allow plastic strain but hold rho fixed')
    parser.add_argument('--disable-wp-gc-coupling', action='store_true',
                        help='Diagnostic ablation: set emergent Gc_local=Gc, no Wp toughening')
    parser.add_argument('--stop-on-invalid', action='store_true',
                        help='Fail fast when rho, Wp/Wext, Gc_local, KJ, or d_frac becomes invalid')
    parser.add_argument('--invalid-wp-wext-pct', type=float, default=None,
                        help='Stop-on-invalid threshold for 100*Wp/Wext')
    parser.add_argument('--invalid-min-step', type=int, default=None,
                        help='Do not apply stop-on-invalid before this step (default 3)')
    parser.add_argument('--invalid-wext-min', type=float, default=None,
                        help='Minimum Wext before Wp/Wext invalid check is active')
    parser.add_argument('--invalid-sigma-eq-GPa', type=float, default=None,
                        help='Stop-on-invalid threshold for local von Mises stress [GPa]')
    parser.add_argument('--invalid-dep-eq-increment', type=float, default=None,
                        help='Stop-on-invalid threshold for accepted equivalent plastic strain increment')

    args = parser.parse_args()

    if args.list_plastic_barrier_systems:
        if args.plastic_barrier_json is None:
            raise ValueError('--plastic-barrier-json is required for --list-plastic-barrier-systems')
        systems = _available_exp_floor_systems(args.plastic_barrier_json, exclude_si=False)
        print('Available EXP_floor barrier systems:')
        for sysname in systems:
            print(f'  - {sysname}')
        return

    # Build config from preset
    if args.preset == 'ceramic':
        cfg = make_ceramic_config()
    elif args.preset == 'cohesive':
        cfg = make_cohesive_dbtt_config()
    elif args.preset == 'dbtt':
        cfg = make_dbtt_config()
    else:
        cfg = make_emergent_config()

    # Override from CLI
    if args.temperatures:
        cfg.T_list = args.temperatures
    if args.steps:
        cfg.loading.n_steps = args.steps
    if args.dU_top is not None:
        cfg.loading.dU_top = float(args.dU_top)
    if args.load_subcycle is not None:
        fac = max(int(args.load_subcycle), 1)
        if fac > 1:
            cfg.loading.n_steps = int(cfg.loading.n_steps) * fac
            cfg.loading.dU_top = float(cfg.loading.dU_top) / fac
            print(f">> Load-path subcycling: n_steps -> {cfg.loading.n_steps}, dU_top -> {cfg.loading.dU_top:.3e} m (factor {fac})")
    if args.nx is not None:
        cfg.mesh.nx = int(args.nx)
    if args.ny is not None:
        cfg.mesh.ny = int(args.ny)
    if args.mesh_jitter is not None:
        cfg.mesh.jitter = float(args.mesh_jitter)
    if args.ell_factor is not None:
        cfg.mesh.ell_factor = float(args.ell_factor)
    if args.ell is not None:
        cfg.mesh.ell_absolute_m = float(args.ell)
    if args.no_plots:
        cfg.diagnostics.make_plots = False
    if args.output_dir is not None:
        cfg.output_dir = args.output_dir
    if args.Gc0 is not None:
        cfg.phase_field.Gc0_athermal = float(args.Gc0)
    if args.wp_gc_coupling_mode is not None:
        cfg.phase_field.wp_gc_coupling_mode = args.wp_gc_coupling_mode
    if args.wp_gc_efficiency is not None:
        cfg.phase_field.plastic_work_to_Gc_efficiency = float(args.wp_gc_efficiency)
    if args.gc_local_cap_factor is not None:
        cfg.phase_field.Gc_local_cap_factor = float(args.gc_local_cap_factor)
    if args.toughening_storage_coeff is not None:
        cfg.phase_field.toughening_storage_coeff = float(args.toughening_storage_coeff)
    if args.toughening_dissipation_coeff is not None:
        cfg.phase_field.toughening_dissipation_coeff = float(args.toughening_dissipation_coeff)
    if args.toughening_relax_per_step is not None:
        cfg.phase_field.toughening_relax_per_step = float(args.toughening_relax_per_step)
    if args.no_toughening_front_only:
        cfg.phase_field.toughening_front_only = False
    if args.no_toughening_energy_audit:
        cfg.phase_field.toughening_include_in_energy_audit = False
    if args.process_zone_mode is not None:
        cfg.process_zone.enabled = (args.process_zone_mode == 'on')
    if args.pz_qgc_driver is not None:
        cfg.process_zone.qgc_driver = args.pz_qgc_driver
    if args.pz_emission_eta0 is not None:
        cfg.process_zone.emission_eta0 = float(args.pz_emission_eta0)
    if args.pz_emission_H_scale is not None:
        cfg.process_zone.emission_H_scale = float(args.pz_emission_H_scale)
    if args.pz_emission_probability_cap is not None:
        cfg.process_zone.emission_probability_cap = float(args.pz_emission_probability_cap)
    if args.pz_mobility_enabled is not None:
        cfg.process_zone.mobility_enabled = bool(args.pz_mobility_enabled)
    if args.pz_mobility_eta0 is not None:
        cfg.process_zone.mobility_eta0 = float(args.pz_mobility_eta0)
    if args.pz_mobility_H_scale is not None:
        cfg.process_zone.mobility_H_scale = float(args.pz_mobility_H_scale)
    if args.pz_mobility_probability_cap is not None:
        cfg.process_zone.mobility_probability_cap = float(args.pz_mobility_probability_cap)
    if args.pz_mobility_backstress_factor is not None:
        cfg.process_zone.mobility_backstress_factor = float(args.pz_mobility_backstress_factor)
    if args.pz_mobility_escape_fraction is not None:
        cfg.process_zone.mobility_escape_fraction = float(args.pz_mobility_escape_fraction)
    if args.pz_storage_min_fraction is not None:
        cfg.process_zone.storage_min_fraction = float(args.pz_storage_min_fraction)
    if args.pz_storage_max_fraction is not None:
        cfg.process_zone.storage_max_fraction = float(args.pz_storage_max_fraction)
    if args.pz_storage_backstress_boost is not None:
        cfg.process_zone.storage_backstress_boost = float(args.pz_storage_backstress_boost)
    if args.pz_source_availability_enabled is not None:
        cfg.process_zone.source_availability_enabled = bool(args.pz_source_availability_enabled)
    if args.pz_source_rho_sat is not None:
        cfg.process_zone.source_rho_sat = float(args.pz_source_rho_sat)
    if args.pz_source_rho_sat_fraction is not None:
        cfg.process_zone.source_rho_sat_fraction = float(args.pz_source_rho_sat_fraction)
    if args.pz_source_availability_power is not None:
        cfg.process_zone.source_availability_power = float(args.pz_source_availability_power)
    if args.pz_source_availability_floor is not None:
        cfg.process_zone.source_availability_floor = float(args.pz_source_availability_floor)
    if args.pz_source_backstress_scale is not None:
        cfg.process_zone.source_backstress_scale = float(args.pz_source_backstress_scale)
    if args.pz_multihit_enabled is not None:
        cfg.process_zone.multihit_enabled = bool(args.pz_multihit_enabled)
    if args.pz_multihit_apply_to is not None:
        cfg.process_zone.multihit_apply_to = str(args.pz_multihit_apply_to)
    if args.pz_multihit_path_length_nm is not None:
        cfg.process_zone.multihit_path_length_nm = float(args.pz_multihit_path_length_nm)
    if args.pz_multihit_path_length_b is not None:
        cfg.process_zone.multihit_path_length_b = float(args.pz_multihit_path_length_b)
    if args.pz_multihit_density_power is not None:
        cfg.process_zone.multihit_density_power = float(args.pz_multihit_density_power)
    if args.pz_multihit_max_hits is not None:
        cfg.process_zone.multihit_max_hits = int(args.pz_multihit_max_hits)
    if args.pz_multihit_min_hits is not None:
        cfg.process_zone.multihit_min_hits = int(args.pz_multihit_min_hits)
    if args.pz_rho_source_sat_cap is not None:
        cfg.process_zone.rho_source_saturation_cap_enabled = bool(args.pz_rho_source_sat_cap)
    if args.pz_rho_increment_per_event is not None:
        cfg.process_zone.rho_increment_per_event = float(args.pz_rho_increment_per_event)
    if args.pz_qgc_from_emission_factor is not None:
        cfg.process_zone.qgc_from_emission_factor = float(args.pz_qgc_from_emission_factor)
    if args.pz_backstress_model is not None:
        cfg.process_zone.backstress_model = args.pz_backstress_model
    if args.pz_backstress_alpha is not None:
        cfg.process_zone.backstress_alpha = float(args.pz_backstress_alpha)
    if args.pz_backstress_rate_ref is not None:
        cfg.process_zone.backstress_rate_ref = float(args.pz_backstress_rate_ref)
    if args.pz_backstress_scale is not None:
        cfg.process_zone.backstress_scale = float(args.pz_backstress_scale)
    if args.pz_memory_backstress_factor is not None:
        cfg.process_zone.memory_backstress_factor = float(args.pz_memory_backstress_factor)
    if args.pz_crack_shielding_coeff is not None:
        cfg.process_zone.crack_shielding_coeff = float(args.pz_crack_shielding_coeff)
    if args.pz_stored_energy_coeff is not None:
        cfg.process_zone.stored_energy_coeff = float(args.pz_stored_energy_coeff)
    if args.pz_stored_release_efficiency is not None:
        cfg.process_zone.stored_energy_release_efficiency = float(args.pz_stored_release_efficiency)
    if args.pz_stored_release_cap_factor is not None:
        cfg.process_zone.stored_energy_release_cap_factor = float(args.pz_stored_release_cap_factor)
    if args.pz_Gc_net_floor_factor is not None:
        cfg.process_zone.Gc_net_floor_factor = float(args.pz_Gc_net_floor_factor)
    if args.no_pz_crack_shielding:
        cfg.process_zone.crack_shielding_enabled = False
    if args.no_pz_stored_energy:
        cfg.process_zone.stored_energy_enabled = False
    if args.pz_recovery_model is not None:
        cfg.process_zone.recovery_model = args.pz_recovery_model
    if args.pz_recovery_eta0 is not None:
        cfg.process_zone.recovery_eta0 = float(args.pz_recovery_eta0)
    if args.pz_recovery_Q_eV is not None:
        cfg.process_zone.recovery_Q_eV = float(args.pz_recovery_Q_eV)
    if args.pz_dynamic_recovery_coeff is not None:
        cfg.process_zone.dynamic_recovery_coeff = float(args.pz_dynamic_recovery_coeff)
    if args.pz_crack_advance_memory_erasure is not None:
        cfg.process_zone.crack_advance_memory_erasure = float(args.pz_crack_advance_memory_erasure)
    if args.pz_crack_hazard is not None:
        cfg.process_zone.crack_hazard_enabled = (args.pz_crack_hazard == 'on')
    if args.pz_crack_hazard_model is not None:
        cfg.process_zone.crack_hazard_model = args.pz_crack_hazard_model
    if args.pz_crack_hazard_drive is not None:
        cfg.process_zone.crack_hazard_drive = args.pz_crack_hazard_drive
    if getattr(args, 'pz_crack_r_pz', None) is not None:
        cfg.process_zone.crack_process_zone_r_pz_m = float(args.pz_crack_r_pz)
    if getattr(args, 'pz_crack_sigma_cap_GPa', None) is not None:
        cfg.process_zone.crack_sigma_cap_Pa = float(args.pz_crack_sigma_cap_GPa) * 1e9
    if getattr(args, 'pz_crack_emission', False):
        cfg.process_zone.crack_emission_enabled = True
    if getattr(args, 'pz_emission_rho_max', None) is not None:
        cfg.process_zone.crack_emission_rho_max = float(args.pz_emission_rho_max)
    if getattr(args, 'pz_emission_energy_scale', None) is not None:
        cfg.process_zone.crack_emission_energy_scale = float(args.pz_emission_energy_scale)
    if getattr(args, 'pz_emission_entropy_scale', None) is not None:
        cfg.process_zone.crack_emission_entropy_scale = float(args.pz_emission_entropy_scale)
    if getattr(args, 'pz_emission_B_target', None) is not None:
        cfg.process_zone.crack_emission_B_target = float(args.pz_emission_B_target)
    if getattr(args, 'frac_S0_neg_kB', None) is not None:
        cfg.fracture_barrier.S0_neg_kB = float(args.frac_S0_neg_kB)
    if getattr(args, 'frac_sigma0_S_GPa', None) is not None:
        cfg.fracture_barrier.sigma0_S_GPa = float(args.frac_sigma0_S_GPa)
    if getattr(args, 'no_frac_monotone_barrier', False):
        cfg.fracture_barrier.monotone_stress = False
    if args.pz_crack_first_passage is not None:
        cfg.process_zone.crack_first_passage = bool(args.pz_crack_first_passage)
    if args.pz_crack_B_target is not None:
        cfg.process_zone.crack_B_target = float(args.pz_crack_B_target)
    if getattr(args, 'no_pz_crack_stored_release_front_mask', False):
        cfg.process_zone.crack_stored_release_front_masked = False
    if args.pz_crack_B_cap is not None:
        cfg.process_zone.crack_B_cap = float(args.pz_crack_B_cap)
    if args.pz_crack_multihit_m is not None:
        cfg.process_zone.crack_multihit_m = float(args.pz_crack_multihit_m)
    if args.pz_crack_multihit_tau is not None:
        cfg.process_zone.crack_multihit_tau_s = float(args.pz_crack_multihit_tau)
    if args.pz_crack_B_relax_time is not None:
        cfg.process_zone.crack_B_relax_time_s = float(args.pz_crack_B_relax_time)
    if getattr(args, 'pz_crack_advance_mode', None) is not None:
        cfg.process_zone.crack_advance_mode = str(args.pz_crack_advance_mode)
    if getattr(args, 'pz_crack_fire_fw_min', None) is not None:
        cfg.process_zone.crack_fire_front_weight_min = float(args.pz_crack_fire_fw_min)
    if getattr(args, 'pz_crack_fire_griffith_frac', None) is not None:
        cfg.process_zone.crack_fire_griffith_frac = float(args.pz_crack_fire_griffith_frac)
    if getattr(args, 'pz_crack_consume_toughening', None) is not None:
        cfg.process_zone.crack_fire_consume_toughening_frac = float(args.pz_crack_consume_toughening)
    if getattr(args, 'pz_crack_consume_radius', None) is not None:
        cfg.process_zone.crack_fire_consume_radius_factor = float(args.pz_crack_consume_radius)
    if getattr(args, 'pz_crack_fired_gc_relief', None) is not None:
        cfg.process_zone.crack_fired_gc_relief = float(args.pz_crack_fired_gc_relief)
    if getattr(args, 'pz_crack_source_fw_min', None) is not None:
        cfg.process_zone.crack_source_fw_min = float(args.pz_crack_source_fw_min)
    if getattr(args, 'plastic_taylor_athermal_alpha', None) is not None:
        cfg.dislocations.taylor_athermal_alpha = float(args.plastic_taylor_athermal_alpha)
    if args.tip_wake_relax is not None:
        cfg.tip_memory.wake_relax = float(args.tip_wake_relax)
    if args.pz_crack_resolved_stress_scale is not None:
        cfg.process_zone.crack_resolved_stress_scale = float(args.pz_crack_resolved_stress_scale)
    if args.pz_crack_eta0 is not None:
        cfg.process_zone.crack_eta0 = float(args.pz_crack_eta0)
    if args.pz_crack_H0_eV is not None:
        cfg.process_zone.crack_H0_eV = float(args.pz_crack_H0_eV)
    if args.pz_crack_drive_exponent is not None:
        cfg.process_zone.crack_drive_exponent = float(args.pz_crack_drive_exponent)
    if args.pz_crack_probability_cap is not None:
        cfg.process_zone.crack_probability_cap = float(args.pz_crack_probability_cap)
    if args.pz_crack_drive_scale is not None:
        cfg.process_zone.crack_drive_scale = float(args.pz_crack_drive_scale)
    if args.no_pz_crack_effective_drive:
        cfg.process_zone.crack_use_effective_drive = False
    if args.pz_crack_effective_drive_mix is not None:
        cfg.process_zone.crack_effective_drive_mix = float(args.pz_crack_effective_drive_mix)
    if args.pz_crack_front_radius_factor is not None:
        cfg.process_zone.crack_front_radius_factor = float(args.pz_crack_front_radius_factor)
    if args.pz_crack_front_grad_threshold is not None:
        cfg.process_zone.crack_front_grad_threshold = float(args.pz_crack_front_grad_threshold)
    if args.pz_crack_front_wake_dmax is not None:
        cfg.process_zone.crack_front_wake_dmax = float(args.pz_crack_front_wake_dmax)
    if args.save_every is not None:
        cfg.diagnostics.save_fields = True
        cfg.diagnostics.save_every = max(args.save_every, 1)
        cfg.diagnostics.save_field_pngs = True
    if getattr(args, 'no_progress', False):
        cfg.diagnostics.progress = False
    if getattr(args, 'progress_interval', None) is not None:
        cfg.diagnostics.progress_interval_s = max(float(args.progress_interval), 0.0)
    if getattr(args, 'progress_every', None) is not None:
        cfg.diagnostics.progress_every = max(int(args.progress_every), 1)
    if args.memory_mode is not None:
        cfg.tip_memory.mode = args.memory_mode
        cfg.tip_memory.enabled = args.memory_mode != 'off'
    if args.tip_memory_gain is not None:
        cfg.tip_memory.state_gain = args.tip_memory_gain
    if args.tip_M_max is not None:
        cfg.tip_memory.M_max = args.tip_M_max
    if args.tip_amp_max is not None:
        cfg.tip_memory.amp_max = args.tip_amp_max
    if args.tip_shield_max is not None:
        cfg.tip_memory.shield_max = args.tip_shield_max
    if args.tip_blunt_work is not None:
        cfg.tip_memory.blunt_per_work = args.tip_blunt_work
    if args.tip_sharpen_damage is not None:
        cfg.tip_memory.sharpen_per_damage = args.tip_sharpen_damage
    if args.tip_drive_exponent is not None:
        cfg.tip_memory.drive_exponent = args.tip_drive_exponent
    if args.no_tip_drive_coupling:
        cfg.tip_memory.couple_to_damage_drive = False
    if args.enable_kinetic_damage_drive:
        cfg.phase_field.use_kinetic_damage_drive = True
    if args.pf_damage_cap is not None:
        cfg.phase_field.max_damage_increment_per_stagger = args.pf_damage_cap
    if args.rho_cap is not None:
        cfg.dislocations.rho_cap = args.rho_cap
    if args.dot_ep_max is not None:
        cfg.dislocations.dot_ep_max = args.dot_ep_max
    if args.plastic_update_mode is not None:
        cfg.dislocations.plastic_update_mode = args.plastic_update_mode
    if args.flow_epsdot_ref is not None:
        cfg.dislocations.flow_epsdot_ref = args.flow_epsdot_ref
    if getattr(args, 'peierls_autocalibrate', False):
        cfg.dislocations.peierls_autocalibrate = True
    if getattr(args, 'peierls_floor_min_MPa', None) is not None:
        cfg.dislocations.peierls_floor_min_MPa = float(args.peierls_floor_min_MPa)
    if getattr(args, 'peierls_cal_T', None) is not None:
        cfg.dislocations.peierls_cal_T_K = float(args.peierls_cal_T)
    if getattr(args, 'peierls_S_kB', None) is not None:
        cfg.dislocations.peierls_S_kB = float(args.peierls_S_kB)
    if getattr(args, 'peierls_v0_b3', None) is not None:
        cfg.dislocations.peierls_v0_b3 = float(args.peierls_v0_b3)
    if getattr(args, 'taylor_multihit', False):
        cfg.dislocations.taylor_multihit = True
    if getattr(args, 'taylor_corr_rho_c', None) is not None:
        cfg.dislocations.taylor_corr_rho_c = float(args.taylor_corr_rho_c)
    if getattr(args, 'taylor_renewal_time', None) is not None:
        cfg.dislocations.taylor_renewal_time_s = float(args.taylor_renewal_time)
    if getattr(args, 'taylor_m_max', None) is not None:
        cfg.dislocations.taylor_m_max = float(args.taylor_m_max)
    if getattr(args, 'taylor_m_exponent', None) is not None:
        cfg.dislocations.taylor_m_exponent = float(args.taylor_m_exponent)
    if args.thermo_consistency_mode is not None:
        cfg.dislocations.thermo_consistency_mode = args.thermo_consistency_mode
    if args.thermo_event_strain is not None:
        cfg.dislocations.thermo_event_strain = args.thermo_event_strain
    if args.thermo_onsager_max_fraction is not None:
        cfg.dislocations.thermo_onsager_max_fraction = args.thermo_onsager_max_fraction
    if args.thermo_work_mode is not None:
        cfg.dislocations.thermo_use_avg_stress_work = (args.thermo_work_mode == 'avg_stress')
    if args.thermo_adaptive_substepping:
        cfg.dislocations.thermo_adaptive_substepping = True
    if args.thermo_max_substeps is not None:
        cfg.dislocations.thermo_max_substeps = int(args.thermo_max_substeps)
    if args.thermo_max_dep_increment is not None:
        cfg.dislocations.thermo_max_dep_increment = float(args.thermo_max_dep_increment)
    if args.thermo_max_hazard_increment is not None:
        cfg.dislocations.thermo_max_hazard_increment = float(args.thermo_max_hazard_increment)
    if args.no_thermo_energy_audit:
        cfg.dislocations.thermo_energy_audit = False
    if args.thermo_energy_rel_tol is not None:
        cfg.dislocations.thermo_energy_rel_tol = float(args.thermo_energy_rel_tol)
    if args.memory_energetics:
        cfg.tip_memory.use_memory_energetics = True
    if args.no_memory_energetics:
        cfg.tip_memory.use_memory_energetics = False
    if args.memory_energy_r_coeff is not None:
        cfg.tip_memory.memory_energy_r_coeff = float(args.memory_energy_r_coeff)
    if args.memory_energy_z_coeff is not None:
        cfg.tip_memory.memory_energy_z_coeff = float(args.memory_energy_z_coeff)
    if args.memory_dissipation_r_coeff is not None:
        cfg.tip_memory.memory_dissipation_r_coeff = float(args.memory_dissipation_r_coeff)
    if args.memory_dissipation_z_coeff is not None:
        cfg.tip_memory.memory_dissipation_z_coeff = float(args.memory_dissipation_z_coeff)
    if args.plastic_H0_eV is not None:
        cfg.plasticity_barrier.H0_J = args.plastic_H0_eV * EV_TO_J
    if args.plastic_vstar_b3 is not None:
        cfg.plasticity_barrier.v0_c = args.plastic_vstar_b3 * cfg.material.b**3
    if args.plastic_barrier_model == 'rational_Hv':
        cfg.plasticity_barrier.model_type = 'rational_Hv'
    elif args.plastic_barrier_model == 'exp_floor':
        if args.plastic_barrier_json is None:
            raise ValueError('--plastic-barrier-json is required for --plastic-barrier-model exp_floor')
        if not _plastic_barrier_system_request_is_all(args.plastic_barrier_system):
            _load_exp_floor_plastic_barrier(
                cfg, args.plastic_barrier_json, args.plastic_barrier_system,
                energy_scale=args.plastic_barrier_scale,
                entropy_scale=args.plastic_exp_entropy_scale,
                stress_scale=args.plastic_exp_stress_scale,
            )
    if args.plastic_exp_vmax_b3 is not None:
        cfg.plasticity_barrier.exp_v_max_b3 = args.plastic_exp_vmax_b3
    if args.max_plastic_strain_increment is not None:
        cfg.dislocations.max_plastic_strain_increment = args.max_plastic_strain_increment
    if args.max_rho_relative_increment is not None:
        cfg.dislocations.max_rho_relative_increment = args.max_rho_relative_increment
    if args.disable_plasticity:
        cfg.dislocations.enable_plasticity = False
    if args.freeze_rho:
        cfg.dislocations.freeze_rho = True
    if args.disable_wp_gc_coupling:
        cfg.phase_field.wp_gc_coupling_mode = 'off'
        cfg.phase_field.plastic_work_to_Gc_efficiency = 0.0
        cfg.phase_field.Gc_local_cap_factor = 1.0
    if args.stop_on_invalid:
        cfg.stop_on_invalid = True
    if args.invalid_wp_wext_pct is not None:
        cfg.invalid_wp_wext_pct = args.invalid_wp_wext_pct
    if args.invalid_min_step is not None:
        cfg.invalid_min_step = args.invalid_min_step
    if args.invalid_wext_min is not None:
        cfg.invalid_Wext_min = args.invalid_wext_min
    if args.invalid_sigma_eq_GPa is not None:
        cfg.invalid_sigma_eq_GPa = float(args.invalid_sigma_eq_GPa)
    if args.invalid_dep_eq_increment is not None:
        cfg.invalid_dep_eq_increment = float(args.invalid_dep_eq_increment)
    if args.allow_invalid_soft_tearing:
        cfg.invalid_wp_wext_pct = np.inf
    if args.no_auto_stop:
        cfg.auto_stop.enabled = False

    # Optional baseline Gc sweep, with optional EXP_floor material-barrier sweep.
    if args.Gc_list is not None:
        base_out = cfg.output_dir
        gc_values = [float(x) for x in args.Gc_list]
        if args.plastic_barrier_model == 'exp_floor' and _plastic_barrier_system_request_is_all(args.plastic_barrier_system):
            systems = _systems_from_request(args.plastic_barrier_json, args.plastic_barrier_system)
        else:
            systems = [None]
        for Gc0 in gc_values:
            for sysname in systems:
                cfg_i = copy.deepcopy(cfg)
                cfg_i.phase_field.Gc0_athermal = float(Gc0)
                parts = [base_out, f"Gc_{_safe_float_tag(Gc0)}Jm2"]
                if sysname is not None:
                    parts.append(_safe_system_tag(sysname))
                    _load_exp_floor_plastic_barrier(
                        cfg_i, args.plastic_barrier_json, sysname,
                        energy_scale=args.plastic_barrier_scale,
                        entropy_scale=args.plastic_exp_entropy_scale,
                        stress_scale=args.plastic_exp_stress_scale,
                    )
                    if args.plastic_exp_vmax_b3 is not None:
                        cfg_i.plasticity_barrier.exp_v_max_b3 = args.plastic_exp_vmax_b3
                cfg_i.output_dir = os.path.join(*parts)
                run_simulation(cfg_i)
        return

    # Optional small material sweep at one temperature.  This is not a general
    # H0/T sweep; it is a diagnostic way to ask whether the EXP_floor parameter
    # set, especially sigc(T) and entropy, is suppressing/activating plasticity
    # at the selected temperature.  Si is excluded by default because it is not
    # a good metallic dislocation-plasticity analog here.
    if (args.plastic_barrier_model == 'exp_floor' and
            _plastic_barrier_system_request_is_all(args.plastic_barrier_system)):
        systems = _systems_from_request(args.plastic_barrier_json, args.plastic_barrier_system)
        base_out = cfg.output_dir
        for sysname in systems:
            cfg_i = copy.deepcopy(cfg)
            cfg_i.output_dir = os.path.join(base_out, _safe_system_tag(sysname))
            _load_exp_floor_plastic_barrier(
                cfg_i, args.plastic_barrier_json, sysname,
                energy_scale=args.plastic_barrier_scale,
                entropy_scale=args.plastic_exp_entropy_scale,
                stress_scale=args.plastic_exp_stress_scale,
            )
            if args.plastic_exp_vmax_b3 is not None:
                cfg_i.plasticity_barrier.exp_v_max_b3 = args.plastic_exp_vmax_b3
            run_simulation(cfg_i)
        return

    run_simulation(cfg)


if __name__ == '__main__':
    main()
