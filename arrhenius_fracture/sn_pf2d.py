"""Full-field 2-D S-N crack-initiation simulation from a blunt surface feature.

This module is deliberately separate from the pre-cracked DeltaK/Paris driver.
The initial condition has d(x)=0 everywhere.  Stress concentration is generated
by a smooth half-elliptical free-surface notch.  Cyclic Arrhenius hazards evolve
spatial plastic-state and localization fields; those fields modify both the
local crack-opening hazard and the phase-field fracture resistance.

The two pilot cases have identical barriers and geometry:
  no_shield : plastic localization/stored-energy term only
  shielded  : same localization plus back-stress/barrier shielding feedback

The code records both a local nucleation-clock first passage and the emergence of
a connected phase-field crack from the blunt notch root.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import gammainc
from scipy.sparse.linalg import spsolve

from .config import KB, EV_TO_J, make_emergent_config
from .fem import (
    plane_strain_D, assemble_mechanics, solve_dirichlet,
    assemble_pf_matrices, project_gp_to_nodes, stress_state,
)
from .materials import FractureModel
from .phase_field import update_phase_field, at2_surface_energy
from .sn_geometry import BluntNotchGeometry, make_blunt_edge_notch_mesh
from .sn_v1 import make_barriers, KBEV


def _write_csv(path: Path, rows):
    if not rows:
        return
    keys = sorted(set().union(*(r.keys() for r in rows)))
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(rows)


def _cycle_hazard_field(barrier, sigma_max_node, T_K, R, frequency_Hz, n_phase,
                        shift_eV=None, tensile_only=False, multihit_m=1.0, tau_c_s=1e-6):
    phase = np.linspace(0.0, 2.0 * np.pi, n_phase, endpoint=False)
    f = 0.5 * (1.0 + R) + 0.5 * (1.0 - R) * np.cos(phase)
    rates = []
    for fac in f:
        sig = fac * sigma_max_node
        if tensile_only:
            sig = np.maximum(sig, 0.0)
        else:
            sig = np.abs(sig)
        G = barrier.deltaG_eV(sig, T_K)
        if shift_eV is not None:
            G = np.maximum(G + shift_eV, 1e-12)
        lam = barrier.rate_prefactor * np.exp(np.clip(-G / max(KBEV * T_K, 1e-30), -700.0, 0.0))
        if multihit_m > 1.0 + 1e-12:
            lam = gammainc(multihit_m, np.minimum(lam * max(tau_c_s, 1e-30), 1e12)) / max(tau_c_s, 1e-30)
        rates.append(lam)
    return np.mean(np.asarray(rates), axis=0) / max(frequency_Hz, 1e-300)




def _solve_symmetric_tension(K, Rint, u, bnd, mesh, Uy_top, Uy_bot):
    """Symmetric displacement control with a single x anchor.

    The general fracture solver fixes x at both bottom corners, which is useful
    for the original notched plate but creates grip-corner stress concentrations
    in an initiation study.  Here top/bottom y displacement is prescribed and a
    single right-mid node fixes rigid x translation, allowing Poisson contraction.
    """
    ndof = mesh.ndof
    prescribed = np.zeros(ndof, dtype=bool)
    up = np.zeros(ndof)
    prescribed[2*bnd.top_nodes + 1] = True
    up[2*bnd.top_nodes + 1] = Uy_top
    prescribed[2*bnd.bot_nodes + 1] = True
    up[2*bnd.bot_nodes + 1] = Uy_bot
    x, y = mesh.nodes[:,0], mesh.nodes[:,1]
    anchor = int(np.argmin((x-mesh.nodes[:,0].max())**2 + y**2))
    prescribed[2*anchor] = True
    up[2*anchor] = 0.0
    free = ~prescribed
    Kc = K.tocsr()
    du_p = up[prescribed] - u[prescribed]
    rhs = -Rint[free] - Kc[np.ix_(free, prescribed)] @ du_p
    un = u.copy()
    un[free] = u[free] + spsolve(Kc[np.ix_(free, free)], rhs)
    un[prescribed] = up[prescribed]
    Rfull = Rint + Kc @ (un-u)
    Ftop = float(np.sum(Rfull[2*bnd.top_nodes + 1]))
    return un, Ftop


def _calibrate_Umax(mesh, bnd, cfg, Dmat, sigma_max_target_Pa):
    u = np.zeros(mesh.ndof)
    ep = np.zeros((3, mesh.ne))
    rho = np.full(mesh.ne, 1e12)
    d = np.zeros(mesh.nn)
    K, Rint, *_ = assemble_mechanics(mesh, u, ep, rho, d, Dmat, cfg.material)
    Uprobe = 1e-8
    u1, Ftop = _solve_symmetric_tension(K, Rint, u, bnd, mesh, Uprobe, -Uprobe)
    sigma_nom_probe = abs(Ftop) / max(cfg.geometry.Lx, 1e-30)
    if sigma_nom_probe <= 0 or not np.isfinite(sigma_nom_probe):
        raise RuntimeError("Failed nominal-stress displacement calibration")
    return Uprobe * sigma_max_target_Pa / sigma_nom_probe


def _adjacency(mesh):
    adj = [set() for _ in range(mesh.nn)]
    for tri in mesh.elems:
        a, b, c = map(int, tri)
        adj[a].update((b, c)); adj[b].update((a, c)); adj[c].update((a, b))
    return adj


def _connected_crack_extent(mesh, d, root_xy, dthr, root_seed_radius, adj):
    root = np.asarray(root_xy)
    active = np.asarray(d >= dthr, bool)
    dist = np.linalg.norm(mesh.nodes - root[None, :], axis=1)
    seeds = np.where(active & (dist <= root_seed_radius))[0]
    if len(seeds) == 0:
        return 0.0, 0
    seen = set(map(int, seeds)); stack = list(map(int, seeds))
    while stack:
        i = stack.pop()
        for j in adj[i]:
            if active[j] and j not in seen:
                seen.add(j); stack.append(j)
    idx = np.fromiter(seen, dtype=int)
    extent = max(float(np.max(mesh.nodes[idx, 0]) - root[0]), 0.0)
    return extent, len(idx)


def _plot_fields(mesh, fields, out_png: Path, root_xy):
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
        ax.set_aspect("equal"); ax.set_title(title); ax.set_xlabel("x (mm)"); ax.set_ylabel("y (mm)")
        fig.colorbar(pc, ax=ax)
    fig.savefig(out_png, dpi=220)
    plt.close(fig)


def run_case_stress(args, case_name: str, sigma_a_MPa: float):
    shield_on = case_name == "shielded"
    cfg = make_emergent_config()
    cfg.geometry.Lx = args.Lx
    cfg.geometry.Ly = args.Ly
    geom = BluntNotchGeometry(args.Lx, args.Ly, args.notch_depth_m, args.notch_half_height_m)
    mesh, bnd, root_xy = make_blunt_edge_notch_mesh(
        geom, nx=args.nx, ny=args.ny, jitter=args.jitter,
        root_h_fine=args.root_h_fine, seed=args.seed,
    )
    Dmat = plane_strain_D(cfg.material)
    Md, Kd = assemble_pf_matrices(mesh)
    ell = args.ell_m if args.ell_m is not None else args.ell_factor * mesh.hbar_tip
    fracture_model = FractureModel(cfg.fracture_barrier, cfg.material, pf=cfg.phase_field)
    Gc_base = float(fracture_model.Gc_of_T(args.T, ell, method=args.Gc_method))

    sigma_max = 2.0 * sigma_a_MPa * 1e6 / max(1.0 - args.R, 1e-30)
    Umax = _calibrate_Umax(mesh, bnd, cfg, Dmat, sigma_max)

    emit, crack = make_barriers(args.S_emit_kB, args.S_crack_kB, args.emit_energy_scale)
    u = np.zeros(mesh.ndof)
    ep_gp = np.zeros((3, mesh.ne))
    rho_gp = np.full(mesh.ne, args.rho0)
    d = np.zeros(mesh.nn)   # critical: no initial crack
    Hhist = np.zeros(mesh.nn)
    B_emit = np.zeros(mesh.nn)
    B_nuc = np.zeros(mesh.nn)
    P = np.zeros(mesh.nn)
    Dloc = np.zeros(mesh.nn)
    Gc_eff = np.full(mesh.nn, Gc_base)
    adj = _adjacency(mesh)

    outdir = Path(args.out) / case_name / (f"sigmaA_{sigma_a_MPa:g}MPa".replace(".", "p"))
    outdir.mkdir(parents=True, exist_ok=True)
    with (outdir / "run_args.json").open("w") as f:
        json.dump(vars(args) | {"case": case_name, "sigma_a_MPa": sigma_a_MPa,
                                "root_radius_m": geom.root_radius, "ell_m": ell,
                                "Umax_m": Umax}, f, indent=2, sort_keys=True)

    cycles = 0.0
    cycles_clock = None
    cycles_pf = None
    rows = []
    last_sigma = np.zeros(mesh.nn)

    for ib in range(args.max_blocks):
        if cycles >= args.cycles_max:
            break
        # Mechanical equilibrium at cycle maximum.
        K, Rint, *_ = assemble_mechanics(mesh, u, ep_gp, rho_gp, d, Dmat, cfg.material)
        u, Ftop = _solve_symmetric_tension(K, Rint, u, bnd, mesh, Umax, -Umax)
        _, sigma_eq, sigma1, psi_gp = stress_state(mesh, u, ep_gp, d, Dmat, cfg.material)
        sigma_emit_node = np.maximum(project_gp_to_nodes(mesh, sigma_eq), 0.0)
        sigma_nuc_node = np.maximum(project_gp_to_nodes(mesh, sigma1), 0.0)
        psi_node = project_gp_to_nodes(mesh, psi_gp)
        last_sigma = sigma_nuc_node.copy()

        # Internal-state coupling.  The two cases have identical plastic-event
        # physics.  They differ only in whether the accumulated plastic state
        # feeds back onto crack opening (chi/Gsh) and PF resistance.
        chi = args.shield_chi if shield_on else 0.0
        Gsh = args.Gshield_eV if shield_on else 0.0

        def hazards_for_state(P_state, D_state):
            sigma_back_state = args.sigma_back_max_GPa * 1e9 * P_state
            sig_emit_eff_state = np.maximum(sigma_emit_node - sigma_back_state, 0.0)
            mu_emit_state = args.emit_site_multiplicity * _cycle_hazard_field(
                emit, sig_emit_eff_state, args.T, args.R, args.frequency_Hz, args.n_phase,
                tensile_only=False,
            )
            sig_nuc_eff_state = np.maximum(
                sigma_nuc_node - chi * sigma_back_state, 0.0
            )
            state_shift_state = Gsh * P_state - args.Gstored_eV * D_state
            mu_nuc_state = _cycle_hazard_field(
                crack, sig_nuc_eff_state, args.T, args.R, args.frequency_Hz, args.n_phase,
                shift_eV=state_shift_state, tensile_only=True,
                multihit_m=args.multihit_m, tau_c_s=args.multihit_tau_s,
            )
            return mu_emit_state, mu_nuc_state, sigma_back_state, state_shift_state

        mu_emit_pre, mu_nuc_pre, sigma_back_pre, state_shift_pre = hazards_for_state(P, Dloc)

        remaining = args.cycles_max - cycles
        dN = min(args.block_cycles, remaining)
        max_mu_e = float(np.max(mu_emit_pre))
        max_mu_n_pre = float(np.max(mu_nuc_pre))
        if max_mu_e > 0:
            sensP = (1.0 - P) / max(args.B_shield_events, 1e-30) / max(args.target_dP, 1e-30)
            sensD = (1.0 - Dloc) / max(args.B_damage_events, 1e-30) / max(args.target_dD, 1e-30)
            state_rate = mu_emit_pre * np.maximum(sensP, sensD)
            max_state_rate = float(np.max(state_rate))
            if max_state_rate > 0:
                dN = min(dN, 1.0 / max_state_rate)
        if max_mu_n_pre > 0:
            dN = min(dN, args.target_dB_nuc / max_mu_n_pre)
        dN = max(min(dN, remaining), args.min_block_cycles)
        dN = min(dN, remaining)
        if dN <= 0:
            break

        # Predictor-corrector integration over the accepted cycle block.  The
        # original pilot used only the PRE-state crack hazard, so a block that
        # generated shielding could not benefit from that shielding until the
        # next block.  Here both plastic and crack clocks use a trapezoidal
        # state integration.  If crack hazard grows strongly inside the block,
        # the block is shortened and retried.
        B_emit_old = B_emit.copy()
        B_nuc_old = B_nuc.copy()
        mu_emit_post = mu_emit_pre
        mu_nuc_post = mu_nuc_pre
        for _pc in range(4):
            # Euler predictor for the plastic clock/state.
            B_emit_pred = B_emit_old + mu_emit_pre * dN
            P_pred = 1.0 - np.exp(-B_emit_pred / max(args.B_shield_events, 1e-30))
            D_pred = 1.0 - np.exp(-B_emit_pred / max(args.B_damage_events, 1e-30))
            mu_emit_pred, _, _, _ = hazards_for_state(P_pred, D_pred)

            # Correct plastic clock, then evaluate crack hazard at corrected state.
            dB_emit = 0.5 * (mu_emit_pre + mu_emit_pred) * dN
            B_emit_trial = B_emit_old + dB_emit
            P_trial = 1.0 - np.exp(-B_emit_trial / max(args.B_shield_events, 1e-30))
            D_trial = 1.0 - np.exp(-B_emit_trial / max(args.B_damage_events, 1e-30))
            mu_emit_post, mu_nuc_post, sigma_back_post, state_shift_post = hazards_for_state(P_trial, D_trial)
            dB_nuc = 0.5 * (mu_nuc_pre + mu_nuc_post) * dN

            max_dB_nuc = float(np.max(dB_nuc))
            if max_dB_nuc <= 1.001 * args.target_dB_nuc or dN <= args.min_block_cycles:
                break
            dN_new = 0.98 * dN * args.target_dB_nuc / max(max_dB_nuc, 1e-300)
            dN = max(min(dN_new, dN, remaining), args.min_block_cycles)

        B_emit = B_emit_trial
        P = P_trial
        Dloc = D_trial
        B_nuc = B_nuc_old + dB_nuc
        cycles_old = cycles
        cycles += dN

        max_mu_n = float(np.max(0.5 * (mu_nuc_pre + mu_nuc_post)))

        # Smooth phase-field resistance evolution.  No hard Gc floor/cap is used;
        # both effects are bounded because P,Dloc in [0,1].
        log_ratio = args.shield_log_gain * P * (1.0 if shield_on else 0.0) - args.damage_log_drop * Dloc
        Gc_eff = Gc_base * np.exp(log_ratio)

        # Hazard-weight the variational PF update.  p_block is the probability of
        # at least one crack-opening event during this accepted cycle block.
        p_block = 1.0 - np.exp(-np.clip(dB_nuc, 0.0, 700.0))
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

        if cycles_clock is None and float(np.max(B_nuc)) >= 1.0:
            crossing = (B_nuc_old < 1.0) & (B_nuc >= 1.0) & (dB_nuc > 0.0)
            if np.any(crossing):
                frac = (1.0 - B_nuc_old[crossing]) / np.maximum(dB_nuc[crossing], 1e-300)
                cycles_clock = cycles_old + float(np.clip(np.min(frac), 0.0, 1.0)) * dN
            else:
                cycles_clock = cycles
        extent, nconn = _connected_crack_extent(
            mesh, d, root_xy, args.pf_damage_threshold,
            args.root_seed_radius_factor * ell, adj,
        )
        if cycles_pf is None and extent >= args.pf_crack_extent_factor * ell:
            cycles_pf = cycles

        rows.append({
            "block": ib, "case": case_name, "sigma_a_MPa": sigma_a_MPa,
            "cycles_total": cycles, "dN": dN, "Ftop_N_per_m": float(Ftop),
            "sigma_max_node_Pa": float(np.max(sigma_nuc_node)),
            "mu_emit_max": max_mu_e, "mu_nuc_max": max_mu_n,
            "mu_nuc_pre_max": float(np.max(mu_nuc_pre)),
            "mu_nuc_post_max": float(np.max(mu_nuc_post)),
            "dB_nuc_max": float(np.max(dB_nuc)),
            "backstress_max_Pa": float(np.max(sigma_back_post)),
            "state_shift_min_eV": float(np.min(state_shift_post)),
            "state_shift_max_eV": float(np.max(state_shift_post)),
            "B_emit_max": float(np.max(B_emit)), "B_nuc_max": float(np.max(B_nuc)),
            "P_max": float(np.max(P)), "Dloc_max": float(np.max(Dloc)),
            "Gc_min_J_m2": float(np.min(Gc_eff)), "Gc_max_J_m2": float(np.max(Gc_eff)),
            "d_max": float(np.max(d)), "connected_crack_extent_m": extent,
            "connected_crack_nodes": nconn,
            "cycles_to_nucleation_clock": cycles_clock if cycles_clock is not None else np.nan,
            "cycles_to_pf_crack": cycles_pf if cycles_pf is not None else np.nan,
            "at2_surface_energy_J_per_m": float(at2_surface_energy(mesh, d, ell, Gc_eff)),
        })

        if args.print_every and ib % args.print_every == 0:
            print(f"2D {case_name} sigma_a={sigma_a_MPa:g} MPa block={ib} N={cycles:.3e} "
                  f"Bnu={np.max(B_nuc):.3g} P={np.max(P):.3g} D={np.max(Dloc):.3g} "
                  f"dmax={np.max(d):.3g} crack={extent*1e6:.2f}um")
        if ib == 0 or (args.snapshot_every > 0 and ib % args.snapshot_every == 0):
            _plot_fields(mesh, [("phase field d", d), ("plastic state P", P),
                                ("localization D", Dloc), ("Gc eff", Gc_eff)],
                         outdir / f"fields_block_{ib:05d}.png", root_xy)
        if args.stop_after_pf_crack and cycles_pf is not None:
            break

    _write_csv(outdir / "sn_pf2d_history.csv", rows)
    _plot_fields(mesh, [("phase field d", d), ("plastic state P", P),
                        ("localization D", Dloc), ("sigma1 max-load (GPa)", last_sigma*1e-9)],
                 outdir / "fields_final.png", root_xy)
    summary = {
        "model": "SN_PF2D_blunt_surface_feature",
        "case": case_name,
        "sigma_a_MPa": sigma_a_MPa,
        "T_K": args.T,
        "R": args.R,
        "frequency_Hz": args.frequency_Hz,
        "cycles_total": cycles,
        "cycles_to_nucleation_clock": cycles_clock,
        "cycles_to_pf_crack": cycles_pf,
        "status": "pf_crack" if cycles_pf is not None else ("clock_first_passage" if cycles_clock is not None else "right_censored"),
        "notch_depth_m": geom.depth_a,
        "notch_half_height_m": geom.half_height_b,
        "notch_root_radius_m": geom.root_radius,
        "ell_m": ell,
        "Gc_base_J_m2": Gc_base,
        "chi_back": chi,
        "Gshield_eV": Gsh,
        "B_emit_final_max": float(np.max(B_emit)),
        "B_nuc_final_max": float(np.max(B_nuc)),
        "P_final_max": float(np.max(P)),
        "Dloc_final_max": float(np.max(Dloc)),
        "d_final_max": float(np.max(d)),
        "Gc_final_min_J_m2": float(np.min(Gc_eff)),
        "history_csv": str(outdir / "sn_pf2d_history.csv"),
    }
    with (outdir / "summary.json").open("w") as f:
        json.dump(summary, f, indent=2, sort_keys=True)
    return summary


def run_sweep(args):
    Path(args.out).mkdir(parents=True, exist_ok=True)
    summaries = []
    for case in ("no_shield", "shielded"):
        for s in args.sigma_a_MPa:
            print(f"=== 2D case={case} sigma_a={s:g} MPa ===")
            summaries.append(run_case_stress(args, case, float(s)))
    _write_csv(Path(args.out) / "sn_pf2d_summary.csv", summaries)
    return summaries


def build_parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", default="runs/sn_pf2d_two_case")
    p.add_argument("--T", type=float, default=300.0)
    p.add_argument("--sigma-a-MPa", nargs="+", type=float, default=[500, 600, 700], dest="sigma_a_MPa")
    p.add_argument("--R", type=float, default=0.1)
    p.add_argument("--frequency-Hz", type=float, default=1000.0, dest="frequency_Hz")
    p.add_argument("--cycles-max", type=float, default=1e9, dest="cycles_max")
    p.add_argument("--block-cycles", type=float, default=1e6, dest="block_cycles")
    p.add_argument("--min-block-cycles", type=float, default=1e-6, dest="min_block_cycles")
    p.add_argument("--max-blocks", type=int, default=3000, dest="max_blocks")
    p.add_argument("--target-dP", type=float, default=0.02, dest="target_dP")
    p.add_argument("--target-dD", type=float, default=0.02, dest="target_dD")
    p.add_argument("--target-dB-nuc", type=float, default=0.05, dest="target_dB_nuc")
    p.add_argument("--n-phase", type=int, default=64, dest="n_phase")

    p.add_argument("--Lx", type=float, default=2e-3)
    p.add_argument("--Ly", type=float, default=4e-3)
    p.add_argument("--notch-depth-m", type=float, default=0.15e-3, dest="notch_depth_m")
    p.add_argument("--notch-half-height-m", type=float, default=0.30e-3, dest="notch_half_height_m")
    p.add_argument("--nx", type=int, default=50)
    p.add_argument("--ny", type=int, default=100)
    p.add_argument("--jitter", type=float, default=0.08)
    p.add_argument("--root-h-fine", type=float, default=20e-6, dest="root_h_fine")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--ell-m", type=float, default=None, dest="ell_m")
    p.add_argument("--ell-factor", type=float, default=3.0, dest="ell_factor")
    p.add_argument("--Gc-method", choices=["lambertw", "hazard"], default="lambertw", dest="Gc_method")
    p.add_argument("--rho0", type=float, default=1e12)

    p.add_argument("--S-emit-kB", type=float, default=-40.0, dest="S_emit_kB")
    p.add_argument("--emit-energy-scale", type=float, default=0.75, dest="emit_energy_scale")
    p.add_argument("--S-crack-kB", type=float, default=0.0, dest="S_crack_kB")
    p.add_argument("--emit-site-multiplicity", type=float, default=5e8, dest="emit_site_multiplicity")
    p.add_argument("--B-shield-events", type=float, default=50.0, dest="B_shield_events")
    p.add_argument("--B-damage-events", type=float, default=500.0, dest="B_damage_events")
    p.add_argument("--sigma-back-max-GPa", type=float, default=1.0, dest="sigma_back_max_GPa")
    p.add_argument("--shield-chi", type=float, default=0.6, dest="shield_chi")
    p.add_argument("--Gshield-eV", type=float, default=0.35, dest="Gshield_eV")
    p.add_argument("--Gstored-eV", type=float, default=0.25, dest="Gstored_eV")
    p.add_argument("--shield-log-gain", type=float, default=0.8, dest="shield_log_gain")
    p.add_argument("--damage-log-drop", type=float, default=1.5, dest="damage_log_drop")
    p.add_argument("--multihit-m", type=float, default=3.0, dest="multihit_m")
    p.add_argument("--multihit-tau-s", type=float, default=1e-6, dest="multihit_tau_s")

    p.add_argument("--n-stagger", type=int, default=2, dest="n_stagger")
    p.add_argument("--max-damage-increment", type=float, default=0.05, dest="max_damage_increment")
    p.add_argument("--damage-drive-cap", type=float, default=20.0, dest="damage_drive_cap")
    p.add_argument("--pf-damage-threshold", type=float, default=0.5, dest="pf_damage_threshold")
    p.add_argument("--root-seed-radius-factor", type=float, default=3.0, dest="root_seed_radius_factor")
    p.add_argument("--pf-crack-extent-factor", type=float, default=3.0, dest="pf_crack_extent_factor")
    p.add_argument("--stop-after-pf-crack", action="store_true", default=True, dest="stop_after_pf_crack")
    p.add_argument("--snapshot-every", type=int, default=100, dest="snapshot_every")
    p.add_argument("--print-every", type=int, default=10, dest="print_every")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    run_sweep(args)


if __name__ == "__main__":
    main()
