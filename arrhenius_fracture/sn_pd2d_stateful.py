"""Blunt-scratch S-N initiation with intact FEM and a stateful PD patch.

Active architecture
-------------------
1. Intact 2-D FEM resolves cyclic stress, Arrhenius plastic eigenstrain,
   dislocation-density evolution, residual stress, and optional ALE scratch
   reshaping.
2. A local nonlocal spring/peridynamic patch receives the FEM shell motion and
   plastic eigenstrain and re-equilibrates as cohesive bonds soften.
3. Candidate sites undergo finite-memory multi-hit completion, reversible
   embryo formation, stabilization/healing, stable-defect growth, and spatial
   bond linkage.
4. Crack formation is a graph-connectivity event in the softened/broken bond
   network.  No phase-field or AT1/AT2 module is imported or called.

This first pilot is one-way coupled: PD bond degradation redistributes the local
patch deformation but does not yet modify the global FEM stiffness.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

from .config import ElasticProperties
from .sn_arrhenius_chain import build_chain_from_namespace
from .sn_geometry import (
    BluntNotchGeometry,
    apply_local_ale_surface_update,
    identify_feature_surface_nodes,
    local_root_radius,
    local_root_xy,
    make_blunt_edge_notch_mesh,
    rebuild_mesh_geometry,
)
from .sn_intact_fem import (
    affine_stress_control_displacements,
    cycle_stress_histories,
    plane_strain_D,
    project_gp_to_nodes,
    project_plastic_state,
    representative_plastic_cycle,
    stress_state_intact,
    surface_morphology_proposal,
)
from .sn_v1 import make_barriers
from .stateful_peridynamics import StatefulPDConfig, StatefulPDPatch


MODEL_ID = "SN_2D_intact_FEM_stateful_local_peridynamics_v3_1_root_localized"


def _write_csv(path: Path, rows):
    if not rows:
        return
    keys = sorted(set().union(*(r.keys() for r in rows)))
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def _plot_fem_fields(mesh, fields, out_png, root_xy, feature_nodes=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri

    n = len(fields)
    fig, axes = plt.subplots(1, n, figsize=(4.3 * n, 4.2), constrained_layout=True)
    axes = np.atleast_1d(axes)
    tri = mtri.Triangulation(mesh.nodes[:, 0] * 1e3, mesh.nodes[:, 1] * 1e3, mesh.elems)
    for ax, (title, data) in zip(axes, fields):
        pc = ax.tripcolor(tri, np.asarray(data), shading="gouraud")
        if feature_nodes is not None:
            q = mesh.nodes[np.asarray(feature_nodes, int)] * 1e3
            ax.plot(q[:, 0], q[:, 1], "k-", lw=0.8)
        ax.plot(root_xy[0] * 1e3, root_xy[1] * 1e3, "rx", ms=7)
        ax.set_aspect("equal")
        ax.set_title(title)
        ax.set_xlabel("x (mm)")
        ax.set_ylabel("y (mm)")
        fig.colorbar(pc, ax=ax, fraction=0.046, pad=0.04)
    fig.savefig(out_png, dpi=220)
    plt.close(fig)


def _pd_config_from_args(args):
    return StatefulPDConfig(
        patch_radius_m=args.pd_patch_radius_m,
        horizon_m=args.pd_horizon_m,
        boundary_shell_m=args.pd_boundary_shell_m,
        residual_bond_stiffness=args.pd_residual_stiffness,
        initiation_radius_m=args.pd_initiation_radius_m,
        initiation_taper_m=args.pd_initiation_taper_m,
        initiation_back_extent_m=args.pd_initiation_back_extent_m,
        site_density_m2=args.site_density_m2,
        hit_count=args.hit_count,
        hit_memory_s=args.hit_memory_s,
        birth_scale=args.birth_scale,
        nu_stabilize_s=args.nu_stabilize_s,
        nu_heal_s=args.nu_heal_s,
        stabilize_stress_Pa=args.stabilize_stress_GPa * 1e9,
        stabilize_width_Pa=args.stabilize_width_GPa * 1e9,
        stabilize_plastic_gain=args.stabilize_plastic_gain,
        heal_return_fraction=args.heal_return_fraction,
        nu_grow_s=args.nu_grow_s,
        grow_stress_Pa=args.grow_stress_GPa * 1e9,
        grow_width_Pa=args.grow_width_GPa * 1e9,
        stable_count_scale=args.stable_count_scale,
        nu_link_s=args.nu_link_s,
        link_stress_Pa=args.link_stress_GPa * 1e9,
        link_width_Pa=args.link_width_GPa * 1e9,
        link_orientation_power=args.link_orientation_power,
        link_orientation_floor=args.link_orientation_floor,
        neighbor_link_gain=args.neighbor_link_gain,
        max_transition_probability=args.max_transition_probability,
        broken_damage=args.broken_damage,
        root_seed_radius_m=args.root_seed_radius_m,
        established_extent_m=args.established_extent_m,
        pd_amplification_cap=args.pd_amplification_cap,
        amplification_damage_scale=args.pd_amplification_damage_scale,
        random_seed=args.pd_seed if args.pd_seed is not None else args.seed,
        softening_damage=args.softening_damage,
    )


def run_case_stress(args, case_name: str, sigma_a_MPa: float):
    shield_on = case_name == "shielded"
    mat = ElasticProperties(E=args.E_GPa * 1e9, nu=args.nu, b=args.b_m, Tm=args.Tm_K)
    Dmat = plane_strain_D(mat)
    geom = BluntNotchGeometry(args.Lx, args.Ly, args.notch_depth_m, args.notch_half_height_m)
    mesh, bnd, _ = make_blunt_edge_notch_mesh(
        geom,
        nx=args.nx,
        ny=args.ny,
        jitter=args.jitter,
        root_h_fine=args.root_h_fine,
        seed=args.seed,
    )
    feature_nodes = identify_feature_surface_nodes(mesh, geom)
    root_xy = local_root_xy(mesh, feature_nodes)
    root_radius0 = local_root_radius(mesh, feature_nodes)
    fixed_mesh_nodes = np.unique(np.r_[bnd.top_nodes, bnd.bot_nodes])

    plast_chain = build_chain_from_namespace(args, mat.b)
    _, crack = make_barriers(-40.0, args.S_crack_kB, args.emit_energy_scale)
    pd_cfg = _pd_config_from_args(args)
    patch = StatefulPDPatch(mesh, geom, root_xy, mat, pd_cfg)
    pd_state = patch.initial_state()

    sigma_max = 2.0 * sigma_a_MPa * 1e6 / max(1.0 - args.R, 1e-30)
    sigma_min = args.R * sigma_max
    u = np.zeros(mesh.ndof)
    ep_gp = np.zeros((3, mesh.ne))
    rho_gp = np.full(mesh.ne, args.rho0)
    epsp_acc_gp = np.zeros(mesh.ne)

    outdir = Path(args.out) / case_name / (f"sigmaA_{sigma_a_MPa:g}MPa".replace(".", "p"))
    outdir.mkdir(parents=True, exist_ok=True)
    with (outdir / "run_args.json").open("w") as f:
        json.dump(
            vars(args)
            | {
                "case": case_name,
                "sigma_a_MPa": sigma_a_MPa,
                "root_radius_initial_m": root_radius0,
                "pd_points": len(patch.xy),
                "pd_bonds": len(patch.bonds),
            },
            f,
            indent=2,
            sort_keys=True,
        )

    cycles = 0.0
    Wp_total = 0.0
    rows = []
    last_residual = np.zeros(mesh.nn)
    last_hist = None
    last_diag = None

    for ib in range(args.max_blocks):
        if cycles >= args.cycles_max or pd_state.cycles_connected is not None:
            break

        Umax, Umin, u_zero, F0, _ = affine_stress_control_displacements(
            mesh, bnd, mat, Dmat, ep_gp, sigma_max, sigma_min, u
        )
        cyc = representative_plastic_cycle(
            mesh,
            bnd,
            mat,
            Dmat,
            ep_gp,
            rho_gp,
            Umax,
            Umin,
            args.T,
            args.frequency_Hz,
            args.plastic_n_phase,
            plast_chain,
            u_zero,
            args.k_store,
            args.k_dyn,
            args.rho_floor,
            args.rho_cap,
            args.max_dep_phase,
            args.max_rho_rel_phase,
        )
        dep_tensor_cycle = cyc["dep_tensor_cycle"]
        dep_eq_cycle = np.maximum(cyc["dep_eq_cycle"], 0.0)
        drho_cycle = cyc["drho_cycle"]

        epsp_node, rho_node, P, Dloc = project_plastic_state(
            mesh, epsp_acc_gp, rho_gp, args.epsp_shield_scale, args.epsp_damage_scale
        )
        chi = args.shield_chi if shield_on else 0.0
        Gsh = args.Gshield_eV if shield_on else 0.0
        sigma_back = args.sigma_back_max_GPa * 1e9 * P
        state_shift = Gsh * P - args.Gstored_eV * Dloc
        hist_pre = cycle_stress_histories(
            mesh, bnd, mat, Dmat, ep_gp, Umax, Umin, args.hazard_n_phase, u_zero
        )
        ep_node_pre = project_gp_to_nodes(mesh, ep_gp)
        _, _, _, bond_amp_pre, point_amp_pre = patch.solve_local_mechanics(
            pd_state, hist_pre["u_max"], ep_node_pre
        )
        pre_rates = patch.preview_rates(
            pd_state,
            crack,
            hist_pre["sigma_node"],
            args.T,
            args.frequency_Hz,
            state_shift,
            sigma_back,
            chi,
            P,
            point_amp_pre,
            bond_amp_pre,
        )

        remaining = args.cycles_max - cycles
        dN = min(args.block_cycles, remaining)
        max_dep_cycle = float(np.max(dep_eq_cycle)) if dep_eq_cycle.size else 0.0
        if max_dep_cycle > 0.0:
            dN = min(dN, args.target_dep_eq_block / max_dep_cycle)
        rel_rho = float(np.max(np.abs(drho_cycle) / np.maximum(rho_gp, args.rho0)))
        if rel_rho > 0.0:
            dN = min(dN, args.target_rho_rel_block / rel_rho)
        max_hit = float(np.max(pre_rates["mu_hit"]))
        if np.isfinite(args.target_hit_memory) and args.target_hit_memory > 0.0 and max_hit > 0.0:
            dN = min(dN, args.target_hit_memory / max_hit)
        if pre_rates["max_rate_per_cycle"] > 0.0:
            target_hazard = -math.log(max(1.0 - args.max_transition_probability, 1e-12))
            dN = min(dN, target_hazard / pre_rates["max_rate_per_cycle"])

        if args.enable_geometry_evolution and max_dep_cycle > 0.0:
            dh_cycle, _, _ = surface_morphology_proposal(
                mesh,
                feature_nodes,
                dep_tensor_cycle,
                args.morph_band_length_m,
                args.morph_normal_weight,
                args.morph_shear_weight,
            )
            max_dh_cycle = float(np.max(np.abs(dh_cycle))) if dh_cycle.size else 0.0
            if max_dh_cycle > 0.0:
                dN = min(dN, args.target_surface_move_fraction * mesh.hbar_tip / max_dh_cycle)

        dN = max(min(dN, remaining), args.min_block_cycles)
        dN = min(dN, remaining)
        if dN <= 0.0:
            break

        dep_tensor_block = dep_tensor_cycle * dN
        dep_eq_block = dep_eq_cycle * dN
        ep_gp = ep_gp + dep_tensor_block
        rho_gp = np.clip(rho_gp + drho_cycle * dN, args.rho_floor, args.rho_cap)
        epsp_acc_gp = epsp_acc_gp + dep_eq_block
        Wp_total += float(np.sum(cyc["Wp_cycle_gp"] * mesh.area_e) * dN)

        mesh_scale = 0.0
        dh_block = np.zeros(len(feature_nodes))
        gamma_nt = np.zeros(len(feature_nodes))
        if args.enable_geometry_evolution:
            dh_block, _, gamma_nt = surface_morphology_proposal(
                mesh,
                feature_nodes,
                dep_tensor_block,
                args.morph_band_length_m,
                args.morph_normal_weight,
                args.morph_shear_weight,
            )
            mesh_scale = apply_local_ale_surface_update(
                mesh,
                feature_nodes,
                dh_block,
                decay_length=args.morph_decay_length_m,
                fixed_nodes=fixed_mesh_nodes,
                max_move=args.max_surface_move_fraction * mesh.hbar_tip,
                min_area_fraction=args.min_area_fraction,
            )
            if mesh_scale > 0.0:
                root_xy = local_root_xy(mesh, feature_nodes)
                rebuild_mesh_geometry(mesh, root_xy)
                patch.update_geometry(mesh, root_xy)

        Umax2, Umin2, u_zero2, F02, _ = affine_stress_control_displacements(
            mesh, bnd, mat, Dmat, ep_gp, sigma_max, sigma_min, u_zero
        )
        _, _, s1_res, _ = stress_state_intact(mesh, u_zero2, ep_gp, Dmat, mat)
        residual_node = project_gp_to_nodes(mesh, s1_res)
        last_residual = residual_node.copy()

        epsp_node, rho_node, P, Dloc = project_plastic_state(
            mesh, epsp_acc_gp, rho_gp, args.epsp_shield_scale, args.epsp_damage_scale
        )
        sigma_back = args.sigma_back_max_GPa * 1e9 * P
        state_shift = Gsh * P - args.Gstored_eV * Dloc
        hist_post = cycle_stress_histories(
            mesh, bnd, mat, Dmat, ep_gp, Umax2, Umin2, args.hazard_n_phase, u_zero2
        )
        ep_node_post = project_gp_to_nodes(mesh, ep_gp)
        _, _, _, bond_amp_post, point_amp_post = patch.solve_local_mechanics(
            pd_state, hist_post["u_max"], ep_node_post
        )
        sigma_combined = np.concatenate([hist_pre["sigma_node"], hist_post["sigma_node"]], axis=0)
        point_amp = 0.5 * (point_amp_pre + point_amp_post)
        bond_amp = 0.5 * (bond_amp_pre + bond_amp_post)
        diag = patch.update(
            pd_state,
            crack,
            sigma_combined,
            args.T,
            args.frequency_Hz,
            dN,
            cycles,
            state_shift,
            sigma_back,
            chi,
            P,
            point_amp,
            bond_amp,
        )
        cycles += dN
        last_hist = hist_post
        last_diag = diag

        root_radius = local_root_radius(mesh, feature_nodes)
        row = {
            "block": ib,
            "case": case_name,
            "sigma_a_MPa": sigma_a_MPa,
            "cycles_total": cycles,
            "dN": dN,
            "Umax_m": Umax2,
            "Umin_m": Umin2,
            "F0_residual_N_per_m": F02,
            "dep_eq_cycle_max": max_dep_cycle,
            "dep_eq_block_max": float(np.max(dep_eq_block)),
            "rho_max_m2": float(np.max(rho_gp)),
            "rho_mean_m2": float(np.mean(rho_gp)),
            "epsp_acc_max": float(np.max(epsp_acc_gp)),
            "residual_sigma1_max_Pa": float(np.max(residual_node)),
            "sigma1_cycle_max_Pa": float(np.max(hist_post["s1_node"])),
            "P_max": float(np.max(P)),
            "Dloc_max": float(np.max(Dloc)),
            "pd_hit_memory_max": diag.max_hit_memory,
            "pd_completion_max": diag.max_completion,
            "pd_embryo_max": diag.max_embryo,
            "pd_stable_max": diag.max_stable,
            "pd_growth_max": diag.max_growth,
            "pd_expected_embryos": diag.expected_embryos,
            "pd_expected_births_cumulative": diag.expected_births_cumulative,
            "pd_expected_stable": diag.expected_stable,
            "pd_realized_embryos": diag.realized_embryos,
            "pd_realized_births_cumulative": diag.realized_births_cumulative,
            "pd_realized_stable": diag.realized_stable,
            "pd_bond_damage_max": diag.max_bond_damage,
            "pd_broken_bonds": diag.broken_bonds,
            "pd_connected_extent_m": diag.connected_extent_m,
            "pd_connected_bonds": diag.connected_bonds,
            "pd_amplification_max": diag.max_pd_amplification,
            "pd_effective_stress_max_Pa": diag.max_effective_stress_Pa,
            "pd_max_rate_per_cycle": diag.max_rate_per_cycle,
            "pd_hit_rate_max_s": diag.max_hit_rate_s,
            "pd_birth_rate_max_per_cycle": diag.max_birth_rate_per_cycle,
            "pd_expected_candidate_sites": diag.expected_candidate_sites,
            "pd_realized_candidate_sites": diag.realized_candidate_sites,
            "root_x_m": root_xy[0],
            "root_y_m": root_xy[1],
            "root_radius_m": root_radius,
            "root_radius_over_initial": root_radius / max(root_radius0, 1e-30) if np.isfinite(root_radius) else np.nan,
            "surface_move_max_m": float(np.max(np.abs(dh_block))) if dh_block.size else 0.0,
            "surface_shear_gamma_nt_max": float(np.max(np.abs(gamma_nt))) if gamma_nt.size else 0.0,
            "ale_accept_scale": mesh_scale,
            "plastic_work_J_per_m": Wp_total,
            "cycles_first_embryo": pd_state.cycles_first_embryo if pd_state.cycles_first_embryo is not None else np.nan,
            "cycles_first_stable": pd_state.cycles_first_stable if pd_state.cycles_first_stable is not None else np.nan,
            "cycles_first_expected_embryo": pd_state.cycles_first_expected_embryo if pd_state.cycles_first_expected_embryo is not None else np.nan,
            "cycles_first_expected_stable": pd_state.cycles_first_expected_stable if pd_state.cycles_first_expected_stable is not None else np.nan,
            "cycles_first_softening": pd_state.cycles_first_softening if pd_state.cycles_first_softening is not None else np.nan,
            "cycles_connected": pd_state.cycles_connected if pd_state.cycles_connected is not None else np.nan,
        }
        rows.append(row)

        if args.print_every and ib % args.print_every == 0:
            print(
                f"STATEFUL_PD {case_name} sigma_a={sigma_a_MPa:g}MPa block={ib} "
                f"N={cycles:.3e} dN={dN:.3g} dep={row['dep_eq_block_max']:.2e} "
                f"Ebirth={diag.expected_births_cumulative:.2g} Estab={diag.expected_stable:.2g} "
                f"Rbirth={diag.realized_births_cumulative:d} Rstab={diag.realized_stable:d} "
                f"omega={diag.max_bond_damage:.3g} broken={diag.broken_bonds} "
                f"a_conn={diag.connected_extent_m*1e6:.1f}um"
            )

        if ib == 0 or (args.snapshot_every > 0 and ib % args.snapshot_every == 0):
            patch.plot_snapshot(
                pd_state,
                outdir / f"pd_patch_block_{ib:05d}.png",
                title=f"{case_name}, sigma_a={sigma_a_MPa:g} MPa, N={cycles:.3e}",
            )
            _plot_fem_fields(
                mesh,
                [
                    ("accumulated eps_p", epsp_node),
                    ("rho (m^-2)", rho_node),
                    ("residual sigma1 (MPa)", residual_node * 1e-6),
                ],
                outdir / f"fem_fields_block_{ib:05d}.png",
                root_xy,
                feature_nodes,
            )
        u = hist_post["u_end"].copy()

    _write_csv(outdir / "sn_stateful_pd_history.csv", rows)
    patch.plot_snapshot(pd_state, outdir / "pd_patch_final.png", title=f"final N={cycles:.3e}")
    patch.plot_initiation_diagnostics(
        pd_state,
        outdir / "pd_initiation_diagnostics_final.png",
        title=f"{case_name}, sigma_a={sigma_a_MPa:g} MPa, final N={cycles:.3e}",
    )
    epsp_node, rho_node, P, Dloc = project_plastic_state(
        mesh, epsp_acc_gp, rho_gp, args.epsp_shield_scale, args.epsp_damage_scale
    )
    _plot_fem_fields(
        mesh,
        [
            ("accumulated eps_p", epsp_node),
            ("rho (m^-2)", rho_node),
            ("residual sigma1 (MPa)", last_residual * 1e-6),
        ],
        outdir / "fem_fields_final.png",
        root_xy,
        feature_nodes,
    )
    np.savez_compressed(
        outdir / "pd_state_final.npz",
        global_nodes=patch.global_nodes,
        xy=patch.xy,
        bonds=patch.bonds,
        boundary=patch.boundary,
        initiation_weight=patch.initiation_weight,
        mean_candidate_sites=patch.mean_candidate_sites,
        available=pd_state.available,
        embryo=pd_state.embryo,
        stable=pd_state.stable,
        inactive=pd_state.inactive,
        candidate_sites=pd_state.candidate_sites,
        available_sites=pd_state.available_sites,
        embryo_sites=pd_state.embryo_sites,
        stable_sites=pd_state.stable_sites,
        inactive_sites=pd_state.inactive_sites,
        born_sites_cumulative=pd_state.born_sites_cumulative,
        healed_sites_cumulative=pd_state.healed_sites_cumulative,
        hit_memory=pd_state.hit_memory,
        completion=pd_state.completion,
        growth=pd_state.growth,
        born_cumulative=pd_state.born_cumulative,
        healed_cumulative=pd_state.healed_cumulative,
        bond_damage=pd_state.bond_damage,
        point_amplification=getattr(patch, "last_point_amp", np.ones(len(patch.xy))),
        bond_amplification=getattr(patch, "last_bond_amp", np.ones(len(patch.bonds))),
        hit_rate_s=getattr(patch, "last_rates", {}).get("hit_rate_s", np.zeros(len(patch.xy))),
        hit_intensity_per_cycle=getattr(patch, "last_rates", {}).get("mu_hit", np.zeros(len(patch.xy))),
        birth_intensity_per_cycle=getattr(patch, "last_rates", {}).get("mu_birth", np.zeros(len(patch.xy))),
        effective_stress_Pa=getattr(patch, "last_rates", {}).get("smax", np.zeros(len(patch.xy))),
    )

    final_connected_extent, final_connected_bonds = patch.connected_crack(pd_state)
    status = "connected_crack" if pd_state.cycles_connected is not None else "right_censored"
    summary = {
        "model": MODEL_ID,
        "coupling": "one_way_FEM_to_PD_with_local_PD_redistribution",
        "case": case_name,
        "sigma_a_MPa": sigma_a_MPa,
        "T_K": args.T,
        "R": args.R,
        "frequency_Hz": args.frequency_Hz,
        "cycles_total": cycles,
        "cycles_first_embryo": pd_state.cycles_first_embryo,
        "cycles_first_stable": pd_state.cycles_first_stable,
        "cycles_first_expected_embryo": pd_state.cycles_first_expected_embryo,
        "cycles_first_expected_stable": pd_state.cycles_first_expected_stable,
        "cycles_first_softening": pd_state.cycles_first_softening,
        "cycles_connected": pd_state.cycles_connected,
        "status": status,
        "pd_points": len(patch.xy),
        "pd_bonds": len(patch.bonds),
        "pd_horizon_m": pd_cfg.horizon_m,
        "pd_patch_radius_m": pd_cfg.patch_radius_m,
        "pd_initiation_radius_m": pd_cfg.initiation_radius_m,
        "pd_initiation_taper_m": pd_cfg.initiation_taper_m,
        "pd_initiation_back_extent_m": pd_cfg.initiation_back_extent_m,
        "pd_established_extent_criterion_m": pd_cfg.established_extent_m,
        "pd_connected_extent_final_m": float(final_connected_extent),
        "pd_connected_bonds_final": int(final_connected_bonds),
        "pd_broken_bonds_final": int(np.count_nonzero(pd_state.bond_damage >= pd_cfg.broken_damage)),
        "pd_bond_damage_final_max": float(np.max(pd_state.bond_damage)),
        "pd_expected_candidate_sites": float(np.sum(patch.mean_candidate_sites)),
        "pd_expected_embryos_final": float(np.sum(pd_state.embryo * patch.mean_candidate_sites)),
        "pd_expected_births_cumulative_final": float(np.sum(pd_state.born_cumulative * patch.mean_candidate_sites)),
        "pd_expected_stable_final": float(np.sum(pd_state.stable * patch.mean_candidate_sites)),
        "pd_candidate_sites_realized": int(np.sum(pd_state.candidate_sites)),
        "pd_realized_births_cumulative_final": int(np.sum(pd_state.born_sites_cumulative)),
        "pd_realized_embryos_final": int(np.sum(pd_state.embryo_sites)),
        "pd_realized_stable_final": int(np.sum(pd_state.stable_sites)),
        "pd_hit_memory_s": float(args.hit_memory_s),
        "pd_hit_memory_cycles": float(args.frequency_Hz * args.hit_memory_s),
        "pd_hit_process": "continuous_time_Poisson",
        "pd_amplification_damage_scale": float(pd_cfg.amplification_damage_scale),
        "pd_softening_damage_threshold": float(pd_cfg.softening_damage),
        "root_radius_initial_m": root_radius0,
        "root_radius_final_m": local_root_radius(mesh, feature_nodes),
        "epsp_acc_final_max": float(np.max(epsp_acc_gp)),
        "rho_final_max_m2": float(np.max(rho_gp)),
        "residual_sigma1_final_max_Pa": float(np.max(last_residual)),
        "plastic_work_final_J_per_m": Wp_total,
        "history_csv": str(outdir / "sn_stateful_pd_history.csv"),
    }
    with (outdir / "summary.json").open("w") as f:
        json.dump(summary, f, indent=2, sort_keys=True)
    return summary


def run_sweep(args):
    Path(args.out).mkdir(parents=True, exist_ok=True)
    summaries = []
    for case in args.cases:
        for s in args.sigma_a_MPa:
            print(f"=== STATEFUL PD case={case} sigma_a={s:g} MPa ===")
            summaries.append(run_case_stress(args, case, float(s)))
    _write_csv(Path(args.out) / "sn_stateful_pd_summary.csv", summaries)
    return summaries


def build_parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", default="runs/sn_stateful_pd_pilot")
    p.add_argument("--cases", nargs="+", choices=["no_shield", "shielded"], default=["no_shield", "shielded"])
    p.add_argument("--T", type=float, default=300.0)
    p.add_argument("--sigma-a-MPa", nargs="+", type=float, default=[500, 600, 700], dest="sigma_a_MPa")
    p.add_argument("--R", type=float, default=0.1)
    p.add_argument("--frequency-Hz", type=float, default=1000.0, dest="frequency_Hz")
    p.add_argument("--cycles-max", type=float, default=1e9, dest="cycles_max")
    p.add_argument("--block-cycles", type=float, default=1e7, dest="block_cycles")
    p.add_argument("--min-block-cycles", type=float, default=1e-6, dest="min_block_cycles")
    p.add_argument("--max-blocks", type=int, default=3000, dest="max_blocks")
    p.add_argument("--target-dep-eq-block", type=float, default=2e-4, dest="target_dep_eq_block")
    p.add_argument("--target-rho-rel-block", type=float, default=0.05, dest="target_rho_rel_block")
    p.add_argument("--target-hit-memory", type=float, default=float("inf"), dest="target_hit_memory")
    p.add_argument("--max-transition-probability", type=float, default=0.08, dest="max_transition_probability")
    p.add_argument("--plastic-n-phase", type=int, default=12, dest="plastic_n_phase")
    p.add_argument("--hazard-n-phase", type=int, default=16, dest="hazard_n_phase")

    p.add_argument("--Lx", type=float, default=2e-3)
    p.add_argument("--Ly", type=float, default=4e-3)
    p.add_argument("--notch-depth-m", type=float, default=0.15e-3, dest="notch_depth_m")
    p.add_argument("--notch-half-height-m", type=float, default=0.30e-3, dest="notch_half_height_m")
    p.add_argument("--nx", type=int, default=36)
    p.add_argument("--ny", type=int, default=72)
    p.add_argument("--jitter", type=float, default=0.08)
    p.add_argument("--root-h-fine", type=float, default=30e-6, dest="root_h_fine")
    p.add_argument("--seed", type=int, default=42)

    p.add_argument("--E-GPa", type=float, default=410.0, dest="E_GPa")
    p.add_argument("--nu", type=float, default=0.28)
    p.add_argument("--b-m", type=float, default=2.74e-10, dest="b_m")
    p.add_argument("--Tm-K", type=float, default=3695.0, dest="Tm_K")
    p.add_argument("--rho0", type=float, default=1e12)
    p.add_argument("--rho-floor", type=float, default=1e8, dest="rho_floor")
    p.add_argument("--rho-cap", type=float, default=1e17, dest="rho_cap")
    p.add_argument("--k-store", type=float, default=np.sqrt(2.0), dest="k_store")
    p.add_argument("--k-dyn", type=float, default=1.0, dest="k_dyn")

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

    p.add_argument("--enable-geometry-evolution", action="store_true", default=True, dest="enable_geometry_evolution")
    p.add_argument("--disable-geometry-evolution", action="store_false", dest="enable_geometry_evolution")
    p.add_argument("--morph-band-length-m", type=float, default=100e-6, dest="morph_band_length_m")
    p.add_argument("--morph-normal-weight", type=float, default=0.25, dest="morph_normal_weight")
    p.add_argument("--morph-shear-weight", type=float, default=1.0, dest="morph_shear_weight")
    p.add_argument("--morph-decay-length-m", type=float, default=150e-6, dest="morph_decay_length_m")
    p.add_argument("--target-surface-move-fraction", type=float, default=0.10, dest="target_surface_move_fraction")
    p.add_argument("--max-surface-move-fraction", type=float, default=0.20, dest="max_surface_move_fraction")
    p.add_argument("--min-area-fraction", type=float, default=0.15, dest="min_area_fraction")

    p.add_argument("--pd-patch-radius-m", type=float, default=0.45e-3, dest="pd_patch_radius_m")
    p.add_argument("--pd-horizon-m", type=float, default=90e-6, dest="pd_horizon_m")
    p.add_argument("--pd-boundary-shell-m", type=float, default=100e-6, dest="pd_boundary_shell_m")
    p.add_argument("--pd-residual-stiffness", type=float, default=1e-7, dest="pd_residual_stiffness")
    p.add_argument("--pd-initiation-radius-m", type=float, default=240e-6, dest="pd_initiation_radius_m")
    p.add_argument("--pd-initiation-taper-m", type=float, default=60e-6, dest="pd_initiation_taper_m")
    p.add_argument("--pd-initiation-back-extent-m", type=float, default=60e-6, dest="pd_initiation_back_extent_m")
    p.add_argument("--pd-amplification-cap", type=float, default=4.0, dest="pd_amplification_cap")
    p.add_argument("--pd-amplification-damage-scale", type=float, default=0.05, dest="pd_amplification_damage_scale")
    p.add_argument("--pd-seed", type=int, default=None, dest="pd_seed",
                   help="seed for the discrete candidate-site realization; defaults to --seed")
    p.add_argument("--site-density-m2", type=float, default=5e10, dest="site_density_m2")
    p.add_argument("--hit-count", type=float, default=3.0, dest="hit_count")
    p.add_argument("--hit-memory-s", type=float, default=1e-6, dest="hit_memory_s")
    p.add_argument("--birth-scale", type=float, default=1.0, dest="birth_scale")
    p.add_argument("--nu-stabilize-s", type=float, default=5e2, dest="nu_stabilize_s")
    p.add_argument("--nu-heal-s", type=float, default=2e2, dest="nu_heal_s")
    p.add_argument("--stabilize-stress-GPa", type=float, default=1.4, dest="stabilize_stress_GPa")
    p.add_argument("--stabilize-width-GPa", type=float, default=0.25, dest="stabilize_width_GPa")
    p.add_argument("--stabilize-plastic-gain", type=float, default=2.0, dest="stabilize_plastic_gain")
    p.add_argument("--heal-return-fraction", type=float, default=0.9, dest="heal_return_fraction")
    p.add_argument("--nu-grow-s", type=float, default=2e-2, dest="nu_grow_s")
    p.add_argument("--grow-stress-GPa", type=float, default=1.2, dest="grow_stress_GPa")
    p.add_argument("--grow-width-GPa", type=float, default=0.25, dest="grow_width_GPa")
    p.add_argument("--stable-count-scale", type=float, default=2.0, dest="stable_count_scale")
    p.add_argument("--nu-link-s", type=float, default=8e-3, dest="nu_link_s")
    p.add_argument("--link-stress-GPa", type=float, default=1.1, dest="link_stress_GPa")
    p.add_argument("--link-width-GPa", type=float, default=0.20, dest="link_width_GPa")
    p.add_argument("--link-orientation-power", type=float, default=2.0, dest="link_orientation_power")
    p.add_argument("--link-orientation-floor", type=float, default=0.05, dest="link_orientation_floor")
    p.add_argument("--neighbor-link-gain", type=float, default=1.5, dest="neighbor_link_gain")
    p.add_argument("--broken-damage", type=float, default=0.95, dest="broken_damage")
    p.add_argument("--softening-damage", type=float, default=1e-3, dest="softening_damage",
                   help="cohesive-damage level used to record first measurable softening")
    p.add_argument("--root-seed-radius-m", type=float, default=120e-6, dest="root_seed_radius_m")
    p.add_argument("--established-extent-m", type=float, default=240e-6, dest="established_extent_m")

    p.add_argument("--snapshot-every", type=int, default=25, dest="snapshot_every")
    p.add_argument("--print-every", type=int, default=1, dest="print_every")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    run_sweep(args)


if __name__ == "__main__":
    main()
