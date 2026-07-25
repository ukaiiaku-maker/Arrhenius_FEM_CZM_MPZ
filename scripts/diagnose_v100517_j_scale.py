#!/usr/bin/env python3
"""Controlled elastic diagnostic for the FEM domain-J toughness scale.

The production solver is not modified.  This script solves the same finite
specimen in a deliberately elastic sharp-wake calibration and compares:

1. the current domain integral using the tension-only energy field returned by
   ``assemble_mechanics``;
2. the same domain integral using the full stored elastic energy density
   ``0.5 * eps_e : sigma``; and
3. the stored-energy drop at fixed imposed displacement for a small virtual
   crack extension.

The result is diagnostic evidence only.  It does not rescale K, change a
parameter row, or alter the production crack-growth law.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

import numpy as np

from arrhenius_fracture.config import make_emergent_config
from arrhenius_fracture.crystal import cubic_plane_strain_D
from arrhenius_fracture.fem import (
    assemble_mechanics,
    elastic_energy_densities,
    solve_dirichlet,
)
from arrhenius_fracture.j_integral import compute_J_integral
from arrhenius_fracture.mesh import make_boundary_data, make_tri_mesh


def _damage_for_tip(mesh, geom, a_tip: float) -> np.ndarray:
    x = mesh.nodes[:, 0]
    y = mesh.nodes[:, 1]
    d = np.zeros(mesh.nn)
    d[(x <= a_tip) & (np.abs(y) <= geom.notch_half_thickness)] = 1.0
    return d


def _solve_elastic_state(mesh, bnd, geom, D, mat, a_tip: float, Uapp: float) -> dict[str, Any]:
    d = _damage_for_tip(mesh, geom, a_tip)
    u = np.zeros(mesh.ndof)
    ep_gp = np.zeros((3, mesh.ne))
    rho_gp = np.full(mesh.ne, 5.0e12)
    Kmat, Rint, sigma_gp, seq_gp, s1_gp, psi_tension = assemble_mechanics(
        mesh, u, ep_gp, rho_gp, d, D, mat
    )
    u, Ftop = solve_dirichlet(Kmat, Rint, u, bnd, 0.5 * Uapp, -0.5 * Uapp)
    Kmat, Rint, sigma_gp, seq_gp, s1_gp, psi_tension = assemble_mechanics(
        mesh, u, ep_gp, rho_gp, d, D, mat
    )
    psi_stored, psi_undegraded = elastic_energy_densities(
        mesh, u, ep_gp, sigma_gp, D
    )
    return {
        "u": u,
        "d": d,
        "sigma_gp": sigma_gp,
        "psi_tension": psi_tension,
        "psi_stored": psi_stored,
        "psi_undegraded": psi_undegraded,
        "stored_energy_J_per_m": float(np.sum(psi_stored * mesh.area_e)),
        "Ftop_N_per_m": float(Ftop),
    }


def _relative_error(value: float, reference: float) -> float:
    return abs(value - reference) / max(abs(reference), 1.0e-30)


def _relative_spread(values: list[float]) -> float:
    array = np.asarray(values, dtype=float)
    mean = float(np.mean(np.abs(array)))
    return float((np.max(array) - np.min(array)) / max(mean, 1.0e-30))


def run_diagnostic(
    out: Path,
    tip_h_um: list[float],
    contour_outer_um: list[float],
    delta_a_um: float,
    Uapp_um: float,
    theta_deg: float,
    nx: int,
    ny: int,
) -> dict[str, Any]:
    cfg = make_emergent_config()
    mat = cfg.material
    geom = cfg.geometry
    D = cubic_plane_strain_D(523.0e9, 203.0e9, 160.0e9, theta_deg)
    delta_a = float(delta_a_um) * 1.0e-6
    Uapp = float(Uapp_um) * 1.0e-6
    rows: list[dict[str, Any]] = []
    mesh_summaries: list[dict[str, Any]] = []

    for h_um in tip_h_um:
        cfg.mesh.nx = int(nx)
        cfg.mesh.ny = int(ny)
        cfg.mesh.tip_h_fine = float(h_um) * 1.0e-6
        cfg.mesh.tip_ratio = 1.15
        mesh = make_tri_mesh(geom, cfg.mesh, seed=42)
        bnd = make_boundary_data(mesh, geom)
        a0 = float(geom.a0)
        state0 = _solve_elastic_state(mesh, bnd, geom, D, mat, a0, Uapp)
        state1 = _solve_elastic_state(mesh, bnd, geom, D, mat, a0 + delta_a, Uapp)
        J_fd = max(
            (state0["stored_energy_J_per_m"] - state1["stored_energy_J_per_m"]) / delta_a,
            0.0,
        )
        K_fd = float(np.sqrt(max(J_fd, 0.0) * mat.Eprime))
        current_values: list[float] = []
        stored_values: list[float] = []
        current_errors: list[float] = []
        stored_errors: list[float] = []
        crack_segments = [
            (np.array([0.0, 0.0]), np.array([a0, 0.0]))
        ]
        h_local = mesh.hbar_tip if mesh.hbar_tip > 0.0 else mesh.hbar
        for outer_um in contour_outer_um:
            ell = float(outer_um) * 1.0e-6 / 8.0
            _, _, info_current = compute_J_integral(
                mesh,
                state0["u"],
                state0["sigma_gp"],
                state0["psi_tension"],
                state0["d"],
                np.array([a0, 0.0]),
                np.array([1.0, 0.0]),
                mat,
                ell=ell,
                crack_segments=crack_segments,
                exclude_radius=2.0 * h_local,
            )
            _, _, info_stored = compute_J_integral(
                mesh,
                state0["u"],
                state0["sigma_gp"],
                state0["psi_stored"],
                state0["d"],
                np.array([a0, 0.0]),
                np.array([1.0, 0.0]),
                mat,
                ell=ell,
                crack_segments=crack_segments,
                exclude_radius=2.0 * h_local,
            )
            J_current = abs(float(info_current.get("J_signed", info_current["J"])))
            J_stored = abs(float(info_stored.get("J_signed", info_stored["J"])))
            K_current = float(np.sqrt(max(J_current, 0.0) * mat.Eprime))
            K_stored = float(np.sqrt(max(J_stored, 0.0) * mat.Eprime))
            current_values.append(J_current)
            stored_values.append(J_stored)
            current_errors.append(_relative_error(J_current, J_fd))
            stored_errors.append(_relative_error(J_stored, J_fd))
            rows.append(
                {
                    "tip_h_um": float(h_um),
                    "mesh_nodes": int(mesh.nn),
                    "mesh_elements": int(mesh.ne),
                    "hbar_tip_um": float(h_local * 1.0e6),
                    "contour_outer_um": float(outer_um),
                    "J_current_tension_only_J_m2": J_current,
                    "K_current_tension_only_MPa_sqrt_m": K_current / 1.0e6,
                    "J_full_stored_J_m2": J_stored,
                    "K_full_stored_MPa_sqrt_m": K_stored / 1.0e6,
                    "J_fixed_displacement_energy_drop_J_m2": J_fd,
                    "K_fixed_displacement_energy_drop_MPa_sqrt_m": K_fd / 1.0e6,
                    "current_over_stored_J": J_current / max(J_stored, 1.0e-30),
                    "current_relative_error_to_energy_drop": current_errors[-1],
                    "stored_relative_error_to_energy_drop": stored_errors[-1],
                    "current_active_elements": int(info_current.get("n_active_elements", 0)),
                    "stored_active_elements": int(info_stored.get("n_active_elements", 0)),
                }
            )
        mesh_summaries.append(
            {
                "tip_h_um": float(h_um),
                "hbar_tip_um": float(h_local * 1.0e6),
                "mesh_nodes": int(mesh.nn),
                "mesh_elements": int(mesh.ne),
                "J_fd_J_m2": J_fd,
                "K_fd_MPa_sqrt_m": K_fd / 1.0e6,
                "J_current_median_J_m2": float(np.median(current_values)),
                "J_stored_median_J_m2": float(np.median(stored_values)),
                "current_contour_relative_spread": _relative_spread(current_values),
                "stored_contour_relative_spread": _relative_spread(stored_values),
                "current_median_relative_error_to_energy_drop": float(np.median(current_errors)),
                "stored_median_relative_error_to_energy_drop": float(np.median(stored_errors)),
            }
        )

    current_error = float(np.median([
        item["current_median_relative_error_to_energy_drop"] for item in mesh_summaries
    ]))
    stored_error = float(np.median([
        item["stored_median_relative_error_to_energy_drop"] for item in mesh_summaries
    ]))
    current_spread = float(max(
        item["current_contour_relative_spread"] for item in mesh_summaries
    ))
    stored_spread = float(max(
        item["stored_contour_relative_spread"] for item in mesh_summaries
    ))
    stored_mesh_values = [item["J_stored_median_J_m2"] for item in mesh_summaries]
    stored_mesh_spread = _relative_spread(stored_mesh_values) if len(stored_mesh_values) > 1 else 0.0

    if stored_error < 0.5 * current_error and stored_error <= 0.25:
        conclusion = "full_stored_energy_materially_improves_J_agreement"
    elif current_error <= 0.20 and stored_error <= 0.20:
        conclusion = "no_large_energy_density_bias_detected"
    else:
        conclusion = "inconclusive_requires_refined_J_calibration"

    acceptance = {
        "stored_median_relative_error_le_0p25": stored_error <= 0.25,
        "stored_contour_relative_spread_le_0p15": stored_spread <= 0.15,
        "stored_mesh_relative_spread_le_0p15": stored_mesh_spread <= 0.15,
    }
    acceptance["pass"] = all(acceptance.values())
    payload = {
        "schema": "v10.0.5.17_J_scale_diagnostic_v1",
        "scope": "controlled elastic sharp-wake calibration; production physics unchanged",
        "crystal_theta_deg": float(theta_deg),
        "crystal_C11_Pa": 523.0e9,
        "crystal_C12_Pa": 203.0e9,
        "crystal_C44_Pa": 160.0e9,
        "Uapp_um": float(Uapp_um),
        "delta_a_um": float(delta_a_um),
        "tip_h_um_requested": [float(value) for value in tip_h_um],
        "contour_outer_um": [float(value) for value in contour_outer_um],
        "current_energy_definition": "tension-only psi returned by assemble_mechanics",
        "alternative_energy_definition": "full stored 0.5*eps_e:sigma",
        "reference_definition": "stored-energy drop at fixed imposed displacement divided by delta_a",
        "current_median_relative_error_to_energy_drop": current_error,
        "stored_median_relative_error_to_energy_drop": stored_error,
        "current_max_contour_relative_spread": current_spread,
        "stored_max_contour_relative_spread": stored_spread,
        "stored_mesh_relative_spread": stored_mesh_spread,
        "conclusion": conclusion,
        "acceptance": acceptance,
        "mesh_summaries": mesh_summaries,
    }

    out.mkdir(parents=True, exist_ok=True)
    csv_path = out / "J_scale_contour_mesh_diagnostic.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    json_path = out / "J_scale_diagnostic_summary.json"
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload, indent=2, sort_keys=True))
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--tip-h-um", type=float, nargs="+", default=[5.0, 2.5])
    parser.add_argument("--contour-outer-um", type=float, nargs="+", default=[100.0, 160.0, 240.0])
    parser.add_argument("--delta-a-um", type=float, default=10.0)
    parser.add_argument("--Uapp-um", type=float, default=20.0)
    parser.add_argument("--theta-deg", type=float, default=30.0)
    parser.add_argument("--nx", type=int, default=36)
    parser.add_argument("--ny", type=int, default=72)
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args()
    payload = run_diagnostic(
        args.out.expanduser().resolve(),
        list(args.tip_h_um),
        list(args.contour_outer_um),
        args.delta_a_um,
        args.Uapp_um,
        args.theta_deg,
        args.nx,
        args.ny,
    )
    if args.strict and not payload["acceptance"]["pass"]:
        raise SystemExit("J-scale diagnostic failed the strict convergence/agreement gate")


if __name__ == "__main__":
    main()
