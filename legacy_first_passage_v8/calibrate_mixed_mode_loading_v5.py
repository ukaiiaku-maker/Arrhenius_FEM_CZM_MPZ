#!/usr/bin/env python3
"""Anisotropic traction-phase and directional-factor calibration for v5."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

MODEL_ID = "FEM_CZM_mixed_mode_calibration_v5_anisotropic_calibrated_tip"


def vals(s):
    return [float(x) for x in str(s).replace(",", " ").split() if x]


def err(v, t):
    return (float(v)-float(t)+180.0) % 360.0 - 180.0


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--out", default="runs/mixed_mode_fem_czm_v5_anisotropic_calibration")
    p.add_argument("--target-psi-deg", default="-60 -45 -30 -15 0 15 30 45 60")
    p.add_argument("--U-cal-m", type=float, default=2e-7)
    p.add_argument("--nx", type=int, default=24)
    p.add_argument("--ny", type=int, default=48)
    p.add_argument("--tip-h-fine", type=float, default=3e-6)
    p.add_argument("--tip-ratio", type=float, default=1.25)
    p.add_argument("--mesh-seed", type=int, default=42)
    p.add_argument("--crystal-theta-deg", type=float, default=45.0)
    p.add_argument("--crystal-C11", type=float, default=523e9)
    p.add_argument("--crystal-C12", type=float, default=203e9)
    p.add_argument("--crystal-C44", type=float, default=160e9)
    p.add_argument("--cleave-gamma-aniso", type=float, default=0.3)
    p.add_argument("--traction-probe-radius-m", type=float, default=10e-6)
    p.add_argument("--traction-annulus-half-width", type=float, default=0.45)
    p.add_argument("--traction-sector-half-angle-deg", type=float, default=40.0)
    p.add_argument("--psi-tol-deg", type=float, default=0.75)
    p.add_argument("--basis-condition-max", type=float, default=1e8)
    p.add_argument("--shear-emission-weight", type=float, default=1.0)
    p.add_argument("--directional-factor-max", type=float, default=5.0)
    a = p.parse_args()

    from arrhenius_fracture.config import make_emergent_config
    from arrhenius_fracture.mesh import make_tri_mesh, make_boundary_data
    from arrhenius_fracture.fem import assemble_mechanics, solve_dirichlet
    from arrhenius_fracture.crystal import cubic_plane_strain_D, zener_ratio
    from arrhenius_fracture.j_integral import compute_J_integral
    from arrhenius_fracture.mixed_mode_first_passage_v5 import (
        AnisotropicCalibratedTipContext, _mixed_solve_factory,
        process_zone_traction_probe, directional_shape_metrics,
        directional_drive_factors, shear_sign_from_basis,
        loading_angle_from_response_basis, traction_phase_deg,
        energy_matrix_from_basis,
    )

    cfg = make_emergent_config()
    cfg.mesh.nx = a.nx
    cfg.mesh.ny = a.ny
    cfg.mesh.tip_h_fine = a.tip_h_fine
    cfg.mesh.tip_ratio = a.tip_ratio
    mesh = make_tri_mesh(cfg.geometry, cfg.mesh, seed=a.mesh_seed)
    bnd = make_boundary_data(mesh, cfg.geometry)
    mat = cfg.material
    D = cubic_plane_strain_D(a.crystal_C11, a.crystal_C12, a.crystal_C44,
                             a.crystal_theta_deg)
    u0 = np.zeros(mesh.ndof)
    ep = np.zeros((3, mesh.ne))
    rho = np.zeros(mesh.ne)
    d = np.zeros(mesh.nn)
    d[bnd.notch_nodes] = 1
    K, R, *_ = assemble_mechanics(mesh, u0, ep, rho, d, D, mat)
    cache = {}

    def evaluate(alpha):
        key = round(float(alpha), 10)
        if key in cache:
            return dict(cache[key])
        ctx = AnisotropicCalibratedTipContext(
            alpha, 0.0, a.crystal_theta_deg, a.cleave_gamma_aniso,
            a.traction_probe_radius_m, a.traction_annulus_half_width,
            a.traction_sector_half_angle_deg)
        solve = _mixed_solve_factory(solve_dirichlet, ctx)
        u, F = solve(K, R, u0.copy(), bnd, 0.5*a.U_cal_m, -0.5*a.U_cal_m)
        _, _, sig, _, _, psi = assemble_mechanics(mesh, u, ep, rho, d, D, mat)
        probe = process_zone_traction_probe(
            mesh, sig, d, np.array([cfg.geometry.a0, 0.0]), np.array([1.0, 0.0]),
            a.traction_probe_radius_m, a.traction_annulus_half_width,
            a.traction_sector_half_angle_deg)
        metrics = directional_shape_metrics(
            probe["stress_tensor"], a.crystal_theta_deg, a.cleave_gamma_aniso,
            np.array([1.0, 0.0])) if probe.get("reliable", False) else {"reliable": False}
        J, KJ, _ = compute_J_integral(
            mesh, u, sig, psi, d, np.array([cfg.geometry.a0, 0.0]),
            np.array([1.0, 0.0]), mat, ell=20e-6)
        row = {
            "loading_angle_deg": float(alpha),
            "sigma_nn_raw_Pa": probe.get("sigma_nn_Pa", np.nan),
            "tau_tn_raw_Pa": probe.get("tau_tn_Pa", np.nan),
            "probe_sigma1_Pa": metrics.get("sigma1_Pa", np.nan),
            "cleavage_shape": metrics.get("cleavage_shape", np.nan),
            "slip_shape": metrics.get("slip_shape", np.nan),
            "candidate_angle_deg": metrics.get("candidate_angle_deg", np.nan),
            "candidate_gamma_rel": metrics.get("candidate_gamma_rel", np.nan),
            "slip_system_name": metrics.get("slip_system_name", "none"),
            "probe_reliable": bool(probe.get("reliable", False) and metrics.get("reliable", False)),
            "probe_n_elements": probe.get("n_elements", 0),
            "J_J_per_m2": J,
            "KJ_reference_Pa_sqrt_m": KJ,
            "generalized_reaction_N": F,
        }
        cache[key] = dict(row)
        return row

    opening = evaluate(0.0)
    sliding = evaluate(90.0)
    Mraw = np.array([
        [opening["sigma_nn_raw_Pa"], sliding["sigma_nn_raw_Pa"]],
        [opening["tau_tn_raw_Pa"], sliding["tau_tn_raw_Pa"]],
    ])
    sign = shear_sign_from_basis(Mraw)
    M = np.diag([1.0, sign]) @ Mraw
    cond = float(np.linalg.cond(M))
    if not np.isfinite(cond) or cond > a.basis_condition_max:
        raise SystemExit(f"traction basis invalid cond={cond}")

    equal = evaluate(45.0)
    G = energy_matrix_from_basis(opening["J_J_per_m2"], sliding["J_J_per_m2"],
                                 equal["J_J_per_m2"], a.U_cal_m)
    ew = np.linalg.eigvalsh(G)
    energy_ok = bool(np.all(np.isfinite(ew)) and np.min(ew) > 0)

    # The zero-phase state is the normalization that recovers the original
    # Mode-I calibrated parameterization exactly.
    alpha0 = loading_angle_from_response_basis(M, 0.0)
    ref = evaluate(alpha0)
    ref_c = float(ref["cleavage_shape"])
    ref_s = max(float(ref["slip_shape"]), 0.0)
    if not np.isfinite(ref_c) or ref_c <= 1e-12:
        raise SystemExit(f"invalid Mode-I reference cleavage shape: {ref_c}")

    rows = []
    for target in vals(a.target_psi_deg):
        try:
            alpha = loading_angle_from_response_basis(M, target)
        except Exception as ex:
            rows.append({"target_psi_deg": target, "phase_converged": False,
                         "reason": f"basis_error:{ex}"})
            continue
        r = evaluate(alpha)
        phase = traction_phase_deg(r["sigma_nn_raw_Pa"], r["tau_tn_raw_Pa"], sign)
        e = err(phase, target)
        factors = directional_drive_factors(
            r["cleavage_shape"], r["slip_shape"], ref_c, ref_s,
            a.shear_emission_weight, a.directional_factor_max)
        row = {
            **r,
            **factors,
            "target_psi_deg": target,
            "loading_angle_deg": alpha,
            "traction_shear_sign": sign,
            "achieved_traction_phase_deg": phase,
            "traction_phase_error_deg": e,
            "phase_converged": bool(r["probe_reliable"] and abs(e) <= a.psi_tol_deg),
            "basis_condition": cond,
            "zener_ratio": zener_ratio(a.crystal_C11, a.crystal_C12, a.crystal_C44),
            "crystal_theta_deg": a.crystal_theta_deg,
            "reference_loading_angle_deg": alpha0,
            "reference_cleavage_shape": ref_c,
            "reference_slip_shape": ref_s,
            "shear_emission_weight": a.shear_emission_weight,
            "directional_factor_max": a.directional_factor_max,
            "energy_matrix_positive_definite": energy_ok,
            "response_11_Pa": M[0, 0], "response_12_Pa": M[0, 1],
            "response_21_Pa": M[1, 0], "response_22_Pa": M[1, 1],
            "energy_G11": G[0, 0], "energy_G12": G[0, 1], "energy_G22": G[1, 1],
        }
        rows.append(row)

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    csv_path = out/"mixed_mode_loading_calibration_v5.csv"
    with csv_path.open("w", newline="") as fp:
        cols = sorted({k for r in rows for k in r})
        w = csv.DictWriter(fp, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    (out/"anisotropic_calibrated_tip_basis.json").write_text(json.dumps({
        "model": MODEL_ID,
        "raw_traction_basis_Pa": Mraw.tolist(),
        "normalized_traction_basis_Pa": M.tolist(),
        "traction_shear_sign": sign,
        "basis_condition": cond,
        "energy_matrix": G.tolist(),
        "energy_eigenvalues": ew.tolist(),
        "energy_matrix_positive_definite": energy_ok,
        "reference_loading_angle_deg": alpha0,
        "reference_cleavage_shape": ref_c,
        "reference_slip_shape": ref_s,
        "shear_emission_weight": a.shear_emission_weight,
        "directional_factor_max": a.directional_factor_max,
        "crystal": {
            "C11": a.crystal_C11, "C12": a.crystal_C12, "C44": a.crystal_C44,
            "theta_deg": a.crystal_theta_deg,
            "zener_ratio": zener_ratio(a.crystal_C11, a.crystal_C12, a.crystal_C44),
        },
        "probe_radius_m": a.traction_probe_radius_m,
    }, indent=2))

    print("raw traction basis [Pa per calibration amplitude]:\n", Mraw)
    print("shear sign:", sign)
    print("normalized traction basis:\n", M)
    print("condition:", cond)
    print("Mode-I reference shapes:", {"cleavage": ref_c, "slip": ref_s})
    for r in rows:
        print({k: r.get(k) for k in (
            "target_psi_deg", "loading_angle_deg", "achieved_traction_phase_deg",
            "traction_phase_error_deg", "phase_converged",
            "cleavage_factor", "emission_factor", "candidate_angle_deg")})
    bad = [r for r in rows if not r.get("phase_converged", False)]
    if bad:
        raise SystemExit("v5 anisotropic calibration failed: " +
                         ", ".join(str(r.get("target_psi_deg")) for r in bad))
    print("wrote", csv_path)


if __name__ == "__main__":
    main()
