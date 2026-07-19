"""Production-path J parity utilities for v10.0.5.9.

The v10.0.5.8 benchmark established an independent fixed-grip elastic reference
for the finite killed notch.  This module records the *actual production*
initialization state immediately after the maximum-load FEM equilibrium solve and
re-evaluates that same state with:

* the exact production J-domain semantics;
* full stored degraded elastic energy;
* the legacy tensile-filtered energy;
* a no-exclusion ablation at identical fields.

No crack kinetics, cohesive state, plasticity law, or material parameter is
modified by these helpers.  The audit entry point disables plastic evolution only
for the one-step elastic probe.
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

from .config import JIntegralConfig
from .fem import elastic_energy_densities
from .j_integral import compute_J_integral

POINT_RELEASE = "10.0.5.9"
PROBE_JSON = "production_j_probe_v10_0_5_9.json"
SUMMARY_JSON = "production_j_parity_v10_0_5_9.json"
PROBE_CSV = "production_j_probe_summary_v10_0_5_9.csv"
CONTOUR_CSV = "production_j_contours_v10_0_5_9.csv"


def _finite(value: Any) -> bool:
    try:
        return bool(np.isfinite(float(value)))
    except (TypeError, ValueError):
        return False


def _relative_difference(a: float, b: float) -> float:
    scale = max(abs(float(a)), abs(float(b)), 1.0e-300)
    return abs(float(a) - float(b)) / scale


def _write_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = [dict(row) for row in rows]
    if not rows:
        path.write_text("")
        return
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _effective_killed_tip(mesh, d: np.ndarray, half_thickness_m: float) -> float:
    mask = (np.asarray(d, dtype=float) > 0.5) & (
        np.abs(mesh.nodes[:, 1]) <= float(half_thickness_m) + 1.0e-15
    )
    if not np.any(mask):
        return math.nan
    return float(np.max(mesh.nodes[mask, 0]))


def _json_segments(segments) -> list[list[list[float]]]:
    output: list[list[list[float]]] = []
    for p0, p1 in list(segments or []):
        output.append([
            [float(p0[0]), float(p0[1])],
            [float(p1[0]), float(p1[1])],
        ])
    return output


def _production_refinement_radius_m(
    width_m: float,
    height_m: float,
    requested_tip_h_m: float,
) -> float:
    h = max(float(requested_tip_h_m), 0.0)
    if h <= 0.0:
        return 0.0
    return min(0.15 * max(float(width_m), float(height_m)), max(40.0 * h, 0.05e-3))


def record_production_j_probe_v10059(
    *,
    path: str | Path,
    mesh,
    u: np.ndarray,
    ep_gp: np.ndarray,
    sigma_gp: np.ndarray,
    psi_tension_gp: np.ndarray,
    d: np.ndarray,
    D: np.ndarray,
    Kmat,
    mat,
    Ftop_N_per_thickness: float,
    Uapp_m: float,
    tip_xy: np.ndarray,
    direction: np.ndarray,
    half_thickness_m: float,
    kill_r_m: float,
    production_ell_m: float,
    production_segments,
    production_exclude_radius_m: float,
    production_path: str,
    contour_outer_m: Iterable[float],
    specimen_width_m: float,
    specimen_height_m: float,
    requested_tip_h_m: float,
    crack_backend_name: str,
    crystal_anisotropic: bool,
    crystal_theta_deg: float,
) -> dict[str, Any]:
    """Persist one post-equilibrium production state and its J ablations."""
    output_path = Path(path).resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    tip = np.asarray(tip_xy, dtype=float).reshape(2)
    growth = np.asarray(direction, dtype=float).reshape(2)
    norm = float(np.linalg.norm(growth))
    growth = growth / norm if norm > 1.0e-30 else np.asarray([1.0, 0.0])
    segments = list(production_segments or [])
    contours = sorted(set(float(value) for value in contour_outer_m))
    if any(value <= 0.0 for value in contours):
        raise ValueError("all parity-audit contour radii must be positive")

    psi_stored_gp, psi_undegraded_gp = elastic_energy_densities(
        mesh, u, ep_gp, sigma_gp, D
    )
    U_density = float(np.sum(psi_stored_gp * mesh.area_e))
    U_matrix = float(0.5 * np.asarray(u) @ (Kmat @ np.asarray(u)))
    energy_closure = abs(U_density - U_matrix) / max(abs(U_matrix), 1.0e-300)

    width = float(specimen_width_m)
    height = float(specimen_height_m)
    sigma_gross = abs(float(Ftop_N_per_thickness)) / max(width, 1.0e-300)
    effective_tip = _effective_killed_tip(mesh, d, half_thickness_m)
    refinement_radius = _production_refinement_radius_m(
        width, height, requested_tip_h_m
    )

    rows: list[dict[str, Any]] = []
    requested_outer = 8.0 * float(production_ell_m)
    all_outer = sorted(set([*contours, requested_outer]))
    for outer in all_outer:
        ell = float(outer) / 8.0
        cfg = JIntegralConfig(
            r_inner_factor=2.0,
            r_outer_factor=8.0,
            q_type="plateau",
        )
        common = dict(
            mesh=mesh,
            u=u,
            sigma_gp=sigma_gp,
            d=d,
            crack_tip=tip,
            crack_direction=growth,
            mat=mat,
            ell=ell,
            cfg=cfg,
            crack_segments=segments,
        )
        J_full, K_full, info_full = compute_J_integral(
            psi_e_gp=psi_stored_gp,
            exclude_radius=float(production_exclude_radius_m),
            **common,
        )
        J_tension, K_tension, info_tension = compute_J_integral(
            psi_e_gp=psi_tension_gp,
            exclude_radius=float(production_exclude_radius_m),
            **common,
        )
        J_full_no_exclusion, K_full_no_exclusion, info_full_no_exclusion = (
            compute_J_integral(
                psi_e_gp=psi_stored_gp,
                exclude_radius=0.0,
                **common,
            )
        )
        J_tension_no_exclusion, K_tension_no_exclusion, _ = compute_J_integral(
            psi_e_gp=psi_tension_gp,
            exclude_radius=0.0,
            **common,
        )
        nearest = min(tip[0], width - tip[0], 0.5 * height)
        within_safety = bool(outer <= 0.80 * nearest)
        row = {
            "point_release": POINT_RELEASE,
            "outer_radius_m": float(outer),
            "outer_radius_um": float(outer) * 1.0e6,
            "is_requested_production_contour": bool(
                math.isclose(float(outer), requested_outer, rel_tol=1.0e-12, abs_tol=1.0e-15)
            ),
            "production_path": str(production_path),
            "production_exclude_radius_m": float(production_exclude_radius_m),
            "production_exclude_radius_um": float(production_exclude_radius_m) * 1.0e6,
            "contour_within_80pct_boundary": within_safety,
            "J_active_elements": int(info_full.get("n_active_elements", 0)),
            "J_active_elements_no_exclusion": int(
                info_full_no_exclusion.get("n_active_elements", 0)
            ),
            "J_full_J_per_m2": float(J_full),
            "J_full_signed_J_per_m2": float(info_full.get("J_signed", math.nan)),
            "J_tension_filtered_J_per_m2": float(J_tension),
            "J_tension_filtered_signed_J_per_m2": float(
                info_tension.get("J_signed", math.nan)
            ),
            "J_full_no_exclusion_J_per_m2": float(J_full_no_exclusion),
            "J_tension_no_exclusion_J_per_m2": float(J_tension_no_exclusion),
            "KJ_full_isotropic_MPa_sqrt_m": float(K_full) / 1.0e6,
            "KJ_tension_isotropic_MPa_sqrt_m": float(K_tension) / 1.0e6,
            "KJ_full_no_exclusion_isotropic_MPa_sqrt_m": (
                float(K_full_no_exclusion) / 1.0e6
            ),
            "KJ_tension_no_exclusion_isotropic_MPa_sqrt_m": (
                float(K_tension_no_exclusion) / 1.0e6
            ),
            "J_full_over_sigma2_m_per_Pa": (
                float(J_full) / sigma_gross**2 if sigma_gross > 0.0 else math.nan
            ),
            "J_tension_over_sigma2_m_per_Pa": (
                float(J_tension) / sigma_gross**2 if sigma_gross > 0.0 else math.nan
            ),
            "J_full_no_exclusion_over_sigma2_m_per_Pa": (
                float(J_full_no_exclusion) / sigma_gross**2
                if sigma_gross > 0.0
                else math.nan
            ),
            "KJ_full_over_sigma_sqrt_m": (
                float(K_full) / sigma_gross if sigma_gross > 0.0 else math.nan
            ),
            "KJ_full_no_exclusion_over_sigma_sqrt_m": (
                float(K_full_no_exclusion) / sigma_gross
                if sigma_gross > 0.0
                else math.nan
            ),
        }
        rows.append(row)

    production_row = min(
        rows,
        key=lambda row: abs(float(row["outer_radius_m"]) - requested_outer),
    )
    payload = {
        "schema": "production_j_probe_v10_0_5_9",
        "point_release": POINT_RELEASE,
        "constitutive_physics_changed": False,
        "plastic_evolution_enabled_in_probe": False,
        "probe_role": "post-equilibrium elastic audit of the actual production initialization path",
        "Uapp_m": float(Uapp_m),
        "Uapp_um": float(Uapp_m) * 1.0e6,
        "Ftop_N_per_thickness": float(Ftop_N_per_thickness),
        "sigma_gross_Pa": sigma_gross,
        "sigma_gross_MPa": sigma_gross / 1.0e6,
        "U_el_density_J_per_m": U_density,
        "U_el_matrix_J_per_m": U_matrix,
        "elastic_energy_closure_relative_error": energy_closure,
        "mesh": {
            "nodes": int(mesh.nn),
            "elements": int(mesh.ne),
            "hbar_tip_m": float(mesh.hbar_tip),
            "hbar_global_m": float(mesh.hbar),
            "requested_tip_h_m": float(requested_tip_h_m),
            "production_refinement_radius_m": refinement_radius,
            "production_refinement_radius_um": refinement_radius * 1.0e6,
        },
        "geometry": {
            "width_m": width,
            "height_m": height,
            "half_thickness_m": float(half_thickness_m),
            "requested_tip_xy_m": [float(tip[0]), float(tip[1])],
            "effective_killed_tip_m": effective_tip,
            "kill_radius_m": float(kill_r_m),
            "crack_backend": str(crack_backend_name),
            "production_segments": _json_segments(segments),
        },
        "elasticity": {
            "crystal_anisotropic": bool(crystal_anisotropic),
            "crystal_theta_deg": float(crystal_theta_deg),
            "K_conversion_note": "K values use isotropic Eprime for display only; J/sigma^2 is primary",
        },
        "production_J_request": {
            "path": str(production_path),
            "ell_m": float(production_ell_m),
            "outer_radius_m": requested_outer,
            "exclude_radius_m": float(production_exclude_radius_m),
            "row": production_row,
        },
        "contours": rows,
    }
    output_path.write_text(json.dumps(payload, indent=2, default=str))
    return payload


def _closest_contour(probe: Mapping[str, Any], outer_radius_m: float) -> dict[str, Any]:
    rows = list(probe.get("contours", []))
    if not rows:
        raise ValueError("production J probe contains no contour rows")
    return dict(
        min(rows, key=lambda row: abs(float(row["outer_radius_m"]) - float(outer_radius_m)))
    )


def analyze_production_j_parity_v10059(
    *,
    reference: Mapping[str, Any],
    probes: Sequence[Mapping[str, Any]],
    selected_outer_radius_m: float = 240.0e-6,
    parity_relative_tolerance: float = 0.10,
    elastic_scaling_relative_tolerance: float = 0.02,
    energy_closure_relative_tolerance: float = 1.0e-6,
) -> dict[str, Any]:
    """Compare production probes with the converged v10.0.5.8 reference."""
    probes = [dict(probe) for probe in probes]
    if len(probes) < 2:
        raise ValueError("at least two production openings are required")
    convergence = dict(reference.get("convergence", {}))
    geometry = dict(reference.get("geometry_factors", {}))
    G_ref = float(convergence["fixed_grip_G_finest_J_per_m2"])
    sigma_ref = float(geometry["sigma_gross_MPa"]) * 1.0e6
    reference_metric = G_ref / max(sigma_ref**2, 1.0e-300)

    case_rows: list[dict[str, Any]] = []
    production_metrics: list[float] = []
    no_exclusion_metrics: list[float] = []
    for probe in sorted(probes, key=lambda item: float(item.get("Uapp_m", 0.0))):
        row = _closest_contour(probe, selected_outer_radius_m)
        production_metric = float(row["J_full_over_sigma2_m_per_Pa"])
        no_exclusion_metric = float(row["J_full_no_exclusion_over_sigma2_m_per_Pa"])
        production_metrics.append(production_metric)
        no_exclusion_metrics.append(no_exclusion_metric)
        tension_difference = _relative_difference(
            float(row["J_full_J_per_m2"]),
            float(row["J_tension_filtered_J_per_m2"]),
        )
        case_rows.append(
            {
                "Uapp_um": float(probe["Uapp_um"]),
                "sigma_gross_MPa": float(probe["sigma_gross_MPa"]),
                "outer_radius_um": float(row["outer_radius_um"]),
                "production_path": row.get("production_path"),
                "production_exclude_radius_um": float(
                    row.get("production_exclude_radius_um", 0.0)
                ),
                "J_full_J_per_m2": float(row["J_full_J_per_m2"]),
                "J_tension_filtered_J_per_m2": float(
                    row["J_tension_filtered_J_per_m2"]
                ),
                "J_full_no_exclusion_J_per_m2": float(
                    row["J_full_no_exclusion_J_per_m2"]
                ),
                "J_full_over_sigma2_m_per_Pa": production_metric,
                "J_full_no_exclusion_over_sigma2_m_per_Pa": no_exclusion_metric,
                "reference_G_over_sigma2_m_per_Pa": reference_metric,
                "production_over_reference": production_metric / reference_metric,
                "no_exclusion_over_reference": no_exclusion_metric / reference_metric,
                "full_vs_tension_relative_difference": tension_difference,
                "energy_closure_relative_error": float(
                    probe["elastic_energy_closure_relative_error"]
                ),
                "production_refinement_radius_um": float(
                    probe["mesh"]["production_refinement_radius_um"]
                ),
                "mesh_hbar_tip_um": float(probe["mesh"]["hbar_tip_m"]) * 1.0e6,
                "effective_killed_tip_um": float(
                    probe["geometry"]["effective_killed_tip_m"]
                ) * 1.0e6,
            }
        )

    production_center = float(np.median(production_metrics))
    no_exclusion_center = float(np.median(no_exclusion_metrics))
    production_scaling_error = max(
        abs(value / production_center - 1.0) for value in production_metrics
    )
    no_exclusion_scaling_error = max(
        abs(value / no_exclusion_center - 1.0) for value in no_exclusion_metrics
    )
    production_ratio = production_center / reference_metric
    no_exclusion_ratio = no_exclusion_center / reference_metric

    closure_passed = all(
        float(row["energy_closure_relative_error"])
        <= float(energy_closure_relative_tolerance)
        for row in case_rows
    )
    elastic_scaling_passed = bool(
        production_scaling_error <= float(elastic_scaling_relative_tolerance)
    )
    production_parity_passed = bool(
        abs(production_ratio - 1.0) <= float(parity_relative_tolerance)
    )
    no_exclusion_parity_passed = bool(
        abs(no_exclusion_ratio - 1.0) <= float(parity_relative_tolerance)
    )

    if not closure_passed:
        status = "production_elastic_energy_closure_failed"
    elif not elastic_scaling_passed:
        status = "production_J_not_quadratic_in_load"
    elif production_parity_passed:
        status = "production_J_parity_passed"
    elif no_exclusion_parity_passed:
        status = "production_exclusion_disk_controls_J_mismatch"
    elif any(
        float(row["production_refinement_radius_um"])
        < float(row["outer_radius_um"])
        for row in case_rows
    ):
        status = "production_refinement_extent_or_mesh_support_mismatch"
    else:
        status = "production_J_fixed_grip_mismatch_unresolved"

    return {
        "schema": "production_j_parity_v10_0_5_9",
        "point_release": POINT_RELEASE,
        "status": status,
        "passed": status == "production_J_parity_passed",
        "constitutive_physics_changed": False,
        "selected_outer_radius_m": float(selected_outer_radius_m),
        "acceptance": {
            "parity_relative_tolerance": float(parity_relative_tolerance),
            "elastic_scaling_relative_tolerance": float(
                elastic_scaling_relative_tolerance
            ),
            "energy_closure_relative_tolerance": float(
                energy_closure_relative_tolerance
            ),
        },
        "reference": {
            "G_fixed_grip_J_per_m2": G_ref,
            "sigma_gross_MPa": sigma_ref / 1.0e6,
            "G_over_sigma2_m_per_Pa": reference_metric,
            "source_schema": reference.get("schema"),
        },
        "production": {
            "median_J_over_sigma2_m_per_Pa": production_center,
            "median_over_reference": production_ratio,
            "maximum_load_scaling_relative_error": production_scaling_error,
            "parity_passed": production_parity_passed,
            "elastic_scaling_passed": elastic_scaling_passed,
        },
        "no_exclusion_ablation": {
            "median_J_over_sigma2_m_per_Pa": no_exclusion_center,
            "median_over_reference": no_exclusion_ratio,
            "maximum_load_scaling_relative_error": no_exclusion_scaling_error,
            "parity_passed": no_exclusion_parity_passed,
        },
        "energy_closure_passed": closure_passed,
        "cases": case_rows,
    }


__all__ = [
    "POINT_RELEASE",
    "PROBE_JSON",
    "SUMMARY_JSON",
    "PROBE_CSV",
    "CONTOUR_CSV",
    "record_production_j_probe_v10059",
    "analyze_production_j_parity_v10059",
]
