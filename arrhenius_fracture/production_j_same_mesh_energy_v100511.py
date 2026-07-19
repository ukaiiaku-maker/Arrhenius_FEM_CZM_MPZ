"""v10.0.5.11 same-production-mesh fixed-grip energy-release audit."""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from .fixed_grip_elastic_audit_v10058 import _solve_elastic_state
from .production_j_refinement_support_v100510 import (
    analyze_refinement_support_v100510,
    record_production_j_refinement_probe_v100510,
)

POINT_RELEASE = "10.0.5.11"
PROBE_JSON = "production_j_same_mesh_probe_v10_0_5_11.json"
SUMMARY_JSON = "production_j_same_mesh_energy_v10_0_5_11.json"
PROBE_CSV = "production_j_same_mesh_cases_v10_0_5_11.csv"
CONTOUR_CSV = "production_j_same_mesh_contours_v10_0_5_11.csv"
RADIAL_CSV = "production_mesh_radial_support_v10_0_5_11.csv"
LAUNCH_FAILURE_JSON = "production_j_same_mesh_launch_failure_v10_0_5_11.json"


def second_order_forward_release_rate_v100511(
    a0_m: float,
    a1_m: float,
    a2_m: float,
    U0_J_per_m: float,
    U1_J_per_m: float,
    U2_J_per_m: float,
) -> dict[str, float]:
    h1 = float(a1_m) - float(a0_m)
    h2 = float(a2_m) - float(a1_m)
    if h1 <= 0.0 or h2 <= 0.0:
        raise ValueError("same-mesh crack stations must increase strictly")
    asymmetry = abs(h1 - h2) / max(0.5 * (h1 + h2), 1.0e-300)
    if asymmetry > 1.0e-10:
        raise ValueError(
            "v10.0.5.11 requires equally spaced first three production corridor stations"
        )
    h = 0.5 * (h1 + h2)
    G_h = (float(U0_J_per_m) - float(U1_J_per_m)) / h
    G_2h = (float(U0_J_per_m) - float(U2_J_per_m)) / (2.0 * h)
    G_second = (
        3.0 * float(U0_J_per_m)
        - 4.0 * float(U1_J_per_m)
        + float(U2_J_per_m)
    ) / (2.0 * h)
    consistency = abs(G_h - G_2h) / max(abs(G_second), 1.0e-300)
    return {
        "station_spacing_m": h,
        "station_spacing_um": h * 1.0e6,
        "G_forward_h_J_per_m2": G_h,
        "G_forward_2h_J_per_m2": G_2h,
        "G_second_order_forward_J_per_m2": G_second,
        "one_step_two_step_relative_difference": consistency,
        "station_spacing_asymmetry": asymmetry,
    }


def _same_mesh_stations(mesh: Any, base_tip_m: float) -> list[float]:
    centers = np.asarray(
        getattr(mesh, "production_refinement_centers_m", []), dtype=float
    )
    if centers.ndim != 2 or centers.shape[1] < 2:
        raise RuntimeError("v10.0.5.11 mesh lacks production corridor centers")
    xs = sorted(
        set(
            float(row[0])
            for row in centers
            if abs(float(row[1])) <= 1.0e-12
            and float(row[0]) >= float(base_tip_m) - 1.0e-15
        )
    )
    if len(xs) < 3:
        raise RuntimeError("v10.0.5.11 requires at least three forward corridor stations")
    if not math.isclose(xs[0], float(base_tip_m), rel_tol=0.0, abs_tol=1.0e-12):
        raise RuntimeError(
            "first production corridor station does not match the effective starter tip"
        )
    return xs[:3]


def record_production_j_same_mesh_probe_v100511(**kwargs) -> dict[str, Any]:
    boundary = kwargs.pop("boundary_data")
    total_opening = float(kwargs.pop("total_grip_opening_m"))
    path = Path(kwargs["path"]).resolve()
    mesh = kwargs["mesh"]
    material = kwargs["mat"]
    D = np.asarray(kwargs["D"], dtype=float)
    half_thickness = float(kwargs["half_thickness_m"])

    payload = record_production_j_refinement_probe_v100510(**kwargs)
    base_tip = float(payload["geometry"]["effective_killed_tip_m"])
    stations = _same_mesh_stations(mesh, base_tip)
    states = [
        _solve_elastic_state(
            mesh,
            boundary,
            material,
            D,
            crack_m=station,
            half_thickness_m=half_thickness,
            total_grip_opening_m=total_opening,
            kappa=1.0e-6,
        )
        for station in stations
    ]
    release = second_order_forward_release_rate_v100511(
        *(state["effective_crack_m"] for state in states),
        *(state["U_el_density_J_per_m"] for state in states),
    )
    sigma_Pa = float(payload["sigma_gross_MPa"]) * 1.0e6
    G = float(release["G_second_order_forward_J_per_m2"])
    release.update(
        {
            "crack_stations_m": [float(state["effective_crack_m"]) for state in states],
            "crack_stations_um": [1.0e6 * float(state["effective_crack_m"]) for state in states],
            "elastic_energies_J_per_m": [float(state["U_el_density_J_per_m"]) for state in states],
            "elastic_energy_matrix_J_per_m": [float(state["U_el_matrix_J_per_m"]) for state in states],
            "energy_closure_relative_errors": [float(state["energy_closure_relative_error"]) for state in states],
            "G_over_sigma2_m_per_Pa": G / max(sigma_Pa * sigma_Pa, 1.0e-300),
            "same_mesh_and_boundary_conditions": True,
            "cohesive_segments_in_energy_perturbation": 0,
        }
    )
    payload.update(
        {
            "schema": "production_j_same_mesh_probe_v10_0_5_11",
            "point_release": POINT_RELEASE,
            "same_mesh_fixed_grip_energy_release": release,
        }
    )
    path.write_text(json.dumps(payload, indent=2, default=str))
    return payload


def analyze_same_mesh_energy_v100511(
    *,
    reference: Mapping[str, Any],
    probes: Sequence[Mapping[str, Any]],
    selected_outer_radius_m: float = 240.0e-6,
    accepted_contours_m: Sequence[float] = (180e-6, 240e-6, 300e-6),
    parity_relative_tolerance: float = 0.10,
    elastic_scaling_relative_tolerance: float = 0.02,
    energy_closure_relative_tolerance: float = 1.0e-6,
    contour_stability_relative_tolerance: float = 0.10,
    derivative_consistency_relative_tolerance: float = 0.10,
) -> dict[str, Any]:
    base = analyze_refinement_support_v100510(
        reference=reference,
        probes=probes,
        selected_outer_radius_m=selected_outer_radius_m,
        accepted_contours_m=accepted_contours_m,
        parity_relative_tolerance=parity_relative_tolerance,
        elastic_scaling_relative_tolerance=elastic_scaling_relative_tolerance,
        energy_closure_relative_tolerance=energy_closure_relative_tolerance,
        contour_stability_relative_tolerance=contour_stability_relative_tolerance,
    )
    accepted = sorted(set(float(v) for v in accepted_contours_m))
    cases = []
    G_metrics = []
    derivative_ok = True
    same_mesh_parity_ok = True
    for probe in sorted(probes, key=lambda p: float(p.get("Uapp_m", 0.0))):
        release = dict(probe["same_mesh_fixed_grip_energy_release"])
        G = float(release["G_second_order_forward_J_per_m2"])
        G_metric = float(release["G_over_sigma2_m_per_Pa"])
        G_metrics.append(G_metric)
        contours = []
        for radius in accepted:
            row = min(
                probe["contours"],
                key=lambda r: abs(float(r["outer_radius_m"]) - radius),
            )
            if not math.isclose(float(row["outer_radius_m"]), radius, rel_tol=1e-10, abs_tol=1e-12):
                raise ValueError(f"probe lacks accepted contour {radius:.16g} m")
            contours.append(float(row["J_full_J_per_m2"]))
        J_median = float(np.median(np.asarray(contours, dtype=float)))
        ratio = J_median / max(G, 1.0e-300)
        derivative_case_ok = (
            float(release["one_step_two_step_relative_difference"])
            <= derivative_consistency_relative_tolerance
        )
        parity_case_ok = abs(ratio - 1.0) <= parity_relative_tolerance
        derivative_ok = bool(derivative_ok and derivative_case_ok)
        same_mesh_parity_ok = bool(same_mesh_parity_ok and parity_case_ok)
        cases.append(
            {
                "Uapp_um": float(probe["Uapp_um"]),
                "G_same_mesh_J_per_m2": G,
                "G_same_mesh_over_sigma2_m_per_Pa": G_metric,
                "J_accepted_contour_median_J_per_m2": J_median,
                "J_median_over_G_same_mesh": ratio,
                "derivative_consistency_relative_error": float(
                    release["one_step_two_step_relative_difference"]
                ),
                "derivative_consistency_passed": derivative_case_ok,
                "same_mesh_J_G_parity_passed": parity_case_ok,
                "crack_stations_um": release["crack_stations_um"],
            }
        )
    G_metrics = np.asarray(G_metrics, dtype=float)
    G_scaling_error = float(
        np.max(np.abs(G_metrics / max(float(np.median(G_metrics)), 1.0e-300) - 1.0))
    )
    G_scaling_ok = G_scaling_error <= elastic_scaling_relative_tolerance
    if not base.get("physical_refinement_support_passed", False):
        status = "production_refinement_support_inadequate"
    elif not base.get("contour_stability_passed", False):
        status = "production_J_contour_instability"
    elif not derivative_ok:
        status = "same_mesh_energy_derivative_not_converged"
    elif not G_scaling_ok:
        status = "same_mesh_G_not_quadratic_in_load"
    elif same_mesh_parity_ok:
        status = "production_J_same_mesh_energy_parity_passed"
    else:
        status = "production_J_same_mesh_energy_mismatch"
    return {
        "schema": "production_j_same_mesh_energy_v10_0_5_11",
        "point_release": POINT_RELEASE,
        "status": status,
        "passed": status == "production_J_same_mesh_energy_parity_passed",
        "constitutive_physics_changed": False,
        "base_v10_0_5_9_analysis": base["base_v10_0_5_9_analysis"],
        "base_v10_0_5_10_analysis": base,
        "same_mesh_cases": cases,
        "same_mesh_derivative_consistency_passed": derivative_ok,
        "same_mesh_G_scaling_relative_error": G_scaling_error,
        "same_mesh_G_scaling_passed": G_scaling_ok,
        "same_mesh_J_G_parity_passed": same_mesh_parity_ok,
        "external_v10_0_5_8_reference_retained_as_diagnostic_only": True,
    }


__all__ = [
    "POINT_RELEASE", "PROBE_JSON", "SUMMARY_JSON", "PROBE_CSV", "CONTOUR_CSV",
    "RADIAL_CSV", "LAUNCH_FAILURE_JSON", "second_order_forward_release_rate_v100511",
    "record_production_j_same_mesh_probe_v100511", "analyze_same_mesh_energy_v100511",
]
