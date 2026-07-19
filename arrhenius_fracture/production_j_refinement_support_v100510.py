"""v10.0.5.10 physical refinement-support diagnostics."""
from __future__ import annotations

import json
import math
import os
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

from .production_j_parity_v10059 import (
    analyze_production_j_parity_v10059,
    record_production_j_probe_v10059,
)

POINT_RELEASE = "10.0.5.10"
PROBE_JSON = "production_j_refinement_probe_v10_0_5_10.json"
SUMMARY_JSON = "production_j_refinement_support_v10_0_5_10.json"
PROBE_CSV = "production_j_refinement_cases_v10_0_5_10.csv"
CONTOUR_CSV = "production_j_refinement_contours_v10_0_5_10.csv"
RADIAL_CSV = "production_mesh_radial_support_v10_0_5_10.csv"


def _values(text: str, scale: float = 1.0) -> list[float]:
    return [float(x) * scale for x in str(text).replace(",", " ").split() if x]


def deduplicate_contours_v100510(
    rows: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    ordered = sorted((dict(r) for r in rows), key=lambda r: float(r["outer_radius_m"]))
    unique: list[dict[str, Any]] = []
    for row in ordered:
        radius = float(row["outer_radius_m"])
        if unique and math.isclose(
            radius,
            float(unique[-1]["outer_radius_m"]),
            rel_tol=1.0e-12,
            abs_tol=1.0e-15,
        ):
            old = unique[-1]
            keep = row if row.get("is_requested_production_contour", False) else old
            keep = dict(keep)
            keep["is_requested_production_contour"] = bool(
                row.get("is_requested_production_contour", False)
                or old.get("is_requested_production_contour", False)
            )
            unique[-1] = keep
        else:
            unique.append(row)
    return unique


def _element_h(mesh) -> np.ndarray:
    x = np.asarray(mesh.nodes)[np.asarray(mesh.elems, dtype=int)]
    return np.mean(
        np.stack(
            [
                np.linalg.norm(x[:, 1] - x[:, 0], axis=1),
                np.linalg.norm(x[:, 2] - x[:, 1], axis=1),
                np.linalg.norm(x[:, 0] - x[:, 2], axis=1),
            ],
            axis=1,
        ),
        axis=1,
    )


def radial_mesh_support_v100510(
    mesh, tip_xy: np.ndarray, radial_edges_m: Iterable[float]
) -> list[dict[str, Any]]:
    edges = sorted(set(float(v) for v in radial_edges_m))
    if len(edges) < 2 or edges[0] < 0.0:
        raise ValueError("radial support edges must be increasing and nonnegative")
    cent = np.asarray(mesh.nodes)[np.asarray(mesh.elems, dtype=int)].mean(axis=1)
    radius = np.linalg.norm(cent - np.asarray(tip_xy).reshape(1, 2), axis=1)
    h = _element_h(mesh)
    configured = float(getattr(mesh, "production_refinement_radius_m", math.nan))
    out: list[dict[str, Any]] = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        if hi <= lo:
            raise ValueError("radial support edges must be strictly increasing")
        values = h[(radius >= lo) & (radius < hi)]
        p95 = float(np.quantile(values, 0.95)) if values.size else math.nan
        out.append(
            {
                "point_release": POINT_RELEASE,
                "r_inner_um": lo * 1.0e6,
                "r_outer_um": hi * 1.0e6,
                "element_count": int(values.size),
                "inside_configured_refinement": bool(
                    np.isfinite(configured) and hi <= configured + 1.0e-15
                ),
                "h_mean_um": float(np.mean(values)) * 1.0e6 if values.size else math.nan,
                "h_median_um": float(np.median(values)) * 1.0e6 if values.size else math.nan,
                "h_p95_um": p95 * 1.0e6 if values.size else math.nan,
                "h_max_um": float(np.max(values)) * 1.0e6 if values.size else math.nan,
                "annulus_width_over_h_p95": (hi - lo) / max(p95, 1.0e-300)
                if values.size
                else math.nan,
            }
        )
    return out


def record_production_j_refinement_probe_v100510(**kwargs) -> dict[str, Any]:
    path = Path(kwargs["path"]).resolve()
    mesh = kwargs["mesh"]
    tip = np.asarray(kwargs["tip_xy"], dtype=float).reshape(2)
    payload = record_production_j_probe_v10059(**kwargs)
    radius = float(getattr(mesh, "production_refinement_radius_m", math.nan))
    policy = getattr(mesh, "production_refinement_policy", None)
    if not np.isfinite(radius) or radius <= 0.0 or policy is None:
        raise RuntimeError("v10.0.5.10 did not receive its fixed physical refinement mesh")

    radial = _values(
        os.environ.get("ARRHENIUS_V100510_RADIAL_EDGES_UM", "0 60 100 140 180 240 300 330"),
        1.0e-6,
    )
    radial.extend([0.0, radius])
    original_count = len(payload.get("contours", []))
    contours = deduplicate_contours_v100510(payload.get("contours", []))
    requested = float(payload["production_J_request"]["outer_radius_m"])
    payload.update(
        {
            "schema": "production_j_refinement_probe_v10_0_5_10",
            "point_release": POINT_RELEASE,
            "base_probe_schema": "production_j_probe_v10_0_5_9",
            "contours": contours,
            "contour_duplicate_count_removed": original_count - len(contours),
            "radial_mesh_support": radial_mesh_support_v100510(mesh, tip, radial),
        }
    )
    payload["mesh"].update(
        {
            "production_refinement_radius_m": radius,
            "production_refinement_radius_um": radius * 1.0e6,
            "production_refinement_policy": str(policy),
            "production_refinement_centers_m": getattr(mesh, "production_refinement_centers_m", None),
            "production_refinement_h_fine_m": float(
                getattr(mesh, "production_refinement_h_fine_m", math.nan)
            ),
            "production_refinement_tip_ratio": float(
                getattr(mesh, "production_refinement_tip_ratio", math.nan)
            ),
        }
    )
    payload["production_J_request"]["row"] = min(
        contours, key=lambda row: abs(float(row["outer_radius_m"]) - requested)
    )
    path.write_text(json.dumps(payload, indent=2, default=str))
    return payload


def _closest(probe: Mapping[str, Any], radius_m: float) -> dict[str, Any]:
    rows = list(probe.get("contours", []))
    row = min(rows, key=lambda r: abs(float(r["outer_radius_m"]) - radius_m))
    if not math.isclose(float(row["outer_radius_m"]), radius_m, rel_tol=1e-10, abs_tol=1e-12):
        raise ValueError(f"probe lacks requested contour {radius_m:.16g} m")
    return dict(row)


def analyze_refinement_support_v100510(
    *,
    reference: Mapping[str, Any],
    probes: Sequence[Mapping[str, Any]],
    selected_outer_radius_m: float = 240.0e-6,
    accepted_contours_m: Sequence[float] = (180e-6, 240e-6, 300e-6),
    parity_relative_tolerance: float = 0.10,
    elastic_scaling_relative_tolerance: float = 0.02,
    energy_closure_relative_tolerance: float = 1.0e-6,
    contour_stability_relative_tolerance: float = 0.10,
) -> dict[str, Any]:
    accepted = sorted(set(float(v) for v in accepted_contours_m))
    base = analyze_production_j_parity_v10059(
        reference=reference,
        probes=probes,
        selected_outer_radius_m=selected_outer_radius_m,
        parity_relative_tolerance=parity_relative_tolerance,
        elastic_scaling_relative_tolerance=elastic_scaling_relative_tolerance,
        energy_closure_relative_tolerance=energy_closure_relative_tolerance,
    )
    cases = []
    max_spread = 0.0
    support = True
    for probe in sorted(probes, key=lambda p: float(p.get("Uapp_m", 0.0))):
        metrics = np.asarray(
            [float(_closest(probe, r)["J_full_over_sigma2_m_per_Pa"]) for r in accepted]
        )
        center = float(np.median(metrics))
        spread = float((np.max(metrics) - np.min(metrics)) / max(abs(center), 1e-300))
        radius = float(probe["mesh"]["production_refinement_radius_m"])
        case_support = radius >= max(accepted) - 1.0e-15
        support = bool(support and case_support)
        max_spread = max(max_spread, spread)
        cases.append(
            {
                "Uapp_um": float(probe["Uapp_um"]),
                "production_refinement_radius_um": radius * 1.0e6,
                "accepted_contours_um": [r * 1.0e6 for r in accepted],
                "J_over_sigma2_m_per_Pa": metrics.tolist(),
                "peak_to_peak_over_median": spread,
                "physical_refinement_support_passed": case_support,
            }
        )
    stable = max_spread <= contour_stability_relative_tolerance
    if not base.get("energy_closure_passed", False):
        status = "production_elastic_energy_closure_failed"
    elif not base["production"].get("elastic_scaling_passed", False):
        status = "production_J_not_quadratic_in_load"
    elif not support:
        status = "production_refinement_support_inadequate"
    elif not stable:
        status = "production_J_contour_instability"
    elif base["production"].get("parity_passed", False):
        status = "production_refinement_support_parity_passed"
    else:
        status = "production_J_parity_failed_with_adequate_support"
    return {
        "schema": "production_j_refinement_support_v10_0_5_10",
        "point_release": POINT_RELEASE,
        "status": status,
        "passed": status == "production_refinement_support_parity_passed",
        "constitutive_physics_changed": False,
        "mesh_policy_changed_for_audit_only": True,
        "base_v10_0_5_9_analysis": base,
        "accepted_contours_m": accepted,
        "contour_stability_relative_tolerance": contour_stability_relative_tolerance,
        "maximum_contour_peak_to_peak_over_median": max_spread,
        "contour_stability_passed": stable,
        "physical_refinement_support_passed": support,
        "contour_cases": cases,
    }


__all__ = [
    "POINT_RELEASE", "PROBE_JSON", "SUMMARY_JSON", "PROBE_CSV", "CONTOUR_CSV",
    "RADIAL_CSV", "deduplicate_contours_v100510", "radial_mesh_support_v100510",
    "record_production_j_refinement_probe_v100510", "analyze_refinement_support_v100510",
]
