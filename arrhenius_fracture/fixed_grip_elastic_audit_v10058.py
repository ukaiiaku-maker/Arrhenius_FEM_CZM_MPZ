"""Fixed-grip elastic FEM convergence audit for the v10.0.5.x geometry.

This module supplies the mechanics-derived reference that the earlier KJ audit was
missing.  It deliberately bypasses Arrhenius hazards, plasticity evolution and
cohesive propagation.  For the actual finite killed-notch geometry it compares:

* the fixed-displacement energy-release rate obtained from nearby crack lengths;
* the domain J integral evaluated with the full stored elastic energy density;
* the legacy domain J integral evaluated with the tensile-filtered energy density;
* the sharp finite-width edge-crack LEFM factor, retained only as a secondary
  geometry comparison.

The primary validation is J_full / G_fixed_grip, not KJ / K_LEFM.  This avoids
mistaking the finite starter notch and fixed-grip compliance for a sharp-crack
infinite-loading calibration error.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from .config import ElasticProperties, GeometryConfig, JIntegralConfig, MeshConfig
from .crystal import W_C11, W_C12, W_C44, cubic_plane_strain_D, zener_ratio
from .fem import (
    assemble_mechanics,
    elastic_energy_densities,
    plane_strain_D,
    solve_dirichlet,
)
from .j_integral import compute_J_integral
from .kj_audit_v10056 import edge_crack_tension_geometry_factor
from .mesh import make_boundary_data, make_tri_mesh

POINT_RELEASE = "10.0.5.8"
ENERGY_CSV = "fixed_grip_energy_release_v10_0_5_8.csv"
J_CSV = "fixed_grip_J_contours_v10_0_5_8.csv"
SUMMARY_JSON = "fixed_grip_elastic_convergence_v10_0_5_8.json"
FIGURE_PNG = "fixed_grip_elastic_convergence_v10_0_5_8.png"


@dataclass(frozen=True)
class FixedGripAuditConfig:
    width_m: float = 2.0e-3
    height_m: float = 4.0e-3
    crack_m: float = 0.5e-3
    notch_half_thickness_m: float = 0.08e-3
    total_grip_opening_m: float = 2.0e-6
    nx: int = 36
    ny: int = 72
    tip_ratio: float = 1.20
    seed: int = 42
    kappa: float = 1.0e-6
    anisotropic: bool = True
    crystal_theta_deg: float = 45.0
    C11_Pa: float = W_C11
    C12_Pa: float = W_C12
    C44_Pa: float = W_C44
    mesh_relative_tolerance: float = 0.10
    crack_increment_relative_tolerance: float = 0.10
    contour_relative_tolerance: float = 0.10
    minimum_plateau_points: int = 3
    minimum_J_active_elements: int = 12
    J_over_G_min: float = 0.80
    J_over_G_max: float = 1.20
    energy_closure_relative_tolerance: float = 1.0e-6

    def validate(self) -> "FixedGripAuditConfig":
        if self.width_m <= 0.0 or self.height_m <= 0.0:
            raise ValueError("specimen dimensions must be positive")
        if not (0.0 < self.crack_m < self.width_m):
            raise ValueError("crack_m must lie inside the specimen")
        if self.notch_half_thickness_m <= 0.0:
            raise ValueError("notch_half_thickness_m must be positive")
        if self.total_grip_opening_m <= 0.0:
            raise ValueError("total_grip_opening_m must be positive")
        if self.nx < 4 or self.ny < 8:
            raise ValueError("mesh background resolution is too small")
        if not (0.0 < self.kappa < 1.0):
            raise ValueError("kappa must be between zero and one")
        if self.minimum_plateau_points < 2:
            raise ValueError("minimum_plateau_points must be at least two")
        if not (0.0 < self.J_over_G_min <= self.J_over_G_max):
            raise ValueError("invalid J/G acceptance interval")
        return self


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
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


def _finite(value: Any) -> bool:
    try:
        return bool(np.isfinite(float(value)))
    except (TypeError, ValueError):
        return False


def _relative_difference(a: float, b: float) -> float:
    scale = max(abs(float(a)), abs(float(b)), 1.0e-300)
    return abs(float(a) - float(b)) / scale


def _damage_for_crack(mesh, crack_m: float, half_thickness_m: float) -> np.ndarray:
    x = mesh.nodes[:, 0]
    y = mesh.nodes[:, 1]
    d = np.zeros(mesh.nn, dtype=float)
    d[(x <= float(crack_m)) & (np.abs(y) <= float(half_thickness_m))] = 1.0
    return d


def _effective_killed_tip(mesh, d: np.ndarray, half_thickness_m: float) -> float:
    mask = (np.asarray(d) > 0.5) & (
        np.abs(mesh.nodes[:, 1]) <= float(half_thickness_m) + 1.0e-15
    )
    if not np.any(mask):
        return float("nan")
    return float(np.max(mesh.nodes[mask, 0]))


def three_point_derivative_at_center(
    x_minus: float,
    x_zero: float,
    x_plus: float,
    f_minus: float,
    f_zero: float,
    f_plus: float,
) -> float:
    """Derivative at x_zero for unequal left/right spacings.

    The expression is the derivative of the quadratic interpolant through the
    three supplied points.  It reduces to the centered difference for equal
    spacings.
    """
    h_left = float(x_zero) - float(x_minus)
    h_right = float(x_plus) - float(x_zero)
    if not (_finite(h_left) and _finite(h_right)) or h_left <= 0.0 or h_right <= 0.0:
        raise ValueError("effective crack lengths must increase strictly")
    c_minus = -h_right / (h_left * (h_left + h_right))
    c_zero = (h_right - h_left) / (h_left * h_right)
    c_plus = h_left / (h_right * (h_left + h_right))
    return c_minus * float(f_minus) + c_zero * float(f_zero) + c_plus * float(f_plus)


def fixed_grip_release_rate(
    *,
    a_minus_m: float,
    a_zero_m: float,
    a_plus_m: float,
    U_minus_J_per_m: float,
    U_zero_J_per_m: float,
    U_plus_J_per_m: float,
) -> dict[str, float]:
    """Compute G=-dU_el/da at fixed imposed displacement."""
    dU_da = three_point_derivative_at_center(
        a_minus_m,
        a_zero_m,
        a_plus_m,
        U_minus_J_per_m,
        U_zero_J_per_m,
        U_plus_J_per_m,
    )
    h_left = float(a_zero_m) - float(a_minus_m)
    h_right = float(a_plus_m) - float(a_zero_m)
    G_backward = -(
        float(U_zero_J_per_m) - float(U_minus_J_per_m)
    ) / h_left
    G_forward = -(
        float(U_plus_J_per_m) - float(U_zero_J_per_m)
    ) / h_right
    G_center = -float(dU_da)
    one_sided_relative_spread = _relative_difference(G_backward, G_forward)
    step_asymmetry = abs(h_left - h_right) / max(0.5 * (h_left + h_right), 1.0e-300)
    return {
        "G_fixed_grip_J_per_m2": G_center,
        "G_backward_J_per_m2": G_backward,
        "G_forward_J_per_m2": G_forward,
        "one_sided_G_relative_spread": one_sided_relative_spread,
        "effective_left_increment_m": h_left,
        "effective_right_increment_m": h_right,
        "effective_increment_asymmetry": step_asymmetry,
    }


def select_contiguous_plateau(
    rows: Iterable[dict[str, Any]],
    *,
    value_key: str,
    relative_tolerance: float,
    minimum_points: int,
    minimum_active_elements: int = 0,
) -> dict[str, Any]:
    """Select the widest contiguous contour interval satisfying a spread bound."""
    valid = [
        dict(row)
        for row in rows
        if _finite(row.get("outer_radius_m"))
        and _finite(row.get(value_key))
        and float(row[value_key]) > 0.0
        and int(float(row.get("J_active_elements", 0) or 0)) >= int(minimum_active_elements)
        and bool(row.get("contour_within_safety_limit", True))
    ]
    valid.sort(key=lambda row: float(row["outer_radius_m"]))
    candidates: list[tuple[int, float, int, int]] = []
    for i in range(len(valid)):
        for j in range(i + int(minimum_points), len(valid) + 1):
            values = np.asarray([float(row[value_key]) for row in valid[i:j]], dtype=float)
            center = float(np.median(values))
            spread = float(np.max(np.abs(values / center - 1.0))) if center > 0.0 else math.inf
            if spread <= float(relative_tolerance):
                candidates.append((j - i, spread, i, j))
    if not candidates:
        return {
            "status": "no_valid_plateau",
            "value_key": value_key,
            "valid_contour_count": len(valid),
            "minimum_points": int(minimum_points),
            "relative_tolerance": float(relative_tolerance),
        }
    n, spread, i, j = sorted(candidates, key=lambda item: (-item[0], item[1], item[2]))[0]
    selected = valid[i:j]
    values = np.asarray([float(row[value_key]) for row in selected], dtype=float)
    center_index = len(selected) // 2
    return {
        "status": "plateau_selected",
        "value_key": value_key,
        "n_points": int(n),
        "maximum_relative_spread": float(spread),
        "median_value": float(np.median(values)),
        "minimum_value": float(np.min(values)),
        "maximum_value": float(np.max(values)),
        "outer_radii_m": [float(row["outer_radius_m"]) for row in selected],
        "selected_outer_radius_m": float(selected[center_index]["outer_radius_m"]),
        "selected_row": dict(selected[center_index]),
    }


def _elasticity(config: FixedGripAuditConfig, material: ElasticProperties) -> np.ndarray:
    if config.anisotropic:
        return cubic_plane_strain_D(
            config.C11_Pa,
            config.C12_Pa,
            config.C44_Pa,
            config.crystal_theta_deg,
        )
    return plane_strain_D(material)


def _solve_elastic_state(
    mesh,
    boundary,
    material: ElasticProperties,
    D: np.ndarray,
    *,
    crack_m: float,
    half_thickness_m: float,
    total_grip_opening_m: float,
    kappa: float,
) -> dict[str, Any]:
    d = _damage_for_crack(mesh, crack_m, half_thickness_m)
    u0 = np.zeros(mesh.ndof, dtype=float)
    ep_gp = np.zeros((3, mesh.ne), dtype=float)
    rho_gp = np.zeros(mesh.ne, dtype=float)
    Kmat, Rint, _, _, _, _ = assemble_mechanics(
        mesh,
        u0,
        ep_gp,
        rho_gp,
        d,
        D,
        material,
        kappa=float(kappa),
    )
    u, Ftop = solve_dirichlet(
        Kmat,
        Rint,
        u0,
        boundary,
        0.5 * float(total_grip_opening_m),
        -0.5 * float(total_grip_opening_m),
    )
    _, _, sigma_gp, sigma_eq_gp, sigma1_gp, psi_tension = assemble_mechanics(
        mesh,
        u,
        ep_gp,
        rho_gp,
        d,
        D,
        material,
        kappa=float(kappa),
    )
    psi_stored, psi_undegraded = elastic_energy_densities(
        mesh,
        u,
        ep_gp,
        sigma_gp,
        D,
    )
    U_density = float(np.sum(psi_stored * mesh.area_e))
    U_matrix = float(0.5 * u @ (Kmat @ u))
    closure = abs(U_density - U_matrix) / max(abs(U_matrix), 1.0e-300)
    return {
        "requested_crack_m": float(crack_m),
        "effective_crack_m": _effective_killed_tip(mesh, d, half_thickness_m),
        "d": d,
        "u": u,
        "sigma_gp": sigma_gp,
        "sigma_eq_gp": sigma_eq_gp,
        "sigma1_gp": sigma1_gp,
        "psi_tension_gp": psi_tension,
        "psi_stored_gp": psi_stored,
        "psi_undegraded_gp": psi_undegraded,
        "Ftop_N_per_thickness": float(Ftop),
        "U_el_density_J_per_m": U_density,
        "U_el_matrix_J_per_m": U_matrix,
        "energy_closure_relative_error": closure,
    }


def _contour_geometry(
    outer_radius_m: float,
    *,
    width_m: float,
    height_m: float,
    tip_x_m: float,
    safety_fraction: float = 0.80,
) -> dict[str, Any]:
    nearest = min(
        float(tip_x_m),
        float(width_m) - float(tip_x_m),
        0.5 * float(height_m),
    )
    safe = float(safety_fraction) * nearest
    return {
        "outer_radius_m": float(outer_radius_m),
        "nearest_tip_boundary_m": nearest,
        "safe_outer_radius_limit_m": safe,
        "contour_within_safety_limit": bool(float(outer_radius_m) <= safe),
    }


def _J_rows_for_state(
    mesh,
    state: dict[str, Any],
    material: ElasticProperties,
    *,
    contour_outer_m: Iterable[float],
    width_m: float,
    height_m: float,
    half_thickness_m: float,
    tip_h_fine_m: float,
) -> list[dict[str, Any]]:
    tip_x = float(state["effective_crack_m"])
    tip = np.asarray([tip_x, 0.0], dtype=float)
    direction = np.asarray([1.0, 0.0], dtype=float)
    segments = [(np.asarray([0.0, 0.0]), np.asarray([tip_x, 0.0]))]
    rows: list[dict[str, Any]] = []
    for outer in sorted(set(float(value) for value in contour_outer_m)):
        geom = _contour_geometry(
            outer,
            width_m=width_m,
            height_m=height_m,
            tip_x_m=tip_x,
        )
        if not geom["contour_within_safety_limit"]:
            rows.append(
                {
                    **geom,
                    "tip_h_fine_m": float(tip_h_fine_m),
                    "audit_status": "rejected_boundary_intersection",
                }
            )
            continue
        cfg = JIntegralConfig(r_inner_factor=2.0, r_outer_factor=8.0, q_type="plateau")
        ell = outer / 8.0
        J_full, K_full, info_full = compute_J_integral(
            mesh,
            state["u"],
            state["sigma_gp"],
            state["psi_stored_gp"],
            state["d"],
            tip,
            direction,
            material,
            ell,
            cfg=cfg,
            crack_segments=segments,
            exclude_radius=0.0,
        )
        J_tension, K_tension, info_tension = compute_J_integral(
            mesh,
            state["u"],
            state["sigma_gp"],
            state["psi_tension_gp"],
            state["d"],
            tip,
            direction,
            material,
            ell,
            cfg=cfg,
            crack_segments=segments,
            exclude_radius=0.0,
        )
        rows.append(
            {
                **geom,
                "tip_h_fine_m": float(tip_h_fine_m),
                "annulus_width_m": 0.75 * outer,
                "annulus_radial_elements_estimate": 0.75 * outer / max(float(tip_h_fine_m), 1.0e-300),
                "J_active_elements": int(info_full.get("n_active_elements", 0)),
                "J_full_J_per_m2": float(J_full),
                "J_full_signed_J_per_m2": float(info_full.get("J_signed", math.nan)),
                "KJ_full_isotropic_MPa_sqrt_m": float(K_full) / 1.0e6,
                "J_tension_filtered_J_per_m2": float(J_tension),
                "J_tension_filtered_signed_J_per_m2": float(info_tension.get("J_signed", math.nan)),
                "KJ_tension_filtered_isotropic_MPa_sqrt_m": float(K_tension) / 1.0e6,
                "audit_status": "complete",
            }
        )
    return rows


def _render_figure(out: Path, energy_rows: list[dict[str, Any]], j_rows: list[dict[str, Any]]) -> str | None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return None

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)
    ax = axes[0]
    deltas = sorted(set(float(row["requested_increment_um"]) for row in energy_rows))
    for delta in deltas:
        selected = sorted(
            [row for row in energy_rows if float(row["requested_increment_um"]) == delta],
            key=lambda row: float(row["tip_h_fine_um"]),
        )
        ax.plot(
            [float(row["tip_h_fine_um"]) for row in selected],
            [float(row["G_fixed_grip_J_per_m2"]) for row in selected],
            marker="o",
            label=f"da={delta:g} um",
        )
    ax.set_xlabel("tip mesh spacing (um)")
    ax.set_ylabel("fixed-grip G (J/m2)")
    ax.invert_xaxis()
    ax.legend()

    ax = axes[1]
    hs = sorted(set(float(row["tip_h_fine_um"]) for row in j_rows))
    for h in hs:
        selected = sorted(
            [
                row
                for row in j_rows
                if float(row["tip_h_fine_um"]) == h
                and _finite(row.get("J_full_over_G_reference"))
            ],
            key=lambda row: float(row["outer_radius_um"]),
        )
        ax.plot(
            [float(row["outer_radius_um"]) for row in selected],
            [float(row["J_full_over_G_reference"]) for row in selected],
            marker="o",
            label=f"h={h:g} um",
        )
    ax.axhline(1.0, linestyle="--")
    ax.set_xlabel("J outer radius (um)")
    ax.set_ylabel("J_full / G_fixed-grip")
    ax.legend()

    path = out / FIGURE_PNG
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return str(path)


def run_audit(
    config: FixedGripAuditConfig,
    *,
    tip_h_fine_m: Iterable[float],
    crack_increment_m: Iterable[float],
    contour_outer_m: Iterable[float],
    out: str | Path,
) -> dict[str, Any]:
    config = config.validate()
    hs = sorted(set(float(value) for value in tip_h_fine_m), reverse=True)
    increments = sorted(set(float(value) for value in crack_increment_m), reverse=True)
    contours = sorted(set(float(value) for value in contour_outer_m))
    if len(hs) < 2:
        raise ValueError("at least two tip mesh spacings are required")
    if len(increments) < 2:
        raise ValueError("at least two crack increments are required")
    if len(contours) < config.minimum_plateau_points:
        raise ValueError("not enough contour radii for the requested plateau")
    if any(value <= 0.0 for value in hs + increments + contours):
        raise ValueError("mesh spacings, crack increments and contours must be positive")
    if config.crack_m - max(increments) <= 0.0:
        raise ValueError("largest crack decrement leaves the specimen")
    if config.crack_m + max(increments) >= config.width_m:
        raise ValueError("largest crack increment leaves the ligament")

    root = Path(out).resolve()
    root.mkdir(parents=True, exist_ok=True)
    material = ElasticProperties()
    D = _elasticity(config, material)
    geom = GeometryConfig(
        Lx=config.width_m,
        Ly=config.height_m,
        a0=config.crack_m,
        notch_half_thickness=config.notch_half_thickness_m,
    )

    energy_rows: list[dict[str, Any]] = []
    all_j_rows: list[dict[str, Any]] = []
    state_zero_by_h: dict[float, dict[str, Any]] = {}

    for h in hs:
        mesh_cfg = MeshConfig(
            nx=config.nx,
            ny=config.ny,
            jitter=0.0,
            tip_h_fine=h,
            tip_ratio=config.tip_ratio,
        )
        mesh = make_tri_mesh(geom, mesh_cfg, seed=config.seed, tip_center=(config.crack_m, 0.0))
        boundary = make_boundary_data(mesh, geom)
        state_cache: dict[float, dict[str, Any]] = {}

        def state(crack: float) -> dict[str, Any]:
            key = float(crack)
            if key not in state_cache:
                state_cache[key] = _solve_elastic_state(
                    mesh,
                    boundary,
                    material,
                    D,
                    crack_m=key,
                    half_thickness_m=config.notch_half_thickness_m,
                    total_grip_opening_m=config.total_grip_opening_m,
                    kappa=config.kappa,
                )
            return state_cache[key]

        zero = state(config.crack_m)
        state_zero_by_h[h] = zero
        for requested_da in increments:
            minus = state(config.crack_m - requested_da)
            plus = state(config.crack_m + requested_da)
            release = fixed_grip_release_rate(
                a_minus_m=minus["effective_crack_m"],
                a_zero_m=zero["effective_crack_m"],
                a_plus_m=plus["effective_crack_m"],
                U_minus_J_per_m=minus["U_el_density_J_per_m"],
                U_zero_J_per_m=zero["U_el_density_J_per_m"],
                U_plus_J_per_m=plus["U_el_density_J_per_m"],
            )
            G = float(release["G_fixed_grip_J_per_m2"])
            Kiso = math.sqrt(max(G, 0.0) * material.Eprime)
            sigma = abs(float(zero["Ftop_N_per_thickness"])) / config.width_m
            denom = sigma * math.sqrt(math.pi * max(float(zero["effective_crack_m"]), 1.0e-300))
            energy_rows.append(
                {
                    "point_release": POINT_RELEASE,
                    "tip_h_fine_m": h,
                    "tip_h_fine_um": h * 1.0e6,
                    "mesh_hbar_tip_m": float(mesh.hbar_tip),
                    "mesh_hbar_global_m": float(mesh.hbar),
                    "mesh_nodes": int(mesh.nn),
                    "mesh_elements": int(mesh.ne),
                    "requested_increment_m": requested_da,
                    "requested_increment_um": requested_da * 1.0e6,
                    "requested_crack_minus_m": config.crack_m - requested_da,
                    "requested_crack_zero_m": config.crack_m,
                    "requested_crack_plus_m": config.crack_m + requested_da,
                    "effective_crack_minus_m": minus["effective_crack_m"],
                    "effective_crack_zero_m": zero["effective_crack_m"],
                    "effective_crack_plus_m": plus["effective_crack_m"],
                    "U_minus_J_per_m": minus["U_el_density_J_per_m"],
                    "U_zero_J_per_m": zero["U_el_density_J_per_m"],
                    "U_plus_J_per_m": plus["U_el_density_J_per_m"],
                    "maximum_energy_closure_relative_error": max(
                        float(minus["energy_closure_relative_error"]),
                        float(zero["energy_closure_relative_error"]),
                        float(plus["energy_closure_relative_error"]),
                    ),
                    **release,
                    "K_fixed_grip_isotropic_MPa_sqrt_m": Kiso / 1.0e6,
                    "sigma_gross_zero_MPa": sigma / 1.0e6,
                    "Y_fixed_grip_energy": Kiso / max(denom, 1.0e-300),
                    "G_positive": bool(G > 0.0),
                }
            )

        all_j_rows.extend(
            _J_rows_for_state(
                mesh,
                zero,
                material,
                contour_outer_m=contours,
                width_m=config.width_m,
                height_m=config.height_m,
                half_thickness_m=config.notch_half_thickness_m,
                tip_h_fine_m=h,
            )
        )

    # Reference G for each mesh is the smallest requested crack increment.
    smallest_da = min(increments)
    G_by_h: dict[float, float] = {}
    for h in hs:
        row = next(
            row
            for row in energy_rows
            if math.isclose(float(row["tip_h_fine_m"]), h)
            and math.isclose(float(row["requested_increment_m"]), smallest_da)
        )
        G_by_h[h] = float(row["G_fixed_grip_J_per_m2"])
    for row in all_j_rows:
        h = float(row["tip_h_fine_m"])
        row["tip_h_fine_um"] = h * 1.0e6
        row["outer_radius_um"] = float(row["outer_radius_m"]) * 1.0e6
        Gref = G_by_h.get(h, math.nan)
        row["G_reference_J_per_m2"] = Gref
        if _finite(row.get("J_full_J_per_m2")) and Gref > 0.0:
            row["J_full_over_G_reference"] = float(row["J_full_J_per_m2"]) / Gref
            row["J_tension_filtered_over_G_reference"] = (
                float(row["J_tension_filtered_J_per_m2"]) / Gref
            )
        else:
            row["J_full_over_G_reference"] = math.nan
            row["J_tension_filtered_over_G_reference"] = math.nan

    finest_h = min(hs)
    second_finest_h = sorted(hs)[1]
    two_smallest_da = sorted(increments)[:2]
    finest_small = next(
        row
        for row in energy_rows
        if math.isclose(float(row["tip_h_fine_m"]), finest_h)
        and math.isclose(float(row["requested_increment_m"]), two_smallest_da[0])
    )
    finest_next_da = next(
        row
        for row in energy_rows
        if math.isclose(float(row["tip_h_fine_m"]), finest_h)
        and math.isclose(float(row["requested_increment_m"]), two_smallest_da[1])
    )
    second_finest_small = next(
        row
        for row in energy_rows
        if math.isclose(float(row["tip_h_fine_m"]), second_finest_h)
        and math.isclose(float(row["requested_increment_m"]), two_smallest_da[0])
    )
    mesh_error = _relative_difference(
        float(finest_small["G_fixed_grip_J_per_m2"]),
        float(second_finest_small["G_fixed_grip_J_per_m2"]),
    )
    increment_error = _relative_difference(
        float(finest_small["G_fixed_grip_J_per_m2"]),
        float(finest_next_da["G_fixed_grip_J_per_m2"]),
    )
    energy_closure_error = max(
        float(row["maximum_energy_closure_relative_error"])
        for row in energy_rows
        if math.isclose(float(row["tip_h_fine_m"]), finest_h)
    )

    finest_J_rows = [
        row for row in all_j_rows if math.isclose(float(row["tip_h_fine_m"]), finest_h)
    ]
    plateau = select_contiguous_plateau(
        finest_J_rows,
        value_key="J_full_J_per_m2",
        relative_tolerance=config.contour_relative_tolerance,
        minimum_points=config.minimum_plateau_points,
        minimum_active_elements=config.minimum_J_active_elements,
    )
    G_finest = float(finest_small["G_fixed_grip_J_per_m2"])
    if plateau.get("status") == "plateau_selected" and G_finest > 0.0:
        J_over_G = float(plateau["median_value"]) / G_finest
    else:
        J_over_G = math.nan

    positive_G = bool(
        all(float(row["G_fixed_grip_J_per_m2"]) > 0.0 for row in energy_rows)
    )
    energy_closure_passed = bool(
        energy_closure_error <= config.energy_closure_relative_tolerance
    )
    crack_increment_converged = bool(
        increment_error <= config.crack_increment_relative_tolerance
    )
    mesh_converged = bool(mesh_error <= config.mesh_relative_tolerance)
    J_plateau_passed = plateau.get("status") == "plateau_selected"
    J_energy_agreement_passed = bool(
        _finite(J_over_G)
        and config.J_over_G_min <= J_over_G <= config.J_over_G_max
    )

    if not positive_G:
        status = "nonpositive_fixed_grip_energy_release"
    elif not energy_closure_passed:
        status = "elastic_energy_closure_failed"
    elif not crack_increment_converged:
        status = "crack_increment_not_converged"
    elif not mesh_converged:
        status = "mesh_not_converged"
    elif not J_plateau_passed:
        status = "full_energy_J_plateau_not_found"
    elif not J_energy_agreement_passed:
        status = "full_energy_J_fixed_grip_mismatch"
    else:
        status = "fixed_grip_elastic_audit_passed"

    sigma_finest = float(finest_small["sigma_gross_zero_MPa"]) * 1.0e6
    K_energy_finest = float(finest_small["K_fixed_grip_isotropic_MPa_sqrt_m"]) * 1.0e6
    a_eff_finest = float(finest_small["effective_crack_zero_m"])
    Y_energy = K_energy_finest / max(
        sigma_finest * math.sqrt(math.pi * a_eff_finest), 1.0e-300
    )
    Y_edge = edge_crack_tension_geometry_factor(a_eff_finest / config.width_m)

    _write_csv(root / ENERGY_CSV, energy_rows)
    _write_csv(root / J_CSV, all_j_rows)
    figure = _render_figure(root, energy_rows, all_j_rows)
    payload = {
        "schema": "fixed_grip_elastic_convergence_v10_0_5_8",
        "point_release": POINT_RELEASE,
        "status": status,
        "passed": status == "fixed_grip_elastic_audit_passed",
        "constitutive_physics_changed": False,
        "benchmark_role": (
            "geometry-specific elastic reference for the actual finite killed notch "
            "under symmetric fixed-grip displacement"
        ),
        "config": asdict(config),
        "mesh_spacings_m": hs,
        "crack_increments_m": increments,
        "contour_outer_radii_m": contours,
        "elasticity": {
            "anisotropic": bool(config.anisotropic),
            "crystal_theta_deg": float(config.crystal_theta_deg),
            "C11_Pa": float(config.C11_Pa),
            "C12_Pa": float(config.C12_Pa),
            "C44_Pa": float(config.C44_Pa),
            "zener_ratio": float(zener_ratio(config.C11_Pa, config.C12_Pa, config.C44_Pa)),
            "K_conversion_note": (
                "K values use isotropic Eprime only for common-unit display; the primary "
                "anisotropy-safe comparison is J_full/G_fixed_grip"
            ),
        },
        "convergence": {
            "finest_tip_h_m": finest_h,
            "second_finest_tip_h_m": second_finest_h,
            "smallest_crack_increment_m": two_smallest_da[0],
            "second_smallest_crack_increment_m": two_smallest_da[1],
            "fixed_grip_G_finest_J_per_m2": G_finest,
            "mesh_relative_difference": mesh_error,
            "mesh_relative_tolerance": config.mesh_relative_tolerance,
            "mesh_converged": mesh_converged,
            "crack_increment_relative_difference": increment_error,
            "crack_increment_relative_tolerance": config.crack_increment_relative_tolerance,
            "crack_increment_converged": crack_increment_converged,
            "elastic_energy_closure_relative_error": energy_closure_error,
            "elastic_energy_closure_relative_tolerance": config.energy_closure_relative_tolerance,
            "elastic_energy_closure_passed": energy_closure_passed,
            "positive_energy_release_all_cases": positive_G,
        },
        "full_energy_J_plateau": plateau,
        "J_full_over_G_fixed_grip": J_over_G,
        "J_over_G_acceptance": [config.J_over_G_min, config.J_over_G_max],
        "J_energy_agreement_passed": J_energy_agreement_passed,
        "geometry_factors": {
            "Y_fixed_grip_energy_actual_notch": Y_energy,
            "Y_sharp_edge_crack_reference": Y_edge,
            "Y_fixed_grip_over_sharp_edge": Y_energy / max(Y_edge, 1.0e-300),
            "effective_crack_m": a_eff_finest,
            "sigma_gross_MPa": sigma_finest / 1.0e6,
            "K_fixed_grip_isotropic_MPa_sqrt_m": K_energy_finest / 1.0e6,
        },
        "outputs": {
            "energy_release_csv": str(root / ENERGY_CSV),
            "J_contours_csv": str(root / J_CSV),
            "figure": figure,
        },
    }
    (root / SUMMARY_JSON).write_text(json.dumps(payload, indent=2, default=str))
    print(f"FIXED-GRIP ELASTIC AUDIT STATUS: {status}")
    print(root / SUMMARY_JSON)
    return payload


def _values(text: str, scale: float = 1.0) -> list[float]:
    return [float(token) * scale for token in str(text).replace(",", " ").split() if token]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("runs/v10_0_5_8_fixed_grip_elastic_convergence_v1"))
    parser.add_argument("--width-mm", type=float, default=2.0)
    parser.add_argument("--height-mm", type=float, default=4.0)
    parser.add_argument("--crack-mm", type=float, default=0.5)
    parser.add_argument("--notch-half-thickness-um", type=float, default=80.0)
    parser.add_argument("--grip-opening-um", type=float, default=2.0)
    parser.add_argument("--tip-h-um", default="10 5 2.5")
    parser.add_argument("--crack-increment-um", default="20 10 5")
    parser.add_argument("--contour-outer-um", default="100 140 180 240 300")
    parser.add_argument("--nx", type=int, default=36)
    parser.add_argument("--ny", type=int, default=72)
    parser.add_argument("--tip-ratio", type=float, default=1.20)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--kappa", type=float, default=1.0e-6)
    parser.add_argument("--anisotropic", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--crystal-theta-deg", type=float, default=45.0)
    parser.add_argument("--C11-GPa", type=float, default=W_C11 / 1.0e9)
    parser.add_argument("--C12-GPa", type=float, default=W_C12 / 1.0e9)
    parser.add_argument("--C44-GPa", type=float, default=W_C44 / 1.0e9)
    parser.add_argument("--mesh-rel-tol", type=float, default=0.10)
    parser.add_argument("--increment-rel-tol", type=float, default=0.10)
    parser.add_argument("--contour-rel-tol", type=float, default=0.10)
    parser.add_argument("--minimum-plateau-points", type=int, default=3)
    parser.add_argument("--minimum-J-active-elements", type=int, default=12)
    parser.add_argument("--J-over-G-min", type=float, default=0.80)
    parser.add_argument("--J-over-G-max", type=float, default=1.20)
    parser.add_argument("--energy-closure-rel-tol", type=float, default=1.0e-6)
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    args = build_parser().parse_args(list(argv) if argv is not None else None)
    config = FixedGripAuditConfig(
        width_m=args.width_mm * 1.0e-3,
        height_m=args.height_mm * 1.0e-3,
        crack_m=args.crack_mm * 1.0e-3,
        notch_half_thickness_m=args.notch_half_thickness_um * 1.0e-6,
        total_grip_opening_m=args.grip_opening_um * 1.0e-6,
        nx=args.nx,
        ny=args.ny,
        tip_ratio=args.tip_ratio,
        seed=args.seed,
        kappa=args.kappa,
        anisotropic=args.anisotropic,
        crystal_theta_deg=args.crystal_theta_deg,
        C11_Pa=args.C11_GPa * 1.0e9,
        C12_Pa=args.C12_GPa * 1.0e9,
        C44_Pa=args.C44_GPa * 1.0e9,
        mesh_relative_tolerance=args.mesh_rel_tol,
        crack_increment_relative_tolerance=args.increment_rel_tol,
        contour_relative_tolerance=args.contour_rel_tol,
        minimum_plateau_points=args.minimum_plateau_points,
        minimum_J_active_elements=args.minimum_J_active_elements,
        J_over_G_min=args.J_over_G_min,
        J_over_G_max=args.J_over_G_max,
        energy_closure_relative_tolerance=args.energy_closure_rel_tol,
    )
    payload = run_audit(
        config,
        tip_h_fine_m=_values(args.tip_h_um, 1.0e-6),
        crack_increment_m=_values(args.crack_increment_um, 1.0e-6),
        contour_outer_m=_values(args.contour_outer_um, 1.0e-6),
        out=args.out,
    )
    return 0 if bool(payload.get("passed")) else 2


__all__ = [
    "POINT_RELEASE",
    "FixedGripAuditConfig",
    "three_point_derivative_at_center",
    "fixed_grip_release_rate",
    "select_contiguous_plateau",
    "run_audit",
    "build_parser",
    "main",
]
