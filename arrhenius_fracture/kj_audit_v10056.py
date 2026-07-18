"""v10.0.5.6 remote-stress/KJ audit and first-passage bracket utilities.

The 2-D solver reports force per unit out-of-plane thickness.  Gross nominal
stress is therefore F_top/Lx.  The audit compares the domain-integral KJ with a
single-edge-crack finite-width reference and, independently, checks whether the
circular J domain can close inside the specimen.

No constitutive kinetics are changed here.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Iterable, Mapping, Sequence


POINT_RELEASE = "10.0.5.6"


@dataclass(frozen=True)
class SpecimenGeometryV10056:
    width_m: float = 2.0e-3
    height_m: float = 4.0e-3
    initial_crack_m: float = 0.5e-3

    def validate(self) -> "SpecimenGeometryV10056":
        if self.width_m <= 0.0 or self.height_m <= 0.0:
            raise ValueError("specimen width and height must be positive")
        if not (0.0 < self.initial_crack_m < self.width_m):
            raise ValueError("initial crack must lie strictly inside specimen width")
        return self

    @property
    def a_over_W(self) -> float:
        return self.initial_crack_m / self.width_m

    @property
    def nearest_tip_boundary_m(self) -> float:
        # Tip is at (a0,0) in a body spanning x=[0,Lx], y=[-Ly/2,Ly/2].
        return min(
            self.initial_crack_m,
            self.width_m - self.initial_crack_m,
            0.5 * self.height_m,
        )


def edge_crack_tension_geometry_factor(a_over_W: float) -> float:
    """Single-edge crack finite-width correction for uniform remote tension.

    Polynomial convention:
      K_I = sigma_gross * sqrt(pi*a) * Y(a/W)
    for 0 < a/W <= 0.6.
    """
    x = float(a_over_W)
    if not (0.0 < x <= 0.6):
        raise ValueError("edge-crack polynomial requires 0 < a/W <= 0.6")
    return (
        1.12
        - 0.231 * x
        + 10.55 * x**2
        - 21.72 * x**3
        + 30.39 * x**4
    )


def gross_nominal_stress_Pa(Ftop_N_per_thickness: float, width_m: float) -> float:
    if width_m <= 0.0:
        raise ValueError("width_m must be positive")
    return abs(float(Ftop_N_per_thickness)) / float(width_m)


def net_section_stress_Pa(
    Ftop_N_per_thickness: float,
    geometry: SpecimenGeometryV10056,
) -> float:
    geometry.validate()
    ligament = geometry.width_m - geometry.initial_crack_m
    return abs(float(Ftop_N_per_thickness)) / ligament


def reference_edge_crack_K_Pa_sqrt_m(
    sigma_gross_Pa: float,
    geometry: SpecimenGeometryV10056,
) -> float:
    geometry.validate()
    Y = edge_crack_tension_geometry_factor(geometry.a_over_W)
    return (
        max(float(sigma_gross_Pa), 0.0)
        * math.sqrt(math.pi * geometry.initial_crack_m)
        * Y
    )


def contour_geometry_audit(
    *,
    outer_radius_m: float,
    geometry: SpecimenGeometryV10056,
    safety_fraction: float = 0.80,
) -> dict[str, Any]:
    geometry.validate()
    outer = float(outer_radius_m)
    if outer <= 0.0:
        raise ValueError("outer_radius_m must be positive")
    nearest = geometry.nearest_tip_boundary_m
    safe_limit = float(safety_fraction) * nearest
    return {
        "outer_radius_m": outer,
        "nearest_tip_boundary_m": nearest,
        "safe_outer_radius_limit_m": safe_limit,
        "outer_over_nearest_boundary": outer / nearest,
        "contour_closes_inside_body": bool(outer < nearest),
        "contour_within_safety_limit": bool(outer <= safe_limit),
    }


def build_kj_audit_row(
    *,
    Ftop_N_per_thickness: float,
    KJ_Pa_sqrt_m: float,
    outer_radius_m: float,
    geometry: SpecimenGeometryV10056,
    n_active_elements: int | None = None,
    safety_fraction: float = 0.80,
) -> dict[str, Any]:
    sigma_gross = gross_nominal_stress_Pa(
        Ftop_N_per_thickness, geometry.width_m
    )
    sigma_net = net_section_stress_Pa(Ftop_N_per_thickness, geometry)
    K_ref = reference_edge_crack_K_Pa_sqrt_m(sigma_gross, geometry)
    geom = contour_geometry_audit(
        outer_radius_m=outer_radius_m,
        geometry=geometry,
        safety_fraction=safety_fraction,
    )
    KJ = max(float(KJ_Pa_sqrt_m), 0.0)
    row: dict[str, Any] = {
        **geom,
        "Ftop_N_per_m_thickness": float(Ftop_N_per_thickness),
        "sigma_gross_MPa": sigma_gross / 1.0e6,
        "sigma_net_MPa": sigma_net / 1.0e6,
        "a_over_W": geometry.a_over_W,
        "edge_geometry_factor_Y": edge_crack_tension_geometry_factor(
            geometry.a_over_W
        ),
        "KJ_MPa_sqrt_m": KJ / 1.0e6,
        "K_LEFM_gross_MPa_sqrt_m": K_ref / 1.0e6,
        "KJ_over_K_LEFM_gross": KJ / K_ref if K_ref > 0.0 else math.nan,
        "KJ_per_sigma_gross_sqrt_m": KJ / sigma_gross
        if sigma_gross > 0.0
        else math.nan,
    }
    if n_active_elements is not None:
        row["J_active_elements"] = int(n_active_elements)
    return row


def select_contour_plateau(
    rows: Sequence[Mapping[str, Any]],
    *,
    relative_tolerance: float = 0.10,
    minimum_points: int = 3,
    minimum_active_elements: int = 12,
) -> dict[str, Any]:
    """Select a stable valid KJ/sigma plateau from consecutive contour radii."""
    valid = []
    for source in rows:
        row = dict(source)
        active = int(row.get("J_active_elements", minimum_active_elements))
        slope = float(row.get("KJ_per_sigma_gross_sqrt_m", math.nan))
        if (
            bool(row.get("contour_within_safety_limit", False))
            and active >= minimum_active_elements
            and math.isfinite(slope)
            and slope > 0.0
        ):
            valid.append(row)
    valid.sort(key=lambda row: float(row["outer_radius_m"]))

    best: tuple[int, int, float] | None = None
    n = len(valid)
    for start in range(n):
        for stop in range(start + minimum_points, n + 1):
            window = valid[start:stop]
            slopes = sorted(
                float(row["KJ_per_sigma_gross_sqrt_m"]) for row in window
            )
            median = slopes[len(slopes) // 2]
            spread = max(abs(value / median - 1.0) for value in slopes)
            if spread <= relative_tolerance:
                score = (stop - start, -spread)
                if best is None or score > (best[1] - best[0], -best[2]):
                    best = (start, stop, spread)

    if best is None:
        return {
            "status": "no_valid_plateau",
            "selected_outer_radius_m": None,
            "valid_contour_count": len(valid),
            "minimum_points": int(minimum_points),
            "relative_tolerance": float(relative_tolerance),
        }

    start, stop, spread = best
    window = valid[start:stop]
    selected = window[len(window) // 2]
    ratios = [float(row["KJ_over_K_LEFM_gross"]) for row in window]
    return {
        "status": "plateau_selected",
        "selected_outer_radius_m": float(selected["outer_radius_m"]),
        "selected_rJ_cluster_ell_m": float(selected["outer_radius_m"]) / 8.0,
        "plateau_outer_radii_m": [float(row["outer_radius_m"]) for row in window],
        "plateau_max_relative_spread": float(spread),
        "plateau_KJ_over_K_LEFM_min": min(ratios),
        "plateau_KJ_over_K_LEFM_max": max(ratios),
        "selected_row": selected,
        "valid_contour_count": len(valid),
        "relative_tolerance": float(relative_tolerance),
    }


LEGACY_LIMITER_LABELS = {
    0: "legacy_or_stochastic_unmapped",
    1: "cleavage_clock",
    2: "stored_process_zone",
    3: "emitted_process_zone",
    4: "mobile_process_zone",
    5: "escape_process_zone",
    6: "peierls_transport",
    7: "taylor_transport",
    8: "cycle_horizon",
}


def enrich_stochastic_block_rows(
    rows: Sequence[Mapping[str, Any]],
    scheduler_records: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    """Merge authoritative textual stochastic scheduler records into block rows."""
    if len(rows) != len(scheduler_records):
        raise ValueError(
            "block/scheduler record count mismatch: "
            f"{len(rows)} != {len(scheduler_records)}"
        )
    output: list[dict[str, Any]] = []
    for raw, sched in zip(rows, scheduler_records):
        row = dict(raw)
        code = int(float(row.get("cycle_limiter_code", -1)))
        mode = str(sched.get("mode", "unknown"))
        base_limiter = str(sched.get("base_limiter", "unknown"))
        final_limiter = str(sched.get("final_limiter", sched.get("limiter", "unknown")))
        row.update(
            {
                "stochastic_scheduler_mode": mode,
                "stochastic_event_rate_per_cycle": float(
                    sched.get("event_rate_per_cycle", 0.0)
                ),
                "stochastic_expected_state_events": float(
                    sched.get("expected_state_events", 0.0)
                ),
                "base_cycle_limiter": base_limiter,
                "final_cycle_limiter": final_limiter,
                "cycle_limiter_label": final_limiter
                if final_limiter != "unknown"
                else LEGACY_LIMITER_LABELS.get(code, f"unknown_code_{code}"),
            }
        )
        output.append(row)
    return output


def classify_first_passage_rows(
    rows: Iterable[Mapping[str, Any]],
) -> dict[str, Any]:
    ordered = sorted(
        (dict(row) for row in rows),
        key=lambda row: float(row["delta_sigma_requested_MPa"]),
    )
    no_fire = [row for row in ordered if not bool(row.get("first_passage_observed"))]
    fired = [row for row in ordered if bool(row.get("first_passage_observed"))]
    lower = max(no_fire, key=lambda row: float(row["delta_sigma_requested_MPa"])) if no_fire else None
    upper = min(fired, key=lambda row: float(row["delta_sigma_requested_MPa"])) if fired else None
    bracketed = bool(
        lower is not None
        and upper is not None
        and float(lower["delta_sigma_requested_MPa"])
        < float(upper["delta_sigma_requested_MPa"])
    )
    return {
        "status": "bracketed" if bracketed else "unbracketed",
        "lower_no_first_passage": lower,
        "upper_first_passage": upper,
        "stress_interval_MPa": (
            [
                float(lower["delta_sigma_requested_MPa"]),
                float(upper["delta_sigma_requested_MPa"]),
            ]
            if bracketed
            else None
        ),
        "n_cases": len(ordered),
        "n_first_passage": len(fired),
        "n_no_first_passage": len(no_fire),
    }


__all__ = [
    "POINT_RELEASE",
    "SpecimenGeometryV10056",
    "edge_crack_tension_geometry_factor",
    "gross_nominal_stress_Pa",
    "net_section_stress_Pa",
    "reference_edge_crack_K_Pa_sqrt_m",
    "contour_geometry_audit",
    "build_kj_audit_row",
    "select_contour_plateau",
    "enrich_stochastic_block_rows",
    "classify_first_passage_rows",
]
