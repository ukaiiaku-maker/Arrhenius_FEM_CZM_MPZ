"""Zero-dimensional persistent-site screening model for v9.13.

This is a mean-field reduction of the validated one-dimensional persistent-site
model.  It deliberately reuses the same candidate objects, EXP-floor surfaces,
persistent source multiplicity, dynamic front width, dynamic tip radius,
implicit backstress-limited emission root, Peierls/Taylor surfaces, and cleavage
multi-hit clock.  It removes only the ahead-of-tip spatial coordinate.

The reduction is a screening fidelity.  Candidate promotion still requires the
full one-dimensional autonomous R-curve calculation.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Any, Sequence

import numpy as np

from .emergent_gnd_rcurve_v913 import RCurveLoadingMap
from .emergent_gnd_state_v913 import (
    effective_front_width_m,
    persistent_site_multiplicity,
    solve_backstress_limited_activations,
)
from .emergent_gnd_types_v912 import EV_TO_J, KB_J_PER_K, ExpFloorSurface
from .emergent_gnd_types_v913 import CandidateParameters, CommonPhysics


@dataclass(frozen=True)
class ZeroDReductionGeometry:
    dx_m: float
    state_strip_width_m: float
    source_bins: int
    tip_weight_norm: float
    source_weight_mean: float
    cell_area_m2: float
    active_arc_factor: float
    density_increment_per_activation_m2: tuple[float, ...]
    slip_count_increment_per_activation: tuple[float, ...]


@dataclass
class ZeroDState:
    mobile_m2: np.ndarray
    retained_m2: np.ndarray
    local_slip_count_by_system: np.ndarray
    extension_m: float = 0.0
    time_s: float = 0.0
    max_backstress_Pa: float = 0.0
    max_tip_radius_m: float = 0.0
    min_front_width_m: float = float("inf")
    max_multiplicity: float = 0.0
    cumulative_activations: np.ndarray = field(default_factory=lambda: np.zeros(0))


@dataclass(frozen=True)
class ZeroDRunSettings:
    target_projected_extension_m: float = 50.0e-6
    load_increment_factor: float = 2.0
    maximum_applied_displacement_m: float = 2.0e-4
    maximum_load_steps: int = 200_000
    translation_profile: str = "source_zone_tip_weighted"

    def validate(self) -> None:
        if self.target_projected_extension_m <= 0.0:
            raise ValueError("target_projected_extension_m must be positive")
        if self.load_increment_factor <= 0.0:
            raise ValueError("load_increment_factor must be positive")
        if self.maximum_applied_displacement_m <= 0.0:
            raise ValueError("maximum_applied_displacement_m must be positive")
        if self.maximum_load_steps < 1:
            raise ValueError("maximum_load_steps must be positive")
        if self.translation_profile != "source_zone_tip_weighted":
            raise ValueError("unsupported zero-D translation profile")


@dataclass
class ZeroDEvent:
    event_index: int
    threshold_action: float
    K_MPa_sqrt_m: float
    applied_displacement_m: float
    cumulative_projected_extension_m: float
    path_advance_m: float
    projected_advance_m: float
    tip_radius_m: float
    front_width_m: float
    backstress_mean_Pa: float
    multiplicity_per_system: float


@dataclass
class ZeroDResult:
    candidate_id: str
    temperature_K: float
    status: str
    events: list[ZeroDEvent]
    max_backstress_Pa: float
    max_tip_radius_m: float
    min_front_width_m: float
    max_multiplicity: float
    numerical_contract: dict[str, Any]

    @property
    def achieved_projected_extension_m(self) -> float:
        return (
            float(self.events[-1].cumulative_projected_extension_m)
            if self.events
            else 0.0
        )

    def checkpoint_K(self, extension_m: float) -> float:
        if not self.events:
            return float("nan")
        target = max(float(extension_m), 0.0)
        x = np.asarray(
            [event.cumulative_projected_extension_m for event in self.events],
            dtype=float,
        )
        index = int(np.searchsorted(x, target, side="left"))
        if index >= len(self.events):
            return float("nan")
        return float(self.events[index].K_MPa_sqrt_m)

    def as_dict(self) -> dict[str, Any]:
        return {
            "candidate_id": self.candidate_id,
            "temperature_K": self.temperature_K,
            "status": self.status,
            "n_events": len(self.events),
            "achieved_projected_extension_um": (
                self.achieved_projected_extension_m * 1.0e6
            ),
            "K_first_MPa_sqrt_m": (
                self.events[0].K_MPa_sqrt_m if self.events else float("nan")
            ),
            "K_25um_MPa_sqrt_m": self.checkpoint_K(25.0e-6),
            "K_50um_MPa_sqrt_m": self.checkpoint_K(50.0e-6),
            "K_75um_MPa_sqrt_m": self.checkpoint_K(75.0e-6),
            "K_100um_MPa_sqrt_m": self.checkpoint_K(100.0e-6),
            "max_backstress_GPa": self.max_backstress_Pa * 1.0e-9,
            "max_tip_radius_um": self.max_tip_radius_m * 1.0e6,
            "min_front_width_um": self.min_front_width_m * 1.0e6,
            "max_source_multiplicity": self.max_multiplicity,
            "events": [vars(event) for event in self.events],
            "numerical_contract": dict(self.numerical_contract),
        }


def reduction_geometry(physics: CommonPhysics) -> ZeroDReductionGeometry:
    physics.validate()
    dx = float(physics.mpz_length_m) / int(physics.n_bins)
    x = (np.arange(physics.n_bins, dtype=float) + 0.5) * dx
    strip = max(
        float(physics.blunting_length_m),
        dx,
        abs(float(physics.b_m)),
        1.0e-30,
    )
    weights = np.exp(-x / strip)
    norm = max(float(np.sum(weights)), 1.0e-30)
    nsrc = max(
        min(
            int(math.ceil(max(float(physics.source_zone_length_m), dx) / dx)),
            int(physics.n_bins),
        ),
        1,
    )
    source_weight_mean = float(np.sum(weights[:nsrc])) / float(nsrc)
    cell_area = dx * strip
    conversion = np.asarray(
        physics.activation_to_line_content_per_system,
        dtype=float,
    )
    density = conversion * source_weight_mean / max(norm * cell_area, 1.0e-30)
    slip = conversion * source_weight_mean
    active_arc = float(physics.reference_source_area_m2) / max(
        float(physics.r0_m) * float(physics.reference_front_width_m),
        1.0e-30,
    )
    return ZeroDReductionGeometry(
        dx_m=dx,
        state_strip_width_m=strip,
        source_bins=nsrc,
        tip_weight_norm=norm,
        source_weight_mean=source_weight_mean,
        cell_area_m2=cell_area,
        active_arc_factor=active_arc,
        density_increment_per_activation_m2=tuple(float(v) for v in density),
        slip_count_increment_per_activation=tuple(float(v) for v in slip),
    )


def _arrhenius_rate(
    surface: ExpFloorSurface,
    stress_Pa: np.ndarray | float,
    temperature_K: float,
    nu0_s: float,
) -> np.ndarray:
    barrier = surface.barrier_eV(stress_Pa, temperature_K)
    exponent = -np.asarray(barrier, dtype=float) * EV_TO_J / (
        KB_J_PER_K * float(temperature_K)
    )
    return float(nu0_s) * np.exp(np.clip(exponent, -700.0, 0.0))


def _cleavage_effective_rate(
    candidate: CandidateParameters,
    physics: CommonPhysics,
    stress_Pa: float,
    temperature_K: float,
) -> float:
    raw = float(
        _arrhenius_rate(
            candidate.cleavage,
            max(float(stress_Pa), 0.0),
            temperature_K,
            physics.cleavage_nu0_s,
        )
    )
    if physics.cleavage_hits <= 1.0 + 1.0e-12:
        return raw
    from scipy.special import gammainc

    exposure = min(raw * physics.cleavage_correlation_time_s, 1.0e12)
    return float(
        gammainc(physics.cleavage_hits, exposure)
        / physics.cleavage_correlation_time_s
    )


def _emission_factors(physics: CommonPhysics, extension_m: float) -> np.ndarray:
    breakpoints = np.asarray(physics.emission_geometry_extension_m, dtype=float)
    if breakpoints.size == 0:
        return np.asarray(physics.emission_schmid_factors, dtype=float)
    factors = np.asarray(physics.emission_geometry_factors, dtype=float)
    index = int(np.searchsorted(breakpoints, max(extension_m, 0.0), side="right") - 1)
    index = min(max(index, 0), breakpoints.size - 1)
    return factors[index].copy()


def _state_geometry(
    candidate: CandidateParameters,
    physics: CommonPhysics,
    reduced: ZeroDReductionGeometry,
    state: ZeroDState,
) -> dict[str, Any]:
    rho = np.maximum(state.mobile_m2, 0.0) + np.maximum(state.retained_m2, 0.0)
    rho_width = max(
        float(physics.rho_forest_floor_m2) + float(np.sum(rho)),
        float(physics.reference_density_m2),
    )
    minimum = max(float(physics.minimum_front_width_m), abs(float(physics.b_m)))
    maximum = (
        max(float(physics.maximum_front_width_m), minimum)
        if float(physics.maximum_front_width_m) > 0.0
        else max(float(physics.reference_front_width_m), minimum)
    )
    width = effective_front_width_m(
        rho_width,
        reference_width_m=physics.reference_front_width_m,
        reference_density_m2=physics.reference_density_m2,
        minimum_width_m=minimum,
        maximum_width_m=maximum,
    )
    slip_count = float(np.sum(np.maximum(state.local_slip_count_by_system, 0.0)))
    radius = max(
        float(physics.r0_m)
        + max(float(candidate.c_blunt), 0.0)
        * abs(float(physics.b_m))
        * slip_count
        * max(float(physics.blunting_slip_fraction), 0.0),
        float(physics.r0_m),
    )
    multiplicity = persistent_site_multiplicity(
        candidate.rho_source0_m2,
        radius,
        width,
        reduced.active_arc_factor,
    )
    resolved = max(abs(float(candidate.taylor.stress_fraction)), 1.0e-6)
    sigma_back = (
        max(float(physics.persistent_backstress_scale), 0.0)
        * float(physics.G_Pa)
        * abs(float(physics.b_m))
        * np.sqrt(np.maximum(rho, 0.0))
        / resolved
    )
    return {
        "rho_by_system_m2": rho,
        "rho_width_m2": rho_width,
        "front_width_m": width,
        "tip_radius_m": radius,
        "multiplicity_per_system": multiplicity,
        "sigma_back_by_system_Pa": sigma_back,
    }


def _translation_retention_factor(
    distance_m: float,
    physics: CommonPhysics,
    reduced: ZeroDReductionGeometry,
) -> float:
    distance = max(float(distance_m), 0.0)
    if distance <= 0.0:
        return 1.0
    dx = reduced.dx_m
    n = int(physics.n_bins)
    source = np.zeros(n, dtype=float)
    source[: reduced.source_bins] = 1.0
    shifted = np.zeros_like(source)
    length = float(physics.mpz_length_m)
    for i in range(n):
        left = i * dx - distance
        right = (i + 1) * dx - distance
        if right <= 0.0 or left >= length:
            continue
        inside_left = max(left, 0.0)
        inside_right = min(right, length)
        if inside_right <= inside_left:
            continue
        j0 = max(int(math.floor(inside_left / dx)), 0)
        j1 = min(
            int(math.floor((inside_right - 1.0e-15 * dx) / dx)),
            n - 1,
        )
        for j in range(j0, j1 + 1):
            overlap_left = max(inside_left, j * dx)
            overlap_right = min(inside_right, (j + 1) * dx)
            shifted[j] += source[i] * max(overlap_right - overlap_left, 0.0) / dx
    x = (np.arange(n, dtype=float) + 0.5) * dx
    weights = np.exp(-x / reduced.state_strip_width_m)
    before = float(np.sum(source * weights))
    after = float(np.sum(shifted * weights))
    return min(max(after / max(before, 1.0e-30), 0.0), 1.0)


def _advance_state(
    candidate: CandidateParameters,
    physics: CommonPhysics,
    reduced: ZeroDReductionGeometry,
    state: ZeroDState,
    *,
    duration_s: float,
    K_MPa_sqrt_m: float,
    temperature_K: float,
) -> float:
    dt = max(float(duration_s), 0.0)
    if dt <= 0.0:
        geometry = _state_geometry(candidate, physics, reduced, state)
        sigma = max(float(K_MPa_sqrt_m), 0.0) * 1.0e6 / math.sqrt(
            2.0 * math.pi * max(float(geometry["tip_radius_m"]), physics.b_m)
        )
        return _cleavage_effective_rate(candidate, physics, sigma, temperature_K)

    geometry = _state_geometry(candidate, physics, reduced, state)
    radius = float(geometry["tip_radius_m"])
    sigma_applied = max(float(K_MPa_sqrt_m), 0.0) * 1.0e6 / math.sqrt(
        2.0 * math.pi * max(radius, physics.b_m)
    )
    factors = np.abs(_emission_factors(physics, state.extension_m))
    drive = factors * sigma_applied
    rho0 = np.asarray(geometry["rho_by_system_m2"], dtype=float)
    sigma_back0 = np.asarray(geometry["sigma_back_by_system_Pa"], dtype=float)
    multiplicity = float(geometry["multiplicity_per_system"])
    rho_per = np.asarray(reduced.density_increment_per_activation_m2, dtype=float)
    slip_per = np.asarray(reduced.slip_count_increment_per_activation, dtype=float)
    resolved = max(abs(float(candidate.taylor.stress_fraction)), 1.0e-6)
    backstress_prefactor = (
        max(float(physics.persistent_backstress_scale), 0.0)
        * float(physics.G_Pa)
        * abs(float(physics.b_m))
        / resolved
    )

    activations = np.zeros(int(physics.n_systems), dtype=float)
    for system in range(int(physics.n_systems)):
        if drive[system] <= sigma_back0[system]:
            continue

        def rate_at(stress: float) -> float:
            if stress <= 0.0:
                return 0.0
            return float(
                _arrhenius_rate(
                    candidate.emission,
                    stress,
                    temperature_K,
                    physics.emission_nu0_s,
                )
            )

        activations[system] = solve_backstress_limited_activations(
            multiplicity=multiplicity,
            dt_s=dt,
            drive_stress_Pa=float(drive[system]),
            rho_initial_m2=float(rho0[system]),
            rho_increment_per_activation_m2=float(rho_per[system]),
            backstress_prefactor_Pa_sqrt_m2=backstress_prefactor,
            rate_function=rate_at,
            tolerance=physics.implicit_tolerance,
            max_iterations=physics.implicit_max_iterations,
        )

    state.mobile_m2 += activations * rho_per
    state.local_slip_count_by_system += activations * slip_per
    state.cumulative_activations += activations

    total_density = np.maximum(state.mobile_m2 + state.retained_m2, 0.0)
    forest = max(
        float(physics.rho_forest_floor_m2) + float(np.sum(total_density)),
        1.0,
    )
    spacing = 1.0 / (2.0 * math.sqrt(forest))
    jump = float(physics.jump_fraction_of_forest_spacing) * spacing
    external = factors * sigma_applied

    p_surface = candidate.peierls.surface(candidate.emission)
    p_stress = candidate.peierls.stress_fraction * external
    p_rate = np.asarray(
        _arrhenius_rate(
            p_surface,
            np.maximum(p_stress, 0.0),
            temperature_K,
            candidate.peierls.nu0_s,
        ),
        dtype=float,
    )
    peierls_velocity = jump * p_rate
    mfp = float(physics.mean_free_path_coefficient) / math.sqrt(forest)
    encounter = (
        max(float(physics.encounter_efficiency), 0.0)
        * np.abs(peierls_velocity)
        / max(mfp, 1.0e-30)
    )

    t_surface = candidate.taylor.surface(candidate.emission)
    phi = spacing / max(float(physics.b_m), 1.0e-30)
    if math.isfinite(float(physics.taylor_phi_max)):
        phi = min(phi, float(physics.taylor_phi_max))
    t_stress = candidate.taylor.stress_fraction * external * phi
    t_single = np.asarray(
        _arrhenius_rate(
            t_surface,
            np.maximum(t_stress, 0.0),
            temperature_K,
            candidate.taylor.nu0_s,
        ),
        dtype=float,
    )
    corr_length = candidate.taylor_corr_scale / (
        2.0 * math.sqrt(max(candidate.taylor_corr_rho_c_m2, 1.0e-300))
    )
    order = 1.0 + 2.0 * corr_length * math.sqrt(forest)
    taylor = t_single / max(order, 1.0)

    total = state.mobile_m2 + state.retained_m2
    exchange = encounter + taylor
    active = exchange > 0.0
    if np.any(active):
        equilibrium = np.zeros_like(total)
        equilibrium[active] = encounter[active] / exchange[active] * total[active]
        decay = np.ones_like(total)
        decay[active] = np.exp(np.clip(-exchange[active] * dt, -700.0, 0.0))
        state.retained_m2[active] = (
            equilibrium[active]
            + (state.retained_m2[active] - equilibrium[active]) * decay[active]
        )
        state.retained_m2 = np.minimum(np.maximum(state.retained_m2, 0.0), total)
        state.mobile_m2 = np.maximum(total - state.retained_m2, 0.0)

    state.time_s += dt
    geometry = _state_geometry(candidate, physics, reduced, state)
    state.max_backstress_Pa = max(
        state.max_backstress_Pa,
        float(np.mean(geometry["sigma_back_by_system_Pa"])),
    )
    state.max_tip_radius_m = max(state.max_tip_radius_m, float(geometry["tip_radius_m"]))
    state.min_front_width_m = min(state.min_front_width_m, float(geometry["front_width_m"]))
    state.max_multiplicity = max(
        state.max_multiplicity,
        float(geometry["multiplicity_per_system"]),
    )
    sigma_c = max(float(K_MPa_sqrt_m), 0.0) * 1.0e6 / math.sqrt(
        2.0 * math.pi * max(float(geometry["tip_radius_m"]), physics.b_m)
    )
    return _cleavage_effective_rate(candidate, physics, sigma_c, temperature_K)


def run_zero_d_rcurve(
    candidate: CandidateParameters,
    physics: CommonPhysics,
    loading_map: RCurveLoadingMap,
    temperature_K: float,
    *,
    settings: ZeroDRunSettings | None = None,
) -> ZeroDResult:
    settings = settings or ZeroDRunSettings()
    settings.validate()
    loading_map.validate()
    physics.validate()
    reduced = reduction_geometry(physics)
    nsys = int(physics.n_systems)
    state = ZeroDState(
        mobile_m2=np.zeros(nsys, dtype=float),
        retained_m2=np.zeros(nsys, dtype=float),
        local_slip_count_by_system=np.zeros(nsys, dtype=float),
        cumulative_activations=np.zeros(nsys, dtype=float),
        max_tip_radius_m=float(physics.r0_m),
    )
    events: list[ZeroDEvent] = []
    displacement = 0.0
    cumulative_projected = 0.0
    load_steps = 0
    dU = float(loading_map.nominal_dU_m) * float(settings.load_increment_factor)
    dt = float(loading_map.nominal_dt_s) * float(settings.load_increment_factor)

    for index in range(loading_map.n_events):
        if cumulative_projected + 1.0e-15 >= settings.target_projected_extension_m:
            break
        threshold = float(loading_map.threshold_actions[index])
        action = 0.0
        geometry_factor = float(loading_map.K_per_U_MPa_sqrt_m_per_m[index])
        while action < threshold:
            load_steps += 1
            if load_steps > settings.maximum_load_steps:
                status = "right_censored_maximum_load_steps"
                return ZeroDResult(
                    candidate.candidate_id,
                    float(temperature_K),
                    status,
                    events,
                    state.max_backstress_Pa,
                    state.max_tip_radius_m,
                    state.min_front_width_m,
                    state.max_multiplicity,
                    _numerical_contract(settings, reduced, load_steps),
                )
            if displacement >= settings.maximum_applied_displacement_m:
                status = "right_censored_maximum_displacement"
                return ZeroDResult(
                    candidate.candidate_id,
                    float(temperature_K),
                    status,
                    events,
                    state.max_backstress_Pa,
                    state.max_tip_radius_m,
                    state.min_front_width_m,
                    state.max_multiplicity,
                    _numerical_contract(settings, reduced, load_steps),
                )
            K0 = geometry_factor * displacement
            K1 = geometry_factor * (displacement + dU)
            Kmid = 0.5 * (K0 + K1)
            rate = _advance_state(
                candidate,
                physics,
                reduced,
                state,
                duration_s=dt,
                K_MPa_sqrt_m=Kmid,
                temperature_K=float(temperature_K),
            )
            action += max(rate, 0.0) * dt
            displacement += dU

        path_advance = float(loading_map.path_advances_m[index])
        projected = float(loading_map.projected_advances_m[index])
        cumulative_projected += projected
        geometry = _state_geometry(candidate, physics, reduced, state)
        events.append(
            ZeroDEvent(
                event_index=index,
                threshold_action=threshold,
                K_MPa_sqrt_m=geometry_factor * displacement,
                applied_displacement_m=displacement,
                cumulative_projected_extension_m=cumulative_projected,
                path_advance_m=path_advance,
                projected_advance_m=projected,
                tip_radius_m=float(geometry["tip_radius_m"]),
                front_width_m=float(geometry["front_width_m"]),
                backstress_mean_Pa=float(
                    np.mean(geometry["sigma_back_by_system_Pa"])
                ),
                multiplicity_per_system=float(geometry["multiplicity_per_system"]),
            )
        )
        keep = _translation_retention_factor(path_advance, physics, reduced)
        state.mobile_m2 *= keep
        state.retained_m2 *= keep
        state.local_slip_count_by_system *= keep
        state.extension_m += path_advance

    status = (
        "complete"
        if cumulative_projected + 1.0e-15 >= settings.target_projected_extension_m
        else "right_censored_loading_map_exhausted"
    )
    return ZeroDResult(
        candidate.candidate_id,
        float(temperature_K),
        status,
        events,
        state.max_backstress_Pa,
        state.max_tip_radius_m,
        state.min_front_width_m,
        state.max_multiplicity,
        _numerical_contract(settings, reduced, load_steps),
    )


def _numerical_contract(
    settings: ZeroDRunSettings,
    reduced: ZeroDReductionGeometry,
    load_steps: int,
) -> dict[str, Any]:
    return {
        "schema": "v9.13_persistent_zero_d_reduction_v1",
        "fidelity": "zero_dimensional_mean_field_screen",
        "candidate_contract": "v9.13_active_fields_only",
        "finite_source_inventory": False,
        "source_depletion_on_emission": False,
        "source_refresh_on_crack_advance": False,
        "explicit_recovery": False,
        "persistent_multiplicity": True,
        "dynamic_tip_radius": True,
        "dynamic_front_width": True,
        "implicit_backstress_limited_emission": True,
        "spatial_state": False,
        "translation_profile": settings.translation_profile,
        "load_increment_factor": settings.load_increment_factor,
        "load_steps": int(load_steps),
        "reduction_geometry": vars(reduced),
    }


def local_peak_metrics(
    temperatures_K: Sequence[float],
    values: Sequence[float],
    *,
    desired_min_K: float = 850.0,
    desired_max_K: float = 1100.0,
) -> dict[str, float | bool]:
    T = np.asarray(temperatures_K, dtype=float)
    y = np.asarray(values, dtype=float)
    finite = np.isfinite(T) & np.isfinite(y)
    if np.sum(finite) < 3:
        return {
            "peak_temperature_K": float("nan"),
            "peak_value": float("nan"),
            "two_sided_prominence": float("nan"),
            "post_peak_drop": float("nan"),
            "high_temperature_rebound": float("nan"),
            "peak_internal": False,
            "peak_in_desired_window": False,
        }
    order = np.argsort(T[finite])
    T = T[finite][order]
    y = y[finite][order]
    internal = np.arange(1, len(T) - 1)
    local = internal[(y[internal] > y[internal - 1]) & (y[internal] > y[internal + 1])]
    if local.size:
        best = int(local[np.argmax(y[local])])
        peak_internal = True
    else:
        best = int(np.argmax(y))
        peak_internal = 0 < best < len(T) - 1
    left = float(y[best] - np.max(y[:best])) if best > 0 else float("-inf")
    right_neighbor = (
        float(y[best] - y[best + 1]) if best + 1 < len(y) else float("-inf")
    )
    prominence = min(left, right_neighbor)
    if best + 1 < len(y):
        post_min = float(np.min(y[best + 1 :]))
        drop = float(y[best] - post_min)
        rebound = float(np.max(y[best + 1 :]) - y[best])
    else:
        drop = float("nan")
        rebound = float("nan")
    return {
        "peak_temperature_K": float(T[best]),
        "peak_value": float(y[best]),
        "two_sided_prominence": float(prominence),
        "post_peak_drop": float(drop),
        "high_temperature_rebound": float(rebound),
        "peak_internal": bool(peak_internal),
        "peak_in_desired_window": bool(
            peak_internal and desired_min_K <= T[best] <= desired_max_K
        ),
    }


__all__ = [name for name in globals() if not name.startswith("_")]
