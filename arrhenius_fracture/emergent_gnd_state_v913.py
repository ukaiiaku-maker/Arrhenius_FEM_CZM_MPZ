"""Persistent-site, backstress-limited v9.13 emergent-GND state.

This module is the one-dimensional homolog of the v10.2.21 two-dimensional
source overlay. Nucleation-site multiplicity is persistent and state dependent;
there is no finite source inventory, crack-advance source refresh, or explicit
mobile/retained recovery. Emission is integrated implicitly against the
evolving Taylor backstress.
"""
from __future__ import annotations

import math
from typing import Callable, Mapping

import numpy as np

from .emergent_gnd_state_v912_energy import EmergentGNDState as _EnergyState
from .emergent_gnd_types_v913 import CandidateParameters, CommonPhysics


MODEL_ID = "v9.13_persistent_site_backstress_blunting"
SOURCE_MODEL = "persistent_areal_sites_backstress_limited_no_inventory"


def effective_front_width_m(
    rho_unsigned_m2: float,
    *,
    reference_width_m: float,
    reference_density_m2: float,
    minimum_width_m: float,
    maximum_width_m: float,
) -> float:
    """Density-limited front correlation width, anchored at the forest floor."""
    rho = max(float(rho_unsigned_m2), float(reference_density_m2), 1.0)
    width = float(reference_width_m) * math.sqrt(
        float(reference_density_m2) / rho
    )
    lower = max(float(minimum_width_m), 1.0e-30)
    upper = max(float(maximum_width_m), lower)
    return min(max(width, lower), upper)


def persistent_site_multiplicity(
    rho_site0_m2: float,
    tip_radius_m: float,
    front_width_m: float,
    active_arc_factor: float,
) -> float:
    area = (
        max(float(active_arc_factor), 0.0)
        * max(float(tip_radius_m), 0.0)
        * max(float(front_width_m), 0.0)
    )
    return max(float(rho_site0_m2), 0.0) * area


def solve_backstress_limited_activations(
    *,
    multiplicity: float,
    dt_s: float,
    drive_stress_Pa: float,
    rho_initial_m2: float,
    rho_increment_per_activation_m2: float,
    backstress_prefactor_Pa_sqrt_m2: float,
    rate_function: Callable[[float], float],
    tolerance: float = 1.0e-10,
    max_iterations: int = 96,
) -> float:
    """Backward-Euler mean activation count with a mechanical blocking root."""
    M = max(float(multiplicity), 0.0)
    dt = max(float(dt_s), 0.0)
    drive = max(float(drive_stress_Pa), 0.0)
    rho0 = max(float(rho_initial_m2), 0.0)
    rho_per = max(float(rho_increment_per_activation_m2), 0.0)
    kback = max(float(backstress_prefactor_Pa_sqrt_m2), 0.0)
    if M <= 0.0 or dt <= 0.0 or drive <= 0.0:
        return 0.0
    if rho_per <= 0.0 or kback <= 0.0:
        raise RuntimeError(
            "persistent-site emission requires positive line conversion "
            "and backstress coupling"
        )

    back0 = kback * math.sqrt(rho0)
    sigma0 = drive - back0
    if sigma0 <= 0.0:
        return 0.0
    rate0 = max(float(rate_function(sigma0)), 0.0)
    if not math.isfinite(rate0) or rate0 <= 0.0:
        return 0.0

    rho_block = (drive / kback) ** 2
    upper = max((rho_block - rho0) / rho_per, 0.0)
    if upper <= 0.0:
        return 0.0

    def residual(value: float) -> float:
        rho = rho0 + rho_per * max(float(value), 0.0)
        sigma_eff = drive - kback * math.sqrt(max(rho, 0.0))
        rate = (
            0.0
            if sigma_eff <= 0.0
            else max(float(rate_function(sigma_eff)), 0.0)
        )
        if not math.isfinite(rate):
            rate = 0.0
        return float(value) - M * rate * dt

    lo = 0.0
    hi = upper
    if residual(hi) < 0.0:
        raise RuntimeError("failed to bracket persistent-site backstress root")
    scale = max(upper, 1.0)
    for _ in range(int(max_iterations)):
        mid = 0.5 * (lo + hi)
        value = residual(mid)
        if (
            abs(value) <= float(tolerance) * scale
            or (hi - lo) <= float(tolerance) * scale
        ):
            return max(mid, 0.0)
        if value > 0.0:
            hi = mid
        else:
            lo = mid
    return max(0.5 * (lo + hi), 0.0)


class EmergentGNDState(_EnergyState):
    """Signed mobile/retained state with persistent crack-tip nucleation sites."""

    def __init__(self, candidate: CandidateParameters, physics: CommonPhysics):
        super().__init__(candidate, physics)
        self.accumulated_slip_m2 = np.zeros_like(self.mobile_m2)
        self._active_arc_factor = self.c.reference_source_area_m2 / (
            max(self.c.r0_m, self.c.b_m, 1.0e-30)
            * self.c.reference_front_width_m
        )
        if not math.isfinite(self._active_arc_factor) or self._active_arc_factor <= 0.0:
            raise RuntimeError("invalid persistent-site active arc factor")
        self.source_available_m2[...] = self.p.rho_source0_m2
        self.source_capacity_m2[...] = self.p.rho_source0_m2
        self.last_source_activations = np.zeros(self.c.n_systems)
        self.last_line_content = np.zeros(self.c.n_systems)
        self.last_emission_drive_Pa = np.zeros(self.c.n_systems)
        self.last_sigma_back_initial_Pa = np.zeros(self.c.n_systems)
        self.last_sigma_effective_final_Pa = np.zeros(self.c.n_systems)
        self.last_rate_initial_s = np.zeros(self.c.n_systems)
        self.last_rate_final_s = np.zeros(self.c.n_systems)
        self.last_tip_radius_before_advance_m = self.tip_radius_m()
        self.last_tip_radius_after_advance_m = self.last_tip_radius_before_advance_m

    @property
    def state_strip_width_m(self) -> float:
        return max(
            float(self.c.blunting_length_m),
            float(self.dx),
            abs(float(self.c.b_m)),
            1.0e-30,
        )

    @property
    def cell_area_m2(self) -> float:
        return self.dx * self.state_strip_width_m

    def integration_metadata(self) -> dict[str, object]:
        metadata = dict(super().integration_metadata())
        metadata.update(
            {
                "model_id": MODEL_ID,
                "source_model": SOURCE_MODEL,
                "finite_source_inventory": False,
                "source_depletion_on_emission": False,
                "source_refresh_on_crack_advance": False,
                "site_multiplicity_in_arrhenius_hazard": True,
                "multiplicity_geometry": "rho_source0*c_arc*r_tip*w_eff",
                "tip_radius_state": "r0+c_blunt*b*local_accumulated_slip",
                "resharpening": "moving_frame_convection_of_accumulated_slip",
                "front_width_state": (
                    "reference_width*sqrt(reference_density/rho_unsigned)"
                ),
                "backstress_population": "unsigned_mobile_plus_retained",
                "shielding_population": "signed_retained",
                "emission_integrator": "implicit_backward_euler_backstress_root",
                "explicit_recovery_active": False,
                "reference_source_area_m2": float(
                    self.c.reference_source_area_m2
                ),
                "reference_front_width_m": float(
                    self.c.reference_front_width_m
                ),
                "reference_density_m2": float(self.c.reference_density_m2),
                "active_arc_factor": float(self._active_arc_factor),
                "activation_to_line_content_per_system": list(
                    self.c.activation_to_line_content_per_system
                ),
            }
        )
        return metadata

    def source_available_fraction(self) -> float:
        return 1.0

    def _source_zone_bin_count(self) -> int:
        length = max(float(self.c.source_zone_length_m), float(self.dx))
        return max(
            min(int(math.ceil(length / float(self.dx))), self.c.n_bins),
            1,
        )

    def _tip_weights(self) -> tuple[np.ndarray, float]:
        length = self.state_strip_width_m
        weights = np.exp(-np.asarray(self.x, dtype=float) / length)
        return weights, max(float(np.sum(weights)), 1.0e-30)

    def unsigned_tip_density_by_system_m2(self) -> np.ndarray:
        weights, norm = self._tip_weights()
        total = np.sum(
            np.maximum(self.mobile_m2, 0.0)
            + np.maximum(self.retained_m2, 0.0),
            axis=1,
        )
        return np.maximum(
            np.sum(total * weights[None, :], axis=1) / norm,
            0.0,
        )

    def local_accumulated_slip_count(self) -> float:
        weights, _ = self._tip_weights()
        line_content = np.maximum(self.accumulated_slip_m2, 0.0) * self.cell_area_m2
        return float(
            np.sum(line_content * weights[None, None, :])
            * max(float(self.c.blunting_slip_fraction), 0.0)
        )

    def tip_radius_m(self) -> float:
        radius = (
            float(self.c.r0_m)
            + max(float(self.p.c_blunt), 0.0)
            * abs(float(self.c.b_m))
            * self.local_accumulated_slip_count()
        )
        return max(radius, float(self.c.r0_m))

    def source_geometry(self) -> dict[str, float]:
        rho_by_system = self.unsigned_tip_density_by_system_m2()
        rho_width = max(
            float(self.c.rho_forest_floor_m2) + float(np.sum(rho_by_system)),
            float(self.c.reference_density_m2),
        )
        minimum = max(
            float(self.c.minimum_front_width_m),
            float(self.dx),
            abs(float(self.c.b_m)),
        )
        maximum = (
            float(self.c.maximum_front_width_m)
            if self.c.maximum_front_width_m > 0.0
            else float(self.c.mpz_length_m)
        )
        width = effective_front_width_m(
            rho_width,
            reference_width_m=self.c.reference_front_width_m,
            reference_density_m2=self.c.reference_density_m2,
            minimum_width_m=minimum,
            maximum_width_m=maximum,
        )
        radius = self.tip_radius_m()
        multiplicity = persistent_site_multiplicity(
            self.p.rho_source0_m2,
            radius,
            width,
            self._active_arc_factor,
        )
        return {
            "tip_radius_m": radius,
            "front_width_m": width,
            "rho_width_m2": rho_width,
            "source_area_m2": self._active_arc_factor * radius * width,
            "multiplicity_per_system": multiplicity,
        }

    def backstress_state(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        rho = self.unsigned_tip_density_by_system_m2()
        tau = (
            max(float(self.c.persistent_backstress_scale), 0.0)
            * float(self.c.G_Pa)
            * abs(float(self.c.b_m))
            * np.sqrt(np.maximum(rho, 0.0))
        )
        resolved = max(abs(float(self.p.taylor.stress_fraction)), 1.0e-6)
        sigma = tau / resolved
        return rho, tau, sigma

    def _density_increment_per_activation(
        self,
        system: int,
        nsrc: int,
    ) -> float:
        conversion = np.asarray(
            self.c.activation_to_line_content_per_system,
            dtype=float,
        )
        weights, norm = self._tip_weights()
        source_weight = float(np.sum(weights[:nsrc])) / max(float(nsrc), 1.0)
        return (
            float(conversion[system])
            * source_weight
            / max(norm * self.cell_area_m2, 1.0e-30)
        )

    def local_rates(
        self,
        K_MPa_sqrt_m: float,
        T_K: float,
    ) -> dict[str, np.ndarray]:
        K_raw = float(K_MPa_sqrt_m)
        K_applied = max(K_raw, 0.0)
        K_eff = max(K_raw - self.K_shield_MPa_sqrt_m(), 0.0)
        radius = self.tip_radius_m()
        stress_scale = 1.0e6 / math.sqrt(2.0 * math.pi * radius)
        sigma_applied = K_applied * stress_scale
        sigma_shielded = K_eff * stress_scale

        schmid = np.asarray(
            self.c.emission_schmid_factors,
            dtype=float,
        )[:, None]
        tau_external = schmid * sigma_shielded
        tau_gnd = self.tau_gnd_Pa()
        tau_eff = tau_external + tau_gnd

        forest = self.forest_density_m2()
        spacing = 1.0 / (2.0 * np.sqrt(forest))
        jump = self.c.jump_fraction_of_forest_spacing * spacing

        p_surface = self.p.peierls.surface(self.p.emission)
        peierls = self._signed_rate(
            p_surface,
            self.p.peierls.stress_fraction * tau_eff,
            T_K,
            self.p.peierls.nu0_s,
        )
        velocity = jump * peierls

        t_surface = self.p.taylor.surface(self.p.emission)
        taylor_stress = (
            self.p.taylor.stress_fraction
            * tau_eff
            * spacing
            / max(self.c.b_m, 1.0e-30)
        )
        taylor_single = np.maximum(
            self._signed_rate(
                t_surface,
                taylor_stress,
                T_K,
                self.p.taylor.nu0_s,
            ),
            0.0,
        )
        corr_length = self.p.taylor_corr_scale / (
            2.0 * math.sqrt(max(self.p.taylor_corr_rho_c_m2, 1.0e-300))
        )
        order = 1.0 + 2.0 * corr_length * np.sqrt(forest)
        taylor = taylor_single / np.maximum(order, 1.0)
        mfp = self.c.mean_free_path_coefficient / np.sqrt(forest)
        encounter = np.abs(velocity) / np.maximum(mfp, 1.0e-30)

        rho_back, tau_back, sigma_back = self.backstress_state()
        drive = np.maximum(
            np.abs(
                np.asarray(self.c.emission_schmid_factors, dtype=float)
            )
            * sigma_shielded,
            0.0,
        )
        sigma_emit = np.maximum(drive - sigma_back, 0.0)
        per_site_rate = self._arrhenius_rate(
            self.p.emission.barrier_eV(sigma_emit, T_K),
            T_K,
            self.c.emission_nu0_s,
        )

        emission = np.zeros(
            (self.c.n_systems, 2, self.c.n_bins),
            dtype=float,
        )
        nsrc = self._source_zone_bin_count()
        for system, sign in enumerate(self.c.emission_signs):
            q = 1 if sign > 0 else 0
            emission[system, q, :nsrc] = per_site_rate[system]

        tau_applied_column = schmid * sigma_applied
        tau_shielded_column = schmid * sigma_shielded
        tau_applied = np.broadcast_to(
            tau_applied_column,
            tau_gnd.shape,
        ).copy()
        tau_shielded = np.broadcast_to(
            tau_shielded_column,
            tau_gnd.shape,
        ).copy()

        return {
            "velocity_m_s": velocity,
            "taylor_completion_s": taylor,
            "encounter_s": encounter,
            "emission_rate_s": emission,
            "recovery_rate_s": np.asarray(0.0),
            "tau_external_Pa": tau_applied,
            "tau_nonlocal_shielding_Pa": tau_shielded - tau_applied,
            "tau_gnd_Pa": tau_gnd,
            "tau_eff_Pa": tau_shielded + tau_gnd,
            "emission_drive_Pa": drive,
            "emission_effective_stress_Pa": sigma_emit,
            "emission_rate_per_site_s": np.asarray(per_site_rate, dtype=float),
            "rho_back_by_system_m2": rho_back,
            "tau_back_by_system_Pa": tau_back,
            "sigma_back_by_system_Pa": sigma_back,
            "temperature_K": np.asarray(float(T_K)),
            "tip_radius_m": np.asarray(radius),
        }

    def _emit_exact(
        self,
        rates: Mapping[str, np.ndarray],
        dt: float,
    ) -> float:
        dt = max(float(dt), 0.0)
        if dt <= 0.0:
            return 0.0

        geometry = self.source_geometry()
        multiplicity = float(geometry["multiplicity_per_system"])
        nsrc = self._source_zone_bin_count()
        drive = np.asarray(rates["emission_drive_Pa"], dtype=float)
        rho0 = np.asarray(rates["rho_back_by_system_m2"], dtype=float)
        sigma_back0 = np.asarray(
            rates["sigma_back_by_system_Pa"],
            dtype=float,
        )
        T_K = float(rates["temperature_K"])
        conversion = np.asarray(
            self.c.activation_to_line_content_per_system,
            dtype=float,
        )
        resolved = max(abs(float(self.p.taylor.stress_fraction)), 1.0e-6)
        backstress_prefactor = (
            max(float(self.c.persistent_backstress_scale), 0.0)
            * float(self.c.G_Pa)
            * abs(float(self.c.b_m))
            / resolved
        )
        if backstress_prefactor <= 0.0:
            raise RuntimeError(
                "persistent-site source requires positive Taylor backstress"
            )

        activations = np.zeros(self.c.n_systems, dtype=float)
        line_by_system = np.zeros(self.c.n_systems, dtype=float)
        rates_initial = np.zeros(self.c.n_systems, dtype=float)
        rates_final = np.zeros(self.c.n_systems, dtype=float)
        sigma_final = np.zeros(self.c.n_systems, dtype=float)

        for system in range(self.c.n_systems):
            if drive[system] <= sigma_back0[system]:
                continue
            rho_per = self._density_increment_per_activation(system, nsrc)

            def rate_at(stress: float) -> float:
                if stress <= 0.0:
                    return 0.0
                barrier = self.p.emission.barrier_eV(stress, T_K)
                return float(
                    self._arrhenius_rate(
                        barrier,
                        T_K,
                        self.c.emission_nu0_s,
                    )
                )

            sigma_initial = max(
                float(drive[system] - sigma_back0[system]),
                0.0,
            )
            rates_initial[system] = rate_at(sigma_initial)
            activations[system] = solve_backstress_limited_activations(
                multiplicity=multiplicity,
                dt_s=dt,
                drive_stress_Pa=float(drive[system]),
                rho_initial_m2=float(rho0[system]),
                rho_increment_per_activation_m2=rho_per,
                backstress_prefactor_Pa_sqrt_m2=backstress_prefactor,
                rate_function=rate_at,
                tolerance=self.c.implicit_tolerance,
                max_iterations=self.c.implicit_max_iterations,
            )
            line_by_system[system] = (
                activations[system] * conversion[system]
            )
            rho_final = (
                float(rho0[system])
                + rho_per * activations[system]
            )
            sigma_final[system] = max(
                float(drive[system])
                - backstress_prefactor * math.sqrt(max(rho_final, 0.0)),
                0.0,
            )
            rates_final[system] = rate_at(sigma_final[system])

        for system, sign in enumerate(self.c.emission_signs):
            q = 1 if sign > 0 else 0
            density_increment = (
                line_by_system[system]
                / float(nsrc)
                / self.cell_area_m2
            )
            self.mobile_m2[system, q, :nsrc] += density_increment
            self.accumulated_slip_m2[system, q, :nsrc] += density_increment

        self.last_source_activations = activations
        self.last_line_content = line_by_system
        self.last_emission_drive_Pa = drive.copy()
        self.last_sigma_back_initial_Pa = sigma_back0.copy()
        self.last_sigma_effective_final_Pa = sigma_final
        self.last_rate_initial_s = rates_initial
        self.last_rate_final_s = rates_final
        self.source_available_m2[...] = self.p.rho_source0_m2

        return float(
            np.sum(line_by_system) / max(self.state_strip_width_m, 1.0e-30)
        )

    @staticmethod
    def _translate_density_field(
        field: np.ndarray,
        distance_m: float,
        dx: float,
        length_m: float,
    ) -> np.ndarray:
        d = max(float(distance_m), 0.0)
        if d <= 0.0:
            return np.asarray(field, dtype=float).copy()
        source = np.asarray(field, dtype=float)
        out = np.zeros_like(source)
        n_bins = source.shape[-1]
        for i in range(n_bins):
            left = i * dx - d
            right = (i + 1) * dx - d
            if right <= 0.0 or left >= length_m:
                continue
            inside_left = max(left, 0.0)
            inside_right = min(right, length_m)
            if inside_right <= inside_left:
                continue
            j0 = max(int(math.floor(inside_left / dx)), 0)
            j1 = min(
                int(math.floor((inside_right - 1.0e-15 * dx) / dx)),
                n_bins - 1,
            )
            for j in range(j0, j1 + 1):
                overlap_left = max(inside_left, j * dx)
                overlap_right = min(inside_right, (j + 1) * dx)
                fraction = max(overlap_right - overlap_left, 0.0) / dx
                if fraction > 0.0:
                    out[..., j] += source[..., i] * fraction
        return out

    def translate_tip(self, da_m: float) -> None:
        da = max(float(da_m), 0.0)
        if da <= 0.0:
            return
        self.last_tip_radius_before_advance_m = self.tip_radius_m()
        for name in (
            "mobile_m2",
            "retained_m2",
            "accumulated_slip_m2",
        ):
            setattr(
                self,
                name,
                self._translate_density_field(
                    getattr(self, name),
                    da,
                    self.dx,
                    self.c.mpz_length_m,
                ),
            )
        self.source_available_m2[...] = self.p.rho_source0_m2
        self.source_capacity_m2[...] = self.p.rho_source0_m2
        self.extension_m += da
        self.last_tip_radius_after_advance_m = self.tip_radius_m()

    def advance_time(
        self,
        duration_s: float,
        K_MPa_sqrt_m: float,
        T_K: float,
    ) -> dict[str, float]:
        remaining = max(float(duration_s), 0.0)
        totals = {"emitted_per_m": 0.0, "annihilated_per_m": 0.0}
        steps = 0
        while remaining > 0.0:
            steps += 1
            if steps > 10_000_000:
                raise RuntimeError(
                    "persistent emergent-GND integration exceeded max feedback steps"
                )
            dt = self._substep({}, remaining)
            increment = self._advance_spatial_step(
                dt,
                K_MPa_sqrt_m,
                T_K,
            )
            totals["emitted_per_m"] += increment["emitted_per_m"]
            totals["annihilated_per_m"] += increment["annihilated_per_m"]
            self.time_s += dt
            remaining -= dt
        return totals

    def cleavage_rate_s(
        self,
        K_MPa_sqrt_m: float,
        T_K: float,
        *,
        neutral: bool = False,
    ) -> float:
        shield = 0.0 if neutral else self.K_shield_MPa_sqrt_m()
        radius = float(self.c.r0_m) if neutral else self.tip_radius_m()
        K_eff = max(float(K_MPa_sqrt_m) - shield, 0.0)
        sigma = K_eff * 1.0e6 / math.sqrt(
            2.0 * math.pi * max(radius, self.c.b_m, 1.0e-30)
        )
        raw = float(
            self._arrhenius_rate(
                self.p.cleavage.barrier_eV(sigma, T_K),
                T_K,
                self.c.cleavage_nu0_s,
            )
        )
        if self.c.cleavage_hits <= 1.0 + 1.0e-12:
            return raw
        from scipy.special import gammainc

        exposure = min(raw * self.c.cleavage_correlation_time_s, 1.0e12)
        return float(
            gammainc(self.c.cleavage_hits, exposure)
            / self.c.cleavage_correlation_time_s
        )

    def diagnostics(
        self,
        residence_time_s: float,
        K_MPa_sqrt_m: float,
        T_K: float,
    ) -> dict[str, float]:
        data = dict(
            super().diagnostics(
                residence_time_s,
                K_MPa_sqrt_m,
                T_K,
            )
        )
        geometry = self.source_geometry()
        rho, tau, sigma = self.backstress_state()
        ratio = np.divide(
            sigma,
            np.maximum(self.last_emission_drive_Pa, 1.0),
            out=np.zeros_like(sigma),
            where=self.last_emission_drive_Pa > 0.0,
        )
        data.update(
            {
                "source_available_fraction": 1.0,
                "persistent_source_inventory_active": 0.0,
                "persistent_source_refresh_active": 0.0,
                "persistent_site_density_m2": float(self.p.rho_source0_m2),
                "persistent_site_multiplicity_per_system": float(
                    geometry["multiplicity_per_system"]
                ),
                "persistent_site_source_area_m2": float(
                    geometry["source_area_m2"]
                ),
                "persistent_site_front_width_m": float(
                    geometry["front_width_m"]
                ),
                "persistent_site_width_density_m2": float(
                    geometry["rho_width_m2"]
                ),
                "persistent_tip_radius_m": float(geometry["tip_radius_m"]),
                "persistent_source_zone_bins": float(
                    self._source_zone_bin_count()
                ),
                "persistent_rho_back_mean_m2": float(np.mean(rho)),
                "persistent_tau_back_mean_Pa": float(np.mean(tau)),
                "persistent_sigma_back_mean_Pa": float(np.mean(sigma)),
                "persistent_backstress_drive_ratio_max": float(np.max(ratio)),
                "persistent_last_source_activations": float(
                    np.sum(self.last_source_activations)
                ),
                "persistent_last_line_content": float(
                    np.sum(self.last_line_content)
                ),
                "persistent_local_accumulated_slip_count": float(
                    self.local_accumulated_slip_count()
                ),
                "tip_radius_before_advance_m": float(
                    self.last_tip_radius_before_advance_m
                ),
                "tip_radius_after_advance_m": float(
                    self.last_tip_radius_after_advance_m
                ),
                "tip_resharpening_by_advance_m": max(
                    float(self.last_tip_radius_before_advance_m)
                    - float(self.last_tip_radius_after_advance_m),
                    0.0,
                ),
            }
        )
        return data


__all__ = [
    "MODEL_ID",
    "SOURCE_MODEL",
    "EmergentGNDState",
    "effective_front_width_m",
    "persistent_site_multiplicity",
    "solve_backstress_limited_activations",
]
