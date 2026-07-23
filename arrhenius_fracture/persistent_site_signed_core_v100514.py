"""Geometry, signed projections, and persistent emission mixin."""
from __future__ import annotations

import copy
import math
from typing import Any, Callable

import numpy as np

from .persistent_site_signed_support_v100514 import (
    effective_front_width_m,
    persistent_site_multiplicity,
    solve_backstress_limited_activations,
)


class PersistentSiteSignedCoreMixin:
    def copy(self) -> "PersistentSiteSignedMPZStateV100514":
        return copy.deepcopy(self)

    @property
    def mobile(self) -> np.ndarray:
        return self.mobile_positive + self.mobile_negative

    @property
    def retained(self) -> np.ndarray:
        return self.retained_positive + self.retained_negative

    @property
    def accumulated_slip(self) -> np.ndarray:
        return self.accumulated_slip_positive + self.accumulated_slip_negative

    @property
    def mobile_count(self) -> float:
        return float(np.sum(self.mobile))

    @property
    def retained_count(self) -> float:
        return float(np.sum(self.retained))

    @property
    def available_site_fraction(self) -> float:
        return 1.0

    def _weights(self) -> tuple[np.ndarray, float, float]:
        Lb = max(self.blunting_length_m, self.dx, self.b_m)
        q = np.exp(-self.x / Lb)
        return (
            q,
            max(float(np.sum(q)), 1.0e-30),
            max(self.blunting_length_m, self.dx),
        )

    def local_unsigned_density_by_system_m2(self) -> np.ndarray:
        q, norm, width = self._weights()
        weighted = np.sum((self.mobile + self.retained) * q[None, :], axis=1) / norm
        return np.maximum(weighted / max(self.dx * width, 1.0e-30), 0.0)

    def forest_density_profile_m2(self) -> np.ndarray:
        width = max(self.blunting_length_m, self.dx)
        unsigned = np.sum(self.mobile + self.retained, axis=0)
        return np.maximum(
            float(self.candidate.rho_forest_floor_m2)
            + unsigned / max(self.dx * width, 1.0e-30),
            1.0,
        )

    def local_slip_count(self) -> float:
        q, _, _ = self._weights()
        return float(np.sum(self.accumulated_slip * q[None, :]))

    def blunted_radius(
        self,
        r0: float | None = None,
        c_blunt: float | None = None,
        b: float | None = None,
    ) -> float:
        radius0 = self.r0_m if r0 is None else max(float(r0), self.b_m)
        coefficient = self.candidate.c_blunt if c_blunt is None else float(c_blunt)
        burgers = self.b_m if b is None else abs(float(b))
        return float(
            max(
                radius0
                + max(coefficient, 0.0) * burgers * self.local_slip_count(),
                radius0,
            )
        )

    def source_geometry(self) -> dict[str, float]:
        rho_by_system = self.local_unsigned_density_by_system_m2()
        rho_width = max(
            float(self.candidate.rho_forest_floor_m2)
            + float(np.sum(rho_by_system)),
            self.reference_density_m2,
        )
        width = effective_front_width_m(
            rho_width,
            reference_width_m=self.reference_width_m,
            reference_density_m2=self.reference_density_m2,
            minimum_physical_width_m=self.minimum_front_width_m,
            burgers_m=self.b_m,
            maximum_width_m=self.maximum_front_width_m,
        )
        radius = self.blunted_radius()
        multiplicity = persistent_site_multiplicity(
            self.candidate.rho_source0_m2,
            self.active_arc_factor,
            radius,
            width,
        )
        return {
            "rho_width_m2": rho_width,
            "front_width_m": width,
            "tip_radius_m": radius,
            "source_area_m2": self.active_arc_factor * radius * width,
            "multiplicity_per_system": multiplicity,
            "active_arc_factor": self.active_arc_factor,
            "front_width_grid_independent": True,
        }

    def backstress(self) -> tuple[np.ndarray, np.ndarray]:
        rho = self.local_unsigned_density_by_system_m2()
        kback = self.G_Pa * self.b_m / max(
            float(self.candidate.taylor_stress_fraction), 1.0e-12
        )
        return rho, kback * np.sqrt(np.maximum(rho, 0.0))

    def shielding_K(
        self,
        G_shear: float | None = None,
        nu: float | None = None,
        b: float | None = None,
    ) -> float:
        signed_active = self.retained_positive - self.retained_negative
        if self.candidate.mobile_shield_fraction != 0.0:
            signed_active = signed_active + self.candidate.mobile_shield_fraction * (
                self.mobile_positive - self.mobile_negative
            )
        K = float(
            np.sum(
                self.kernel.active_kernel_Pa_sqrt_m_per_signed_line
                * signed_active
            )
        )
        if self.wake_shielding:
            signed_wake = self.wake_retained_positive - self.wake_retained_negative
            if self.candidate.mobile_shield_fraction != 0.0:
                signed_wake = (
                    signed_wake
                    + self.candidate.mobile_shield_fraction
                    * (self.wake_mobile_positive - self.wake_mobile_negative)
                )
            wake_kernel = self.kernel.wake_kernel_Pa_sqrt_m_per_signed_line
            n = min(wake_kernel.shape[1], signed_wake.shape[1])
            K += float(np.sum(wake_kernel[:, :n] * signed_wake[:, :n]))
        return K

    def _source_bin_count(self) -> int:
        return max(
            min(int(math.ceil(self.source_zone_length_m / self.dx)), self.n_bins),
            1,
        )

    def _density_increment_per_activation(self, system: int) -> float:
        q, norm, width = self._weights()
        nsrc = self._source_bin_count()
        source_weight = float(np.sum(q[:nsrc])) / float(nsrc)
        return (
            float(self.kernel.activation_to_line_content_by_system[int(system)])
            * source_weight
            / max(norm * self.dx * width, 1.0e-30)
        )

    def emit_persistent(
        self,
        *,
        dt_s: float,
        T_K: float,
        opening_stress_Pa: float,
        drive_factors: np.ndarray,
        tau_signed_Pa: np.ndarray,
        rate_function: Callable[[float, float], float],
        tolerance: float = 1.0e-10,
        max_iterations: int = 96,
    ) -> dict[str, Any]:
        dt = max(float(dt_s), 0.0)
        factors = np.asarray(drive_factors, dtype=float).reshape(-1)
        tau = np.asarray(tau_signed_Pa, dtype=float).reshape(-1)
        if factors.shape != (self.n_systems,) or tau.shape != (self.n_systems,):
            raise ValueError("two signed anisotropic drive channels are required")
        if not np.all(np.isfinite(factors)) or not np.all(np.isfinite(tau)):
            raise ValueError("anisotropic drive channels must be finite")
        geometry = self.source_geometry()
        multiplicity = float(geometry["multiplicity_per_system"])
        rho0, sigma_back0 = self.backstress()
        kback = self.G_Pa * self.b_m / max(
            float(self.candidate.taylor_stress_fraction), 1.0e-12
        )
        drive = np.maximum(
            factors * max(float(opening_stress_Pa), 0.0), 0.0
        )
        signs = np.sign(tau)
        activations = np.zeros(self.n_systems)
        blocking = np.zeros(self.n_systems)
        blocked = np.zeros(self.n_systems, dtype=bool)
        rates_initial = np.zeros(self.n_systems)
        rates_final = np.zeros(self.n_systems)
        sigma_initial = np.zeros(self.n_systems)
        sigma_final = np.zeros(self.n_systems)
        line_content = np.zeros(self.n_systems)
        nsrc = self._source_bin_count()
        for system in range(self.n_systems):
            sigma_initial[system] = max(
                drive[system] - sigma_back0[system], 0.0
            )
            if (
                signs[system] == 0.0
                or sigma_initial[system] <= 0.0
                or dt <= 0.0
            ):
                continue
            rates_initial[system] = max(
                float(rate_function(sigma_initial[system], T_K)), 0.0
            )
            rho_per = self._density_increment_per_activation(system)
            activation, block, is_blocked = solve_backstress_limited_activations(
                multiplicity=multiplicity,
                dt_s=dt,
                drive_stress_Pa=float(drive[system]),
                rho_initial_m2=float(rho0[system]),
                rho_increment_per_activation_m2=rho_per,
                backstress_prefactor_Pa_sqrt_m2=kback,
                rate_function=lambda stress, T=T_K: rate_function(stress, T),
                tolerance=tolerance,
                max_iterations=max_iterations,
            )
            activations[system] = activation
            blocking[system] = block
            blocked[system] = is_blocked
            line_content[system] = activation * float(
                self.kernel.activation_to_line_content_by_system[system]
            )
            rho_final = rho0[system] + rho_per * activation
            sigma_final[system] = max(
                drive[system] - kback * math.sqrt(max(rho_final, 0.0)), 0.0
            )
            rates_final[system] = max(
                float(rate_function(sigma_final[system], T_K)), 0.0
            )
            amount = line_content[system] / float(nsrc)
            if signs[system] > 0.0:
                self.mobile_positive[system, :nsrc] += amount
                self.accumulated_slip_positive[system, :nsrc] += amount
            else:
                self.mobile_negative[system, :nsrc] += amount
                self.accumulated_slip_negative[system, :nsrc] += amount
        self.available_sites = self.site_capacity.copy()
        self.tip_source_activity = np.ones(self.n_systems)
        emitted = float(np.sum(line_content))
        self.emitted_total += emitted
        out = {
            "dN_emit": emitted,
            "source_activations": float(np.sum(activations)),
            "activations_by_system": activations,
            "blocking_activations_by_system": blocking,
            "mechanical_blocking_active_by_system": blocked,
            "line_content_by_system": line_content,
            "burgers_sign_by_system": signs,
            "drive_stress_by_system_Pa": drive,
            "rho_back_initial_by_system_m2": rho0,
            "sigma_back_initial_by_system_Pa": sigma_back0,
            "sigma_emit_initial_by_system_Pa": sigma_initial,
            "sigma_emit_final_by_system_Pa": sigma_final,
            "rate_initial_by_system_s": rates_initial,
            "rate_final_by_system_s": rates_final,
            "aggregate_hazard_initial_by_system_s": multiplicity * rates_initial,
            "persistent_site_geometry": geometry,
            "source_sites_refreshed": 0.0,
            "available_site_fraction": 1.0,
            "finite_source_inventory_active": False,
            "source_depletion_active": False,
            "source_refresh_active": False,
        }
        self.last_emission = copy.deepcopy(out)
        return out


__all__ = ["PersistentSiteSignedCoreMixin"]
