"""Signed mobile/retained state for the v9.12 emergent-GND campaign."""
from __future__ import annotations

import math
from typing import Mapping

import numpy as np
from scipy.optimize import brentq
from scipy.special import gammainc

from .emergent_gnd_types_v912 import (
    CandidateParameters,
    CommonPhysics,
    EV_TO_J,
    KB_J_PER_K,
    KB_EV_PER_K,
    ExpFloorSurface,
    positive,
)


class EmergentGNDState:
    """Front-attached physical-density MPZ with no fitted saturation law."""

    def __init__(self, candidate: CandidateParameters, physics: CommonPhysics):
        physics.validate()
        self.p = candidate
        self.c = physics
        self.dx = physics.mpz_length_m / physics.n_bins
        self.x = (np.arange(physics.n_bins, dtype=float) + 0.5) * self.dx
        shape = (physics.n_systems, 2, physics.n_bins)
        self.mobile_m2 = np.zeros(shape, dtype=float)
        self.retained_m2 = np.zeros(shape, dtype=float)
        self.source_available_m2 = np.full(
            (physics.n_systems, physics.n_bins),
            candidate.rho_source0_m2,
            dtype=float,
        )
        self.source_capacity_m2 = self.source_available_m2.copy()
        self.time_s = 0.0
        self.extension_m = 0.0
        self._gnd_kernel = self._analytical_gnd_kernel()
        self._shield_kernel = self._analytical_shield_kernel()

    @property
    def cell_area_m2(self) -> float:
        return self.dx * self.c.active_strip_width_m

    def _analytical_gnd_kernel(self) -> np.ndarray:
        delta = self.x[:, None] - self.x[None, :]
        rc = self.c.core_regularization_b * self.c.b_m
        prefactor = self.c.G_Pa * self.c.b_m / (
            2.0 * math.pi * max(1.0 - self.c.nu, 1.0e-12)
        )
        return (
            prefactor
            * delta
            / (delta * delta + rc * rc)
            * self.cell_area_m2
        )

    def _analytical_shield_kernel(self) -> np.ndarray:
        x = np.maximum(self.x, self.c.core_regularization_b * self.c.b_m)
        per_line = (
            self.c.G_Pa
            * self.c.b_m
            / max(1.0 - self.c.nu, 1.0e-12)
            / np.sqrt(2.0 * math.pi * x)
        )
        return per_line * self.cell_area_m2 / 1.0e6

    def set_mechanical_kernels(
        self,
        *,
        gnd_stress_kernel_Pa_per_m2: np.ndarray | None = None,
        shielding_kernel_MPa_sqrt_m_per_m2: np.ndarray | None = None,
    ) -> None:
        if gnd_stress_kernel_Pa_per_m2 is not None:
            kernel = np.asarray(gnd_stress_kernel_Pa_per_m2, dtype=float)
            if kernel.shape != (self.c.n_bins, self.c.n_bins):
                raise ValueError("GND kernel must have shape (n_bins, n_bins)")
            self._gnd_kernel = kernel.copy()
        if shielding_kernel_MPa_sqrt_m_per_m2 is not None:
            kernel = np.asarray(
                shielding_kernel_MPa_sqrt_m_per_m2,
                dtype=float,
            )
            if kernel.shape != (self.c.n_bins,):
                raise ValueError("shielding kernel must have shape (n_bins,)")
            self._shield_kernel = kernel.copy()

    def signed_gnd_m2(self) -> np.ndarray:
        return self.retained_m2[:, 1, :] - self.retained_m2[:, 0, :]

    def forest_density_m2(self) -> np.ndarray:
        total = np.sum(self.mobile_m2 + self.retained_m2, axis=1)
        interaction = np.asarray(self.c.forest_interaction_matrix, dtype=float)
        return np.maximum(
            self.c.rho_forest_floor_m2 + interaction @ total,
            1.0,
        )

    def tau_gnd_Pa(self) -> np.ndarray:
        own = np.asarray(
            [self._gnd_kernel @ row for row in self.signed_gnd_m2()]
        )
        return np.asarray(self.c.gnd_stress_projection_matrix, dtype=float) @ own

    def K_shield_MPa_sqrt_m(self) -> float:
        factors = np.asarray(self.c.shielding_orientation_factors, dtype=float)
        return float(
            np.sum(
                factors[:, None]
                * self.signed_gnd_m2()
                * self._shield_kernel[None, :]
            )
        )

    def source_available_fraction(self) -> float:
        nsrc = max(
            1,
            min(
                self.c.n_bins,
                int(math.ceil(self.c.source_zone_length_m / self.dx)),
            ),
        )
        cap = float(np.sum(self.source_capacity_m2[:, :nsrc]))
        return 0.0 if cap <= 0.0 else float(
            np.sum(self.source_available_m2[:, :nsrc]) / cap
        )

    @staticmethod
    def _arrhenius_rate(
        barrier_eV: np.ndarray | float,
        T_K: float,
        nu0_s: float,
    ) -> np.ndarray:
        exponent = -np.asarray(barrier_eV, dtype=float) * EV_TO_J / (
            KB_J_PER_K * float(T_K)
        )
        return positive(nu0_s) * np.exp(np.clip(exponent, -700.0, 0.0))

    @classmethod
    def _signed_rate(
        cls,
        surface: ExpFloorSurface,
        signed_stress_Pa: np.ndarray,
        T_K: float,
        nu0_s: float,
    ) -> np.ndarray:
        forward = cls._arrhenius_rate(
            surface.barrier_eV(np.maximum(signed_stress_Pa, 0.0), T_K),
            T_K,
            nu0_s,
        )
        reverse = cls._arrhenius_rate(
            surface.barrier_eV(np.maximum(-signed_stress_Pa, 0.0), T_K),
            T_K,
            nu0_s,
        )
        return forward - reverse

    def _recovery_rate_s(self, T_K: float) -> float:
        if self.p.recovery_nu0_s <= 0.0:
            return 0.0
        G = max(
            self.p.recovery_H0_eV
            - self.p.recovery_activation_entropy_kB
            * KB_EV_PER_K
            * (float(T_K) - self.p.cleavage.Tref_K),
            0.0,
        )
        return float(self._arrhenius_rate(G, T_K, self.p.recovery_nu0_s))

    def local_rates(
        self,
        K_MPa_sqrt_m: float,
        T_K: float,
    ) -> dict[str, np.ndarray]:
        K_eff = max(K_MPa_sqrt_m - self.K_shield_MPa_sqrt_m(), 0.0)
        sigma_open = K_eff * 1.0e6 / math.sqrt(
            2.0 * math.pi * self.c.r0_m
        )
        tau_external = (
            np.asarray(self.c.emission_schmid_factors, dtype=float)[:, None]
            * sigma_open
        )
        tau_eff = tau_external + self.tau_gnd_Pa()
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

        emission = np.zeros(
            (self.c.n_systems, 2, self.c.n_bins),
            dtype=float,
        )
        nsrc = max(
            1,
            min(
                self.c.n_bins,
                int(math.ceil(self.c.source_zone_length_m / self.dx)),
            ),
        )
        for system, sign in enumerate(self.c.emission_signs):
            q = 1 if sign > 0 else 0
            stress = float(sign) * tau_eff[system, :nsrc]
            emission[system, q, :nsrc] = self._arrhenius_rate(
                self.p.emission.barrier_eV(np.maximum(stress, 0.0), T_K),
                T_K,
                self.c.emission_nu0_s,
            )
        return {
            "velocity_m_s": velocity,
            "taylor_completion_s": taylor,
            "encounter_s": encounter,
            "emission_rate_s": emission,
            "recovery_rate_s": np.asarray(self._recovery_rate_s(T_K)),
        }

    def _substep(self, rates: Mapping[str, np.ndarray], remaining_s: float) -> float:
        maximum = max(
            float(np.max(np.abs(rates["emission_rate_s"]))),
            float(np.max(np.abs(rates["taylor_completion_s"]))),
            float(np.max(np.abs(rates["encounter_s"]))),
            float(rates["recovery_rate_s"]),
        )
        vmax = float(np.max(np.abs(rates["velocity_m_s"])))
        if vmax > 0.0:
            maximum = max(maximum, vmax / max(self.dx, 1.0e-30))
        if maximum <= 0.0:
            return remaining_s
        return max(
            min(remaining_s, self.c.max_fractional_state_change / maximum),
            self.c.min_substep_s,
        )

    @staticmethod
    def _advect(
        field: np.ndarray,
        velocity: np.ndarray,
        dx: float,
        dt: float,
    ) -> np.ndarray:
        if velocity.shape != field.shape:
            raise ValueError("velocity field must match mobile density shape")
        out = field.copy()
        for system in range(field.shape[0]):
            for q in range(2):
                f = field[system, q]
                v = velocity[system, q]
                face_v = np.empty(f.size + 1)
                face_v[1:-1] = 0.5 * (v[:-1] + v[1:])
                face_v[0], face_v[-1] = v[0], v[-1]
                flux = np.zeros(f.size + 1)
                for j in range(1, f.size):
                    flux[j] = face_v[j] * (
                        f[j - 1] if face_v[j] >= 0.0 else f[j]
                    )
                flux[-1] = max(face_v[-1], 0.0) * f[-1]
                out[system, q] = np.maximum(
                    f - dt * (flux[1:] - flux[:-1]) / dx,
                    0.0,
                )
        return out

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
            if steps > 1_000_000:
                raise RuntimeError("emergent-GND integration exceeded max substeps")
            rates = self.local_rates(K_MPa_sqrt_m, T_K)
            dt = min(self._substep(rates, remaining), remaining)

            emission = rates["emission_rate_s"]
            for system, sign in enumerate(self.c.emission_signs):
                q = 1 if sign > 0 else 0
                fraction = 1.0 - np.exp(
                    -np.minimum(emission[system, q] * dt, 700.0)
                )
                emitted = self.source_available_m2[system] * fraction
                self.source_available_m2[system] -= emitted
                self.mobile_m2[system, q] += emitted
                totals["emitted_per_m"] += float(np.sum(emitted) * self.dx)

            ke = rates["encounter_s"][:, None, :]
            kt = rates["taylor_completion_s"][:, None, :]
            total = self.mobile_m2 + self.retained_m2
            exchange = ke + kt
            frac_r = np.divide(
                ke,
                exchange,
                out=np.zeros_like(exchange),
                where=exchange > 0.0,
            )
            r_eq = frac_r * total
            decay = np.exp(-np.minimum(exchange * dt, 700.0))
            self.retained_m2 = np.clip(
                r_eq + (self.retained_m2 - r_eq) * decay,
                0.0,
                total,
            )
            self.mobile_m2 = total - self.retained_m2

            recovery = float(rates["recovery_rate_s"])
            if recovery > 0.0:
                self.retained_m2 *= math.exp(-min(recovery * dt, 700.0))

            capture = self.c.annihilation_capture_radius_b * self.c.b_m
            for system in range(self.c.n_systems):
                relative_velocity = 2.0 * np.abs(rates["velocity_m_s"][system])
                pair_rate = (
                    2.0
                    * capture
                    * relative_velocity
                    * self.retained_m2[system, 0]
                    * self.retained_m2[system, 1]
                )
                removed = np.minimum(
                    pair_rate * dt,
                    np.minimum(
                        self.retained_m2[system, 0],
                        self.retained_m2[system, 1],
                    ),
                )
                self.retained_m2[system, 0] -= removed
                self.retained_m2[system, 1] -= removed
                totals["annihilated_per_m"] += float(
                    2.0 * np.sum(removed) * self.dx
                )

            velocity = np.repeat(
                rates["velocity_m_s"][:, None, :],
                2,
                axis=1,
            )
            velocity[:, 0, :] *= -1.0
            self.mobile_m2 = self._advect(
                self.mobile_m2,
                velocity,
                self.dx,
                dt,
            )
            self.time_s += dt
            remaining -= dt
        return totals

    def translate_tip(self, da_m: float) -> None:
        da = max(float(da_m), 0.0)
        if da <= 0.0:
            return
        if self.c.n_bins == 1:
            refresh = 1.0 - math.exp(
                -da / max(self.p.source_refresh_length_m, 1.0e-30)
            )
            self.source_available_m2 = np.minimum(
                self.source_available_m2
                + (self.source_capacity_m2 - self.source_available_m2) * refresh,
                self.source_capacity_m2,
            )
            self.extension_m += da
            return
        sample = self.x + da
        for array in (self.mobile_m2, self.retained_m2):
            for system in range(self.c.n_systems):
                for q in range(2):
                    array[system, q] = np.interp(
                        sample,
                        self.x,
                        array[system, q],
                        right=0.0,
                    )
        refresh = 1.0 - math.exp(
            -da / max(self.p.source_refresh_length_m, 1.0e-30)
        )
        for system in range(self.c.n_systems):
            shifted = np.interp(
                sample,
                self.x,
                self.source_available_m2[system],
                right=self.p.rho_source0_m2,
            )
            self.source_available_m2[system] = np.minimum(
                shifted
                + (self.source_capacity_m2[system] - shifted) * refresh,
                self.source_capacity_m2[system],
            )
        self.extension_m += da

    def cleavage_rate_s(
        self,
        K_MPa_sqrt_m: float,
        T_K: float,
        *,
        neutral: bool = False,
    ) -> float:
        shield = 0.0 if neutral else self.K_shield_MPa_sqrt_m()
        K_eff = max(K_MPa_sqrt_m - shield, 0.0)
        sigma = K_eff * 1.0e6 / math.sqrt(2.0 * math.pi * self.c.r0_m)
        raw = float(
            self._arrhenius_rate(
                self.p.cleavage.barrier_eV(sigma, T_K),
                T_K,
                self.c.cleavage_nu0_s,
            )
        )
        if self.c.cleavage_hits <= 1.0 + 1.0e-12:
            return raw
        exposure = min(raw * self.c.cleavage_correlation_time_s, 1.0e12)
        return float(
            gammainc(self.c.cleavage_hits, exposure)
            / self.c.cleavage_correlation_time_s
        )

    def required_K_for_rate(
        self,
        target_rate_s: float,
        T_K: float,
        *,
        neutral: bool,
    ) -> float:
        target = max(float(target_rate_s), 0.0)
        if target <= 0.0:
            return 0.0

        def residual(K: float) -> float:
            return self.cleavage_rate_s(K, T_K, neutral=neutral) - target

        if residual(0.0) >= 0.0:
            return 0.0
        upper = 300.0
        while residual(upper) < 0.0 and upper < 1.0e5:
            upper *= 2.0
        if residual(upper) < 0.0:
            raise RuntimeError("target cleavage rate cannot be bracketed")
        return float(brentq(residual, 0.0, upper))

    def delta_K_micro_MPa_sqrt_m(
        self,
        T_K: float,
        target_rate_s: float = 1.0e-3,
    ) -> float:
        neutral = self.required_K_for_rate(target_rate_s, T_K, neutral=True)
        stateful = self.required_K_for_rate(target_rate_s, T_K, neutral=False)
        return stateful - neutral

    def diagnostics(
        self,
        residence_time_s: float,
        K_MPa_sqrt_m: float,
        T_K: float,
    ) -> dict[str, float]:
        rates = self.local_rates(K_MPa_sqrt_m, T_K)
        area = self.dx * self.c.active_strip_width_m
        return {
            "retained_line_count_per_unit_thickness": float(
                np.sum(self.retained_m2) * area
            ),
            "gnd_abs_line_count_per_unit_thickness": float(
                np.sum(np.abs(self.signed_gnd_m2())) * area
            ),
            "tau_gnd_tip_MPa": float(
                np.max(np.abs(self.tau_gnd_Pa()[:, 0])) / 1.0e6
            ),
            "K_shield_MPa_sqrt_m": self.K_shield_MPa_sqrt_m(),
            "source_available_fraction": self.source_available_fraction(),
            "pi_store_max": float(np.max(rates["encounter_s"]))
            * residence_time_s,
            "pi_release_max": float(np.max(rates["taylor_completion_s"]))
            * residence_time_s,
        }
