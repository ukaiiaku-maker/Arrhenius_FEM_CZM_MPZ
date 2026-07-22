"""Diagnostic energy bookkeeping for the corrected v9.12 emergent-GND state.

The bookkeeping in this module does not alter any constitutive rate or state
update.  It integrates a reduced 1-D Orowan resolved-power measure per unit
crack-front thickness and reports a logarithmic dislocation line-energy proxy.
These quantities are useful for distinguishing a loss of stationary-crack
shielding from increasing ductile dissipation, but they are not by themselves a
full FEM J-integral, a validated J-R curve, or a Charpy impact energy.
"""
from __future__ import annotations

import math
from typing import Mapping

import numpy as np

from .emergent_gnd_state_v912_stiff import EmergentGNDState as _StiffState
from .emergent_gnd_types_v912 import CandidateParameters, CommonPhysics


class EmergentGNDState(_StiffState):
    """Corrected stiff state with non-coupled work/energy diagnostics."""

    def __init__(self, candidate: CandidateParameters, physics: CommonPhysics):
        super().__init__(candidate, physics)
        self.external_plastic_work_J_per_m = 0.0
        self.nonlocal_shielding_work_J_per_m = 0.0
        self.internal_stress_work_J_per_m = 0.0
        self.effective_plastic_work_J_per_m = 0.0
        self.effective_plastic_dissipation_J_per_m = 0.0

    def integration_metadata(self) -> dict[str, float | int | str | bool]:
        metadata = dict(super().integration_metadata())
        metadata.update(
            {
                "energy_bookkeeping": (
                    "reduced_1d_orowan_power_and_log_line_energy_v1"
                ),
                "energy_bookkeeping_active": bool(self.c.n_bins > 1),
                "energy_bookkeeping_feedback_active": False,
                "energy_units_scope": "per_unit_crack_front_thickness",
                "external_work_sign_convention": (
                    "signed_tau_from_K_applied_times_gamma_dot"
                ),
                "feedback_work_decomposition": (
                    "nonlocal_K_shield_projection_plus_local_tau_gnd"
                ),
                "effective_dissipation_convention": (
                    "positive_part_tau_effective_gamma_dot"
                ),
            }
        )
        return metadata

    def local_rates(
        self,
        K_MPa_sqrt_m: float,
        T_K: float,
    ) -> dict[str, np.ndarray]:
        rates = dict(super().local_rates(K_MPa_sqrt_m, T_K))
        stress_scale = 1.0e6 / math.sqrt(2.0 * math.pi * self.c.r0_m)
        K_raw = float(K_MPa_sqrt_m)
        K_applied = max(K_raw, 0.0)
        K_eff = max(K_raw - self.K_shield_MPa_sqrt_m(), 0.0)
        sigma_applied = K_applied * stress_scale
        sigma_shielded = K_eff * stress_scale

        tau_gnd = self.tau_gnd_Pa()
        schmid = np.asarray(
            self.c.emission_schmid_factors,
            dtype=float,
        )[:, None]
        tau_applied_column = schmid * sigma_applied
        tau_shielded_column = schmid * sigma_shielded
        tau_applied = np.broadcast_to(tau_applied_column, tau_gnd.shape).copy()
        tau_shielded = np.broadcast_to(
            tau_shielded_column,
            tau_gnd.shape,
        ).copy()
        rates["tau_external_Pa"] = tau_applied
        rates["tau_nonlocal_shielding_Pa"] = tau_shielded - tau_applied
        rates["tau_gnd_Pa"] = tau_gnd
        rates["tau_eff_Pa"] = tau_shielded + tau_gnd
        return rates

    def _plastic_power_snapshot(
        self,
        rates: Mapping[str, np.ndarray],
    ) -> dict[str, float]:
        """Return resolved Orowan power per unit crack-front thickness.

        Opposite dislocation signs move in opposite spatial directions under the
        same signed Peierls velocity.  Their Burgers-sign/velocity products
        therefore add in the net plastic shear rate, giving

            gamma_dot_alpha = b * (rho_m,+ + rho_m,-) * v_alpha.
        """
        velocity = np.asarray(rates["velocity_m_s"], dtype=float)
        mobile_total = np.sum(np.maximum(self.mobile_m2, 0.0), axis=1)
        gamma_dot = self.c.b_m * mobile_total * velocity

        tau_external = np.asarray(rates["tau_external_Pa"], dtype=float)
        tau_shielding = np.asarray(
            rates["tau_nonlocal_shielding_Pa"],
            dtype=float,
        )
        tau_gnd = np.asarray(rates["tau_gnd_Pa"], dtype=float)
        tau_eff = np.asarray(rates["tau_eff_Pa"], dtype=float)
        expected = (self.c.n_systems, self.c.n_bins)
        for array, name in (
            (velocity, "velocity_m_s"),
            (tau_external, "tau_external_Pa"),
            (tau_shielding, "tau_nonlocal_shielding_Pa"),
            (tau_gnd, "tau_gnd_Pa"),
            (tau_eff, "tau_eff_Pa"),
        ):
            if array.shape != expected:
                raise ValueError(f"{name} must have shape {expected}")

        area_per_unit_thickness = self.cell_area_m2
        external_density = tau_external * gamma_dot
        shielding_density = tau_shielding * gamma_dot
        internal_density = tau_gnd * gamma_dot
        effective_density = tau_eff * gamma_dot
        return {
            "external_plastic_power_J_per_m_s": float(
                np.sum(external_density) * area_per_unit_thickness
            ),
            "nonlocal_shielding_power_J_per_m_s": float(
                np.sum(shielding_density) * area_per_unit_thickness
            ),
            "internal_stress_power_J_per_m_s": float(
                np.sum(internal_density) * area_per_unit_thickness
            ),
            "effective_plastic_power_J_per_m_s": float(
                np.sum(effective_density) * area_per_unit_thickness
            ),
            "effective_plastic_dissipation_J_per_m_s": float(
                np.sum(np.maximum(effective_density, 0.0))
                * area_per_unit_thickness
            ),
        }

    def _accumulate_work_interval(
        self,
        start: Mapping[str, float],
        end: Mapping[str, float],
        duration_s: float,
    ) -> None:
        duration = max(float(duration_s), 0.0)
        if duration <= 0.0:
            return
        mappings = (
            (
                "external_plastic_work_J_per_m",
                "external_plastic_power_J_per_m_s",
            ),
            (
                "nonlocal_shielding_work_J_per_m",
                "nonlocal_shielding_power_J_per_m_s",
            ),
            (
                "internal_stress_work_J_per_m",
                "internal_stress_power_J_per_m_s",
            ),
            (
                "effective_plastic_work_J_per_m",
                "effective_plastic_power_J_per_m_s",
            ),
            (
                "effective_plastic_dissipation_J_per_m",
                "effective_plastic_dissipation_J_per_m_s",
            ),
        )
        for attribute, key in mappings:
            increment = 0.5 * duration * (
                float(start[key]) + float(end[key])
            )
            setattr(self, attribute, float(getattr(self, attribute)) + increment)

    def _line_energy_terms_J_per_m(self) -> dict[str, float]:
        """Return a fixed-physics logarithmic line-energy proxy.

        The outer cutoff is the local forest spacing and the inner cutoff is the
        common regularized core radius.  This is a reduced stored-line-energy
        diagnostic, not the complete elastic energy of the signed GND field.
        """
        forest = self.forest_density_m2()
        outer = 1.0 / (2.0 * np.sqrt(np.maximum(forest, 1.0)))
        core = max(
            self.c.core_regularization_b * self.c.b_m,
            self.c.b_m,
            1.0e-30,
        )
        log_factor = np.log(np.maximum(outer / core, 1.0))
        line_tension_J_per_m = (
            self.c.G_Pa
            * self.c.b_m
            * self.c.b_m
            / (4.0 * math.pi * max(1.0 - self.c.nu, 1.0e-12))
            * log_factor
        )
        mobile = np.sum(np.maximum(self.mobile_m2, 0.0), axis=1)
        retained = np.sum(np.maximum(self.retained_m2, 0.0), axis=1)
        area = self.cell_area_m2
        mobile_energy = float(np.sum(mobile * line_tension_J_per_m) * area)
        retained_energy = float(
            np.sum(retained * line_tension_J_per_m) * area
        )
        return {
            "mobile_line_energy_J_per_m": mobile_energy,
            "retained_line_energy_J_per_m": retained_energy,
            "total_line_energy_J_per_m": mobile_energy + retained_energy,
        }

    def _advance_spatial_step(
        self,
        dt: float,
        K_MPa_sqrt_m: float,
        T_K: float,
    ) -> dict[str, float]:
        """Advance the v2 state and integrate diagnostic power by trapezoids."""
        totals = {"emitted_per_m": 0.0, "annihilated_per_m": 0.0}

        rates_start = self.local_rates(K_MPa_sqrt_m, T_K)
        power_start = self._plastic_power_snapshot(rates_start)
        totals["annihilated_per_m"] += self._annihilate_exact(
            rates_start, 0.5 * dt
        )
        self._coupled_mobile_retained(rates_start, 0.5 * dt)

        rates_mid = self.local_rates(K_MPa_sqrt_m, T_K)
        power_mid = self._plastic_power_snapshot(rates_mid)
        self._accumulate_work_interval(power_start, power_mid, 0.5 * dt)

        totals["emitted_per_m"] += self._emit_exact(rates_mid, dt)

        rates_post_emit = self.local_rates(K_MPa_sqrt_m, T_K)
        power_post_emit = self._plastic_power_snapshot(rates_post_emit)
        self._coupled_mobile_retained(rates_post_emit, 0.5 * dt)

        rates_end = self.local_rates(K_MPa_sqrt_m, T_K)
        power_end = self._plastic_power_snapshot(rates_end)
        self._accumulate_work_interval(
            power_post_emit, power_end, 0.5 * dt
        )
        totals["annihilated_per_m"] += self._annihilate_exact(
            rates_end, 0.5 * dt
        )
        return totals

    def diagnostics(
        self,
        residence_time_s: float,
        K_MPa_sqrt_m: float,
        T_K: float,
    ) -> dict[str, float]:
        diagnostics = dict(
            super().diagnostics(residence_time_s, K_MPa_sqrt_m, T_K)
        )
        diagnostics.update(self._line_energy_terms_J_per_m())
        diagnostics.update(
            {
                "external_plastic_work_J_per_m": float(
                    self.external_plastic_work_J_per_m
                ),
                "nonlocal_shielding_work_J_per_m": float(
                    self.nonlocal_shielding_work_J_per_m
                ),
                "internal_stress_work_J_per_m": float(
                    self.internal_stress_work_J_per_m
                ),
                "effective_plastic_work_J_per_m": float(
                    self.effective_plastic_work_J_per_m
                ),
                "effective_plastic_dissipation_J_per_m": float(
                    self.effective_plastic_dissipation_J_per_m
                ),
            }
        )
        extension = max(float(self.extension_m), 0.0)
        if extension > 0.0:
            diagnostics["external_plastic_work_per_crack_area_J_m2"] = (
                float(self.external_plastic_work_J_per_m) / extension
            )
            diagnostics[
                "effective_plastic_dissipation_per_crack_area_J_m2"
            ] = float(self.effective_plastic_dissipation_J_per_m) / extension
        else:
            diagnostics["external_plastic_work_per_crack_area_J_m2"] = 0.0
            diagnostics[
                "effective_plastic_dissipation_per_crack_area_J_m2"
            ] = 0.0
        return diagnostics


__all__ = ["EmergentGNDState"]
