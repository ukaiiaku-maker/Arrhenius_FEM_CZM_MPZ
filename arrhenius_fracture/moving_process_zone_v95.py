"""Spatial local-density moving-process-zone kinetics for v9.5.

The v9.4 developed-state screen exposed a structural inconsistency: the
Taylor forest density was calculated from total retained count divided by the
square of the full process-zone length, while shielding was calculated from
individual near-tip retained lines.  For a 50--200 micrometre MPZ this kept the
Taylor law at its density floor even when a physically meaningful near-tip
retained population existed.

This implementation preserves the conservative source, mobile, retained,
shielding, blunting, moving-frame, and detailed-balance PT architecture but
calculates forest density and stress locally in every MPZ bin.
"""
from __future__ import annotations

import math
from typing import Any

import numpy as np

from .config import EV_TO_J, KB
from .moving_process_zone import MovingProcessZoneState as _BaseState


class MovingProcessZoneState(_BaseState):
    """v9.5 MPZ state with bin-local forest density and stress."""

    def local_forest_density_m2(self) -> np.ndarray:
        """Forest density from retained lines in each physical MPZ bin.

        The effective out-of-plane/process-zone width uses the existing
        blunting length.  One retained line in one bin therefore contributes
        ``1/(dx*width)`` rather than ``1/L_pz**2`` to the local density.
        """
        width = max(float(self.cfg.blunting_length_m), self.dx, 1.0e-12)
        local_count = np.sum(np.maximum(self.retained, 0.0), axis=0)
        return np.maximum(
            float(self.cfg.pt_forest_density_floor_m2)
            + local_count / max(self.dx * width, 1.0e-30),
            1.0,
        )

    def local_stress_profile_Pa(self, tip_stress_Pa: float) -> np.ndarray:
        """LEFM-like decay of the supplied effective tip stress over the MPZ."""
        ref = max(float(self.cfg.blunting_length_m), self.dx, 1.0e-12)
        return max(float(tip_stress_Pa), 0.0) * np.sqrt(
            ref / np.maximum(ref + self.x, ref)
        )

    @staticmethod
    def _activated_rate_array(
        nu0: float,
        barrier_eV: float,
        activation_volume_b3: float,
        stress_Pa: np.ndarray,
        T_K: float,
        b: float,
    ) -> np.ndarray:
        H = max(float(barrier_eV), 0.0) * EV_TO_J
        V = max(float(activation_volume_b3), 0.0) * abs(float(b)) ** 3
        G = np.maximum(H - np.maximum(stress_Pa, 0.0) * V, 0.0)
        exponent = -G / max(KB * float(T_K), 1.0e-30)
        return max(float(nu0), 0.0) * np.exp(
            np.clip(exponent, -700.0, 0.0)
        )

    def initialize_forest_profile(
        self,
        rho_tip_m2: float,
        decay_length_m: float | None = None,
        available_site_fraction: float = 1.0,
        slip_per_retained: float = 1.0,
    ) -> None:
        """Initialize a physically dimensioned retained-density branch.

        This is used only by developed-state continuation audits.  Virgin
        production runs retain the default zero-state initialization.
        """
        floor = float(self.cfg.pt_forest_density_floor_m2)
        excess = max(float(rho_tip_m2) - floor, 0.0)
        length = max(
            float(decay_length_m)
            if decay_length_m is not None
            else 0.15 * self.length_m,
            self.dx,
        )
        profile = excess * np.exp(-self.x / length)
        width = max(float(self.cfg.blunting_length_m), self.dx, 1.0e-12)
        counts = profile * self.dx * width
        system_weight = np.maximum(np.abs(self.orientation_factors), 1.0e-12)
        system_weight /= np.sum(system_weight)
        self.retained = system_weight[:, None] * counts[None, :]
        self.mobile.fill(0.0)
        self.accumulated_slip = (
            max(float(slip_per_retained), 0.0) * self.retained.copy()
        )
        frac = float(np.clip(available_site_fraction, 0.0, 1.0))
        self.available_sites = self.site_capacity * frac

    def evolve(
        self,
        dt_s: float,
        T_K: float,
        stress_Pa: float,
        b: float,
        emission_hazard_integral: float = 0.0,
        system_weights: np.ndarray | None = None,
    ) -> dict[str, float]:
        """Advance the spatial MPZ with local detailed-balance PT kinetics."""
        dt_s = max(float(dt_s), 0.0)
        emitted = self._source_commit_from_hazard(
            emission_hazard_integral, system_weights
        )
        source_recovered = self._recover_source_sites(dt_s)

        stress_profile = self.local_stress_profile_Pa(stress_Pa)
        rho_profile = self.local_forest_density_m2()
        pt_diag: dict[str, float] = {}

        if bool(getattr(self.cfg, "use_emission_derived_pt", True)):
            from .emission_derived_plasticity import (
                CorrelatedTaylorConfig,
                EmissionDerivedPeierlsTaylorConfig,
                EmissionDerivedPeierlsTaylorModel,
                ExpFloorSurface,
                MechanismScale,
            )

            pt_cfg = EmissionDerivedPeierlsTaylorConfig(
                parent=ExpFloorSurface(
                    G00_eV=self.cfg.pt_emit_G00_eV,
                    gT_eV_per_K=self.cfg.pt_emit_gT_eV_per_K,
                    sigc0_Pa=self.cfg.pt_emit_sigc0_Pa,
                    sT_Pa_per_K=self.cfg.pt_emit_sT_Pa_per_K,
                    Tref_K=self.cfg.pt_emit_Tref_K,
                    a=self.cfg.pt_emit_exp_a,
                    n=self.cfg.pt_emit_exp_n,
                    floor_fraction=self.cfg.pt_emit_floor_frac,
                    floor_min_eV=self.cfg.pt_emit_floor_min_eV,
                    floor_max_fraction=self.cfg.pt_emit_floor_max_frac,
                ),
                peierls=MechanismScale(
                    self.cfg.pt_peierls_energy_ratio,
                    self.cfg.pt_peierls_entropy_ratio,
                    self.cfg.pt_peierls_stress_ratio,
                    self.cfg.pt_peierls_nu0_s,
                ),
                taylor=MechanismScale(
                    self.cfg.pt_taylor_energy_ratio,
                    self.cfg.pt_taylor_entropy_ratio,
                    self.cfg.pt_taylor_stress_ratio,
                    self.cfg.pt_taylor_nu0_s,
                ),
                correlated_taylor=CorrelatedTaylorConfig(
                    self.cfg.pt_taylor_corr_rho_c,
                    self.cfg.pt_taylor_renewal_time_s,
                    self.cfg.pt_taylor_m_exponent,
                    self.cfg.pt_taylor_m_scale,
                    self.cfg.pt_taylor_m_cap,
                ),
                peierls_stress_fraction=self.cfg.pt_peierls_stress_fraction,
                taylor_stress_fraction=self.cfg.pt_taylor_stress_fraction,
                taylor_phi_max=self.cfg.pt_taylor_phi_max,
                mobile_fraction_low_density=self.cfg.pt_mobile_fraction,
                mobile_saturation_density_m2=(
                    self.cfg.pt_mobile_saturation_density_m2
                ),
                mobile_density_floor_m2=self.cfg.pt_mobile_density_floor_m2,
                jump_fraction_of_forest_spacing=self.cfg.pt_jump_fraction,
                jump_length_min_m=self.cfg.pt_jump_length_min_m,
            )
            pt_model = EmissionDerivedPeierlsTaylorModel(pt_cfg)
            rates = pt_model.rates(stress_profile, rho_profile, T_K, b)
            peierls_profile = np.asarray(
                rates["peierls_rate_s"], dtype=float
            ).reshape(-1)
            detrap_profile = np.asarray(
                rates["taylor_completion_rate_s"], dtype=float
            ).reshape(-1)
            series_profile = np.asarray(
                rates["series_rate_s"], dtype=float
            ).reshape(-1)

            mobile_by_bin = np.sum(np.maximum(self.mobile, 0.0), axis=0)
            if float(np.sum(mobile_by_bin)) > 0.0:
                glide_rate = float(
                    np.sum(peierls_profile * mobile_by_bin)
                    / np.sum(mobile_by_bin)
                )
            else:
                nsrc = max(
                    min(int(self.cfg.source_bin_count), self.n_bins), 1
                )
                glide_rate = float(np.mean(peierls_profile[:nsrc]))
            detrap_rate = float(np.max(detrap_profile))
            series_rate = float(np.max(series_profile))
            pt_diag = {
                "peierls_rate_s": glide_rate,
                "peierls_rate_min_s": float(np.min(peierls_profile)),
                "peierls_rate_max_s": float(np.max(peierls_profile)),
                "taylor_completion_rate_s": detrap_rate,
                "taylor_completion_rate_min_s": float(
                    np.min(detrap_profile)
                ),
                "taylor_completion_rate_max_s": float(
                    np.max(detrap_profile)
                ),
                "peierls_taylor_series_rate_s": series_rate,
                "taylor_m_eff": float(
                    np.max(np.asarray(rates["taylor_m_eff"], dtype=float))
                ),
                "G_peierls_eV": float(
                    np.min(np.asarray(rates["G_peierls_eV"], dtype=float))
                ),
                "G_taylor_eV": float(
                    np.min(np.asarray(rates["G_taylor_eV"], dtype=float))
                ),
                "rho_forest_m2": float(np.max(rho_profile)),
                "rho_forest_min_m2": float(np.min(rho_profile)),
                "rho_forest_median_m2": float(np.median(rho_profile)),
                "rho_forest_max_m2": float(np.max(rho_profile)),
                "local_stress_min_Pa": float(np.min(stress_profile)),
                "local_stress_max_Pa": float(np.max(stress_profile)),
            }
        else:
            peierls_profile = np.full(
                self.n_bins,
                self._activated_rate(
                    self.cfg.glide_nu0_s,
                    self.cfg.glide_barrier_eV,
                    self.cfg.glide_activation_volume_b3,
                    self.cfg.glide_stress_fraction * max(float(stress_Pa), 0.0),
                    T_K,
                    b,
                ),
            )
            detrap_profile = np.full(
                self.n_bins,
                self._activated_rate(
                    self.cfg.detrap_nu0_s,
                    self.cfg.detrap_barrier_eV,
                    self.cfg.detrap_activation_volume_b3,
                    max(float(stress_Pa), 0.0),
                    T_K,
                    b,
                ),
            )
            glide_rate = float(peierls_profile[0])
            detrap_rate = float(detrap_profile[0])
            series_rate = 0.0 if glide_rate <= 0 or detrap_rate <= 0 else (
                1.0 / (1.0 / glide_rate + 1.0 / detrap_rate)
            )

        velocity = max(float(self.cfg.glide_step_m), 0.0) * glide_rate
        trap_profile = self._activated_rate_array(
            self.cfg.trap_nu0_s,
            self.cfg.trap_barrier_eV,
            self.cfg.trap_activation_volume_b3,
            stress_profile,
            T_K,
            b,
        )
        retained_recovery_profile = self._activated_rate_array(
            self.cfg.retained_recovery_nu0_s,
            self.cfg.retained_recovery_barrier_eV,
            self.cfg.retained_recovery_activation_volume_b3,
            stress_profile,
            T_K,
            b,
        )

        if dt_s > 0.0:
            ftrap = 1.0 - np.exp(-np.minimum(trap_profile * dt_s, 700.0))
            fdetrap = 1.0 - np.exp(
                -np.minimum(detrap_profile * dt_s, 700.0)
            )
            frec_r = 1.0 - np.exp(
                -np.minimum(retained_recovery_profile * dt_s, 700.0)
            )
            frec_m = 1.0 - math.exp(
                -min(
                    max(self.cfg.mobile_recovery_rate_s, 0.0) * dt_s,
                    700.0,
                )
            )
        else:
            ftrap = np.zeros(self.n_bins)
            fdetrap = np.zeros(self.n_bins)
            frec_r = np.zeros(self.n_bins)
            frec_m = 0.0

        dm_to_r = self.mobile * ftrap[None, :]
        dr_to_m = self.retained * fdetrap[None, :]
        dr_rec = np.maximum(self.retained - dr_to_m, 0.0) * frec_r[None, :]
        dm_rec = np.maximum(self.mobile - dm_to_r, 0.0) * frec_m
        self.mobile += dr_to_m - dm_to_r - dm_rec
        self.retained += dm_to_r - dr_to_m - dr_rec
        trapped = float(np.sum(dm_to_r))
        detrapped = float(np.sum(dr_to_m))
        recovered = float(np.sum(dr_rec) + np.sum(dm_rec))

        annihilated = 0.0
        kpair = max(float(self.cfg.pair_annihilation_rate_per_count_s), 0.0)
        if kpair > 0.0 and self.n_systems >= 2 and dt_s > 0.0:
            for i in range(0, self.n_systems - 1, 2):
                pair = np.minimum(self.retained[i], self.retained[i + 1])
                frac = 1.0 - np.exp(-np.minimum(kpair * pair * dt_s, 700.0))
                removed = pair * frac
                self.retained[i] -= removed
                self.retained[i + 1] -= removed
                annihilated += 2.0 * float(np.sum(removed))

        self.mobile, escaped = self._advect_forward_field(
            self.mobile, velocity * dt_s
        )
        self.mobile = np.maximum(self.mobile, 0.0)
        self.retained = np.maximum(self.retained, 0.0)
        self.escaped_total += escaped
        self.recovered_total += recovered + annihilated
        self.time_s += dt_s
        return {
            "dN_emit": float(np.sum(emitted)),
            "dN_source_recovered": source_recovered,
            "dN_trapped": trapped,
            "dN_detrapped": detrapped,
            "dN_escaped": escaped,
            "dN_recovered": recovered,
            "dN_annihilated": annihilated,
            "glide_rate_s": glide_rate,
            "glide_velocity_m_s": velocity,
            "trap_rate_s": float(np.max(trap_profile)),
            "detrap_rate_s": detrap_rate,
            "retained_recovery_rate_s": float(
                np.max(retained_recovery_profile)
            ),
            "peierls_taylor_series_rate_s": series_rate,
            "transport_substeps": 1.0,
            "available_site_fraction": self.available_site_fraction,
            **pt_diag,
        }


__all__ = ["MovingProcessZoneState"]
