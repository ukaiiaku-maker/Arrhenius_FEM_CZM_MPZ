"""Unified Peierls-transport/Taylor-retention moving process zone (v9.10).

The v9.9 validation exposed a mismatch between the analytical retention proxy
and the spatial MPZ: the proxy treated slow PT escape as retained shielding,
while the spatial model required a separate legacy Arrhenius trap law.  This
module removes that independent trap barrier.

The state equations are now

    dN_m/dt = R_emit - k_enc N_m + k_T N_r - k_esc N_m - k_mrec N_m
    dN_r/dt =          k_enc N_m - k_T N_r             - k_rrec N_r

with

    v_P   = jump_length * lambda_P
    k_enc = eta_enc * v_P * sqrt(rho_f)
    k_T   = lambda_T_completion
    k_esc = v_P / L_MPZ.

Thus Peierls motion transports mobile dislocations and creates obstacle
encounters; Taylor completion releases retained dislocations.  A frozen
Peierls branch cannot manufacture retained shielding.  There is no independent
trap activation barrier and no constitutive cap.
"""
from __future__ import annotations

import math
from typing import Any

import numpy as np

from .emission_derived_plasticity import (
    CorrelatedTaylorConfig,
    EmissionDerivedPeierlsTaylorConfig,
    ExpFloorSurface,
)
from .emission_derived_plasticity_v97 import (
    EmissionDerivedPeierlsTaylorModel,
    IndependentEntropyMechanismScale,
)
from .moving_process_zone_v95 import MovingProcessZoneState as _SpatialStateV95


class MovingProcessZoneState(_SpatialStateV95):
    """v9.5 spatial state with unified encounter/retention kinetics."""

    state_model = "moving_pz_v910_unified_transport_retention"

    def _pt_model(self) -> EmissionDerivedPeierlsTaylorModel:
        cfg = self.cfg
        parent = ExpFloorSurface(
            G00_eV=float(cfg.pt_emit_G00_eV),
            gT_eV_per_K=float(cfg.pt_emit_gT_eV_per_K),
            sigc0_Pa=float(cfg.pt_emit_sigc0_Pa),
            sT_Pa_per_K=float(cfg.pt_emit_sT_Pa_per_K),
            Tref_K=float(cfg.pt_emit_Tref_K),
            a=float(cfg.pt_emit_exp_a),
            n=float(cfg.pt_emit_exp_n),
            floor_fraction=float(cfg.pt_emit_floor_frac),
            floor_min_eV=float(cfg.pt_emit_floor_min_eV),
            floor_max_fraction=float(cfg.pt_emit_floor_max_frac),
        )
        return EmissionDerivedPeierlsTaylorModel(
            EmissionDerivedPeierlsTaylorConfig(
                parent=parent,
                peierls=IndependentEntropyMechanismScale(
                    energy_ratio=float(cfg.pt_peierls_energy_ratio),
                    activation_entropy_kB=float(cfg.pt_peierls_entropy_ratio),
                    stress_ratio=float(cfg.pt_peierls_stress_ratio),
                    rate_prefactor_s=float(cfg.pt_peierls_nu0_s),
                ),
                taylor=IndependentEntropyMechanismScale(
                    energy_ratio=float(cfg.pt_taylor_energy_ratio),
                    activation_entropy_kB=float(cfg.pt_taylor_entropy_ratio),
                    stress_ratio=float(cfg.pt_taylor_stress_ratio),
                    rate_prefactor_s=float(cfg.pt_taylor_nu0_s),
                ),
                correlated_taylor=CorrelatedTaylorConfig(
                    rho_c_m2=float(cfg.pt_taylor_corr_rho_c),
                    renewal_time_s=1.0,
                    m_exponent=float(cfg.pt_taylor_m_exponent),
                    m_scale=float(cfg.pt_taylor_m_scale),
                    m_cap=float("inf"),
                ),
                peierls_stress_fraction=float(cfg.pt_peierls_stress_fraction),
                taylor_stress_fraction=float(cfg.pt_taylor_stress_fraction),
                taylor_phi_max=float("inf"),
                mobile_fraction_low_density=float(cfg.pt_mobile_fraction),
                mobile_saturation_density_m2=float("inf"),
                mobile_density_floor_m2=0.0,
                jump_fraction_of_forest_spacing=float(cfg.pt_jump_fraction),
                jump_length_min_m=0.0,
                rate_cap_s=float("inf"),
            )
        )

    @staticmethod
    def encounter_rate_s(
        peierls_rate_s: np.ndarray | float,
        jump_length_m: np.ndarray | float,
        rho_forest_m2: np.ndarray | float,
        encounter_efficiency: float,
    ) -> np.ndarray:
        """Geometric mobile-obstacle encounter rate v_P / l_f."""
        p = np.maximum(np.asarray(peierls_rate_s, dtype=float), 0.0)
        jump = np.maximum(np.asarray(jump_length_m, dtype=float), 0.0)
        rho = np.maximum(np.asarray(rho_forest_m2, dtype=float), 0.0)
        eta = max(float(encounter_efficiency), 0.0)
        return eta * jump * p * np.sqrt(rho)

    @staticmethod
    def _exchange_mobile_retained(
        mobile: np.ndarray,
        retained: np.ndarray,
        encounter_rate_s: np.ndarray,
        taylor_release_rate_s: np.ndarray,
        dt_s: float,
    ) -> tuple[np.ndarray, np.ndarray, float, float]:
        """Exact two-state exchange for each system/bin over one accepted step."""
        if dt_s <= 0.0:
            return mobile, retained, 0.0, 0.0
        ke = np.maximum(np.asarray(encounter_rate_s, dtype=float), 0.0)[None, :]
        kt = np.maximum(np.asarray(taylor_release_rate_s, dtype=float), 0.0)[None, :]
        total = np.maximum(mobile, 0.0) + np.maximum(retained, 0.0)
        rate = ke + kt
        frac_r_eq = np.divide(ke, rate, out=np.zeros_like(rate), where=rate > 0.0)
        r_eq = frac_r_eq * total
        decay = np.exp(-np.minimum(rate * float(dt_s), 700.0))
        new_r = r_eq + (np.maximum(retained, 0.0) - r_eq) * decay
        new_r = np.clip(new_r, 0.0, total)
        new_m = total - new_r
        trapped = float(np.sum(np.maximum(new_r - retained, 0.0)))
        released = float(np.sum(np.maximum(retained - new_r, 0.0)))
        return new_m, new_r, trapped, released

    def evolve(
        self,
        dt_s: float,
        T_K: float,
        stress_Pa: float,
        b: float,
        emission_hazard_integral: float = 0.0,
        system_weights=None,
    ) -> dict[str, float]:
        dt_s = max(float(dt_s), 0.0)
        emitted = self._source_commit_from_hazard(
            emission_hazard_integral, system_weights
        )
        source_recovered = self._recover_source_sites(dt_s)

        stress_profile = self.local_stress_profile_Pa(stress_Pa)
        rho_profile = self.local_forest_density_m2()
        model = self._pt_model()
        rates = model.rates(stress_profile, rho_profile, T_K, b)
        peierls_profile = np.asarray(rates["peierls_rate_s"], dtype=float).reshape(-1)
        taylor_profile = np.asarray(
            rates["taylor_completion_rate_s"], dtype=float
        ).reshape(-1)
        jump_profile = np.asarray(rates["jump_length_m"], dtype=float).reshape(-1)
        velocity_profile = jump_profile * peierls_profile
        eta = float(getattr(self.cfg, "pt_encounter_efficiency", 1.0))
        encounter_profile = self.encounter_rate_s(
            peierls_profile, jump_profile, rho_profile, eta
        )

        mobile_before = self.mobile.copy()
        retained_before = self.retained.copy()
        self.mobile, self.retained, trapped, released = (
            self._exchange_mobile_retained(
                self.mobile,
                self.retained,
                encounter_profile,
                taylor_profile,
                dt_s,
            )
        )

        retained_recovery_rate = max(
            float(getattr(self.cfg, "retained_recovery_nu0_s", 0.0)), 0.0
        )
        mobile_recovery_rate = max(
            float(getattr(self.cfg, "mobile_recovery_rate_s", 0.0)), 0.0
        )
        if dt_s > 0.0:
            fr = 1.0 - math.exp(-min(retained_recovery_rate * dt_s, 700.0))
            fm = 1.0 - math.exp(-min(mobile_recovery_rate * dt_s, 700.0))
        else:
            fr = fm = 0.0
        rec_r = self.retained * fr
        rec_m = self.mobile * fm
        self.retained -= rec_r
        self.mobile -= rec_m
        recovered = float(np.sum(rec_r) + np.sum(rec_m))

        # Conservative scalar advection speed weighted by the mobile population.
        mobile_by_bin = np.sum(np.maximum(self.mobile, 0.0), axis=0)
        if float(np.sum(mobile_by_bin)) > 0.0:
            velocity = float(
                np.sum(velocity_profile * mobile_by_bin)
                / np.sum(mobile_by_bin)
            )
            glide_rate = float(
                np.sum(peierls_profile * mobile_by_bin)
                / np.sum(mobile_by_bin)
            )
        else:
            nsrc = max(min(int(self.cfg.source_bin_count), self.n_bins), 1)
            velocity = float(np.mean(velocity_profile[:nsrc]))
            glide_rate = float(np.mean(peierls_profile[:nsrc]))
        self.mobile, escaped = self._advect_forward_field(
            self.mobile, max(velocity, 0.0) * dt_s
        )
        self.mobile = np.maximum(self.mobile, 0.0)
        self.retained = np.maximum(self.retained, 0.0)

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

        self.escaped_total += escaped
        self.recovered_total += recovered + annihilated
        self.time_s += dt_s
        series_profile = np.asarray(rates["series_rate_s"], dtype=float).reshape(-1)
        return {
            "dN_emit": float(np.sum(emitted)),
            "dN_source_recovered": float(source_recovered),
            "dN_trapped": float(trapped),
            "dN_detrapped": float(released),
            "dN_escaped": float(escaped),
            "dN_recovered": float(recovered),
            "dN_annihilated": float(annihilated),
            "mobile_before_exchange": float(np.sum(mobile_before)),
            "retained_before_exchange": float(np.sum(retained_before)),
            "glide_rate_s": float(glide_rate),
            "glide_velocity_m_s": float(velocity),
            "peierls_rate_s": float(glide_rate),
            "peierls_rate_min_s": float(np.min(peierls_profile)),
            "peierls_rate_max_s": float(np.max(peierls_profile)),
            "taylor_completion_rate_s": float(np.max(taylor_profile)),
            "taylor_completion_rate_min_s": float(np.min(taylor_profile)),
            "taylor_completion_rate_max_s": float(np.max(taylor_profile)),
            "encounter_rate_s": float(np.max(encounter_profile)),
            "encounter_rate_min_s": float(np.min(encounter_profile)),
            "encounter_rate_max_s": float(np.max(encounter_profile)),
            "pt_encounter_efficiency": float(eta),
            "peierls_taylor_series_rate_s": float(np.max(series_profile)),
            "taylor_m_eff": float(
                np.max(np.asarray(rates["taylor_m_eff"], dtype=float))
            ),
            "G_peierls_eV": float(
                np.min(np.asarray(rates["G_peierls_eV"], dtype=float))
            ),
            "G_taylor_eV": float(
                np.min(np.asarray(rates["G_taylor_eV"], dtype=float))
            ),
            "jump_length_m": float(np.max(jump_profile)),
            "rho_forest_m2": float(np.max(rho_profile)),
            "rho_forest_min_m2": float(np.min(rho_profile)),
            "rho_forest_median_m2": float(np.median(rho_profile)),
            "rho_forest_max_m2": float(np.max(rho_profile)),
            "local_stress_min_Pa": float(np.min(stress_profile)),
            "local_stress_max_Pa": float(np.max(stress_profile)),
            "transport_substeps": 1.0,
            "available_site_fraction": float(self.available_site_fraction),
            "unified_transport_retention_active": 1.0,
            "legacy_trap_barrier_active": 0.0,
            "pt_independent_entropy_active": 1.0,
        }


__all__ = ["MovingProcessZoneState"]
