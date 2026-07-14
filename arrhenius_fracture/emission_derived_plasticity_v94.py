"""Signed detailed-balance completion of emission-derived PT kinetics.

This module preserves the v9.3 emission-derived EXP-floor surfaces and
correlated Taylor multi-hit closure, but restores the signed forward-minus-
reverse kinetic law used by the prior DDD Peierls implementation.  Net plastic
flow is therefore exactly zero at zero effective stress.
"""
from __future__ import annotations

import numpy as np
from scipy.special import gammainc

from .emission_derived_plasticity import (
    CorrelatedTaylorConfig,
    EmissionDerivedPeierlsTaylorConfig,
    EmissionDerivedPeierlsTaylorModel as _BaseModel,
    ExpFloorSurface,
    MechanismScale,
    config_from_dislocation_config,
)


class EmissionDerivedPeierlsTaylorModel(_BaseModel):
    """v9.4 PT model with zero-stress detailed balance."""

    def raw_zero_stress_barrier_eV(
        self, mechanism: str, T_K: float
    ) -> float:
        """Return the unclamped scaled zero-stress free energy."""
        key = str(mechanism).lower()
        if key in {"peierls", "p"}:
            scale = self.cfg.peierls
        elif key in {"taylor", "t"}:
            scale = self.cfg.taylor
        elif key in {"emission", "emit", "e"}:
            scale = MechanismScale(1.0, 1.0, 1.0, 1.0)
        else:
            raise ValueError(f"unknown mechanism: {mechanism}")
        dT = float(T_K) - float(self.cfg.parent.Tref_K)
        return (
            float(scale.energy_ratio) * float(self.cfg.parent.G00_eV)
            + float(scale.entropy_ratio)
            * float(self.cfg.parent.gT_eV_per_K) * dT
        )

    def peierls_rates(
        self, sigma_eq_Pa: np.ndarray | float, T_K: float
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return net, forward and zero-stress reverse Peierls rates."""
        tau = max(float(self.cfg.peierls_stress_fraction), 0.0) * np.maximum(
            np.asarray(sigma_eq_Pa, dtype=float), 0.0
        )
        forward = self._arrhenius_rate(
            self.barrier_eV("peierls", tau, T_K),
            T_K,
            self.cfg.peierls.rate_prefactor_s,
        )
        reverse_scalar = self._arrhenius_rate(
            self.barrier_eV("peierls", 0.0, T_K),
            T_K,
            self.cfg.peierls.rate_prefactor_s,
        )
        reverse = np.broadcast_to(reverse_scalar, np.shape(forward))
        return np.maximum(forward - reverse, 0.0), forward, reverse

    def peierls_rate(
        self, sigma_eq_Pa: np.ndarray | float, T_K: float
    ) -> np.ndarray:
        return self.peierls_rates(sigma_eq_Pa, T_K)[0]

    def taylor_single_hit_rates(
        self,
        sigma_eq_Pa: np.ndarray | float,
        rho_forest_m2: np.ndarray | float,
        T_K: float,
        b_m: float,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return net, forward and zero-stress reverse Taylor hit rates."""
        local_stress = self.taylor_local_stress(
            sigma_eq_Pa, rho_forest_m2, b_m
        )
        forward = self._arrhenius_rate(
            self.barrier_eV("taylor", local_stress, T_K),
            T_K,
            self.cfg.taylor.rate_prefactor_s,
        )
        reverse_scalar = self._arrhenius_rate(
            self.barrier_eV("taylor", 0.0, T_K),
            T_K,
            self.cfg.taylor.rate_prefactor_s,
        )
        reverse = np.broadcast_to(reverse_scalar, np.shape(forward))
        return np.maximum(forward - reverse, 0.0), forward, reverse

    def taylor_single_hit_rate(
        self,
        sigma_eq_Pa: np.ndarray | float,
        rho_forest_m2: np.ndarray | float,
        T_K: float,
        b_m: float,
    ) -> np.ndarray:
        return self.taylor_single_hit_rates(
            sigma_eq_Pa, rho_forest_m2, T_K, b_m
        )[0]

    def taylor_completion_rates(
        self,
        sigma_eq_Pa: np.ndarray | float,
        rho_forest_m2: np.ndarray | float,
        T_K: float,
        b_m: float,
    ) -> dict[str, np.ndarray]:
        """Apply the correlated multi-hit renewal to both directions."""
        single_net, single_forward, single_reverse = (
            self.taylor_single_hit_rates(
                sigma_eq_Pa, rho_forest_m2, T_K, b_m
            )
        )
        order = self.cfg.correlated_taylor.hit_order(rho_forest_m2)
        renewal = max(
            float(self.cfg.correlated_taylor.renewal_time_s), 1.0e-30
        )
        exposure_forward = np.minimum(
            np.maximum(single_forward * renewal, 0.0), 1.0e12
        )
        exposure_reverse = np.minimum(
            np.maximum(single_reverse * renewal, 0.0), 1.0e12
        )
        completion_forward = gammainc(order, exposure_forward) / renewal
        completion_reverse = gammainc(order, exposure_reverse) / renewal
        completion_net = np.maximum(
            completion_forward - completion_reverse, 0.0
        )
        return {
            "net": np.asarray(completion_net, dtype=float),
            "forward": np.asarray(completion_forward, dtype=float),
            "reverse": np.asarray(completion_reverse, dtype=float),
            "order": np.asarray(order, dtype=float),
            "single_net": np.asarray(single_net, dtype=float),
            "single_forward": np.asarray(single_forward, dtype=float),
            "single_reverse": np.asarray(single_reverse, dtype=float),
        }

    def taylor_completion_rate(
        self,
        sigma_eq_Pa: np.ndarray | float,
        rho_forest_m2: np.ndarray | float,
        T_K: float,
        b_m: float,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        out = self.taylor_completion_rates(
            sigma_eq_Pa, rho_forest_m2, T_K, b_m
        )
        return out["net"], out["order"], out["single_net"]

    def rates(
        self,
        sigma_eq_Pa: np.ndarray | float,
        rho_forest_m2: np.ndarray | float,
        T_K: float,
        b_m: float,
    ) -> dict[str, np.ndarray]:
        p_net, p_forward, p_reverse = self.peierls_rates(
            sigma_eq_Pa, T_K
        )
        taylor = self.taylor_completion_rates(
            sigma_eq_Pa, rho_forest_m2, T_K, b_m
        )
        series = self.series_rate(p_net, taylor["net"])
        rho_mobile = self.mobile_density(rho_forest_m2)
        jump = self.jump_length(rho_forest_m2)
        equivalent_rate = (
            max(float(self.cfg.equivalent_strain_factor), 0.0)
            * rho_mobile
            * abs(float(b_m))
            * jump
            * series
        )
        cap = float(self.cfg.rate_cap_s)
        if np.isfinite(cap):
            equivalent_rate = np.minimum(
                equivalent_rate, max(cap, 0.0)
            )
        local_stress = self.taylor_local_stress(
            sigma_eq_Pa, rho_forest_m2, b_m
        )
        return {
            "peierls_rate_s": np.asarray(p_net, dtype=float),
            "peierls_forward_rate_s": np.asarray(p_forward, dtype=float),
            "peierls_reverse_rate_s": np.asarray(p_reverse, dtype=float),
            "taylor_single_hit_rate_s": taylor["single_net"],
            "taylor_single_hit_forward_rate_s": taylor["single_forward"],
            "taylor_single_hit_reverse_rate_s": taylor["single_reverse"],
            "taylor_completion_rate_s": taylor["net"],
            "taylor_completion_forward_rate_s": taylor["forward"],
            "taylor_completion_reverse_rate_s": taylor["reverse"],
            "series_rate_s": np.asarray(series, dtype=float),
            "equivalent_plastic_rate_s": np.asarray(
                equivalent_rate, dtype=float
            ),
            "taylor_m_eff": taylor["order"],
            "rho_mobile_m2": np.asarray(rho_mobile, dtype=float),
            "jump_length_m": np.asarray(jump, dtype=float),
            "taylor_local_stress_Pa": np.asarray(local_stress, dtype=float),
            "G_peierls_eV": np.asarray(
                self.barrier_eV(
                    "peierls",
                    max(float(self.cfg.peierls_stress_fraction), 0.0)
                    * np.maximum(
                        np.asarray(sigma_eq_Pa, dtype=float), 0.0
                    ),
                    T_K,
                ),
                dtype=float,
            ),
            "G_taylor_eV": np.asarray(
                self.barrier_eV("taylor", local_stress, T_K), dtype=float
            ),
            "raw_G0_peierls_eV": np.asarray(
                self.raw_zero_stress_barrier_eV("peierls", T_K)
            ),
            "raw_G0_taylor_eV": np.asarray(
                self.raw_zero_stress_barrier_eV("taylor", T_K)
            ),
        }


__all__ = [
    "CorrelatedTaylorConfig",
    "EmissionDerivedPeierlsTaylorConfig",
    "EmissionDerivedPeierlsTaylorModel",
    "ExpFloorSurface",
    "MechanismScale",
    "config_from_dislocation_config",
]
