"""Uncapped emission-derived Peierls--Taylor kinetics (MPZ v9.6).

This production closure removes the exploratory algebraic caps and saturation
functions introduced in v9.3--v9.5.  Peierls and Taylor remain scaled
EXP-floor descendants of the active crack-tip emission surface and retain exact
forward-minus-reverse detailed balance.

Density enters through the physical forest spacing

    delta = 1 / (2 sqrt(rho_f))

and through the number of correlated Taylor obstacles in a correlation length.
There is no ``phi_max``, ``m_max``, mobile-density saturation, jump-length
floor, renewal-time conversion, constitutive density cap, or plastic-rate cap.
For constant conditions the Taylor completion rate is the inverse mean gamma
waiting time, lambda_T / m.  Transient solvers may additionally integrate the
same gamma hazard explicitly without changing this mean-rate constitutive
screen.
"""
from __future__ import annotations

import numpy as np

from .emission_derived_plasticity import (
    CorrelatedTaylorConfig,
    EmissionDerivedPeierlsTaylorConfig,
    ExpFloorSurface,
    MechanismScale,
    config_from_dislocation_config,
)
from .emission_derived_plasticity_v94 import (
    EmissionDerivedPeierlsTaylorModel as _DetailedBalanceModel,
)


class EmissionDerivedPeierlsTaylorModel(_DetailedBalanceModel):
    """Detailed-balance Peierls--Taylor model with no artificial saturation."""

    def forest_spacing_m(
        self, rho_forest_m2: np.ndarray | float
    ) -> np.ndarray:
        rho = np.maximum(np.asarray(rho_forest_m2, dtype=float), 1.0e-300)
        return 1.0 / (2.0 * np.sqrt(rho))

    def taylor_amplification(
        self, rho_forest_m2: np.ndarray | float, b_m: float
    ) -> np.ndarray:
        b = max(abs(float(b_m)), 1.0e-300)
        return self.forest_spacing_m(rho_forest_m2) / b

    def taylor_local_stress(
        self,
        sigma_eq_Pa: np.ndarray | float,
        rho_forest_m2: np.ndarray | float,
        b_m: float,
    ) -> np.ndarray:
        tau = max(float(self.cfg.taylor_stress_fraction), 0.0) * np.maximum(
            np.asarray(sigma_eq_Pa, dtype=float), 0.0
        )
        return tau * self.taylor_amplification(rho_forest_m2, b_m)

    def correlation_length_m(self) -> float:
        """Physical correlation length implied by the legacy rho_c coordinate."""
        rho_c = max(
            float(self.cfg.correlated_taylor.rho_c_m2), 1.0e-300
        )
        scale = max(float(self.cfg.correlated_taylor.m_scale), 0.0)
        return scale / (2.0 * np.sqrt(rho_c))

    def natural_hit_order(
        self, rho_forest_m2: np.ndarray | float
    ) -> np.ndarray:
        """Unbounded obstacle count in one physical correlation length.

        For the default exponent p=1,

            m = 1 + 2 L_corr sqrt(rho_f).

        The exponent is retained as a correlation-geometry parameter, not as a
        saturation mechanism.  ``m_cap`` is intentionally ignored.
        """
        rho = np.maximum(np.asarray(rho_forest_m2, dtype=float), 0.0)
        p = max(float(self.cfg.correlated_taylor.m_exponent), 1.0e-12)
        obstacle_count = 2.0 * self.correlation_length_m() * np.sqrt(rho)
        return 1.0 + np.power(np.maximum(obstacle_count, 0.0), p)

    def taylor_completion_rates(
        self,
        sigma_eq_Pa: np.ndarray | float,
        rho_forest_m2: np.ndarray | float,
        T_K: float,
        b_m: float,
    ) -> dict[str, np.ndarray]:
        """Mean gamma-completion rates with no renewal window or hit cap."""
        single_net, single_forward, single_reverse = (
            self.taylor_single_hit_rates(
                sigma_eq_Pa, rho_forest_m2, T_K, b_m
            )
        )
        order = self.natural_hit_order(rho_forest_m2)
        completion_forward = single_forward / order
        completion_reverse = single_reverse / order
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

    @staticmethod
    def series_rate(
        rate_a: np.ndarray | float, rate_b: np.ndarray | float
    ) -> np.ndarray:
        a = np.maximum(np.asarray(rate_a, dtype=float), 0.0)
        b = np.maximum(np.asarray(rate_b, dtype=float), 0.0)
        out = np.zeros(np.broadcast_shapes(a.shape, b.shape), dtype=float)
        aa, bb = np.broadcast_arrays(a, b)
        active = (aa > 0.0) & (bb > 0.0)
        out[active] = aa[active] * bb[active] / (aa[active] + bb[active])
        return out

    def mobile_density(
        self,
        rho_forest_m2: np.ndarray | float,
        rho_mobile_m2: np.ndarray | float | None = None,
    ) -> np.ndarray:
        """Return an explicit mobile state or an unsaturated bulk fallback.

        The moving process zone carries mobile and retained populations
        separately and does not need the fallback.  The current bulk FEM stores
        only one density field, so until that field is split it uses the
        documented linear carrier fraction without an algebraic saturation.
        """
        if rho_mobile_m2 is not None:
            return np.maximum(np.asarray(rho_mobile_m2, dtype=float), 0.0)
        rho = np.maximum(np.asarray(rho_forest_m2, dtype=float), 0.0)
        return max(float(self.cfg.mobile_fraction_low_density), 0.0) * rho

    def jump_length(
        self, rho_forest_m2: np.ndarray | float
    ) -> np.ndarray:
        return max(
            float(self.cfg.jump_fraction_of_forest_spacing), 0.0
        ) * self.forest_spacing_m(rho_forest_m2)

    def rates(
        self,
        sigma_eq_Pa: np.ndarray | float,
        rho_forest_m2: np.ndarray | float,
        T_K: float,
        b_m: float,
        rho_mobile_m2: np.ndarray | float | None = None,
    ) -> dict[str, np.ndarray]:
        p_net, p_forward, p_reverse = self.peierls_rates(
            sigma_eq_Pa, T_K
        )
        taylor = self.taylor_completion_rates(
            sigma_eq_Pa, rho_forest_m2, T_K, b_m
        )
        series = self.series_rate(p_net, taylor["net"])
        rho_mobile = self.mobile_density(rho_forest_m2, rho_mobile_m2)
        jump = self.jump_length(rho_forest_m2)
        equivalent_rate = (
            max(float(self.cfg.equivalent_strain_factor), 0.0)
            * rho_mobile
            * abs(float(b_m))
            * jump
            * series
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
            "forest_spacing_m": np.asarray(
                self.forest_spacing_m(rho_forest_m2), dtype=float
            ),
            "taylor_amplification": np.asarray(
                self.taylor_amplification(rho_forest_m2, b_m), dtype=float
            ),
            "correlation_length_m": np.asarray(self.correlation_length_m()),
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
            "constitutive_caps_active": np.asarray(False),
        }


__all__ = [
    "CorrelatedTaylorConfig",
    "EmissionDerivedPeierlsTaylorConfig",
    "EmissionDerivedPeierlsTaylorModel",
    "ExpFloorSurface",
    "MechanismScale",
    "config_from_dislocation_config",
]
