"""Independent-shape, entropy-decoupled Peierls--Taylor kinetics (v9.10.2).

The earlier calibration models inherited the emission EXP-floor shape for the
Peierls and Taylor barriers.  This module keeps the same uncapped,
detailed-balance kinetics but permits each transport mechanism to carry its own
EXP-floor alpha and n values:

    G_j(sigma,T) = G_floor,j + (G0,j-G_floor,j)
                   exp[-alpha_j (sigma/sigma_c,j)^n_j].

Barrier enthalpies, activation entropies, alpha, and n are therefore all
independently searchable, while the optimizer continues to enforce H_P < H_T
and the full G_P(sigma,T) <= G_T(sigma,T) surface ordering.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .emission_derived_plasticity import (
    CorrelatedTaylorConfig,
    EmissionDerivedPeierlsTaylorConfig,
    ExpFloorSurface,
)
from .emission_derived_plasticity_v97 import (
    EmissionDerivedPeierlsTaylorModel as _EntropyModel,
)


@dataclass(frozen=True)
class IndependentShapeEntropyMechanismScale:
    """Mechanism scale with independent entropy and EXP-floor shape."""

    energy_ratio: float
    activation_entropy_kB: float
    exp_a: float
    exp_n: float
    stress_ratio: float = 1.0
    rate_prefactor_s: float = 1.0e12
    entropy_ratio: float = 0.0


class EmissionDerivedPeierlsTaylorModel(_EntropyModel):
    """v9.7/v9.6 kinetics with mechanism-specific alpha and n."""

    @classmethod
    def _scaled_surface_values(
        cls,
        parent: ExpFloorSurface,
        scale,
        sigma_Pa: np.ndarray | float,
        T_K: float,
    ) -> np.ndarray:
        sigma = np.maximum(np.asarray(sigma_Pa, dtype=float), 0.0)
        dT = float(T_K) - float(parent.Tref_K)
        energy_ratio = max(float(scale.energy_ratio), 0.0)
        gT = cls.mechanism_gT_eV_per_K(parent, scale)
        G0 = max(
            energy_ratio * float(parent.G00_eV) + gT * dT,
            1.0e-12,
        )
        sigc_parent = float(parent.sigc0_Pa) + float(parent.sT_Pa_per_K) * dT
        sigc = max(float(scale.stress_ratio) * sigc_parent, 1.0)
        raw_floor = max(
            float(parent.floor_min_eV) * max(energy_ratio, 1.0e-12),
            float(parent.floor_fraction) * G0,
        )
        floor = min(float(parent.floor_max_fraction) * G0, raw_floor)
        alpha = max(float(getattr(scale, "exp_a", parent.a)), 0.0)
        exponent = max(float(getattr(scale, "exp_n", parent.n)), 1.0e-9)
        x = sigma / sigc
        return np.maximum(
            floor + (G0 - floor) * np.exp(-alpha * np.power(x, exponent)),
            0.0,
        )

    def rates(
        self,
        sigma_eq_Pa: np.ndarray | float,
        rho_forest_m2: np.ndarray | float,
        T_K: float,
        b_m: float,
        rho_mobile_m2: np.ndarray | float | None = None,
    ) -> dict[str, np.ndarray]:
        out = super().rates(
            sigma_eq_Pa,
            rho_forest_m2,
            T_K,
            b_m,
            rho_mobile_m2=rho_mobile_m2,
        )
        out.update(
            {
                "peierls_exp_a": np.asarray(
                    float(getattr(self.cfg.peierls, "exp_a", self.cfg.parent.a))
                ),
                "peierls_exp_n": np.asarray(
                    float(getattr(self.cfg.peierls, "exp_n", self.cfg.parent.n))
                ),
                "taylor_exp_a": np.asarray(
                    float(getattr(self.cfg.taylor, "exp_a", self.cfg.parent.a))
                ),
                "taylor_exp_n": np.asarray(
                    float(getattr(self.cfg.taylor, "exp_n", self.cfg.parent.n))
                ),
                "independent_pt_shape_active": np.asarray(True),
            }
        )
        return out


__all__ = [
    "CorrelatedTaylorConfig",
    "EmissionDerivedPeierlsTaylorConfig",
    "EmissionDerivedPeierlsTaylorModel",
    "ExpFloorSurface",
    "IndependentShapeEntropyMechanismScale",
]
