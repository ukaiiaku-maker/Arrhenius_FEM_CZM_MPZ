"""Entropy-decoupled uncapped Peierls--Taylor calibration model (v9.7).

The v9.6 audit showed that the historical 0.005/0.02 energy ratios produce
nearly zero macroscopic flow stresses.  It also showed that multiplying each
mechanism's temperature slope by the active emission slope couples unrelated
thermal responses and can cancel an already-small reference barrier.

This module keeps the v9.6 uncapped, detailed-balance, emission-shaped
Peierls--Taylor kinetics, but permits the Peierls and Taylor activation
entropies to be specified independently in units of k_B:

    G0_j(T) = r_H,j G00_emit - S*_j k_B (T - Tref).

The stress-shape parameters remain inherited from the emission EXP-floor
surface.  This is a calibration model until a physically admissible parameter
family has been selected; it is not activated globally by package import.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .config import KB, EV_TO_J
from .emission_derived_plasticity import (
    CorrelatedTaylorConfig,
    EmissionDerivedPeierlsTaylorConfig,
    ExpFloorSurface,
    MechanismScale,
)
from .emission_derived_plasticity_v96 import (
    EmissionDerivedPeierlsTaylorModel as _UncappedModel,
)

KB_EV_PER_K = KB / EV_TO_J


@dataclass(frozen=True)
class IndependentEntropyMechanismScale:
    """Mechanism scale with an independently specified activation entropy.

    ``activation_entropy_kB`` is S*/k_B.  Positive S* lowers the activation
    free energy as temperature rises; negative S* raises it.  The legacy
    ``entropy_ratio`` field is retained only for serialization compatibility
    and is ignored by the v9.7 model when ``activation_entropy_kB`` is finite.
    """

    energy_ratio: float
    activation_entropy_kB: float
    stress_ratio: float = 1.0
    rate_prefactor_s: float = 1.0e12
    entropy_ratio: float = 0.0


class EmissionDerivedPeierlsTaylorModel(_UncappedModel):
    """v9.6 uncapped PT kinetics with independently tunable entropy."""

    @staticmethod
    def mechanism_gT_eV_per_K(
        parent: ExpFloorSurface,
        scale: MechanismScale | IndependentEntropyMechanismScale,
    ) -> float:
        entropy = getattr(scale, "activation_entropy_kB", None)
        if entropy is not None and np.isfinite(float(entropy)):
            return -float(entropy) * KB_EV_PER_K
        return float(getattr(scale, "entropy_ratio", 0.0)) * float(
            parent.gT_eV_per_K
        )

    @classmethod
    def _scaled_surface_values(
        cls,
        parent: ExpFloorSurface,
        scale: MechanismScale | IndependentEntropyMechanismScale,
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
        x = sigma / sigc
        return np.maximum(
            floor
            + (G0 - floor)
            * np.exp(
                -max(float(parent.a), 0.0)
                * np.power(x, max(float(parent.n), 1.0e-9))
            ),
            0.0,
        )

    def raw_zero_stress_barrier_eV(
        self, mechanism: str, T_K: float
    ) -> float:
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
            + self.mechanism_gT_eV_per_K(self.cfg.parent, scale) * dT
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
        pS = getattr(self.cfg.peierls, "activation_entropy_kB", np.nan)
        tS = getattr(self.cfg.taylor, "activation_entropy_kB", np.nan)
        out.update(
            {
                "peierls_activation_entropy_kB": np.asarray(float(pS)),
                "taylor_activation_entropy_kB": np.asarray(float(tS)),
                "peierls_gT_eV_per_K": np.asarray(
                    self.mechanism_gT_eV_per_K(
                        self.cfg.parent, self.cfg.peierls
                    )
                ),
                "taylor_gT_eV_per_K": np.asarray(
                    self.mechanism_gT_eV_per_K(
                        self.cfg.parent, self.cfg.taylor
                    )
                ),
                "entropy_decoupled_from_emission": np.asarray(True),
            }
        )
        return out


__all__ = [
    "CorrelatedTaylorConfig",
    "EmissionDerivedPeierlsTaylorConfig",
    "EmissionDerivedPeierlsTaylorModel",
    "ExpFloorSurface",
    "IndependentEntropyMechanismScale",
    "KB_EV_PER_K",
]
