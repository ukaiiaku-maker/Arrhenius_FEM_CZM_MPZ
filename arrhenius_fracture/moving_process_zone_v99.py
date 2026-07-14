"""Spatial MPZ adapter for independently calibrated PT barriers (v9.9).

The active production state remains v9.5/v9.6.  This module is used only by the
v9.9 promotion workflow so that the spatial moving-process-zone calculation can
consume the absolute Peierls/Taylor barriers, independent activation entropies,
and prefactors identified by the joint-response/continuation searches.

The v9.5 state constructs ``MechanismScale`` objects internally.  For this
adapter only, the legacy ``entropy_ratio`` slots are interpreted as activation
entropies in units of k_B and translated to the v9.7 independent-entropy model.
No constitutive caps or saturation functions are reintroduced.
"""
from __future__ import annotations

from dataclasses import replace
from typing import Any

from . import emission_derived_plasticity as _pt_base
from .emission_derived_plasticity_v97 import (
    EmissionDerivedPeierlsTaylorModel as _IndependentEntropyModel,
    IndependentEntropyMechanismScale,
)
from .moving_process_zone_v95 import MovingProcessZoneState as _SpatialStateV95


class _IndependentEntropyPTAdapter(_IndependentEntropyModel):
    """Translate v9.5 ``MechanismScale`` fields into v9.7 semantics."""

    def __init__(self, cfg: Any):
        peierls = IndependentEntropyMechanismScale(
            energy_ratio=float(cfg.peierls.energy_ratio),
            activation_entropy_kB=float(cfg.peierls.entropy_ratio),
            stress_ratio=float(cfg.peierls.stress_ratio),
            rate_prefactor_s=float(cfg.peierls.rate_prefactor_s),
        )
        taylor = IndependentEntropyMechanismScale(
            energy_ratio=float(cfg.taylor.energy_ratio),
            activation_entropy_kB=float(cfg.taylor.entropy_ratio),
            stress_ratio=float(cfg.taylor.stress_ratio),
            rate_prefactor_s=float(cfg.taylor.rate_prefactor_s),
        )
        super().__init__(replace(cfg, peierls=peierls, taylor=taylor))


class MovingProcessZoneState(_SpatialStateV95):
    """v9.5 local-density MPZ using v9.7 independent-entropy PT kinetics.

    The temporary module substitution is confined to one synchronous call to
    ``evolve`` and is restored in ``finally``.  Promotion jobs are therefore run
    in separate processes rather than multithreaded inside one interpreter.
    """

    def evolve(
        self,
        dt_s: float,
        T_K: float,
        stress_Pa: float,
        b: float,
        emission_hazard_integral: float = 0.0,
        system_weights=None,
    ) -> dict[str, float]:
        original = _pt_base.EmissionDerivedPeierlsTaylorModel
        _pt_base.EmissionDerivedPeierlsTaylorModel = _IndependentEntropyPTAdapter
        try:
            out = super().evolve(
                dt_s,
                T_K,
                stress_Pa,
                b,
                emission_hazard_integral=emission_hazard_integral,
                system_weights=system_weights,
            )
        finally:
            _pt_base.EmissionDerivedPeierlsTaylorModel = original
        out.update(
            {
                "pt_independent_entropy_active": 1.0,
                "pt_peierls_activation_entropy_kB": float(
                    self.cfg.pt_peierls_entropy_ratio
                ),
                "pt_taylor_activation_entropy_kB": float(
                    self.cfg.pt_taylor_entropy_ratio
                ),
                "pt_peierls_nu0_s": float(self.cfg.pt_peierls_nu0_s),
                "pt_taylor_nu0_s": float(self.cfg.pt_taylor_nu0_s),
            }
        )
        return out


__all__ = ["MovingProcessZoneState"]
