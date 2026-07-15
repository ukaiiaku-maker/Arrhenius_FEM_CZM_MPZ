"""Activate the independently shaped v9.10.2 Peierls/Taylor law in 2-D bulk FEM.

The legacy 2-D path constructed a base emission-derived model that inherited the
emission alpha/n for Peierls and Taylor. This context manager replaces only the
two symbols imported inside ``plasticity.update_plasticity`` so the bulk solver
uses the exact independent H0, S*, alpha, and n values selected by v9.10.2/3.
"""
from __future__ import annotations

from contextlib import contextmanager
from typing import Any

from .emission_derived_plasticity import (
    CorrelatedTaylorConfig,
    EmissionDerivedPeierlsTaylorConfig,
    ExpFloorSurface,
)
from .emission_derived_plasticity_v9102 import (
    EmissionDerivedPeierlsTaylorModel,
    IndependentShapeEntropyMechanismScale,
)


def independent_config_from_dislocation_config(disl_cfg: Any) -> EmissionDerivedPeierlsTaylorConfig:
    parent = ExpFloorSurface(
        G00_eV=float(disl_cfg.pt_emit_G00_eV),
        gT_eV_per_K=float(disl_cfg.pt_emit_gT_eV_per_K),
        sigc0_Pa=float(disl_cfg.pt_emit_sigc0_Pa),
        sT_Pa_per_K=float(disl_cfg.pt_emit_sT_Pa_per_K),
        Tref_K=float(disl_cfg.pt_emit_Tref_K),
        a=float(disl_cfg.pt_emit_exp_a),
        n=float(disl_cfg.pt_emit_exp_n),
        floor_fraction=float(disl_cfg.pt_emit_floor_frac),
        floor_min_eV=float(disl_cfg.pt_emit_floor_min_eV),
        floor_max_fraction=float(disl_cfg.pt_emit_floor_max_frac),
    )
    return EmissionDerivedPeierlsTaylorConfig(
        parent=parent,
        peierls=IndependentShapeEntropyMechanismScale(
            energy_ratio=float(disl_cfg.pt_peierls_energy_ratio),
            activation_entropy_kB=float(getattr(
                disl_cfg, "pt_peierls_activation_entropy_kB",
                disl_cfg.pt_peierls_entropy_ratio,
            )),
            exp_a=float(disl_cfg.pt_peierls_exp_a),
            exp_n=float(disl_cfg.pt_peierls_exp_n),
            stress_ratio=float(disl_cfg.pt_peierls_stress_ratio),
            rate_prefactor_s=float(disl_cfg.pt_peierls_nu0_s),
        ),
        taylor=IndependentShapeEntropyMechanismScale(
            energy_ratio=float(disl_cfg.pt_taylor_energy_ratio),
            activation_entropy_kB=float(getattr(
                disl_cfg, "pt_taylor_activation_entropy_kB",
                disl_cfg.pt_taylor_entropy_ratio,
            )),
            exp_a=float(disl_cfg.pt_taylor_exp_a),
            exp_n=float(disl_cfg.pt_taylor_exp_n),
            stress_ratio=float(disl_cfg.pt_taylor_stress_ratio),
            rate_prefactor_s=float(disl_cfg.pt_taylor_nu0_s),
        ),
        correlated_taylor=CorrelatedTaylorConfig(
            rho_c_m2=float(disl_cfg.pt_taylor_corr_rho_c),
            renewal_time_s=1.0,
            m_exponent=1.0,
            m_scale=float(disl_cfg.pt_taylor_m_scale),
            m_cap=float("inf"),
        ),
        peierls_stress_fraction=float(disl_cfg.pt_peierls_stress_fraction),
        taylor_stress_fraction=float(disl_cfg.pt_taylor_stress_fraction),
        taylor_phi_max=float("inf"),
        mobile_fraction_low_density=float(disl_cfg.pt_mobile_fraction),
        mobile_saturation_density_m2=float("inf"),
        mobile_density_floor_m2=0.0,
        jump_fraction_of_forest_spacing=float(disl_cfg.pt_jump_fraction),
        jump_length_min_m=0.0,
        rate_cap_s=float("inf"),
    )


@contextmanager
def independent_bulk_pt_active():
    """Temporarily make the 2-D bulk plasticity import the v9.10.2 model."""
    from . import emission_derived_plasticity as base

    old_model = base.EmissionDerivedPeierlsTaylorModel
    old_builder = base.config_from_dislocation_config
    base.EmissionDerivedPeierlsTaylorModel = EmissionDerivedPeierlsTaylorModel
    base.config_from_dislocation_config = independent_config_from_dislocation_config
    try:
        yield
    finally:
        base.EmissionDerivedPeierlsTaylorModel = old_model
        base.config_from_dislocation_config = old_builder


__all__ = [
    "EmissionDerivedPeierlsTaylorModel",
    "independent_config_from_dislocation_config",
    "independent_bulk_pt_active",
]
