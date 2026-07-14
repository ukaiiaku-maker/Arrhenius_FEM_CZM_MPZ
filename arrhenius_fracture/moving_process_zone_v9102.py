"""Unified MPZ state with independent Peierls and Taylor EXP-floor shapes.

This is the v9.10 transport/retention state with the same mobile--retained
kinetics, but its PT model reads independent alpha and n values from the
promotion configuration.  It remains isolated from the package-level active
production state until spatial and 2-D validation are complete.
"""
from __future__ import annotations

from .emission_derived_plasticity import (
    CorrelatedTaylorConfig,
    EmissionDerivedPeierlsTaylorConfig,
    ExpFloorSurface,
)
from .emission_derived_plasticity_v9102 import (
    EmissionDerivedPeierlsTaylorModel,
    IndependentShapeEntropyMechanismScale,
)
from .moving_process_zone_v910 import MovingProcessZoneState as _UnifiedStateV910


class MovingProcessZoneState(_UnifiedStateV910):
    """v9.10 unified state with four independently shaped barriers."""

    state_model = "moving_pz_v9102_independent_four_barrier_shapes"

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
                peierls=IndependentShapeEntropyMechanismScale(
                    energy_ratio=float(cfg.pt_peierls_energy_ratio),
                    activation_entropy_kB=float(cfg.pt_peierls_entropy_ratio),
                    exp_a=float(
                        getattr(cfg, "pt_peierls_exp_a", cfg.pt_emit_exp_a)
                    ),
                    exp_n=float(
                        getattr(cfg, "pt_peierls_exp_n", cfg.pt_emit_exp_n)
                    ),
                    stress_ratio=float(cfg.pt_peierls_stress_ratio),
                    rate_prefactor_s=float(cfg.pt_peierls_nu0_s),
                ),
                taylor=IndependentShapeEntropyMechanismScale(
                    energy_ratio=float(cfg.pt_taylor_energy_ratio),
                    activation_entropy_kB=float(cfg.pt_taylor_entropy_ratio),
                    exp_a=float(
                        getattr(cfg, "pt_taylor_exp_a", cfg.pt_emit_exp_a)
                    ),
                    exp_n=float(
                        getattr(cfg, "pt_taylor_exp_n", cfg.pt_emit_exp_n)
                    ),
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


__all__ = ["MovingProcessZoneState"]
