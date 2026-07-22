"""Campaign wrapper that records reduced work/energy diagnostics.

All candidate parsing, scoring, and registry utilities are inherited unchanged
from ``emergent_gnd_campaign_v912``.  Only the temperature-protocol driver is
overridden so the corrected energy-bookkeeping state and its checkpoint fields
are used.
"""
from __future__ import annotations

from typing import Sequence

from .emergent_gnd_campaign_v912 import *  # noqa: F401,F403
from .emergent_gnd_state_v912_energy import EmergentGNDState
from .emergent_gnd_types_v912 import (
    CandidateParameters,
    CommonPhysics,
    ProtocolSegment,
    TemperatureResult,
)


ENERGY_RESULT_FIELDS = (
    "external_plastic_work_J_per_m",
    "nonlocal_shielding_work_J_per_m",
    "internal_stress_work_J_per_m",
    "effective_plastic_work_J_per_m",
    "effective_plastic_dissipation_J_per_m",
    "external_plastic_work_per_crack_area_J_m2",
    "effective_plastic_dissipation_per_crack_area_J_m2",
    "mobile_line_energy_J_per_m",
    "retained_line_energy_J_per_m",
    "total_line_energy_J_per_m",
)


def run_temperature_protocol(
    candidate: CandidateParameters,
    physics: CommonPhysics,
    protocol: Sequence[ProtocolSegment],
    T_K: float,
    *,
    target_cleavage_rate_s: float = 1.0e-3,
) -> TemperatureResult:
    state = EmergentGNDState(candidate, physics)
    result = TemperatureResult(candidate.candidate_id, float(T_K))
    result.numerical_integration = state.integration_metadata()
    for field_name in ENERGY_RESULT_FIELDS:
        setattr(result, field_name, [])

    for segment in protocol:
        midpoint_K = 0.5 * (
            segment.K_start_MPa_sqrt_m + segment.K_end_MPa_sqrt_m
        )
        state.advance_time(segment.duration_s, midpoint_K, T_K)

        # Measure depletion after the dwell and before moving-tip translation
        # introduces fresh source-bearing material.
        source_fraction_pre_advance = state.source_available_fraction()

        state.translate_tip(segment.da_m)
        if segment.da_m > 0.0:
            tip_speed = segment.da_m / max(segment.duration_s, 1.0e-30)
            residence = physics.mpz_length_m / max(tip_speed, 1.0e-30)
        else:
            residence = segment.duration_s
        diag = state.diagnostics(residence, midpoint_K, T_K)

        result.extensions_um.append(segment.extension_end_m * 1.0e6)
        result.K_applied_MPa_sqrt_m.append(segment.K_end_MPa_sqrt_m)
        result.delta_K_micro_MPa_sqrt_m.append(
            state.delta_K_micro_MPa_sqrt_m(T_K, target_cleavage_rate_s)
        )
        result.K_shield_MPa_sqrt_m.append(diag["K_shield_MPa_sqrt_m"])
        result.tau_gnd_tip_MPa.append(diag["tau_gnd_tip_MPa"])
        result.retained_line_count_per_unit_thickness.append(
            diag["retained_line_count_per_unit_thickness"]
        )
        result.gnd_abs_line_count_per_unit_thickness.append(
            diag["gnd_abs_line_count_per_unit_thickness"]
        )
        result.source_available_fraction.append(
            diag["source_available_fraction"]
        )
        result.source_available_fraction_pre_advance.append(
            source_fraction_pre_advance
        )
        result.pi_store_max.append(diag["pi_store_max"])
        result.pi_release_max.append(diag["pi_release_max"])
        for field_name in ENERGY_RESULT_FIELDS:
            getattr(result, field_name).append(float(diag[field_name]))
    return result


__all__ = [
    name for name in globals() if not name.startswith("_")
]
