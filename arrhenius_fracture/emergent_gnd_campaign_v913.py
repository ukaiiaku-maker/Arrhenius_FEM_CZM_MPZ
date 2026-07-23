"""Campaign utilities for the v9.13 persistent-site 1-D transfer."""
from __future__ import annotations

from typing import Any, Mapping, Sequence

from .emergent_gnd_campaign_v912 import (
    developed_delta_K,
    dump_result_json,
    load_protocol_csv,
    score_microstructural_transition,
)
from .emergent_gnd_state_v913 import EmergentGNDState
from .emergent_gnd_types_v912 import ExpFloorSurface, PTMechanism, ProtocolSegment, TemperatureResult
from .emergent_gnd_types_v913 import CandidateParameters, CommonPhysics


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

PERSISTENT_RESULT_FIELDS = (
    "persistent_site_multiplicity_per_system",
    "persistent_site_source_area_m2",
    "persistent_site_front_width_m",
    "persistent_site_width_density_m2",
    "persistent_tip_radius_m",
    "persistent_rho_back_mean_m2",
    "persistent_tau_back_mean_Pa",
    "persistent_sigma_back_mean_Pa",
    "persistent_backstress_drive_ratio_max",
    "persistent_last_source_activations",
    "persistent_last_line_content",
    "persistent_local_accumulated_slip_count",
    "tip_radius_before_advance_m",
    "tip_radius_after_advance_m",
    "tip_resharpening_by_advance_m",
)


def _float_or_default(row: Mapping[str, Any], key: str, default: float) -> float:
    value = row.get(key, default)
    if value in (None, ""):
        value = default
    return float(value)


def candidate_from_registry_row(row: Mapping[str, Any]) -> CandidateParameters:
    """Parse one top-five row under the persistent-site/no-recovery contract."""
    Tref = _float_or_default(row, "Tref_K", 481.33)

    def surface(prefix: str) -> ExpFloorSurface:
        return ExpFloorSurface(
            G00_eV=float(row[f"{prefix}_G00_eV"]),
            gT_eV_per_K=float(row[f"{prefix}_gT_eV_per_K"]),
            sigc0_Pa=float(row[f"{prefix}_sigc0_GPa"]) * 1.0e9,
            sT_Pa_per_K=float(row[f"{prefix}_sT_GPa_per_K"]) * 1.0e9,
            exp_a=float(row[f"{prefix}_exp_a"]),
            exp_n=float(row[f"{prefix}_exp_n"]),
            floor_fraction=float(row[f"{prefix}_floor_frac"]),
            Tref_K=Tref,
        )

    rho_source = row.get("rho_source0_m2")
    if rho_source in (None, ""):
        raise KeyError("v9.13 requires rho_source0_m2")
    c_blunt = row.get("c_blunt")
    if c_blunt in (None, ""):
        raise KeyError("v9.13 requires c_blunt")

    return CandidateParameters(
        candidate_id=str(row["candidate_id"]),
        cleavage=surface("cleave"),
        emission=surface("emit"),
        peierls=PTMechanism(
            float(row["peierls_H0_eV"]),
            float(row["peierls_activation_entropy_kB"]),
            float(row["peierls_exp_a"]),
            float(row["peierls_exp_n"]),
            _float_or_default(row, "peierls_nu0_s", 1.0e12),
        ),
        taylor=PTMechanism(
            float(row["taylor_H0_eV"]),
            float(row["taylor_activation_entropy_kB"]),
            float(row["taylor_exp_a"]),
            float(row["taylor_exp_n"]),
            _float_or_default(row, "taylor_nu0_s", 1.0e11),
        ),
        rho_source0_m2=float(rho_source),
        source_refresh_length_m=_float_or_default(
            row, "source_refresh_length_um", 0.0
        ) * 1.0e-6,
        taylor_corr_rho_c_m2=float(row["taylor_corr_rho_c_m2"]),
        taylor_corr_scale=float(row["taylor_corr_scale"]),
        recovery_nu0_s=0.0,
        recovery_H0_eV=_float_or_default(row, "recovery_H0_eV", 0.0),
        recovery_activation_entropy_kB=_float_or_default(
            row, "recovery_activation_entropy_kB", 0.0
        ),
        c_blunt=float(c_blunt),
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
    for field_name in ENERGY_RESULT_FIELDS + PERSISTENT_RESULT_FIELDS:
        setattr(result, field_name, [])

    for segment in protocol:
        midpoint_K = 0.5 * (
            segment.K_start_MPa_sqrt_m + segment.K_end_MPa_sqrt_m
        )
        state.advance_time(segment.duration_s, midpoint_K, T_K)
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
        result.source_available_fraction.append(1.0)
        result.source_available_fraction_pre_advance.append(
            source_fraction_pre_advance
        )
        result.pi_store_max.append(diag["pi_store_max"])
        result.pi_release_max.append(diag["pi_release_max"])
        for field_name in ENERGY_RESULT_FIELDS + PERSISTENT_RESULT_FIELDS:
            getattr(result, field_name).append(float(diag[field_name]))
    return result


__all__ = [name for name in globals() if not name.startswith("_")]
