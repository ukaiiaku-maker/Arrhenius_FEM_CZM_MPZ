"""Campaign I/O and scoring for the v9.12 emergent-GND formulation."""
from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from .emergent_gnd_state_v912_stiff import EmergentGNDState
from .emergent_gnd_types_v912 import (
    CandidateParameters,
    CommonPhysics,
    ExpFloorSurface,
    PTMechanism,
    ProtocolSegment,
    TemperatureResult,
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
    for segment in protocol:
        midpoint_K = 0.5 * (
            segment.K_start_MPa_sqrt_m + segment.K_end_MPa_sqrt_m
        )
        state.advance_time(segment.duration_s, midpoint_K, T_K)

        # Measure source depletion at the end of the physical dwell, before
        # moving-tip translation introduces fresh material into the source zone.
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
    return result


def developed_delta_K(
    result: TemperatureResult,
    window_um: tuple[float, float] = (10.0, 30.0),
) -> float:
    extension = np.asarray(result.extensions_um, dtype=float)
    values = np.asarray(result.delta_K_micro_MPa_sqrt_m, dtype=float)
    mask = (extension >= window_um[0]) & (extension <= window_um[1])
    if not np.any(mask):
        raise ValueError(f"no checkpoints inside developed window {window_um}")
    return float(np.median(values[mask]))


def score_microstructural_transition(
    temperatures_K: Sequence[float],
    developed_delta_K_MPa_sqrt_m: Sequence[float],
    *,
    min_amplitude: float = 8.0,
    target_localization: float = 0.50,
    max_width_K: float = 200.0,
) -> dict[str, float | bool]:
    """Score only Delta_K_micro(T); K0(T) is absent by construction."""
    T = np.asarray(temperatures_K, dtype=float)
    y = np.asarray(developed_delta_K_MPa_sqrt_m, dtype=float)
    order = np.argsort(T)
    T, y = T[order], y[order]
    if T.size < 3:
        raise ValueError("at least three temperatures are required")
    positive = np.maximum(np.diff(y), 0.0)
    total = float(np.sum(positive))
    localization = float(np.max(positive) / total) if total > 0.0 else 0.0
    amplitude = float(np.max(y) - np.min(y))
    if amplitude > 0.0:
        normalized = (y - np.min(y)) / amplitude
        envelope = np.maximum.accumulate(normalized)
        width = max(
            float(np.interp(0.90, envelope, T))
            - float(np.interp(0.10, envelope, T)),
            0.0,
        )
    else:
        width = float("inf")
    fit = np.polyfit(T, y, 1)
    residual = float(np.sum((y - np.polyval(fit, T)) ** 2))
    total_variance = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - residual / total_variance if total_variance > 0.0 else 1.0
    passed = (
        amplitude >= min_amplitude
        and localization >= target_localization
        and width <= max_width_K
    )
    score = (
        max(amplitude - min_amplitude, 0.0)
        + 10.0 * max(localization - target_localization, 0.0)
        + max(max_width_K - width, 0.0) / max(max_width_K, 1.0)
        + max(0.95 - r2, 0.0)
    )
    return {
        "amplitude_MPa_sqrt_m": amplitude,
        "largest_jump_localization": localization,
        "transition_width_10_90_K": width,
        "linear_r2": r2,
        "pass": bool(passed),
        "score": float(score),
    }


def candidate_from_registry_row(row: Mapping[str, Any]) -> CandidateParameters:
    Tref = float(row.get("Tref_K", 481.33))

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
        raise KeyError(
            "v9.12 requires physical rho_source0_m2; "
            "source_sites_per_system is not a constitutive substitute"
        )
    return CandidateParameters(
        candidate_id=str(row["candidate_id"]),
        cleavage=surface("cleave"),
        emission=surface("emit"),
        peierls=PTMechanism(
            float(row["peierls_H0_eV"]),
            float(row["peierls_activation_entropy_kB"]),
            float(row["peierls_exp_a"]),
            float(row["peierls_exp_n"]),
            float(row.get("peierls_nu0_s", 1.0e12)),
        ),
        taylor=PTMechanism(
            float(row["taylor_H0_eV"]),
            float(row["taylor_activation_entropy_kB"]),
            float(row["taylor_exp_a"]),
            float(row["taylor_exp_n"]),
            float(row.get("taylor_nu0_s", 1.0e11)),
        ),
        rho_source0_m2=float(rho_source),
        source_refresh_length_m=float(row["source_refresh_length_um"]) * 1.0e-6,
        taylor_corr_rho_c_m2=float(row["taylor_corr_rho_c_m2"]),
        taylor_corr_scale=float(row["taylor_corr_scale"]),
        recovery_nu0_s=float(row.get("recovery_nu0_s", 0.0) or 0.0),
        recovery_H0_eV=float(row.get("recovery_H0_eV", 0.0) or 0.0),
        recovery_activation_entropy_kB=float(
            row.get("recovery_activation_entropy_kB", 0.0) or 0.0
        ),
    )


def load_protocol_csv(path: str | Path) -> list[ProtocolSegment]:
    segments: list[ProtocolSegment] = []
    with Path(path).open(newline="") as fp:
        for row in csv.DictReader(fp):
            segments.append(
                ProtocolSegment(
                    float(row["extension_start_um"]) * 1.0e-6,
                    float(row["extension_end_um"]) * 1.0e-6,
                    float(row["K_start_MPa_sqrt_m"]),
                    float(row["K_end_MPa_sqrt_m"]),
                    float(row["duration_s"]),
                )
            )
    if not segments:
        raise RuntimeError(f"protocol is empty: {path}")
    return segments


def dump_result_json(path: str | Path, payload: Mapping[str, Any]) -> None:
    Path(path).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
