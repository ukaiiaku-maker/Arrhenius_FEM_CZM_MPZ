"""Deterministic temperature and cleavage-shelf transforms for v9.13 candidates.

The transforms operate on registry rows and never modify shared geometry,
loading, stochastic thresholds, source density, blunting, or state-translation
constants.
"""
from __future__ import annotations

import math
from typing import Any, Iterable, Mapping

from .emergent_gnd_contract_v913 import (
    ACTIVE_CANDIDATE_PARAMETER_FIELDS,
    effective_candidate_parameters,
)


def _normalized_row(row: Mapping[str, Any]) -> dict[str, Any]:
    """Return a mutable row containing every active field as a finite float."""
    out = dict(row)
    active = effective_candidate_parameters(row)
    out.update(active)
    return out


def _new_candidate_id(row: Mapping[str, Any], suffix: str) -> str:
    candidate_id = str(row.get("candidate_id", "candidate"))
    return f"{candidate_id}__{suffix}"


def surface_linear_values(
    row: Mapping[str, Any],
    prefix: str,
    temperature_K: float,
) -> tuple[float, float]:
    """Return zero-stress energy (eV) and characteristic stress (GPa)."""
    active = effective_candidate_parameters(row)
    if prefix not in ("cleave", "emit"):
        raise ValueError("prefix must be 'cleave' or 'emit'")
    temperature = float(temperature_K)
    tref = active["Tref_K"]
    energy = active[f"{prefix}_G00_eV"] + active[
        f"{prefix}_gT_eV_per_K"
    ] * (temperature - tref)
    stress = active[f"{prefix}_sigc0_GPa"] + active[
        f"{prefix}_sT_GPa_per_K"
    ] * (temperature - tref)
    return float(energy), float(stress)


def validate_positive_barrier_domain(
    row: Mapping[str, Any],
    temperatures_K: Iterable[float],
    *,
    minimum_zero_stress_energy_eV: float = 0.0,
    minimum_characteristic_stress_GPa: float = 0.0,
) -> None:
    """Reject transforms that rely on the constitutive positivity clamps."""
    minimum_energy = float(minimum_zero_stress_energy_eV)
    minimum_stress = float(minimum_characteristic_stress_GPa)
    if minimum_energy < 0.0 or minimum_stress < 0.0:
        raise ValueError("minimum barrier-domain values must be nonnegative")
    failures: list[str] = []
    for temperature in temperatures_K:
        temperature = float(temperature)
        for prefix in ("cleave", "emit"):
            energy, stress = surface_linear_values(row, prefix, temperature)
            if not math.isfinite(energy) or energy <= minimum_energy:
                failures.append(
                    f"{prefix} G0({temperature:g} K)={energy:.9g} eV"
                )
            if not math.isfinite(stress) or stress <= minimum_stress:
                failures.append(
                    f"{prefix} sigc({temperature:g} K)={stress:.9g} GPa"
                )
    if failures:
        raise ValueError(
            "transformed barrier leaves the validated positive domain: "
            + "; ".join(failures)
        )


def temperature_scale_candidate_row(
    row: Mapping[str, Any],
    scale: float,
    *,
    candidate_id: str | None = None,
) -> dict[str, Any]:
    """Translate the complete kinetic response from T to ``scale*T``.

    The exact identity applies while the fractional EXP-floor remains active.
    The absolute ``floor_min_eV`` is a shared constitutive constant rather than a
    candidate coordinate, so callers should verify event identity numerically.
    """
    lam = float(scale)
    if not math.isfinite(lam) or lam <= 0.0:
        raise ValueError("temperature scale must be positive and finite")
    out = _normalized_row(row)
    out["Tref_K"] *= lam
    for field in (
        "cleave_G00_eV",
        "emit_G00_eV",
        "peierls_H0_eV",
        "taylor_H0_eV",
    ):
        out[field] *= lam
    for field in ("cleave_sT_GPa_per_K", "emit_sT_GPa_per_K"):
        out[field] /= lam
    out["candidate_id"] = (
        str(candidate_id)
        if candidate_id is not None
        else _new_candidate_id(row, f"temperature_scale_{lam:.8g}")
    )
    out["transform_kind"] = "exact_temperature_axis_scale"
    out["temperature_scale"] = lam
    return out


def anchored_cleavage_pivot_row(
    row: Mapping[str, Any],
    *,
    shelf_temperature_K: float,
    anchor_temperature_K: float,
    shelf_energy_factor: float = 1.0,
    shelf_stress_factor: float = 1.0,
    candidate_id: str | None = None,
) -> dict[str, Any]:
    """Tilt only the linear cleavage surface about a fixed high-T anchor.

    At ``anchor_temperature_K`` the complete cleavage barrier is unchanged for
    every applied stress.  At ``shelf_temperature_K`` the zero-stress energy and
    characteristic stress are multiplied by the requested factors.  Emission,
    Peierls, Taylor, source, shielding, backstress, and geometry parameters are
    untouched.
    """
    shelf_temperature = float(shelf_temperature_K)
    anchor_temperature = float(anchor_temperature_K)
    energy_factor = float(shelf_energy_factor)
    stress_factor = float(shelf_stress_factor)
    if not math.isfinite(shelf_temperature) or not math.isfinite(anchor_temperature):
        raise ValueError("pivot temperatures must be finite")
    if shelf_temperature >= anchor_temperature:
        raise ValueError("shelf temperature must be below anchor temperature")
    if energy_factor <= 0.0 or stress_factor <= 0.0:
        raise ValueError("shelf factors must be positive")

    out = _normalized_row(row)
    tref = float(out["Tref_K"])
    energy_shelf, stress_shelf = surface_linear_values(out, "cleave", shelf_temperature)
    energy_anchor, stress_anchor = surface_linear_values(
        out, "cleave", anchor_temperature
    )
    delta_temperature = anchor_temperature - shelf_temperature

    new_gT = (energy_anchor - energy_factor * energy_shelf) / delta_temperature
    new_G00 = energy_anchor - new_gT * (anchor_temperature - tref)
    new_sT = (stress_anchor - stress_factor * stress_shelf) / delta_temperature
    new_sigc0 = stress_anchor - new_sT * (anchor_temperature - tref)

    out["cleave_G00_eV"] = float(new_G00)
    out["cleave_gT_eV_per_K"] = float(new_gT)
    out["cleave_sigc0_GPa"] = float(new_sigc0)
    out["cleave_sT_GPa_per_K"] = float(new_sT)
    out["candidate_id"] = (
        str(candidate_id)
        if candidate_id is not None
        else _new_candidate_id(
            row,
            "cleavage_pivot_"
            f"T{shelf_temperature:g}_to_T{anchor_temperature:g}_"
            f"G{energy_factor:.6g}_S{stress_factor:.6g}",
        )
    )
    out["transform_kind"] = "anchored_linear_cleavage_pivot"
    out["shelf_temperature_K"] = shelf_temperature
    out["anchor_temperature_K"] = anchor_temperature
    out["shelf_energy_factor"] = energy_factor
    out["shelf_stress_factor"] = stress_factor
    return out


def scale_cleavage_stress_axis_row(
    row: Mapping[str, Any],
    scale: float,
    *,
    candidate_id: str | None = None,
) -> dict[str, Any]:
    """Scale the cleavage characteristic-stress axis at every temperature."""
    factor = float(scale)
    if not math.isfinite(factor) or factor <= 0.0:
        raise ValueError("cleavage stress scale must be positive and finite")
    out = _normalized_row(row)
    out["cleave_sigc0_GPa"] *= factor
    out["cleave_sT_GPa_per_K"] *= factor
    out["candidate_id"] = (
        str(candidate_id)
        if candidate_id is not None
        else _new_candidate_id(row, f"cleavage_stress_scale_{factor:.8g}")
    )
    out["transform_kind"] = "global_cleavage_stress_axis_scale"
    out["cleavage_stress_scale"] = factor
    return out


def active_parameter_row(row: Mapping[str, Any]) -> dict[str, Any]:
    """Return candidate ID plus the canonical active parameter coordinates."""
    out = {"candidate_id": str(row["candidate_id"])}
    active = effective_candidate_parameters(row)
    out.update({field: active[field] for field in ACTIVE_CANDIDATE_PARAMETER_FIELDS})
    for key in (
        "transform_kind",
        "temperature_scale",
        "shelf_temperature_K",
        "anchor_temperature_K",
        "shelf_energy_factor",
        "shelf_stress_factor",
        "cleavage_stress_scale",
    ):
        if key in row:
            out[key] = row[key]
    return out


__all__ = [
    "active_parameter_row",
    "anchored_cleavage_pivot_row",
    "scale_cleavage_stress_axis_row",
    "surface_linear_values",
    "temperature_scale_candidate_row",
    "validate_positive_barrier_domain",
]
