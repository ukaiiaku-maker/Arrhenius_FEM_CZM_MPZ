#!/usr/bin/env python3
"""Export the final ceramic, weakT, and DBTT MPZ parameterizations.

The exporter reads the exact candidate rows selected after reduced spatial MPZ
validation and writes implementation-ready CSV, JSON, and Markdown records.
It separates optimized material coordinates from fixed constitutive constants
and records the exact EXP-floor equations used by v9.10.2/v9.10.3.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from arrhenius_fracture.config import EV_TO_J, KB


KB_EV_PER_K = KB / EV_TO_J

FINAL_CANDIDATES = {
    "ceramic": "ceramic_restart02_candidate00",
    "weakT": "weakT_restart00_candidate00",
    "DBTT": "DBTT_restart01_candidate05",
}

DEFAULT_MANIFESTS = {
    "ceramic": [
        Path("runs/mpz_v9_10_2_selected_materials_v1/ceramic/spatial_promotion_manifest.csv"),
        Path("runs/mpz_v9_10_2_independent_shape_global_search_v1/ceramic/spatial_promotion_manifest.csv"),
        Path("runs/mpz_v9_10_2_independent_shape_global_search_v1/ceramic/unified_global_all_candidates.csv"),
    ],
    "weakT": [
        Path("runs/mpz_v9_10_2_selected_materials_v1/weakT/spatial_promotion_manifest.csv"),
        Path("runs/mpz_v9_10_2_independent_shape_global_search_v1/weakT/spatial_promotion_manifest.csv"),
        Path("runs/mpz_v9_10_2_independent_shape_global_search_v1/weakT/unified_global_all_candidates.csv"),
    ],
    "DBTT": [
        Path("runs/mpz_v9_10_3_dbtt_1000um_manifest_v1/DBTT/spatial_promotion_manifest.csv"),
        Path("runs/mpz_v9_10_3_dbtt_targeted_global_search_v1/DBTT/spatial_promotion_manifest.csv"),
        Path("runs/mpz_v9_10_3_dbtt_targeted_global_search_v1/DBTT/unified_global_all_candidates.csv"),
    ],
}

VALIDATION_METRICS = {
    "ceramic": Path("runs/mpz_v9_10_2_ceramic_spatial_500um_v1/unified_spatial_temperature_metrics.csv"),
    "weakT": Path("runs/mpz_v9_10_2_weakT_spatial_500um_v1/unified_spatial_temperature_metrics.csv"),
    "DBTT": Path("runs/mpz_v9_10_3_dbtt_spatial_1000um_v1/unified_spatial_temperature_metrics.csv"),
}

OPTIMIZED_FIELDS = [
    "cleave_G00_eV", "cleave_gT_eV_per_K", "cleave_sigc0_GPa",
    "cleave_sT_GPa_per_K", "cleave_exp_a", "cleave_exp_n",
    "cleave_floor_frac",
    "emit_G00_eV", "emit_gT_eV_per_K", "emit_sigc0_GPa",
    "emit_sT_GPa_per_K", "emit_exp_a", "emit_exp_n",
    "emit_floor_frac",
    "peierls_H0_eV", "peierls_exp_a", "peierls_exp_n",
    "delta_H_PT_eV", "taylor_H0_eV", "taylor_exp_a",
    "taylor_exp_n", "peierls_activation_entropy_kB",
    "taylor_activation_entropy_kB", "taylor_corr_rho_c_m2",
    "taylor_corr_scale", "source_sites_per_system",
    "encounter_efficiency", "retained_recovery_rate_s",
    "source_refresh_length_um", "c_blunt",
]

LOG_DERIVATIONS = {
    "taylor_corr_rho_c_m2": "log10_taylor_corr_rho_c_m2",
    "taylor_corr_scale": "log10_taylor_corr_scale",
    "source_sites_per_system": "log10_source_sites_per_system",
    "encounter_efficiency": "log10_encounter_efficiency",
    "retained_recovery_rate_s": "log10_retained_recovery_rate_s",
    "source_refresh_length_um": "log10_source_refresh_length_um",
}

FIXED_CONSTANTS: dict[str, Any] = {
    "Tref_K": 481.33,
    "floor_min_eV_cleave_emit": 1.0e-4,
    "floor_max_fraction": 0.95,
    "cleavage_attempt_frequency_s-1": 1.0e12,
    "emission_attempt_frequency_s-1": 1.0e11,
    "peierls_attempt_frequency_s-1": 1.0e12,
    "taylor_attempt_frequency_s-1": 1.0e11,
    "cleavage_hit_order": 3.0,
    "cleavage_correlation_time_s": 1.0e-6,
    "peierls_stress_fraction": 1.0 / math.sqrt(3.0),
    "taylor_stress_fraction": 1.0 / math.sqrt(3.0),
    "taylor_m_exponent": 1.0,
    "taylor_renewal_time_s_serialized_only": 1.0,
    "n_slip_system_populations_in_reduced_MPZ": 2,
    "initial_forest_density_m-2": 5.0e12,
    "burgers_vector_m": 2.74e-10,
    "shear_modulus_Pa": 160.0e9,
    "poisson_ratio": 0.28,
    "zeroD_process_zone_radius_m": 1.0e-6,
    "zeroD_escape_length_m": 50.0e-6,
    "source_recovery_rate_s-1": 0.0,
    "mobile_recovery_rate_s-1": 0.0,
    "pair_annihilation_rate_per_count_s-1": 0.0,
    "mobile_shield_fraction": 0.0,
    "mobile_density_saturation_m-2": "infinity",
    "mobile_density_floor_m-2": 0.0,
    "jump_length_floor_m": 0.0,
    "taylor_amplification_cap": "infinity",
    "taylor_hit_order_cap": "infinity",
    "plastic_rate_cap_s-1": "infinity",
}


def finite_float(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if np.isfinite(number) else None


def locate_candidate(material_class: str, candidate_id: str, override: Path | None) -> tuple[pd.Series, Path]:
    paths = [override] if override is not None else DEFAULT_MANIFESTS[material_class]
    attempted: list[str] = []
    for path in paths:
        if path is None:
            continue
        attempted.append(str(path))
        if not path.exists():
            continue
        table = pd.read_csv(path)
        if "candidate_id" not in table.columns:
            continue
        match = table[table["candidate_id"].astype(str) == candidate_id]
        if len(match) == 1:
            return match.iloc[0].copy(), path
        if len(match) > 1:
            raise RuntimeError(f"Duplicate candidate {candidate_id!r} in {path}")
    raise FileNotFoundError(
        f"Could not find {candidate_id!r} for {material_class}. Tried: {attempted}"
    )


def normalize_candidate(material_class: str, row: pd.Series, source: Path) -> dict[str, Any]:
    record: dict[str, Any] = {
        "material_class": material_class,
        "candidate_id": str(row["candidate_id"]),
        "source_manifest": str(source),
        "calibration_version": "v9.10.3" if material_class == "DBTT" else "v9.10.2",
    }
    for field in OPTIMIZED_FIELDS:
        value = finite_float(row.get(field))
        if value is not None:
            record[field] = value

    for derived, log_name in LOG_DERIVATIONS.items():
        if derived not in record:
            log_value = finite_float(row.get(log_name))
            if log_value is not None:
                record[derived] = 10.0 ** log_value

    if "taylor_H0_eV" not in record:
        hp = record.get("peierls_H0_eV")
        delta = record.get("delta_H_PT_eV")
        if hp is not None and delta is not None:
            record["taylor_H0_eV"] = hp + delta

    required = [name for name in OPTIMIZED_FIELDS if name not in record]
    if required:
        raise RuntimeError(
            f"Candidate {record['candidate_id']} is missing implementation fields: {required}"
        )

    record["peierls_gT_eV_per_K"] = (
        -record["peierls_activation_entropy_kB"] * KB_EV_PER_K
    )
    record["taylor_gT_eV_per_K"] = (
        -record["taylor_activation_entropy_kB"] * KB_EV_PER_K
    )
    record["peierls_energy_ratio_to_emission"] = (
        record["peierls_H0_eV"] / record["emit_G00_eV"]
    )
    record["taylor_energy_ratio_to_emission"] = (
        record["taylor_H0_eV"] / record["emit_G00_eV"]
    )
    record["source_capacity_total"] = (
        FIXED_CONSTANTS["n_slip_system_populations_in_reduced_MPZ"]
        * record["source_sites_per_system"]
    )
    return record


def mechanism_rows(record: dict[str, Any]) -> list[dict[str, Any]]:
    material_class = record["material_class"]
    candidate_id = record["candidate_id"]
    rows: list[dict[str, Any]] = []

    for prefix, mechanism, nu0 in [
        ("cleave", "cleavage", 1.0e12),
        ("emit", "emission", 1.0e11),
    ]:
        rows.append({
            "material_class": material_class,
            "candidate_id": candidate_id,
            "mechanism": mechanism,
            "reference_barrier_eV": record[f"{prefix}_G00_eV"],
            "thermal_slope_eV_per_K": record[f"{prefix}_gT_eV_per_K"],
            "activation_entropy_kB": np.nan,
            "critical_stress_ref_GPa": record[f"{prefix}_sigc0_GPa"],
            "critical_stress_slope_GPa_per_K": record[f"{prefix}_sT_GPa_per_K"],
            "alpha": record[f"{prefix}_exp_a"],
            "n": record[f"{prefix}_exp_n"],
            "floor_fraction": record[f"{prefix}_floor_frac"],
            "floor_min_eV": 1.0e-4,
            "floor_max_fraction": 0.95,
            "attempt_frequency_s-1": nu0,
            "stress_driver": "crack-tip effective stress",
        })

    for prefix, mechanism, nu0 in [
        ("peierls", "Peierls", 1.0e12),
        ("taylor", "Taylor", 1.0e11),
    ]:
        energy_ratio = record[f"{prefix}_energy_ratio_to_emission"]
        rows.append({
            "material_class": material_class,
            "candidate_id": candidate_id,
            "mechanism": mechanism,
            "reference_barrier_eV": record[f"{prefix}_H0_eV"],
            "thermal_slope_eV_per_K": record[f"{prefix}_gT_eV_per_K"],
            "activation_entropy_kB": record[f"{prefix}_activation_entropy_kB"],
            "critical_stress_ref_GPa": record["emit_sigc0_GPa"],
            "critical_stress_slope_GPa_per_K": record["emit_sT_GPa_per_K"],
            "alpha": record[f"{prefix}_exp_a"],
            "n": record[f"{prefix}_exp_n"],
            "floor_fraction": record["emit_floor_frac"],
            "floor_min_eV": 1.0e-4 * energy_ratio,
            "floor_max_fraction": 0.95,
            "attempt_frequency_s-1": nu0,
            "stress_driver": (
                "resolved stress = sigma_eq/sqrt(3)"
                if mechanism == "Peierls"
                else "local Taylor stress = sigma_eq/sqrt(3)/(2*b*sqrt(rho_f))"
            ),
        })
    return rows


def validation_rows(material_class: str, candidate_id: str) -> list[dict[str, Any]]:
    path = VALIDATION_METRICS[material_class]
    if not path.exists():
        return []
    table = pd.read_csv(path)
    if "candidate_id" not in table.columns:
        return []
    selected = table[table["candidate_id"].astype(str) == candidate_id].copy()
    keep = [
        "candidate_id", "target_class", "T_K", "completed", "K_init",
        "K_plateau", "delta_KR", "final_K_shield_MPa_sqrt_m",
        "final_retained_count", "final_mobile_count",
        "final_available_site_fraction",
    ]
    return selected[[name for name in keep if name in selected.columns]].to_dict(orient="records")


def json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, list):
        return [json_safe(v) for v in value]
    if isinstance(value, (np.floating, float)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    return value


def write_notes(out: Path, records: list[dict[str, Any]]) -> None:
    candidate_lines = "\n".join(
        f"- **{r['material_class']}**: `{r['candidate_id']}` ({r['calibration_version']})"
        for r in records
    )
    text = f"""# Final MPZ barrier parameter sets

{candidate_lines}

## Barrier equation

For cleavage and emission,

```text
G0_j(T) = G00_j + gT_j (T - Tref)
sigma_c,j(T) = [sigma_c0,j + sT_j (T - Tref)] x 1e9 Pa
Gfloor,j(T) = min[0.95 G0_j, max(1e-4 eV, f_j G0_j)]
G_j(sigma,T) = Gfloor,j + [G0_j-Gfloor,j]
                 exp[-alpha_j (sigma/sigma_c,j)^n_j]
lambda_raw,j = nu0_j exp[-G_j/(k_B T)]
```

`Tref = {FIXED_CONSTANTS['Tref_K']} K`. Cleavage uses a three-hit gamma completion
rate with correlation time `1e-6 s`; emission uses the raw Arrhenius rate times
the currently available source population.

For Peierls and Taylor,

```text
G0_j(T) = H0_j - (S*_j/k_B) k_B[eV/K] (T - Tref)
```

They inherit the emission critical-stress function and floor fraction, but use
their own `H0`, activation entropy, `alpha`, and `n`. The Peierls and Taylor
attempt frequencies are `1e12 s^-1` and `1e11 s^-1`, respectively.

Taylor uses

```text
delta = 1/(2 sqrt(rho_f))
tau_T,local = (sigma_eq/sqrt(3)) delta/b
L_corr = m_scale/(2 sqrt(rho_c))
m = 1 + [2 L_corr sqrt(rho_f)]^1
lambda_T,completion = lambda_T,single / m
```

with forward-minus-reverse detailed balance. The Peierls and Taylor rates are
combined in series. No density, stress, hit-order, jump-length, mobile-density,
or plastic-rate cap is active.

## Files

- `final_mpz_parameter_sets_wide.csv`: one complete row per material class.
- `final_mpz_barrier_parameters_long.csv`: one row per class and barrier.
- `final_mpz_parameter_sets.json`: structured portable record.
- `final_mpz_fixed_constants.json`: non-optimized constitutive constants.
- `final_mpz_validation_metrics.csv`: available reduced-spatial validation data.
"""
    (out / "IMPLEMENTATION_NOTES.md").write_text(text)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ceramic-manifest", type=Path)
    parser.add_argument("--weakt-manifest", type=Path)
    parser.add_argument("--dbtt-manifest", type=Path)
    parser.add_argument(
        "--out", type=Path,
        default=Path("runs/final_mpz_material_parameters_v1"),
    )
    args = parser.parse_args()

    overrides = {
        "ceramic": args.ceramic_manifest,
        "weakT": args.weakt_manifest,
        "DBTT": args.dbtt_manifest,
    }
    records: list[dict[str, Any]] = []
    barriers: list[dict[str, Any]] = []
    validations: list[dict[str, Any]] = []

    for material_class in ("ceramic", "weakT", "DBTT"):
        candidate_id = FINAL_CANDIDATES[material_class]
        row, source = locate_candidate(
            material_class, candidate_id, overrides[material_class]
        )
        record = normalize_candidate(material_class, row, source)
        records.append(record)
        barriers.extend(mechanism_rows(record))
        validations.extend(validation_rows(material_class, candidate_id))

    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)

    pd.DataFrame(records).to_csv(
        out / "final_mpz_parameter_sets_wide.csv", index=False
    )
    pd.DataFrame(barriers).to_csv(
        out / "final_mpz_barrier_parameters_long.csv", index=False
    )
    pd.DataFrame(validations).to_csv(
        out / "final_mpz_validation_metrics.csv", index=False
    )
    (out / "final_mpz_parameter_sets.json").write_text(
        json.dumps(json_safe({r["material_class"]: r for r in records}), indent=2)
    )
    (out / "final_mpz_fixed_constants.json").write_text(
        json.dumps(json_safe(FIXED_CONSTANTS), indent=2)
    )
    write_notes(out, records)

    display = pd.DataFrame(records)[
        [
            "material_class", "candidate_id", "cleave_G00_eV",
            "emit_G00_eV", "peierls_H0_eV", "taylor_H0_eV",
            "source_sites_per_system", "encounter_efficiency",
            "retained_recovery_rate_s", "source_refresh_length_um",
            "c_blunt",
        ]
    ]
    print(display.to_string(index=False))
    print(f"\nWrote implementation-ready parameter records to: {out}")


if __name__ == "__main__":
    main()
