#!/usr/bin/env python3
"""Export the final ceramic, weakT, and DBTT model parameterizations.

The exporter reads the selected candidate rows from the completed v9.10.2 and
v9.10.3 searches and writes portable CSV, JSON, and Markdown records containing:

* the four independent EXP-floor barrier surfaces;
* the Peierls/Taylor Arrhenius prefactors and entropy slopes;
* correlated Taylor and moving-process-zone state parameters;
* fixed constitutive constants required to reproduce the calibrated model;
* the reduced-spatial validation response used to select each class.

No optimization is performed.  The script is a traceability/export utility for
reimplementing the selected model in another code base.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import subprocess
from typing import Any

import numpy as np
import pandas as pd

from arrhenius_fracture.config import EV_TO_J, KB
import optimize_mpz_v9_10_2_independent_shape_global as v102


KB_EV_PER_K = KB / EV_TO_J
TREF_K = 481.33

SELECTIONS = {
    "ceramic": {
        "candidate_id": "ceramic_restart02_candidate00",
        "candidate_csv": Path(
            "runs/mpz_v9_10_2_independent_shape_global_search_v1/"
            "ceramic/unified_global_all_candidates.csv"
        ),
        "validation_csv": Path(
            "runs/mpz_v9_10_2_ceramic_spatial_500um_v1/"
            "unified_spatial_temperature_metrics.csv"
        ),
        "validation_extension_um": 500.0,
        "validation_status": "PASSED_CERAMIC_500UM_REDUCED_SPATIAL_GATE",
    },
    "weakT": {
        "candidate_id": "weakT_restart00_candidate00",
        "candidate_csv": Path(
            "runs/mpz_v9_10_2_independent_shape_global_search_v1/"
            "weakT/unified_global_all_candidates.csv"
        ),
        "validation_csv": Path(
            "runs/mpz_v9_10_2_weakT_spatial_500um_v1/"
            "unified_spatial_temperature_metrics.csv"
        ),
        "validation_extension_um": 500.0,
        "validation_status": "PASSED_WEAKT_500UM_REDUCED_SPATIAL_GATE",
    },
    "DBTT": {
        "candidate_id": "DBTT_restart01_candidate05",
        "candidate_csv": Path(
            "runs/mpz_v9_10_3_dbtt_targeted_global_search_v1/"
            "DBTT/unified_global_all_candidates.csv"
        ),
        "validation_csv": Path(
            "runs/mpz_v9_10_3_dbtt_spatial_1000um_v1/"
            "unified_spatial_temperature_metrics.csv"
        ),
        "validation_extension_um": 1000.0,
        "validation_status": "PASSED_DBTT_1000UM_TARGET_AWARE_REDUCED_SPATIAL_GATE",
    },
}


def finite_float(value: Any, name: str) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} is not numeric: {value!r}") from exc
    if not np.isfinite(out):
        raise ValueError(f"{name} is not finite: {value!r}")
    return out


def load_candidate(path: Path, candidate_id: str) -> pd.Series:
    if not path.exists():
        raise FileNotFoundError(f"Candidate file not found: {path}")
    frame = pd.read_csv(path)
    if "candidate_id" not in frame.columns:
        raise KeyError(f"candidate_id column missing from {path}")
    selected = frame[frame["candidate_id"].astype(str) == candidate_id]
    if len(selected) != 1:
        raise RuntimeError(
            f"Expected exactly one row for {candidate_id} in {path}; "
            f"found {len(selected)}"
        )
    return selected.iloc[0].copy()


def decode_candidate(row: pd.Series) -> dict[str, float]:
    missing = [name for name in v102.PARAMETER_NAMES if name not in row.index]
    if missing:
        raise KeyError(f"Candidate row is missing optimized parameters: {missing}")
    vector = np.asarray(
        [finite_float(row[name], name) for name in v102.PARAMETER_NAMES],
        dtype=float,
    )
    return v102.decode(vector)


def mechanism_barrier_rows(
    class_name: str,
    candidate_id: str,
    p: dict[str, float],
) -> list[dict[str, Any]]:
    emit0 = max(float(p["emit_G00_eV"]), 1.0e-12)

    definitions = [
        {
            "mechanism": "cleavage",
            "G_ref_eV": p["cleave_G00_eV"],
            "gT_eV_per_K": p["cleave_gT_eV_per_K"],
            "sigma_c_ref_GPa": p["cleave_sigc0_GPa"],
            "sigma_c_slope_GPa_per_K": p["cleave_sT_GPa_per_K"],
            "alpha": p["cleave_exp_a"],
            "n": p["cleave_exp_n"],
            "floor_fraction": p["cleave_floor_frac"],
            "floor_min_eV": 1.0e-4,
            "floor_max_fraction": 0.95,
            "attempt_frequency_s^-1": 1.0e12,
            "thermal_parameterization": "G0=Gref+gT*(T-Tref)",
            "stress_scale_source": "independent cleavage sigma_c(T)",
            "floor_source": "independent cleavage floor fraction",
        },
        {
            "mechanism": "emission",
            "G_ref_eV": p["emit_G00_eV"],
            "gT_eV_per_K": p["emit_gT_eV_per_K"],
            "sigma_c_ref_GPa": p["emit_sigc0_GPa"],
            "sigma_c_slope_GPa_per_K": p["emit_sT_GPa_per_K"],
            "alpha": p["emit_exp_a"],
            "n": p["emit_exp_n"],
            "floor_fraction": p["emit_floor_frac"],
            "floor_min_eV": 1.0e-4,
            "floor_max_fraction": 0.95,
            "attempt_frequency_s^-1": 1.0e11,
            "thermal_parameterization": "G0=Gref+gT*(T-Tref)",
            "stress_scale_source": "independent emission sigma_c(T)",
            "floor_source": "independent emission floor fraction",
        },
        {
            "mechanism": "Peierls",
            "G_ref_eV": p["peierls_H0_eV"],
            "gT_eV_per_K": -p["peierls_activation_entropy_kB"] * KB_EV_PER_K,
            "sigma_c_ref_GPa": p["emit_sigc0_GPa"],
            "sigma_c_slope_GPa_per_K": p["emit_sT_GPa_per_K"],
            "alpha": p["peierls_exp_a"],
            "n": p["peierls_exp_n"],
            "floor_fraction": p["emit_floor_frac"],
            "floor_min_eV": 1.0e-4 * p["peierls_H0_eV"] / emit0,
            "floor_max_fraction": 0.95,
            "attempt_frequency_s^-1": 1.0e12,
            "thermal_parameterization": "G0=H0-S*kB*(T-Tref)",
            "stress_scale_source": "inherits emission sigma_c(T); stress_ratio=1",
            "floor_source": "inherits emission floor fraction",
        },
        {
            "mechanism": "Taylor",
            "G_ref_eV": p["taylor_H0_eV"],
            "gT_eV_per_K": -p["taylor_activation_entropy_kB"] * KB_EV_PER_K,
            "sigma_c_ref_GPa": p["emit_sigc0_GPa"],
            "sigma_c_slope_GPa_per_K": p["emit_sT_GPa_per_K"],
            "alpha": p["taylor_exp_a"],
            "n": p["taylor_exp_n"],
            "floor_fraction": p["emit_floor_frac"],
            "floor_min_eV": 1.0e-4 * p["taylor_H0_eV"] / emit0,
            "floor_max_fraction": 0.95,
            "attempt_frequency_s^-1": 1.0e11,
            "thermal_parameterization": "G0=H0-S*kB*(T-Tref)",
            "stress_scale_source": "inherits emission sigma_c(T); stress_ratio=1",
            "floor_source": "inherits emission floor fraction",
        },
    ]

    rows: list[dict[str, Any]] = []
    for record in definitions:
        gT = float(record["gT_eV_per_K"])
        record.update(
            {
                "class": class_name,
                "candidate_id": candidate_id,
                "Tref_K": TREF_K,
                "equivalent_activation_entropy_kB": -gT / KB_EV_PER_K,
            }
        )
        rows.append(record)
    return rows


def state_parameter_row(
    class_name: str,
    candidate_id: str,
    p: dict[str, float],
) -> dict[str, Any]:
    return {
        "class": class_name,
        "candidate_id": candidate_id,
        "taylor_corr_rho_c_m^-2": p["taylor_corr_rho_c_m2"],
        "taylor_corr_scale": p["taylor_corr_scale"],
        "taylor_renewal_time_s": 1.0,
        "taylor_hit_order_exponent": 1.0,
        "taylor_hit_order_cap": "inf",
        "source_sites_per_system": p["source_sites_per_system"],
        "number_of_source_systems": 2,
        "source_recovery_rate_s^-1": 0.0,
        "source_refresh_length_um": p["source_refresh_length_um"],
        "encounter_efficiency": p["encounter_efficiency"],
        "retained_recovery_rate_s^-1": p["retained_recovery_rate_s"],
        "mobile_recovery_rate_s^-1": 0.0,
        "pair_annihilation_rate_per_count_s^-1": 0.0,
        "mobile_shield_fraction": 0.0,
        "c_blunt": p["c_blunt"],
        "Peierls_attempt_frequency_s^-1": 1.0e12,
        "Taylor_attempt_frequency_s^-1": 1.0e11,
        "Peierls_stress_ratio": 1.0,
        "Taylor_stress_ratio": 1.0,
        "mobile_fraction_low_density": 0.01,
        "mobile_saturation_density_m^-2": "inf",
        "mobile_density_floor_m^-2": 0.0,
        "jump_length_min_m": 0.0,
        "Taylor_stress_amplification_cap": "inf",
        "constitutive_rate_cap_s^-1": "inf",
    }


def validation_rows(path: Path, candidate_id: str, class_name: str) -> list[dict[str, Any]]:
    if not path.exists():
        raise FileNotFoundError(f"Validation file not found: {path}")
    frame = pd.read_csv(path)
    selected = frame[frame["candidate_id"].astype(str) == candidate_id].copy()
    if selected.empty:
        raise RuntimeError(f"No validation rows for {candidate_id} in {path}")
    selected["class"] = class_name
    columns = [
        "class",
        "candidate_id",
        "T_K",
        "completed",
        "K_init",
        "K_plateau",
        "delta_KR",
        "final_peierls_traverse_number",
        "final_K_shield_MPa_sqrt_m",
        "final_retained_count",
        "final_mobile_count",
    ]
    return selected[[name for name in columns if name in selected.columns]].to_dict(
        orient="records"
    )


def git_head() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True, stderr=subprocess.DEVNULL
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return "UNKNOWN"


def json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    if isinstance(value, np.generic):
        return json_safe(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return "inf" if value > 0 else "-inf"
    return value


def markdown_table(frame: pd.DataFrame, columns: list[str]) -> str:
    view = frame[columns].copy()
    headers = "| " + " | ".join(columns) + " |"
    separator = "| " + " | ".join(["---"] * len(columns)) + " |"
    lines = [headers, separator]
    for _, row in view.iterrows():
        values = []
        for value in row:
            if isinstance(value, (float, np.floating)):
                values.append(f"{float(value):.8g}")
            else:
                values.append(str(value))
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def write_markdown(
    path: Path,
    candidate_frame: pd.DataFrame,
    barrier_frame: pd.DataFrame,
    state_frame: pd.DataFrame,
    validation_frame: pd.DataFrame,
    head: str,
) -> None:
    sections = [
        "# Final three-class Arrhenius MPZ parameterizations",
        "",
        f"Repository commit: `{head}`",
        "",
        "These are the reduced-spatially validated ceramic, weak-temperature/FCC-like, and DBTT parameterizations. The DBTT candidate passed the target-aware 1000 µm reduced-spatial gate; ceramic and weakT passed their 500 µm gates.",
        "",
        "## Selected candidates",
        "",
        markdown_table(
            candidate_frame,
            [
                "class",
                "candidate_id",
                "validation_extension_um",
                "validation_status",
            ],
        ),
        "",
        "## Barrier law",
        "",
        "For each mechanism j:",
        "",
        "`G_j(sigma,T) = G_floor,j(T) + [G0,j(T)-G_floor,j(T)] exp[-alpha_j (sigma/sigma_c,j(T))^n_j]`",
        "",
        "`G_floor,j(T) = min(0.95 G0,j, max(G_floor,min,j, f_floor,j G0,j))`.",
        "",
        "Cleavage and emission use `G0(T)=Gref+gT(T-Tref)`. Peierls and Taylor use `G0(T)=H0-S*kB(T-Tref)`. Peierls and Taylor inherit the emission stress scale and floor fraction but have independent H0, entropy, alpha, and n.",
        "",
        "Arrhenius single-hit rate: `lambda_j = nu0_j exp[-G_j/(kB T)]`. Cleavage uses the correlated three-hit completion law with `m=3` and `tau=1e-6 s`. Taylor uses the correlated completion law with density-dependent hit order, `renewal_time=1 s`, and no hit-order cap.",
        "",
        "## Barrier parameters",
        "",
        markdown_table(
            barrier_frame,
            [
                "class",
                "mechanism",
                "G_ref_eV",
                "gT_eV_per_K",
                "equivalent_activation_entropy_kB",
                "sigma_c_ref_GPa",
                "sigma_c_slope_GPa_per_K",
                "alpha",
                "n",
                "floor_fraction",
                "attempt_frequency_s^-1",
            ],
        ),
        "",
        "## Process-zone and correlated-Taylor parameters",
        "",
        markdown_table(
            state_frame,
            [
                "class",
                "taylor_corr_rho_c_m^-2",
                "taylor_corr_scale",
                "source_sites_per_system",
                "source_refresh_length_um",
                "encounter_efficiency",
                "retained_recovery_rate_s^-1",
                "c_blunt",
            ],
        ),
        "",
        "## State equations",
        "",
        "`dN_m/dt = R_emit - k_enc N_m + k_T N_r - k_esc N_m - k_mrec N_m`",
        "",
        "`dN_r/dt = k_enc N_m - k_T N_r - k_rrec N_r`",
        "",
        "with `v_P=jump*lambda_P`, `k_enc=eta_enc v_P sqrt(rho_f)`, `k_T=lambda_T,completion`, and `k_esc=v_P/L_MPZ`.",
        "",
        "## Reduced-spatial validation response",
        "",
        markdown_table(
            validation_frame.sort_values(["class", "T_K"]),
            [
                "class",
                "T_K",
                "K_init",
                "K_plateau",
                "delta_KR",
            ],
        ),
        "",
        "## Fixed implementation constants",
        "",
        "- `Tref = 481.33 K`",
        "- `b = 2.74e-10 m`, `G = 160 GPa`, `nu = 0.28` in the reduced calibration closure",
        "- cleavage `nu0 = 1e12 s^-1`, emission `nu0 = 1e11 s^-1`",
        "- Peierls `nu0 = 1e12 s^-1`, Taylor `nu0 = 1e11 s^-1`",
        "- two source/slip systems; source recovery, mobile recovery, pair annihilation, and mobile shielding are zero",
        "- no mobile-density, jump-length, Taylor-amplification, hit-order, or constitutive-rate cap",
        "- spatial validation used `r_pz=1 µm`, `L_MPZ=100 µm`, and 200 MPZ bins",
        "",
        "The JSON export is the authoritative machine-readable record. The CSV files are flattened views for tables and transfer to other codes.",
    ]
    path.write_text("\n".join(sections) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("runs/final_three_class_model_parameters_v1"),
    )
    args = parser.parse_args()

    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)

    candidate_records: list[dict[str, Any]] = []
    barrier_records: list[dict[str, Any]] = []
    state_records: list[dict[str, Any]] = []
    validation_records: list[dict[str, Any]] = []
    nested: dict[str, Any] = {}

    for class_name, selection in SELECTIONS.items():
        candidate_id = str(selection["candidate_id"])
        row = load_candidate(Path(selection["candidate_csv"]), candidate_id)
        parameters = decode_candidate(row)
        barriers = mechanism_barrier_rows(class_name, candidate_id, parameters)
        state = state_parameter_row(class_name, candidate_id, parameters)
        validation = validation_rows(
            Path(selection["validation_csv"]), candidate_id, class_name
        )

        candidate_record = {
            "class": class_name,
            "candidate_id": candidate_id,
            "candidate_source": str(selection["candidate_csv"]),
            "validation_source": str(selection["validation_csv"]),
            "validation_extension_um": selection["validation_extension_um"],
            "validation_status": selection["validation_status"],
            "objective": finite_float(row.get("objective", np.nan), "objective"),
        }
        candidate_records.append(candidate_record)
        barrier_records.extend(barriers)
        state_records.append(state)
        validation_records.extend(validation)

        nested[class_name] = {
            "selection": candidate_record,
            "optimized_and_derived_parameters": parameters,
            "barriers": {record["mechanism"]: record for record in barriers},
            "state_parameters": state,
            "validation_response": validation,
        }

    head = git_head()
    fixed = {
        "repository_commit": head,
        "model_family": "v9.10.2 independent four-shape unified MPZ; v9.10.3 target-aware DBTT selection",
        "Tref_K": TREF_K,
        "KB_eV_per_K": KB_EV_PER_K,
        "b_m": 2.74e-10,
        "G_shear_Pa": 160.0e9,
        "poisson_ratio": 0.28,
        "rho0_m^-2": 5.0e12,
        "r_pz_m": 1.0e-6,
        "validation_L_MPZ_m": 100.0e-6,
        "validation_MPZ_bins": 200,
        "cleavage_attempt_frequency_s^-1": 1.0e12,
        "emission_attempt_frequency_s^-1": 1.0e11,
        "Peierls_attempt_frequency_s^-1": 1.0e12,
        "Taylor_attempt_frequency_s^-1": 1.0e11,
        "cleavage_hit_order": 3.0,
        "cleavage_correlation_time_s": 1.0e-6,
        "Taylor_renewal_time_s": 1.0,
        "Taylor_hit_order_exponent": 1.0,
        "Taylor_hit_order_cap": "inf",
        "number_of_source_systems": 2,
        "source_recovery_rate_s^-1": 0.0,
        "mobile_recovery_rate_s^-1": 0.0,
        "pair_annihilation_rate_per_count_s^-1": 0.0,
        "mobile_shield_fraction": 0.0,
        "constitutive_caps_active": False,
    }

    candidate_frame = pd.DataFrame(candidate_records)
    barrier_frame = pd.DataFrame(barrier_records)
    state_frame = pd.DataFrame(state_records)
    validation_frame = pd.DataFrame(validation_records)

    candidate_frame.to_csv(out / "selected_candidates.csv", index=False)
    barrier_frame.to_csv(out / "final_barrier_parameters_long.csv", index=False)
    state_frame.to_csv(out / "final_process_zone_parameters.csv", index=False)
    validation_frame.to_csv(out / "final_validation_response.csv", index=False)

    wide_rows = []
    for class_name, group in barrier_frame.groupby("class", sort=False):
        row: dict[str, Any] = {"class": class_name}
        for _, mechanism in group.iterrows():
            prefix = str(mechanism["mechanism"]).lower()
            for name in [
                "G_ref_eV",
                "gT_eV_per_K",
                "equivalent_activation_entropy_kB",
                "sigma_c_ref_GPa",
                "sigma_c_slope_GPa_per_K",
                "alpha",
                "n",
                "floor_fraction",
                "floor_min_eV",
                "floor_max_fraction",
                "attempt_frequency_s^-1",
            ]:
                row[f"{prefix}_{name}"] = mechanism[name]
        wide_rows.append(row)
    pd.DataFrame(wide_rows).to_csv(
        out / "final_barrier_parameters_wide.csv", index=False
    )

    payload = {
        "schema": "arrhenius_mpz_final_three_class_parameters_v1",
        "fixed_constants": fixed,
        "classes": nested,
    }
    (out / "final_three_class_model_parameters.json").write_text(
        json.dumps(json_safe(payload), indent=2)
    )

    write_markdown(
        out / "FINAL_THREE_CLASS_MODEL_PARAMETERS.md",
        candidate_frame,
        barrier_frame,
        state_frame,
        validation_frame,
        head,
    )

    print(f"Wrote final parameter package to: {out}")
    print("\nSelected candidates:")
    print(candidate_frame.to_string(index=False))
    print("\nBarrier parameters:")
    print(
        barrier_frame[
            [
                "class",
                "mechanism",
                "G_ref_eV",
                "gT_eV_per_K",
                "sigma_c_ref_GPa",
                "sigma_c_slope_GPa_per_K",
                "alpha",
                "n",
                "floor_fraction",
                "attempt_frequency_s^-1",
            ]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
