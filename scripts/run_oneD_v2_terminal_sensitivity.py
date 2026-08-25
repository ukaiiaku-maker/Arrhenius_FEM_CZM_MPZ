#!/usr/bin/env python3
"""Deterministic two-provider material and lifecycle sensitivity study."""
from __future__ import annotations

from dataclasses import replace
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_outputs" / "oneD_v2_terminal_predictive_program"
sys.path.insert(0, str(ROOT))

from reduced_fracture_v2.predictive import default_lifecycle_reduction
from scripts.run_oneD_v2_predictive_campaign import IDS, inputs, run_case, summary

TEMPERATURES = (600.0, 1000.0, 1200.0)
LEVELS = (-0.25, -0.10, 0.10, 0.25)
REFERENCE_IDS = ("v913_zeroD_sobol_0242980", "v913_zeroD_sobol_0202500")

# Material coordinates are shared across providers.  Physics coordinates are
# controlled common-law sensitivity probes, not candidate-row changes.
MATERIAL_GROUPS = {
    "cleavage_stress_free_barrier": ("cleave_G00_eV", "cleave_gT_eV_per_K"),
    "cleavage_activation_stress_shape": (
        "cleave_sigc0_GPa", "cleave_sT_GPa_per_K", "cleave_exp_a", "cleave_exp_n",
    ),
    "emission_stress_free_barrier": ("emit_G00_eV", "emit_gT_eV_per_K"),
    "emission_activation_stress_shape": (
        "emit_sigc0_GPa", "emit_sT_GPa_per_K", "emit_exp_a", "emit_exp_n",
    ),
    "peierls_barrier_entropy": ("peierls_H0_eV", "peierls_activation_entropy_kB"),
    "peierls_shape": ("peierls_exp_a", "peierls_exp_n"),
    "taylor_barrier_entropy": ("taylor_H0_eV", "taylor_activation_entropy_kB"),
    "taylor_shape": ("taylor_exp_a", "taylor_exp_n"),
    "initial_source_density": ("rho_source0_m2",),
    "taylor_correlation": ("taylor_corr_rho_c_m2", "taylor_corr_scale"),
    "blunting_coefficient": ("c_blunt",),
}
PHYSICS_GROUPS = {
    "backstress_coefficient": ("persistent_backstress_scale",),
    "process_zone_length": ("mpz_length_m",),
    "source_zone_length": ("source_zone_length_m",),
    "front_width_scale": ("reference_front_width_m",),
}
LIFECYCLE_GROUPS = {
    "lifecycle_hazard_progress": ("hazard_progress_scale", "minimum_hazard_progress_scale"),
    "lifecycle_event_length": ("nominal_advance_m",),
    "lifecycle_translation_length": ("translation_length_scale",),
    "lifecycle_reload_threshold": ("reload_gap_threshold_m",),
}
METRICS = (
    "first_event_native_KJ_MPa_sqrt_m",
    "onset_envelope_max_MPa_sqrt_m",
    "precursor_reinitiation_count",
    "physical_avalanche_count",
    "largest_avalanche_fraction",
    "mean_event_size_um",
    "max_tip_radius_um",
    "max_backstress_GPa",
    "minimum_front_width_um",
    "max_source_multiplicity",
)
SENSITIVITY_MAXIMUM_INTERVALS = 10_000


def _scaled(value, delta):
    return float(value) * (1.0 + float(delta))


def _material_row(row: pd.Series, group: str, delta: float) -> pd.Series:
    varied = row.copy()
    varied["candidate_id"] = f"{row.candidate_id}__sens__{group}__{delta:+.2f}"
    for field in MATERIAL_GROUPS[group]:
        varied[field] = _scaled(row[field], delta)
    return varied


def _physics_state(physics, group: str, delta: float):
    updates = {field: _scaled(getattr(physics, field), delta) for field in PHYSICS_GROUPS[group]}
    return replace(physics, **updates)


def _lifecycle_state(backend: str, group: str, delta: float):
    base = default_lifecycle_reduction(backend)
    updates = {}
    for field in LIFECYCLE_GROUPS[group]:
        value = float(getattr(base, field))
        # A zero lower plateau remains zero; all active coordinates scale.
        updates[field] = _scaled(value, delta) if value != 0.0 else 0.0
    return replace(base, **updates)


def _record(result, group: str, owner: str, delta: float, baseline: dict | None = None):
    record = summary(result)
    record.update({
        "sensitivity_group": group,
        "parameter_owner": owner,
        "delta_fraction": float(delta),
        "is_baseline": group == "BASELINE",
    })
    for metric in METRICS:
        base = np.nan if baseline is None else float(baseline[metric])
        value = float(record[metric])
        record[f"baseline__{metric}"] = base
        record[f"delta__{metric}"] = value - base if np.isfinite(base) else np.nan
        record[f"relative__{metric}"] = (
            (value - base) / max(abs(base), 1.0e-12) if np.isfinite(base) else np.nan
        )
    return record


def _fingerprint(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    checkpoint = OUT / "oneD_v2_parameter_sensitivity_checkpoint.parquet"
    physics, rows, providers = inputs()
    rows = rows[rows.candidate_id.isin(REFERENCE_IDS)]
    records = []
    for _, source_row in rows.iterrows():
        material = IDS[str(source_row.candidate_id)]
        for backend, (mechanics, drive) in providers.items():
            for temperature in TEMPERATURES:
                base_result = run_case(
                    source_row, material, temperature, backend, mechanics, drive,
                    physics, 100.0,
                    maximum_intervals=SENSITIVITY_MAXIMUM_INTERVALS,
                )
                base_summary = summary(base_result)
                records.append(_record(base_result, "BASELINE", "NONE", 0.0))
                for group in MATERIAL_GROUPS:
                    for delta in LEVELS:
                        result = run_case(
                            _material_row(source_row, group, delta), material,
                            temperature, backend, mechanics, drive, physics, 100.0,
                            maximum_intervals=SENSITIVITY_MAXIMUM_INTERVALS,
                        )
                        records.append(_record(result, group, "SHARED_MATERIAL_ROW", delta, base_summary))
                for group in PHYSICS_GROUPS:
                    for delta in LEVELS:
                        result = run_case(
                            source_row, material, temperature, backend, mechanics,
                            drive, _physics_state(physics, group, delta), 100.0,
                            maximum_intervals=SENSITIVITY_MAXIMUM_INTERVALS,
                        )
                        records.append(_record(result, group, "CONTROLLED_COMMON_PHYSICS", delta, base_summary))
                for group in LIFECYCLE_GROUPS:
                    for delta in LEVELS:
                        result = run_case(
                            source_row, material, temperature, backend, mechanics,
                            drive, physics, 100.0,
                            lifecycle=_lifecycle_state(backend, group, delta),
                            maximum_intervals=SENSITIVITY_MAXIMUM_INTERVALS,
                        )
                        records.append(_record(result, group, "BACKEND_REDUCTION", delta, base_summary))
                # The accepted common physics has zero direct shielding
                # orientation factors.  Record that coordinate explicitly as
                # structurally inactive rather than inventing a material law.
                for delta in LEVELS:
                    result = run_case(
                        source_row, material, temperature, backend, mechanics,
                        drive, physics, 100.0,
                        maximum_intervals=SENSITIVITY_MAXIMUM_INTERVALS,
                    )
                    records.append(_record(
                        result, "direct_shielding_coefficient",
                        "INACTIVE_ACCEPTED_COMMON_PHYSICS", delta, base_summary,
                    ))
                print(
                    f"SENSITIVITY_PROGRESS material={material} provider={backend} "
                    f"temperature_K={temperature:g}",
                    flush=True,
                )
                pd.DataFrame(records).to_parquet(checkpoint, index=False)
    frame = pd.DataFrame(records)
    path = OUT / "oneD_v2_parameter_sensitivity_results.parquet"
    frame.to_parquet(path, index=False)
    checkpoint.unlink(missing_ok=True)

    rows_out = []
    perturbed = frame[~frame.is_baseline]
    for (material, temperature, group, owner), local in perturbed.groupby(
        ["material_class", "temperature_K", "sensitivity_group", "parameter_owner"]
    ):
        item = {
            "material_class": material,
            "temperature_K": temperature,
            "sensitivity_group": group,
            "parameter_owner": owner,
        }
        for metric in METRICS:
            signs = {}
            magnitudes = {}
            for backend in ("PF", "FEMCZM"):
                q = local[local.provider == backend].set_index("delta_fraction")
                lo = float(q.loc[-0.25, metric])
                hi = float(q.loc[0.25, metric])
                signs[backend] = int(np.sign(hi - lo))
                base = float(q.iloc[0][f"baseline__{metric}"])
                magnitudes[backend] = abs(hi - lo) / max(abs(base), 1.0e-12)
            item[f"PF_sign__{metric}"] = signs["PF"]
            item[f"FEMCZM_sign__{metric}"] = signs["FEMCZM"]
            item[f"sign_agreement__{metric}"] = signs["PF"] == signs["FEMCZM"]
            item[f"mean_magnitude__{metric}"] = np.mean(list(magnitudes.values()))
        rows_out.append(item)
    summary_frame = pd.DataFrame(rows_out)
    summary_path = OUT / "oneD_v2_parameter_sensitivity_summary.csv"
    summary_frame.to_csv(summary_path, index=False)
    manifest = {
        "schema": "oneD_v2_parameter_sensitivity_v1",
        "deterministic": True,
        "materials": [IDS[x] for x in REFERENCE_IDS],
        "temperatures_K": list(TEMPERATURES),
        "providers": list(providers),
        "perturbation_levels": list(LEVELS),
        "material_groups": MATERIAL_GROUPS,
        "controlled_common_physics_groups": PHYSICS_GROUPS,
        "backend_reduction_groups": LIFECYCLE_GROUPS,
        "direct_shielding_policy": "RECORDED_INACTIVE_ZERO_ACCEPTED_PHYSICS",
        "case_count": int(len(frame)),
        "maximum_intervals_per_case": SENSITIVITY_MAXIMUM_INTERVALS,
        "results_sha256": _fingerprint(path),
        "summary_sha256": _fingerprint(summary_path),
        "new_2D_PF_runs": 0,
        "new_2D_FEMCZM_runs": 0,
    }
    (OUT / "oneD_v2_parameter_sensitivity_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    print(f"SENSITIVITY_COMPLETE cases={len(frame)} sha256={manifest['results_sha256']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
