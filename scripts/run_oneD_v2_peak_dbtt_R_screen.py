#!/usr/bin/env python3
"""Checkpointable two-provider Peak-R/DBTT-R reduced-model screen."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.stats import qmc

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_outputs" / "oneD_v2_peak_dbtt_rcurve_search"
sys.path.insert(0, str(ROOT))

from arrhenius_fracture.emergent_gnd_campaign_v913 import candidate_from_registry_row
from arrhenius_fracture.emergent_gnd_contract_v913 import ACTIVE_CANDIDATE_PARAMETER_FIELDS
from reduced_fracture_v2.rcurve_propensity import summarize_rcurve_propensity
from scripts.run_oneD_v2_predictive_campaign import inputs, run_case, summary
from scripts.search_oneD_v2_provider_robust import source_rows


CONTROL_IDS = {
    "Peak": "v913_zeroD_sobol_0242980",
    "DBTT": "v913_zeroD_sobol_0202500",
}
CLASS_SEEDS = {"Peak": 8666, "DBTT": 1008666}
TEMPERATURES = {
    "Peak": (600.0, 900.0, 1000.0, 1200.0),
    "DBTT": (300.0, 600.0, 900.0, 1000.0, 1100.0, 1200.0),
}
SEARCH_FIELDS = (
    "cleave_G00_eV", "cleave_gT_eV_per_K", "cleave_sigc0_GPa",
    "cleave_sT_GPa_per_K", "cleave_exp_a", "cleave_exp_n",
    "emit_G00_eV", "emit_gT_eV_per_K", "emit_sigc0_GPa",
    "emit_sT_GPa_per_K", "emit_exp_a", "emit_exp_n",
    "peierls_H0_eV", "peierls_activation_entropy_kB",
    "taylor_H0_eV", "taylor_activation_entropy_kB",
    "rho_source0_m2", "taylor_corr_rho_c_m2", "taylor_corr_scale", "c_blunt",
)
SIGNED_FIELDS = {
    "cleave_gT_eV_per_K", "cleave_sT_GPa_per_K",
    "emit_gT_eV_per_K", "emit_sT_GPa_per_K",
    "peierls_activation_entropy_kB", "taylor_activation_entropy_kB",
}
LOG_SPAN_DECADES = {
    "cleave_G00_eV": 0.24, "cleave_sigc0_GPa": 0.22,
    "cleave_exp_a": 0.18, "cleave_exp_n": 0.18,
    "emit_G00_eV": 0.28, "emit_sigc0_GPa": 0.34,
    "emit_exp_a": 0.22, "emit_exp_n": 0.22,
    "peierls_H0_eV": 0.22, "taylor_H0_eV": 0.22,
    "rho_source0_m2": 0.55, "taylor_corr_rho_c_m2": 0.50,
    "taylor_corr_scale": 0.32, "c_blunt": 0.40,
}
DIRECTED_GROUPS = {
    "cleavage_barrier": ("cleave_G00_eV", "cleave_gT_eV_per_K"),
    "cleavage_stress_shape": (
        "cleave_sigc0_GPa", "cleave_sT_GPa_per_K", "cleave_exp_a", "cleave_exp_n",
    ),
    "emission_barrier": ("emit_G00_eV", "emit_gT_eV_per_K"),
    "emission_stress_shape": (
        "emit_sigc0_GPa", "emit_sT_GPa_per_K", "emit_exp_a", "emit_exp_n",
    ),
    "peierls_barrier_entropy": ("peierls_H0_eV", "peierls_activation_entropy_kB"),
    "taylor_barrier_entropy": ("taylor_H0_eV", "taylor_activation_entropy_kB"),
    "source_correlation": ("rho_source0_m2", "taylor_corr_rho_c_m2", "taylor_corr_scale"),
    "blunting": ("c_blunt",),
}
OBJECTIVES = (
    "completion_domain_objective",
    "temperature_class_objective",
    "reinitiation_objective",
    "precursor_topology_objective",
    "final_avalanche_objective",
    "control_onset_deviation_objective",
    "provider_robustness_objective",
    "process_zone_state_objective",
)


def canonical_parameter_json(row: pd.Series) -> str:
    payload = {
        field: float(row[field])
        for field in ACTIVE_CANDIDATE_PARAMETER_FIELDS
    }
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def parameter_sha256(row: pd.Series) -> str:
    return hashlib.sha256(canonical_parameter_json(row).encode("utf-8")).hexdigest()


def _valid(row: pd.Series) -> bool:
    try:
        candidate = candidate_from_registry_row(row)
        for temperature in (300.0, 1200.0):
            if candidate.cleavage.zero_stress_eV(temperature) <= 1.0e-6:
                return False
            if candidate.emission.zero_stress_eV(temperature) <= 1.0e-6:
                return False
            if candidate.cleavage.characteristic_stress_Pa(temperature) <= 1.0e6:
                return False
            if candidate.emission.characteristic_stress_Pa(temperature) <= 1.0e6:
                return False
    except (ValueError, OverflowError, FloatingPointError):
        return False
    return True


def _signed_scales(control: pd.Series, archived: pd.DataFrame) -> dict[str, float]:
    scales = {}
    for field in SIGNED_FIELDS:
        values = archived[field].astype(float).replace([np.inf, -np.inf], np.nan).dropna()
        spread = float(values.quantile(0.90) - values.quantile(0.10)) if len(values) else 0.0
        scales[field] = max(abs(float(control[field])) * 0.60, spread * 0.25, 1.0e-12)
    return scales


def _identity(row: pd.Series, material: str, method: str, anchor: str) -> pd.Series:
    result = row.copy()
    digest = parameter_sha256(result)
    if method != "HISTORICAL_CONTROL":
        result["candidate_id"] = f"oneD_v2_{material.lower()}_R_{digest[:16]}"
    result["target_response_class"] = material
    result["search_campaign_id"] = "oneD_v2_peak_dbtt_R_v1"
    result["parent_or_anchor_id"] = anchor
    result["search_generation_method"] = method
    result["active_parameter_schema"] = "V913_ACTIVE_CANDIDATE_29"
    result["active_parameter_count"] = len(ACTIVE_CANDIDATE_PARAMETER_FIELDS)
    result["canonical_parameter_json"] = canonical_parameter_json(result)
    result["parameter_sha256"] = digest
    return result


def population(material: str, sobol_power: int) -> tuple[pd.DataFrame, dict]:
    _, controls, _ = inputs()
    control = controls[controls.candidate_id == CONTROL_IDS[material]].iloc[0].copy()
    archived = source_rows().copy()
    signed_scales = _signed_scales(control, archived)
    records = [_identity(control, material, "HISTORICAL_CONTROL", str(control.candidate_id))]

    sampler = qmc.Sobol(
        len(SEARCH_FIELDS), scramble=True,
        seed=240219 if material == "Peak" else 200250,
    )
    for vector in sampler.random_base2(int(sobol_power)):
        row = control.copy()
        for field, coordinate in zip(SEARCH_FIELDS, vector):
            centered = 2.0 * float(coordinate) - 1.0
            if field in SIGNED_FIELDS:
                row[field] = float(control[field]) + centered * signed_scales[field]
            else:
                row[field] = float(control[field]) * 10.0 ** (
                    centered * LOG_SPAN_DECADES[field]
                )
        if _valid(row):
            records.append(_identity(row, material, "SOBOL_LOCAL_MATERIAL", str(control.candidate_id)))

    for group, fields in DIRECTED_GROUPS.items():
        for delta in (-0.25, -0.10, 0.10, 0.25):
            row = control.copy()
            for field in fields:
                row[field] = float(control[field]) * (1.0 + delta)
            if _valid(row):
                records.append(_identity(
                    row, material, f"DIRECTED_{group.upper()}_{delta:+.2f}", str(control.candidate_id)
                ))

    for _, endpoint in archived.iterrows():
        for alpha in (0.25, 0.50):
            row = control.copy()
            for field in SEARCH_FIELDS:
                a, b = float(control[field]), float(endpoint[field])
                if field in SIGNED_FIELDS or a <= 0.0 or b <= 0.0:
                    row[field] = (1.0 - alpha) * a + alpha * b
                else:
                    row[field] = float(np.exp((1.0 - alpha) * np.log(a) + alpha * np.log(b)))
            if _valid(row):
                records.append(_identity(
                    row, material, f"ARCHIVED_MORPH_{alpha:.2f}",
                    f"{control.candidate_id}|{endpoint.candidate_id}",
                ))
    frame = pd.DataFrame(records).drop_duplicates("parameter_sha256", keep="first")
    frame = frame.sort_values(["search_generation_method", "candidate_id"], kind="stable").reset_index(drop=True)
    bounds = {
        "schema": "oneD_v2_peak_dbtt_R_search_bounds_v1",
        "material": material,
        "control_id": CONTROL_IDS[material],
        "sobol_power": int(sobol_power),
        "candidate_count": int(len(frame)),
        "searched_fields": list(SEARCH_FIELDS),
        "fixed_active_fields": sorted(set(ACTIVE_CANDIDATE_PARAMETER_FIELDS) - set(SEARCH_FIELDS)),
        "log_span_decades": LOG_SPAN_DECADES,
        "signed_additive_half_spans": signed_scales,
        "excluded_common_physics": [
            "persistent_backstress_scale", "mpz_length_m", "source_zone_length_m",
            "shielding_orientation_factors",
        ],
        "exclusion_reason": "qualified common physics or inactive zero coordinate; not a class material row",
    }
    return frame, bounds


def _base_record(row: pd.Series, material: str, provider: str, temperature: float) -> dict:
    record = {
        "candidate_id": str(row.candidate_id),
        "target_response_class": material,
        "provider": provider,
        "temperature_K": float(temperature),
        "hazard_seed": CLASS_SEEDS[material],
        "target_um": 100.0,
        "search_stage": "FAST_TWO_PROVIDER_SCREEN",
        "search_campaign_id": row.search_campaign_id,
        "parent_or_anchor_id": row.parent_or_anchor_id,
        "search_generation_method": row.search_generation_method,
        "active_parameter_schema": row.active_parameter_schema,
        "active_parameter_count": int(row.active_parameter_count),
        "canonical_parameter_json": row.canonical_parameter_json,
        "parameter_sha256": row.parameter_sha256,
    }
    for field in ACTIVE_CANDIDATE_PARAMETER_FIELDS:
        record[field] = float(row[field])
    return record


def case_record(row, material, temperature, provider, mechanics, drive, physics):
    record = _base_record(row, material, provider, temperature)
    try:
        result = run_case(
            row, material, temperature, provider, mechanics, drive, physics, 100.0,
            seed=CLASS_SEEDS[material], maximum_intervals=10_000,
        )
        record.update(summary(result))
        record.update(summarize_rcurve_propensity(result))
        onsets = json.loads(record["onset_candidates_json"])
        record.update({
            "maximum_onset_tip_radius_um": max((x["tip_radius_m"] * 1.0e6 for x in onsets), default=np.nan),
            "maximum_onset_backstress_GPa": max((x["backstress_Pa"] * 1.0e-9 for x in onsets), default=np.nan),
            "maximum_onset_mobile_density_m2": max((x["mobile_density_m2"] for x in onsets), default=np.nan),
            "maximum_onset_retained_density_m2": max((x["retained_density_m2"] for x in onsets), default=np.nan),
            "maximum_onset_source_multiplicity": max((x["source_multiplicity"] for x in onsets), default=np.nan),
            "signed_shielding_MPa_sqrt_m": 0.0,
            "signed_shielding_policy": "INACTIVE_ACCEPTED_COMMON_PHYSICS",
            "map_oracle_fallback_count": int(result["status"] == "RIGHT_CENSORED_DRIVE_MAP_BOUND"),
            "terminal_reason": result["status"],
            "case_error": "",
        })
    except Exception as exc:  # preserve every evaluated numerical/parameter failure
        record.update({
            "status": "EVALUATION_FAILURE",
            "terminal_reason": type(exc).__name__,
            "case_error": str(exc),
            "physical_avalanche_count": 0,
            "precursor_reinitiation_count": 0,
            "N_reinit": 0,
            "largest_avalanche_fraction": 0.0,
            "map_oracle_fallback_count": 0,
        })
    return record


def _nondominated(values: np.ndarray) -> np.ndarray:
    finite = np.nan_to_num(values, nan=1.0e12, posinf=1.0e12, neginf=-1.0e12)
    keep = np.ones(len(finite), dtype=bool)
    for index, row in enumerate(finite):
        dominated = np.all(finite <= row, axis=1) & np.any(finite < row, axis=1)
        dominated[index] = False
        keep[index] = not bool(np.any(dominated))
    return keep


def _case(cases: pd.DataFrame, provider: str, temperature: float) -> pd.Series:
    local = cases[(cases.provider == provider) & (cases.temperature_K == temperature)]
    if len(local) != 1:
        raise ValueError(f"missing unique case for {provider} at {temperature:g} K")
    return local.iloc[0]


def candidate_features(cases: pd.DataFrame, material: str, control: pd.DataFrame) -> dict:
    complete = cases.status.eq("TARGET_RIGHT_CENSORED")
    completion = float((~complete).sum())
    state = 0.0
    for row in cases.itertuples():
        state += max(float(getattr(row, "max_tip_radius_um", 0.0)) - 300.0, 0.0) / 300.0
        state += max(float(getattr(row, "max_backstress_GPa", 0.0)) - 50.0, 0.0) / 50.0
        state += max(0.001 - float(getattr(row, "minimum_front_width_um", 1.0)), 0.0) / 0.001
    onset_deviation = 0.0
    provider_curve = {}
    provider_reinit = {}
    provider_counts = {}
    temperature_penalty = reinitiation = topology = final_penalty = 0.0
    for provider in ("PF", "FEMCZM"):
        local = cases[cases.provider == provider].sort_values("temperature_K")
        base = control[control.provider == provider].sort_values("temperature_K")
        initial = local.initial_onset_native_KJ_MPa_sqrt_m.to_numpy(float)
        envelope = local.maximum_onset_native_KJ_MPa_sqrt_m.to_numpy(float)
        relative = local.relative_deltaK_reinit.to_numpy(float)
        counts = local.N_reinit.to_numpy(float)
        fractions = local.largest_avalanche_fraction.to_numpy(float)
        base_initial = base.initial_onset_native_KJ_MPa_sqrt_m.to_numpy(float)
        onset_deviation += float(np.nanmean(np.abs(np.log(np.maximum(initial, 1.0e-9) / np.maximum(base_initial, 1.0e-9)))))
        provider_curve[provider] = initial / max(float(initial[0]), 1.0e-9)
        provider_reinit[provider] = relative
        provider_counts[provider] = counts
        final_penalty += float(np.nanmean(np.maximum(0.60 - fractions, 0.0)))
        if material == "Peak":
            peak_rows = local[local.temperature_K.isin((900.0, 1000.0))]
            edge_rows = local[local.temperature_K.isin((600.0, 1200.0))]
            peak_envelope = float(peak_rows.maximum_onset_native_KJ_MPa_sqrt_m.max())
            edge_envelope = float(edge_rows.maximum_onset_native_KJ_MPa_sqrt_m.max())
            prominence = (peak_envelope - edge_envelope) / max(abs(peak_envelope), 1.0e-9)
            maximum_temperature = float(local.loc[local.maximum_onset_native_KJ_MPa_sqrt_m.idxmax(), "temperature_K"])
            temperature_penalty += max(0.08 - prominence, 0.0) + (0.0 if maximum_temperature in (900.0, 1000.0) else 1.0)
            peak_count = peak_rows.N_reinit.to_numpy(float)
            peak_relative = peak_rows.relative_deltaK_reinit.to_numpy(float)
            reinitiation += float(np.nanmean(np.maximum(0.10 - peak_relative, 0.0)))
            topology += float(np.nanmean(np.maximum(1.0 - peak_count, 0.0) + np.maximum(peak_count - 3.0, 0.0)))
            topology += float(np.maximum(_case(local, provider, 600.0).N_reinit - 1.0, 0.0))
            topology += float(np.maximum(_case(local, provider, 1200.0).N_reinit - 2.0, 0.0))
        else:
            low = local[local.temperature_K.isin((300.0, 600.0))]
            upper = local[local.temperature_K.isin((1000.0, 1100.0))]
            low_envelope = low.maximum_onset_native_KJ_MPa_sqrt_m.to_numpy(float)
            upper_envelope = upper.maximum_onset_native_KJ_MPa_sqrt_m.to_numpy(float)
            shelf_variation = float(np.ptp(low_envelope) / max(float(np.mean(low_envelope)), 1.0e-9))
            transition = (float(np.mean(upper_envelope)) - float(np.mean(low_envelope))) / max(float(np.mean(low_envelope)), 1.0e-9)
            temperature_penalty += max(shelf_variation - 0.20, 0.0) + max(0.20 - transition, 0.0)
            upper_relative = upper.relative_deltaK_reinit.to_numpy(float)
            reinitiation += float(np.nanmean(np.maximum(0.15 - upper_relative, 0.0)))
            upper_count = upper.N_reinit.to_numpy(float)
            low_count = low.N_reinit.to_numpy(float)
            topology += float(np.nanmean(np.maximum(1.0 - upper_count, 0.0) + np.maximum(upper_count - 3.0, 0.0)))
            topology += float(np.nanmean(np.maximum(low_count - 1.0, 0.0)))
            high = _case(local, provider, 1200.0)
            topology += float(max(float(high.N_reinit) - 3.0, 0.0))
    robustness = (
        float(np.mean(np.abs(provider_curve["PF"] - provider_curve["FEMCZM"])))
        + float(np.mean(np.abs(provider_reinit["PF"] - provider_reinit["FEMCZM"])))
        + 0.25 * float(np.mean(np.abs(provider_counts["PF"] - provider_counts["FEMCZM"])))
    )
    return {
        "completion_domain_objective": completion,
        "temperature_class_objective": float(temperature_penalty),
        "reinitiation_objective": float(reinitiation),
        "precursor_topology_objective": float(topology),
        "final_avalanche_objective": float(final_penalty),
        "control_onset_deviation_objective": float(onset_deviation),
        "provider_robustness_objective": float(robustness),
        "process_zone_state_objective": float(state),
        "mean_relative_deltaK_reinit": float(cases.relative_deltaK_reinit.mean()),
        "maximum_relative_deltaK_reinit": float(cases.relative_deltaK_reinit.max()),
        "mean_N_reinit": float(cases.N_reinit.mean()),
        "minimum_largest_avalanche_fraction": float(cases.largest_avalanche_fraction.min()),
        "all_cases_target_complete": bool(complete.all()),
        "status_counts_json": json.dumps(cases.status.value_counts().to_dict(), sort_keys=True),
    }


def aggregate_population(cases: pd.DataFrame, candidates: pd.DataFrame, material: str) -> pd.DataFrame:
    control_id = CONTROL_IDS[material]
    control = cases[cases.candidate_id == control_id]
    records = []
    by_id = candidates.set_index("candidate_id", drop=False)
    for candidate_id, local in cases.groupby("candidate_id", sort=False):
        row = by_id.loc[candidate_id]
        item = {
            "candidate_id": candidate_id,
            "target_response_class": material,
            "parameter_sha256": row.parameter_sha256,
            "canonical_parameter_json": row.canonical_parameter_json,
            "search_generation_method": row.search_generation_method,
            "parent_or_anchor_id": row.parent_or_anchor_id,
            **candidate_features(local, material, control),
        }
        for field in ACTIVE_CANDIDATE_PARAMETER_FIELDS:
            item[field] = float(row[field])
        records.append(item)
    frame = pd.DataFrame(records)
    frame["pareto_nondominated"] = _nondominated(frame.loc[:, OBJECTIVES].to_numpy(float))
    values = frame.loc[:, OBJECTIVES].to_numpy(float)
    lo, hi = np.nanmin(values, axis=0), np.nanmax(values, axis=0)
    normalized = (values - lo) / np.maximum(hi - lo, 1.0e-12)
    front = normalized[frame.pareto_nondominated.to_numpy(bool)]
    frame["distance_to_pareto_front"] = np.min(
        np.linalg.norm(normalized[:, None, :] - front[None, :, :], axis=2), axis=1
    )
    frame["pareto_near_front"] = frame.distance_to_pareto_front <= 0.20
    frame["option_status"] = np.where(
        frame.candidate_id.eq(control_id), "HISTORICAL_CONTROL",
        np.where(frame.pareto_nondominated, "PARETO_NONDOMINATED",
                 np.where(frame.pareto_near_front, "PARETO_NEAR_FRONT", "SCREENED")),
    )
    return frame.sort_values([*OBJECTIVES, "candidate_id"], kind="stable")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--class", dest="material", choices=("Peak", "DBTT"), required=True)
    parser.add_argument("--sobol-power", type=int, default=7)
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    material = args.material
    slug = "peak" if material == "Peak" else "dbtt"
    candidates, bounds = population(material, args.sobol_power)
    candidates_path = OUT / f"oneD_v2_{slug}_R_candidates.csv"
    candidates.to_csv(candidates_path, index=False)
    (OUT / f"oneD_v2_{slug}_R_search_bounds.json").write_text(
        json.dumps(bounds, indent=2, sort_keys=True) + "\n"
    )
    checkpoint = OUT / f"oneD_v2_{slug}_R_screen_checkpoint.parquet"
    records = pd.read_parquet(checkpoint).to_dict("records") if checkpoint.exists() else []
    done = {
        (str(row["candidate_id"]), str(row["provider"]), float(row["temperature_K"]))
        for row in records
    }
    physics, _, providers = inputs()
    for index, row in candidates.iterrows():
        for provider, (mechanics, drive) in providers.items():
            for temperature in TEMPERATURES[material]:
                key = (str(row.candidate_id), provider, float(temperature))
                if key in done:
                    continue
                records.append(case_record(
                    row, material, temperature, provider, mechanics, drive, physics
                ))
                done.add(key)
        if (index + 1) % 4 == 0 or index + 1 == len(candidates):
            pd.DataFrame(records).to_parquet(checkpoint, index=False)
            print(
                f"RCURVE_SCREEN_PROGRESS class={material} candidate={index + 1}/{len(candidates)} "
                f"cases={len(records)}", flush=True,
            )
    cases = pd.DataFrame(records)
    population_path = OUT / f"oneD_v2_{slug}_R_search_population.parquet"
    cases.to_parquet(population_path, index=False)
    aggregate = aggregate_population(cases, candidates, material)
    pareto_path = OUT / f"oneD_v2_{slug}_R_pareto.csv"
    aggregate[aggregate.pareto_nondominated | aggregate.pareto_near_front].to_csv(
        pareto_path, index=False
    )
    checkpoint.unlink(missing_ok=True)
    manifest = {
        "schema": "oneD_v2_peak_dbtt_R_fast_screen_v1",
        "material": material,
        "control_id": CONTROL_IDS[material],
        "candidate_count": int(len(candidates)),
        "evaluated_case_count": int(len(cases)),
        "temperatures_K": list(TEMPERATURES[material]),
        "providers": list(providers),
        "seed": CLASS_SEEDS[material],
        "target_um": 100.0,
        "objectives": list(OBJECTIVES),
        "objective_policy": "SEPARATE_PARETO_DIMENSIONS_NO_OPAQUE_WEIGHTED_SCORE",
        "parameter_identity": "SHA256_CANONICAL_FULL_PRECISION_V913_ACTIVE_29",
        "population_sha256": hashlib.sha256(population_path.read_bytes()).hexdigest(),
        "pareto_sha256": hashlib.sha256(pareto_path.read_bytes()).hexdigest(),
        "new_2D_PF_runs": 0,
        "new_2D_FEMCZM_runs": 0,
    }
    (OUT / f"oneD_v2_{slug}_R_screen_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    print(
        f"RCURVE_SCREEN_COMPLETE class={material} candidates={len(candidates)} "
        f"pareto={int(aggregate.pareto_nondominated.sum())}", flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
