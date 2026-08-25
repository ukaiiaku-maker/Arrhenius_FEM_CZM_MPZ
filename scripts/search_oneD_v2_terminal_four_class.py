#!/usr/bin/env python3
"""Search shared weak-T and ceramic rows across both predictive providers."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.stats import qmc

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_outputs" / "oneD_v2_terminal_predictive_program"
PF_MATERIALS = Path(
    "/Volumes/Data/Data/Nanopillar_calculation/"
    "PF-fracture-fatigue_v10_2_21_persistent_sites_top1/"
    "arrhenius_fracture/data/materials"
)
sys.path.insert(0, str(ROOT))

from arrhenius_fracture.emergent_gnd_contract_v913 import ACTIVE_CANDIDATE_PARAMETER_FIELDS
from scripts.run_oneD_v2_predictive_campaign import inputs, run_case, summary
from scripts.search_oneD_v2_provider_robust import source_rows

SCREEN_TEMPERATURES = (300.0, 600.0, 1000.0, 1200.0)
VALIDATION_TEMPERATURES = (300.0, 600.0, 900.0, 1000.0, 1100.0, 1200.0)
VALIDATION_SEEDS = (2008666, 2008667, 2008668)
ANCHORS = {
    "weak-T": "v913_zeroD_sobol_0129902",
    "ceramic-like": "v913_zeroD_sobol_0077080",
}
CLASS_SEEDS = {"weak-T": 2008666, "ceramic-like": 3008666}
MAXIMUM_INTERVALS = 10_000

SEARCH_SPAN_DECADES = {
    "cleave_G00_eV": 0.30,
    "cleave_gT_eV_per_K": 0.30,
    "cleave_sigc0_GPa": 0.30,
    "cleave_sT_GPa_per_K": 0.30,
    "cleave_exp_a": 0.25,
    "cleave_exp_n": 0.25,
    "emit_G00_eV": 0.35,
    "emit_gT_eV_per_K": 0.35,
    "emit_sigc0_GPa": 0.35,
    "emit_sT_GPa_per_K": 0.35,
    "emit_exp_a": 0.30,
    "emit_exp_n": 0.30,
    "peierls_H0_eV": 0.30,
    "peierls_activation_entropy_kB": 0.30,
    "taylor_H0_eV": 0.30,
    "taylor_activation_entropy_kB": 0.30,
    "rho_source0_m2": 0.80,
    "taylor_corr_rho_c_m2": 0.70,
    "taylor_corr_scale": 0.40,
    "c_blunt": 0.50,
}
OBJECTIVES = (
    "completion_objective",
    "class_trend_objective",
    "physical_topology_objective",
    "precursor_objective",
    "avalanche_size_objective",
    "onset_scale_objective",
    "provider_robustness_objective",
    "state_realism_objective",
)


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def population(material: str) -> pd.DataFrame:
    sources = source_rows().copy()
    sources["search_source"] = "ARCHIVED_V913_POOL"
    sources["search_anchor"] = sources.candidate_id
    anchor = sources[sources.candidate_id == ANCHORS[material]].iloc[0]
    fields = tuple(SEARCH_SPAN_DECADES)
    unit = qmc.Sobol(len(fields), scramble=True, seed=314159 if material == "weak-T" else 271828).random_base2(5)
    local = []
    slug = material.replace("-", "_")
    for index, vector in enumerate(unit):
        row = anchor.copy()
        row["candidate_id"] = f"oneD_v2_terminal_{slug}_{index:04d}"
        row["search_source"] = "SOBOL_LOCAL_SHARED_TWO_PROVIDER"
        row["search_anchor"] = anchor.candidate_id
        for field, coordinate in zip(fields, vector):
            value = float(anchor[field])
            factor = 10.0 ** ((2.0 * float(coordinate) - 1.0) * SEARCH_SPAN_DECADES[field])
            row[field] = value * factor
        local.append(row)
    return pd.concat([sources, pd.DataFrame(local)], ignore_index=True).drop_duplicates("candidate_id")


def case_record(row, material, temperature, provider, mechanics, drive, physics, seed, stage):
    result = run_case(
        row, material, temperature, provider, mechanics, drive, physics, 100.0,
        seed=seed, maximum_intervals=MAXIMUM_INTERVALS,
    )
    record = summary(result)
    record.update({
        "search_class": material,
        "search_stage": stage,
        "hazard_seed": int(seed),
        "search_source": row.get("search_source", ""),
        "search_anchor": row.get("search_anchor", ""),
    })
    return record


def objectives(cases: pd.DataFrame, material: str) -> dict:
    completed = cases.status.eq("TARGET_RIGHT_CENSORED")
    completion = int((~completed).sum())
    topology = precursor = largest = trend = onset_scale = state = 0.0
    curves = {}
    topology_curves = {}
    enhancement = 0.0
    for provider in ("PF", "FEMCZM"):
        q = cases[cases.provider == provider].groupby("temperature_K").mean(numeric_only=True).sort_index()
        initial = q.first_event_native_KJ_MPa_sqrt_m.to_numpy(float)
        envelope = q.onset_envelope_max_MPa_sqrt_m.to_numpy(float)
        avalanches = q.physical_avalanche_count.to_numpy(float)
        precursors = q.precursor_reinitiation_count.to_numpy(float)
        fractions = q.largest_avalanche_fraction.to_numpy(float)
        curves[provider] = initial - initial[0]
        topology_curves[provider] = avalanches
        if material == "weak-T":
            trend += np.ptp(initial) / 8.0 + np.ptp(envelope) / 12.0
            topology += float(np.mean(np.abs(avalanches - 1.0)))
            precursor += float(np.mean(precursors))
            largest += float(np.mean(np.abs(fractions - 1.0)))
            enhancement += float(np.mean(np.maximum(envelope / np.maximum(initial, 1.0e-9) - 1.5, 0.0)))
            onset_scale += float(np.mean(np.maximum(initial - 60.0, 0.0) / 30.0 + np.maximum(8.0 - initial, 0.0) / 8.0))
        else:
            # Ceramic-like: low scale, single cascade, no upper-temperature
            # enhancement, and flat-to-declining initiation.
            trend += max(initial[-1] - initial[0], 0.0) / 5.0 + np.sum(np.maximum(np.diff(initial), 0.0)) / 10.0
            topology += float(np.mean(np.abs(avalanches - 1.0)))
            precursor += float(np.mean(precursors))
            largest += float(np.mean(np.abs(fractions - 1.0)))
            enhancement += float(np.mean(np.maximum(envelope / np.maximum(initial, 1.0e-9) - 1.25, 0.0)))
            onset_scale += float(np.mean(np.maximum(initial - 35.0, 0.0) / 20.0))
        state += (
            max(float(q.max_tip_radius_um.max()) - 1000.0, 0.0) / 1000.0
            + max(float(q.max_backstress_GPa.max()) - 50.0, 0.0) / 50.0
            + max(0.001 - float(q.minimum_front_width_um.min()), 0.0) / 0.001
        )
    robustness = (
        float(np.mean(np.abs(curves["PF"] - curves["FEMCZM"]))) / 20.0
        + float(np.mean(np.abs(topology_curves["PF"] - topology_curves["FEMCZM"])))
    )
    physical_topology = topology + enhancement
    contract = bool(
        completion == 0
        and physical_topology <= (1.0 if material == "weak-T" else 0.5)
        and precursor <= (1.0 if material == "weak-T" else 0.5)
        and trend <= (3.0 if material == "weak-T" else 1.0)
        and state == 0.0
    )
    return {
        "completion_objective": completion,
        "class_trend_objective": float(trend),
        "physical_topology_objective": float(physical_topology),
        "precursor_objective": float(precursor),
        "avalanche_size_objective": float(largest),
        "onset_scale_objective": float(onset_scale),
        "provider_robustness_objective": float(robustness),
        "state_realism_objective": float(state),
        "full_class_contract_pass": contract,
    }


def nondominated(frame: pd.DataFrame) -> np.ndarray:
    values = frame.loc[:, OBJECTIVES].to_numpy(float)
    keep = np.ones(len(frame), dtype=bool)
    for i, row in enumerate(values):
        dominated = np.all(values <= row, axis=1) & np.any(values < row, axis=1)
        dominated[i] = False
        keep[i] = not bool(np.any(dominated))
    return keep


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    physics, _, providers = inputs()
    all_cases = []
    all_validated = []
    all_metrics = []
    checkpoint = OUT / "oneD_v2_terminal_search_checkpoint.parquet"
    for material in ANCHORS:
        candidates = population(material)
        screen_records = []
        for index, (_, row) in enumerate(candidates.iterrows(), start=1):
            for temperature in SCREEN_TEMPERATURES:
                for provider, (mechanics, drive) in providers.items():
                    screen_records.append(case_record(
                        row, material, temperature, provider, mechanics, drive,
                        physics, CLASS_SEEDS[material], "SCREEN",
                    ))
            if index % 4 == 0 or index == len(candidates):
                pd.DataFrame(all_cases + screen_records).to_parquet(checkpoint, index=False)
                print(f"SEARCH_PROGRESS class={material} {index}/{len(candidates)}", flush=True)
        screen = pd.DataFrame(screen_records)
        all_cases.extend(screen_records)
        by_id = candidates.set_index("candidate_id")
        metric_rows = []
        for candidate_id, local in screen.groupby("candidate_id", sort=False):
            metric_rows.append({
                "search_class": material,
                "candidate_id": candidate_id,
                "search_source": by_id.loc[candidate_id].search_source,
                "search_anchor": by_id.loc[candidate_id].search_anchor,
                **objectives(local, material),
            })
        metric = pd.DataFrame(metric_rows)
        metric["pareto_nondominated_screen"] = nondominated(metric)
        ranking = metric.sort_values(list(OBJECTIVES), kind="stable")
        shortlist = list(ranking.head(4).candidate_id)
        for candidate_id in metric[metric.pareto_nondominated_screen].candidate_id:
            if candidate_id not in shortlist and len(shortlist) < 6:
                shortlist.append(candidate_id)
        validation_records = []
        for candidate_id in shortlist:
            row = by_id.loc[candidate_id].copy()
            row["candidate_id"] = candidate_id
            seeds = VALIDATION_SEEDS if material == "weak-T" else tuple(x + 1000000 for x in VALIDATION_SEEDS)
            for seed in seeds:
                for temperature in VALIDATION_TEMPERATURES:
                    for provider, (mechanics, drive) in providers.items():
                        validation_records.append(case_record(
                            row, material, temperature, provider, mechanics, drive,
                            physics, seed, "REPEATED_VALIDATION",
                        ))
            print(f"VALIDATION_PROGRESS class={material} candidate={candidate_id}", flush=True)
        validation = pd.DataFrame(validation_records)
        all_cases.extend(validation_records)
        all_validated.extend(validation_records)
        for candidate_id, local in validation.groupby("candidate_id", sort=False):
            result = {
                "search_class": material,
                "candidate_id": candidate_id,
                "search_source": by_id.loc[candidate_id].search_source,
                "search_anchor": by_id.loc[candidate_id].search_anchor,
                **objectives(local, material),
            }
            for field in ACTIVE_CANDIDATE_PARAMETER_FIELDS:
                result[field] = float(by_id.loc[candidate_id][field])
            all_metrics.append(result)
    cases = pd.DataFrame(all_cases)
    cases.to_parquet(OUT / "oneD_v2_search_population.parquet", index=False)
    pareto = pd.DataFrame(all_metrics)
    pareto["pareto_nondominated_validation"] = False
    for material, indices in pareto.groupby("search_class").groups.items():
        pareto.loc[list(indices), "pareto_nondominated_validation"] = nondominated(pareto.loc[list(indices)])
    pareto["selection_status"] = np.where(
        pareto.full_class_contract_pass,
        "ELIGIBLE_SHARED_CLASS_ROW",
        "PARETO_DIAGNOSTIC_NOT_FULL_CONTRACT",
    )
    pareto.sort_values(["search_class", *OBJECTIVES], kind="stable").to_csv(
        OUT / "oneD_v2_pareto_candidates.csv", index=False,
    )
    checkpoint.unlink(missing_ok=True)
    manifest = {
        "schema": "oneD_v2_terminal_shared_four_class_search_v1",
        "searched_classes": list(ANCHORS),
        "retained_without_search": ["Peak", "DBTT"],
        "screen_temperatures_K": list(SCREEN_TEMPERATURES),
        "validation_temperatures_K": list(VALIDATION_TEMPERATURES),
        "validation_seed_count": len(VALIDATION_SEEDS),
        "objectives": list(OBJECTIVES),
        "objective_policy": "PARETO_PLUS_LEXICOGRAPHIC_NO_OPAQUE_WEIGHTED_SCALAR",
        "screen_case_count": int(sum(x["search_stage"] == "SCREEN" for x in all_cases)),
        "validation_case_count": int(sum(x["search_stage"] == "REPEATED_VALIDATION" for x in all_cases)),
        "eligible_rows": pareto.loc[pareto.full_class_contract_pass, ["search_class", "candidate_id"]].to_dict("records"),
        "population_sha256": _sha(OUT / "oneD_v2_search_population.parquet"),
        "pareto_sha256": _sha(OUT / "oneD_v2_pareto_candidates.csv"),
        "new_2D_PF_runs": 0,
        "new_2D_FEMCZM_runs": 0,
    }
    (OUT / "oneD_v2_parameter_search_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    print(f"SEARCH_COMPLETE cases={len(cases)} eligible={int(pareto.full_class_contract_pass.sum())}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
