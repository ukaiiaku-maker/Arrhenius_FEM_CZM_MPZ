#!/usr/bin/env python3
"""Focused Pareto refinement between topology-good and trend-good anchors."""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.stats import qmc

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_outputs" / "oneD_v2_terminal_predictive_program"
sys.path.insert(0, str(ROOT))

from arrhenius_fracture.emergent_gnd_campaign_v913 import candidate_from_registry_row
from arrhenius_fracture.emergent_gnd_contract_v913 import ACTIVE_CANDIDATE_PARAMETER_FIELDS
from scripts.run_oneD_v2_predictive_campaign import inputs
from scripts.search_oneD_v2_provider_robust import source_rows
from scripts.search_oneD_v2_terminal_four_class import (
    CLASS_SEEDS,
    MAXIMUM_INTERVALS,
    OBJECTIVES,
    SCREEN_TEMPERATURES,
    VALIDATION_SEEDS,
    VALIDATION_TEMPERATURES,
    case_record,
    nondominated,
    objectives,
)

ANCHORS = {
    "weak-T": ("v913_zeroD_sobol_0116955", "v913_zeroD_sobol_0129902"),
    "ceramic-like": ("v913_zeroD_sobol_0134035", "v913_zeroD_sobol_0077080"),
}
INTERPOLATED_FIELDS = (
    "cleave_G00_eV", "cleave_gT_eV_per_K", "cleave_sigc0_GPa",
    "cleave_sT_GPa_per_K", "cleave_exp_a", "cleave_exp_n",
    "emit_G00_eV", "emit_gT_eV_per_K", "emit_sigc0_GPa",
    "emit_sT_GPa_per_K", "emit_exp_a", "emit_exp_n",
    "peierls_H0_eV", "peierls_activation_entropy_kB",
    "taylor_H0_eV", "taylor_activation_entropy_kB",
    "rho_source0_m2", "taylor_corr_rho_c_m2", "taylor_corr_scale", "c_blunt",
)
SIGNED = {
    "cleave_gT_eV_per_K", "cleave_sT_GPa_per_K",
    "emit_gT_eV_per_K", "emit_sT_GPa_per_K",
    "peierls_activation_entropy_kB", "taylor_activation_entropy_kB",
}


def focused_population(material: str) -> pd.DataFrame:
    source = source_rows().set_index("candidate_id")
    topology = source.loc[ANCHORS[material][0]]
    trend = source.loc[ANCHORS[material][1]]
    sampler = qmc.Sobol(len(INTERPOLATED_FIELDS), scramble=True,
                        seed=161803 if material == "weak-T" else 141421)
    rows = []
    for index, vector in enumerate(sampler.random_base2(6)):
        row = topology.copy()
        row["candidate_id"] = f"oneD_v2_focused_{material.replace('-', '_')}_{index:04d}"
        row["search_source"] = "FOCUSED_TOPOLOGY_TREND_HYPERRECTANGLE"
        row["search_anchor"] = "|".join(ANCHORS[material])
        for field, alpha in zip(INTERPOLATED_FIELDS, vector):
            a, b = float(topology[field]), float(trend[field])
            if field in SIGNED or a <= 0.0 or b <= 0.0:
                row[field] = (1.0 - float(alpha)) * a + float(alpha) * b
            else:
                row[field] = np.exp((1.0 - float(alpha)) * np.log(a) + float(alpha) * np.log(b))
        try:
            candidate = candidate_from_registry_row(row)
            for temperature in (300.0, 1200.0):
                assert candidate.cleavage.zero_stress_eV(temperature) > 1.0e-6
                assert candidate.emission.zero_stress_eV(temperature) > 1.0e-6
                assert candidate.cleavage.characteristic_stress_Pa(temperature) > 1.0e6
                assert candidate.emission.characteristic_stress_Pa(temperature) > 1.0e6
        except (AssertionError, ValueError):
            continue
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> int:
    physics, _, providers = inputs()
    existing_cases_path = OUT / "oneD_v2_search_population.parquet"
    existing_pareto_path = OUT / "oneD_v2_pareto_candidates.csv"
    existing_cases = pd.read_parquet(existing_cases_path)
    existing_pareto = pd.read_csv(existing_pareto_path)
    new_cases = []
    new_metrics = []
    checkpoint = OUT / "oneD_v2_terminal_refinement_checkpoint.parquet"
    for material in ANCHORS:
        candidates = focused_population(material)
        screen_records = []
        for index, (_, row) in enumerate(candidates.iterrows(), start=1):
            for temperature in SCREEN_TEMPERATURES:
                for provider, (mechanics, drive) in providers.items():
                    screen_records.append(case_record(
                        row, material, temperature, provider, mechanics, drive,
                        physics, CLASS_SEEDS[material], "FOCUSED_SCREEN",
                    ))
            if index % 4 == 0 or index == len(candidates):
                pd.DataFrame(new_cases + screen_records).to_parquet(checkpoint, index=False)
                print(f"REFINE_PROGRESS class={material} {index}/{len(candidates)}", flush=True)
        new_cases.extend(screen_records)
        screen = pd.DataFrame(screen_records)
        by_id = candidates.set_index("candidate_id")
        metric_rows = []
        for candidate_id, local in screen.groupby("candidate_id", sort=False):
            metric_rows.append({"candidate_id": candidate_id, **objectives(local, material)})
        metric = pd.DataFrame(metric_rows)
        metric["pareto"] = nondominated(metric)
        primary = metric.sort_values([
            "completion_objective", "physical_topology_objective",
            "precursor_objective", "class_trend_objective",
            "onset_scale_objective", "provider_robustness_objective",
        ], kind="stable")
        trend_first = metric.sort_values([
            "completion_objective", "class_trend_objective",
            "physical_topology_objective", "precursor_objective",
            "onset_scale_objective", "provider_robustness_objective",
        ], kind="stable")
        shortlist = list(primary.head(4).candidate_id)
        for candidate_id in trend_first.head(3).candidate_id:
            if candidate_id not in shortlist:
                shortlist.append(candidate_id)
        for candidate_id in metric[metric.pareto].candidate_id:
            if candidate_id not in shortlist and len(shortlist) < 8:
                shortlist.append(candidate_id)
        validation = []
        seeds = VALIDATION_SEEDS if material == "weak-T" else tuple(x + 1000000 for x in VALIDATION_SEEDS)
        for candidate_id in shortlist:
            row = by_id.loc[candidate_id].copy()
            row["candidate_id"] = candidate_id
            for seed in seeds:
                for temperature in VALIDATION_TEMPERATURES:
                    for provider, (mechanics, drive) in providers.items():
                        validation.append(case_record(
                            row, material, temperature, provider, mechanics, drive,
                            physics, seed, "FOCUSED_REPEATED_VALIDATION",
                        ))
            print(f"REFINE_VALIDATION class={material} candidate={candidate_id}", flush=True)
        new_cases.extend(validation)
        validation = pd.DataFrame(validation)
        for candidate_id, local in validation.groupby("candidate_id", sort=False):
            record = {
                "search_class": material,
                "candidate_id": candidate_id,
                "search_source": by_id.loc[candidate_id].search_source,
                "search_anchor": by_id.loc[candidate_id].search_anchor,
                **objectives(local, material),
            }
            for field in ACTIVE_CANDIDATE_PARAMETER_FIELDS:
                record[field] = float(by_id.loc[candidate_id][field])
            new_metrics.append(record)
    combined_cases = pd.concat([existing_cases, pd.DataFrame(new_cases)], ignore_index=True)
    combined_cases.to_parquet(existing_cases_path, index=False)
    combined = pd.concat([existing_pareto, pd.DataFrame(new_metrics)], ignore_index=True)
    combined["pareto_nondominated_validation"] = False
    for material, indices in combined.groupby("search_class").groups.items():
        combined.loc[list(indices), "pareto_nondominated_validation"] = nondominated(combined.loc[list(indices)])
    combined["selection_status"] = np.where(
        combined.full_class_contract_pass,
        "ELIGIBLE_SHARED_CLASS_ROW",
        "PARETO_DIAGNOSTIC_NOT_FULL_CONTRACT",
    )
    combined.sort_values(["search_class", *OBJECTIVES], kind="stable").to_csv(existing_pareto_path, index=False)
    checkpoint.unlink(missing_ok=True)
    manifest_path = OUT / "oneD_v2_parameter_search_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest.update({
        "focused_refinement": True,
        "focused_anchors": ANCHORS,
        "focused_candidate_count": int(sum(len(focused_population(x)) for x in ANCHORS)),
        "focused_validation_candidate_count": int(len(new_metrics)),
        "total_case_count": int(len(combined_cases)),
        "eligible_rows": combined.loc[combined.full_class_contract_pass, ["search_class", "candidate_id"]].to_dict("records"),
    })
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(f"REFINEMENT_COMPLETE cases={len(combined_cases)} eligible={int(combined.full_class_contract_pass.sum())}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
