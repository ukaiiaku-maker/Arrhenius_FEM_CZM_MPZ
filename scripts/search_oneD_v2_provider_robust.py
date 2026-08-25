#!/usr/bin/env python3
"""Provider-robust V2 DBTT search with separate Pareto objectives.

The search is deliberately confined to the reduced model.  It reads the
qualified, candidate-independent mechanics and tensor-drive maps and never
launches a two-dimensional trajectory.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.stats import qmc

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_outputs" / "oneD_v2_predictive_model"
PF_MATERIALS = Path(
    "/Volumes/Data/Data/Nanopillar_calculation/"
    "PF-fracture-fatigue_v10_2_21_persistent_sites_top1/"
    "arrhenius_fracture/data/materials"
)
sys.path.insert(0, str(ROOT))

from arrhenius_fracture.emergent_gnd_contract_v913 import (  # noqa: E402
    ACTIVE_CANDIDATE_PARAMETER_FIELDS,
)
from scripts.run_oneD_v2_predictive_campaign import (  # noqa: E402
    inputs,
    run_case,
    summary,
)

SCREEN_TEMPERATURES = (300.0, 600.0, 1000.0, 1200.0)
VALIDATION_TEMPERATURES = (300.0, 600.0, 900.0, 1000.0, 1100.0, 1200.0)
BASE_SEED = 1008666
VALIDATION_SEEDS = (1008666, 1008667, 1008668)
DBTT_CURRENT = "v913_zeroD_sobol_0202500"
DBTT_ALTERNATE = "v913_zeroD_sobol_0116955"

LOG_PERTURB = {
    "cleave_G00_eV": 0.16,
    "cleave_sigc0_GPa": 0.18,
    "cleave_exp_a": 0.20,
    "cleave_exp_n": 0.18,
    "emit_G00_eV": 0.20,
    "emit_sigc0_GPa": 0.28,
    "emit_exp_a": 0.25,
    "emit_exp_n": 0.25,
    "rho_source0_m2": 0.60,
    "taylor_corr_rho_c_m2": 0.60,
    "taylor_corr_scale": 0.30,
    "c_blunt": 0.45,
}
SIGNED_FIELDS = {
    "cleave_gT_eV_per_K",
    "cleave_sT_GPa_per_K",
    "emit_gT_eV_per_K",
    "emit_sT_GPa_per_K",
    "peierls_activation_entropy_kB",
    "taylor_activation_entropy_kB",
}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_rows() -> pd.DataFrame:
    names = (
        "v10_2_23_v913_top10_persistent_site_registry.csv",
        "v10_2_24_v913_top10_upper_shelf_registry.csv",
        "v10_2_25_v913_paper_campaign_registry.csv",
        "v10_2_27_v913_four_class_paper_registry.csv",
    )
    frames = []
    for name in names:
        frame = pd.read_csv(PF_MATERIALS / name)
        frame["search_source"] = name
        frames.append(frame)
    return pd.concat(frames, ignore_index=True).drop_duplicates("candidate_id")


def interpolate(a: pd.Series, b: pd.Series, alpha: float, candidate_id: str) -> pd.Series:
    row = a.copy()
    row["candidate_id"] = candidate_id
    row["search_source"] = "LOG_PARAMETER_MORPH"
    row["search_anchor"] = f"{a.candidate_id}|{b.candidate_id}"
    for field in ACTIVE_CANDIDATE_PARAMETER_FIELDS:
        av, bv = float(a[field]), float(b[field])
        if field in SIGNED_FIELDS or av <= 0.0 or bv <= 0.0:
            row[field] = (1.0 - alpha) * av + alpha * bv
        else:
            row[field] = np.exp((1.0 - alpha) * np.log(av) + alpha * np.log(bv))
    return row


def population() -> pd.DataFrame:
    sources = source_rows()
    archival = sources[
        sources.material_class.astype(str).str.lower().isin(("dbtt", "peak"))
    ].copy()
    archival["search_source"] = "ARCHIVED_V913_POOL"
    archival["search_anchor"] = archival.candidate_id
    a = sources[sources.candidate_id == DBTT_CURRENT].iloc[0]
    b = sources[sources.candidate_id == DBTT_ALTERNATE].iloc[0]
    morphs = [
        interpolate(a, b, float(alpha), f"oneD_v2_dbtt_morph_{index:03d}")
        for index, alpha in enumerate(np.linspace(0.025, 0.975, 39), start=1)
    ]
    anchors = (a, b)
    sampler = qmc.Sobol(len(LOG_PERTURB), scramble=True, seed=271828)
    unit = sampler.random_base2(6)  # 64 local candidates per anchor.
    locals_: list[pd.Series] = []
    fields = tuple(LOG_PERTURB)
    for anchor_index, anchor in enumerate(anchors):
        for sample_index, vector in enumerate(unit):
            row = anchor.copy()
            row["candidate_id"] = (
                f"oneD_v2_dbtt_local_a{anchor_index}_{sample_index:03d}"
            )
            row["search_source"] = "SOBOL_LOCAL_REDUCED_SCREEN"
            row["search_anchor"] = anchor.candidate_id
            for field, u in zip(fields, vector):
                row[field] = float(anchor[field]) * 10.0 ** (
                    (2.0 * float(u) - 1.0) * LOG_PERTURB[field]
                )
            locals_.append(row)
    combined = pd.concat(
        [archival, pd.DataFrame(morphs), pd.DataFrame(locals_)],
        ignore_index=True,
    ).drop_duplicates("candidate_id")
    return combined


def case_record(row, temperature, provider, mechanics, drive, physics, seed, stage):
    result = run_case(
        row, "DBTT", temperature, provider, mechanics, drive, physics, 100.0,
        seed=seed,
    )
    record = summary(result)
    record.update(
        {
            "search_stage": stage,
            "hazard_seed": int(seed),
            "search_source": row.get("search_source", ""),
            "search_anchor": row.get("search_anchor", ""),
        }
    )
    return record


def candidate_objectives(cases: pd.DataFrame) -> dict[str, float | int | bool]:
    def value(provider, temperature, field):
        local = cases[(cases.provider == provider) & (cases.temperature_K == temperature)]
        return float(local[field].mean())

    status = int((cases.status != "TARGET_RIGHT_CENSORED").sum())
    shelf = {}
    gain = {}
    drop = {}
    topology_error = 0.0
    precursor_error = 0.0
    largest_error = 0.0
    event_size_error = 0.0
    normalized_curves = {}
    for provider, target_avalanches, target_largest in (
        ("PF", 2.0, 0.55),
        ("FEMCZM", 3.0, 0.88),
    ):
        low = np.array(
            [value(provider, t, "onset_envelope_max_MPa_sqrt_m") for t in (300.0, 600.0)]
        )
        middle = np.array(
            [
                value(provider, t, "onset_envelope_max_MPa_sqrt_m")
                for t in sorted(set(cases.temperature_K) & {900.0, 1000.0, 1100.0})
            ]
        )
        high = value(provider, 1200.0, "onset_envelope_max_MPa_sqrt_m")
        shelf[provider] = float(np.ptp(low))
        gain[provider] = float(np.max(middle) - np.mean(low))
        drop[provider] = float(np.max(middle) - high)
        at_1000 = cases[(cases.provider == provider) & (cases.temperature_K == 1000.0)]
        avalanches = float(at_1000.physical_avalanche_count.mean())
        precursors = float(at_1000.precursor_reinitiation_count.mean())
        topology_error += abs(avalanches - target_avalanches)
        precursor_error += abs(precursors - (target_avalanches - 1.0))
        largest_error += abs(float(at_1000.largest_avalanche_fraction.mean()) - target_largest)
        event_size_error += abs(float(at_1000.mean_event_size_um.mean()) - 5.0) / 5.0
        ordered = cases[cases.provider == provider].groupby("temperature_K")[
            "onset_envelope_max_MPa_sqrt_m"
        ].mean().sort_index().to_numpy(float)
        normalized_curves[provider] = ordered - float(np.mean(low))
    provider_robustness = float(
        np.mean(np.abs(normalized_curves["PF"] - normalized_curves["FEMCZM"])) / 20.0
    )
    onset_objective = (
        sum(max(v - 10.0, 0.0) / 10.0 for v in shelf.values())
        + sum(max(8.0 - v, 0.0) / 8.0 for v in gain.values())
        + sum(max(-v, 0.0) / 5.0 for v in drop.values())
    )
    radius = float(cases.max_tip_radius_um.max())
    backstress = float(cases.max_backstress_GPa.max())
    width = float(cases.minimum_front_width_um.min())
    state_realism = (
        max(radius - 100.0, 0.0) / 100.0
        + max(backstress - 50.0, 0.0) / 50.0
        + max(0.01 - width, 0.0) / 0.01
    )
    return {
        "status_objective": status,
        "class_topology_objective": topology_error,
        "onset_envelope_objective": onset_objective,
        "precursor_count_objective": precursor_error,
        "avalanche_topology_objective": topology_error + largest_error,
        "event_size_objective": event_size_error,
        "provider_robustness_objective": provider_robustness,
        "state_realism_objective": state_realism,
        "PF_low_shelf_span": shelf["PF"],
        "FEMCZM_low_shelf_span": shelf["FEMCZM"],
        "PF_transition_gain": gain["PF"],
        "FEMCZM_transition_gain": gain["FEMCZM"],
        "PF_high_temperature_drop": drop["PF"],
        "FEMCZM_high_temperature_drop": drop["FEMCZM"],
        "all_cases_complete": status == 0,
        "full_dbtt_contract_pass": bool(
            status == 0
            and topology_error <= 1.0
            and precursor_error <= 1.0
            and all(v <= 10.0 for v in shelf.values())
            and all(v >= 8.0 for v in gain.values())
            and all(v >= 0.0 for v in drop.values())
            and state_realism == 0.0
        ),
    }


OBJECTIVES = (
    "status_objective",
    "class_topology_objective",
    "onset_envelope_objective",
    "precursor_count_objective",
    "avalanche_topology_objective",
    "event_size_objective",
    "provider_robustness_objective",
    "state_realism_objective",
)


def nondominated(frame: pd.DataFrame) -> np.ndarray:
    values = frame.loc[:, OBJECTIVES].to_numpy(float)
    keep = np.ones(len(frame), dtype=bool)
    for index, row in enumerate(values):
        dominated = np.all(values <= row, axis=1) & np.any(values < row, axis=1)
        dominated[index] = False
        keep[index] = not bool(np.any(dominated))
    return keep


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    physics, _, providers = inputs()
    candidates = population()
    checkpoint = OUT / "oneD_v2_search_screen_checkpoint.parquet"
    if checkpoint.exists():
        cached = pd.read_parquet(checkpoint)
        records = cached.to_dict("records")
        completed_ids = set(cached.candidate_id.astype(str))
        print(f"SCREEN_RESUME {len(completed_ids)}/{len(candidates)}", flush=True)
    else:
        records = []
        completed_ids = set()
    total = len(candidates)
    for index, (_, row) in enumerate(candidates.iterrows(), start=1):
        if str(row.candidate_id) in completed_ids:
            continue
        for temperature in SCREEN_TEMPERATURES:
            for provider, (mechanics, drive) in providers.items():
                records.append(
                    case_record(
                        row, temperature, provider, mechanics, drive, physics,
                        BASE_SEED, "SCREEN",
                    )
                )
        if index % 10 == 0 or index == total:
            pd.DataFrame(records).to_parquet(checkpoint, index=False)
            print(f"SCREEN_PROGRESS {index}/{total}", flush=True)
    screen = pd.DataFrame(records)
    metrics = []
    by_id = candidates.set_index("candidate_id")
    for candidate_id, local in screen.groupby("candidate_id", sort=False):
        row = {
            "candidate_id": candidate_id,
            "search_source": by_id.loc[candidate_id].get("search_source", ""),
            "search_anchor": by_id.loc[candidate_id].get("search_anchor", ""),
            **candidate_objectives(local),
        }
        metrics.append(row)
    metrics = pd.DataFrame(metrics)
    metrics["pareto_nondominated_screen"] = nondominated(metrics)
    ranking = metrics.sort_values(
        [
            "status_objective",
            "class_topology_objective",
            "onset_envelope_objective",
            "provider_robustness_objective",
            "state_realism_objective",
        ]
    )
    shortlist = list(ranking.head(10).candidate_id)
    for candidate_id in metrics[metrics.pareto_nondominated_screen].candidate_id:
        if candidate_id not in shortlist and len(shortlist) < 16:
            shortlist.append(candidate_id)
    validation_records = []
    for index, candidate_id in enumerate(shortlist, start=1):
        row = by_id.loc[candidate_id].copy()
        row["candidate_id"] = candidate_id
        for seed in VALIDATION_SEEDS:
            for temperature in VALIDATION_TEMPERATURES:
                for provider, (mechanics, drive) in providers.items():
                    validation_records.append(
                        case_record(
                            row, temperature, provider, mechanics, drive, physics,
                            seed, "REPEATED_VALIDATION",
                        )
                    )
        print(f"VALIDATION_PROGRESS {index}/{len(shortlist)}", flush=True)
    validation = pd.DataFrame(validation_records)
    validated_metrics = []
    for candidate_id, local in validation.groupby("candidate_id", sort=False):
        row = {
            "candidate_id": candidate_id,
            "search_source": by_id.loc[candidate_id].get("search_source", ""),
            "search_anchor": by_id.loc[candidate_id].get("search_anchor", ""),
            **candidate_objectives(local),
        }
        for field in ACTIVE_CANDIDATE_PARAMETER_FIELDS:
            row[field] = float(by_id.loc[candidate_id][field])
        validated_metrics.append(row)
    pareto = pd.DataFrame(validated_metrics)
    pareto["pareto_nondominated_validation"] = nondominated(pareto)
    pareto["selection_status"] = np.where(
        pareto.full_dbtt_contract_pass,
        "ELIGIBLE_SHARED_DBTT_ROW",
        "REJECTED_FULL_DBTT_CONTRACT",
    )
    combined_cases = pd.concat([screen, validation], ignore_index=True)
    combined_cases.to_parquet(OUT / "oneD_v2_search_population.parquet", index=False)
    checkpoint.unlink(missing_ok=True)
    pareto.sort_values(
        ["full_dbtt_contract_pass", "class_topology_objective", "onset_envelope_objective"],
        ascending=[False, True, True],
    ).to_csv(OUT / "oneD_v2_pareto_candidates.csv", index=False)
    manifest = {
        "schema": "oneD_v2_provider_robust_search_v1",
        "screen_candidate_count": int(len(candidates)),
        "screen_case_count": int(len(screen)),
        "validation_candidate_count": int(len(shortlist)),
        "validation_case_count": int(len(validation)),
        "screen_temperatures_K": list(SCREEN_TEMPERATURES),
        "validation_temperatures_K": list(VALIDATION_TEMPERATURES),
        "validation_seeds": list(VALIDATION_SEEDS),
        "objectives": list(OBJECTIVES),
        "eligible_shared_dbtt_rows": pareto.loc[
            pareto.full_dbtt_contract_pass, "candidate_id"
        ].tolist(),
        "new_2D_PF_runs": 0,
        "new_2D_FEMCZM_runs": 0,
        "source_registry_hashes": {
            name: sha(PF_MATERIALS / name)
            for name in (
                "v10_2_23_v913_top10_persistent_site_registry.csv",
                "v10_2_24_v913_top10_upper_shelf_registry.csv",
                "v10_2_25_v913_paper_campaign_registry.csv",
                "v10_2_27_v913_four_class_paper_registry.csv",
            )
        },
    }
    (OUT / "oneD_v2_parameter_search_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    print(
        "SEARCH_COMPLETE "
        f"screen={len(candidates)} validation={len(shortlist)} "
        f"eligible={int(pareto.full_dbtt_contract_pass.sum())}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
