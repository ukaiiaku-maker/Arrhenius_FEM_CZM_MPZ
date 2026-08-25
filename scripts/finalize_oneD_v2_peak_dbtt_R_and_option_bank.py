#!/usr/bin/env python3
"""Finalize the focused R search and the versioned four-class option bank."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_outputs" / "oneD_v2_peak_dbtt_rcurve_search"
OLD = ROOT / "analysis_outputs" / "oneD_v2_terminal_predictive_program"
sys.path.insert(0, str(ROOT))

from arrhenius_fracture.emergent_gnd_contract_v913 import ACTIVE_CANDIDATE_PARAMETER_FIELDS
from scripts.run_oneD_v2_peak_dbtt_R_screen import OBJECTIVES as R_OBJECTIVES
from scripts.search_oneD_v2_terminal_four_class import (
    OBJECTIVES as OLD_OBJECTIVES, nondominated, objectives, population as old_population,
)


BANK_VERSION = "oneD_v2_fracture_option_bank_v1"
CLASS_ORDER = ("Peak", "DBTT", "weak-T", "ceramic-like")
CONTROLS = {
    "Peak": "v913_zeroD_sobol_0242980",
    "DBTT": "v913_zeroD_sobol_0202500",
    "weak-T": "oneD_v2_focused_weak_T_0016",
    "ceramic-like": "oneD_v2_focused_ceramic_like_0018",
}
FOCUSED = {
    "Peak": [
        "oneD_v2_peak_R_41f8789bcbc1f097",
        "oneD_v2_peak_R_a2e923a7587543db",
        "oneD_v2_peak_R_7834115ae79cd559",
    ],
    "DBTT": [
        "oneD_v2_dbtt_R_9f5160f509e713e2",
        "oneD_v2_dbtt_R_c4859a34963f15af",
        "oneD_v2_dbtt_R_7b3b5b45354e64ef",
        "oneD_v2_dbtt_R_3e168d9381eafa7b",
        "oneD_v2_dbtt_R_41a3f2096b3c4c54",
        "oneD_v2_dbtt_R_c2dabb42101cf812",
    ],
}
TARGET_BANK_COUNTS = {"Peak": 16, "DBTT": 16, "weak-T": 14, "ceramic-like": 12}
SHORTLIST_COUNTS = {"Peak": 6, "DBTT": 6, "weak-T": 4, "ceramic-like": 4}
FORBIDDEN_MATERIAL_FIELDS = {
    "hazard_progress_scale", "reload_gap_threshold_m", "nominal_advance_m",
    "translation_length_scale", "event_length_mode", "threshold_mode",
    "map_bounds", "exact_oracle_cache", "provider_correction_factor",
}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def git_head(path: Path = ROOT) -> str:
    return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=path, text=True).strip()


def canonical_json(row: pd.Series) -> str:
    return json.dumps(
        {field: float(row[field]) for field in ACTIVE_CANDIDATE_PARAMETER_FIELDS},
        sort_keys=True, separators=(",", ":"), allow_nan=False,
    )


def material_hash(row: pd.Series) -> str:
    return hashlib.sha256(canonical_json(row).encode()).hexdigest()


def units() -> dict[str, str]:
    result = {}
    for field in ACTIVE_CANDIDATE_PARAMETER_FIELDS:
        if field.endswith("_eV"):
            result[field] = "eV"
        elif field.endswith("_eV_per_K"):
            result[field] = "eV/K"
        elif field.endswith("_GPa"):
            result[field] = "GPa"
        elif field.endswith("_GPa_per_K"):
            result[field] = "GPa/K"
        elif field.endswith("_m2"):
            result[field] = "m^-2"
        elif field.endswith("_s"):
            result[field] = "s^-1"
        else:
            result[field] = "dimensionless"
    return result


def candidate_frames() -> dict[str, pd.DataFrame]:
    frames = {
        "Peak": pd.read_csv(OUT / "oneD_v2_peak_R_candidates.csv"),
        "DBTT": pd.read_csv(OUT / "oneD_v2_dbtt_R_candidates.csv"),
    }
    prior_registry = pd.read_csv(OLD / "oneD_v2_new_four_class_registry.csv")
    for material in ("weak-T", "ceramic-like"):
        frame = old_population(material).copy()
        control = prior_registry[
            (prior_registry.material_class == material)
            & (prior_registry.candidate_id == CONTROLS[material])
        ].copy()
        if len(control) != 1:
            raise RuntimeError(f"missing unique prior selected row for {material}")
        control["search_anchor"] = "PRIOR_QUALIFIED_SELECTED_ROW"
        control["search_source"] = "PRIOR_QUALIFIED_SELECTED_ROW"
        frame = pd.concat([frame, control], ignore_index=True, sort=False)
        frame["target_response_class"] = material
        frame["search_campaign_id"] = "oneD_v2_terminal_shared_four_class_search_v1"
        frame["parent_or_anchor_id"] = frame.search_anchor
        frame["search_generation_method"] = frame.search_source
        frame["active_parameter_schema"] = "V913_ACTIVE_CANDIDATE_29"
        frame["active_parameter_count"] = len(ACTIVE_CANDIDATE_PARAMETER_FIELDS)
        frame["canonical_parameter_json"] = frame.apply(canonical_json, axis=1)
        frame["parameter_sha256"] = frame.apply(material_hash, axis=1)
        frames[material] = frame
    return frames


def response_frames() -> dict[str, pd.DataFrame]:
    result = {}
    for material, slug in (("Peak", "peak"), ("DBTT", "dbtt")):
        screen = pd.read_parquet(OUT / f"oneD_v2_{slug}_R_search_population.parquet").copy()
        screen["response_source"] = "FOCUSED_FAST_SCREEN"
        validation = pd.read_csv(OUT / "oneD_v2_peak_dbtt_R_multiseed_validation.csv")
        validation = validation[validation.target_response_class == material].copy()
        validation["response_source"] = "FOCUSED_MULTISEED_VALIDATION"
        result[material] = pd.concat([screen, validation], ignore_index=True, sort=False)
    old = pd.read_parquet(OLD / "oneD_v2_search_population.parquet")
    for material in ("weak-T", "ceramic-like"):
        local = old[old.material_class == material].copy()
        local["target_response_class"] = material
        local["N_reinit"] = local.precursor_reinitiation_count
        local["initial_onset_native_KJ_MPa_sqrt_m"] = local.first_event_native_KJ_MPa_sqrt_m
        local["maximum_onset_native_KJ_MPa_sqrt_m"] = local.onset_envelope_max_MPa_sqrt_m
        local["relative_deltaK_reinit"] = (
            local.maximum_onset_native_KJ_MPa_sqrt_m - local.initial_onset_native_KJ_MPa_sqrt_m
        ) / local.initial_onset_native_KJ_MPa_sqrt_m.clip(lower=1.0e-12)
        local["response_source"] = "PRIOR_COMPLETED_SHARED_FOUR_CLASS_SEARCH"
        result[material] = local
    return result


def _objective_library(material: str, candidates: pd.DataFrame,
                       responses: pd.DataFrame) -> pd.DataFrame:
    if material in ("Peak", "DBTT"):
        slug = material.lower()
        return pd.read_csv(OUT / f"oneD_v2_{slug}_R_pareto.csv")
    by_id = candidates.set_index("candidate_id")
    rows = []
    for candidate_id, local in responses[responses.search_stage == "SCREEN"].groupby("candidate_id"):
        row = {
            "candidate_id": candidate_id,
            "target_response_class": material,
            **objectives(local, material),
        }
        source = by_id.loc[candidate_id]
        row.update({field: float(source[field]) for field in ACTIVE_CANDIDATE_PARAMETER_FIELDS})
        row["parameter_sha256"] = source.parameter_sha256
        row["canonical_parameter_json"] = source.canonical_parameter_json
        rows.append(row)
    frame = pd.DataFrame(rows)
    frame["pareto_nondominated"] = nondominated(frame.rename(
        columns={name: name for name in OLD_OBJECTIVES}
    ))
    values = frame[list(OLD_OBJECTIVES)].to_numpy(float)
    lo, hi = np.nanmin(values, axis=0), np.nanmax(values, axis=0)
    normalized = (values - lo) / np.maximum(hi - lo, 1.0e-12)
    front = normalized[frame.pareto_nondominated]
    frame["distance_to_pareto_front"] = np.min(
        np.linalg.norm(normalized[:, None] - front[None], axis=2), axis=1
    )
    frame["pareto_near_front"] = frame.distance_to_pareto_front <= 0.20
    frame["option_status"] = np.where(
        frame.pareto_nondominated, "PARETO_NONDOMINATED",
        np.where(frame.pareto_near_front, "PARETO_NEAR_FRONT", "SCREENED"),
    )
    return frame[frame.pareto_nondominated | frame.pareto_near_front].copy()


def _candidate_features(material: str, library: pd.DataFrame,
                        responses: pd.DataFrame) -> pd.DataFrame:
    objective_names = list(R_OBJECTIVES if material in ("Peak", "DBTT") else OLD_OBJECTIVES)
    available = [name for name in objective_names if name in library]
    features = library[["candidate_id", *available]].copy()
    aggregate = responses.groupby("candidate_id").agg(
        onset_mean=("initial_onset_native_KJ_MPa_sqrt_m", "mean"),
        onset_min=("initial_onset_native_KJ_MPa_sqrt_m", "min"),
        onset_max=("initial_onset_native_KJ_MPa_sqrt_m", "max"),
        reinit_mean=("relative_deltaK_reinit", "mean"),
        reinit_max=("relative_deltaK_reinit", "max"),
        precursor_mean=("N_reinit", "mean"),
        avalanche_fraction_mean=("largest_avalanche_fraction", "mean"),
        radius_max=("max_tip_radius_um", "max"),
        backstress_max=("max_backstress_GPa", "max"),
    ).reset_index()
    return features.merge(aggregate, on="candidate_id", how="left")


def maximin_select(material: str, library: pd.DataFrame, candidates: pd.DataFrame,
                   responses: pd.DataFrame, count: int, required: list[str]) -> list[str]:
    pool = library.drop_duplicates("candidate_id").copy()
    missing_required = [item for item in required if item not in set(pool.candidate_id)]
    if missing_required:
        pool = pd.concat([
            pool,
            pd.DataFrame({"candidate_id": missing_required}),
        ], ignore_index=True, sort=False)
    params = candidates.drop_duplicates("candidate_id").set_index("candidate_id")
    pool = pool[pool.candidate_id.isin(params.index)].copy()
    feature = _candidate_features(material, pool, responses).set_index("candidate_id")
    matrix_parts = [params.loc[pool.candidate_id, list(ACTIVE_CANDIDATE_PARAMETER_FIELDS)].to_numpy(float)]
    matrix_parts.append(feature.loc[pool.candidate_id].select_dtypes(include=[np.number]).to_numpy(float))
    values = np.concatenate(matrix_parts, axis=1)
    values = np.nan_to_num(values, nan=np.nanmedian(values, axis=0), posinf=0.0, neginf=0.0)
    lo, hi = np.min(values, axis=0), np.max(values, axis=0)
    values = (values - lo) / np.maximum(hi - lo, 1.0e-30)
    ids = list(pool.candidate_id.astype(str))
    selected = []
    for candidate_id in required:
        if candidate_id in ids and candidate_id not in selected:
            selected.append(candidate_id)
    while len(selected) < min(count, len(ids)):
        if not selected:
            selected.append(ids[0]); continue
        chosen = [ids.index(item) for item in selected]
        distances = np.min(
            np.linalg.norm(values[:, None, :] - values[chosen][None, :, :], axis=2), axis=1
        )
        for index in chosen:
            distances[index] = -1.0
        selected.append(ids[int(np.argmax(distances))])
    return selected


def _response_json(records: pd.DataFrame) -> str:
    keep = [
        "response_source", "provider", "temperature_K", "hazard_seed", "target_um",
        "status", "terminal_reason", "initial_onset_native_KJ_MPa_sqrt_m",
        "maximum_onset_native_KJ_MPa_sqrt_m", "deltaK_reinit_MPa_sqrt_m",
        "relative_deltaK_reinit", "N_reinit", "physical_avalanche_count",
        "largest_avalanche_fraction", "mean_event_size_um", "median_event_size_um",
        "reload_increment_min_um", "reload_increment_max_um", "max_tip_radius_um",
        "minimum_front_width_um", "max_backstress_GPa", "max_source_multiplicity",
        "maximum_onset_mobile_density_m2", "maximum_onset_retained_density_m2",
        "signed_shielding_MPa_sqrt_m", "map_oracle_fallback_count",
    ]
    present = [name for name in keep if name in records]
    local = records[present].replace({np.nan: None})
    return json.dumps(local.to_dict("records"), sort_keys=True, separators=(",", ":"), default=str)


def _roles(material: str, selected: list[str], candidates: pd.DataFrame,
           responses: pd.DataFrame) -> dict[str, list[str]]:
    local = responses[responses.candidate_id.isin(selected)].groupby("candidate_id").agg(
        onset=("initial_onset_native_KJ_MPa_sqrt_m", "mean"),
        reinit=("relative_deltaK_reinit", "max"),
        precursors=("N_reinit", "mean"),
        fraction=("largest_avalanche_fraction", "mean"),
        radius=("max_tip_radius_um", "max"),
        backstress=("max_backstress_GPa", "max"),
    )
    params = candidates.drop_duplicates("candidate_id").set_index("candidate_id")
    roles = {candidate: ["CURATED_OPTION", "DIVERSE_PARETO"] for candidate in selected}
    roles[CONTROLS[material]].extend(["CONTROL", "BALANCED"])
    for field, low_role, high_role in (
        ("onset", "LOW_ONSET", "HIGH_ONSET"),
        ("reinit", "LOW_PRECURSOR_COUNT", "R_ENRICHED"),
        ("precursors", "LOW_PRECURSOR_COUNT", "HIGH_PRECURSOR_COUNT"),
        ("fraction", "LOW_FINAL_AVALANCHE_FRACTION", "HIGH_FINAL_AVALANCHE_FRACTION"),
        ("radius", "LOW_BACKSTRESS", "HIGH_BACKSTRESS"),
        ("backstress", "LOW_BACKSTRESS", "HIGH_BACKSTRESS"),
    ):
        series = local[field].dropna()
        if len(series):
            roles[str(series.idxmin())].append(low_role)
            roles[str(series.idxmax())].append(high_role)
    blunt = params.loc[selected, "c_blunt"].astype(float)
    roles[str(blunt.idxmin())].append("LOW_BLUNTING")
    roles[str(blunt.idxmax())].append("HIGH_BLUNTING")
    if material == "DBTT":
        roles[FOCUSED[material][0]].extend(["LOW_T_SHELF", "SHIFTED_TRANSITION", "PROVIDER_SENSITIVE"])
    if material == "Peak":
        roles[FOCUSED[material][0]].extend(["HIGH_T_SOFTENING", "PROVIDER_SENSITIVE"])
    return {key: sorted(set(value)) for key, value in roles.items()}


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    producer_commit = git_head()
    source_commit = subprocess.check_output(
        ["git", "rev-parse", "5045fee"], cwd=ROOT, text=True
    ).strip()
    domain_path = OLD / "oneD_v2_domain_of_usefulness.json"
    domain_hash = sha(domain_path)
    candidates_by_class = candidate_frames()
    responses_by_class = response_frames()
    libraries = {
        material: _objective_library(material, candidates_by_class[material], responses_by_class[material])
        for material in CLASS_ORDER
    }

    all_material_entries = []
    all_response = []
    for material in CLASS_ORDER:
        candidates = candidates_by_class[material]
        responses = responses_by_class[material].copy()
        hashes = candidates.drop_duplicates("candidate_id").set_index("candidate_id").parameter_sha256
        responses["parameter_sha256"] = responses.candidate_id.map(hashes)
        all_response.append(responses)
        for row in candidates.itertuples(index=False):
            series = pd.Series(row._asdict())
            all_material_entries.append({
                "candidate_id": str(series.candidate_id),
                "target_response_class": material,
                "search_campaign_id": str(series.get("search_campaign_id", "")),
                "parent_or_anchor_id": str(series.get("parent_or_anchor_id", "")),
                "search_generation_method": str(series.get("search_generation_method", "")),
                "parameter_sha256": str(series.parameter_sha256),
                "canonical_parameter_json": str(series.canonical_parameter_json),
                **{field: float(series[field]) for field in ACTIVE_CANDIDATE_PARAMETER_FIELDS},
            })
    response_all = pd.concat(all_response, ignore_index=True, sort=False)
    raw_materials = pd.DataFrame(all_material_entries)
    population_rows = []
    for digest, group in raw_materials.groupby("parameter_sha256", sort=False):
        aliases = sorted(set(group.candidate_id.astype(str)))
        candidate_id = aliases[0]
        if any(alias.startswith("v913_") for alias in aliases):
            candidate_id = sorted(alias for alias in aliases if alias.startswith("v913_"))[0]
        row = group.iloc[0]
        records = response_all[response_all.parameter_sha256 == digest]
        population_rows.append({
            "candidate_id": candidate_id,
            "candidate_version": BANK_VERSION,
            "target_response_class": "|".join(sorted(set(group.target_response_class), key=CLASS_ORDER.index)),
            "search_campaign_id": "|".join(sorted(set(group.search_campaign_id))),
            "parent_or_anchor_id": "|".join(sorted(set(group.parent_or_anchor_id))),
            "search_generation_method": "|".join(sorted(set(group.search_generation_method))),
            "active_parameter_schema": "V913_ACTIVE_CANDIDATE_29",
            "active_parameter_count": len(ACTIVE_CANDIDATE_PARAMETER_FIELDS),
            "canonical_parameter_json": row.canonical_parameter_json,
            "parameter_sha256": digest,
            "candidate_aliases_json": json.dumps(aliases),
            "source_commit": source_commit,
            "analysis_commit": producer_commit,
            "registry_commit": producer_commit,
            "qualified_domain_sha256": domain_hash,
            "response_record_count": len(records),
            "response_records_json": _response_json(records),
            "fatigue_evaluated": False,
            "fatigue_validation_status": "NOT_EVALUATED",
            **{field: float(row[field]) for field in ACTIVE_CANDIDATE_PARAMETER_FIELDS},
        })
    complete = pd.DataFrame(population_rows).sort_values("candidate_id")
    complete_path = OUT / "oneD_v2_complete_material_option_population.parquet"
    complete.to_parquet(complete_path, index=False)

    pareto_parts = []
    for material, library in libraries.items():
        local = library.copy()
        local["library_target_response_class"] = material
        pareto_parts.append(local)
    pareto = pd.concat(pareto_parts, ignore_index=True, sort=False)
    pareto["parameter_sha256"] = pareto.apply(
        lambda row: material_hash(row), axis=1
    )
    pareto["candidate_version"] = BANK_VERSION
    pareto["qualified_domain_sha256"] = domain_hash
    pareto_path = OUT / "oneD_v2_pareto_material_option_library.parquet"
    pareto.drop_duplicates(["parameter_sha256", "library_target_response_class"]).to_parquet(
        pareto_path, index=False
    )

    selected_by_class = {}
    for material in CLASS_ORDER:
        required = [CONTROLS[material], *FOCUSED.get(material, [])]
        selected_by_class[material] = maximin_select(
            material, libraries[material], candidates_by_class[material],
            responses_by_class[material], TARGET_BANK_COUNTS[material], required,
        )
    direct = pd.read_csv(OUT / "pf_2d_peak_dbtt_R_transfer_summary.csv")
    old_direct = pd.read_csv(OLD / "oneD_v2_pf_transfer_results.csv")
    bank_rows = []
    for material in CLASS_ORDER:
        candidates = candidates_by_class[material].drop_duplicates("candidate_id").set_index("candidate_id")
        roles = _roles(material, selected_by_class[material], candidates.reset_index(), responses_by_class[material])
        for candidate_id in selected_by_class[material]:
            row = candidates.loc[candidate_id]
            digest = str(row.parameter_sha256)
            records = responses_by_class[material][responses_by_class[material].candidate_id == candidate_id]
            status = {"CURATED_OPTION"}
            if candidate_id == CONTROLS[material]:
                status.update({"HISTORICAL_CONTROL", "CURRENT_V2_SELECTED"})
            library_row = libraries[material][libraries[material].candidate_id == candidate_id]
            if len(library_row):
                if bool(library_row.iloc[0].get("pareto_nondominated", library_row.iloc[0].get("pareto_nondominated_validation", False))):
                    status.add("PARETO_NONDOMINATED")
                else:
                    status.add("PARETO_NEAR_FRONT")
            direct_status = "NOT_RUN"
            artifact_json = "[]"
            if candidate_id in set(direct.candidate_id) and bool(
                direct[direct.candidate_id == candidate_id].candidate_role.eq("R_FINALIST").any()
            ):
                q = direct[direct.candidate_id == candidate_id]
                direct_status = "FAILED_NO_POSITIVE_RELOAD_SEPARATED_RESISTANCE"
                status.update({"R_ENRICHED_OPTION", "DIRECT_PF_REJECTED", "PROVIDER_SENSITIVE_OPTION"})
                artifact_json = json.dumps(q[["temperature_K", "hazard_seed", "source_steps_file", "source_steps_sha256"]].to_dict("records"), sort_keys=True)
            elif candidate_id in set(direct.candidate_id):
                q = direct[direct.candidate_id == candidate_id]
                direct_status = "AVAILABLE_AUTHORITATIVE_TEMPERATURE_MATCHED_UNPAIRED_SEEDS"
                status.add("DIRECT_PF_VALIDATED")
                artifact_json = json.dumps(q[["temperature_K", "hazard_seed", "source_steps_file", "source_steps_sha256"]].to_dict("records"), sort_keys=True)
            elif candidate_id in set(old_direct.candidate_id):
                q = old_direct[old_direct.candidate_id == candidate_id]
                direct_status = "PASSED_PRIOR_BOUNDED_PF_TRANSFER"
                status.add("DIRECT_PF_VALIDATED")
                artifact_json = json.dumps(q[["temperature_K", "pf_2D_seed", "source_steps_file", "source_steps_sha256"]].to_dict("records"), sort_keys=True)
            aggregate = {
                "initial_onset_mean_MPa_sqrt_m": float(records.initial_onset_native_KJ_MPa_sqrt_m.mean()),
                "initial_onset_min_MPa_sqrt_m": float(records.initial_onset_native_KJ_MPa_sqrt_m.min()),
                "initial_onset_max_MPa_sqrt_m": float(records.initial_onset_native_KJ_MPa_sqrt_m.max()),
                "maximum_relative_deltaK_reinit": float(records.relative_deltaK_reinit.max()),
                "mean_N_reinit": float(records.N_reinit.mean()),
                "minimum_largest_avalanche_fraction": float(records.largest_avalanche_fraction.min()),
                "maximum_tip_radius_um": float(records.max_tip_radius_um.max()),
                "maximum_backstress_GPa": float(records.max_backstress_GPa.max()),
                "provider_robustness_status": "PROVIDER_SENSITIVE" if "PROVIDER_SENSITIVE" in roles[candidate_id] else "SCREENED_TRADEOFF",
            }
            bank_rows.append({
                "candidate_id": candidate_id,
                "candidate_version": BANK_VERSION,
                "target_response_class": material,
                "search_campaign_id": str(row.get("search_campaign_id", "")),
                "parent_or_anchor_id": str(row.get("parent_or_anchor_id", "")),
                "search_generation_method": str(row.get("search_generation_method", "")),
                "active_parameter_schema": "V913_ACTIVE_CANDIDATE_29",
                "active_parameter_count": len(ACTIVE_CANDIDATE_PARAMETER_FIELDS),
                "canonical_parameter_json": str(row.canonical_parameter_json),
                "parameter_sha256": digest,
                "source_commit": source_commit,
                "analysis_commit": producer_commit,
                "registry_commit": producer_commit,
                "qualified_domain_sha256": domain_hash,
                "option_status": "|".join(sorted(status)),
                "option_role": "|".join(roles[candidate_id]),
                "direct_PF_validation_status": direct_status,
                "direct_PF_artifacts_json": artifact_json,
                "fatigue_evaluated": False,
                "fatigue_validation_status": "NOT_EVALUATED",
                **aggregate,
                **{field: float(row[field]) for field in ACTIVE_CANDIDATE_PARAMETER_FIELDS},
            })
    bank = pd.DataFrame(bank_rows).sort_values(["target_response_class", "candidate_id"])
    material_columns = [
        "candidate_id", "candidate_version", "target_response_class",
        "search_campaign_id", "parent_or_anchor_id", "search_generation_method",
        "active_parameter_schema", "active_parameter_count", "canonical_parameter_json",
        "parameter_sha256", "source_commit", "analysis_commit", "registry_commit",
        *ACTIVE_CANDIDATE_PARAMETER_FIELDS,
    ]
    material_registry = bank[material_columns].drop_duplicates("parameter_sha256")
    material_registry.to_csv(OUT / "oneD_v2_fracture_option_bank_material_registry.csv", index=False)
    material_registry.to_parquet(OUT / "oneD_v2_fracture_option_bank_material_registry.parquet", index=False)
    bank.to_csv(OUT / "oneD_v2_fracture_option_bank.csv", index=False)
    bank.assign(response_records_json=bank.parameter_sha256.map(
        complete.set_index("parameter_sha256").response_records_json
    )).to_parquet(OUT / "oneD_v2_fracture_option_bank.parquet", index=False)
    pd.DataFrame([
        {"candidate_id": row.candidate_id, "parameter_sha256": row.parameter_sha256,
         "parameter_name": field, "value": float(row[field]), "units": units()[field]}
        for _, row in bank.iterrows() for field in ACTIVE_CANDIDATE_PARAMETER_FIELDS
    ]).to_csv(OUT / "oneD_v2_fracture_option_bank_long_parameters.csv", index=False)
    curated_ids = set(bank.candidate_id)
    response_features = response_all[response_all.candidate_id.isin(curated_ids)].copy()
    response_features.to_csv(OUT / "oneD_v2_fracture_option_bank_response_features.csv", index=False)

    shortlist_rows = []
    for material in CLASS_ORDER:
        chosen = selected_by_class[material][:SHORTLIST_COUNTS[material]]
        shortlist_rows.append(bank[(bank.target_response_class == material) & bank.candidate_id.isin(chosen)])
    shortlist = pd.concat(shortlist_rows, ignore_index=True)
    shortlist.to_csv(OUT / "oneD_v2_future_joint_search_shortlist.csv", index=False)

    final_registry = bank[bank.candidate_id.isin((CONTROLS["Peak"], CONTROLS["DBTT"]))].copy()
    final_registry["monotonic_fracture_decision"] = "RETAIN_CONTROL"
    final_registry.to_csv(OUT / "oneD_v2_peak_dbtt_rcurve_registry.csv", index=False)
    decision = {
        "schema": "oneD_v2_peak_dbtt_R_final_decision_v1",
        "Peak": {
            "decision": "RETAIN_CONTROL",
            "candidate_id": CONTROLS["Peak"],
            "variant_status": "NO_CREDIBLE_R_ENRICHED_ROW",
            "reason": "all reduced Peak-R signals were FEMCZM-only and direct PF remained one avalanche",
        },
        "DBTT": {
            "decision": "RETAIN_CONTROL",
            "candidate_id": CONTROLS["DBTT"],
            "variant_status": "NO_CREDIBLE_R_ENRICHED_ROW",
            "reason": "direct PF reinitiation onsets declined by about 9.35 MPa sqrt(m) at 1100-1200 K",
        },
        "FEMCZM_decision": "No FEM/CZM material row is changed; no new FEM/CZM simulation was run.",
        "focused_variants_preserved_in_option_bank": FOCUSED,
        "fatigue_evaluated": False,
        "fatigue_validation_status": "NOT_EVALUATED",
    }
    (OUT / "oneD_v2_peak_dbtt_R_final_decision.json").write_text(
        json.dumps(decision, indent=2, sort_keys=True) + "\n"
    )

    manifest_paths = [
        complete_path, pareto_path, OUT / "oneD_v2_fracture_option_bank.csv",
        OUT / "oneD_v2_fracture_option_bank.parquet",
        OUT / "oneD_v2_fracture_option_bank_material_registry.csv",
        OUT / "oneD_v2_fracture_option_bank_material_registry.parquet",
        OUT / "oneD_v2_fracture_option_bank_long_parameters.csv",
        OUT / "oneD_v2_fracture_option_bank_response_features.csv",
        OUT / "oneD_v2_future_joint_search_shortlist.csv",
        OUT / "oneD_v2_peak_dbtt_rcurve_registry.csv",
        OUT / "oneD_v2_peak_dbtt_R_final_decision.json",
    ]
    manifest = {
        "schema": "oneD_v2_fracture_option_bank_manifest_v1",
        "bank_version": BANK_VERSION,
        "source_commit": source_commit,
        "analysis_commit": producer_commit,
        "registry_commit": producer_commit,
        "qualified_domain_path": str(domain_path.relative_to(ROOT)),
        "qualified_domain_sha256": domain_hash,
        "candidate_counts": {
            "complete_unique_materials": len(complete),
            "pareto_and_near": len(pareto),
            "curated_by_class": bank.target_response_class.value_counts().to_dict(),
            "shortlist_by_class": shortlist.target_response_class.value_counts().to_dict(),
        },
        "material_parameter_fields": list(ACTIVE_CANDIDATE_PARAMETER_FIELDS),
        "material_parameter_units": units(),
        "forbidden_backend_fields": sorted(FORBIDDEN_MATERIAL_FIELDS),
        "backend_reduction_configuration": "SEPARATE_EXISTING_V2_PROVIDER_CONFIGURATION",
        "pareto_near_front_distance": 0.20,
        "diversity_method": "MAXIMIN_NORMALIZED_ACTIVE_PARAMETERS_PLUS_RESPONSE_FEATURES",
        "fatigue_evaluated": False,
        "fatigue_validation_status": "NOT_EVALUATED",
        "new_FEMCZM_runs": 0,
        "new_PF_runs": 6,
        "artifacts": {path.name: sha(path) for path in manifest_paths},
    }
    manifest_path = OUT / "oneD_v2_fracture_option_bank_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    provenance = {
        "schema": "oneD_v2_peak_dbtt_R_provenance_manifest_v1",
        **{key: manifest[key] for key in (
            "source_commit", "analysis_commit", "registry_commit",
            "qualified_domain_path", "qualified_domain_sha256",
            "new_FEMCZM_runs", "new_PF_runs",
        )},
        "focused_branch": "codex/oneD-v2-peak-dbtt-rcurve-search",
        "pf_branch": "codex/oneD-v2-peak-dbtt-R-pf-transfer",
        "pf_runner_commit": "130b305",
        "direct_pf_manifest_sha256": sha(OUT / "pf_2d_peak_dbtt_R_transfer_manifest.json"),
        "scientific_fingerprint": hashlib.sha256(json.dumps({
            "material_hashes": sorted(bank.parameter_sha256),
            "decision": decision,
            "domain": domain_hash,
        }, sort_keys=True, separators=(",", ":")).encode()).hexdigest(),
    }
    (OUT / "oneD_v2_peak_dbtt_R_provenance_manifest.json").write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n"
    )
    write_reports(bank, pareto, decision, manifest)
    write_figures(bank)
    print(
        f"FINALIZE_COMPLETE complete={len(complete)} pareto={len(pareto)} "
        f"bank={len(bank)} shortlist={len(shortlist)}"
    )
    return 0


def write_reports(bank: pd.DataFrame, pareto: pd.DataFrame, decision: dict, manifest: dict) -> None:
    reports = {
        "ONE_D_V2_PEAK_R_SEARCH.md": """# Peak R-propensity search\n\nThe focused screen evaluated 204 shared material vectors under both providers. No candidate produced reload-separated reinitiation in PF near 900–1000 K. Apparent enrichment was FEM/CZM-only and therefore failed the shared-row provider gate. The control is retained; provider-sensitive alternatives are preserved in the option bank.\n""",
        "ONE_D_V2_DBTT_R_SEARCH.md": """# DBTT R-propensity search\n\nThe screen found five strict canonical-seed candidates with a low-temperature shelf and upper-temperature reinitiation under both reduced providers. Three-seed validation showed residual low-temperature FEM/CZM seed sensitivity. The physically conservative finalist `oneD_v2_dbtt_R_9f5160f509e713e2` was advanced to direct PF.\n""",
        "ONE_D_V2_PEAK_DBTT_R_FINALIST_VALIDATION.md": """# Peak/DBTT R finalist validation\n\nAll 510 reduced validation cases reached their 100 or 300 µm target without map fallback. Peak enrichment remained backend-specific. DBTT upper-temperature enrichment survived all three reduced seeds but retained provider-scale and low-temperature seed sensitivity. Exact PF onset mechanics matched the native map; the DBTT wake-transition tensor factor differed by about 15%, requiring direct PF. No new FEM/CZM solve was run.\n""",
        "PF_2D_PEAK_R_TRANSFER_VALIDATION.md": """# Direct PF Peak-R transfer\n\nThe Peak-R finalist reached the 100 µm right-censor target at 600, 900, and 1200 K. Each case is one physical avalanche with no reload-separated reinitiation candidate. Event-wise PF native KJ variation is a model-native driving trajectory, not an R-curve. Transfer status: **failed for R enrichment**.\n""",
        "PF_2D_DBTT_R_TRANSFER_VALIDATION.md": """# Direct PF DBTT-R transfer\n\nThe DBTT-R finalist reached the 100 µm right-censor target at 600, 1100, and 1200 K. Each case has two reload-separated physical avalanches, but the second onset is lower: −0.44 MPa√m at 600 K and approximately −9.35 MPa√m at 1100–1200 K. Transfer status: **failed for rising resistance**.\n""",
        "ONE_D_V2_PEAK_DBTT_R_FINAL_DECISION.md": """# Final Peak/DBTT R decision\n\n**Peak: RETAIN_CONTROL; NO_CREDIBLE_R_ENRICHED_ROW.**\n\n**DBTT: RETAIN_CONTROL; NO_CREDIBLE_R_ENRICHED_ROW.**\n\nThe tested variants remain fracture-side diagnostic options, not promoted production rows. No FEM/CZM material row is changed; no new FEM/CZM simulation was run. Weak-T and ceramic-like selected rows remain unchanged. Fatigue was not evaluated.\n""",
        "ONE_D_V2_FRACTURE_OPTION_BANK.md": f"""# V2 fracture material option bank\n\nVersion: `{BANK_VERSION}`. The bank contains {len(bank)} curated full-precision shared material vectors: {bank.target_response_class.value_counts().to_dict()}. It preserves controls, non-dominated/near-front tradeoffs, provider-sensitive diagnostic variants, and prior weak-T/ceramic alternatives. Selection uses maximin distance in normalized material-plus-fracture-response space.\n\nThis is a **fatigue-ready material option** library only: every row has `fatigue_evaluated=false` and `fatigue_validation_status=NOT_EVALUATED`. No fatigue code or cyclic response was used. Backend lifecycle constants are excluded from the material registry.\n""",
        "ONE_D_V2_PEAK_DBTT_R_PROVENANCE.md": f"""# Peak/DBTT R provenance\n\nProducer/analysis/registry code commit: `{manifest['analysis_commit']}`. Predictive source baseline: `{manifest['source_commit']}`. Qualified domain SHA-256: `{manifest['qualified_domain_sha256']}`; the same hash is used in every bank row and manifest. Direct work comprised six PF finalist cases with at most two workers and zero new FEM/CZM runs.\n""",
    }
    for name, text in reports.items():
        (ROOT / name).write_text(text)


def write_figures(bank: pd.DataFrame) -> None:
    plt.rcParams.update({"figure.dpi": 140, "font.size": 8.5})
    reduced = pd.read_csv(OUT / "oneD_v2_peak_dbtt_R_multiseed_validation.csv")
    direct = pd.read_csv(OUT / "pf_2d_peak_dbtt_R_transfer_summary.csv")
    for material, stem in (("Peak", "PEAK"), ("DBTT", "DBTT")):
        ids = [CONTROLS[material], FOCUSED[material][0]]
        fig, ax = plt.subplots(figsize=(6.5, 4.2), constrained_layout=True)
        for candidate, style in zip(ids, ("--", "-")):
            q = reduced[(reduced.target_response_class == material) & (reduced.candidate_id == candidate) & (reduced.target_um == 100)]
            for provider, marker in (("PF", "o"), ("FEMCZM", "s")):
                z = q[q.provider == provider].groupby("temperature_K").agg(
                    initial=("initial_onset_native_KJ_MPa_sqrt_m", "mean"),
                    maximum=("maximum_onset_native_KJ_MPa_sqrt_m", "mean"),
                ).reset_index()
                ax.plot(z.temperature_K, z.initial, marker=marker, ls=style, label=f"{candidate.split('_')[-1]} {provider}")
                ax.fill_between(z.temperature_K, z.initial, z.maximum, alpha=.12)
        ax.set(xlabel="Temperature (K)", ylabel="Reload-separated onset native KJ (MPa√m)", title=f"{material}: onset and reinitiation envelope")
        ax.grid(alpha=.25); ax.legend(fontsize=6)
        fig.savefig(ROOT / f"{stem}_R_ONSET_AND_REINITIATION_ENVELOPE.png", dpi=180); plt.close(fig)

        fig, ax = plt.subplots(figsize=(6.5, 4.2), constrained_layout=True)
        for candidate, style in zip(ids, ("--", "-")):
            q = reduced[(reduced.target_response_class == material) & (reduced.candidate_id == candidate) & (reduced.target_um == 100)]
            z = q.groupby(["provider", "temperature_K"]).physical_avalanche_count.mean().reset_index()
            for provider, marker in (("PF", "o"), ("FEMCZM", "s")):
                a = z[z.provider == provider]
                ax.plot(a.temperature_K, a.physical_avalanche_count, marker=marker, ls=style, label=f"{candidate.split('_')[-1]} {provider}")
        ax.set(xlabel="Temperature (K)", ylabel="Physical-avalanche count", title=f"{material}: reduced avalanche topology")
        ax.grid(alpha=.25); ax.legend(fontsize=6)
        fig.savefig(ROOT / f"{stem}_R_AVALANCHE_TOPOLOGY.png", dpi=180); plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 4.5), constrained_layout=True)
    q = reduced[(reduced.target_um == 100) & reduced.candidate_id.isin([*CONTROLS.values(), FOCUSED["Peak"][0], FOCUSED["DBTT"][0]])]
    z = q.groupby(["target_response_class", "candidate_id", "provider"]).relative_deltaK_reinit.mean().reset_index()
    for key, group in z.groupby(["target_response_class", "candidate_id"]):
        ax.plot(group.provider, group.relative_deltaK_reinit, "o-", label=f"{key[0]} {key[1].split('_')[-1]}")
    ax.set(ylabel="Mean reduced relative ΔK reinit", title="Control versus R-screen variants (reduced)"); ax.grid(alpha=.25); ax.legend(fontsize=6)
    fig.savefig(ROOT / "CONTROL_VS_R_ENRICHED_REDUCED.png", dpi=180); plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 4.5), constrained_layout=True)
    for (material, role), group in direct.groupby(["material_class", "candidate_role"]):
        ax.plot(group.temperature_K, group.signed_max_reinitiation_minus_initial_K_MPa_sqrt_m.fillna(0), "o-", label=f"{material} {role}")
    ax.axhline(0, color="black", lw=.8); ax.set(xlabel="Temperature (K)", ylabel="Signed max reinitiation − initial K (MPa√m)", title="Direct PF reload-separated resistance candidates")
    ax.grid(alpha=.25); ax.legend(fontsize=7)
    fig.savefig(ROOT / "CONTROL_VS_R_ENRICHED_PF.png", dpi=180); plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 4.5), constrained_layout=True)
    q = bank[bank.target_response_class.isin(["Peak", "DBTT"])]
    for material, marker in (("Peak", "o"), ("DBTT", "s")):
        z = q[q.target_response_class == material]
        ax.scatter(z.maximum_tip_radius_um, z.maximum_backstress_GPa, marker=marker, label=material, alpha=.75)
    ax.set(xlabel="Maximum reduced tip radius (µm)", ylabel="Maximum reduced backstress (GPa)", title="Peak/DBTT curated process-zone states")
    ax.grid(alpha=.25); ax.legend()
    fig.savefig(ROOT / "PEAK_DBTT_R_PROCESS_ZONE_STATE.png", dpi=180); plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(9, 4), constrained_layout=True)
    counts = bank.target_response_class.value_counts().reindex(CLASS_ORDER)
    axes[0].bar(counts.index, counts.values); axes[0].tick_params(axis="x", rotation=25); axes[0].set(ylabel="Curated options", title="Fracture option bank")
    decisions = ["RETAIN\nCONTROL", "RETAIN\nCONTROL"]
    axes[1].bar(["Peak", "DBTT"], [1, 1], color=["#4C78A8", "#F58518"])
    for i, text in enumerate(decisions): axes[1].text(i, .5, text, ha="center", va="center", color="white", weight="bold")
    axes[1].set(yticks=[], title="Final monotonic-fracture decision")
    fig.savefig(ROOT / "PEAK_DBTT_R_FINAL_SUMMARY.png", dpi=180); plt.close(fig)


if __name__ == "__main__":
    raise SystemExit(main())
