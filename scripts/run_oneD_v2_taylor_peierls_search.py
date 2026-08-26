#!/usr/bin/env python3
"""Staged Taylor/Peierls-only kinetic screen for Peak and DBTT.

Stage 0 generates an empirical-bounds Sobol regime map.  Stage 1 evaluates a
diverse 512-row subset plus the immutable class control under both qualified
reduced mechanics providers.  All artifacts are checkpointed and fail closed
against the class control material contract.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
from typing import Any, Mapping

import numpy as np
import pandas as pd
from scipy.stats import qmc

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_outputs" / "oneD_v2_taylor_peierls_rcurve_search"
sys.path.insert(0, str(ROOT))

from arrhenius_fracture.emergent_gnd_campaign_v913 import candidate_from_registry_row
from arrhenius_fracture.emergent_gnd_contract_v913 import ACTIVE_CANDIDATE_PARAMETER_FIELDS
from arrhenius_fracture.zero_d_persistent_v913 import (
    ZeroDState,
    reduction_geometry,
    source_kinetic_diagnostics,
)
from reduced_fracture_v2.rcurve_propensity import summarize_rcurve_propensity
from reduced_fracture_v2.taylor_peierls_contract import (
    SEARCH_WHITELIST,
    barrier_hashes,
    full_material_hash,
    full_material_json,
    taylor_peierls_hash,
    validate_candidate,
)
from scripts.run_oneD_v2_predictive_campaign import inputs, run_case, summary
from scripts.search_oneD_v2_provider_robust import source_rows


CONTROL_IDS = {
    "Peak": "v913_zeroD_sobol_0242980",
    "DBTT": "v913_zeroD_sobol_0202500",
}
CLASS_SEEDS = {"Peak": 8666, "DBTT": 1008666}
TEMPERATURES = {
    "Peak": (600.0, 800.0, 900.0, 1000.0, 1100.0, 1200.0),
    "DBTT": (300.0, 600.0, 800.0, 900.0, 1000.0, 1100.0, 1200.0),
}
REPRESENTATIVE_TEMPERATURES = {
    "Peak": (600.0, 1000.0, 1200.0),
    "DBTT": (300.0, 1000.0, 1200.0),
}
SIGNED_FIELDS = {
    "peierls_activation_entropy_kB", "taylor_activation_entropy_kB",
}
CAMPAIGN_ID = "oneD_v2_taylor_peierls_rcurve_search_v1"
SOURCE_POPULATIONS = (
    ROOT / "candidates" / "v9_13_persistent_sites_top5_registry.csv",
    ROOT / "analysis_outputs" / "oneD_v2_peak_dbtt_rcurve_search"
    / "oneD_v2_complete_material_option_population.parquet",
)
CONTROL_REGISTRY = (
    ROOT / "analysis_outputs" / "oneD_v2_peak_dbtt_rcurve_search"
    / "oneD_v2_fracture_option_bank_material_registry.csv"
)


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _head() -> str:
    return subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()


def _controls() -> dict[str, pd.Series]:
    registry = pd.read_csv(CONTROL_REGISTRY)
    controls: dict[str, pd.Series] = {}
    for material, candidate_id in CONTROL_IDS.items():
        local = registry[registry.candidate_id.astype(str) == candidate_id]
        if len(local) != 1:
            raise RuntimeError(f"missing unique authoritative control {candidate_id}")
        controls[material] = local.iloc[0].copy()
    return controls


def _bound_source() -> pd.DataFrame:
    frames = [source_rows()]
    for path in SOURCE_POPULATIONS:
        frame = pd.read_csv(path) if path.suffix == ".csv" else pd.read_parquet(path)
        frames.append(frame)
    result = pd.concat(frames, ignore_index=True, sort=False)
    return result.dropna(subset=list(SEARCH_WHITELIST)).copy()


def transformed_bounds() -> dict[str, Any]:
    archived = _bound_source()
    controls = _controls()
    fields: dict[str, Any] = {}
    for field in SEARCH_WHITELIST:
        values = pd.to_numeric(archived[field], errors="coerce")
        values = values[np.isfinite(values)]
        if field not in SIGNED_FIELDS:
            values = values[values > 0.0]
        observed_min = float(values.min())
        observed_max = float(values.max())
        lower = float(values.quantile(0.05))
        upper = float(values.quantile(0.95))
        control_values = [float(controls[name][field]) for name in ("Peak", "DBTT")]
        lower = max(observed_min, min(lower, *control_values))
        upper = min(observed_max, max(upper, *control_values))
        if not lower < upper:
            raise RuntimeError(f"degenerate empirical bound for {field}")
        fields[field] = {
            "lower": lower,
            "upper": upper,
            "transform": "LINEAR_ADDITIVE" if field in SIGNED_FIELDS else "LOG10_POSITIVE",
            "observed_min": observed_min,
            "observed_max": observed_max,
            "empirical_quantiles": [0.05, 0.95],
        }
    return {
        "schema": "oneD_v2_taylor_peierls_empirical_bounds_v1",
        "fields": fields,
        "sources": [
            {"path": str(path.relative_to(ROOT)), "sha256": _sha(path)}
            for path in SOURCE_POPULATIONS
        ] + [{
            "path": "PF authoritative V9.13 material registries via source_rows()",
            "candidate_count_after_identity_deduplication": int(len(source_rows())),
        }],
        "policy": (
            "POOLED_OBSERVED_P05_P95_EXPANDED_TO_INCLUDE_BOTH_CONTROLS_"
            "AND_NEVER_BEYOND_OBSERVED_EXTREMA"
        ),
    }


def _coordinate(field: str, unit: float, definition: Mapping[str, Any]) -> float:
    lower, upper = float(definition["lower"]), float(definition["upper"])
    if definition["transform"] == "LINEAR_ADDITIVE":
        return lower + float(unit) * (upper - lower)
    return 10.0 ** (math.log10(lower) + float(unit) * (math.log10(upper) - math.log10(lower)))


def _candidate_identity(row: pd.Series, material: str, method: str) -> pd.Series:
    control_id = CONTROL_IDS[material]
    result = row.copy()
    digest = full_material_hash(result)
    if method != "HISTORICAL_CONTROL":
        result["candidate_id"] = f"oneD_v2_{material.lower()}_TP_{digest[:16]}"
    result["target_class"] = material
    result["parent_control_id"] = control_id
    result["search_stage"] = "STAGE0_KINETIC_REGIME_MAP"
    result["search_method"] = method
    result["search_campaign_id"] = CAMPAIGN_ID
    result["canonical_parameter_json"] = full_material_json(result)
    result["full_material_sha256"] = digest
    result["taylor_peierls_subvector_sha256"] = taylor_peierls_hash(result)
    for field, value in barrier_hashes(result).items():
        result[field] = value
    return result


def sobol_population(material: str, power: int, bounds: Mapping[str, Any]) -> pd.DataFrame:
    controls = _controls()
    control = controls[material]
    records = [_candidate_identity(control, material, "HISTORICAL_CONTROL")]
    sampler = qmc.Sobol(
        d=len(SEARCH_WHITELIST), scramble=True,
        seed=551903 if material == "Peak" else 551907,
    )
    for vector in sampler.random_base2(power):
        row = control.copy()
        for field, unit in zip(SEARCH_WHITELIST, vector):
            row[field] = _coordinate(field, float(unit), bounds["fields"][field])
        validate_candidate(row, control)
        records.append(_candidate_identity(row, material, "SOBOL_EMPIRICAL_BOUNDS"))
    frame = pd.DataFrame(records).drop_duplicates("full_material_sha256", keep="first")
    for _, row in frame.iterrows():
        validate_candidate(row, control)
    return frame.sort_values(["search_method", "candidate_id"], kind="stable").reset_index(drop=True)


def _representative_states(material: str) -> list[dict[str, Any]]:
    physics, _, providers = inputs()
    control = _controls()[material]
    mechanics, drive = providers["PF"]
    records: list[dict[str, Any]] = []
    for temperature in REPRESENTATIVE_TEMPERATURES[material]:
        result = run_case(
            control, material, temperature, "PF", mechanics, drive, physics,
            100.0, seed=CLASS_SEEDS[material], maximum_intervals=10_000,
        )
        if not result["events"]:
            raise RuntimeError(f"control has no representative event at {temperature:g} K")
        event = result["events"][0]
        nsys = int(physics.n_systems)
        slip = max(
            (float(event["tip_radius_m"]) - float(physics.r0_m))
            / max(
                float(control["c_blunt"]) * abs(float(physics.b_m))
                * float(physics.blunting_slip_fraction),
                1.0e-300,
            ),
            0.0,
        )
        records.append({
            "temperature_K": float(temperature),
            "K_MPa_sqrt_m": float(event["native_KJ_MPa_sqrt_m"]),
            "drive_factors": [float(value) for value in event["drive_factors"]],
            "extension_m": float(event["extension_before_m"]),
            "mobile_density_m2": float(event["mobile_density_m2"]),
            "retained_density_m2": float(event["retained_density_m2"]),
            "slip_count": slip,
            "state_source": "PF_REDUCED_CONTROL_FIRST_PRE_EVENT_STATE",
        })
    return records


def regime_map(candidates: pd.DataFrame, material: str) -> pd.DataFrame:
    physics, _, _ = inputs()
    reduced = reduction_geometry(physics)
    control = _controls()[material]
    states = _representative_states(material)
    records: list[dict[str, Any]] = []
    for index, row in candidates.iterrows():
        validate_candidate(row, control)
        model = candidate_from_registry_row(row)
        local: list[dict[str, float]] = []
        for representative in states:
            nsys = int(physics.n_systems)
            state = ZeroDState(
                mobile_m2=np.full(nsys, representative["mobile_density_m2"] / nsys),
                retained_m2=np.full(nsys, representative["retained_density_m2"] / nsys),
                local_slip_count_by_system=np.full(nsys, representative["slip_count"] / nsys),
                extension_m=representative["extension_m"],
                cumulative_activations=np.zeros(nsys),
            )
            diagnostic = source_kinetic_diagnostics(
                model, physics, reduced, state,
                K_MPa_sqrt_m=representative["K_MPa_sqrt_m"],
                temperature_K=representative["temperature_K"],
                drive_factors_override=representative["drive_factors"],
            )
            def median(name: str) -> float:
                return float(np.median(np.asarray(diagnostic[name], dtype=float)))
            local.append({
                "temperature_K": representative["temperature_K"],
                "log10_peierls_rate_s": math.log10(max(median("peierls_rate_s_by_system"), 1e-300)),
                "log10_peierls_velocity_m_s": math.log10(max(median("peierls_velocity_m_s_by_system"), 1e-300)),
                "log10_encounter_rate_s": math.log10(max(median("encounter_rate_s_by_system"), 1e-300)),
                "log10_taylor_completion_rate_s": math.log10(max(median("taylor_completion_rate_s_by_system"), 1e-300)),
                "log10_chi_ret": math.log10(max(median("chi_ret_by_system"), 1e-300)),
                "log10_chi_taylor_completion": math.log10(max(median("chi_taylor_completion_by_system"), 1e-300)),
                "retained_equilibrium_fraction": median("retained_equilibrium_fraction_by_system"),
            })
        item = {
            "candidate_id": str(row.candidate_id),
            "target_class": material,
            "full_material_sha256": str(row.full_material_sha256),
            "taylor_peierls_subvector_sha256": str(row.taylor_peierls_subvector_sha256),
            "diagnostic_states_json": json.dumps(local, sort_keys=True, separators=(",", ":")),
        }
        for field in local[0]:
            if field == "temperature_K":
                continue
            values = np.asarray([record[field] for record in local], dtype=float)
            item[f"{field}_median"] = float(np.median(values))
            item[f"{field}_min"] = float(np.min(values))
            item[f"{field}_max"] = float(np.max(values))
        log_chi = float(item["log10_chi_ret_median"])
        equilibrium = float(item["retained_equilibrium_fraction_median"])
        if log_chi < -1.0 or equilibrium < 0.10:
            regime = "TRANSPORT_DOMINATED"
        elif log_chi > 1.0 and equilibrium > 0.90:
            regime = "RETENTION_DOMINATED"
        else:
            regime = "BALANCED_TRANSPORT_RETENTION"
        item["kinetic_regime"] = regime
        records.append(item)
        if (index + 1) % 512 == 0:
            print(f"STAGE0_PROGRESS class={material} candidates={index + 1}/{len(candidates)}", flush=True)
    return pd.DataFrame(records)


def _normalized_selection_features(
    candidates: pd.DataFrame, diagnostics: pd.DataFrame, bounds: Mapping[str, Any]
) -> np.ndarray:
    columns: list[np.ndarray] = []
    for field in SEARCH_WHITELIST:
        values = candidates[field].to_numpy(float)
        definition = bounds["fields"][field]
        if definition["transform"] == "LOG10_POSITIVE":
            values = np.log10(values)
            lower, upper = math.log10(definition["lower"]), math.log10(definition["upper"])
        else:
            lower, upper = definition["lower"], definition["upper"]
        columns.append((values - lower) / max(upper - lower, 1.0e-300))
    aligned = diagnostics.set_index("candidate_id").loc[candidates.candidate_id]
    for field in (
        "log10_peierls_rate_s_median", "log10_taylor_completion_rate_s_median",
        "log10_chi_ret_median", "log10_chi_taylor_completion_median",
        "retained_equilibrium_fraction_median",
    ):
        values = aligned[field].to_numpy(float)
        lower, upper = np.nanpercentile(values, [1.0, 99.0])
        columns.append(np.clip((values - lower) / max(upper - lower, 1.0e-300), 0.0, 1.0))
    return np.column_stack(columns)


def diverse_selection(
    candidates: pd.DataFrame, diagnostics: pd.DataFrame, bounds: Mapping[str, Any], count: int
) -> pd.DataFrame:
    control_mask = candidates.search_method.eq("HISTORICAL_CONTROL").to_numpy()
    pool = candidates[~control_mask].reset_index(drop=True)
    features = _normalized_selection_features(pool, diagnostics, bounds)
    aligned = diagnostics.set_index("candidate_id").loc[pool.candidate_id]
    forced: set[int] = set()
    for field in (
        "log10_peierls_rate_s_median", "log10_taylor_completion_rate_s_median",
        "log10_chi_ret_median", "log10_chi_taylor_completion_median",
        "retained_equilibrium_fraction_median",
    ):
        values = aligned[field].to_numpy(float)
        forced.update((int(np.argmin(values)), int(np.argmax(values))))
    for regime in sorted(aligned.kinetic_regime.unique()):
        indices = np.flatnonzero(aligned.kinetic_regime.to_numpy() == regime)
        if len(indices):
            forced.add(int(indices[len(indices) // 2]))
    selected = sorted(forced)
    distance = np.full(len(pool), np.inf)
    if selected:
        distance = np.min(
            np.sum((features[:, None, :] - features[np.asarray(selected)][None, :, :]) ** 2, axis=2),
            axis=1,
        )
    distance[np.asarray(selected, dtype=int)] = -np.inf
    while len(selected) < count:
        index = int(np.argmax(distance))
        selected.append(index)
        candidate_distance = np.sum((features - features[index]) ** 2, axis=1)
        distance = np.minimum(distance, candidate_distance)
        distance[np.asarray(selected, dtype=int)] = -np.inf
    selected_frame = pool.iloc[selected].copy()
    selected_frame["stage1_selection_rank"] = np.arange(1, len(selected_frame) + 1)
    control = candidates[control_mask].copy()
    control["stage1_selection_rank"] = 0
    result = pd.concat([control, selected_frame], ignore_index=True)
    return result.sort_values("stage1_selection_rank", kind="stable").reset_index(drop=True)


def prepare(power: int, selected_count: int) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    bounds = transformed_bounds()
    (OUT / "oneD_v2_taylor_peierls_search_bounds.json").write_text(
        json.dumps(bounds, indent=2, sort_keys=True) + "\n"
    )
    all_candidates: list[pd.DataFrame] = []
    all_diagnostics: list[pd.DataFrame] = []
    all_selected: list[pd.DataFrame] = []
    for material in ("Peak", "DBTT"):
        candidates = sobol_population(material, power, bounds)
        diagnostics = regime_map(candidates, material)
        selected = diverse_selection(candidates, diagnostics, bounds, selected_count)
        all_candidates.append(candidates)
        all_diagnostics.append(diagnostics)
        all_selected.append(selected)
        print(
            f"STAGE0_COMPLETE class={material} generated={len(candidates)} selected={len(selected)}",
            flush=True,
        )
    population = pd.concat(all_candidates, ignore_index=True)
    diagnostics = pd.concat(all_diagnostics, ignore_index=True)
    selected = pd.concat(all_selected, ignore_index=True)
    population.to_parquet(OUT / "oneD_v2_taylor_peierls_stage0_population.parquet", index=False)
    diagnostics.to_csv(OUT / "oneD_v2_taylor_peierls_timescale_diagnostics.csv", index=False)
    selected.to_csv(OUT / "oneD_v2_taylor_peierls_stage1_selection.csv", index=False)
    manifest = {
        "schema": "oneD_v2_taylor_peierls_stage0_v1",
        "producer_commit": _head(),
        "campaign_id": CAMPAIGN_ID,
        "sobol_power": power,
        "sobol_rows_per_class": 2 ** power,
        "control_rows_per_class": 1,
        "selected_sobol_rows_per_class": selected_count,
        "selected_total_rows_per_class": selected_count + 1,
        "whitelist": list(SEARCH_WHITELIST),
        "bounds_sha256": _sha(OUT / "oneD_v2_taylor_peierls_search_bounds.json"),
        "population_sha256": _sha(OUT / "oneD_v2_taylor_peierls_stage0_population.parquet"),
        "diagnostics_sha256": _sha(OUT / "oneD_v2_taylor_peierls_timescale_diagnostics.csv"),
        "selection_sha256": _sha(OUT / "oneD_v2_taylor_peierls_stage1_selection.csv"),
        "new_FEMCZM_runs": 0,
    }
    (OUT / "oneD_v2_taylor_peierls_stage0_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )


_WORKER_INPUTS: tuple[Any, Any, Any] | None = None


def _worker_initialize() -> None:
    global _WORKER_INPUTS
    physics, _, providers = inputs()
    _WORKER_INPUTS = (physics, providers, _controls())


def _case_task(task: tuple[dict[str, Any], str, str, float, float, int, str]) -> dict[str, Any]:
    if _WORKER_INPUTS is None:
        _worker_initialize()
    assert _WORKER_INPUTS is not None
    physics, providers, controls = _WORKER_INPUTS
    row_dict, material, provider, temperature, target_um, seed, stage = task
    row = pd.Series(row_dict)
    validate_candidate(row, controls[material])
    base = {
        "candidate_id": str(row.candidate_id), "target_class": material,
        "parent_control_id": str(row.parent_control_id), "provider": provider,
        "temperature_K": float(temperature), "hazard_seed": int(seed),
        "target_um": float(target_um), "search_stage": stage,
        "search_method": str(row.search_method),
        "canonical_parameter_json": str(row.canonical_parameter_json),
        "full_material_sha256": str(row.full_material_sha256),
        "cleavage_barrier_sha256": str(row.cleavage_barrier_sha256),
        "emission_barrier_sha256": str(row.emission_barrier_sha256),
        "taylor_peierls_subvector_sha256": str(row.taylor_peierls_subvector_sha256),
    }
    for field in ACTIVE_CANDIDATE_PARAMETER_FIELDS:
        base[field] = float(row[field])
    mechanics, drive = providers[provider]
    try:
        result = run_case(
            row, material, temperature, provider, mechanics, drive, physics,
            target_um, seed=seed, maximum_intervals=25_000,
        )
        base.update(summary(result))
        base.update(summarize_rcurve_propensity(result))
        onsets = json.loads(base["onset_candidates_json"])
        base.update({
            "maximum_onset_tip_radius_um": max((x["tip_radius_m"] * 1e6 for x in onsets), default=np.nan),
            "maximum_onset_backstress_GPa": max((x["backstress_Pa"] * 1e-9 for x in onsets), default=np.nan),
            "maximum_onset_mobile_density_m2": max((x["mobile_density_m2"] for x in onsets), default=np.nan),
            "maximum_onset_retained_density_m2": max((x["retained_density_m2"] for x in onsets), default=np.nan),
            "maximum_onset_chi_ret": max((x.get("chi_ret_max", np.nan) for x in onsets), default=np.nan),
            "maximum_onset_chi_taylor_completion": max((x.get("chi_taylor_completion_max", np.nan) for x in onsets), default=np.nan),
            "maximum_onset_retained_equilibrium_fraction": max((x.get("retained_equilibrium_fraction_max", np.nan) for x in onsets), default=np.nan),
            "terminal_reason": str(result["status"]), "case_error": "",
            "map_oracle_fallback_count": int(result["status"] == "RIGHT_CENSORED_DRIVE_MAP_BOUND"),
        })
    except Exception as exc:
        base.update({
            "status": "EVALUATION_FAILURE", "terminal_reason": type(exc).__name__,
            "case_error": str(exc), "event_count": 0, "physical_avalanche_count": 0,
            "precursor_reinitiation_count": 0, "N_reinit": 0,
            "largest_avalanche_fraction": 0.0, "map_oracle_fallback_count": 0,
        })
    return base


def screen(workers: int, checkpoint_every: int, executor_kind: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    selection = pd.read_csv(OUT / "oneD_v2_taylor_peierls_stage1_selection.csv")
    checkpoint = OUT / "oneD_v2_taylor_peierls_stage1_100um_checkpoint.parquet"
    records = pd.read_parquet(checkpoint).to_dict("records") if checkpoint.exists() else []
    done = {
        (str(x["candidate_id"]), str(x["provider"]), float(x["temperature_K"]))
        for x in records
    }
    tasks = []
    for _, row in selection.iterrows():
        material = str(row.target_class)
        for provider in ("PF", "FEMCZM"):
            for temperature in TEMPERATURES[material]:
                key = (str(row.candidate_id), provider, float(temperature))
                if key not in done:
                    tasks.append((row.to_dict(), material, provider, temperature, 100.0,
                                  CLASS_SEEDS[material], "STAGE1_TWO_PROVIDER_100UM"))
    print(
        f"STAGE1_START pending={len(tasks)} existing={len(records)} "
        f"workers={workers} executor={executor_kind}", flush=True,
    )
    if executor_kind == "process":
        executor_context = ProcessPoolExecutor(
            max_workers=workers, initializer=_worker_initialize
        )
    else:
        _worker_initialize()
        executor_context = ThreadPoolExecutor(max_workers=workers)
    with executor_context as executor:
        futures = [executor.submit(_case_task, task) for task in tasks]
        for completed, future in enumerate(as_completed(futures), start=1):
            records.append(future.result())
            if completed % checkpoint_every == 0 or completed == len(futures):
                frame = pd.DataFrame(records).sort_values(
                    ["target_class", "candidate_id", "provider", "temperature_K"], kind="stable"
                )
                frame.to_parquet(checkpoint, index=False)
                print(f"STAGE1_PROGRESS completed={completed}/{len(futures)} total={len(frame)}", flush=True)
    cases = pd.DataFrame(records).sort_values(
        ["target_class", "candidate_id", "provider", "temperature_K"], kind="stable"
    ).reset_index(drop=True)
    target = OUT / "oneD_v2_taylor_peierls_stage1_100um_cases.parquet"
    cases.to_parquet(target, index=False)
    manifest = {
        "schema": "oneD_v2_taylor_peierls_stage1_100um_v1",
        "producer_commit": _head(), "campaign_id": CAMPAIGN_ID,
        "case_count": int(len(cases)), "candidate_count": int(cases.candidate_id.nunique()),
        "providers": ["PF", "FEMCZM"], "temperatures_K": TEMPERATURES,
        "canonical_seeds": CLASS_SEEDS, "target_um": 100.0,
        "reduced_worker_count": int(workers), "reduced_executor": executor_kind,
        "heavy_PF_worker_count": 0,
        "results_sha256": _sha(target), "new_FEMCZM_runs": 0,
        "status_counts": cases.status.value_counts(dropna=False).to_dict(),
    }
    (OUT / "oneD_v2_taylor_peierls_stage1_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    checkpoint.unlink(missing_ok=True)
    print(f"STAGE1_COMPLETE cases={len(cases)}", flush=True)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("prepare", "screen", "all"))
    parser.add_argument("--sobol-power", type=int, default=12)
    parser.add_argument("--selected-per-class", type=int, default=512)
    parser.add_argument("--workers", type=int, default=min(max(os.cpu_count() or 2, 2), 8))
    parser.add_argument("--executor", choices=("thread", "process"), default="thread")
    parser.add_argument("--checkpoint-every", type=int, default=128)
    args = parser.parse_args()
    if args.mode in {"prepare", "all"}:
        prepare(args.sobol_power, args.selected_per_class)
    if args.mode in {"screen", "all"}:
        screen(args.workers, args.checkpoint_every, args.executor)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
