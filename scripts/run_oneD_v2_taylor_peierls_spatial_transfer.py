#!/usr/bin/env python3
"""Screen DBTT Taylor/Peierls rows with the spatial v9.13 state engine."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
from typing import Any

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SOURCE_OUT = ROOT / "analysis_outputs" / "oneD_v2_taylor_peierls_rcurve_search"
OUT = ROOT / "analysis_outputs" / "oneD_v2_taylor_peierls_spatial_transfer"
PF_TRANSFER_BASE = SOURCE_OUT / "oneD_v2_taylor_peierls_pf_transfer_registry.csv"
sys.path.insert(0, str(ROOT))

from arrhenius_fracture.emergent_gnd_campaign_v913 import candidate_from_registry_row
from arrhenius_fracture.emergent_gnd_contract_v913 import ACTIVE_CANDIDATE_PARAMETER_FIELDS
from arrhenius_fracture.emergent_gnd_rcurve_v913 import run_autonomous_rcurve
from reduced_fracture_v2.predictive import (
    PF_PREDICTIVE_LIFECYCLE_V2,
    stochastic_loading_map,
)
from reduced_fracture_v2.spatial_transfer import (
    ObservedEmergentGNDState,
    physical_onset_event_indices,
    profile_rows,
    sampled_history,
    state_snapshot,
)
from reduced_fracture_v2.taylor_peierls_contract import (
    SEARCH_WHITELIST,
    barrier_hashes,
    full_material_hash,
    validate_candidate,
)
from scripts.run_oneD_v2_predictive_campaign import inputs


CONTROL_ID = "v913_zeroD_sobol_0202500"
TEMPERATURE_K = 1100.0
HAZARD_SEED = 1008666
TARGET_UM = 300.0
CHECKPOINT_SUMMARY = OUT / "oneD_v2_spatial_1d_candidate_summary_checkpoint.parquet"
CHECKPOINT_HISTORY = OUT / "oneD_v2_spatial_1d_state_history_checkpoint.parquet"
SUMMARY_PATH = OUT / "oneD_v2_spatial_1d_candidate_summary.parquet"
HISTORY_PATH = OUT / "oneD_v2_spatial_1d_state_history.parquet"


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _head() -> str:
    return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()


_WORKER: dict[str, Any] | None = None


def _worker_init() -> None:
    global _WORKER
    physics, _, providers = inputs()
    mechanics, drive = providers["PF"]
    population = pd.read_parquet(SOURCE_OUT / "oneD_v2_taylor_peierls_complete_population.parquet")
    control = population[population.candidate_id.eq(CONTROL_ID)].iloc[0]
    loading = stochastic_loading_map(
        mechanics,
        seed=HAZARD_SEED,
        target_extension_m=TARGET_UM * 1.0e-6,
    )
    _WORKER = {
        "physics": physics,
        "mechanics": mechanics,
        "drive": drive,
        "loading": loading,
        "control": control,
    }


def _run_candidate(row_data: dict[str, Any], *, save_profiles: bool = False) -> dict[str, Any]:
    if _WORKER is None:
        _worker_init()
    assert _WORKER is not None
    row = pd.Series(row_data)
    validate_candidate(row, _WORKER["control"])
    snapshots: list[dict[str, Any]] = []

    def observer(phase, context, state):
        record = state_snapshot(
            phase, context, state, temperature_K=TEMPERATURE_K
        )
        if save_profiles:
            record["_state_payload"] = {
                "x_m": np.asarray(state.x, dtype=float).copy(),
                "mobile_m2": np.asarray(state.mobile_m2, dtype=float).copy(),
                "retained_m2": np.asarray(state.retained_m2, dtype=float).copy(),
                "accumulated_slip_m2": np.asarray(
                    state.accumulated_slip_m2, dtype=float
                ).copy(),
            }
        snapshots.append(record)

    drive = _WORKER["drive"]
    result = run_autonomous_rcurve(
        candidate_from_registry_row(row),
        _WORKER["physics"],
        _WORKER["loading"],
        TEMPERATURE_K,
        target_projected_extension_m=TARGET_UM * 1.0e-6,
        maximum_applied_displacement_m=1.0e-3,
        maximum_integration_substeps=750_000,
        translation_mode="hazard_coupled",
        emission_drive_provider=lambda extension, state: drive.evaluate(
            extension, state.tip_radius_m()
        ),
        state_factory=ObservedEmergentGNDState,
        state_observer=observer,
    )
    onsets, avalanches = physical_onset_event_indices(
        result.events,
        reload_gap_threshold_m=PF_PREDICTIVE_LIFECYCLE_V2.reload_gap_threshold_m,
    )
    sampled = sampled_history(snapshots, onsets) if result.status == "complete" else []
    onset_K = [float(result.events[index].K_MPa_sqrt_m) for index in onsets]
    counts = np.bincount(np.asarray(avalanches, dtype=int)) if avalanches else np.zeros(0)
    summary = {
        "candidate_id": str(row.candidate_id),
        "full_material_sha256": str(row.full_material_sha256),
        "taylor_peierls_subvector_sha256": str(row.taylor_peierls_subvector_sha256),
        "cleavage_barrier_sha256": str(row.cleavage_barrier_sha256),
        "emission_barrier_sha256": str(row.emission_barrier_sha256),
        "temperature_K": TEMPERATURE_K,
        "hazard_seed": HAZARD_SEED,
        "target_um": TARGET_UM,
        "status": result.status,
        "event_count": len(result.events),
        "physical_avalanche_count": int(len(counts)),
        "N_reinit": max(len(onsets) - 1, 0),
        "largest_avalanche_fraction": (
            float(np.max(counts) / max(len(result.events), 1)) if len(counts) else 0.0
        ),
        "first_onset_K_MPa_sqrt_m": onset_K[0] if onset_K else np.nan,
        "maximum_onset_K_MPa_sqrt_m": max(onset_K, default=np.nan),
        "deltaK_reinit_MPa_sqrt_m": (
            max(onset_K[1:]) - onset_K[0] if len(onset_K) > 1 else np.nan
        ),
        "achieved_projected_extension_um": result.achieved_projected_extension_m * 1.0e6,
        "maximum_tip_radius_um": result.max_tip_radius_m * 1.0e6,
        "maximum_backstress_GPa": result.max_backstress_Pa * 1.0e-9,
        "minimum_front_width_um": result.min_front_width_m * 1.0e6,
        "maximum_source_multiplicity": result.max_source_multiplicity,
        "final_applied_opening_um": result.final_applied_displacement_m * 1.0e6,
        "final_time_s": result.final_elapsed_time_s,
        "spatial_state_model": result.numerical_integration["model_id"],
        "observer_feedback": False,
        "new_FEMCZM_runs": 0,
    }
    history: list[dict[str, Any]] = []
    profiles: list[dict[str, Any]] = []
    for sample in sampled:
        payload = sample.pop("_state_payload", None)
        clean = {"candidate_id": str(row.candidate_id), **sample}
        history.append(clean)
        if payload is not None:
            profiles.extend(profile_rows(str(row.candidate_id), clean, payload))
    return {"summary": summary, "history": history, "profiles": profiles}


def _task(row_data: dict[str, Any]) -> dict[str, Any]:
    try:
        return _run_candidate(row_data)
    except Exception as exc:
        return {
            "summary": {
                "candidate_id": str(row_data["candidate_id"]),
                "status": "EVALUATION_FAILURE",
                "case_error": f"{type(exc).__name__}: {exc}",
                "temperature_K": TEMPERATURE_K,
                "hazard_seed": HAZARD_SEED,
                "target_um": TARGET_UM,
            },
            "history": [],
            "profiles": [],
        }


def _screen_pool() -> pd.DataFrame:
    pareto = pd.read_parquet(SOURCE_OUT / "oneD_v2_taylor_peierls_pareto_library.parquet")
    pool = pareto[pareto.target_class.eq("DBTT")].copy()
    if len(pool) != 513 or not pool.candidate_id.is_unique:
        raise RuntimeError("expected 513 unique fully evaluated DBTT rows")
    return pool.sort_values("candidate_id", kind="stable").reset_index(drop=True)


def run_screen(workers: int, checkpoint_every: int) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    pool = _screen_pool()
    summaries = pd.read_parquet(CHECKPOINT_SUMMARY).to_dict("records") if CHECKPOINT_SUMMARY.exists() else []
    histories = pd.read_parquet(CHECKPOINT_HISTORY).to_dict("records") if CHECKPOINT_HISTORY.exists() else []
    done = {str(row["candidate_id"]) for row in summaries}
    tasks = [row.to_dict() for _, row in pool.iterrows() if str(row.candidate_id) not in done]
    print(f"SPATIAL_SCREEN_START pending={len(tasks)} existing={len(done)} workers={workers}", flush=True)
    with ProcessPoolExecutor(max_workers=workers, initializer=_worker_init) as executor:
        futures = [executor.submit(_task, task) for task in tasks]
        for completed, future in enumerate(as_completed(futures), start=1):
            item = future.result()
            summaries.append(item["summary"])
            histories.extend(item["history"])
            if completed % checkpoint_every == 0 or completed == len(futures):
                pd.DataFrame(summaries).sort_values("candidate_id", kind="stable").to_parquet(
                    CHECKPOINT_SUMMARY, index=False
                )
                pd.DataFrame(histories).sort_values(
                    ["candidate_id", "sample_role", "requested_extension_um"], kind="stable"
                ).to_parquet(CHECKPOINT_HISTORY, index=False)
                print(f"SPATIAL_SCREEN_PROGRESS completed={completed}/{len(futures)}", flush=True)
    summary = pd.DataFrame(summaries).sort_values("candidate_id", kind="stable").reset_index(drop=True)
    history = pd.DataFrame(histories).sort_values(
        ["candidate_id", "sample_role", "requested_extension_um"], kind="stable"
    ).reset_index(drop=True)
    if len(summary) != len(pool) or not summary.candidate_id.is_unique:
        raise RuntimeError("spatial screen is incomplete or has duplicate identities")
    summary.to_parquet(SUMMARY_PATH, index=False)
    history.to_parquet(HISTORY_PATH, index=False)
    CHECKPOINT_SUMMARY.unlink(missing_ok=True)
    CHECKPOINT_HISTORY.unlink(missing_ok=True)
    print(f"SPATIAL_SCREEN_COMPLETE candidates={len(summary)} history_rows={len(history)}", flush=True)


FEATURE_METRICS = (
    "mobile_line_content", "retained_line_content", "retained_fraction",
    "mobile_centroid_tip_relative_m", "mobile_width_m",
    "retained_centroid_tip_relative_m", "retained_width_m",
    "near_tip_mobile_line_content", "near_tip_retained_line_content",
    "wake_mobile_departed_line_content", "wake_retained_departed_line_content",
    "wake_retained_centroid_laboratory_m", "wake_retained_width_m",
    "tip_radius_m", "front_width_m", "backstress_Pa", "K_shield_MPa_sqrt_m",
    "source_multiplicity", "peierls_transport_velocity_abs_max_m_s",
    "peierls_transport_velocity_abs_mean_m_s", "encounter_retention_rate_max_s",
    "encounter_retention_rate_mean_s", "taylor_completion_rate_max_s",
    "taylor_completion_rate_mean_s",
)


def build_feature_vectors(history: pd.DataFrame, summary: pd.DataFrame) -> pd.DataFrame:
    checkpoint = history[history.sample_role.eq("EXTENSION_CHECKPOINT")]
    records = []
    for candidate_id, local in checkpoint.groupby("candidate_id", sort=False):
        record: dict[str, Any] = {"candidate_id": candidate_id}
        for _, row in local.iterrows():
            label = f"x{int(round(float(row.requested_extension_um))):03d}um"
            for metric in FEATURE_METRICS:
                record[f"{label}__{metric}"] = float(row[metric])
        onsets = history[
            history.candidate_id.eq(candidate_id) & history.sample_role.eq("PHYSICAL_ONSET")
        ]
        for metric in FEATURE_METRICS:
            record[f"onset_max__{metric}"] = float(onsets[metric].max())
            record[f"onset_final__{metric}"] = float(onsets.iloc[-1][metric])
        records.append(record)
    features = pd.DataFrame(records)
    features = features.merge(summary, on="candidate_id", how="left", validate="one_to_one")
    control_onset = float(features.loc[features.candidate_id.eq(CONTROL_ID), "first_onset_K_MPa_sqrt_m"].iloc[0])
    features["first_onset_relative_error_to_control"] = (
        features.first_onset_K_MPa_sqrt_m - control_onset
    ) / control_onset
    features["selection_eligible"] = (
        features.status.eq("complete")
        & features.achieved_projected_extension_um.ge(TARGET_UM - 1.0e-9)
        & features.first_onset_relative_error_to_control.abs().le(0.20)
    )
    return features


def _standardized_matrix(features: pd.DataFrame) -> tuple[np.ndarray, list[str]]:
    columns = [column for column in features if "__" in column]
    values = features[columns].to_numpy(float).copy()
    for index, column in enumerate(columns):
        if (
            "line_content" in column or "multiplicity" in column
            or "rate_" in column or "velocity_" in column
        ):
            values[:, index] = np.sign(values[:, index]) * np.log10(1.0 + np.abs(values[:, index]))
    median = np.nanmedian(values, axis=0)
    q25, q75 = np.nanpercentile(values, [25.0, 75.0], axis=0)
    scale = np.where(q75 > q25, q75 - q25, np.nanstd(values, axis=0))
    scale = np.where(scale > 0.0, scale, 1.0)
    matrix = (values - median) / scale
    if np.any(~np.isfinite(matrix)):
        raise RuntimeError("nonfinite spatial response feature matrix")
    return matrix, columns


def _pick_extreme(
    eligible: pd.DataFrame, field: str, selected: set[str], *, absolute: bool = False
) -> str:
    values = eligible[field].abs() if absolute else eligible[field]
    order = values.sort_values(ascending=False, kind="stable").index
    for index in order:
        candidate = str(eligible.loc[index, "candidate_id"])
        if candidate not in selected:
            return candidate
    raise RuntimeError(f"no unselected candidate remains for {field}")


def select_candidates(features: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    eligible = features[features.selection_eligible].copy().reset_index(drop=True)
    if CONTROL_ID not in set(eligible.candidate_id):
        raise RuntimeError("control failed the spatial selection gates")
    matrix, columns = _standardized_matrix(eligible)
    distances = np.sqrt(np.maximum(
        np.sum((matrix[:, None, :] - matrix[None, :, :]) ** 2, axis=2), 0.0
    ))
    selected: set[str] = {CONTROL_ID}
    roles = [("CONTROL", CONTROL_ID)]
    role_fields = (
        ("TRANSPORT_EXTREME", "onset_max__peierls_transport_velocity_abs_max_m_s", False),
        ("RETENTION_EXTREME", "x300um__retained_line_content", False),
        ("NEAR_TIP_RETENTION_EXTREME", "onset_max__near_tip_retained_line_content", False),
        ("PERSISTENT_WAKE_EXTREME", "x300um__wake_retained_departed_line_content", False),
        ("MOBILE_TRANSPORT_EXTREME", "x300um__mobile_centroid_tip_relative_m", False),
    )
    for role, field, absolute in role_fields:
        candidate = _pick_extreme(eligible, field, selected, absolute=absolute)
        selected.add(candidate); roles.append((role, candidate))
    backstress_score = (
        eligible["onset_max__backstress_Pa"] / max(float(eligible["onset_max__backstress_Pa"].median()), 1.0)
        + eligible["onset_max__K_shield_MPa_sqrt_m"].abs()
        / max(float(eligible["onset_max__K_shield_MPa_sqrt_m"].abs().median()), 1.0e-30)
    )
    eligible = eligible.assign(_backstress_shielding_score=backstress_score)
    candidate = _pick_extreme(eligible, "_backstress_shielding_score", selected)
    selected.add(candidate); roles.append(("BACKSTRESS_SHIELDING_EXTREME", candidate))
    selected_indices = [int(eligible.index[eligible.candidate_id.eq(candidate_id)][0]) for candidate_id in selected]
    median_distance = np.sqrt(np.sum(matrix * matrix, axis=1))
    central = np.flatnonzero(median_distance <= np.median(median_distance))
    min_to_selected = np.min(distances[:, selected_indices], axis=1)
    objective = min_to_selected - 0.25 * median_distance
    objective[[i for i, cid in enumerate(eligible.candidate_id) if cid in selected]] = -np.inf
    central_objective = objective[central]
    medoid_index = int(central[int(np.argmax(central_objective))])
    medoid = str(eligible.iloc[medoid_index].candidate_id)
    selected.add(medoid); roles.append(("BALANCED_MAXIMIN_MEDOID", medoid))
    selection = pd.DataFrame(roles, columns=["selection_role", "candidate_id"])
    selection["selection_order"] = np.arange(len(selection))
    selection = selection.merge(features, on="candidate_id", how="left", validate="one_to_one")
    if len(selection) != 8 or selection.candidate_id.nunique() != 8:
        raise RuntimeError("spatial transfer selection must contain eight unique candidates")
    ids = eligible.candidate_id.astype(str).tolist()
    pairwise = pd.DataFrame(
        {
            "candidate_id_a": np.repeat(ids, len(ids)),
            "candidate_id_b": np.tile(ids, len(ids)),
            "response_feature_distance": distances.reshape(-1),
            "feature_schema_json": json.dumps(columns, separators=(",", ":")),
        }
    )
    return selection.sort_values("selection_order"), pairwise


def _pf_registry(selection: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, Any]]:
    pool = _screen_pool().set_index("candidate_id", drop=False)
    base_registry = pd.read_csv(PF_TRANSFER_BASE)
    template = base_registry[base_registry.candidate_id.eq(CONTROL_ID)].iloc[0]
    records = []
    primary = []
    for _, selected in selection.sort_values("selection_order").iterrows():
        source = pool.loc[selected.candidate_id]
        row = template.copy()
        role_slug = str(selected.selection_role).lower()
        option_key = f"oneD_v2_DBTT_spatial_{int(selected.selection_order):02d}_{role_slug}"
        row["option_key"] = option_key
        row["candidate_id"] = str(source.candidate_id)
        row["role"] = str(selected.selection_role)
        row["mechanism_summary"] = "Taylor/Peierls-only spatial-transfer state extreme"
        row["validation_status"] = "BOUNDED_SPATIAL_TRANSFER_PF_300UM"
        for field in ACTIVE_CANDIDATE_PARAMETER_FIELDS:
            row[field] = float(source[field])
        records.append(row)
        primary.append({
            "option_key": option_key,
            "candidate_id": str(source.candidate_id),
            "material_class": "DBTT",
            "selection_role": str(selected.selection_role),
            "full_material_sha256": str(source.full_material_sha256),
            "taylor_peierls_subvector_sha256": str(source.taylor_peierls_subvector_sha256),
            "cleavage_barrier_sha256": str(source.cleavage_barrier_sha256),
            "emission_barrier_sha256": str(source.emission_barrier_sha256),
        })
    registry = pd.DataFrame(records)
    payload = {
        "schema": "oneD_v2_taylor_peierls_spatial_transfer_selection_v1",
        "analysis_only": True,
        "canonical_PF_registry_modified": False,
        "canonical_option_order": registry.option_key.astype(str).tolist(),
        "primary_candidates": primary,
        "temperature_K": TEMPERATURE_K,
        "hazard_seed": HAZARD_SEED,
        "target_extension_um": TARGET_UM,
        "theta_deg": 0.0,
        "bulk_plasticity_mode": "tip_only",
        "maximum_concurrent_heavy_workers": 2,
        "new_FEMCZM_runs": 0,
        "producer_commit": _head(),
    }
    return registry, payload


def finalize_selection() -> None:
    summary = pd.read_parquet(SUMMARY_PATH)
    history = pd.read_parquet(HISTORY_PATH)
    features = build_feature_vectors(history, summary)
    selection, pairwise = select_candidates(features)
    features.to_parquet(OUT / "oneD_v2_spatial_1d_response_feature_vectors.parquet", index=False)
    pairwise.to_parquet(OUT / "oneD_v2_spatial_1d_pairwise_distances.parquet", index=False)
    selection.to_csv(OUT / "oneD_v2_spatial_transfer_selected_candidates.csv", index=False)
    selected_ids = selection.candidate_id.astype(str).tolist()
    selected_pairwise = pairwise[
        pairwise.candidate_id_a.isin(selected_ids) & pairwise.candidate_id_b.isin(selected_ids)
    ]
    selected_pairwise.to_csv(OUT / "oneD_v2_spatial_transfer_selected_pairwise_distances.csv", index=False)
    registry, payload = _pf_registry(selection)
    registry_path = OUT / "oneD_v2_spatial_transfer_pf_registry.csv"
    registry.to_csv(registry_path, index=False)
    payload["installed_registry_sha256"] = _sha(registry_path)
    selection_path = OUT / "oneD_v2_spatial_transfer_pf_selection.json"
    selection_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

    _worker_init()
    profile_records = []
    selected_history = []
    for _, selected in selection.sort_values("selection_order").iterrows():
        source = _screen_pool().set_index("candidate_id", drop=False).loc[
            selected.candidate_id
        ]
        item = _run_candidate(source.to_dict(), save_profiles=True)
        profile_records.extend(item["profiles"])
        selected_history.extend(item["history"])
        print(f"SELECTED_PROFILE_COMPLETE role={selected.selection_role} candidate={selected.candidate_id}", flush=True)
    pd.DataFrame(selected_history).to_parquet(
        OUT / "oneD_v2_spatial_1d_selected_state_history.parquet", index=False
    )
    pd.DataFrame(profile_records).to_parquet(
        OUT / "oneD_v2_spatial_1d_selected_profiles.parquet", index=False
    )
    manifest = {
        "schema": "oneD_v2_taylor_peierls_spatial_1d_screen_v1",
        "producer_commit": _head(),
        "source_population_sha256": _sha(
            SOURCE_OUT / "oneD_v2_taylor_peierls_pareto_library.parquet"
        ),
        "candidate_count": int(len(summary)),
        "eligible_candidate_count": int(features.selection_eligible.sum()),
        "selected_candidate_count": int(len(selection)),
        "temperature_K": TEMPERATURE_K,
        "hazard_seed": HAZARD_SEED,
        "target_um": TARGET_UM,
        "state_model": "v9.13_persistent_site_backstress_blunting",
        "observer_feedback": False,
        "barriers_changed": False,
        "non_taylor_peierls_fields_changed": False,
        "new_FEMCZM_runs": 0,
        "artifacts": {},
    }
    for path in sorted(OUT.iterdir()):
        if path.is_file() and path.name != "oneD_v2_spatial_1d_manifest.json":
            manifest["artifacts"][path.name] = _sha(path)
    (OUT / "oneD_v2_spatial_1d_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    print(f"SPATIAL_SELECTION_COMPLETE eligible={int(features.selection_eligible.sum())} selected=8", flush=True)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("screen", "select", "all"))
    parser.add_argument("--workers", type=int, default=min(max(os.cpu_count() or 2, 2), 12))
    parser.add_argument("--checkpoint-every", type=int, default=8)
    args = parser.parse_args()
    if args.mode in {"screen", "all"}:
        run_screen(args.workers, args.checkpoint_every)
    if args.mode in {"select", "all"}:
        finalize_selection()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
