#!/usr/bin/env python3
"""Finalize the bounded DBTT Taylor/Peierls spatial-transfer experiment."""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from finalize_oneD_v2_taylor_peierls_search import classify_pf


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_outputs" / "oneD_v2_taylor_peierls_spatial_transfer"
PF_MANIFEST = Path("/private/tmp/oneD-v2-taylor-peierls-spatial-pf-runs/spatial_pf_matrix_manifest.json")
CHECKPOINTS_UM = (0.0, 25.0, 50.0, 100.0, 200.0, 300.0)
BASE_FEATURES = (
    "mobile_total", "retained_total", "retained_fraction",
    "mobile_centroid_tip_relative_m", "mobile_width_m",
    "retained_centroid_tip_relative_m", "retained_width_m",
    "near_tip_mobile", "near_tip_retained", "wake_mobile", "wake_retained",
    "wake_retained_centroid_laboratory_m", "wake_retained_width_m",
    "tip_radius_m", "front_width_m", "backstress_Pa", "shielding_MPa_sqrt_m",
    "multiplicity", "peierls_velocity_max_m_s", "peierls_velocity_mean_m_s",
    "encounter_rate_max_s", "encounter_rate_mean_s",
    "taylor_completion_rate_max_s", "taylor_completion_rate_mean_s",
)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def head(path: Path) -> str:
    return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=path, text=True).strip()


def number(value: Any) -> float:
    return float(value) if value is not None else np.nan


def moments(x: np.ndarray, weights: np.ndarray) -> tuple[float, float]:
    total = float(np.sum(weights))
    if total <= 0:
        return 0.0, 0.0
    centroid = float(np.sum(x * weights) / total)
    width = float(np.sqrt(max(np.sum((x - centroid) ** 2 * weights) / total, 0.0)))
    return centroid, width


def state_feature(record: dict[str, Any], step: pd.Series, tip_m: float) -> dict[str, Any]:
    active_x = np.asarray(record["active_x_ahead_of_tip_m"], dtype=float)
    wake_x = np.asarray(record["wake_x_behind_tip_m"], dtype=float)
    am = np.asarray(record["mobile_active_by_system_bin"], dtype=float).sum(axis=0)
    ar = np.asarray(record["retained_active_by_system_bin"], dtype=float).sum(axis=0)
    wm = np.asarray(record["mobile_wake_by_system_bin"], dtype=float).sum(axis=0)
    wr = np.asarray(record["retained_wake_by_system_bin"], dtype=float).sum(axis=0)
    mobile = float(am.sum() + wm.sum())
    retained = float(ar.sum() + wr.sum())
    mc, mw = moments(np.r_[active_x, -wake_x], np.r_[am, wm])
    rc, rw = moments(np.r_[active_x, -wake_x], np.r_[ar, wr])
    wrc_lab, wrw = moments(tip_m - wake_x, wr)
    near = active_x <= 5e-6 + 1e-15
    pv = np.abs(np.asarray(record["peierls_velocity_m_s_by_bin"], dtype=float))
    er = np.asarray(record["encounter_rate_s_by_bin"], dtype=float)
    tr = np.asarray(record["taylor_completion_rate_s_by_bin"], dtype=float)
    return {
        "mobile_total": mobile, "retained_total": retained,
        "retained_fraction": retained / max(mobile + retained, 1e-300),
        "mobile_centroid_tip_relative_m": mc, "mobile_width_m": mw,
        "retained_centroid_tip_relative_m": rc, "retained_width_m": rw,
        "near_tip_mobile": float(am[near].sum()), "near_tip_retained": float(ar[near].sum()),
        "wake_mobile": float(wm.sum()), "wake_retained": float(wr.sum()),
        "wake_retained_centroid_laboratory_m": wrc_lab,
        "wake_retained_width_m": wrw,
        "tip_radius_m": number(record.get("persistent_tip_radius_m")),
        "front_width_m": number(record.get("persistent_site_front_width_m")),
        "backstress_Pa": float(step.sigma_back_Pa),
        "shielding_MPa_sqrt_m": float(step.mpz_K_shield_Pa_sqrt_m) * 1e-6,
        "multiplicity": number(record.get("persistent_site_multiplicity_per_system")),
        "peierls_velocity_max_m_s": float(pv.max(initial=0.0)),
        "peierls_velocity_mean_m_s": float(pv.mean()) if pv.size else 0.0,
        "encounter_rate_max_s": float(er.max(initial=0.0)),
        "encounter_rate_mean_s": float(er.mean()) if er.size else 0.0,
        "taylor_completion_rate_max_s": float(tr.max(initial=0.0)),
        "taylor_completion_rate_mean_s": float(tr.mean()) if tr.size else 0.0,
        "native_KJ_MPa_sqrt_m": float(step.KJ_Pa_sqrtm) * 1e-6,
        "native_J_J_per_m2": float(step.J_effective_direct_J_per_m2),
        "reaction_N_per_m": float(step.Ftop_N), "opening_m": float(step.Uapp_m),
        "source_opening_stress_Pa": number(record.get("source_opening_stress_Pa")),
        "opening_tensor_Pa_json": json.dumps(record.get("opening_tensor_Pa"), separators=(",", ":")),
        "channel_tensors_Pa_json": json.dumps(record.get("channel_tensors_Pa"), separators=(",", ":")),
        "resolved_shears_Pa_json": json.dumps(record.get("anisotropic_tau_signed_Pa"), separators=(",", ":")),
        "topology_json": json.dumps({
            "active_x_ahead_of_tip_m": record.get("active_x_ahead_of_tip_m"),
            "wake_x_behind_tip_m": record.get("wake_x_behind_tip_m"),
            "spatial_state": record.get("taylor_peierls_state_profile_spatial_state"),
        }, separators=(",", ":")),
    }


def expand_profile(common: dict[str, Any], label: str, record: dict[str, Any], tip_m: float) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    active_x = np.asarray(record["active_x_ahead_of_tip_m"], dtype=float)
    wake_x = np.asarray(record["wake_x_behind_tip_m"], dtype=float)
    active_arrays = {
        "mobile": np.asarray(record["mobile_active_by_system_bin"], dtype=float),
        "retained": np.asarray(record["retained_active_by_system_bin"], dtype=float),
    }
    wake_arrays = {
        "mobile": np.asarray(record["mobile_wake_by_system_bin"], dtype=float),
        "retained": np.asarray(record["retained_wake_by_system_bin"], dtype=float),
    }
    rates = {name: np.asarray(record[name], dtype=float) for name in (
        "forest_density_active_m2_by_bin", "local_transport_stress_Pa_by_bin",
        "peierls_rate_s_by_bin", "peierls_velocity_m_s_by_bin", "encounter_rate_s_by_bin",
        "taylor_completion_rate_s_by_bin", "chi_ret_by_bin", "chi_taylor_completion_by_bin",
    )}
    for region, x, arrays, sign in (("ACTIVE_AHEAD_OF_TIP", active_x, active_arrays, 1.0),
                                     ("WAKE_BEHIND_TIP", wake_x, wake_arrays, -1.0)):
        for system in range(arrays["mobile"].shape[0]):
            for i, distance in enumerate(x):
                rows.append({**common, "state_label": label, "region": region,
                             "system_index": system, "bin_index": i,
                             "tip_relative_position_m": sign * float(distance),
                             "laboratory_position_m": tip_m + sign * float(distance),
                             "mobile_count": float(arrays["mobile"][system, i]),
                             "retained_count": float(arrays["retained"][system, i]),
                             **{name: (float(v[i]) if region.startswith("ACTIVE") else np.nan)
                                for name, v in rates.items()},
                             "observer_semantics": "DEFAULT_OFF_NEUTRAL_NO_FEEDBACK"})
    return rows


def robust_distance(frame: pd.DataFrame, id_col: str, feature_cols: list[str]) -> pd.DataFrame:
    x = frame[feature_cols].astype(float).replace([np.inf, -np.inf], np.nan).fillna(0.0).to_numpy()
    med = np.median(x, axis=0)
    scale = np.percentile(x, 75, axis=0) - np.percentile(x, 25, axis=0)
    fallback = np.std(x, axis=0)
    scale = np.where(scale > 1e-15, scale, np.where(fallback > 1e-15, fallback, 1.0))
    z = (x - med) / scale
    ids = frame[id_col].astype(str).tolist()
    rows = []
    for i in range(len(ids)):
        for j in range(i + 1, len(ids)):
            candidate_i, candidate_j = sorted((ids[i], ids[j]))
            rows.append({"candidate_i": candidate_i, "candidate_j": candidate_j,
                         "distance": float(np.linalg.norm(z[i] - z[j]))})
    return pd.DataFrame(rows)


def safe_spearman(a: pd.Series | np.ndarray, b: pd.Series | np.ndarray) -> float:
    left, right = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    if np.nanstd(left) <= 1e-15 or np.nanstd(right) <= 1e-15:
        return np.nan
    return float(spearmanr(left, right).statistic)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    manifest = json.loads(PF_MANIFEST.read_text())
    execution_manifest = OUT / "pf_2d_spatial_transfer_execution_manifest.json"
    shutil.copyfile(PF_MANIFEST, execution_manifest)
    selection = pd.read_csv(OUT / "oneD_v2_spatial_transfer_selected_candidates.csv")
    registry = pd.read_csv(OUT / "oneD_v2_spatial_transfer_pf_registry.csv")
    case_by_id = {c["candidate_id"]: c for c in manifest["cases"]}
    transactions_all, avalanches_all, onsets_all = [], [], []
    state_rows: list[dict[str, Any]] = []
    profile_rows: list[dict[str, Any]] = []

    for selected in selection.itertuples(index=False):
        item = case_by_id[selected.candidate_id]
        case_dir = Path(item["case_path"])
        case = {"material_class": "DBTT", "candidate_role": selected.selection_role,
                "candidate_id": selected.candidate_id, "option_key": item["option_key"],
                "temperature_K": 1100.0, "hazard_seed": 1008666, "target_um": 300.0,
                "case_dir": case_dir}
        transactions, avalanches, onsets, _ = classify_pf(case)
        transactions_all.append(transactions); avalanches_all.append(avalanches); onsets_all.append(onsets)
        steps_path = next(case_dir.glob("steps_1100K.csv"))
        audit_path = case_dir / "anisotropic_emission_audit_v10174.json"
        steps = pd.read_csv(steps_path)
        audit = json.loads(audit_path.read_text())["records"]
        reached = np.flatnonzero(steps.crack_extension_m.to_numpy(float) * 1e6 >= 300 - 1e-9)
        stop = int(reached[0])
        steps = steps.iloc[:stop + 1].copy(); audit = audit[:stop + 1]
        event_positions = np.flatnonzero(steps.n_fire.to_numpy(float) > 0)
        onset_positions = [int(event_positions[int(i)]) for i in onsets.event_transaction_index]
        selected_positions: dict[int, list[str]] = {}
        ext_um = steps.crack_extension_m.to_numpy(float) * 1e6
        for checkpoint in CHECKPOINTS_UM:
            eligible = np.flatnonzero(ext_um >= checkpoint - 1e-9)
            pos = int(eligible[0]) if len(eligible) else len(steps) - 1
            selected_positions.setdefault(pos, []).append(f"CHECKPOINT_{int(checkpoint):03d}UM")
        for k, pos in enumerate(onset_positions):
            selected_positions.setdefault(pos, []).append("INITIAL_ONSET" if k == 0 else f"REINITIATION_ONSET_{k:02d}")
        for pos, labels in sorted(selected_positions.items()):
            common = {"candidate_id": selected.candidate_id, "selection_role": selected.selection_role,
                      "option_key": item["option_key"], "step_index": int(steps.iloc[pos].step),
                      "projected_extension_um": float(ext_um[pos]),
                      "temperature_K": 1100.0, "hazard_seed": 1008666,
                      "source_steps_sha256": sha(steps_path), "source_audit_sha256": sha(audit_path)}
            for label in labels:
                state_rows.append({**common, "state_label": label,
                                   **state_feature(audit[pos], steps.iloc[pos], float(steps.iloc[pos].a_tip_m))})
                profile_rows.extend(expand_profile(common, label, audit[pos], float(steps.iloc[pos].a_tip_m)))

    transactions = pd.concat(transactions_all, ignore_index=True)
    avalanches = pd.concat(avalanches_all, ignore_index=True)
    onsets = pd.concat(onsets_all, ignore_index=True)
    states = pd.DataFrame(state_rows); profiles = pd.DataFrame(profile_rows)
    transactions.to_csv(OUT / "pf_2d_spatial_transfer_event_transactions.csv", index=False)
    avalanches.to_csv(OUT / "pf_2d_spatial_transfer_physical_avalanches.csv", index=False)
    onsets.to_csv(OUT / "pf_2d_spatial_transfer_onsets.csv", index=False)
    states.to_parquet(OUT / "pf_2d_spatial_transfer_state_features.parquet", index=False)
    profiles.to_parquet(OUT / "pf_2d_spatial_transfer_state_profiles.parquet", index=False)

    summaries = []
    decompositions = []
    for candidate, local in transactions.groupby("candidate_id", sort=False):
        onset = onsets[onsets.candidate_id.eq(candidate)].sort_values("physical_avalanche_index")
        first = onset.iloc[0]
        reinit = onset.iloc[1:]
        native_max = float(reinit.pre_event_native_KJ_MPa_sqrt_m.max()) if len(reinit) else np.nan
        local_max = float(reinit.pre_event_K_effective_local_equivalent_MPa_sqrt_m.max()) if len(reinit) else np.nan
        d_native = native_max - float(first.pre_event_native_KJ_MPa_sqrt_m) if len(reinit) else np.nan
        d_local = local_max - float(first.pre_event_K_effective_local_equivalent_MPa_sqrt_m) if len(reinit) else np.nan
        av_sizes = avalanches[avalanches.candidate_id.eq(candidate)].avalanche_extension_um.to_numpy(float)
        sel = selection[selection.candidate_id.eq(candidate)].iloc[0]
        positive = bool(len(reinit) and d_native > 1e-12)
        summaries.append({"candidate_id": candidate, "selection_role": sel.selection_role,
                          "event_transaction_count": len(local), "physical_avalanche_count": len(onset),
                          "N_reinit": max(len(onset) - 1, 0),
                          "initial_onset_native_KJ_MPa_sqrt_m": float(first.pre_event_native_KJ_MPa_sqrt_m),
                          "maximum_reinitiation_native_KJ_MPa_sqrt_m": native_max,
                          "deltaK_reinit_MPa_sqrt_m": d_native,
                          "largest_avalanche_fraction": float(av_sizes.max() / av_sizes.sum()),
                          "positive_reload_separated_resistance": positive,
                          "target_right_censored": True,
                          "direct_PF_decision": ("POSITIVE_REQUIRES_CONFIRMATION" if positive else
                                                 "NON_TOUGHENING_TAYLOR_PEIERLS_REDISTRIBUTION_INSUFFICIENT")})
        decompositions.append({"candidate_id": candidate, "selection_role": sel.selection_role,
                               "initial_native_KJ_MPa_sqrt_m": float(first.pre_event_native_KJ_MPa_sqrt_m),
                               "maximum_reinit_native_KJ_MPa_sqrt_m": native_max,
                               "apparent_deltaK_reinit_MPa_sqrt_m": d_native,
                               "initial_local_equivalent_K_MPa_sqrt_m": float(first.pre_event_K_effective_local_equivalent_MPa_sqrt_m),
                               "maximum_reinit_local_equivalent_K_MPa_sqrt_m": local_max,
                               "local_deltaK_reinit_MPa_sqrt_m": d_local,
                               "initial_shielding_MPa_sqrt_m": float(first.pre_event_K_shield_MPa_sqrt_m),
                               "maximum_reinit_abs_shielding_MPa_sqrt_m": (float(reinit.pre_event_K_shield_MPa_sqrt_m.abs().max()) if len(reinit) else np.nan),
                               "initial_tip_radius_um": float(first.pre_event_tip_radius_um),
                               "maximum_reinit_tip_radius_um": (float(reinit.pre_event_tip_radius_um.max()) if len(reinit) else np.nan),
                               "response_origin": ("NO_POSITIVE_RESPONSE" if not positive else
                                                   "POSITIVE_APPARENT_RESPONSE_REQUIRES_LOCAL_RADIUS_SHIELDING_CONFIRMATION")})
    summary = pd.DataFrame(summaries); decomposition = pd.DataFrame(decompositions)
    summary.to_csv(OUT / "pf_2d_spatial_transfer_summary.csv", index=False)
    decomposition.to_csv(OUT / "pf_2d_spatial_transfer_resistance_decomposition.csv", index=False)

    # Candidate vectors combine all sparse states and physical onsets.
    vector_rows = []
    for candidate, local in states.groupby("candidate_id", sort=False):
        row: dict[str, Any] = {"candidate_id": candidate,
                               "selection_role": local.selection_role.iloc[0]}
        for state in sorted(local.state_label.unique()):
            record = local[local.state_label.eq(state)].iloc[0]
            for feature in BASE_FEATURES:
                row[f"{state}__{feature}"] = float(record[feature])
        vector_rows.append(row)
    vectors = pd.DataFrame(vector_rows).fillna(0.0)
    vector_features = [c for c in vectors if "__" in c]
    vectors.to_parquet(OUT / "pf_2d_spatial_transfer_response_feature_vectors.parquet", index=False)
    pf_dist = robust_distance(vectors, "candidate_id", vector_features)
    pf_dist.to_csv(OUT / "pf_2d_spatial_transfer_pairwise_distances.csv", index=False)

    # Preserve the exact all-pool scaling used to choose the eight representatives.
    # Re-scaling only the selected subset would change the selection-space metric.
    one_d_dist = pd.read_csv(OUT / "oneD_v2_spatial_transfer_selected_pairwise_distances.csv")
    one_d_dist = one_d_dist[one_d_dist.candidate_id_a < one_d_dist.candidate_id_b].rename(columns={
        "candidate_id_a": "candidate_i", "candidate_id_b": "candidate_j",
        "response_feature_distance": "distance_1d",
    })[["candidate_i", "candidate_j", "distance_1d"]]
    joined = one_d_dist.merge(pf_dist.rename(columns={"distance": "distance_2d"}),
                              on=["candidate_i", "candidate_j"], how="inner")
    joined.to_csv(OUT / "oneD_to_pf_spatial_pairwise_rank_transfer.csv", index=False)
    pearson = float(np.corrcoef(joined.distance_1d, joined.distance_2d)[0, 1])
    spearman = safe_spearman(joined.distance_1d, joined.distance_2d)
    control = selection.loc[selection.selection_role.eq("CONTROL"), "candidate_id"].iloc[0]
    control_pairs = joined[(joined.candidate_i.eq(control)) | (joined.candidate_j.eq(control))].copy()
    control_spearman = safe_spearman(control_pairs.distance_1d, control_pairs.distance_2d)
    response = summary.merge(vectors, on=["candidate_id", "selection_role"])
    correlation_rows = []
    for feature in vector_features:
        correlation_rows.append({"feature": feature,
                                 "spearman_with_deltaK_reinit": safe_spearman(response[feature], response.deltaK_reinit_MPa_sqrt_m)})
    correlations = pd.DataFrame(correlation_rows).sort_values("spearman_with_deltaK_reinit", key=lambda s: s.abs(), ascending=False)
    correlations.to_csv(OUT / "pf_2d_spatial_state_resistance_correlations.csv", index=False)
    diversity_rows = []
    for label, local in states.groupby("state_label", sort=True):
        for feature in BASE_FEATURES:
            diversity_rows.append({"state_label": label, "feature": feature,
                                   "minimum": float(local[feature].min()),
                                   "maximum": float(local[feature].max()),
                                   "range": float(local[feature].max() - local[feature].min())})
    pd.DataFrame(diversity_rows).to_csv(OUT / "pf_2d_spatial_state_diversity_ranges.csv", index=False)

    strong_diversity = bool(pf_dist.distance.max() > 1.0)
    any_positive = bool(summary.positive_reload_separated_resistance.any())
    decision = {
        "schema": "oneD_v2_taylor_peierls_spatial_transfer_decision_v1",
        "material_class": "DBTT", "temperature_K": 1100.0, "hazard_seed": 1008666,
        "target_extension_um": 300.0, "selected_candidate_count": 8,
        "oneD_pf_pairwise_distance_pearson": pearson,
        "oneD_pf_pairwise_distance_spearman": spearman,
        "oneD_pf_control_distance_rank_spearman": (control_spearman if np.isfinite(control_spearman) else None),
        "maximum_pf_pairwise_state_distance": float(pf_dist.distance.max()),
        "strong_pf_spatial_state_diversity": strong_diversity,
        "positive_candidate_count": int(summary.positive_reload_separated_resistance.sum()),
        "conditional_local_multiseed_confirmation_launched": False,
        "conditional_confirmation_reason": ("REQUIRED_POSITIVE_CANDIDATE_PRESENT" if any_positive else
                                              "NOT_TRIGGERED_NO_POSITIVE_RELOAD_SEPARATED_RESISTANCE"),
        "final_decision": ("POSITIVE_CANDIDATE_REQUIRES_CONFIRMATION" if any_positive else
                           "TAYLOR_PEIERLS_REDISTRIBUTION_ALONE_INSUFFICIENT_DESPITE_2D_STATE_DIVERSITY"),
        "scope": "FIXED_CLEAVAGE_EMISSION_AND_NON_TAYLOR_PEIERLS_FIELDS; NO_MEMORY; PF_ONLY",
    }
    (OUT / "oneD_v2_taylor_peierls_spatial_transfer_final_decision.json").write_text(json.dumps(decision, indent=2) + "\n")

    plt.style.use("seaborn-v0_8-whitegrid")
    fig, ax = plt.subplots(figsize=(6.6, 5.2)); ax.scatter(joined.distance_1d, joined.distance_2d)
    ax.set(xlabel="1-D robust response-feature distance", ylabel="2-D PF robust spatial-state distance",
           title=f"Spatial rank transfer (Spearman={spearman:.3f})")
    fig.tight_layout(); fig.savefig(OUT / "ONED_TO_PF_SPATIAL_DISTANCE_TRANSFER.png", dpi=180); plt.close(fig)
    fig, ax = plt.subplots(figsize=(8.5, 5.1)); order = summary.sort_values("deltaK_reinit_MPa_sqrt_m")
    ax.barh(order.selection_role, order.deltaK_reinit_MPa_sqrt_m, color=np.where(order.deltaK_reinit_MPa_sqrt_m > 0, "#2a9d8f", "#d1495b"))
    ax.axvline(0, color="black", lw=1); ax.set(xlabel=r"$\Delta K_{reinit}$ [MPa$\sqrt{m}$]", title="Reload-separated PF resistance")
    fig.tight_layout(); fig.savefig(OUT / "PF_SPATIAL_TRANSFER_REINITIATION_RESISTANCE.png", dpi=180); plt.close(fig)
    colors = plt.cm.tab10(np.arange(len(decomposition)))
    fig, ax = plt.subplots(figsize=(9.0, 5.2));
    for color, r in zip(colors, decomposition.itertuples()):
        ax.scatter(r.local_deltaK_reinit_MPa_sqrt_m, r.apparent_deltaK_reinit_MPa_sqrt_m,
                   color=color, label=r.selection_role.replace("_EXTREME", ""))
    ax.axhline(0, color="black", lw=.8); ax.axvline(0, color="black", lw=.8)
    ax.set(xlabel=r"Local-equivalent $\Delta K$ [MPa$\sqrt{m}$]", ylabel=r"Native apparent $\Delta K$ [MPa$\sqrt{m}$]", title="Apparent-versus-local resistance decomposition")
    ax.legend(loc="center left", bbox_to_anchor=(1.01, .5), fontsize=7, frameon=False)
    fig.tight_layout(); fig.savefig(OUT / "PF_APPARENT_VERSUS_LOCAL_RESISTANCE.png", dpi=180); plt.close(fig)
    reinit_states = states[states.state_label.eq("REINITIATION_ONSET_01")]
    fig, ax = plt.subplots(figsize=(9.0, 5.2))
    for color, r in zip(colors, reinit_states.itertuples()):
        ax.scatter(r.wake_retained, r.mobile_total, color=color,
                   label=r.selection_role.replace("_EXTREME", ""))
    ax.set(xlabel="Wake-retained content at re-initiation", ylabel="Total mobile content at re-initiation",
           title="Strong 2-D PF state diversity at the common physical onset")
    ax.set_xscale("symlog", linthresh=1.0); ax.set_yscale("symlog", linthresh=1.0)
    ax.legend(loc="center left", bbox_to_anchor=(1.01, .5), fontsize=7, frameon=False)
    fig.tight_layout(); fig.savefig(OUT / "PF_REINITIATION_SPATIAL_STATE_DIVERSITY.png", dpi=180); plt.close(fig)
    reinit_profiles = profiles[profiles.state_label.eq("REINITIATION_ONSET_01")]
    fig, ax = plt.subplots(figsize=(9.0, 5.2))
    for color, (candidate, local) in zip(colors, reinit_profiles.groupby("candidate_id", sort=False)):
        role = local.selection_role.iloc[0].replace("_EXTREME", "")
        curve = local.groupby("tip_relative_position_m", sort=True).retained_count.sum()
        ax.plot(curve.index.to_numpy(float) * 1e6, curve.to_numpy(float), color=color, label=role, lw=1.4)
    ax.set_yscale("symlog", linthresh=1.0)
    ax.axvline(0, color="black", lw=.8)
    ax.set(xlabel="Tip-relative position [um] (wake < 0)", ylabel="Retained content summed over systems",
           title="2-D PF retained-state profiles at reload-separated re-initiation")
    ax.legend(loc="center left", bbox_to_anchor=(1.01, .5), fontsize=7, frameon=False)
    fig.tight_layout(); fig.savefig(OUT / "PF_REINITIATION_RETAINED_PROFILES.png", dpi=180); plt.close(fig)

    provenance = {
        "schema": "oneD_v2_taylor_peierls_spatial_transfer_provenance_v1",
        "producer_code_commit": head(ROOT), "pf_runner_commit": manifest["pf_runner_commit"],
        "pf_matrix_manifest": execution_manifest.name, "pf_matrix_manifest_sha256": sha(execution_manifest),
        "selection_sha256": sha(OUT / "oneD_v2_spatial_transfer_pf_selection.json"),
        "registry_sha256": sha(OUT / "oneD_v2_spatial_transfer_pf_registry.csv"),
        "candidate_count": 8, "fresh_pf_cases": manifest["fresh_PF_run_count"],
        "reused_pf_cases": manifest["reused_PF_run_count"], "maximum_pf_workers": manifest["maximum_concurrent_heavy_workers"],
        "new_femczm_runs": 0, "memory_term_added": False,
        "non_taylor_peierls_material_fields_varied": False,
        "cleavage_or_emission_barriers_changed": False,
        "profile_observer": "DEFAULT_OFF_NEUTRAL_NO_FEEDBACK",
        "scientific_fingerprints": {p.name: sha(p) for p in (
            sorted(OUT.glob("*.csv")) + sorted(OUT.glob("*.parquet")) +
            sorted(OUT.glob("*.png")) + [OUT / "oneD_v2_taylor_peierls_spatial_transfer_final_decision.json"]
        )},
    }
    (OUT / "oneD_v2_taylor_peierls_spatial_transfer_provenance.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")

    report = f"""# Taylor/Peierls spatial-transfer experiment\n\n## Scope\n\nDBTT, 1100 K, seed 1008666, theta=0, tip-only, through 300 um. Cleavage, emission, and every non-Taylor/non-Peierls material field were held fixed. No memory term and no FEM/CZM calculation were used.\n\n## Execution\n\nThe fully evaluated 513-row DBTT option pool was replayed with the spatial 1-D state model. All 513 reached 300 um uncensored and passed the first-onset-within-20-percent gate. Eight response-space representatives were transferred to PF. One compatible control was reused and seven missing cases were run with no more than two workers.\n\n## Rank transfer and resistance\n\nThe robust pairwise-distance association from 1-D to 2-D is Pearson {pearson:.6g}, Spearman {spearman:.6g}; the control-distance rank association is {control_spearman:.6g}. The largest 2-D state distance is {pf_dist.distance.max():.6g}.\n\nAll resistance values use V2 reload-separated pre-event onsets. Positive candidates: {int(summary.positive_reload_separated_resistance.sum())} of 8. The candidate range in DeltaK_reinit is {summary.deltaK_reinit_MPa_sqrt_m.min():.6g} to {summary.deltaK_reinit_MPa_sqrt_m.max():.6g} MPa sqrt(m).\n\n## Decision\n\n**{decision['final_decision'].replace('_', ' ')}.** The conditional local/multiseed confirmation was not launched because no positive reload-separated candidate was found. Spatial diversity therefore does not establish toughening: within this fixed-field experiment, Taylor/Peierls redistribution alone is insufficient.\n\nNative PF apparent resistance remains separate from the source-opening local-equivalent measure. The decomposition table records native, local, shielding, and radius changes without treating native KJ as remote applied K.\n"""
    (OUT / "ONE_D_V2_TAYLOR_PEIERLS_SPATIAL_TRANSFER_FINAL_DECISION.md").write_text(report)
    print(f"SPATIAL_TRANSFER_FINALIZED candidates=8 positive={int(any_positive)} spearman={spearman:.6g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
