#!/usr/bin/env python3
"""Select a diverse v9.12 batch for directional DBTT, peak, and exploration."""
from __future__ import annotations

import argparse
from pathlib import Path

import joblib
import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--model", required=True)
    p.add_argument("--pool-table", required=True)
    p.add_argument("--pool-registry", required=True)
    p.add_argument("--exclude-registry", action="append", default=[])
    p.add_argument("--batch-size", type=int, required=True)
    p.add_argument("--directional-fraction", type=float, default=0.52)
    p.add_argument("--peak-fraction", type=float, default=0.34)
    p.add_argument("--beta", type=float, default=1.5)
    p.add_argument("--preselect-factor", type=int, default=12)
    p.add_argument("--diversity-weight", type=float, default=0.30)
    p.add_argument("--out", required=True)
    return p.parse_args()


def robust_z(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    median = float(np.nanmedian(values))
    scale = float(np.nanmedian(np.abs(values - median))) * 1.4826
    if not np.isfinite(scale) or scale <= 1.0e-12:
        scale = float(np.nanstd(values))
    return (values - median) / max(scale, 1.0e-12)


def tree_mean_std(model, X: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    prediction = np.asarray(
        [tree.predict(X) for tree in model.estimators_], dtype=float
    )
    return np.mean(prediction, axis=0), np.std(prediction, axis=0)


def probability_true(model, X: np.ndarray) -> np.ndarray:
    classes = list(model.classes_)
    index = classes.index(True) if True in classes else classes.index(1)
    return model.predict_proba(X)[:, index]


def probability(bundle: dict[str, object], target: str, X: np.ndarray) -> np.ndarray:
    models = bundle["models"]
    constants = bundle.get("constants", {})
    if target in models:
        return probability_true(models[target], X)
    if target in constants:
        return np.full(len(X), float(bool(constants[target])))
    return np.ones(len(X), dtype=float)


def diverse_select(
    candidates: np.ndarray,
    acquisition: np.ndarray,
    Xscaled: np.ndarray,
    count: int,
    *,
    preselect_factor: int,
    diversity_weight: float,
) -> list[int]:
    if count <= 0 or candidates.size == 0:
        return []
    count = min(int(count), int(candidates.size))
    order = candidates[np.argsort(acquisition[candidates])[::-1]]
    pre_n = min(len(order), max(count, int(preselect_factor) * count))
    pre = order[:pre_n]
    chosen = [int(pre[0])]
    remaining = set(map(int, pre[1:]))
    amin = float(np.min(acquisition[pre]))
    amax = float(np.max(acquisition[pre]))

    while len(chosen) < count and remaining:
        options = np.fromiter(remaining, dtype=int)
        min_distance = np.min(
            np.linalg.norm(
                Xscaled[options, None, :]
                - Xscaled[np.asarray(chosen)][None, :, :],
                axis=2,
            ),
            axis=1,
        )
        distance_scale = max(float(np.max(min_distance)), 1.0e-12)
        acquisition_norm = (
            acquisition[options] - amin
        ) / max(amax - amin, 1.0e-12)
        merit = (
            (1.0 - diversity_weight) * acquisition_norm
            + diversity_weight * (min_distance / distance_scale)
        )
        pick = int(options[int(np.argmax(merit))])
        chosen.append(pick)
        remaining.remove(pick)
    return chosen


def main() -> int:
    a = parse_args()
    bundle = joblib.load(a.model)
    table = pd.read_csv(a.pool_table, low_memory=False)
    registry = pd.read_csv(a.pool_registry, low_memory=False)
    if not table["candidate_id"].is_unique or not registry["candidate_id"].is_unique:
        raise RuntimeError("candidate_id must be unique")

    excluded: set[str] = set(bundle.get("training_candidate_ids", []))
    for path in a.exclude_registry:
        other = pd.read_csv(path, usecols=["candidate_id"])
        excluded.update(other["candidate_id"].astype(str))

    table["candidate_id"] = table["candidate_id"].astype(str)
    registry["candidate_id"] = registry["candidate_id"].astype(str)
    keep = ~table["candidate_id"].isin(excluded)
    table = table.loc[keep].reset_index(drop=True)
    if table.empty:
        raise RuntimeError("no candidates remain after exclusions")

    features = list(bundle["feature_names"])
    medians = pd.Series(bundle["feature_medians"])
    Xdf = table.reindex(columns=features).apply(pd.to_numeric, errors="coerce")
    X = Xdf.fillna(medians).fillna(0.0).to_numpy(float)
    center = np.median(X, axis=0)
    scale = np.median(np.abs(X - center), axis=0) * 1.4826
    scale[~np.isfinite(scale) | (scale <= 1.0e-12)] = 1.0
    Xscaled = (X - center) / scale

    models = bundle["models"]
    direction_model = models["y__log1p_directional_dbtt_gain_positive"]
    peak_model = models["y__log1p_peak_prominence"]
    direction_mean, direction_std = tree_mean_std(direction_model, X)
    peak_mean, peak_std = tree_mean_std(peak_model, X)
    p_direction = probability(bundle, "y__direction_correct_1d", X)
    p_direction5 = probability(
        bundle, "y__directional_dbtt_ge_threshold_1d", X
    )
    p_peak = probability(bundle, "y__peak_like_1d", X)

    direction_acq = (
        (0.20 + 0.40 * p_direction + 0.40 * p_direction5)
        * (
            robust_z(direction_mean)
            + float(a.beta) * np.maximum(robust_z(direction_std), 0.0)
            + 3.0
        )
    )
    peak_acq = (
        (0.20 + 0.80 * p_peak)
        * (
            robust_z(peak_mean)
            + float(a.beta) * np.maximum(robust_z(peak_std), 0.0)
            + 3.0
        )
    )
    exploration_acq = np.maximum(
        robust_z(direction_std), robust_z(peak_std)
    ) + 2.0

    batch = min(max(int(a.batch_size), 1), len(table))
    n_direction = int(round(batch * float(a.directional_fraction)))
    n_peak = int(round(batch * float(a.peak_fraction)))
    n_direction = min(n_direction, batch)
    n_peak = min(n_peak, batch - n_direction)
    n_explore = batch - n_direction - n_peak

    all_indices = np.arange(len(table), dtype=int)
    family = table.get(
        "campaign_parent_family", pd.Series([""] * len(table))
    ).fillna("").astype(str).str.lower()
    peak_family = family.eq("peak").to_numpy()

    chosen_direction = diverse_select(
        all_indices,
        direction_acq,
        Xscaled,
        n_direction,
        preselect_factor=a.preselect_factor,
        diversity_weight=a.diversity_weight,
    )
    used = set(chosen_direction)

    peak_candidates = all_indices[~np.isin(all_indices, list(used))]
    if np.any(peak_family):
        restricted = peak_candidates[peak_family[peak_candidates]]
        if restricted.size >= n_peak:
            peak_candidates = restricted
    chosen_peak = diverse_select(
        peak_candidates,
        peak_acq,
        Xscaled,
        n_peak,
        preselect_factor=a.preselect_factor,
        diversity_weight=a.diversity_weight,
    )
    used.update(chosen_peak)

    explore_candidates = all_indices[~np.isin(all_indices, list(used))]
    chosen_explore = diverse_select(
        explore_candidates,
        exploration_acq,
        Xscaled,
        n_explore,
        preselect_factor=a.preselect_factor,
        diversity_weight=max(a.diversity_weight, 0.50),
    )

    labels = (
        [(index, "directional") for index in chosen_direction]
        + [(index, "peak") for index in chosen_peak]
        + [(index, "exploration") for index in chosen_explore]
    )
    selected_indices = [index for index, _ in labels]
    selected = table.iloc[selected_indices][["candidate_id"]].copy()
    selected["acquisition_role"] = [label for _, label in labels]
    selected["pred_direction_log1p_mean"] = direction_mean[selected_indices]
    selected["pred_direction_log1p_std"] = direction_std[selected_indices]
    selected["pred_peak_log1p_mean"] = peak_mean[selected_indices]
    selected["pred_peak_log1p_std"] = peak_std[selected_indices]
    selected["pred_p_direction_correct"] = p_direction[selected_indices]
    selected["pred_p_direction_ge_threshold"] = p_direction5[selected_indices]
    selected["pred_p_peak_like"] = p_peak[selected_indices]
    selected["active_rank"] = np.arange(1, len(selected) + 1)

    output = selected.merge(
        registry, on="candidate_id", how="left", validate="one_to_one"
    )
    if output["candidate_id"].isna().any() or len(output) != batch:
        raise RuntimeError("failed to join selected candidates to registry")
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(out, index=False)
    print(
        "DIRECTIONAL_PEAK_BATCH "
        f"selected={len(output)} directional={len(chosen_direction)} "
        f"peak={len(chosen_peak)} exploration={len(chosen_explore)} "
        f"pool={len(table)} out={out}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
