#!/usr/bin/env python3
"""Select a diverse active-learning batch from a v9.12 candidate pool."""
from __future__ import annotations

import argparse
from pathlib import Path

import joblib
import numpy as np
import pandas as pd


def args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--model", required=True)
    p.add_argument("--pool-ml-table", required=True)
    p.add_argument("--pool-registry", required=True)
    p.add_argument("--batch-size", type=int, default=256)
    p.add_argument("--beta", type=float, default=1.5)
    p.add_argument("--preselect-factor", type=int, default=20)
    p.add_argument("--diversity-weight", type=float, default=0.35)
    p.add_argument("--out", required=True)
    return p.parse_args()


def probability_true(model, X: np.ndarray) -> np.ndarray:
    classes = list(model.classes_)
    target = classes.index(True) if True in classes else classes.index(1)
    return model.predict_proba(X)[:, target]


def tree_mean_std(model, X: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    pred = np.asarray([tree.predict(X) for tree in model.estimators_], dtype=float)
    return np.mean(pred, axis=0), np.std(pred, axis=0)


def robust_z(values: np.ndarray) -> np.ndarray:
    med = float(np.median(values))
    scale = float(np.median(np.abs(values - med))) * 1.4826
    if not np.isfinite(scale) or scale <= 1.0e-12:
        scale = float(np.std(values))
    return (values - med) / max(scale, 1.0e-12)


def main() -> int:
    a = args()
    bundle = joblib.load(a.model)
    table = pd.read_csv(a.pool_ml_table)
    registry = pd.read_csv(a.pool_registry)
    if not table["candidate_id"].is_unique or not registry["candidate_id"].is_unique:
        raise RuntimeError("candidate_id must be unique in pool table and registry")

    features = list(bundle["feature_names"])
    medians = pd.Series(bundle["feature_medians"])
    Xdf = table.reindex(columns=features).apply(pd.to_numeric, errors="coerce")
    Xdf = Xdf.fillna(medians).fillna(0.0)
    X = Xdf.to_numpy(float)
    models = bundle["models"]

    p_complete = (
        probability_true(models["complete_classifier"], X)
        if "complete_classifier" in models else np.ones(len(X))
    )
    p_pass = (
        probability_true(models["pass_classifier"], X)
        if "pass_classifier" in models else np.ones(len(X))
    )
    score_mean, score_std = tree_mean_std(models["y__score"], X)
    acquisition = (
        p_complete
        * (0.10 + 0.90 * p_pass)
        * (robust_z(score_mean) + a.beta * np.maximum(robust_z(score_std), 0.0) + 2.0)
    )

    batch = min(max(a.batch_size, 1), len(table))
    pre_n = min(len(table), max(batch, a.preselect_factor * batch))
    pre = np.argsort(acquisition)[::-1][:pre_n]

    Xs = X.copy()
    center = np.median(Xs, axis=0)
    scale = np.median(np.abs(Xs - center), axis=0) * 1.4826
    scale[~np.isfinite(scale) | (scale <= 1.0e-12)] = 1.0
    Xs = (Xs - center) / scale

    chosen = [int(pre[0])]
    remaining = set(map(int, pre[1:]))
    acq_pre = acquisition[pre]
    amin, amax = float(np.min(acq_pre)), float(np.max(acq_pre))
    while len(chosen) < batch and remaining:
        candidates = np.fromiter(remaining, dtype=int)
        min_dist = np.min(
            np.linalg.norm(
                Xs[candidates, None, :] - Xs[np.asarray(chosen)][None, :, :], axis=2
            ),
            axis=1,
        )
        dscale = max(float(np.max(min_dist)), 1.0e-12)
        anorm = (acquisition[candidates] - amin) / max(amax - amin, 1.0e-12)
        merit = (
            (1.0 - a.diversity_weight) * anorm
            + a.diversity_weight * (min_dist / dscale)
        )
        pick = int(candidates[int(np.argmax(merit))])
        chosen.append(pick)
        remaining.remove(pick)

    selected = table.iloc[chosen][["candidate_id"]].copy()
    selected["active_acquisition"] = acquisition[chosen]
    selected["active_p_complete"] = p_complete[chosen]
    selected["active_p_pass"] = p_pass[chosen]
    selected["active_score_mean"] = score_mean[chosen]
    selected["active_score_std"] = score_std[chosen]
    selected["active_rank"] = np.arange(1, len(selected) + 1)
    output = selected.merge(registry, on="candidate_id", how="left", validate="one_to_one")
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(out, index=False)
    print(f"ACTIVE_BATCH selected={len(output)} pool={len(table)} out={out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
