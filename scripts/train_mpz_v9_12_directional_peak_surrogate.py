#!/usr/bin/env python3
"""Train direction-aware DBTT and genuine-peak v9.12 transfer surrogates."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import ExtraTreesClassifier, ExtraTreesRegressor
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    mean_absolute_error,
    r2_score,
    roc_auc_score,
)
from sklearn.model_selection import KFold, StratifiedKFold, cross_val_predict


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--table", required=True)
    p.add_argument("--out-model", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--trees", type=int, default=1200)
    p.add_argument("--folds", type=int, default=5)
    p.add_argument("--seed", type=int, default=9122)
    return p.parse_args()


def truth(series: pd.Series) -> np.ndarray:
    return series.astype(str).str.lower().isin(["1", "true", "yes"]).to_numpy()


def tree_cv_regressor(
    X: np.ndarray,
    y: np.ndarray,
    *,
    trees: int,
    folds: int,
    seed: int,
) -> tuple[ExtraTreesRegressor, dict[str, float]]:
    n_splits = max(2, min(int(folds), len(y)))
    model = ExtraTreesRegressor(
        n_estimators=trees,
        min_samples_leaf=2,
        max_features=0.7,
        n_jobs=-1,
        random_state=seed,
    )
    cv = KFold(n_splits=n_splits, shuffle=True, random_state=seed)
    pred = cross_val_predict(model, X, y, cv=cv, n_jobs=-1)
    metrics = {
        "n": int(len(y)),
        "folds": int(n_splits),
        "cv_r2": float(r2_score(y, pred)),
        "cv_mae": float(mean_absolute_error(y, pred)),
    }
    model.fit(X, y)
    return model, metrics


def tree_cv_classifier(
    X: np.ndarray,
    y: np.ndarray,
    *,
    trees: int,
    folds: int,
    seed: int,
) -> tuple[ExtraTreesClassifier | None, dict[str, object]]:
    classes, counts = np.unique(y, return_counts=True)
    if classes.size < 2:
        return None, {"constant": bool(classes[0]), "n": int(len(y))}
    n_splits = max(2, min(int(folds), int(counts.min())))
    model = ExtraTreesClassifier(
        n_estimators=trees,
        min_samples_leaf=2,
        max_features="sqrt",
        class_weight="balanced",
        n_jobs=-1,
        random_state=seed,
    )
    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    prob = cross_val_predict(
        model, X, y, cv=cv, method="predict_proba", n_jobs=-1
    )[:, 1]
    pred = prob >= 0.5
    metrics: dict[str, object] = {
        "n": int(len(y)),
        "positive": int(np.sum(y)),
        "folds": int(n_splits),
        "balanced_accuracy": float(balanced_accuracy_score(y, pred)),
        "roc_auc": float(roc_auc_score(y, prob)),
        "average_precision": float(average_precision_score(y, prob)),
    }
    model.fit(X, y)
    return model, metrics


def main() -> int:
    a = parse_args()
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(a.table, low_memory=False)
    features = [column for column in df.columns if column.startswith("x_")]
    if not features:
        raise RuntimeError("table contains no x_ features")

    Xdf = df[features].apply(pd.to_numeric, errors="coerce")
    medians = Xdf.median(axis=0).fillna(0.0)
    Xall = Xdf.fillna(medians).fillna(0.0).to_numpy(float)
    complete = (
        truth(df["y__complete_1d"])
        if "y__complete_1d" in df
        else np.ones(len(df), dtype=bool)
    )
    if int(complete.sum()) < 20:
        raise RuntimeError("too few completed 1-D rows")

    bundle: dict[str, object] = {
        "schema_version": 1,
        "model_kind": "mpz_v9_12_directional_peak_transfer",
        "feature_names": features,
        "feature_medians": medians.to_dict(),
        "training_candidate_ids": df.loc[complete, "candidate_id"].astype(str).tolist(),
        "models": {},
        "constants": {},
    }
    metrics: dict[str, object] = {
        "rows": int(len(df)),
        "complete_1d": int(complete.sum()),
        "features": int(len(features)),
    }
    importance_rows: list[dict[str, object]] = []

    regressions = [
        "y__log1p_directional_dbtt_gain_positive",
        "y__log1p_peak_prominence",
        "y__log1p_amplitude_1d",
        "y__min_source_available_fraction_pre_advance_developed_window_1d",
    ]
    for index, target in enumerate(regressions):
        if target not in df:
            continue
        yall = pd.to_numeric(df[target], errors="coerce").to_numpy(float)
        mask = complete & np.isfinite(yall)
        if int(mask.sum()) < 20:
            continue
        model, result = tree_cv_regressor(
            Xall[mask],
            yall[mask],
            trees=a.trees,
            folds=a.folds,
            seed=a.seed + index,
        )
        bundle["models"][target] = model
        metrics[target] = result
        for feature, value in zip(features, model.feature_importances_):
            importance_rows.append(
                {"model": target, "feature": feature, "importance": float(value)}
            )

    classifications = [
        "y__direction_correct_1d",
        "y__directional_dbtt_ge_threshold_1d",
        "y__peak_like_1d",
    ]
    for index, target in enumerate(classifications):
        if target not in df:
            continue
        valid = complete & df[target].notna().to_numpy()
        y = truth(df.loc[valid, target])
        if int(valid.sum()) < 20:
            continue
        model, result = tree_cv_classifier(
            Xall[valid],
            y,
            trees=a.trees,
            folds=a.folds,
            seed=a.seed + 100 + index,
        )
        metrics[target] = result
        if model is None:
            bundle["constants"][target] = bool(np.unique(y)[0])
        else:
            bundle["models"][target] = model
            for feature, value in zip(features, model.feature_importances_):
                importance_rows.append(
                    {"model": target, "feature": feature, "importance": float(value)}
                )

    bundle["training_metrics"] = metrics
    model_path = Path(a.out_model)
    model_path.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump(bundle, model_path)

    importance = pd.DataFrame(importance_rows)
    if not importance.empty:
        importance.to_csv(out / "feature_importance_by_model.csv", index=False)
        (
            importance.groupby("feature", as_index=False)["importance"]
            .mean()
            .sort_values("importance", ascending=False)
            .to_csv(out / "feature_importance_mean.csv", index=False)
        )

    (out / "cross_validation_metrics.json").write_text(
        json.dumps(metrics, indent=2, sort_keys=True) + "\n"
    )
    print(
        "DIRECTIONAL_PEAK_SURROGATE "
        f"rows={len(df)} complete={int(complete.sum())} "
        f"models={len(bundle['models'])} out={model_path}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
