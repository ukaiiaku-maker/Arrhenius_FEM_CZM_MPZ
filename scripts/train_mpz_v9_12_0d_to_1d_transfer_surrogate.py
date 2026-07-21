#!/usr/bin/env python3
"""Train the MPZ v9.12 multi-fidelity 0-D -> 1-D transfer surrogate.

The model learns spatial attenuation/amplification, 1-D amplitude, transition
location, persistence, shielding, and the probability of retaining at least
50 MPa sqrt(m).  Large 0-D responses remain in the training data and are handled
with logarithmic targets rather than removed by an upper-amplitude cutoff.
"""
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
    p.add_argument("--transfer-table", required=True)
    p.add_argument("--out-model", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--trees", type=int, default=1000)
    p.add_argument("--folds", type=int, default=5)
    p.add_argument("--seed", type=int, default=9121)
    return p.parse_args()


def truth(series: pd.Series) -> np.ndarray:
    return series.astype(str).str.lower().isin(["1", "true", "yes"]).to_numpy()


def safe_folds(requested: int, n: int) -> int:
    return max(2, min(int(requested), int(n)))


def main() -> int:
    a = parse_args()
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(a.transfer_table, low_memory=False)
    features = [column for column in df.columns if column.startswith("x_")]
    if not features:
        raise RuntimeError("transfer table contains no x_ features")

    Xdf = df[features].apply(pd.to_numeric, errors="coerce")
    medians = Xdf.median(axis=0).fillna(0.0)
    X = Xdf.fillna(medians).to_numpy(float)
    complete = truth(df["y__complete_1d"]) if "y__complete_1d" in df else np.zeros(len(df), dtype=bool)
    if not np.any(complete):
        raise RuntimeError("no completed 1-D rows in transfer table")

    bundle: dict[str, object] = {
        "schema_version": 1,
        "model_kind": "mpz_v9_12_0d_to_1d_transfer",
        "feature_names": features,
        "feature_medians": medians.to_dict(),
        "models": {},
    }
    metrics: dict[str, object] = {
        "rows": int(len(df)),
        "features": int(len(features)),
        "complete_1d": int(complete.sum()),
    }
    importance = np.zeros(len(features), dtype=float)
    importance_models = 0

    # Numerical feasibility classifier.  It is fitted only when both classes
    # exist; otherwise the observed constant class is recorded in the bundle.
    unique_complete = np.unique(complete)
    if unique_complete.size == 2:
        counts = np.bincount(complete.astype(int), minlength=2)
        n_splits = max(2, min(safe_folds(a.folds, len(df)), int(counts.min())))
        classifier = ExtraTreesClassifier(
            n_estimators=a.trees,
            min_samples_leaf=2,
            max_features="sqrt",
            class_weight="balanced",
            n_jobs=-1,
            random_state=a.seed,
        )
        cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=a.seed)
        prob = cross_val_predict(classifier, X, complete, cv=cv, method="predict_proba", n_jobs=-1)[:, 1]
        pred = prob >= 0.5
        metrics["complete_classifier_cv"] = {
            "folds": n_splits,
            "balanced_accuracy": float(balanced_accuracy_score(complete, pred)),
            "roc_auc": float(roc_auc_score(complete, prob)),
            "average_precision": float(average_precision_score(complete, prob)),
        }
        classifier.fit(X, complete)
        bundle["models"]["complete_classifier"] = classifier
        importance += classifier.feature_importances_
        importance_models += 1
    else:
        bundle["complete_constant"] = bool(unique_complete[0])
        metrics["complete_constant"] = bool(unique_complete[0])

    # Primary scientific classifier: whether the 1-D response retains at least
    # 50 MPa sqrt(m).  Train only on completed 1-D cases.
    threshold_target = "y__retains_50_MPa_sqrt_m_1d"
    if threshold_target in df:
        yclass = truth(df[threshold_target])
        mask = complete & df[threshold_target].notna().to_numpy()
        classes = np.unique(yclass[mask])
        if mask.sum() >= 10 and classes.size == 2:
            counts = np.bincount(yclass[mask].astype(int), minlength=2)
            n_splits = max(2, min(safe_folds(a.folds, int(mask.sum())), int(counts.min())))
            classifier = ExtraTreesClassifier(
                n_estimators=a.trees,
                min_samples_leaf=2,
                max_features="sqrt",
                class_weight="balanced",
                n_jobs=-1,
                random_state=a.seed + 1,
            )
            cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=a.seed + 1)
            prob = cross_val_predict(
                classifier,
                X[mask],
                yclass[mask],
                cv=cv,
                method="predict_proba",
                n_jobs=-1,
            )[:, 1]
            pred = prob >= 0.5
            metrics["retains_50_classifier_cv"] = {
                "n": int(mask.sum()),
                "positive": int(yclass[mask].sum()),
                "folds": n_splits,
                "balanced_accuracy": float(balanced_accuracy_score(yclass[mask], pred)),
                "roc_auc": float(roc_auc_score(yclass[mask], prob)),
                "average_precision": float(average_precision_score(yclass[mask], prob)),
            }
            classifier.fit(X[mask], yclass[mask])
            bundle["models"]["retains_50_classifier"] = classifier
            importance += classifier.feature_importances_
            importance_models += 1
        elif classes.size == 1 and mask.sum() > 0:
            bundle["retains_50_constant"] = bool(classes[0])
            metrics["retains_50_constant"] = bool(classes[0])

    regression_targets = [
        "y__log1p_amplitude_1d",
        "y__log10_transfer_ratio_1d_over_0d",
        "y__jump_temperature_K_1d",
        "y__localization_1d",
        "y__persistence_1d",
        "y__post_peak_drop_fraction_1d",
        "y__log1p_max_shield_1d",
        "y__max_tau_gnd_tip_MPa_1d",
        "y__min_source_available_fraction_1d",
    ]
    for index, target in enumerate(regression_targets):
        if target not in df:
            continue
        y = pd.to_numeric(df[target], errors="coerce").to_numpy(float)
        mask = complete & np.isfinite(y)
        n = int(mask.sum())
        if n < 10:
            continue
        n_splits = safe_folds(a.folds, n)
        regressor = ExtraTreesRegressor(
            n_estimators=a.trees,
            min_samples_leaf=2,
            max_features=0.7,
            n_jobs=-1,
            random_state=a.seed + 10 + index,
        )
        cv = KFold(n_splits=n_splits, shuffle=True, random_state=a.seed + 10 + index)
        pred = cross_val_predict(regressor, X[mask], y[mask], cv=cv, n_jobs=-1)
        metrics[target] = {
            "n": n,
            "folds": n_splits,
            "cv_r2": float(r2_score(y[mask], pred)),
            "cv_mae": float(mean_absolute_error(y[mask], pred)),
        }
        regressor.fit(X[mask], y[mask])
        bundle["models"][target] = regressor
        importance += regressor.feature_importances_
        importance_models += 1

    bundle["training_metrics"] = metrics
    model_path = Path(a.out_model)
    model_path.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump(bundle, model_path)

    if importance_models:
        mean_importance = importance / float(importance_models)
        order = np.argsort(mean_importance)[::-1]
        pd.DataFrame(
            {
                "feature": np.asarray(features)[order],
                "importance": mean_importance[order],
            }
        ).to_csv(out / "feature_importance.csv", index=False)

    (out / "cross_validation_metrics.json").write_text(
        json.dumps(metrics, indent=2, sort_keys=True) + "\n"
    )
    print(
        "TRANSFER_SURROGATE_TRAINED "
        f"rows={len(df)} complete_1d={int(complete.sum())} "
        f"features={len(features)} models={len(bundle['models'])} out={model_path}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
