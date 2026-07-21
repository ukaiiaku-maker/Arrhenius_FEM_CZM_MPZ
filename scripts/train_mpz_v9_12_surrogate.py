#!/usr/bin/env python3
"""Train an Extra-Trees active-learning surrogate for v9.12."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import ExtraTreesClassifier, ExtraTreesRegressor
from sklearn.metrics import accuracy_score, mean_absolute_error, r2_score


def args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--ml-table", required=True)
    p.add_argument("--out-model", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--trees", type=int, default=500)
    p.add_argument("--seed", type=int, default=912)
    return p.parse_args()


def truth(series: pd.Series) -> np.ndarray:
    return series.astype(str).str.lower().isin(["1", "true", "yes"]).to_numpy()


def main() -> int:
    a = args()
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(a.ml_table)
    features = [c for c in df.columns if c.startswith("x_")]
    Xdf = df[features].apply(pd.to_numeric, errors="coerce")
    medians = Xdf.median(axis=0).fillna(0.0)
    X = Xdf.fillna(medians).to_numpy(float)

    if "y__status" in df:
        complete = df["y__status"].fillna("").astype(str).eq("complete").to_numpy()
    else:
        complete = np.isfinite(pd.to_numeric(df.get("y__score"), errors="coerce"))
    if not np.any(complete):
        raise RuntimeError("no completed campaign rows in ML table")

    bundle: dict[str, object] = {
        "feature_names": features,
        "feature_medians": medians.to_dict(),
        "models": {},
        "schema_version": 1,
    }
    metrics: dict[str, object] = {
        "rows": len(df),
        "features": len(features),
        "complete": int(complete.sum()),
    }
    importances = np.zeros(len(features), dtype=float)
    importance_models = 0

    if np.unique(complete).size == 2:
        clf = ExtraTreesClassifier(
            n_estimators=a.trees,
            min_samples_leaf=2,
            max_features="sqrt",
            class_weight="balanced",
            n_jobs=-1,
            random_state=a.seed,
        ).fit(X, complete)
        bundle["models"]["complete_classifier"] = clf
        importances += clf.feature_importances_
        importance_models += 1
        metrics["complete_training_accuracy"] = float(
            accuracy_score(complete, clf.predict(X))
        )

    if "y__pass" in df:
        passed = truth(df["y__pass"])
        mask = complete & df["y__pass"].notna().to_numpy()
        if np.unique(passed[mask]).size == 2:
            pclf = ExtraTreesClassifier(
                n_estimators=a.trees,
                min_samples_leaf=2,
                max_features="sqrt",
                class_weight="balanced",
                n_jobs=-1,
                random_state=a.seed + 1,
            ).fit(X[mask], passed[mask])
            bundle["models"]["pass_classifier"] = pclf
            importances += pclf.feature_importances_
            importance_models += 1
            metrics["pass_training_accuracy"] = float(
                accuracy_score(passed[mask], pclf.predict(X[mask]))
            )

    targets = [
        "y__score",
        "y__amplitude_MPa_sqrt_m",
        "y__largest_jump_localization",
        "y__transition_width_10_90_K",
        "y__max_abs_K_shield_MPa_sqrt_m",
        "y__max_gnd_abs_line_count_per_unit_thickness",
    ]
    for index, target in enumerate(targets):
        if target not in df:
            continue
        y = pd.to_numeric(df[target], errors="coerce").to_numpy(float)
        mask = complete & np.isfinite(y)
        if mask.sum() < 20:
            continue
        model = ExtraTreesRegressor(
            n_estimators=a.trees,
            min_samples_leaf=2,
            max_features=0.7,
            n_jobs=-1,
            random_state=a.seed + 10 + index,
        ).fit(X[mask], y[mask])
        bundle["models"][target] = model
        pred = model.predict(X[mask])
        metrics[target] = {
            "n": int(mask.sum()),
            "training_r2": float(r2_score(y[mask], pred)),
            "training_mae": float(mean_absolute_error(y[mask], pred)),
        }
        importances += model.feature_importances_
        importance_models += 1

    if "y__score" not in bundle["models"]:
        raise RuntimeError("at least 20 finite completed score rows are required")

    bundle["training_metrics"] = metrics
    model_path = Path(a.out_model)
    model_path.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump(bundle, model_path)

    if importance_models:
        importance = importances / importance_models
        order = np.argsort(importance)[::-1]
        pd.DataFrame(
            {"feature": np.asarray(features)[order], "importance": importance[order]}
        ).to_csv(out / "feature_importance.csv", index=False)
    (out / "training_metrics.json").write_text(
        json.dumps(metrics, indent=2, sort_keys=True) + "\n"
    )
    print(
        f"SURROGATE_TRAINED rows={len(df)} features={len(features)} "
        f"models={len(bundle['models'])} out={model_path}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
