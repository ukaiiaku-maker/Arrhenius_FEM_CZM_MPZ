#!/usr/bin/env python3
"""Quality-diversity promotion for the existing v9.10.3 isotropic MPZ search.

This script does not modify constitutive physics or refit candidates. It consumes
v9.10.3 global-search outputs and replaces objective-only shortlist truncation
with an auditable quality-diversity selection for spatial promotion.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

MODEL_ID = "v9.10.4_isotropic_mpz_quality_diversity_promotion"
RESPONSE_METRICS = (
    "K_init_proxy",
    "K_plateau_proxy",
    "delta_KR_proxy",
    "early_rise_per_100um_proxy",
    "plateau_rise_per_100um_proxy",
)
EXCLUDED_PARAMETER_FIELDS = {
    "candidate_id",
    "target_class",
    "restart",
    "candidate_source",
    "objective",
    "accepted_for_spatial_promotion",
    "acceptance_reason",
    "de_success",
    "local_success",
}


@dataclass(frozen=True)
class SelectionConfig:
    count: int = 10
    quality_reserve_fraction: float = 0.30
    quality_weight: float = 0.35
    parameter_weight: float = 0.45
    response_weight: float = 0.55
    pool_factor: int = 12
    preserve_restart_lineages: bool = True

    def validate(self) -> "SelectionConfig":
        if self.count < 1:
            raise ValueError("count must be at least one")
        if not 0.0 <= self.quality_reserve_fraction <= 1.0:
            raise ValueError("quality_reserve_fraction must lie in [0,1]")
        if self.quality_weight < 0.0:
            raise ValueError("quality_weight must be non-negative")
        if self.parameter_weight < 0.0 or self.response_weight < 0.0:
            raise ValueError("diversity weights must be non-negative")
        if self.parameter_weight + self.response_weight <= 0.0:
            raise ValueError("at least one diversity weight must be positive")
        if self.pool_factor < 1:
            raise ValueError("pool_factor must be at least one")
        return self


def _bool_series(values: pd.Series) -> pd.Series:
    if values.dtype == bool:
        return values.fillna(False)
    normalized = values.astype(str).str.strip().str.lower()
    return normalized.isin({"1", "true", "yes", "y", "pass", "passed"})


def _finite_numeric(frame: pd.DataFrame, columns: Iterable[str]) -> list[str]:
    kept: list[str] = []
    for name in columns:
        if name not in frame:
            continue
        values = pd.to_numeric(frame[name], errors="coerce")
        if values.notna().sum() >= 2 and float(values.max() - values.min()) > 0.0:
            kept.append(name)
    return kept


def infer_parameter_columns(candidates: pd.DataFrame) -> list[str]:
    preferred = [
        name
        for name in candidates.columns
        if name not in EXCLUDED_PARAMETER_FIELDS
        and not name.endswith("_loss")
        and not name.startswith("K_")
        and name
        not in {
            "status",
            "search_initialization",
            "objective_mode",
            "barrier_order_margin_eV",
            "min_raw_barrier_eV",
            "min_peierls_traverse_number",
            "plateau_temperature_rise",
            "max_K_shield_MPa_sqrt_m",
            "parameter_count",
            "full_search_space",
            "prior_shortlist_used",
        }
    ]
    return _finite_numeric(candidates, preferred)


def build_response_table(
    candidates: pd.DataFrame, temperature_detail: pd.DataFrame
) -> tuple[pd.DataFrame, list[str]]:
    base = pd.DataFrame({"candidate_id": candidates["candidate_id"].astype(str)})
    if temperature_detail.empty or "candidate_id" not in temperature_detail:
        return base, []
    detail = temperature_detail.copy()
    detail["candidate_id"] = detail["candidate_id"].astype(str)
    detail["T_K"] = pd.to_numeric(detail.get("T_K"), errors="coerce")
    parts = [base.set_index("candidate_id")]
    feature_names: list[str] = []
    for metric in RESPONSE_METRICS:
        if metric not in detail:
            continue
        values = pd.to_numeric(detail[metric], errors="coerce")
        table = detail.assign(_value=values).pivot_table(
            index="candidate_id", columns="T_K", values="_value", aggfunc="mean"
        )
        renamed = {
            temperature: f"response_{metric}_{int(round(float(temperature)))}K"
            for temperature in table.columns
            if math.isfinite(float(temperature))
        }
        table = table.rename(columns=renamed)
        feature_names.extend(str(name) for name in table.columns)
        parts.append(table)
    joined = pd.concat(parts, axis=1).reset_index()
    feature_names = _finite_numeric(joined, feature_names)
    return joined, feature_names


def _robust_scaled(frame: pd.DataFrame, columns: list[str]) -> np.ndarray:
    if not columns:
        return np.zeros((len(frame), 0), dtype=float)
    data = np.column_stack(
        [pd.to_numeric(frame[name], errors="coerce").to_numpy(dtype=float) for name in columns]
    )
    median = np.nanmedian(data, axis=0)
    q25 = np.nanpercentile(data, 25.0, axis=0)
    q75 = np.nanpercentile(data, 75.0, axis=0)
    scale = q75 - q25
    span = np.nanmax(data, axis=0) - np.nanmin(data, axis=0)
    scale = np.where(scale > 1.0e-12, scale, span)
    scale = np.where(scale > 1.0e-12, scale, 1.0)
    data = np.where(np.isfinite(data), data, median[None, :])
    return (data - median[None, :]) / scale[None, :]


def _pairwise_distance(values: np.ndarray) -> np.ndarray:
    if values.shape[1] == 0:
        return np.zeros((values.shape[0], values.shape[0]), dtype=float)
    diff = values[:, None, :] - values[None, :, :]
    distance = np.sqrt(np.mean(diff * diff, axis=2))
    finite = distance[np.isfinite(distance) & (distance > 0.0)]
    if finite.size:
        distance /= float(np.median(finite))
    return np.where(np.isfinite(distance), distance, 0.0)


def _quality_scores(objective: pd.Series) -> np.ndarray:
    values = pd.to_numeric(objective, errors="coerce").to_numpy(dtype=float)
    values = np.where(np.isfinite(values), values, np.inf)
    order = np.argsort(values, kind="mergesort")
    scores = np.zeros(len(values), dtype=float)
    if len(values) == 1:
        scores[0] = 1.0
    elif len(values) > 1:
        scores[order] = 1.0 - np.arange(len(values), dtype=float) / (len(values) - 1)
    return scores


def select_quality_diverse(
    candidates: pd.DataFrame,
    temperature_detail: pd.DataFrame,
    config: SelectionConfig,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    cfg = config.validate()
    required = {"candidate_id", "objective"}
    missing = sorted(required.difference(candidates.columns))
    if missing:
        raise ValueError(f"candidate table is missing {missing}")
    work = candidates.copy().reset_index(drop=True)
    work["candidate_id"] = work["candidate_id"].astype(str)
    work = work.sort_values(["objective", "candidate_id"], kind="mergesort")
    work = work.drop_duplicates("candidate_id", keep="first").reset_index(drop=True)
    if work.empty:
        raise ValueError("candidate table is empty")
    if "accepted_for_spatial_promotion" in work:
        passed = _bool_series(work["accepted_for_spatial_promotion"])
    else:
        passed = pd.Series(False, index=work.index)
    work["_passed"] = passed.to_numpy(dtype=bool)

    pass_indices = list(work.index[work["_passed"]])
    fail_indices = list(work.index[~work["_passed"]])
    ordered_indices = pass_indices + fail_indices
    pool_size = min(len(work), max(cfg.count, cfg.pool_factor * cfg.count))
    pool_indices = ordered_indices[:pool_size]
    pool = work.loc[pool_indices].copy().reset_index().rename(columns={"index": "_source_index"})

    parameter_columns = infer_parameter_columns(pool)
    response_table, response_columns = build_response_table(pool, temperature_detail)
    pool = pool.merge(response_table, on="candidate_id", how="left")
    parameter_distance = _pairwise_distance(_robust_scaled(pool, parameter_columns))
    response_distance = _pairwise_distance(_robust_scaled(pool, response_columns))
    diversity_weight = cfg.parameter_weight + cfg.response_weight
    combined_distance = (
        cfg.parameter_weight * parameter_distance
        + cfg.response_weight * response_distance
    ) / diversity_weight
    quality = _quality_scores(pool["objective"])

    selected: list[int] = []
    reasons: dict[int, str] = {}
    utilities: dict[int, float] = {}

    def add(index: int, reason: str, utility: float) -> None:
        if index not in selected and len(selected) < cfg.count:
            selected.append(index)
            reasons[index] = reason
            utilities[index] = float(utility)

    pool_passers = [i for i in pool.index if bool(pool.loc[i, "_passed"])]
    pool_passers.sort(key=lambda i: (float(pool.loc[i, "objective"]), str(pool.loc[i, "candidate_id"])))

    # Hard guarantee: when all passers fit within the promotion budget, retain every
    # one before considering any failed or near-pass candidate.
    if len(pass_indices) <= cfg.count:
        for index in pool_passers:
            add(index, "all_passers_reserve", quality[index])
    else:
        reserve_count = max(1, int(math.ceil(cfg.count * cfg.quality_reserve_fraction)))
        for index in pool_passers[:reserve_count]:
            add(index, "quality_reserve", quality[index])

    if cfg.preserve_restart_lineages and len(selected) < cfg.count and "restart" in pool:
        restart_values = pd.to_numeric(pool["restart"], errors="coerce")
        lineages = sorted({int(value) for value in restart_values if math.isfinite(value)})
        for lineage in lineages:
            eligible = [
                i
                for i in pool_passers
                if i not in selected and int(float(pool.loc[i, "restart"])) == lineage
            ]
            if eligible:
                best = min(eligible, key=lambda i: float(pool.loc[i, "objective"]))
                add(best, "restart_lineage_reserve", quality[best])
            if len(selected) >= cfg.count:
                break

    while len(selected) < min(cfg.count, len(pool)):
        eligible = [i for i in pool.index if i not in selected]
        if not eligible:
            break
        if selected:
            min_distance = np.min(combined_distance[:, selected], axis=1)
        else:
            min_distance = np.ones(len(pool), dtype=float)
        best_index = None
        best_key = None
        for index in eligible:
            utility = cfg.quality_weight * quality[index] + min_distance[index]
            # Passing status is a strict tie-breaker after utility, not a substitute
            # for the hard all-passers guarantee above.
            key = (
                float(utility),
                int(bool(pool.loc[index, "_passed"])),
                float(quality[index]),
                -float(pool.loc[index, "objective"]),
                str(pool.loc[index, "candidate_id"]),
            )
            if best_key is None or key > best_key:
                best_key = key
                best_index = index
        assert best_index is not None
        add(best_index, "quality_diversity_fill", best_key[0])

    selected_rows = pool.loc[selected].copy()
    selected_rows["selection_rank"] = np.arange(1, len(selected_rows) + 1)
    selected_rows["selection_reason"] = [reasons[index] for index in selected]
    selected_rows["selection_utility"] = [utilities[index] for index in selected]
    selected_rows["selection_quality_score"] = [quality[index] for index in selected]
    if selected:
        selected_rows["selection_min_distance"] = [
            0.0
            if rank == 0
            else float(np.min(combined_distance[index, selected[:rank]]))
            for rank, index in enumerate(selected)
        ]
    else:
        selected_rows["selection_min_distance"] = []
    drop_columns = [name for name in selected_rows if name.startswith("response_")]
    drop_columns += ["_source_index", "_passed"]
    selected_rows = selected_rows.drop(columns=drop_columns, errors="ignore")

    selected_ids = set(selected_rows["candidate_id"].astype(str))
    all_passers_retained = set(work.loc[work["_passed"], "candidate_id"].astype(str)).issubset(selected_ids)
    audit = {
        "schema": MODEL_ID,
        "config": asdict(cfg),
        "n_candidates": int(len(work)),
        "n_passers": int(work["_passed"].sum()),
        "n_pool": int(len(pool)),
        "n_selected": int(len(selected_rows)),
        "all_passers_fit_in_budget": bool(int(work["_passed"].sum()) <= cfg.count),
        "all_passers_retained": bool(all_passers_retained),
        "parameter_columns": parameter_columns,
        "response_columns": response_columns,
        "selected": selected_rows[
            [
                "candidate_id",
                "objective",
                "selection_rank",
                "selection_reason",
                "selection_utility",
                "selection_quality_score",
                "selection_min_distance",
            ]
        ].to_dict(orient="records"),
        "constitutive_physics_modified": False,
        "mechanics_closure": "v9_isotropic_moving_process_zone",
        "requires_spatial_promotion": True,
        "requires_2d_validation_after_spatial_promotion": True,
    }
    return selected_rows.reset_index(drop=True), audit


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--search-root", type=Path, required=True)
    parser.add_argument("--target-class", default="DBTT")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--count", type=int, default=10)
    parser.add_argument("--quality-reserve-fraction", type=float, default=0.30)
    parser.add_argument("--quality-weight", type=float, default=0.35)
    parser.add_argument("--parameter-weight", type=float, default=0.45)
    parser.add_argument("--response-weight", type=float, default=0.55)
    parser.add_argument("--pool-factor", type=int, default=12)
    parser.add_argument("--no-preserve-restart-lineages", action="store_true")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    source = args.search_root / args.target_class
    candidate_path = source / "unified_global_all_candidates.csv"
    detail_path = source / "unified_global_temperature_detail.csv"
    if not candidate_path.is_file():
        raise SystemExit(f"candidate table not found: {candidate_path}")
    if not detail_path.is_file():
        raise SystemExit(f"temperature detail not found: {detail_path}")
    out = args.out / args.target_class
    if out.exists() and not args.force:
        raise SystemExit(f"refusing to overwrite {out}; pass --force explicitly")
    out.mkdir(parents=True, exist_ok=True)

    candidates = pd.read_csv(candidate_path)
    details = pd.read_csv(detail_path)
    selected, audit = select_quality_diverse(
        candidates,
        details,
        SelectionConfig(
            count=args.count,
            quality_reserve_fraction=args.quality_reserve_fraction,
            quality_weight=args.quality_weight,
            parameter_weight=args.parameter_weight,
            response_weight=args.response_weight,
            pool_factor=args.pool_factor,
            preserve_restart_lineages=not args.no_preserve_restart_lineages,
        ),
    )
    selected.to_csv(out / "spatial_promotion_manifest.csv", index=False)
    selected.to_csv(out / "quality_diversity_selected.csv", index=False)
    selected_ids = set(selected["candidate_id"].astype(str))
    details[details["candidate_id"].astype(str).isin(selected_ids)].to_csv(
        out / "quality_diversity_temperature_detail.csv", index=False
    )
    (out / "quality_diversity_selection.json").write_text(
        json.dumps(audit, indent=2, allow_nan=False)
    )
    print(json.dumps(audit, indent=2), flush=True)


if __name__ == "__main__":
    main()
