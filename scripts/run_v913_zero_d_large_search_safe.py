#!/usr/bin/env python3
"""Safety wrapper for the v9.13 persistent zero-D large search.

The vectorized proxy deliberately explores parameter regions that may produce no
finite toughness values. Those candidates must receive a finite rejection score,
not emit NumPy all-NaN warnings or leak NaN/Inf into strict JSON outputs.

The exact-stage gate also requires completion to the requested checkpoint. The
promotion selector preserves every complete strict pass and fills the remaining
pool tier-by-tier with diverse, complete, finite local-peak candidates.
"""
from __future__ import annotations

import copy
import math
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

import arrhenius_fracture.zero_d_search_v913 as search_helpers
import scripts.run_v913_zero_d_large_search as base


SAFETY_SCHEMA = "v9.13_persistent_zero_d_nonfinite_completion_promotion_v2"


def json_safe(value: Any) -> Any:
    """Recursively convert NumPy values and nonfinite floats to JSON null."""
    if isinstance(value, np.generic):
        return json_safe(value.item())
    if isinstance(value, Mapping):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def local_peak_metrics_safe(
    temperatures_K: Any,
    values: Any,
    *,
    desired_min_K: float = 850.0,
    desired_max_K: float = 1100.0,
) -> dict[str, float | bool]:
    """Return finite-or-NaN metrics without sentinel infinities."""
    temperatures = np.asarray(temperatures_K, dtype=float)
    response = np.asarray(values, dtype=float)
    finite = np.isfinite(temperatures) & np.isfinite(response)
    if int(np.sum(finite)) < 3:
        return {
            "peak_temperature_K": float("nan"),
            "peak_value": float("nan"),
            "two_sided_prominence": float("nan"),
            "post_peak_drop": float("nan"),
            "high_temperature_rebound": float("nan"),
            "peak_internal": False,
            "peak_in_desired_window": False,
        }

    order = np.argsort(temperatures[finite])
    temperatures = temperatures[finite][order]
    response = response[finite][order]
    internal_indices = np.arange(1, len(temperatures) - 1)
    local_indices = internal_indices[
        (response[internal_indices] > response[internal_indices - 1])
        & (response[internal_indices] > response[internal_indices + 1])
    ]
    if local_indices.size:
        prominence_candidates = np.minimum(
            response[local_indices] - response[local_indices - 1],
            response[local_indices] - response[local_indices + 1],
        )
        best = int(local_indices[int(np.argmax(prominence_candidates))])
        peak_internal = True
    else:
        best = int(np.argmax(response))
        peak_internal = 0 < best < len(response) - 1

    if peak_internal:
        prominence = float(
            min(
                response[best] - response[best - 1],
                response[best] - response[best + 1],
            )
        )
    else:
        prominence = float("nan")

    if best + 1 < len(response):
        post = response[best + 1 :]
        post_minimum = float(np.min(post))
        post_maximum = float(np.max(post))
        drop = float(response[best] - post_minimum)
        rebound = float(post_maximum - response[best])
    else:
        drop = float("nan")
        rebound = float("nan")

    return {
        "peak_temperature_K": float(temperatures[best]),
        "peak_value": float(response[best]),
        "two_sided_prominence": prominence,
        "post_peak_drop": drop,
        "high_temperature_rebound": rebound,
        "peak_internal": bool(peak_internal),
        "peak_in_desired_window": bool(
            peak_internal
            and desired_min_K <= temperatures[best] <= desired_max_K
        ),
    }


def curve_metrics_matrix_safe(
    temperatures: np.ndarray,
    curves: np.ndarray,
    *,
    peak_min: float,
    peak_max: float,
) -> dict[str, np.ndarray]:
    """Vectorized local-peak metrics with explicit all-nonfinite handling."""
    temperatures = np.asarray(temperatures, dtype=float)
    curves = np.asarray(curves, dtype=float)
    if curves.ndim != 2 or temperatures.ndim != 1:
        raise ValueError("temperature grid must be 1-D and curves must be 2-D")
    if curves.shape[1] != temperatures.size:
        raise ValueError("curve width must match temperature-grid length")

    nrows, ntemperatures = curves.shape
    finite = np.isfinite(curves)
    has_finite = np.any(finite, axis=1)
    peak_index = np.zeros(nrows, dtype=int)
    prominence = np.full(nrows, np.nan, dtype=float)
    best_local_prominence = np.full(nrows, -np.inf, dtype=float)
    local_found = np.zeros(nrows, dtype=bool)

    for index in range(1, ntemperatures - 1):
        triplet_finite = (
            finite[:, index - 1] & finite[:, index] & finite[:, index + 1]
        )
        local = (
            triplet_finite
            & (curves[:, index] > curves[:, index - 1])
            & (curves[:, index] > curves[:, index + 1])
        )
        candidate_prominence = np.minimum(
            curves[:, index] - curves[:, index - 1],
            curves[:, index] - curves[:, index + 1],
        )
        replace = local & (candidate_prominence > best_local_prominence)
        peak_index[replace] = index
        best_local_prominence[replace] = candidate_prominence[replace]
        local_found[replace] = True

    fallback = has_finite & ~local_found
    if np.any(fallback):
        sanitized = np.where(finite[fallback], curves[fallback], -np.inf)
        peak_index[fallback] = np.argmax(sanitized, axis=1)

    rows = np.arange(nrows)
    peak_value = np.full(nrows, np.nan, dtype=float)
    peak_temperature = np.full(nrows, np.nan, dtype=float)
    peak_value[has_finite] = curves[rows[has_finite], peak_index[has_finite]]
    peak_temperature[has_finite] = temperatures[peak_index[has_finite]]
    peak_internal = (
        has_finite
        & (peak_index > 0)
        & (peak_index < ntemperatures - 1)
    )
    prominence[local_found] = best_local_prominence[local_found]

    post_minimum = np.full(nrows, np.nan, dtype=float)
    post_maximum = np.full(nrows, np.nan, dtype=float)
    for row, index in enumerate(peak_index):
        if not has_finite[row] or index + 1 >= ntemperatures:
            continue
        post = curves[row, index + 1 :]
        finite_post = post[np.isfinite(post)]
        if finite_post.size:
            post_minimum[row] = float(np.min(finite_post))
            post_maximum[row] = float(np.max(finite_post))

    drop = peak_value - post_minimum
    rebound = post_maximum - peak_value
    desired = (
        peak_internal
        & (peak_temperature >= float(peak_min))
        & (peak_temperature <= float(peak_max))
    )
    return {
        "peak_temperature_K": peak_temperature,
        "peak_value_MPa_sqrt_m": peak_value,
        "two_sided_prominence_MPa_sqrt_m": prominence,
        "post_peak_drop_MPa_sqrt_m": drop,
        "high_temperature_rebound_MPa_sqrt_m": rebound,
        "peak_internal": peak_internal.astype(int),
        "peak_in_desired_window": desired.astype(int),
    }


_original_exact_candidate_worker = base._exact_candidate_worker
_original_run_contract = base._run_contract
_original_score_frame = base._score_frame


def exact_candidate_worker_safe(row: dict[str, Any]) -> dict[str, Any]:
    """Normalize nonfinite rejection diagnostics before parent serialization."""
    return json_safe(_original_exact_candidate_worker(row))


def score_frame_safe(
    frame: pd.DataFrame,
    prefix: str,
    *,
    minimum_prominence: float,
    minimum_drop: float,
    maximum_rebound: float,
    peak_min: float,
    peak_max: float,
) -> pd.DataFrame:
    """Require exact-stage completion and replace nonfinite objectives."""
    out = _original_score_frame(
        frame,
        prefix,
        minimum_prominence=minimum_prominence,
        minimum_drop=minimum_drop,
        maximum_rebound=maximum_rebound,
        peak_min=peak_min,
        peak_max=peak_max,
    )
    objective_name = f"{prefix}_objective"
    gate_name = f"{prefix}_gate_pass"
    objective = pd.to_numeric(out[objective_name], errors="coerce").to_numpy(float)
    objective = np.where(np.isfinite(objective), objective, 1000.0)

    completion_name = f"{prefix}_complete"
    if completion_name in out:
        complete = (
            pd.to_numeric(out[completion_name], errors="coerce")
            .fillna(0)
            .astype(bool)
            .to_numpy()
        )
        objective = objective + np.where(complete, 0.0, 500.0)
        out[gate_name] = out[gate_name].astype(bool).to_numpy() & complete

    out[objective_name] = objective
    return out


def _finite_columns(frame: pd.DataFrame, names: Sequence[str]) -> np.ndarray:
    mask = np.ones(len(frame), dtype=bool)
    for name in names:
        values = pd.to_numeric(frame[name], errors="coerce").to_numpy(float)
        mask &= np.isfinite(values)
    return mask


def _greedy_quality_diverse(
    features: np.ndarray,
    objective: np.ndarray,
    candidates: np.ndarray,
    selected: list[int],
    count: int,
) -> list[int]:
    """Add quality-weighted farthest points from one promotion tier."""
    remaining = int(count)
    available = [int(index) for index in candidates if int(index) not in selected]
    if remaining <= 0 or not available:
        return selected

    if not selected:
        best = min(
            available,
            key=lambda index: (
                float(objective[index])
                if np.isfinite(objective[index])
                else float("inf"),
                index,
            ),
        )
        selected.append(best)
        available.remove(best)
        remaining -= 1
        if remaining <= 0 or not available:
            return selected

    min_distance = np.full(len(features), np.inf, dtype=float)
    for index in selected:
        distance = np.sum(np.square(features - features[index]), axis=1)
        min_distance = np.minimum(min_distance, distance)

    available_array = np.asarray(available, dtype=int)
    available_objective = objective[available_array]
    order = np.argsort(np.argsort(available_objective, kind="stable"), kind="stable")
    percentile = (order + 1.0) / max(len(order), 1)
    quality = np.exp(-3.0 * percentile)

    for _ in range(min(remaining, len(available_array))):
        distance = np.sqrt(np.maximum(min_distance[available_array], 0.0))
        acquisition = distance * (0.10 + 0.90 * quality)
        pick_position = int(np.argmax(acquisition))
        picked = int(available_array[pick_position])
        selected.append(picked)
        distance_to_pick = np.sum(
            np.square(features - features[picked]),
            axis=1,
        )
        min_distance = np.minimum(min_distance, distance_to_pick)
        available_array = np.delete(available_array, pick_position)
        quality = np.delete(quality, pick_position)
        if available_array.size == 0:
            break
    return selected


def diverse_selection_safe(
    ranked: pd.DataFrame,
    policy: Mapping[str, Any],
    count: int,
) -> pd.DataFrame:
    """Preserve complete strict passes, then fill ordered quality tiers."""
    if count < 1:
        raise ValueError("promotion count must be positive")
    frame = ranked.reset_index(drop=True).copy()

    parameter_finite = _finite_columns(
        frame,
        tuple(base.ACTIVE_CANDIDATE_PARAMETER_FIELDS),
    )
    metric_finite = _finite_columns(
        frame,
        (
            "zeroD_objective",
            "zeroD_peak_temperature_K",
            "zeroD_peak_value_MPa_sqrt_m",
            "zeroD_two_sided_prominence_MPa_sqrt_m",
            "zeroD_post_peak_drop_MPa_sqrt_m",
            "zeroD_high_temperature_rebound_MPa_sqrt_m",
        ),
    )
    complete = (
        pd.to_numeric(frame["zeroD_complete"], errors="coerce")
        .fillna(0)
        .astype(bool)
        .to_numpy()
    )
    internal = frame["zeroD_peak_internal"].astype(bool).to_numpy()
    desired = frame["zeroD_peak_in_desired_window"].astype(bool).to_numpy()
    strict = (
        frame["zeroD_gate_pass"].astype(bool).to_numpy()
        & complete
        & parameter_finite
        & metric_finite
    )

    tier = np.full(len(frame), "", dtype=object)
    tier[strict] = "strict_gate"
    mask = (
        (tier == "")
        & complete
        & parameter_finite
        & metric_finite
        & internal
        & desired
    )
    tier[mask] = "relaxed_desired_peak"
    mask = (
        (tier == "")
        & complete
        & parameter_finite
        & metric_finite
        & internal
    )
    tier[mask] = "internal_peak_outside_window"
    mask = (
        (tier == "")
        & complete
        & parameter_finite
        & metric_finite
    )
    tier[mask] = "finite_complete_boundary"

    eligible = np.flatnonzero(tier != "")
    if eligible.size == 0:
        raise RuntimeError("no complete finite zero-D candidates are promotable")

    features = base._normalize_features(frame, policy)
    objective = pd.to_numeric(
        frame["zeroD_objective"], errors="coerce"
    ).to_numpy(float)
    selected: list[int] = []

    for tier_name in (
        "strict_gate",
        "relaxed_desired_peak",
        "internal_peak_outside_window",
        "finite_complete_boundary",
    ):
        remaining = int(count) - len(selected)
        if remaining <= 0:
            break
        candidates = np.flatnonzero(tier == tier_name)
        if tier_name == "strict_gate" and len(candidates) <= remaining:
            ordered = sorted(
                (int(index) for index in candidates),
                key=lambda index: (objective[index], index),
            )
            selected.extend(ordered)
        else:
            selected = _greedy_quality_diverse(
                features,
                objective,
                candidates,
                selected,
                remaining,
            )

    if len(selected) < min(int(count), len(eligible)):
        raise RuntimeError(
            "promotion selector failed to fill the available complete finite pool"
        )

    result = frame.iloc[selected[:count]].copy()
    result["promotion_tier"] = tier[np.asarray(selected[:count], dtype=int)]
    result["diversity_rank"] = np.arange(1, len(result) + 1)

    strict_ids = set(frame.loc[strict, "candidate_id"].astype(str))
    promoted_ids = set(result["candidate_id"].astype(str))
    if len(strict_ids) <= count and not strict_ids.issubset(promoted_ids):
        missing = sorted(strict_ids - promoted_ids)
        raise RuntimeError(f"promotion omitted complete strict passes: {missing}")
    if not result["zeroD_complete"].astype(bool).all():
        raise RuntimeError("promotion contains incomplete exact candidates")
    return result


def run_contract_safe(
    args: Any,
    policy: Mapping[str, Any],
    anchor_rows: Any,
) -> dict[str, Any]:
    """Version safety and promotion behavior for resume protection."""
    original = _original_run_contract(args, policy, anchor_rows)
    contract = copy.deepcopy(original["contract"])
    contract["nonfinite_safety"] = {
        "schema": SAFETY_SCHEMA,
        "json_nonfinite_encoding": "null",
        "all_nonfinite_curve_policy": "invalid_candidate_finite_penalty",
        "nan_slice_warnings_suppressed_by_explicit_finite_mask": True,
        "exact_gate_requires_completion": True,
        "promotion_requires_complete_finite_metrics": True,
        "all_complete_strict_passes_preserved": True,
        "promotion_tiers": [
            "strict_gate",
            "relaxed_desired_peak",
            "internal_peak_outside_window",
            "finite_complete_boundary",
        ],
    }
    stable = copy.deepcopy(contract)
    stable.pop("created_at_utc", None)
    return {
        "sha256": base._canonical_sha256(stable),
        "contract": contract,
    }


def install_safety_patch() -> None:
    base._json_safe = json_safe
    base.local_peak_metrics = local_peak_metrics_safe
    base._curve_metrics_matrix = curve_metrics_matrix_safe
    search_helpers._curve_metrics_matrix = curve_metrics_matrix_safe
    base._exact_candidate_worker = exact_candidate_worker_safe
    base._score_frame = score_frame_safe
    base._diverse_selection = diverse_selection_safe
    base._run_contract = run_contract_safe


def main() -> int:
    install_safety_patch()
    return int(base.main())


if __name__ == "__main__":
    raise SystemExit(main())
