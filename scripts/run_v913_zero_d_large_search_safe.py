#!/usr/bin/env python3
"""Safety wrapper for the v9.13 persistent zero-D large search.

The vectorized proxy deliberately explores parameter regions that may produce no
finite toughness values.  Those candidates must receive a finite rejection score,
not emit NumPy all-NaN warnings or leak NaN/Inf into strict JSON outputs.
"""
from __future__ import annotations

import copy
import math
from typing import Any, Mapping

import numpy as np

import scripts.run_v913_zero_d_large_search as base


SAFETY_SCHEMA = "v9.13_persistent_zero_d_nonfinite_safety_v1"


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
        triplet_finite = finite[:, index - 1] & finite[:, index] & finite[:, index + 1]
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


def exact_candidate_worker_safe(row: dict[str, Any]) -> dict[str, Any]:
    """Normalize nonfinite rejection diagnostics before parent serialization."""
    return json_safe(_original_exact_candidate_worker(row))


def run_contract_safe(
    args: Any,
    policy: Mapping[str, Any],
    anchor_rows: Any,
) -> dict[str, Any]:
    """Version the safety behavior so unsafe output directories cannot resume."""
    original = _original_run_contract(args, policy, anchor_rows)
    contract = copy.deepcopy(original["contract"])
    contract["nonfinite_safety"] = {
        "schema": SAFETY_SCHEMA,
        "json_nonfinite_encoding": "null",
        "all_nonfinite_curve_policy": "invalid_candidate_finite_penalty",
        "nan_slice_warnings_suppressed_by_explicit_finite_mask": True,
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
    base._exact_candidate_worker = exact_candidate_worker_safe
    base._run_contract = run_contract_safe


def main() -> int:
    install_safety_patch()
    return int(base.main())


if __name__ == "__main__":
    raise SystemExit(main())
