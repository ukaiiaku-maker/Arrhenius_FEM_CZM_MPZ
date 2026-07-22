#!/usr/bin/env python3
"""Add direction-aware DBTT and genuine-peak targets to a v9.12 transfer table."""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--low-max-K", type=float, default=700.0)
    p.add_argument("--high-min-K", type=float, default=1000.0)
    p.add_argument("--peak-min-K", type=float, default=800.0)
    p.add_argument("--peak-max-K", type=float, default=1000.0)
    p.add_argument("--direction-threshold", type=float, default=5.0)
    p.add_argument("--peak-threshold", type=float, default=1.0)
    return p.parse_args()


def temperature_columns(df: pd.DataFrame, prefix: str) -> list[tuple[float, str]]:
    pattern = re.compile(re.escape(prefix) + r"T([0-9]+)K$")
    found: list[tuple[float, str]] = []
    for column in df.columns:
        match = pattern.match(column)
        if match:
            found.append((float(match.group(1)), column))
    return sorted(found)


def add_trajectory_metrics(
    df: pd.DataFrame,
    *,
    prefix: str,
    out_prefix: str,
    low_max_K: float,
    high_min_K: float,
    peak_min_K: float,
    peak_max_K: float,
) -> None:
    columns = temperature_columns(df, prefix)
    if not columns:
        return
    temperatures = np.asarray([item[0] for item in columns], dtype=float)
    values = df[[item[1] for item in columns]].apply(
        pd.to_numeric, errors="coerce"
    ).to_numpy(float)

    low_mask = temperatures <= float(low_max_K)
    high_mask = temperatures >= float(high_min_K)
    prepeak_mask = temperatures < float(peak_min_K)
    if not np.any(low_mask) or not np.any(high_mask):
        raise RuntimeError("trajectory grid does not span requested low/high bands")
    if not np.any(prepeak_mask):
        prepeak_mask = low_mask

    low = np.nanmedian(values[:, low_mask], axis=1)
    high = np.nanmedian(values[:, high_mask], axis=1)
    gain = high - low
    gain_positive = np.maximum(gain, 0.0)

    safe = np.where(np.isfinite(values), values, -np.inf)
    peak_index = np.argmax(safe, axis=1)
    peak = safe[np.arange(len(df)), peak_index]
    peak[~np.isfinite(peak)] = np.nan
    peak_temperature = temperatures[peak_index]
    final = values[:, -1]
    prepeak_baseline = np.nanmedian(values[:, prepeak_mask], axis=1)
    peak_rise = peak - prepeak_baseline
    peak_drop = peak - final
    peak_in_window = (
        (peak_temperature >= float(peak_min_K))
        & (peak_temperature <= float(peak_max_K))
    )
    peak_prominence = np.where(
        peak_in_window,
        np.maximum(np.minimum(peak_rise, peak_drop), 0.0),
        0.0,
    )

    differences = np.diff(values, axis=1)
    positive = np.maximum(differences, 0.0)
    jump_index = np.argmax(np.where(np.isfinite(positive), positive, -np.inf), axis=1)
    jump_high_temperature = temperatures[np.minimum(jump_index + 1, len(temperatures) - 1)]
    persistence = np.divide(
        final,
        peak,
        out=np.zeros_like(final),
        where=np.isfinite(peak) & (peak > 0.0),
    )

    df[f"{out_prefix}low_temperature_baseline"] = low
    df[f"{out_prefix}high_temperature_plateau"] = high
    df[f"{out_prefix}directional_dbtt_gain"] = gain
    df[f"{out_prefix}directional_dbtt_gain_positive"] = gain_positive
    df[f"{out_prefix}log1p_directional_dbtt_gain_positive"] = np.log1p(
        gain_positive
    )
    df[f"{out_prefix}peak_temperature_K"] = peak_temperature
    df[f"{out_prefix}peak_rise"] = peak_rise
    df[f"{out_prefix}peak_drop"] = peak_drop
    df[f"{out_prefix}peak_prominence"] = peak_prominence
    df[f"{out_prefix}log1p_peak_prominence"] = np.log1p(peak_prominence)
    df[f"{out_prefix}largest_positive_jump_high_temperature_K"] = (
        jump_high_temperature
    )
    df[f"{out_prefix}persistence_from_trajectory"] = persistence


def main() -> int:
    a = parse_args()
    df = pd.read_csv(a.input, low_memory=False)

    add_trajectory_metrics(
        df,
        prefix="x_0d__deltaK_",
        out_prefix="x_0d__",
        low_max_K=a.low_max_K,
        high_min_K=a.high_min_K,
        peak_min_K=a.peak_min_K,
        peak_max_K=a.peak_max_K,
    )
    add_trajectory_metrics(
        df,
        prefix="y__deltaK_1d_",
        out_prefix="y__",
        low_max_K=a.low_max_K,
        high_min_K=a.high_min_K,
        peak_min_K=a.peak_min_K,
        peak_max_K=a.peak_max_K,
    )

    if "y__directional_dbtt_gain" in df:
        direction_correct = (
            pd.to_numeric(df["y__directional_dbtt_gain"], errors="coerce") > 0.0
        )
        if "y__peak_temperature_K" in df:
            direction_correct &= (
                pd.to_numeric(df["y__peak_temperature_K"], errors="coerce")
                >= a.peak_min_K
            )
        if "y__largest_positive_jump_high_temperature_K" in df:
            direction_correct &= (
                pd.to_numeric(
                    df["y__largest_positive_jump_high_temperature_K"],
                    errors="coerce",
                )
                >= a.peak_min_K
            )
        if "y__persistence_from_trajectory" in df:
            direction_correct &= (
                pd.to_numeric(
                    df["y__persistence_from_trajectory"], errors="coerce"
                )
                >= 0.70
            )
        df["y__direction_correct_1d"] = direction_correct
        df["y__directional_dbtt_ge_threshold_1d"] = (
            direction_correct
            & (
                pd.to_numeric(
                    df["y__directional_dbtt_gain"], errors="coerce"
                )
                >= a.direction_threshold
            )
        )

    if "y__peak_prominence" in df:
        df["y__peak_like_1d"] = (
            pd.to_numeric(df["y__peak_prominence"], errors="coerce")
            >= a.peak_threshold
        )

    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    print(
        "DIRECTIONAL_PEAK_TABLE "
        f"rows={len(df)} columns={len(df.columns)} out={out}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
