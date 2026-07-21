#!/usr/bin/env python3
"""Build a multi-fidelity 0-D -> 1-D training table for MPZ v9.12.

The table preserves the full 0-D temperature trajectory as input and learns the
spatial transfer rather than fitting a standalone 1-D response from raw
parameters alone.  Extreme 0-D amplitudes are retained and represented through
log-transformed features; no upper-amplitude rejection is applied.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--candidate-registry", required=True)
    p.add_argument("--zero-d-root", required=True)
    p.add_argument("--one-d-root", required=True)
    p.add_argument("--bounds-json", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--transfer-epsilon", type=float, default=1.0e-6)
    p.add_argument("--one-d-threshold", type=float, default=50.0)
    return p.parse_args()


def read_csv(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open(newline="") as fp:
        return list(csv.DictReader(fp))


def read_status(root: Path) -> dict[str, dict[str, str]]:
    for name in ("resilient_records_checkpoint.csv", "ranking.csv"):
        path = root / name
        if path.exists():
            return {row["candidate_id"]: row for row in read_csv(path)}
    return {}


def read_summary(root: Path, candidate_id: str) -> dict[str, Any] | None:
    path = root / candidate_id / "candidate_summary.json"
    if not path.exists():
        return None
    return json.loads(path.read_text())


def trajectory(summary: dict[str, Any] | None) -> tuple[np.ndarray, np.ndarray]:
    if summary is None:
        return np.asarray([], dtype=float), np.asarray([], dtype=float)
    temperatures = np.asarray(summary.get("temperatures_K", []), dtype=float)
    values = np.asarray(
        summary.get("developed_delta_K_micro_MPa_sqrt_m", []), dtype=float
    )
    if temperatures.size != values.size:
        raise RuntimeError("temperature and developed-response lengths differ")
    order = np.argsort(temperatures)
    return temperatures[order], values[order]


def trajectory_metrics(T: np.ndarray, y: np.ndarray) -> dict[str, float]:
    finite = np.isfinite(T) & np.isfinite(y)
    T = T[finite]
    y = y[finite]
    if y.size == 0:
        return {
            "amplitude": math.nan,
            "jump_temperature_K": math.nan,
            "largest_positive_jump": math.nan,
            "localization": math.nan,
            "persistence": math.nan,
            "post_peak_drop_fraction": math.nan,
            "final": math.nan,
            "peak": math.nan,
        }
    ymin = float(np.min(y))
    ymax = float(np.max(y))
    amplitude = ymax - ymin
    positive = np.maximum(np.diff(y), 0.0)
    if positive.size and float(np.sum(positive)) > 0.0:
        j = int(np.argmax(positive))
        jump_temperature = float(T[j + 1])
        largest_positive_jump = float(positive[j])
        localization = float(positive[j] / np.sum(positive))
    else:
        jump_temperature = math.nan
        largest_positive_jump = 0.0
        localization = 0.0
    peak = max(ymax, 0.0)
    final = float(y[-1])
    persistence = final / peak if peak > 0.0 else 0.0
    post_peak_drop = max(ymax - final, 0.0)
    post_peak_drop_fraction = post_peak_drop / peak if peak > 0.0 else 0.0
    return {
        "amplitude": amplitude,
        "jump_temperature_K": jump_temperature,
        "largest_positive_jump": largest_positive_jump,
        "localization": localization,
        "persistence": persistence,
        "post_peak_drop_fraction": post_peak_drop_fraction,
        "final": final,
        "peak": ymax,
    }


def objective_value(summary: dict[str, Any] | None, key: str) -> float:
    if summary is None:
        return math.nan
    value = summary.get("objective", {}).get(key, math.nan)
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def log1p_nonnegative(value: float) -> float:
    return math.log1p(max(float(value), 0.0)) if math.isfinite(value) else math.nan


def main() -> int:
    a = parse_args()
    rows = read_csv(a.candidate_registry)
    if not rows:
        raise RuntimeError("candidate registry is empty")

    bounds_payload = json.loads(Path(a.bounds_json).read_text())
    bound_names = list(bounds_payload.get("search_bounds", bounds_payload))
    root0 = Path(a.zero_d_root)
    root1 = Path(a.one_d_root)
    status0 = read_status(root0)
    status1 = read_status(root1)

    records: list[dict[str, Any]] = []
    all_temperatures: set[float] = set()
    for row in rows:
        candidate_id = str(row["candidate_id"])
        summary0 = read_summary(root0, candidate_id)
        summary1 = read_summary(root1, candidate_id)
        T0, y0 = trajectory(summary0)
        T1, y1 = trajectory(summary1)
        all_temperatures.update(float(v) for v in T0)
        all_temperatures.update(float(v) for v in T1)
        metrics0 = trajectory_metrics(T0, y0)
        metrics1 = trajectory_metrics(T1, y1)

        rec: dict[str, Any] = {"candidate_id": candidate_id}
        for name in bound_names:
            rec[f"x_raw__{name}"] = pd.to_numeric(row.get(name), errors="coerce")

        for T, value in zip(T0, y0):
            tag = int(round(float(T)))
            rec[f"x_0d__deltaK_T{tag}K"] = float(value)
            rec[f"x_0d__log1p_deltaK_T{tag}K"] = log1p_nonnegative(float(value))

        for key, value in metrics0.items():
            rec[f"x_0d__{key}"] = value
        rec["x_0d__log1p_amplitude"] = log1p_nonnegative(metrics0["amplitude"])
        rec["x_0d__log1p_max_shield"] = log1p_nonnegative(
            objective_value(summary0, "max_abs_K_shield_MPa_sqrt_m")
        )
        rec["x_0d__log1p_max_gnd"] = log1p_nonnegative(
            objective_value(summary0, "max_gnd_abs_line_count_per_unit_thickness")
        )
        rec["x_0d__min_source_available_fraction"] = objective_value(
            summary0, "min_source_available_fraction"
        )

        status_one = status1.get(candidate_id, {}).get(
            "status", "complete" if summary1 is not None else "missing"
        )
        complete_one = status_one == "complete" and summary1 is not None
        rec["y__status_1d"] = status_one
        rec["y__complete_1d"] = bool(complete_one)

        if complete_one:
            for T, value in zip(T1, y1):
                tag = int(round(float(T)))
                rec[f"y__deltaK_1d_T{tag}K"] = float(value)
                rec[f"y__log1p_deltaK_1d_T{tag}K"] = log1p_nonnegative(float(value))
            for key, value in metrics1.items():
                rec[f"y__{key}_1d"] = value
            amp0 = max(float(metrics0["amplitude"]), 0.0)
            amp1 = max(float(metrics1["amplitude"]), 0.0)
            eps = max(float(a.transfer_epsilon), 1.0e-300)
            rec["y__log1p_amplitude_1d"] = math.log1p(amp1)
            rec["y__transfer_ratio_1d_over_0d"] = (amp1 + eps) / (amp0 + eps)
            rec["y__log10_transfer_ratio_1d_over_0d"] = math.log10(
                (amp1 + eps) / (amp0 + eps)
            )
            rec["y__retains_50_MPa_sqrt_m_1d"] = bool(
                amp1 >= float(a.one_d_threshold)
            )
            rec["y__log1p_max_shield_1d"] = log1p_nonnegative(
                objective_value(summary1, "max_abs_K_shield_MPa_sqrt_m")
            )
            rec["y__max_tau_gnd_tip_MPa_1d"] = objective_value(
                summary1, "max_tau_gnd_tip_MPa"
            )
            rec["y__min_source_available_fraction_1d"] = objective_value(
                summary1, "min_source_available_fraction"
            )
        records.append(rec)

    temperatures = sorted(all_temperatures)
    fields: list[str] = []
    for rec in records:
        for key in rec:
            if key not in fields:
                fields.append(key)
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(records, columns=fields).to_csv(out, index=False)
    complete_count = sum(bool(r.get("y__complete_1d", False)) for r in records)
    print(
        "TRANSFER_TABLE "
        f"rows={len(records)} complete_1d={complete_count} "
        f"temperatures={','.join(str(int(t)) for t in temperatures)} out={out}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
