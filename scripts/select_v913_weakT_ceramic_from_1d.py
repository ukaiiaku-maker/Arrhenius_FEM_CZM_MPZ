#!/usr/bin/env python3
"""Select weak-temperature/FCC-like and ceramic-like v9.13 candidates.

The selector operates on the completed 384-candidate autonomous 1-D screen.  It
uses the complete event histories, not only the peak-search labels:

* weak-T/FCC-like: small K50 temperature span with a small positive R-curve rise;
* ceramic-like: nearly flat lower/mid-temperature response, weak toughness loss
  concentrated at high temperature, and negligible R-curve development.

The input candidate rows are immutable.  The selected registry copies the exact
active rows and adds only selection metadata columns.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import sys
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import numpy as np
import pandas as pd

from arrhenius_fracture.dbtt_long_alignment_v913 import (
    checkpoint_from_events,
    checkpoint_reached,
    peak_metrics,
)

EXPECTED_TEMPERATURES_K = (700, 800, 900, 950, 1000, 1050, 1100, 1200, 1300, 1400)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--case-root", type=Path)
    parser.add_argument("--candidate-registry", type=Path)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--top-count", type=int, default=20)
    return parser.parse_args()


def read_registry(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(path)
    frame = pd.read_csv(path)
    if "candidate_id" not in frame:
        raise RuntimeError(f"candidate registry lacks candidate_id: {path}")
    if frame["candidate_id"].astype(str).duplicated().any():
        raise RuntimeError("candidate registry contains duplicate candidate_id values")
    return frame


def locate_inputs(args: argparse.Namespace) -> tuple[Path, Path]:
    run_dir = args.run_dir.expanduser().resolve()
    case_root = (
        args.case_root.expanduser().resolve()
        if args.case_root is not None
        else run_dir / "one_d_screen" / "cases"
    )
    registry = (
        args.candidate_registry.expanduser().resolve()
        if args.candidate_registry is not None
        else run_dir / "preparation" / "selected_384_registry.csv"
    )
    if not case_root.is_dir():
        raise FileNotFoundError(case_root)
    if not registry.is_file():
        raise FileNotFoundError(registry)
    return case_root, registry


def load_payloads(case_root: Path) -> list[dict[str, Any]]:
    paths = sorted(case_root.glob("*/T*K.json"))
    if not paths:
        raise RuntimeError(f"no candidate/T*K.json payloads under {case_root}")
    payloads: list[dict[str, Any]] = []
    for path in paths:
        payload = json.loads(path.read_text())
        payload["_source_path"] = str(path)
        payloads.append(payload)
    return payloads


def finite_median(values: list[float]) -> float:
    clean = [float(value) for value in values if math.isfinite(float(value))]
    return float(np.median(clean)) if clean else float("nan")


def finite_std(values: list[float]) -> float:
    clean = [float(value) for value in values if math.isfinite(float(value))]
    return float(np.std(clean)) if clean else float("nan")


def first_event_K(events: list[dict[str, Any]]) -> float:
    if not events:
        return float("nan")
    return float(events[0]["K_MPa_sqrt_m"])


def event_range_K(events: list[dict[str, Any]]) -> float:
    values = np.asarray([float(row["K_MPa_sqrt_m"]) for row in events], dtype=float)
    values = values[np.isfinite(values)]
    return float(np.max(values) - np.min(values)) if values.size else float("nan")


def case_table(payloads: list[dict[str, Any]]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for payload in payloads:
        events = list(payload.get("events", []))
        k_first = first_event_K(events)
        k25 = checkpoint_from_events(events, 25.0, strict=True)
        k50 = checkpoint_from_events(events, 50.0, strict=True)
        rows.append(
            {
                "candidate_id": str(payload["candidate_id"]),
                "temperature_K": int(round(float(payload["temperature_K"]))),
                "status": str(payload.get("status", "")),
                "seed": int(payload.get("seed", -1)),
                "reached_25um": checkpoint_reached(events, 25.0),
                "reached_50um": checkpoint_reached(events, 50.0),
                "K_first_MPa_sqrt_m": k_first,
                "K_25um_MPa_sqrt_m": k25,
                "K_50um_MPa_sqrt_m": k50,
                "R_rise_first_to_50_MPa_sqrt_m": k50 - k_first,
                "R_rise_25_to_50_MPa_sqrt_m": k50 - k25,
                "event_K_range_MPa_sqrt_m": event_range_K(events),
                "source_case_json": payload.get("_source_path", ""),
            }
        )
    return pd.DataFrame(rows).sort_values(["candidate_id", "temperature_K"])


def ordered_value(local: pd.DataFrame, temperature: int, column: str) -> float:
    row = local[local["temperature_K"] == int(temperature)]
    if len(row) != 1:
        return float("nan")
    return float(row.iloc[0][column])


def candidate_metrics(cases: pd.DataFrame, registry_ids: set[str]) -> pd.DataFrame:
    records: list[dict[str, Any]] = []
    expected = set(EXPECTED_TEMPERATURES_K)
    for candidate_id, local in cases.groupby("candidate_id", sort=True):
        if candidate_id not in registry_ids:
            raise RuntimeError(f"case payload candidate absent from registry: {candidate_id}")
        local = local.sort_values("temperature_K")
        temperatures = local["temperature_K"].astype(int).tolist()
        complete_grid = (
            len(local) == len(EXPECTED_TEMPERATURES_K)
            and set(temperatures) == expected
            and bool((local["status"].astype(str) == "complete").all())
            and bool(local[["reached_25um", "reached_50um"]].to_numpy(dtype=bool).all())
        )
        k50 = local["K_50um_MPa_sqrt_m"].to_numpy(dtype=float)
        finite_grid = bool(np.isfinite(k50).all())
        complete_grid = complete_grid and finite_grid

        low_mid = local[local["temperature_K"] <= 1100]
        high = local[local["temperature_K"] >= 1200]
        k50_low_mid = low_mid["K_50um_MPa_sqrt_m"].to_numpy(dtype=float)
        k50_high = high["K_50um_MPa_sqrt_m"].to_numpy(dtype=float)
        total_rise = local["R_rise_first_to_50_MPa_sqrt_m"].to_numpy(dtype=float)
        late_rise = local["R_rise_25_to_50_MPa_sqrt_m"].to_numpy(dtype=float)

        response_peak = peak_metrics(
            local["temperature_K"].to_numpy(dtype=float),
            k50,
        )
        k700 = ordered_value(local, 700, "K_50um_MPa_sqrt_m")
        k1100 = ordered_value(local, 1100, "K_50um_MPa_sqrt_m")
        k1200 = ordered_value(local, 1200, "K_50um_MPa_sqrt_m")
        k1300 = ordered_value(local, 1300, "K_50um_MPa_sqrt_m")
        k1400 = ordered_value(local, 1400, "K_50um_MPa_sqrt_m")
        high_rebound = max(0.0, k1400 - min(k1200, k1300))
        low_mid_median = finite_median(k50_low_mid.tolist())
        high_median = finite_median(k50_high.tolist())
        high_loss = low_mid_median - high_median
        pre_high_change = k1100 - k700
        terminal_high_change = k1400 - k1100

        median_total_rise = finite_median(total_rise.tolist())
        median_late_rise = finite_median(late_rise.tolist())
        median_abs_total_rise = finite_median(np.abs(total_rise).tolist())
        median_abs_late_rise = finite_median(np.abs(late_rise).tolist())
        positive_total_fraction = float(np.mean(total_rise > 0.25))

        weakT_gate = bool(
            complete_grid
            and float(np.ptp(k50)) <= 8.0
            and abs(high_loss) <= 4.0
            and 0.5 <= median_total_rise <= 6.0
            and -0.5 <= median_late_rise <= 3.0
            and positive_total_fraction >= 0.60
            and response_peak.peak_prominence < 4.0
        )
        ceramic_gate = bool(
            complete_grid
            and float(np.ptp(k50_low_mid)) <= 5.0
            and 1.0 <= high_loss <= 8.0
            and abs(pre_high_change) <= 4.0
            and -10.0 <= terminal_high_change <= -1.0
            and median_abs_total_rise <= 2.5
            and median_abs_late_rise <= 1.5
            and high_rebound <= 1.0
            and response_peak.peak_prominence < 3.0
        )

        weakT_score = (
            float(np.ptp(k50)) / 8.0
            + abs(high_loss) / 4.0
            + abs(median_total_rise - 2.5) / 3.0
            + abs(median_late_rise - 0.75) / 2.0
            + finite_std(total_rise.tolist()) / 4.0
            + response_peak.peak_prominence / 4.0
            + (0.0 if complete_grid else 100.0)
        )
        ceramic_score = (
            float(np.ptp(k50_low_mid)) / 5.0
            + abs(high_loss - 3.5) / 4.0
            + abs(pre_high_change) / 4.0
            + abs(median_total_rise) / 2.5
            + median_abs_late_rise / 1.5
            + high_rebound / 1.0
            + response_peak.peak_prominence / 3.0
            + (0.0 if complete_grid else 100.0)
        )

        record: dict[str, Any] = {
            "candidate_id": candidate_id,
            "complete_temperature_grid": complete_grid,
            "case_count": len(local),
            "K50_mean_MPa_sqrt_m": float(np.mean(k50)) if finite_grid else float("nan"),
            "K50_median_MPa_sqrt_m": finite_median(k50.tolist()),
            "K50_temperature_span_MPa_sqrt_m": float(np.ptp(k50)) if finite_grid else float("nan"),
            "K50_low_mid_span_700_1100_MPa_sqrt_m": float(np.ptp(k50_low_mid)),
            "K50_low_mid_median_MPa_sqrt_m": low_mid_median,
            "K50_high_median_1200_1400_MPa_sqrt_m": high_median,
            "high_temperature_toughness_loss_MPa_sqrt_m": high_loss,
            "pre_high_temperature_change_K1100_minus_K700_MPa_sqrt_m": pre_high_change,
            "terminal_change_K1400_minus_K1100_MPa_sqrt_m": terminal_high_change,
            "high_temperature_rebound_MPa_sqrt_m": high_rebound,
            "median_R_rise_first_to_50_MPa_sqrt_m": median_total_rise,
            "median_abs_R_rise_first_to_50_MPa_sqrt_m": median_abs_total_rise,
            "median_R_rise_25_to_50_MPa_sqrt_m": median_late_rise,
            "median_abs_R_rise_25_to_50_MPa_sqrt_m": median_abs_late_rise,
            "positive_R_rise_temperature_fraction": positive_total_fraction,
            "median_event_K_range_MPa_sqrt_m": finite_median(
                local["event_K_range_MPa_sqrt_m"].astype(float).tolist()
            ),
            "K50_peak_temperature_K": response_peak.peak_temperature_K,
            "K50_peak_prominence_MPa_sqrt_m": response_peak.peak_prominence,
            "weakT_gate": weakT_gate,
            "weakT_score": weakT_score,
            "ceramic_gate": ceramic_gate,
            "ceramic_score": ceramic_score,
        }
        for temperature in EXPECTED_TEMPERATURES_K:
            record[f"K50_T{temperature}K_MPa_sqrt_m"] = ordered_value(
                local, temperature, "K_50um_MPa_sqrt_m"
            )
            record[f"Rrise_T{temperature}K_MPa_sqrt_m"] = ordered_value(
                local, temperature, "R_rise_first_to_50_MPa_sqrt_m"
            )
        records.append(record)
    return pd.DataFrame(records)


def ranked(frame: pd.DataFrame, gate: str, score: str) -> pd.DataFrame:
    return frame.sort_values(
        [gate, score, "candidate_id"],
        ascending=[False, True, True],
        kind="stable",
    ).reset_index(drop=True)


def json_safe(value: Any) -> Any:
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return float(value) if math.isfinite(float(value)) else None
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    return value


def select_distinct(weak: pd.DataFrame, ceramic: pd.DataFrame) -> tuple[pd.Series, pd.Series]:
    weak_row = weak.iloc[0]
    ceramic_candidates = ceramic[ceramic["candidate_id"] != weak_row["candidate_id"]]
    if ceramic_candidates.empty:
        raise RuntimeError("unable to select distinct weak-T and ceramic candidates")
    return weak_row, ceramic_candidates.iloc[0]


def write_selected_registry(
    registry: pd.DataFrame,
    weak_row: pd.Series,
    ceramic_row: pd.Series,
    out_path: Path,
) -> None:
    selected_ids = [str(weak_row["candidate_id"]), str(ceramic_row["candidate_id"])]
    source = registry.set_index(registry["candidate_id"].astype(str), drop=False)
    missing = [candidate for candidate in selected_ids if candidate not in source.index]
    if missing:
        raise RuntimeError(f"selected candidate missing from registry: {missing}")
    rows = source.loc[selected_ids].copy().reset_index(drop=True)
    rows.insert(0, "paper_material_class", ["weakT_FCC_like", "ceramic_like"])
    rows.insert(1, "selection_role", ["primary", "primary"])
    rows.insert(2, "oneD_selection_score", [weak_row["weakT_score"], ceramic_row["ceramic_score"]])
    rows.insert(3, "oneD_strict_gate_passed", [weak_row["weakT_gate"], ceramic_row["ceramic_gate"]])
    rows.to_csv(out_path, index=False)


def main() -> int:
    args = parse_args()
    if args.top_count < 1:
        raise ValueError("top-count must be positive")
    case_root, registry_path = locate_inputs(args)
    registry = read_registry(registry_path)
    payloads = load_payloads(case_root)
    cases = case_table(payloads)
    metrics = candidate_metrics(cases, set(registry["candidate_id"].astype(str)))
    weak = ranked(metrics, "weakT_gate", "weakT_score")
    ceramic = ranked(metrics, "ceramic_gate", "ceramic_score")
    weak_row, ceramic_row = select_distinct(weak, ceramic)

    out = args.out_dir.expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)
    cases.to_csv(out / "case_temperature_rcurve_metrics.csv", index=False)
    weak.to_csv(out / "weakT_FCC_like_ranked.csv", index=False)
    ceramic.to_csv(out / "ceramic_like_ranked.csv", index=False)
    write_selected_registry(
        registry,
        weak_row,
        ceramic_row,
        out / "selected_weakT_ceramic_registry.csv",
    )

    summary = {
        "schema": "v9.13_weakT_ceramic_selection_from_384_1d_v1",
        "run_dir": str(args.run_dir.expanduser().resolve()),
        "case_root": str(case_root),
        "candidate_registry": str(registry_path),
        "candidate_count": int(len(metrics)),
        "expected_temperatures_K": list(EXPECTED_TEMPERATURES_K),
        "complete_candidate_count": int(metrics["complete_temperature_grid"].sum()),
        "weakT_strict_gate_count": int(metrics["weakT_gate"].sum()),
        "ceramic_strict_gate_count": int(metrics["ceramic_gate"].sum()),
        "selected": {
            "weakT_FCC_like": weak_row.to_dict(),
            "ceramic_like": ceramic_row.to_dict(),
        },
        "selection_policy": {
            "weakT_FCC_like": (
                "minimum weakT score, preferring strict gate: small K50 temperature span, "
                "small positive first-to-50um R-curve rise, and no pronounced peak"
            ),
            "ceramic_like": (
                "minimum ceramic score, preferring strict gate: flat 700-1100K response, "
                "weak loss concentrated at 1200-1400K, negligible R-curve rise, and no rebound"
            ),
        },
    }
    (out / "selection_summary.json").write_text(
        json.dumps(json_safe(summary), indent=2, sort_keys=True) + "\n"
    )

    columns = [
        "candidate_id",
        "complete_temperature_grid",
        "K50_temperature_span_MPa_sqrt_m",
        "high_temperature_toughness_loss_MPa_sqrt_m",
        "median_R_rise_first_to_50_MPa_sqrt_m",
        "median_R_rise_25_to_50_MPa_sqrt_m",
        "K50_peak_prominence_MPa_sqrt_m",
    ]
    print("V913_WEAKT_FCC_LIKE_SELECTION")
    print(weak.head(args.top_count)[columns + ["weakT_gate", "weakT_score"]].to_string(index=False))
    print("\nV913_CERAMIC_LIKE_SELECTION")
    print(ceramic.head(args.top_count)[columns + ["ceramic_gate", "ceramic_score"]].to_string(index=False))
    print("\nSELECTED")
    print(
        f"weakT_FCC_like={weak_row['candidate_id']} strict={bool(weak_row['weakT_gate'])} "
        f"score={float(weak_row['weakT_score']):.6g}"
    )
    print(
        f"ceramic_like={ceramic_row['candidate_id']} strict={bool(ceramic_row['ceramic_gate'])} "
        f"score={float(ceramic_row['ceramic_score']):.6g}"
    )
    print(f"registry={out / 'selected_weakT_ceramic_registry.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
