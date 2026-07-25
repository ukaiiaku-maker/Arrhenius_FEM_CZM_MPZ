from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "select_v913_weakT_ceramic_from_1d.py"
)
SPEC = importlib.util.spec_from_file_location("select_v913_weakT_ceramic", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def add_candidate(rows, candidate_id, k50_values, total_rise, late_rise):
    for temperature, k50 in zip(MODULE.EXPECTED_TEMPERATURES_K, k50_values):
        rows.append(
            {
                "candidate_id": candidate_id,
                "temperature_K": temperature,
                "status": "complete",
                "seed": 3621,
                "reached_25um": True,
                "reached_50um": True,
                "K_first_MPa_sqrt_m": k50 - total_rise,
                "K_25um_MPa_sqrt_m": k50 - late_rise,
                "K_50um_MPa_sqrt_m": k50,
                "R_rise_first_to_50_MPa_sqrt_m": total_rise,
                "R_rise_25_to_50_MPa_sqrt_m": late_rise,
                "event_K_range_MPa_sqrt_m": abs(total_rise),
                "source_case_json": "synthetic",
            }
        )


def synthetic_cases():
    rows = []
    add_candidate(
        rows,
        "weak_candidate",
        [40.0, 40.5, 41.0, 41.0, 40.8, 40.5, 40.2, 40.0, 39.8, 39.5],
        2.0,
        1.0,
    )
    add_candidate(
        rows,
        "ceramic_candidate",
        [32.0, 32.0, 32.0, 32.0, 32.0, 32.0, 32.0, 30.5, 29.5, 28.5],
        0.0,
        0.0,
    )
    add_candidate(
        rows,
        "peak_candidate",
        [25.0, 30.0, 40.0, 50.0, 60.0, 55.0, 45.0, 35.0, 30.0, 28.0],
        8.0,
        4.0,
    )
    return pd.DataFrame(rows)


def test_selection_distinguishes_weakT_and_ceramic_shapes():
    cases = synthetic_cases()
    metrics = MODULE.candidate_metrics(
        cases,
        {"weak_candidate", "ceramic_candidate", "peak_candidate"},
    )
    by_id = metrics.set_index("candidate_id")
    assert bool(by_id.loc["weak_candidate", "weakT_gate"]) is True
    assert bool(by_id.loc["ceramic_candidate", "ceramic_gate"]) is True
    assert bool(by_id.loc["peak_candidate", "weakT_gate"]) is False
    assert bool(by_id.loc["peak_candidate", "ceramic_gate"]) is False

    weak = MODULE.ranked(metrics, "weakT_gate", "weakT_score")
    ceramic = MODULE.ranked(metrics, "ceramic_gate", "ceramic_score")
    weak_row, ceramic_row = MODULE.select_distinct(weak, ceramic)
    assert weak_row["candidate_id"] == "weak_candidate"
    assert ceramic_row["candidate_id"] == "ceramic_candidate"


def test_monotonic_ceramic_curve_has_zero_peak_prominence():
    metrics = MODULE.candidate_metrics(
        synthetic_cases(),
        {"weak_candidate", "ceramic_candidate", "peak_candidate"},
    )
    ceramic = metrics.set_index("candidate_id").loc["ceramic_candidate"]
    assert ceramic["K50_peak_prominence_MPa_sqrt_m"] == 0.0
    assert ceramic["high_temperature_toughness_loss_MPa_sqrt_m"] > 0.0
