from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from arrhenius_fracture.emergent_gnd_contract_v913 import ACTIVE_CANDIDATE_PARAMETER_FIELDS
from arrhenius_fracture.weakT_ceramic_objectives_v913 import (
    ceramic_score,
    diverse_select,
    response_metrics,
    weakT_score,
)


TEMPERATURES = [
    300,
    400,
    500,
    600,
    700,
    800,
    900,
    950,
    1000,
    1050,
    1100,
    1150,
    1200,
    1250,
    1300,
]

COARSE_TEMPERATURES = [300, 600, 900, 1100, 1200]


def test_weakT_gate_requires_flat_initiation_and_small_positive_rcurve():
    initial = [20.0, 20.2, 20.1, 20.3, 20.1, 20.2, 20.0, 20.1, 20.2, 20.1, 20.0, 19.9, 19.8, 19.7, 19.6]
    developed = [24.0, 24.2, 24.0, 24.3, 24.1, 24.2, 24.0, 24.2, 24.1, 24.0, 23.9, 23.8, 23.7, 23.6, 23.5]
    metrics = response_metrics(TEMPERATURES, initial, developed)
    gate, score = weakT_score(metrics)
    assert gate is True
    assert np.isfinite(score)
    assert metrics["low_temperature_initial_rise_MPa_sqrt_m"] < 1.0
    assert 3.0 < metrics["median_R_rise_MPa_sqrt_m"] < 5.0


def test_ceramic_gate_requires_high_temperature_loss_and_negligible_rcurve():
    initial = [25.0, 25.1, 25.0, 25.1, 25.0, 24.9, 25.0, 24.9, 24.8, 24.8, 24.7, 24.3, 23.5, 22.8, 22.2]
    developed = [25.4, 25.5, 25.4, 25.5, 25.4, 25.3, 25.4, 25.3, 25.2, 25.2, 25.1, 24.7, 23.9, 23.2, 22.6]
    metrics = response_metrics(TEMPERATURES, initial, developed)
    gate, score = ceramic_score(metrics)
    assert gate is True
    assert np.isfinite(score)
    assert metrics["initial_high_temperature_loss_MPa_sqrt_m"] > 0.5
    assert metrics["median_abs_R_rise_MPa_sqrt_m"] < 1.0


def test_five_temperature_grid_supports_both_class_gates():
    weak_initial = [20.0, 20.2, 20.1, 20.0, 19.8]
    weak_developed = [24.0, 24.2, 24.1, 24.0, 23.8]
    weak_metrics = response_metrics(
        COARSE_TEMPERATURES,
        weak_initial,
        weak_developed,
    )
    weak_gate, _ = weakT_score(weak_metrics)
    assert weak_gate is True

    ceramic_initial = [25.0, 25.1, 25.0, 24.8, 22.5]
    ceramic_developed = [25.4, 25.5, 25.4, 25.2, 22.9]
    ceramic_metrics = response_metrics(
        COARSE_TEMPERATURES,
        ceramic_initial,
        ceramic_developed,
    )
    ceramic_gate, _ = ceramic_score(ceramic_metrics)
    assert ceramic_gate is True


def test_failed_transfers_are_rejected_by_low_temperature_and_rcurve_terms():
    initial = [1.0, 2.0, 3.0, 5.0, 8.0, 10.0, 11.0, 11.5, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0]
    developed = [1.0, 4.0, 8.0, 15.0, 22.0, 28.0, 32.0, 35.0, 31.0, 31.0, 30.0, 28.0, 29.0, 30.0, 24.0]
    metrics = response_metrics(TEMPERATURES, initial, developed)
    weak_gate, _ = weakT_score(metrics)
    ceramic_gate, _ = ceramic_score(metrics)
    assert weak_gate is False
    assert ceramic_gate is False
    assert metrics["low_temperature_initial_rise_MPa_sqrt_m"] >= 7.0
    assert metrics["max_abs_R_rise_MPa_sqrt_m"] >= 20.0


def test_diverse_selection_preserves_count_and_unique_candidates():
    policy_path = Path(__file__).resolve().parents[1] / "mpz_v9_13_zero_d_weakT_ceramic_search_policy.json"
    policy = json.loads(policy_path.read_text())
    rows = []
    for index in range(12):
        row = {
            "candidate_id": f"candidate_{index:02d}",
            "gate": index < 6,
            "score": float(index),
        }
        for field in ACTIVE_CANDIDATE_PARAMETER_FIELDS:
            if field in policy["search_dimensions"]:
                spec = policy["search_dimensions"][field]
                fraction = (index + 1) / 13.0
                if spec["mode"] == "log10_delta":
                    value = 10 ** (
                        np.log10(float(spec["low"]))
                        + fraction
                        * (np.log10(float(spec["high"])) - np.log10(float(spec["low"])))
                    )
                else:
                    value = float(spec["low"]) + fraction * (
                        float(spec["high"]) - float(spec["low"])
                    )
                row[field] = value
        rows.append(row)
    frame = pd.DataFrame(rows)
    selected = diverse_select(
        frame,
        policy,
        count=5,
        score_column="score",
        gate_column="gate",
    )
    assert len(selected) == 5
    assert selected["candidate_id"].is_unique
    assert selected.iloc[0]["candidate_id"] == "candidate_00"
