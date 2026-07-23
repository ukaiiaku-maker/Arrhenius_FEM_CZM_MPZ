from __future__ import annotations

import pandas as pd

from scripts.augment_mpz_v9_12_directional_peak_targets import (
    add_directional_peak_classifications,
    add_trajectory_metrics,
)
from scripts.run_v913_autonomous_dbtt_search import _select_rows


def test_existing_peak_objective_applies_to_autonomous_K_trajectory():
    table = pd.DataFrame(
        [
            {
                "candidate_id": "peak",
                "y__K25_1d_T700K": 20.0,
                "y__K25_1d_T800K": 30.0,
                "y__K25_1d_T900K": 50.0,
                "y__K25_1d_T1000K": 35.0,
                "y__K25_1d_T1100K": 25.0,
                "y__K25_1d_T1200K": 20.0,
            }
        ]
    )
    add_trajectory_metrics(
        table,
        prefix="y__K25_1d_",
        out_prefix="y__",
        low_max_K=700.0,
        high_min_K=1000.0,
        peak_min_K=800.0,
        peak_max_K=1000.0,
    )
    add_directional_peak_classifications(
        table,
        out_prefix="y__",
        peak_min_K=800.0,
        direction_threshold=5.0,
        peak_threshold=1.0,
    )

    row = table.iloc[0]
    assert row["y__peak_temperature_K"] == 900.0
    assert row["y__peak_rise"] == 30.0
    assert row["y__peak_drop"] == 30.0
    assert row["y__peak_prominence"] == 30.0
    assert bool(row["y__peak_like_1d"]) is True
    assert bool(row["y__direction_correct_1d"]) is False


def test_nested_sobol_prefix_is_balanced_across_selected_peak_parents():
    rows = []
    for family, parent in (
        ("peak", "parent_a"),
        ("peak", "parent_b"),
        ("plateau", "parent_c"),
    ):
        for index in range(4):
            rows.append(
                {
                    "candidate_id": f"{parent}_{index:04d}",
                    "campaign_parent_id": parent,
                    "campaign_parent_family": family,
                }
            )

    selected = _select_rows(
        rows,
        families=("peak",),
        candidate_ids=(),
        per_parent=2,
        parent_offset=1,
    )
    assert [row["candidate_id"] for row in selected] == [
        "parent_a_0001",
        "parent_a_0002",
        "parent_b_0001",
        "parent_b_0002",
    ]


def test_explicit_candidate_selection_preserves_requested_order():
    rows = [
        {"candidate_id": "a", "campaign_parent_family": "peak"},
        {"candidate_id": "b", "campaign_parent_family": "plateau"},
    ]
    selected = _select_rows(
        rows,
        families=("peak",),
        candidate_ids=("b", "a"),
        per_parent=1,
        parent_offset=0,
    )
    assert [row["candidate_id"] for row in selected] == ["b", "a"]
