from __future__ import annotations

import pandas as pd
import pytest

from arrhenius_fracture.emergent_gnd_contract_v913 import (
    ACTIVE_CANDIDATE_PARAMETER_FIELDS,
)
from scripts.prepare_v913_zero_d_to_1d_screen import select_rows


def candidate(index: int, tier: str, diversity: int, objective: float) -> dict:
    row = {
        "candidate_id": f"candidate_{index}",
        "promotion_tier": tier,
        "zeroD_complete": 1,
        "zeroD_objective": objective,
        "zeroD_rank": index + 1,
        "diversity_rank": diversity,
    }
    for field in ACTIVE_CANDIDATE_PARAMETER_FIELDS:
        row[field] = 1.0 + 0.01 * index
    row["Tref_K"] = 481.33
    row["peierls_nu0_s"] = 1.0e12
    row["taylor_nu0_s"] = 1.0e11
    return row


def test_selection_preserves_all_strict_and_fills_by_diversity() -> None:
    source = pd.DataFrame(
        [
            candidate(0, "strict_gate", 8, 2.0),
            candidate(1, "strict_gate", 7, 1.0),
            candidate(2, "relaxed_desired_peak", 4, 0.2),
            candidate(3, "relaxed_desired_peak", 1, 5.0),
            candidate(4, "relaxed_desired_peak", 2, 4.0),
            candidate(5, "relaxed_desired_peak", 3, 3.0),
        ]
    )
    selected = select_rows(source, 4)
    assert selected["candidate_id"].tolist() == [
        "candidate_1",
        "candidate_0",
        "candidate_3",
        "candidate_4",
    ]
    assert selected["oneD_selection_tier"].tolist() == [
        "strict_zeroD",
        "strict_zeroD",
        "relaxed_diverse_zeroD",
        "relaxed_diverse_zeroD",
    ]


def test_selection_rejects_incomplete_strict_row() -> None:
    rows = [
        candidate(0, "strict_gate", 1, 0.0),
        candidate(1, "relaxed_desired_peak", 2, 1.0),
        candidate(2, "relaxed_desired_peak", 3, 2.0),
    ]
    rows[0]["zeroD_complete"] = 0
    selected = select_rows(pd.DataFrame(rows), 2)
    assert "candidate_0" not in set(selected["candidate_id"])
    assert len(selected) == 2


def test_corrected_registry_schema_without_completion_column_is_accepted() -> None:
    rows = [
        candidate(0, "strict_gate", 3, 0.0),
        candidate(1, "relaxed_desired_peak", 1, 1.0),
        candidate(2, "relaxed_desired_peak", 2, 2.0),
    ]
    for row in rows:
        row.pop("zeroD_complete")
    selected = select_rows(pd.DataFrame(rows), 2)
    assert selected["candidate_id"].tolist() == ["candidate_0", "candidate_1"]
    assert selected["oneD_selection_tier"].tolist() == [
        "strict_zeroD",
        "relaxed_diverse_zeroD",
    ]


def test_registry_without_completion_column_rejects_unknown_tier() -> None:
    rows = [candidate(0, "finite_complete_boundary", 1, 0.0)]
    rows[0].pop("zeroD_complete")
    with pytest.raises(RuntimeError, match="non-promotable tiers"):
        select_rows(pd.DataFrame(rows), 1)
