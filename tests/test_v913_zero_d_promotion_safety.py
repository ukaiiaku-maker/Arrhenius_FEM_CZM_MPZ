from __future__ import annotations

import numpy as np
import pandas as pd

import scripts.run_v913_zero_d_large_search as base
from scripts.run_v913_zero_d_large_search_safe import (
    diverse_selection_safe,
    score_frame_safe,
)


def test_exact_gate_requires_completion() -> None:
    frame = pd.DataFrame(
        {
            "zeroD_two_sided_prominence_MPa_sqrt_m": [8.0, 8.0],
            "zeroD_post_peak_drop_MPa_sqrt_m": [9.0, 9.0],
            "zeroD_high_temperature_rebound_MPa_sqrt_m": [1.0, 1.0],
            "zeroD_peak_temperature_K": [1000.0, 1000.0],
            "zeroD_peak_value_MPa_sqrt_m": [60.0, 60.0],
            "zeroD_peak_internal": [1, 1],
            "zeroD_complete": [1, 0],
        }
    )
    scored = score_frame_safe(
        frame,
        "zeroD",
        minimum_prominence=5.0,
        minimum_drop=5.0,
        maximum_rebound=3.0,
        peak_min=850.0,
        peak_max=1100.0,
    )
    assert scored["zeroD_gate_pass"].tolist() == [True, False]
    assert scored.loc[1, "zeroD_objective"] >= 500.0


def test_diverse_selection_preserves_complete_strict_passes(monkeypatch) -> None:
    rows = []
    for index in range(6):
        row = {
            "candidate_id": f"candidate_{index}",
            "zeroD_objective": float(index),
            "zeroD_peak_temperature_K": 1000.0,
            "zeroD_peak_value_MPa_sqrt_m": 60.0 + index,
            "zeroD_two_sided_prominence_MPa_sqrt_m": 8.0,
            "zeroD_post_peak_drop_MPa_sqrt_m": 9.0,
            "zeroD_high_temperature_rebound_MPa_sqrt_m": 1.0,
            "zeroD_peak_internal": 1,
            "zeroD_peak_in_desired_window": 1,
            "zeroD_complete": 1,
            "zeroD_gate_pass": index in (0, 1, 2),
        }
        for field in base.ACTIVE_CANDIDATE_PARAMETER_FIELDS:
            row[field] = 1.0 + 0.01 * index
        rows.append(row)
    rows[2]["zeroD_complete"] = 0
    frame = pd.DataFrame(rows)

    monkeypatch.setattr(
        base,
        "_normalize_features",
        lambda candidate_frame, _policy: np.column_stack(
            (
                np.arange(len(candidate_frame), dtype=float),
                np.square(np.arange(len(candidate_frame), dtype=float)),
            )
        ),
    )
    promoted = diverse_selection_safe(frame, {}, 4)
    ids = set(promoted["candidate_id"])
    assert {"candidate_0", "candidate_1"}.issubset(ids)
    assert "candidate_2" not in ids
    assert promoted["zeroD_complete"].astype(bool).all()
    assert len(promoted) == 4
