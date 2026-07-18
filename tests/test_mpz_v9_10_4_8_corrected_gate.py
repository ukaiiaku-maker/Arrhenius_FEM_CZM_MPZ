from __future__ import annotations

import numpy as np
import pandas as pd

from evaluate_dynamic_1d_mpz_v9_10_4_8 import corrected_transition_metrics
from rebuild_analytical_promotion_mpz_v9_10_4_8 import rebuild


def test_sharp_final_subinterval_transition_is_not_rejected_as_missing_high_shelf():
    T = np.asarray([700.0, 733.3333333333, 766.6666666667, 800.0])
    full = np.asarray([8.4293442952, 8.3111411957, 9.3962013085, 19.2702586616])
    off = np.asarray([8.4281806481, 8.2636570959, 8.1027496765, 7.9450221666])
    result = corrected_transition_metrics(T, full, off)
    assert result["moving_1d_accept"] is True
    assert result["moving_1d_edge_ratio"] > 2.0
    assert result["moving_1d_midpoint_separation_ratio"] < 1.5
    assert result["moving_1d_transition_width_K"] < 100.0


def test_factor_two_candidate_below_physical_low_toughness_floor_is_rejected():
    T = np.asarray([500.0, 533.3333333333, 566.6666666667, 600.0])
    full = np.asarray([2.6476, 2.6436, 3.8022, 9.1905])
    off = np.asarray([2.6445, 2.5906, 2.5390, 2.4894])
    result = corrected_transition_metrics(T, full, off)
    assert result["moving_1d_accept"] is False
    assert result["moving_1d_reason"] == "moving_1d_low_toughness_below_floor"


def test_rebuild_fills_bracket_after_screen_passes():
    rows = []
    for i, passed in enumerate([True, True, False, False, False]):
        rows.append(
            {
                "candidate_id": f"c{i}",
                "analysis_valid": True,
                "screen_pass": passed,
                "objective": float(i),
                "transition_bracket": "T0700_0800K",
                "coarse_transition_low_T_K": 700.0,
            }
        )
    selected = rebuild(pd.DataFrame(rows), per_bracket_keep=4)
    assert selected.candidate_id.tolist() == ["c0", "c1", "c2", "c3"]
    assert selected.selection_basis.tolist() == [
        "analytical_screen_pass",
        "analytical_screen_pass",
        "best_valid_fill",
        "best_valid_fill",
    ]
