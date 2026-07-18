from __future__ import annotations

import json

import numpy as np
import pandas as pd

import analytical_downselect_mpz_v9_10_4_7 as analytical
import evaluate_dynamic_1d_mpz_v9_10_4_7 as dynamic


def test_fixed_zero_screen_removes_direct_cleavage_temperature_slopes() -> None:
    bounds = analytical.search_bounds("fixed_zero")
    assert bounds["cleave_gT_eV_per_K"] == (0.0, 0.0)
    assert bounds["cleave_sT_GPa_per_K"] == (0.0, 0.0)


def test_sobol_vectors_have_exact_requested_count_and_fixed_dimensions() -> None:
    bounds = analytical.search_bounds("fixed_zero")
    vectors = analytical.sobol_parameter_vectors(37, bounds, seed=123)
    assert vectors.shape == (37, len(analytical.PARAMETER_NAMES))
    g_index = analytical.PARAMETER_NAMES.index("cleave_gT_eV_per_K")
    s_index = analytical.PARAMETER_NAMES.index("cleave_sT_GPa_per_K")
    assert np.all(vectors[:, g_index] == 0.0)
    assert np.all(vectors[:, s_index] == 0.0)


def test_bracket_selection_preserves_independent_transition_windows() -> None:
    rows = []
    for bracket, low, high in (("T0500_0600K", 500.0, 600.0), ("T0700_0800K", 700.0, 800.0)):
        for rank, objective in enumerate((3.0, 1.0, 2.0)):
            rows.append(
                {
                    "candidate_id": f"{bracket}_{rank}",
                    "analysis_valid": True,
                    "screen_pass": True,
                    "transition_bracket": bracket,
                    "coarse_transition_low_T_K": low,
                    "coarse_transition_high_T_K": high,
                    "objective": objective,
                }
            )
    selected = analytical.select_by_bracket(pd.DataFrame(rows), per_bracket_keep=2)
    assert selected.groupby("transition_bracket").size().to_dict() == {
        "T0500_0600K": 2,
        "T0700_0800K": 2,
    }
    assert set(selected[selected.transition_bracket == "T0500_0600K"].objective) == {1.0, 2.0}


def test_dynamic_stage_scores_exactly_four_transition_temperatures() -> None:
    temperatures = np.asarray([600.0, 633.3333333333, 666.6666666667, 700.0])
    full = np.asarray([10.0, 11.0, 19.0, 21.0])
    off = np.asarray([10.0, 10.1, 10.3, 10.5])
    metrics = dynamic._transition_metrics(temperatures, full, off)
    assert metrics["moving_1d_valid"]
    assert metrics["moving_1d_accept"]
    assert metrics["moving_1d_edge_ratio"] == 2.1
    assert metrics["moving_1d_transition_width_K"] <= 100.0


def test_manifest_schedule_requires_four_points() -> None:
    row = pd.Series(
        {
            "candidate_id": "candidate",
            "refinement_transition_temperatures_K": json.dumps([600.0, 633.3333333333, 666.6666666667, 700.0]),
        }
    )
    assert dynamic._schedule_from_row(row) == [600.0, 633.3333333333, 666.6666666667, 700.0]
