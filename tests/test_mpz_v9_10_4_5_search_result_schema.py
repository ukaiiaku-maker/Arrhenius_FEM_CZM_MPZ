from __future__ import annotations

import numpy as np

import optimize_mpz_v9_10_4_5_narrow_dbtt as patched


def test_incomplete_detailed_result_supplies_empty_event_detail() -> None:
    x = np.zeros(len(patched._base.PARAMETER_NAMES), dtype=float)
    result = patched.stabilize_detailed_result(
        {
            "objective": 1.1e6,
            "completion_loss": 1.0e5,
            "parameters": {"source_sites_per_system": 1.0},
            "temperature_detail": [{"T_K": 300.0, "completed": False}],
        },
        x,
        details=True,
    )
    assert result["event_detail"] == []
    assert result["temperature_detail"][0]["T_K"] == 300.0
    assert result["evaluation_status"] == "INCOMPLETE_CANDIDATE"


def test_early_rejection_has_stable_detailed_schema() -> None:
    x = np.zeros(len(patched._base.PARAMETER_NAMES), dtype=float)
    result = patched.stabilize_detailed_result(
        {
            "objective": 1.0e8,
            "parameters": {},
        },
        x,
        details=True,
    )
    assert result["parameters"] == {}
    assert result["temperature_detail"] == []
    assert result["event_detail"] == []
    assert result["evaluation_status"] == "EARLY_REJECTED_CANDIDATE"


def test_non_detailed_result_is_not_expanded() -> None:
    x = np.zeros(len(patched._base.PARAMETER_NAMES), dtype=float)
    original = {"objective": 42.0}
    result = patched.stabilize_detailed_result(original, x, details=False)
    assert result == original
    assert "event_detail" not in result
