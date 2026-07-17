"""Regression tests for v9.10.4.2 invalid-transition audit handling."""
from __future__ import annotations

import math

from audit_mpz_v9_10_4_2_current_dbtt import (
    guarded_best_adjacent_transition,
    transition_summary_fields,
)


def test_guard_adds_nan_schema_for_nonfinite_toughness() -> None:
    result = guarded_best_adjacent_transition(
        [300.0, 400.0, 500.0, 600.0],
        [15.0, float("nan"), 30.0, 31.0],
    )
    assert result["valid"] is False
    assert result["reason"] == "nonfinite_toughness"
    assert math.isnan(float(result["shelf_ratio"]))
    assert math.isnan(float(result["jump_concentration"]))
    assert math.isnan(float(result["transition_low_K"]))


def test_valid_transition_is_not_modified() -> None:
    result = guarded_best_adjacent_transition(
        [300.0, 400.0, 500.0, 600.0],
        [15.0, 15.2, 31.0, 31.5],
    )
    assert result["valid"] is True
    assert result["shelf_ratio"] > 2.0
    assert result["transition_low_K"] == 400.0
    assert result["transition_high_K"] == 500.0


def test_summary_fields_preserve_invalid_reason() -> None:
    result = {
        "valid": False,
        "loss": 1.0e12,
        "reason": "nonfinite_toughness",
        "shelf_ratio": float("nan"),
        "jump_concentration": float("nan"),
        "transition_low_K": float("nan"),
        "transition_high_K": float("nan"),
        "plasticity_off_ratio": float("nan"),
    }
    out = transition_summary_fields("new", result)
    assert out["new_transition_valid"] is False
    assert out["new_transition_reason"] == "nonfinite_toughness"
    assert out["new_transition_loss"] == 1.0e12
    assert math.isnan(out["new_shelf_ratio"])
