import csv
from pathlib import Path

import pytest

import run_v10_0_5_6_kj_audit_bracket_audited as bracket_audited
import run_v10_0_5_6_stochastic_delta_sigma_audited as campaign_audited


def test_audited_reader_normalizes_boolean_strings(tmp_path: Path):
    path = tmp_path / "rows.csv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "first_passage_observed",
                "right_censored",
                "reached_cycle_horizon",
                "reached_target_extension",
            ],
        )
        writer.writeheader()
        writer.writerow(
            {
                "first_passage_observed": "False",
                "right_censored": "0",
                "reached_cycle_horizon": "True",
                "reached_target_extension": "no",
            }
        )
    row = bracket_audited._read_csv_with_booleans(path)[0]
    assert row["first_passage_observed"] is False
    assert row["right_censored"] is False
    assert row["reached_cycle_horizon"] is True
    assert row["reached_target_extension"] is False


def test_bracket_rejects_numerically_censored_case():
    rows = [
        {
            "delta_sigma_requested_MPa": 100.0,
            "first_passage_observed": False,
            "right_censored": True,
        }
    ]
    with pytest.raises(RuntimeError, match="numerically censored"):
        bracket_audited._classify_without_numerical_censoring(rows)


def test_deterministic_calibration_gets_explicit_text_label():
    rows = [
        {
            "step": 1,
            "cycle_limiter_code": 8,
            "cycles_block": 1.0e-9,
        }
    ]
    out = campaign_audited._enrich_with_calibration_fallback(rows, [])
    assert out[0]["stochastic_scheduler_mode"] == "deterministic_calibration"
    assert out[0]["cycle_limiter_label"] == "cycle_horizon"
    assert out[0]["stochastic_expected_state_events"] == 0.0
