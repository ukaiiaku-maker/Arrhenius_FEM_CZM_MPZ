import csv
from pathlib import Path
from types import SimpleNamespace

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


def test_single_front_missing_active_count_uses_radial_resolution(monkeypatch):
    args = SimpleNamespace(
        tip_h_fine_m=2.5e-6,
        minimum_J_active_elements=12,
        plateau_relative_tolerance=0.10,
        plateau_minimum_points=3,
    )
    monkeypatch.setattr(bracket_audited, "_single_front_table_exists", lambda row: False)
    rows = []
    for outer_um, slope, ratio in [
        (140.0, 0.023465, 0.394449),
        (180.0, 0.023865, 0.401169),
        (240.0, 0.023418, 0.393668),
    ]:
        rows.append(
            {
                "outer_radius_m": outer_um * 1e-6,
                "outer_radius_um": outer_um,
                "contour_within_safety_limit": True,
                "J_active_elements": 0,
                "KJ_per_sigma_gross_sqrt_m": slope,
                "KJ_over_K_LEFM_gross": ratio,
                "case_root": f"unused_{outer_um}",
            }
        )
    repaired = bracket_audited._repair_resolution_rows(args, rows)
    assert all(row["J_active_elements"] is None for row in repaired)
    assert all(row["J_active_elements_available"] is False for row in repaired)
    assert all(row["J_resolution_passed"] is True for row in repaired)
    assert repaired[0]["J_annulus_radial_elements_estimate"] == pytest.approx(42.0)

    selected = bracket_audited._select_with_resolution_fallback(args, repaired)
    assert selected["status"] == "plateau_selected"
    assert selected["selected_outer_radius_m"] == pytest.approx(180e-6)
    assert selected["plateau_max_relative_spread"] < 0.02
    assert selected["selected_row"]["J_active_elements"] is None


def test_true_zero_active_count_still_fails_resolution(monkeypatch):
    args = SimpleNamespace(
        tip_h_fine_m=2.5e-6,
        minimum_J_active_elements=12,
        plateau_relative_tolerance=0.10,
        plateau_minimum_points=3,
    )
    monkeypatch.setattr(bracket_audited, "_single_front_table_exists", lambda row: True)
    rows = [
        {
            "outer_radius_m": outer_um * 1e-6,
            "contour_within_safety_limit": True,
            "J_active_elements": 0,
            "KJ_per_sigma_gross_sqrt_m": 0.0235,
            "KJ_over_K_LEFM_gross": 0.40,
            "case_root": f"unused_{outer_um}",
        }
        for outer_um in (140.0, 180.0, 240.0)
    ]
    repaired = bracket_audited._repair_resolution_rows(args, rows)
    assert all(row["J_active_elements_available"] is True for row in repaired)
    assert all(row["J_resolution_passed"] is False for row in repaired)
    selected = bracket_audited._select_with_resolution_fallback(args, repaired)
    assert selected["status"] == "no_valid_plateau"
