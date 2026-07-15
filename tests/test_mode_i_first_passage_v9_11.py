from __future__ import annotations

import csv

import pytest

from arrhenius_fracture.mode_i_first_passage_v9_11 import (
    _clock_summary,
    _derived_mpz_length_args,
    unit_mode_i_directional_factors,
    validate_direct_mode_args,
)


def test_mode_i_directional_factors_are_exactly_unity():
    out = unit_mode_i_directional_factors(0.2, 0.8, 0.3, 0.1)
    assert out["cleavage_factor"] == 1.0
    assert out["emission_factor"] == 1.0
    assert out["directional_factor_cap_active"] is False
    assert out["mode_I_direct_material_calibration"] is True


def test_mode_i_entry_point_rejects_mixed_mode_overrides():
    with pytest.raises(SystemExit):
        validate_direct_mode_args(["--mixity-shear-coeff", "0.2"])
    with pytest.raises(SystemExit):
        validate_direct_mode_args(["--reference-cleavage-shape=0.7"])


def test_mode_i_entry_point_allows_full_solver_controls():
    validate_direct_mode_args([
        "--mode", "2d",
        "--crystal-aniso",
        "--crystal-compete",
        "--crystal-theta-deg", "45",
        "--crack-backend", "adaptive_czm",
    ])


def test_mpz_length_is_propagated_to_legacy_resolution_fields():
    derived, length_m = _derived_mpz_length_args(["--mpz-length-um", "100"])
    assert length_m == pytest.approx(100.0e-6)
    assert derived == [
        "--L-pz", "0.0001",
        "--mpz-length-m", "0.0001",
    ]


def test_conflicting_legacy_process_zone_length_is_rejected():
    with pytest.raises(SystemExit):
        _derived_mpz_length_args([
            "--mpz-length-um", "100",
            "--L-pz", "1e-6",
        ])


def test_clock_summary_reads_final_and_first_fire_residual(tmp_path):
    path = tmp_path / "steps_0700K.csv"
    with path.open("w", newline="") as fp:
        writer = csv.DictWriter(fp, fieldnames=["step", "B", "n_fire"])
        writer.writeheader()
        writer.writerows([
            {"step": 1, "B": 0.2, "n_fire": 0},
            {"step": 2, "B": 0.06, "n_fire": 1},
            {"step": 3, "B": 0.4, "n_fire": 0},
        ])
    out = _clock_summary(tmp_path, 700.0)
    assert out["B_final"] == pytest.approx(0.4)
    assert out["B_first_fire_residual"] == pytest.approx(0.06)
    assert out["B_final_step"] == pytest.approx(3.0)
    assert out["B_first_fire_step"] == pytest.approx(2.0)
    assert out["B_summary_source"] == "steps_0700K.csv"
