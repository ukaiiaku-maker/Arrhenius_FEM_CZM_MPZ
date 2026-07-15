from __future__ import annotations

import pytest

from arrhenius_fracture.mode_i_first_passage_v9_11 import (
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
