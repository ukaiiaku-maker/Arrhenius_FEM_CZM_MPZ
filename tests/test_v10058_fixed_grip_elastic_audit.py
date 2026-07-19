from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from arrhenius_fracture.fixed_grip_elastic_audit_v10058 import (
    FixedGripAuditConfig,
    fixed_grip_release_rate,
    run_audit,
    select_contiguous_plateau,
    three_point_derivative_at_center,
)


def test_three_point_derivative_is_exact_for_quadratic_on_unequal_spacing():
    x0 = 0.5
    xm = 0.47
    xp = 0.52

    def f(x):
        return 3.0 * x * x - 2.0 * x + 7.0

    derivative = three_point_derivative_at_center(xm, x0, xp, f(xm), f(x0), f(xp))
    assert derivative == pytest.approx(6.0 * x0 - 2.0, rel=1.0e-13, abs=1.0e-13)


def test_fixed_grip_release_rate_uses_negative_energy_derivative():
    out = fixed_grip_release_rate(
        a_minus_m=0.49,
        a_zero_m=0.50,
        a_plus_m=0.515,
        U_minus_J_per_m=10.2,
        U_zero_J_per_m=10.0,
        U_plus_J_per_m=9.7,
    )
    assert out["G_fixed_grip_J_per_m2"] > 0.0
    assert out["G_backward_J_per_m2"] > 0.0
    assert out["G_forward_J_per_m2"] > 0.0
    assert out["effective_increment_asymmetry"] > 0.0


def test_contiguous_plateau_prefers_widest_valid_window():
    rows = [
        {
            "outer_radius_m": radius,
            "J_active_elements": 20,
            "contour_within_safety_limit": True,
            "J_full_J_per_m2": value,
        }
        for radius, value in zip(
            [80e-6, 100e-6, 140e-6, 180e-6, 240e-6, 300e-6],
            [8.0, 9.7, 10.0, 10.1, 9.9, 12.0],
        )
    ]
    selected = select_contiguous_plateau(
        rows,
        value_key="J_full_J_per_m2",
        relative_tolerance=0.04,
        minimum_points=3,
        minimum_active_elements=12,
    )
    assert selected["status"] == "plateau_selected"
    assert selected["outer_radii_m"] == pytest.approx([100e-6, 140e-6, 180e-6, 240e-6])
    assert selected["maximum_relative_spread"] <= 0.04


def test_config_rejects_invalid_geometry():
    with pytest.raises(ValueError):
        FixedGripAuditConfig(crack_m=3.0e-3).validate()


def test_small_elastic_audit_writes_artifacts(tmp_path: Path):
    config = FixedGripAuditConfig(
        width_m=1.0e-3,
        height_m=2.0e-3,
        crack_m=0.25e-3,
        notch_half_thickness_m=40.0e-6,
        total_grip_opening_m=0.5e-6,
        nx=16,
        ny=32,
        tip_ratio=1.20,
        anisotropic=False,
        minimum_plateau_points=2,
        minimum_J_active_elements=1,
        mesh_relative_tolerance=10.0,
        crack_increment_relative_tolerance=10.0,
        contour_relative_tolerance=10.0,
        J_over_G_min=1.0e-12,
        J_over_G_max=1.0e12,
        energy_closure_relative_tolerance=1.0e-4,
    )
    payload = run_audit(
        config,
        tip_h_fine_m=(40.0e-6, 20.0e-6),
        crack_increment_m=(40.0e-6, 20.0e-6),
        contour_outer_m=(60.0e-6, 80.0e-6, 100.0e-6),
        out=tmp_path,
    )
    assert payload["schema"] == "fixed_grip_elastic_convergence_v10_0_5_8"
    assert np.isfinite(payload["convergence"]["fixed_grip_G_finest_J_per_m2"])
    assert (tmp_path / "fixed_grip_energy_release_v10_0_5_8.csv").exists()
    assert (tmp_path / "fixed_grip_J_contours_v10_0_5_8.csv").exists()
    summary = tmp_path / "fixed_grip_elastic_convergence_v10_0_5_8.json"
    assert summary.exists()
    saved = json.loads(summary.read_text())
    assert saved["status"] == payload["status"]
