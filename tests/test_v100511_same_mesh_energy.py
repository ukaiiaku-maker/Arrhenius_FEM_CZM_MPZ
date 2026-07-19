from __future__ import annotations

import pytest

from arrhenius_fracture.mode_i_first_passage_v10_0_5_11_same_mesh_probe import (
    validate_source_transform_v100511,
)
from arrhenius_fracture.production_j_same_mesh_energy_v100511 import (
    analyze_same_mesh_energy_v100511,
    second_order_forward_release_rate_v100511,
)


def _reference():
    return {
        "schema": "fixed_grip_elastic_convergence_v10_0_5_8",
        "passed": True,
        "convergence": {"fixed_grip_G_finest_J_per_m2": 100.0},
        "geometry_factors": {"sigma_gross_MPa": 200.0},
    }


def _probe(opening_um: float, *, same_mesh_G_scale: float = 1.0):
    sigma_MPa = 100.0 * opening_um
    sigma_Pa = sigma_MPa * 1.0e6
    reference_metric = 100.0 / (200.0e6) ** 2
    ratios = (0.88, 0.89, 0.90)
    contours = []
    for outer_um, ratio in zip((180.0, 240.0, 300.0), ratios):
        metric = ratio * reference_metric
        J = metric * sigma_Pa**2
        contours.append(
            {
                "outer_radius_m": outer_um * 1.0e-6,
                "outer_radius_um": outer_um,
                "production_path": "straight_progressive_cluster_no_exclusion",
                "production_exclude_radius_um": 0.0,
                "J_full_J_per_m2": J,
                "J_tension_filtered_J_per_m2": J,
                "J_full_no_exclusion_J_per_m2": J,
                "J_full_over_sigma2_m_per_Pa": metric,
                "J_full_no_exclusion_over_sigma2_m_per_Pa": metric,
            }
        )
    median_J = contours[1]["J_full_J_per_m2"]
    G = same_mesh_G_scale * median_J
    return {
        "Uapp_m": opening_um * 1.0e-6,
        "Uapp_um": opening_um,
        "sigma_gross_MPa": sigma_MPa,
        "elastic_energy_closure_relative_error": 1.0e-12,
        "mesh": {
            "production_refinement_radius_m": 330.0e-6,
            "production_refinement_radius_um": 330.0,
            "hbar_tip_m": 2.5e-6,
        },
        "geometry": {"effective_killed_tip_m": 0.5e-3},
        "contours": contours,
        "same_mesh_fixed_grip_energy_release": {
            "G_second_order_forward_J_per_m2": G,
            "G_over_sigma2_m_per_Pa": G / sigma_Pa**2,
            "one_step_two_step_relative_difference": 0.02,
            "crack_stations_um": [500.0, 522.0, 544.0],
        },
    }


def test_second_order_forward_release_is_exact_for_quadratic_energy():
    # U(a) = 10 - 4a - a^2, so G(0) = -dU/da = 4.
    result = second_order_forward_release_rate_v100511(0.0, 1.0, 2.0, 10.0, 5.0, -2.0)
    assert result["G_second_order_forward_J_per_m2"] == pytest.approx(4.0)
    assert result["station_spacing_um"] == pytest.approx(1.0e6)


def test_same_mesh_analysis_can_pass_when_external_cross_mesh_reference_is_low():
    result = analyze_same_mesh_energy_v100511(
        reference=_reference(),
        probes=[_probe(1.0), _probe(2.0)],
    )
    assert result["base_v10_0_5_10_analysis"]["status"] == "production_J_parity_failed_with_adequate_support"
    assert result["status"] == "production_J_same_mesh_energy_parity_passed"
    assert result["passed"] is True


def test_same_mesh_analysis_rejects_internal_J_G_mismatch():
    result = analyze_same_mesh_energy_v100511(
        reference=_reference(),
        probes=[_probe(1.0, same_mesh_G_scale=1.25), _probe(2.0, same_mesh_G_scale=1.25)],
    )
    assert result["status"] == "production_J_same_mesh_energy_mismatch"
    assert result["passed"] is False


def test_v100511_source_transform_compiles_and_supplies_same_mesh_state():
    audit = validate_source_transform_v100511()
    assert audit["source_transform_preflight_passed"] is True
    assert audit["same_mesh_recorder"] is True
    assert audit["boundary_data_supplied"] is True
    assert audit["fixed_grip_opening_supplied"] is True
    assert audit["v10_0_5_9_production_path_preserved"] is True
    assert audit["full_audited_v10055_stack"] is True
    assert audit["constitutive_physics_changed"] is False
