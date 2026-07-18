import math

import pytest

from arrhenius_fracture.kj_audit_v10056 import (
    SpecimenGeometryV10056,
    build_kj_audit_row,
    classify_first_passage_rows,
    contour_geometry_audit,
    edge_crack_tension_geometry_factor,
    enrich_stochastic_block_rows,
    select_contour_plateau,
)


def test_edge_crack_reference_for_active_geometry():
    geometry = SpecimenGeometryV10056()
    assert geometry.a_over_W == pytest.approx(0.25)
    assert edge_crack_tension_geometry_factor(0.25) == pytest.approx(
        1.5009609375
    )
    row = build_kj_audit_row(
        Ftop_N_per_thickness=388.8888889e6 * geometry.width_m,
        KJ_Pa_sqrt_m=23.1342201e6,
        outer_radius_m=250.0e-6,
        geometry=geometry,
        n_active_elements=20,
    )
    assert row["sigma_gross_MPa"] == pytest.approx(388.8888889)
    assert row["K_LEFM_gross_MPa_sqrt_m"] == pytest.approx(23.1342201)
    assert row["KJ_over_K_LEFM_gross"] == pytest.approx(1.0)


def test_legacy_eight_mm_cluster_contour_is_rejected():
    geometry = SpecimenGeometryV10056()
    audit = contour_geometry_audit(
        outer_radius_m=8.0e-3,
        geometry=geometry,
        safety_fraction=0.8,
    )
    assert geometry.nearest_tip_boundary_m == pytest.approx(0.5e-3)
    assert audit["outer_over_nearest_boundary"] == pytest.approx(16.0)
    assert audit["contour_closes_inside_body"] is False
    assert audit["contour_within_safety_limit"] is False


def test_plateau_selection_uses_only_safe_consecutive_radii():
    rows = []
    for outer_um, slope in [
        (60, 0.050),
        (100, 0.061),
        (140, 0.062),
        (180, 0.0615),
        (240, 0.0618),
        (450, 0.090),
    ]:
        rows.append(
            {
                "outer_radius_m": outer_um * 1e-6,
                "contour_within_safety_limit": outer_um <= 400,
                "J_active_elements": 30,
                "KJ_per_sigma_gross_sqrt_m": slope,
                "KJ_over_K_LEFM_gross": slope / 0.0617,
            }
        )
    selected = select_contour_plateau(
        rows, relative_tolerance=0.03, minimum_points=3
    )
    assert selected["status"] == "plateau_selected"
    assert selected["selected_outer_radius_m"] in pytest.approx(
        [140e-6, 180e-6]
    )
    assert 450e-6 not in selected["plateau_outer_radii_m"]


def test_scheduler_text_labels_replace_ambiguous_numeric_code():
    rows = [
        {"step": 1, "cycle_limiter_code": 0},
        {"step": 2, "cycle_limiter_code": 8},
    ]
    records = [
        {
            "mode": "tau_leap",
            "event_rate_per_cycle": 0.2,
            "expected_state_events": 3.0,
            "limiter": "stochastic_tau_leap",
        },
        {
            "mode": "quiet",
            "event_rate_per_cycle": 0.0,
            "expected_state_events": 0.0,
            "limiter": "cycle_horizon",
        },
    ]
    out = enrich_stochastic_block_rows(rows, records)
    assert out[0]["cycle_limiter_label"] == "stochastic_tau_leap"
    assert out[0]["stochastic_scheduler_mode"] == "tau_leap"
    assert out[1]["cycle_limiter_label"] == "cycle_horizon"


def test_first_passage_bracket_classification():
    rows = [
        {"delta_sigma_requested_MPa": 100.0, "first_passage_observed": False},
        {"delta_sigma_requested_MPa": 150.0, "first_passage_observed": False},
        {"delta_sigma_requested_MPa": 200.0, "first_passage_observed": True},
        {"delta_sigma_requested_MPa": 250.0, "first_passage_observed": True},
    ]
    bracket = classify_first_passage_rows(rows)
    assert bracket["status"] == "bracketed"
    assert bracket["stress_interval_MPa"] == [150.0, 200.0]
