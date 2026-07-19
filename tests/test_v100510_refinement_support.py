from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from arrhenius_fracture.config import GeometryConfig, MeshConfig
from arrhenius_fracture.mode_i_first_passage_v10_0_5_10_refinement_probe import (
    validate_source_transform_v100510,
)
from arrhenius_fracture.physical_refinement_mesh_v100510 import (
    clear_physical_refinement_v100510,
    configure_physical_refinement_v100510,
    make_physical_refinement_mesh_v100510,
)
from arrhenius_fracture.production_j_refinement_support_v100510 import (
    analyze_refinement_support_v100510,
    deduplicate_contours_v100510,
    radial_mesh_support_v100510,
)
from run_v10_0_5_10_refinement_support_audit import (
    ENTRY_MODULE,
    _probe_command,
    build_parser,
)


def test_fixed_physical_refinement_mesh_records_actual_radius():
    geom = GeometryConfig(Lx=2e-3, Ly=4e-3, a0=0.5e-3)
    cfg = MeshConfig(tip_h_fine=20e-6, tip_ratio=1.15)
    configure_physical_refinement_v100510(200e-6)
    try:
        mesh = make_physical_refinement_mesh_v100510(geom, cfg, seed=42)
    finally:
        clear_physical_refinement_v100510()
    assert mesh.production_refinement_radius_m == pytest.approx(200e-6)
    assert mesh.production_refinement_policy == "fixed_physical_radius_same_radial_ring_law"
    assert np.min(np.linalg.norm(mesh.nodes - np.array([[geom.a0, 0.0]]), axis=1)) < 1e-12
    rows = radial_mesh_support_v100510(mesh, np.array([geom.a0, 0.0]), [0, 100e-6, 200e-6])
    assert len(rows) == 2
    assert all(row["element_count"] > 0 for row in rows)
    assert all(row["inside_configured_refinement"] for row in rows)


def test_contour_deduplication_prefers_requested_row():
    rows = [
        {"outer_radius_m": 240e-6, "value": 1, "is_requested_production_contour": False},
        {"outer_radius_m": 0.00023999999999999998, "value": 2, "is_requested_production_contour": True},
        {"outer_radius_m": 300e-6, "value": 3, "is_requested_production_contour": False},
    ]
    out = deduplicate_contours_v100510(rows)
    assert len(out) == 2
    assert out[0]["value"] == 2
    assert out[0]["is_requested_production_contour"] is True


def _reference():
    return {
        "schema": "fixed_grip_elastic_convergence_v10_0_5_8",
        "passed": True,
        "convergence": {"fixed_grip_G_finest_J_per_m2": 100.0},
        "geometry_factors": {"sigma_gross_MPa": 200.0},
    }


def _probe(opening_um: float, ratios, radius_um: float = 330.0):
    sigma_MPa = 100.0 * opening_um
    sigma_Pa = sigma_MPa * 1e6
    ref_metric = 100.0 / (200e6) ** 2
    contours = []
    for outer_um, ratio in zip((180.0, 240.0, 300.0), ratios):
        metric = ratio * ref_metric
        J = metric * sigma_Pa**2
        contours.append(
            {
                "outer_radius_m": outer_um * 1e-6,
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
    return {
        "Uapp_m": opening_um * 1e-6,
        "Uapp_um": opening_um,
        "sigma_gross_MPa": sigma_MPa,
        "elastic_energy_closure_relative_error": 1e-12,
        "mesh": {
            "production_refinement_radius_m": radius_um * 1e-6,
            "production_refinement_radius_um": radius_um,
            "hbar_tip_m": 2.5e-6,
        },
        "geometry": {"effective_killed_tip_m": 0.5e-3},
        "contours": contours,
    }


def test_refinement_analysis_passes_parity_and_plateau():
    result = analyze_refinement_support_v100510(
        reference=_reference(),
        probes=[_probe(1.0, [1.02, 1.01, 1.00]), _probe(2.0, [1.02, 1.01, 1.00])],
    )
    assert result["status"] == "production_refinement_support_parity_passed"
    assert result["passed"] is True
    assert result["physical_refinement_support_passed"] is True
    assert result["contour_stability_passed"] is True


def test_refinement_analysis_rejects_contour_instability_before_parity_claim():
    result = analyze_refinement_support_v100510(
        reference=_reference(),
        probes=[_probe(1.0, [1.20, 1.00, 0.80]), _probe(2.0, [1.20, 1.00, 0.80])],
    )
    assert result["status"] == "production_J_contour_instability"
    assert result["passed"] is False


def test_refinement_analysis_rejects_insufficient_physical_radius():
    result = analyze_refinement_support_v100510(
        reference=_reference(),
        probes=[
            _probe(1.0, [1.0, 1.0, 1.0], radius_um=250.0),
            _probe(2.0, [1.0, 1.0, 1.0], radius_um=250.0),
        ],
    )
    assert result["status"] == "production_refinement_support_inadequate"


def test_v100510_source_transform_compiles_and_replaces_only_recorder():
    audit = validate_source_transform_v100510()
    assert audit["source_transform_preflight_passed"] is True
    assert audit["refinement_recorder"] is True
    assert audit["v10_0_5_9_production_path_preserved"] is True
    assert audit["full_audited_v10055_stack"] is True
    assert audit["constitutive_physics_changed"] is False


def test_v100510_runner_uses_new_probe_module(tmp_path: Path):
    args = build_parser().parse_args(["--reference-json", str(tmp_path / "reference.json")])
    command = _probe_command(args, tmp_path / "case", 1e-6)
    assert command[command.index("-m") + 1] == ENTRY_MODULE
    assert command[command.index("--cycle-block-mode") + 1] == "hazard_limited"
    assert args.tip_refinement_radius_um == pytest.approx(330.0)
    assert args.contour_outer_um == "100 140 180 240 300"
