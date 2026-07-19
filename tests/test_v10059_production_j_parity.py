from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from arrhenius_fracture.mode_i_first_passage_v10_0_5_9_production_j_probe import (
    _elastic_update_plasticity,
    validate_source_transform_v10059,
)
from arrhenius_fracture.production_j_parity_v10059 import (
    analyze_production_j_parity_v10059,
)
from run_v10_0_5_9_production_j_parity import (
    _log_tail,
    _probe_command,
    build_parser,
)


def _reference():
    return {
        "schema": "fixed_grip_elastic_convergence_v10_0_5_8",
        "passed": True,
        "convergence": {"fixed_grip_G_finest_J_per_m2": 100.0},
        "geometry_factors": {"sigma_gross_MPa": 200.0},
    }


def _probe(opening_um, ratio, no_exclusion_ratio=None):
    sigma_MPa = 100.0 * opening_um
    sigma_Pa = sigma_MPa * 1.0e6
    ref_metric = 100.0 / (200.0e6) ** 2
    metric = ratio * ref_metric
    no_metric = metric if no_exclusion_ratio is None else no_exclusion_ratio * ref_metric
    J = metric * sigma_Pa**2
    J_no = no_metric * sigma_Pa**2
    return {
        "Uapp_m": opening_um * 1.0e-6,
        "Uapp_um": opening_um,
        "sigma_gross_MPa": sigma_MPa,
        "elastic_energy_closure_relative_error": 1.0e-12,
        "mesh": {"production_refinement_radius_um": 100.0, "hbar_tip_m": 2.5e-6},
        "geometry": {"effective_killed_tip_m": 0.5e-3},
        "contours": [{
            "outer_radius_m": 240.0e-6,
            "outer_radius_um": 240.0,
            "production_path": "anisotropic_root_cluster_with_2killr_exclusion",
            "production_exclude_radius_um": 5.0,
            "J_full_J_per_m2": J,
            "J_tension_filtered_J_per_m2": 1.001 * J,
            "J_full_no_exclusion_J_per_m2": J_no,
            "J_full_over_sigma2_m_per_Pa": metric,
            "J_full_no_exclusion_over_sigma2_m_per_Pa": no_metric,
        }],
    }


def test_parity_analysis_passes_quadratic_response():
    result = analyze_production_j_parity_v10059(
        reference=_reference(), probes=[_probe(1.0, 1.02), _probe(2.0, 1.02)]
    )
    assert result["status"] == "production_J_parity_passed"
    assert result["passed"] is True
    assert result["production"]["median_over_reference"] == pytest.approx(1.02)


def test_parity_analysis_identifies_exclusion_control():
    result = analyze_production_j_parity_v10059(
        reference=_reference(),
        probes=[
            _probe(1.0, 0.65, no_exclusion_ratio=1.01),
            _probe(2.0, 0.65, no_exclusion_ratio=1.01),
        ],
    )
    assert result["status"] == "production_exclusion_disk_controls_J_mismatch"
    assert result["no_exclusion_ablation"]["parity_passed"] is True


def test_parity_analysis_identifies_mesh_support_mismatch():
    result = analyze_production_j_parity_v10059(
        reference=_reference(), probes=[_probe(1.0, 0.65), _probe(2.0, 0.65)]
    )
    assert result["status"] == "production_refinement_extent_or_mesh_support_mismatch"


def test_parity_analysis_rejects_nonquadratic_scaling():
    result = analyze_production_j_parity_v10059(
        reference=_reference(),
        probes=[_probe(1.0, 1.0), _probe(2.0, 1.08)],
        elastic_scaling_relative_tolerance=0.02,
    )
    assert result["status"] == "production_J_not_quadratic_in_load"


def test_elastic_update_plasticity_preserves_state():
    ep = np.arange(12, dtype=float).reshape(3, 4)
    rho = np.arange(4, dtype=float) + 5.0
    ep_out, rho_out, dot = _elastic_update_plasticity(
        ep, rho, np.zeros((3, 4)), object(), 700.0, 1.0, object(), object()
    )
    assert np.array_equal(ep_out, ep)
    assert np.array_equal(rho_out, rho)
    assert np.array_equal(dot, np.zeros_like(rho))
    assert ep_out is not ep
    assert rho_out is not rho


def test_source_transform_compiles_and_preserves_stack():
    result = validate_source_transform_v10059()
    assert result["source_transform_preflight_passed"] is True
    assert result["production_recorder"] is True
    assert result["root_front_production_exclusion"] is True
    assert result["straight_progressive_no_exclusion"] is True
    assert result["no_unconditional_kill_r_read"] is True
    assert result["full_audited_v10055_stack"] is True
    assert result["constitutive_physics_changed"] is False


def test_probe_command_satisfies_vhcf_cycle_mode_contract(tmp_path: Path):
    args = build_parser().parse_args(["--reference-json", str(tmp_path / "reference.json")])
    command = _probe_command(args, tmp_path / "case", 1.0e-6)
    mode_index = command.index("--cycle-block-mode")
    assert command[mode_index + 1] == "hazard_limited"
    assert "requested_cap" not in command
    assert command[command.index("--cycles-max") + 1] == "1"
    assert command[command.index("--max-block-cycles") + 1] == "1"


def test_probe_command_uses_one_authoritative_mpz_length(tmp_path: Path):
    args = build_parser().parse_args([
        "--reference-json",
        str(tmp_path / "reference.json"),
        "--L-pz-um",
        "20",
        "--mpz-length-um",
        "100",
    ])
    command = _probe_command(args, tmp_path / "case", 1.0e-6)
    legacy_m = float(command[command.index("--L-pz") + 1])
    modern_um = float(command[command.index("--mpz-length-um") + 1])
    assert legacy_m == pytest.approx(100.0e-6)
    assert modern_um == pytest.approx(100.0)
    assert legacy_m * 1.0e6 == pytest.approx(modern_um)


def test_log_tail_reports_only_requested_lines(tmp_path: Path):
    path = tmp_path / "probe.log"
    path.write_text("\n".join(f"line {index}" for index in range(10)))
    assert _log_tail(path, max_lines=3) == "line 7\nline 8\nline 9"
