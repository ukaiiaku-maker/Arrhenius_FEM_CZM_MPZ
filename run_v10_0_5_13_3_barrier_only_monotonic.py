#!/usr/bin/env python3
"""v10.0.5.13.3 runner: tip/source-only MPZ with rate-preserving macro-steps."""
from __future__ import annotations

import math
from pathlib import Path
import shlex
import subprocess

import run_v10_0_5_13_2_barrier_only_monotonic as _base

POINT_RELEASE = "10.0.5.13.3"
ENTRY_MODULE = (
    "arrhenius_fracture.mode_i_first_passage_v10_0_5_13_3_barrier_only"
)
TIP_ONLY_POLICY = {
    "policy_id": "preserve_existing_tip_only_moving_mpz_v1005133",
    "bulk_plasticity_mode": "tip_only",
    "mpz_length_um": 100.0,
    "mpz_n_bins": 80,
    "candidate_source_inventory_applied": False,
    "candidate_source_refresh_applied": False,
    "candidate_encounter_recovery_applied": False,
    "candidate_shielding_blunting_applied": False,
    "candidate_initial_state_applied": False,
    "state_configuration_source": "existing_tip_only_2d_solver_and_explicit_cli",
    "state_evolution_source": "existing_moving_crack_tip_MPZ",
    "continuum_bulk_role": "elastic_fem_only",
    "uniform_bulk_mobile_retained_state_active": False,
}
_ORIGINAL_BUILD = _base._build_command
_ORIGINAL_PREFLIGHT = _base._preflight
_ULTIMATE = _base._base._base
_ORIGINAL_COMPLETE = _ULTIMATE._case_is_complete
_ORIGINAL_SUMMARIZE = _ULTIMATE._summarize


def _option_after(cmd: list[str], name: str) -> str | None:
    try:
        index = cmd.index(name)
    except ValueError:
        return None
    return cmd[index + 1] if index + 1 < len(cmd) else None


def _build_command(py, args, option_key, T_K, target_um, case_dir):
    cmd = _ORIGINAL_BUILD(py, args, option_key, T_K, target_um, case_dir)
    old = "arrhenius_fracture.mode_i_first_passage_v10_0_5_13_2_barrier_only"
    try:
        cmd[cmd.index(old)] = ENTRY_MODULE
    except ValueError as exc:
        raise RuntimeError(f"v10.0.5.13.2 command lacks expected entry {old}") from exc
    try:
        mode_index = cmd.index("--bulk-plasticity-mode") + 1
    except ValueError as exc:
        raise RuntimeError("barrier campaign command lacks --bulk-plasticity-mode") from exc
    cmd[mode_index] = "tip_only"
    if "bulk_same_pt_km" in cmd:
        raise RuntimeError("v10.0.5.13.3 command retained bulk_same_pt_km")
    return cmd


def _tip_only_case_is_complete(case_dir: Path, option_key: str, target_um: float) -> bool:
    status = _ULTIMATE._read_json(case_dir / _ULTIMATE.STATUS_FILE)
    production = _ULTIMATE._read_json(case_dir / _ULTIMATE.PRODUCTION_MANIFEST)
    integration = _ULTIMATE._read_json(case_dir / "mpz_v9_11_integration_audit.json")
    complete, extension_um = _ULTIMATE.completion_status(case_dir, target_um)
    bulk_pt = integration.get("bulk_PT", {})
    scope = production.get("plasticity_scope", {})
    return bool(
        status.get("status") == "complete"
        and status.get("option_key") == option_key
        and math.isclose(
            float(status.get("target_extension_um", -1.0)),
            target_um,
            abs_tol=1.0e-9,
        )
        and complete
        and extension_um is not None
        and production.get("run_completed_without_exception") is True
        and production.get("candidate_state_fields_applied") is False
        and production.get("mesh_refinement_runtime", {}).get("actual_radius_verified") is True
        and production.get("barrier_option", {}).get("option_key") == option_key
        and scope.get("bulk_plasticity_mode") == "tip_only"
        and scope.get("moving_crack_tip_mpz_active") is True
        and scope.get("uniform_bulk_mobile_retained_state_active") is False
        and bulk_pt.get("mode") == "tip_only"
        and bulk_pt.get("explicit_mobile_retained_state") is False
        and bulk_pt.get("source_interpretation") == "moving_crack_tip_MPZ_only"
    )


def _tip_only_summarize(case_dir, option_key, T_K, target_um, returncode, reused):
    row = _ORIGINAL_SUMMARIZE(
        case_dir, option_key, T_K, target_um, returncode, reused
    )
    integration = _ULTIMATE._read_json(case_dir / "mpz_v9_11_integration_audit.json")
    production = _ULTIMATE._read_json(case_dir / _ULTIMATE.PRODUCTION_MANIFEST)
    bulk_pt = integration.get("bulk_PT", {})
    scope = production.get("plasticity_scope", {})
    row.update(
        {
            "point_release": POINT_RELEASE,
            "bulk_plasticity_mode": bulk_pt.get("mode"),
            "tip_only_verified": bool(
                bulk_pt.get("mode") == "tip_only"
                and bulk_pt.get("explicit_mobile_retained_state") is False
                and bulk_pt.get("source_interpretation")
                == "moving_crack_tip_MPZ_only"
                and scope.get("bulk_plasticity_mode") == "tip_only"
                and scope.get("moving_crack_tip_mpz_active") is True
            ),
            "moving_crack_tip_mpz_active": scope.get(
                "moving_crack_tip_mpz_active"
            ),
            "uniform_bulk_mobile_retained_state_active": scope.get(
                "uniform_bulk_mobile_retained_state_active"
            ),
        }
    )
    _ULTIMATE._write_json(case_dir / _ULTIMATE.STATUS_FILE, row)
    return row


def _preflight(py: str, run_tests: bool):
    _ORIGINAL_PREFLIGHT(py, False)
    compile_cmd = [
        py,
        "-m",
        "py_compile",
        "arrhenius_fracture/mode_i_first_passage_v10_0_5_13_3_barrier_only.py",
        "run_v10_0_5_13_3_barrier_only_monotonic.py",
    ]
    cp = subprocess.run(compile_cmd, text=True)
    if cp.returncode != 0:
        raise SystemExit(
            f"v10.0.5.13.3 compile failed: {shlex.join(compile_cmd)}"
        )
    if run_tests:
        test_cmd = [
            py,
            "-m",
            "pytest",
            "-q",
            "tests/test_v100513_barrier_only.py",
            "tests/test_v1005131_preserved_state.py",
            "tests/test_v1005132_startup_resolution_warning.py",
            "tests/test_v1005133_tip_only_ramp.py",
            "tests/test_v1005123_phase_c_repairs.py",
        ]
        cp = subprocess.run(test_cmd, text=True)
        if cp.returncode != 0:
            raise SystemExit(
                f"v10.0.5.13.3 tests failed: {shlex.join(test_cmd)}"
            )


def main():
    saved_build = _base._build_command
    saved_preflight = _base._preflight
    saved_release = _base.POINT_RELEASE
    saved_complete = _ULTIMATE._case_is_complete
    saved_summarize = _ULTIMATE._summarize
    saved_policy = _ULTIMATE.TWO_D_STATE_POLICY
    _base._build_command = _build_command
    _base._preflight = _preflight
    _base.POINT_RELEASE = POINT_RELEASE
    _ULTIMATE._case_is_complete = _tip_only_case_is_complete
    _ULTIMATE._summarize = _tip_only_summarize
    # Scope the tip-only campaign metadata and MPZ resolution to this release.
    # The shared v10.0.5.13 registry remains unchanged for legacy contracts.
    _ULTIMATE.TWO_D_STATE_POLICY = dict(TIP_ONLY_POLICY)
    try:
        return _base.main()
    finally:
        _base._build_command = saved_build
        _base._preflight = saved_preflight
        _base.POINT_RELEASE = saved_release
        _ULTIMATE._case_is_complete = saved_complete
        _ULTIMATE._summarize = saved_summarize
        _ULTIMATE.TWO_D_STATE_POLICY = saved_policy


if __name__ == "__main__":
    main()


__all__ = [
    "POINT_RELEASE",
    "ENTRY_MODULE",
    "TIP_ONLY_POLICY",
    "_build_command",
    "_tip_only_case_is_complete",
    "_tip_only_summarize",
    "_preflight",
    "main",
]
