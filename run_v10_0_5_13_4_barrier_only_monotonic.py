#!/usr/bin/env python3
"""v10.0.5.13.4 campaign wrapper for tip-only policy propagation."""
from __future__ import annotations

import shlex
import subprocess

import run_v10_0_5_13_3_barrier_only_monotonic as _base

POINT_RELEASE = "10.0.5.13.4"
ENTRY_MODULE = (
    "arrhenius_fracture.mode_i_first_passage_v10_0_5_13_4_barrier_only"
)
_ORIGINAL_BUILD = _base._build_command
_ORIGINAL_PREFLIGHT = _base._preflight


def _build_command(py, args, option_key, T_K, target_um, case_dir):
    cmd = _ORIGINAL_BUILD(py, args, option_key, T_K, target_um, case_dir)
    old = "arrhenius_fracture.mode_i_first_passage_v10_0_5_13_3_barrier_only"
    try:
        cmd[cmd.index(old)] = ENTRY_MODULE
    except ValueError as exc:
        raise RuntimeError(f"v10.0.5.13.3 command lacks expected entry {old}") from exc
    mode_index = cmd.index("--bulk-plasticity-mode") + 1
    if cmd[mode_index] != "tip_only":
        raise RuntimeError(
            "v10.0.5.13.4 command did not retain tip_only plasticity scope"
        )
    return cmd


def _preflight(py: str, run_tests: bool):
    _ORIGINAL_PREFLIGHT(py, False)
    compile_cmd = [
        py,
        "-m",
        "py_compile",
        "arrhenius_fracture/mode_i_first_passage_v10_0_5_13_4_barrier_only.py",
        "run_v10_0_5_13_4_barrier_only_monotonic.py",
    ]
    cp = subprocess.run(compile_cmd, text=True)
    if cp.returncode != 0:
        raise SystemExit(
            f"v10.0.5.13.4 compile failed: {shlex.join(compile_cmd)}"
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
            "tests/test_v1005134_tip_only_policy_propagation.py",
            "tests/test_v1005123_phase_c_repairs.py",
        ]
        cp = subprocess.run(test_cmd, text=True)
        if cp.returncode != 0:
            raise SystemExit(
                f"v10.0.5.13.4 tests failed: {shlex.join(test_cmd)}"
            )


def main():
    saved_build = _base._build_command
    saved_preflight = _base._preflight
    saved_release = _base.POINT_RELEASE
    _base._build_command = _build_command
    _base._preflight = _preflight
    _base.POINT_RELEASE = POINT_RELEASE
    try:
        return _base.main()
    finally:
        _base._build_command = saved_build
        _base._preflight = saved_preflight
        _base.POINT_RELEASE = saved_release


if __name__ == "__main__":
    main()


__all__ = ["POINT_RELEASE", "ENTRY_MODULE", "_build_command", "_preflight", "main"]
