#!/usr/bin/env python3
"""Phase-C v10.0.5.12.3 live runner.

Repairs two post-run integration defects without changing constitutive or FEM
physics:

* invokes the v10.0.5.12.3 entry that preserves refinement audit metadata
  through CZM topology rebuilds;
* supplies the removed ``numpy.trapz`` name from ``numpy.trapezoid`` when
  required by the inherited campaign summary routine.
"""
from __future__ import annotations

import shlex
import subprocess

import numpy as np

import run_v10_0_5_12_phase_c_monotonic as _base
import run_v10_0_5_12_1_phase_c_monotonic as _radius
import run_v10_0_5_12_2_phase_c_monotonic as _live

POINT_RELEASE = "10.0.5.12.3"
ENTRY_MODULE = "arrhenius_fracture.mode_i_first_passage_v10_0_5_12_3_phase_c"
_ORIGINAL_RADIUS_BUILD = _radius.build_command
_ORIGINAL_PREFLIGHT = _base.preflight


def ensure_numpy_trapz_compat() -> None:
    """Provide the legacy analysis name only when NumPy removed it."""
    if hasattr(np, "trapz"):
        return
    trapezoid = getattr(np, "trapezoid", None)
    if trapezoid is None:
        raise RuntimeError("NumPy provides neither trapezoid nor trapz integration")
    np.trapz = trapezoid


def build_command(py, args, option_key, T_K, target_um, case_dir):
    """Build the validated command and select the v10.0.5.12.3 entry."""
    cmd = _ORIGINAL_RADIUS_BUILD(py, args, option_key, T_K, target_um, case_dir)
    old_module = "arrhenius_fracture.mode_i_first_passage_v10_0_5_12_phase_c"
    try:
        index = cmd.index(old_module)
    except ValueError as exc:
        raise RuntimeError(f"Phase-C command lacks expected entry module {old_module}") from exc
    cmd[index] = ENTRY_MODULE
    return cmd


def preflight(py, run_tests):
    """Compile the inherited stack and run both original and repair tests."""
    _ORIGINAL_PREFLIGHT(py, False)
    compile_command = [
        py,
        "-m",
        "py_compile",
        "arrhenius_fracture/mode_i_first_passage_v10_0_5_12_3_phase_c.py",
        "run_v10_0_5_12_2_phase_c_monotonic.py",
        "run_v10_0_5_12_3_phase_c_monotonic.py",
    ]
    cp = subprocess.run(compile_command, text=True)
    if cp.returncode != 0:
        raise SystemExit(f"Phase-C repair compile failed: {shlex.join(compile_command)}")
    if run_tests:
        test_command = [
            py,
            "-m",
            "pytest",
            "-q",
            "tests/test_v100512_phase_c.py",
            "tests/test_v1005123_phase_c_repairs.py",
        ]
        cp = subprocess.run(test_command, text=True)
        if cp.returncode != 0:
            raise SystemExit(f"Phase-C repair tests failed: {shlex.join(test_command)}")


def main() -> None:
    ensure_numpy_trapz_compat()
    saved_live_release = _live.POINT_RELEASE
    saved_radius_build = _radius.build_command
    saved_preflight = _base.preflight
    _live.POINT_RELEASE = POINT_RELEASE
    _radius.build_command = build_command
    _base.preflight = preflight
    try:
        _live.main()
    finally:
        _live.POINT_RELEASE = saved_live_release
        _radius.build_command = saved_radius_build
        _base.preflight = saved_preflight


if __name__ == "__main__":
    main()


__all__ = [
    "POINT_RELEASE",
    "ENTRY_MODULE",
    "ensure_numpy_trapz_compat",
    "build_command",
    "preflight",
    "main",
    "_base",
]
