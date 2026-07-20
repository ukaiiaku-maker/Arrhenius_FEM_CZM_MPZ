#!/usr/bin/env python3
"""Corrected Phase-C campaign wrapper for cluster-J radius semantics.

``sharp_front --rJ-cluster`` accepts the legacy domain length ``ell`` whose
actual outer radius is approximately ``8*ell``.  The original v10.0.5.12
campaign passed the requested physical 240 um outer radius directly as ell,
which produced a 1.92 mm contour.  This wrapper converts the requested physical
outer radius to the required legacy length while preserving the remainder of
the v10.0.5.12 campaign implementation.
"""
from __future__ import annotations

import math

import run_v10_0_5_12_phase_c_monotonic as _base

POINT_RELEASE = "10.0.5.12.1"
CLUSTER_J_RADIUS_TO_LEGACY_ELL = 8.0
_ORIGINAL_BUILD_COMMAND = _base.build_command


def cluster_j_legacy_length_m(outer_radius_um: float) -> float:
    """Convert requested physical cluster-J outer radius to legacy ell [m]."""
    radius_um = float(outer_radius_um)
    if not math.isfinite(radius_um) or radius_um <= 0.0:
        raise ValueError("cluster-J physical outer radius must be finite and positive")
    return radius_um * 1.0e-6 / CLUSTER_J_RADIUS_TO_LEGACY_ELL


def build_command(py, args, option_key, T_K, target_um, case_dir):
    """Build the original command and correct only ``--rJ-cluster``."""
    cmd = _ORIGINAL_BUILD_COMMAND(py, args, option_key, T_K, target_um, case_dir)
    try:
        index = cmd.index("--rJ-cluster")
    except ValueError as exc:
        raise RuntimeError("Phase-C command lacks --rJ-cluster") from exc
    if index + 1 >= len(cmd):
        raise RuntimeError("Phase-C command has no value after --rJ-cluster")

    ell_m = cluster_j_legacy_length_m(args.cluster_J_outer_um)
    cmd[index + 1] = _base.fstr(ell_m)
    return cmd


def main() -> None:
    saved = _base.build_command
    _base.build_command = build_command
    try:
        _base.main()
    finally:
        _base.build_command = saved


if __name__ == "__main__":
    main()


__all__ = [
    "POINT_RELEASE",
    "CLUSTER_J_RADIUS_TO_LEGACY_ELL",
    "cluster_j_legacy_length_m",
    "build_command",
    "main",
]
