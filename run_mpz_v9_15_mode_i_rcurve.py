#!/usr/bin/env python3
"""Select the v9.15 coupled cohesive-MPZ Mode-I entry point."""
from __future__ import annotations

import run_mpz_v9_11_mode_i_rcurve_3T as _base

_original_build_command = _base.build_command


def _build_command_v915(py, args, class_name, manifest, T_K, case_dir):
    cmd = _original_build_command(py, args, class_name, manifest, T_K, case_dir)
    old = "arrhenius_fracture.mode_i_first_passage_v9_11"
    new = "arrhenius_fracture.mode_i_first_passage_v9_15"
    try:
        cmd[cmd.index(old)] = new
    except ValueError as exc:
        raise RuntimeError(f"v9.11 Mode-I module token not found in command: {cmd}") from exc
    return cmd


def main():
    original = _base.build_command
    _base.build_command = _build_command_v915
    try:
        return _base.main()
    finally:
        _base.build_command = original


if __name__ == "__main__":
    main()
