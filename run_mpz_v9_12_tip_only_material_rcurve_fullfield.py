#!/usr/bin/env python3
"""Launch the v9.12 campaign through the v9.12 Mode-I field renderer."""
from __future__ import annotations

import run_mpz_v9_12_tip_only_material_rcurve as _base


_original_build_command = _base.build_command


def _build_command_fullfield(args, class_name, run_root, force_rerun):
    cmd = _original_build_command(args, class_name, run_root, force_rerun)
    old = "run_mpz_v9_11_mode_i_rcurve_3T.py"
    new = "run_mpz_v9_12_mode_i_rcurve.py"
    try:
        cmd[cmd.index(old)] = new
    except ValueError as exc:
        raise RuntimeError(f"v9.11 R-curve driver token not found in command: {cmd}") from exc
    return cmd


def main():
    original = _base.build_command
    _base.build_command = _build_command_fullfield
    try:
        return _base.main()
    finally:
        _base.build_command = original


if __name__ == "__main__":
    main()
