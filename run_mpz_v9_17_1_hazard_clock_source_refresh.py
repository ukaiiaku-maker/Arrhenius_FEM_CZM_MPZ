#!/usr/bin/env python3
"""Run the v9.17.1 one-fire routing-fixed hazard-clock material gate."""
from __future__ import annotations

import run_mpz_v9_17_hazard_clock_source_refresh as _v917

_original_build = _v917._build_command_v917


def _build_command_v9171(args, class_name, run_root, force_rerun):
    cmd = _original_build(args, class_name, run_root, force_rerun)
    old = "run_mpz_v9_17_mode_i_rcurve.py"
    new = "run_mpz_v9_17_1_mode_i_rcurve.py"
    try:
        cmd[cmd.index(old)] = new
    except ValueError as exc:
        raise RuntimeError(f"v9.17 driver token not found in command: {cmd}") from exc
    return cmd


def main():
    original = _v917._build_command_v917
    _v917._build_command_v917 = _build_command_v9171
    try:
        return _v917.main()
    finally:
        _v917._build_command_v917 = original


if __name__ == "__main__":
    main()
