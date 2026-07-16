#!/usr/bin/env python3
"""Run the v9.18.1 persistent-wake gate with active-event renewal rollback."""
from __future__ import annotations

import run_mpz_v9_18_persistent_plastic_wake as _v918

_original_build = _v918._build_command_v918


def _build_command_v9181(args, class_name, run_root, force_rerun):
    cmd = _original_build(args, class_name, run_root, force_rerun)
    old = "run_mpz_v9_18_mode_i_rcurve.py"
    new = "run_mpz_v9_18_1_mode_i_rcurve.py"
    try:
        cmd[cmd.index(old)] = new
    except ValueError as exc:
        raise RuntimeError(f"v9.18 driver token not found in command: {cmd}") from exc
    return cmd


def main():
    original = _v918._build_command_v918
    _v918._build_command_v918 = _build_command_v9181
    try:
        return _v918.main()
    finally:
        _v918._build_command_v918 = original


if __name__ == "__main__":
    main()
