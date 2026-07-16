#!/usr/bin/env python3
"""v9.18.5 campaign: v9.18.2 completion plus target-stop quality corridor."""
from __future__ import annotations

import run_mpz_v9_18_1_persistent_plastic_wake as _v9181
import run_mpz_v9_18_2_persistent_plastic_wake as _v9182


def _build_command_v9185(args, class_name, run_root, force_rerun):
    cmd = _v9181._original_build(args, class_name, run_root, force_rerun)
    old = "run_mpz_v9_18_mode_i_rcurve.py"
    new = "run_mpz_v9_18_5_mode_i_rcurve.py"
    try:
        cmd[cmd.index(old)] = new
    except ValueError as exc:
        raise RuntimeError(f"v9.18 driver token not found in command: {cmd}") from exc
    return cmd


def main():
    original = _v9181._build_command_v9181
    _v9181._build_command_v9181 = _build_command_v9185
    try:
        return _v9182.main()
    finally:
        _v9181._build_command_v9181 = original


if __name__ == "__main__":
    main()
