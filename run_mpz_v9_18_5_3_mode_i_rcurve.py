#!/usr/bin/env python3
from __future__ import annotations
import run_mpz_v9_11_mode_i_rcurve_3T as _base

_original_build_command = _base.build_command

def _build_command_v91853(py, args, class_name, manifest, T_K, case_dir):
    cmd = _original_build_command(py, args, class_name, manifest, T_K, case_dir)
    old = "arrhenius_fracture.mode_i_first_passage_v9_11"
    new = "arrhenius_fracture.mode_i_first_passage_v9_18_5_3"
    cmd[cmd.index(old)] = new
    return cmd

def main():
    original = _base.build_command
    _base.build_command = _build_command_v91853
    try:
        return _base.main()
    finally:
        _base.build_command = original

if __name__ == "__main__":
    main()
