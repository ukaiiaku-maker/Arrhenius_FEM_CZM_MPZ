#!/usr/bin/env python3
"""Run v9.14 through the conservative adaptive-CZM event path."""
from __future__ import annotations

import json
from pathlib import Path

import run_mpz_v9_13_deterministic_material_transfer as _base
from arrhenius_fracture.remesh_audit_v914 import audit_case


_original_build_command = _base.build_command
_original_run_case = _base.run_case


def _build_command_v914(args, class_name, run_root, force_rerun):
    cmd = _original_build_command(args, class_name, run_root, force_rerun)
    old = "run_mpz_v9_13_mode_i_rcurve.py"
    new = "run_mpz_v9_14_mode_i_rcurve.py"
    try:
        cmd[cmd.index(old)] = new
    except ValueError as exc:
        raise RuntimeError(f"v9.13 driver token not found in command: {cmd}") from exc
    return cmd


def _run_case_v914(args, base_seed, class_name, root):
    row = _original_run_case(args, base_seed, class_name, root)
    case_dir = Path(row["case_dir"])
    audit = audit_case(case_dir, args.T_K)
    row.update({f"v914_{k}": v for k, v in audit.items() if k not in {"case_dir", "T_K"}})
    (case_dir / "v9_14_case_summary.json").write_text(json.dumps(row, indent=2, default=str))
    return row


def main():
    original_build = _base.build_command
    original_run = _base.run_case
    _base.build_command = _build_command_v914
    _base.run_case = _run_case_v914
    try:
        return _base.main()
    finally:
        _base.build_command = original_build
        _base.run_case = original_run


if __name__ == "__main__":
    main()
