#!/usr/bin/env python3
"""Run the v9.18 persistent-plastic-wake material-transfer gate."""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

import run_mpz_v9_17_hazard_clock_source_refresh as _v917
from arrhenius_fracture.coupled_event_audit_v918 import audit_case_v918

_original_build = _v917._build_command_v917
_original_run = _v917._run_case_v917


def _build_command_v918(args, class_name, run_root, force_rerun):
    cmd = _original_build(args, class_name, run_root, force_rerun)
    old = "run_mpz_v9_17_mode_i_rcurve.py"
    new = "run_mpz_v9_18_mode_i_rcurve.py"
    try:
        cmd[cmd.index(old)] = new
    except ValueError as exc:
        raise RuntimeError(f"v9.17 driver token not found in command: {cmd}") from exc
    return cmd


def _run_case_v918(args, base_seed, class_name, root):
    row = _original_run(args, base_seed, class_name, root)
    case_dir = Path(row["case_dir"])
    audit = audit_case_v918(case_dir, args.T_K, args.target_extension_um)
    row.update({f"v918_{k}": v for k, v in audit.items() if k not in {"case_dir", "T_K"}})
    row["analysis_committed_extension_um"] = float(
        audit.get("committed_extension_um", 0.0) or 0.0
    )
    row["target_completed"] = bool(audit.get("target_committed", False))
    row["final_extension_um"] = row["analysis_committed_extension_um"]
    rc = int(row.get("subprocess_returncode", row.get("returncode", 0)) or 0)
    if rc == 0:
        row["status"] = (
            "complete"
            if bool(audit.get("target_committed", False))
            else "right_censored_persistent_wake_uncommitted"
        )
    (case_dir / "v9_18_case_summary.json").write_text(
        json.dumps(row, indent=2, default=str)
    )
    pd.DataFrame([row]).to_csv(case_dir / "v9_18_case_summary.csv", index=False)
    return row


def main():
    original_build = _v917._build_command_v917
    original_run = _v917._run_case_v917
    _v917._build_command_v917 = _build_command_v918
    _v917._run_case_v917 = _run_case_v918
    try:
        return _v917.main()
    finally:
        _v917._build_command_v917 = original_build
        _v917._run_case_v917 = original_run


if __name__ == "__main__":
    main()
