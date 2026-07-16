#!/usr/bin/env python3
"""Run the v9.17 absolute-hazard/source-refresh material-transfer gate."""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

import run_mpz_v9_16_kinetic_trial_opening as _v916
from arrhenius_fracture.coupled_event_audit_v917 import audit_case_v917

_original_build = _v916._build_command_v916
_original_run = _v916._run_case_v916


def _build_command_v917(args, class_name, run_root, force_rerun):
    cmd = _original_build(args, class_name, run_root, force_rerun)
    old = "run_mpz_v9_16_mode_i_rcurve.py"
    new = "run_mpz_v9_17_mode_i_rcurve.py"
    try:
        cmd[cmd.index(old)] = new
    except ValueError as exc:
        raise RuntimeError(f"v9.16 driver token not found in command: {cmd}") from exc
    return cmd


def _run_case_v917(args, base_seed, class_name, root):
    row = _original_run(args, base_seed, class_name, root)
    case_dir = Path(row["case_dir"])
    audit = audit_case_v917(case_dir, args.T_K, args.target_extension_um)
    row.update({f"v917_{k}": v for k, v in audit.items() if k not in {"case_dir", "T_K"}})
    row["analysis_committed_extension_um"] = float(audit.get("committed_extension_um", 0.0) or 0.0)
    row["target_completed"] = bool(audit.get("target_committed", False))
    row["final_extension_um"] = row["analysis_committed_extension_um"]
    rc = int(row.get("subprocess_returncode", row.get("returncode", 0)) or 0)
    if rc == 0:
        row["status"] = (
            "complete"
            if bool(audit.get("target_committed", False))
            else "right_censored_hazard_opening_uncommitted"
        )
    (case_dir / "v9_17_case_summary.json").write_text(
        json.dumps(row, indent=2, default=str)
    )
    pd.DataFrame([row]).to_csv(case_dir / "v9_17_case_summary.csv", index=False)
    return row


def main():
    original_build = _v916._build_command_v916
    original_run = _v916._run_case_v916
    _v916._build_command_v916 = _build_command_v917
    _v916._run_case_v916 = _run_case_v917
    try:
        return _v916.main()
    finally:
        _v916._build_command_v916 = original_build
        _v916._run_case_v916 = original_run


if __name__ == "__main__":
    main()
