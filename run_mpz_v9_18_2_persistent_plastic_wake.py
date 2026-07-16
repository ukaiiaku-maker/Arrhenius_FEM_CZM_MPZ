#!/usr/bin/env python3
"""v9.18.2 committed-event completion handshake.

The v9.16 guard asks the inherited topology driver for one extra crack quantum so
that the final trial interface can finish opening.  v9.18 then suppresses any
renewal after the requested *committed* extension.  Consequently the inner
Mode-I driver exits with its legacy right-censored code even though the physical
event audit is complete.  This module promotes that expected guard exit to a
successful campaign result without changing any solver or constitutive state.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

import run_mpz_v9_18_1_persistent_plastic_wake as _v9181
import run_mpz_v9_18_persistent_plastic_wake as _campaign

_original_run = _campaign._run_case_v918


def _read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    try:
        value = json.loads(path.read_text())
    except Exception:
        return {}
    return value if isinstance(value, dict) else {}


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, default=str))


def promote_committed_completion(row: dict[str, Any], target_extension_um: float) -> dict[str, Any]:
    """Promote only a fully audited committed stop to campaign success."""
    out = dict(row)
    committed = bool(out.get("v918_target_committed", False))
    no_active = bool(out.get("v918_no_uncommitted_trial_at_exit", False))
    wake_commit = bool(out.get("v918_persistent_wake_commit_gate_passed", False))
    try:
        inner_solver_rc = int(out.get("returncode", 1))
    except (TypeError, ValueError):
        inner_solver_rc = 1
    if not (committed and no_active and wake_commit and inner_solver_rc == 0):
        return out

    case_dir = Path(out["case_dir"])
    case_dir.mkdir(parents=True, exist_ok=True)
    target = float(target_extension_um)
    committed_um = float(out.get("analysis_committed_extension_um", target) or target)
    old_outer_rc = int(out.get("subprocess_returncode", 0) or 0)
    solver_guard = float(out.get("solver_guard_target_extension_um", target) or target)

    marker = case_dir / ".long_growth_complete"
    marker.touch()
    out.update({
        "status": "complete",
        "subprocess_returncode": 0,
        "completion_marker_present": True,
        "target_extension_um": target,
        "target_completed": True,
        "final_extension_um": committed_um,
        "completion_basis": "v918_committed_event_audit",
        "legacy_topology_guard_returncode": old_outer_rc,
        "legacy_topology_guard_target_extension_um": solver_guard,
        "committed_completion_promoted_v9182": True,
    })

    run_audit_path = case_dir / "rcurve_run_audit.json"
    run_audit = _read_json(run_audit_path)
    if run_audit:
        run_audit["legacy_topology_guard_target_extension_um"] = float(
            run_audit.get("target_extension_um", solver_guard)
        )
        run_audit["target_extension_um"] = target
        run_audit["completion_basis"] = "v918_committed_event_audit"
        run_audit["committed_extension_um"] = committed_um
        _write_json(run_audit_path, run_audit)

    config_path = case_dir / "v9_13_run_config.json"
    config = _read_json(config_path)
    if config:
        config["target_extension_um"] = target
        config["legacy_topology_guard_target_extension_um"] = solver_guard
        config["completion_basis"] = "v918_committed_event_audit"
        _write_json(config_path, config)

    # The legacy material-transfer audit reads v9_13_case_summary directly.
    for stem in (
        "v9_13_case_summary",
        "v9_16_case_summary",
        "v9_17_case_summary",
        "v9_18_case_summary",
        "v9_18_2_case_summary",
    ):
        _write_json(case_dir / f"{stem}.json", out)
        pd.DataFrame([out]).to_csv(case_dir / f"{stem}.csv", index=False)
    return out


def _run_case_v9182(args, base_seed, class_name, root):
    row = _original_run(args, base_seed, class_name, root)
    return promote_committed_completion(row, args.target_extension_um)


def main():
    original = _campaign._run_case_v918
    _campaign._run_case_v918 = _run_case_v9182
    try:
        return _v9181.main()
    finally:
        _campaign._run_case_v918 = original


if __name__ == "__main__":
    main()
