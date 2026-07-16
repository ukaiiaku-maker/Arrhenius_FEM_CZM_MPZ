from __future__ import annotations

import json
from pathlib import Path

from run_mpz_v9_18_2_persistent_plastic_wake import promote_committed_completion


def test_expected_topology_guard_exit_is_promoted(tmp_path: Path):
    (tmp_path / "rcurve_run_audit.json").write_text(json.dumps({
        "target_extension_um": 15.0,
        "driver": "legacy_topology_guard",
    }))
    (tmp_path / "v9_13_run_config.json").write_text(json.dumps({
        "target_extension_um": 10.0,
    }))
    row = {
        "case_dir": str(tmp_path),
        "status": "right_censored",
        "returncode": 0,
        "subprocess_returncode": 1,
        "completion_marker_present": False,
        "solver_guard_target_extension_um": 15.0,
        "analysis_committed_extension_um": 10.0,
        "v918_target_committed": True,
        "v918_no_uncommitted_trial_at_exit": True,
        "v918_persistent_wake_commit_gate_passed": True,
    }

    out = promote_committed_completion(row, 10.0)

    assert out["status"] == "complete"
    assert out["subprocess_returncode"] == 0
    assert out["legacy_topology_guard_returncode"] == 1
    assert out["committed_completion_promoted_v9182"]
    assert (tmp_path / ".long_growth_complete").exists()
    audit = json.loads((tmp_path / "rcurve_run_audit.json").read_text())
    assert audit["target_extension_um"] == 10.0
    assert audit["legacy_topology_guard_target_extension_um"] == 15.0
    summary = json.loads((tmp_path / "v9_13_case_summary.json").read_text())
    assert summary["subprocess_returncode"] == 0
    assert summary["completion_basis"] == "v918_committed_event_audit"


def test_uncommitted_or_active_event_is_not_promoted(tmp_path: Path):
    row = {
        "case_dir": str(tmp_path),
        "returncode": 0,
        "subprocess_returncode": 1,
        "v918_target_committed": True,
        "v918_no_uncommitted_trial_at_exit": False,
        "v918_persistent_wake_commit_gate_passed": True,
    }
    out = promote_committed_completion(row, 10.0)
    assert out["subprocess_returncode"] == 1
    assert not (tmp_path / ".long_growth_complete").exists()
