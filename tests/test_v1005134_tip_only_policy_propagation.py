from __future__ import annotations

from pathlib import Path
import subprocess
from types import SimpleNamespace

import pytest

from arrhenius_fracture import barrier_only_response_registry_v100513 as registry
from arrhenius_fracture import mode_i_first_passage_v10_0_5_13_4_barrier_only as entry
import run_v10_0_5_13_4_barrier_only_monotonic as runner


def test_full_wrapper_chain_receives_tip_only_policy_and_restores_it(monkeypatch, tmp_path):
    original_core_policy = entry._CORE_ENTRY.TWO_D_STATE_POLICY
    original_preserved_policy = entry._PRESERVED_ENTRY.TWO_D_STATE_POLICY
    original_registry_policy = registry.TWO_D_STATE_POLICY
    original_compact = entry._CORE_ENTRY._compact_audit
    observed = {}

    def fake_core_main(argv):
        observed["core_mode"] = entry._CORE_ENTRY.TWO_D_STATE_POLICY[
            "bulk_plasticity_mode"
        ]
        observed["preserved_mode"] = entry._PRESERVED_ENTRY.TWO_D_STATE_POLICY[
            "bulk_plasticity_mode"
        ]
        observed["registry_mode"] = registry.TWO_D_STATE_POLICY[
            "bulk_plasticity_mode"
        ]
        remaining = list(argv)
        entry._CORE_ENTRY._require_text_or_set(
            remaining,
            "--bulk-plasticity-mode",
            entry._CORE_ENTRY.TWO_D_STATE_POLICY["bulk_plasticity_mode"],
        )
        observed["validated_mode"] = remaining[
            remaining.index("--bulk-plasticity-mode") + 1
        ]
        audit = entry._CORE_ENTRY._compact_audit(
            {
                "candidate_id": "candidate",
                "target_class": "DBTT",
                "barrier_fingerprint_sha256": "barrier",
                "parameter_fingerprint_sha256": "parameter",
                "barrier_fields_transferred": [],
                "candidate_state_fields_ignored": {},
            }
        )
        observed["audit"] = audit
        return [{"T_K": 700.0}]

    monkeypatch.setattr(entry._CORE_ENTRY, "main", fake_core_main)
    result = entry.main(
        [
            "--out",
            str(tmp_path),
            "--bulk-plasticity-mode",
            "tip_only",
        ]
    )
    assert result == [{"T_K": 700.0}]
    assert observed["core_mode"] == "tip_only"
    assert observed["preserved_mode"] == "tip_only"
    assert observed["registry_mode"] == "tip_only"
    assert observed["validated_mode"] == "tip_only"
    assert observed["audit"]["bulk_state_evolves_in_fem"] is False
    assert observed["audit"]["moving_crack_tip_mpz_active"] is True
    assert observed["audit"]["uniform_bulk_mobile_retained_state_active"] is False

    assert entry._CORE_ENTRY.TWO_D_STATE_POLICY is original_core_policy
    assert entry._PRESERVED_ENTRY.TWO_D_STATE_POLICY is original_preserved_policy
    assert registry.TWO_D_STATE_POLICY is original_registry_policy
    assert entry._CORE_ENTRY._compact_audit is original_compact


def test_runner_routes_to_v1005134_entry(monkeypatch, tmp_path: Path):
    base_cmd = [
        "python",
        "-m",
        "arrhenius_fracture.mode_i_first_passage_v10_0_5_13_3_barrier_only",
        "--bulk-plasticity-mode",
        "tip_only",
    ]
    monkeypatch.setattr(runner, "_ORIGINAL_BUILD", lambda *args: list(base_cmd))
    cmd = runner._build_command(
        "python", SimpleNamespace(), "dbtt_primary", 700, 20.0, tmp_path
    )
    assert runner.ENTRY_MODULE in cmd
    assert (
        "arrhenius_fracture.mode_i_first_passage_v10_0_5_13_3_barrier_only"
        not in cmd
    )
    assert cmd[cmd.index("--bulk-plasticity-mode") + 1] == "tip_only"


def test_entry_rejects_non_tip_only_mode(tmp_path: Path):
    with pytest.raises(SystemExit, match="requires --bulk-plasticity-mode tip_only"):
        entry.main(
            [
                "--out",
                str(tmp_path),
                "--bulk-plasticity-mode",
                "bulk_same_pt_km",
            ]
        )


def test_shell_launcher_version_and_syntax():
    path = Path("run_v10_0_5_13_4_barrier_only_monotonic.sh")
    text = path.read_text()
    assert "Point release: 10.0.5.13.4" in text
    assert "Policy propagation: tip_only enforced through v10.0.5.13 core entry" in text
    assert "run_v10_0_5_13_4_barrier_only_monotonic.py" in text
    assert "DU=${DU:-2e-5}" in text
    assert "DT=${DT:-840}" in text
    cp = subprocess.run(["bash", "-n", str(path)], capture_output=True, text=True)
    assert cp.returncode == 0, cp.stderr
