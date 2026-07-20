from __future__ import annotations

from pathlib import Path
import subprocess
from types import SimpleNamespace

import pytest

from arrhenius_fracture.barrier_only_response_registry_v100513 import (
    TWO_D_STATE_POLICY as LEGACY_POLICY,
)
from arrhenius_fracture import mode_i_first_passage_v10_0_5_13_3_barrier_only as entry
import run_v10_0_5_13_3_barrier_only_monotonic as runner


def test_campaign_policy_is_tip_only_and_legacy_policy_is_unchanged():
    policy = runner.TIP_ONLY_POLICY
    assert policy["bulk_plasticity_mode"] == "tip_only"
    assert policy["state_evolution_source"] == "existing_moving_crack_tip_MPZ"
    assert policy["continuum_bulk_role"] == "elastic_fem_only"
    assert policy["uniform_bulk_mobile_retained_state_active"] is False
    # v10.0.5.13.3 must not retroactively alter older release contracts.
    assert LEGACY_POLICY["bulk_plasticity_mode"] == "bulk_same_pt_km"
    assert LEGACY_POLICY["policy_id"] == "preserve_existing_full_2d_state_v100513"


def test_command_forces_tip_only_and_new_entry(monkeypatch, tmp_path: Path):
    base_cmd = [
        "python",
        "-m",
        "arrhenius_fracture.mode_i_first_passage_v10_0_5_13_2_barrier_only",
        "--bulk-plasticity-mode",
        "bulk_same_pt_km",
        "--dU",
        "2e-5",
        "--dt",
        "840",
    ]
    monkeypatch.setattr(runner, "_ORIGINAL_BUILD", lambda *args: list(base_cmd))
    cmd = runner._build_command(
        "python", SimpleNamespace(), "dbtt_primary", 700, 20.0, tmp_path
    )
    assert runner.ENTRY_MODULE in cmd
    mode_index = cmd.index("--bulk-plasticity-mode") + 1
    assert cmd[mode_index] == "tip_only"
    assert "bulk_same_pt_km" not in cmd


def test_entry_rejects_bulk_mode(tmp_path: Path):
    with pytest.raises(SystemExit, match="requires --bulk-plasticity-mode tip_only"):
        entry.main(
            [
                "--out",
                str(tmp_path),
                "--bulk-plasticity-mode",
                "bulk_same_pt_km",
            ]
        )


def test_tip_only_restart_gate(monkeypatch, tmp_path: Path):
    status = {
        "status": "complete",
        "option_key": "dbtt_primary",
        "target_extension_um": 20.0,
    }
    production = {
        "run_completed_without_exception": True,
        "candidate_state_fields_applied": False,
        "mesh_refinement_runtime": {"actual_radius_verified": True},
        "barrier_option": {"option_key": "dbtt_primary"},
        "plasticity_scope": {
            "bulk_plasticity_mode": "tip_only",
            "moving_crack_tip_mpz_active": True,
            "uniform_bulk_mobile_retained_state_active": False,
        },
    }
    integration = {
        "bulk_PT": {
            "mode": "tip_only",
            "explicit_mobile_retained_state": False,
            "source_interpretation": "moving_crack_tip_MPZ_only",
        }
    }
    mapping = {
        runner._ULTIMATE.STATUS_FILE: status,
        runner._ULTIMATE.PRODUCTION_MANIFEST: production,
        "mpz_v9_11_integration_audit.json": integration,
    }
    monkeypatch.setattr(
        runner._ULTIMATE,
        "_read_json",
        lambda path: mapping.get(Path(path).name, {}),
    )
    monkeypatch.setattr(
        runner._ULTIMATE,
        "completion_status",
        lambda case, target: (True, target),
    )
    assert runner._tip_only_case_is_complete(
        tmp_path, "dbtt_primary", 20.0
    )
    integration["bulk_PT"]["mode"] = "bulk_same_pt_km"
    assert not runner._tip_only_case_is_complete(
        tmp_path, "dbtt_primary", 20.0
    )


def test_shell_has_rate_preserving_macro_step_and_valid_syntax():
    path = Path("run_v10_0_5_13_3_barrier_only_monotonic.sh")
    text = path.read_text()
    assert "DU=${DU:-2e-5}" in text
    assert "DT=${DT:-840}" in text
    assert "BASE_DU=${BASE_DU:-2e-7}" in text
    assert "BASE_DT=${BASE_DT:-8.4}" in text
    assert "ALLOW_RATE_CHANGE" in text
    assert "Plasticity scope: tip_only" in text
    cp = subprocess.run(["bash", "-n", str(path)], capture_output=True, text=True)
    assert cp.returncode == 0, cp.stderr


def test_default_macro_step_preserves_historical_rate():
    assert (2.0e-5 / 840.0) == pytest.approx(2.0e-7 / 8.4, rel=1e-15)
