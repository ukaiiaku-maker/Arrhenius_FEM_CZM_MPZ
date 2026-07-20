from __future__ import annotations

from argparse import Namespace
import json
from pathlib import Path
import subprocess

import numpy as np
import pytest

from arrhenius_fracture import crack_backend
from arrhenius_fracture import mode_i_first_passage_v10_0_5_13_barrier_only as entry
from arrhenius_fracture import mode_i_first_passage_v9_18_5 as v9185
from arrhenius_fracture.barrier_only_response_registry_v100513 import (
    BARRIER_FIELDS,
    IGNORED_CANDIDATE_STATE_FIELDS,
    PRIMARY_OPTION_KEYS,
    TWO_D_STATE_POLICY,
    load_barrier_option,
)
import run_v10_0_5_13_barrier_only_monotonic as campaign


def _campaign_args() -> Namespace:
    return Namespace(
        registry=None,
        tip_refinement_radius_um=330.0,
        cluster_J_outer_um=240.0,
        local_J_outer_um=100.0,
        steps=50000,
        nx=36,
        ny=72,
        tip_h_fine=2.5e-6,
        tip_ratio=1.15,
        dU=2.0e-7,
        dt=8.4,
        n_stagger=2,
        print_every=25,
        adaptive_event_target=0.15,
        da_um=5.0,
        theta_deg=45.0,
        save_snapshots=3,
        snapshot_cols=3,
        snapshot_interval_um=50.0,
    )


def test_four_options_transfer_only_barrier_fields():
    assert tuple(PRIMARY_OPTION_KEYS) == (
        "ceramic_primary",
        "weakT_primary",
        "dbtt_primary",
        "peak_primary",
    )
    options = [load_barrier_option(key) for key in PRIMARY_OPTION_KEYS]
    assert len({opt.barrier_fingerprint_sha256 for opt in options}) == 4
    for option in options:
        assert set(option.barrier_row) == set(BARRIER_FIELDS)
        assert set(option.ignored_candidate_state).issubset(
            set(IGNORED_CANDIDATE_STATE_FIELDS)
        )
        assert "source_sites_per_system" not in option.barrier_row
        assert "source_refresh_length_um" not in option.barrier_row
        assert "c_blunt" not in option.barrier_row


def test_legacy_compatibility_row_uses_common_2d_state_not_candidate_state():
    option = load_barrier_option("dbtt_primary")
    assert option.ignored_candidate_state["source_sites_per_system"] == pytest.approx(
        141.0590567476921
    )
    row = option.legacy_row()
    assert row["source_sites_per_system"] == pytest.approx(
        TWO_D_STATE_POLICY["source_sites_per_system"]
    )
    assert row["source_sites_per_system"] != pytest.approx(
        option.ignored_candidate_state["source_sites_per_system"]
    )
    assert row["source_refresh_length_um"] == pytest.approx(
        TWO_D_STATE_POLICY["source_refresh_length_um"]
    )
    assert row["c_blunt"] == pytest.approx(TWO_D_STATE_POLICY["c_blunt"])
    assert row["two_d_state_policy_id"] == TWO_D_STATE_POLICY["policy_id"]


def test_full_command_uses_state_coupled_2d_policy(tmp_path: Path):
    command = campaign._build_command(
        "/example/python",
        _campaign_args(),
        "peak_primary",
        900,
        100.0,
        tmp_path / "case",
    )
    assert command[2] == "arrhenius_fracture.mode_i_first_passage_v10_0_5_13_barrier_only"
    assert command[command.index("--bulk-plasticity-mode") + 1] == "bulk_same_pt_km"
    assert float(command[command.index("--mpz-length-um") + 1]) == pytest.approx(100.0)
    assert int(command[command.index("--mpz-n-bins") + 1]) == 80
    assert float(command[command.index("--target-crack-extension-um") + 1]) == pytest.approx(100.0)
    assert float(command[command.index("--tip-refinement-radius-um") + 1]) == pytest.approx(330.0)
    assert float(command[command.index("--selected-cluster-J-outer-um") + 1]) == pytest.approx(240.0)
    assert float(entry.cluster_j_legacy_length_m(240.0)) == pytest.approx(30.0e-6)


def test_entry_preserves_refinement_metadata_and_records_barrier_only_contract(
    monkeypatch, tmp_path: Path
):
    original_rebuild = crack_backend.rebuild_tri_mesh

    def fake_solver(_argv):
        nodes = np.array(
            [
                [0.0, 0.0],
                [1.0e-3, 0.0],
                [0.0, 1.0e-3],
                [1.0e-3, 1.0e-3],
            ],
            dtype=float,
        )
        elems = np.array([[0, 1, 2], [1, 3, 2]], dtype=int)
        mesh = crack_backend.rebuild_tri_mesh(
            nodes,
            elems,
            tip_centers=np.array([[5.0e-4, 0.0]]),
        )
        v9185._RUNTIME["mesh"] = mesh
        return [{"T_K": 700.0}]

    monkeypatch.setattr(entry._solver, "main", fake_solver)
    out = tmp_path / "case"
    result = entry.main(
        [
            "--barrier-option",
            "dbtt_primary",
            "--tip-refinement-radius-um",
            "330",
            "--selected-cluster-J-outer-um",
            "240",
            "--local-J-outer-um",
            "100",
            "--mode",
            "2d",
            "--bulk-plasticity-mode",
            "bulk_same_pt_km",
            "--mpz-length-um",
            "100",
            "--mpz-n-bins",
            "80",
            "--max-fronts",
            "1",
            "--rJ-cluster",
            "30e-6",
            "--rJ-outer",
            "100e-6",
            "--out",
            str(out),
        ]
    )
    assert result == [{"T_K": 700.0}]
    assert crack_backend.rebuild_tri_mesh is original_rebuild
    payload = json.loads((out / entry.PRODUCTION_MANIFEST).read_text())
    assert payload["point_release"] == "10.0.5.13"
    assert payload["run_completed_without_exception"] is True
    assert payload["candidate_state_fields_applied"] is False
    assert payload["bulk_plasticity_mode"] == "bulk_same_pt_km"
    assert payload["mesh_refinement_runtime"]["actual_radius_verified"] is True
    ignored = payload["barrier_option"]["candidate_state_fields_ignored"]
    assert "source_sites_per_system" in ignored
    assert "source_refresh_length_um" in ignored
    assert "c_blunt" in ignored


def test_restart_gate_skips_only_verified_complete_state_coupled_case(
    monkeypatch, tmp_path: Path
):
    case = tmp_path / "dbtt_primary" / "T0700"
    case.mkdir(parents=True)
    (case / campaign.STATUS_FILE).write_text(
        json.dumps(
            {
                "status": "complete",
                "option_key": "dbtt_primary",
                "target_extension_um": 100.0,
            }
        )
    )
    (case / entry.PRODUCTION_MANIFEST).write_text(
        json.dumps(
            {
                "run_completed_without_exception": True,
                "candidate_state_fields_applied": False,
                "mesh_refinement_runtime": {"actual_radius_verified": True},
                "barrier_option": {"option_key": "dbtt_primary"},
            }
        )
    )
    bulk_path = case / "bulk_state_v9_11_summary.json"
    bulk_path.write_text(
        json.dumps(
            {
                "bulk_explicit_mobile_retained_state": True,
                "bulk_state_update_calls": 12,
            }
        )
    )
    monkeypatch.setattr(campaign, "completion_status", lambda *_: (True, 100.0))
    assert campaign._case_is_complete(case, "dbtt_primary", 100.0)

    production = json.loads((case / entry.PRODUCTION_MANIFEST).read_text())
    production["candidate_state_fields_applied"] = True
    (case / entry.PRODUCTION_MANIFEST).write_text(json.dumps(production))
    assert not campaign._case_is_complete(case, "dbtt_primary", 100.0)

    production["candidate_state_fields_applied"] = False
    (case / entry.PRODUCTION_MANIFEST).write_text(json.dumps(production))
    bulk = json.loads(bulk_path.read_text())
    bulk["bulk_state_update_calls"] = 0
    bulk_path.write_text(json.dumps(bulk))
    assert not campaign._case_is_complete(case, "dbtt_primary", 100.0)


def test_shell_launcher_has_valid_syntax_and_root_import_guard():
    path = "run_v10_0_5_13_barrier_only_monotonic.sh"
    completed = subprocess.run(["bash", "-n", path], capture_output=True, text=True)
    assert completed.returncode == 0, completed.stderr
    text = Path(path).read_text()
    assert 'ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)' in text
    assert "ACTUAL_PACKAGE" in text
    assert "EXPECTED_PACKAGE" in text
    assert "run_v10_0_5_13_barrier_only_monotonic.py" in text
