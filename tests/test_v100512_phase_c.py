from __future__ import annotations

from argparse import Namespace
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from arrhenius_fracture import mode_i_first_passage_v10_0_3_progressive as v1003
from arrhenius_fracture import mode_i_first_passage_v10_0_5_12_phase_c as entry
from arrhenius_fracture import mode_i_first_passage_v9_18_5 as v9185
from arrhenius_fracture.mpz_response_registry_v100512 import (
    EXPECTED_CANDIDATE_IDS,
    PARAMETER_SOURCE,
    PRIMARY_OPTION_KEYS,
    load_option,
    load_registry,
)
import run_v10_0_5_12_phase_c_monotonic as campaign


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
        save_snapshots=11,
        snapshot_cols=6,
        snapshot_interval_um=50.0,
    )


def test_authoritative_four_option_registry():
    registry = load_registry()
    assert tuple(registry) == PRIMARY_OPTION_KEYS
    assert {key: value.candidate_id for key, value in registry.items()} == EXPECTED_CANDIDATE_IDS
    assert (registry["ceramic_primary"].mpz_length_um, registry["ceramic_primary"].mpz_n_bins) == (100.0, 200)
    assert (registry["weakT_primary"].mpz_length_um, registry["weakT_primary"].mpz_n_bins) == (100.0, 200)
    assert (registry["dbtt_primary"].mpz_length_um, registry["dbtt_primary"].mpz_n_bins) == (50.0, 80)
    assert (registry["peak_primary"].mpz_length_um, registry["peak_primary"].mpz_n_bins) == (50.0, 80)
    assert all(float(option.row["mobile_shield_fraction"]) == 0.0 for option in registry.values())


def test_full_matrix_is_exactly_40_cases():
    options, temperatures, target = campaign.default_matrix("full")
    assert options == PRIMARY_OPTION_KEYS
    assert temperatures == tuple(range(300, 1201, 100))
    assert target == 500.0
    assert len(options) * len(temperatures) == 40


def test_command_uses_exact_option_grid_and_validated_mechanics(tmp_path: Path):
    args = _campaign_args()
    command = campaign.build_command(
        "/example/python",
        args,
        "peak_primary",
        900,
        500.0,
        tmp_path / "peak" / "T0900",
    )
    text = " ".join(command)
    assert "--phase-c-option peak_primary" in text
    assert "--mpz-length-um 50" in text
    assert "--mpz-n-bins 80" in text
    assert "--tip-refinement-radius-um 330" in text
    assert "--selected-cluster-J-outer-um 240" in text
    assert "--local-J-outer-um 100" in text
    assert "--tip-h-fine 2.5e-06" in text
    assert "--tip-ratio 1.15" in text
    assert "--max-fronts 1" in text
    assert "--bulk-plasticity-mode tip_only" in text
    assert "--adaptive-event-target 0.15" in text


def test_phase_c_entry_composes_registry_and_refinement(monkeypatch, tmp_path: Path):
    original_source = v1003.PF_SOURCE
    observed = {}

    def fake_solver(argv):
        observed["source"] = v1003.PF_SOURCE
        manifest = v1003.load_material_manifest(
            "peak_primary", parameter_source=PARAMETER_SOURCE
        )
        observed["candidate_id"] = manifest.candidate_id
        observed["mobile_fraction"] = v1003.KineticCampaignCZMConfig(
            mobile_shield_fraction=1.0
        ).mobile_shield_fraction
        v9185._RUNTIME["mesh"] = SimpleNamespace(
            production_refinement_radius_m=330.0e-6,
            production_refinement_centers_m=[[0.0005, 0.0]],
            production_refinement_policy="fixed_physical_radius_same_radial_ring_law",
        )
        return [{"T_K": 900.0}]

    monkeypatch.setattr(entry._v10052, "main", fake_solver)
    out = tmp_path / "case"
    result = entry.main(
        [
            "--phase-c-option", "peak_primary",
            "--tip-refinement-radius-um", "330",
            "--selected-cluster-J-outer-um", "240",
            "--local-J-outer-um", "100",
            "--v10-material-source", PARAMETER_SOURCE,
            "--czm-opening-coupling", "clock_linear",
            "--mpz-length-um", "50",
            "--mpz-n-bins", "80",
            "--max-fronts", "1",
            "--rJ-cluster", "240e-6",
            "--rJ-outer", "100e-6",
            "--out", str(out),
        ]
    )
    assert result == [{"T_K": 900.0}]
    assert observed == {
        "source": PARAMETER_SOURCE,
        "candidate_id": "DBTT_restart05_candidate61",
        "mobile_fraction": 0.0,
    }
    payload = json.loads((out / entry.PRODUCTION_MANIFEST).read_text())
    assert payload["run_completed_without_exception"] is True
    assert payload["option"]["option_key"] == "peak_primary"
    assert payload["mesh_refinement_runtime"]["actual_radius_verified"] is True
    assert v1003.PF_SOURCE == original_source


def test_unsupported_contour_fails_before_solver(monkeypatch, tmp_path: Path):
    called = False

    def fake_solver(_):
        nonlocal called
        called = True

    monkeypatch.setattr(entry._v10052, "main", fake_solver)
    with pytest.raises(SystemExit, match="not supported"):
        entry.main(
            [
                "--phase-c-option", "dbtt_primary",
                "--tip-refinement-radius-um", "200",
                "--selected-cluster-J-outer-um", "240",
                "--local-J-outer-um", "100",
                "--out", str(tmp_path / "unsupported"),
            ]
        )
    assert called is False


def test_registry_manifest_is_one_exact_row(tmp_path: Path):
    option = load_option("dbtt_primary")
    path = option.write_selected_csv(tmp_path / "selected.csv")
    lines = path.read_text().splitlines()
    assert len(lines) == 2
    assert option.candidate_id in lines[1]
