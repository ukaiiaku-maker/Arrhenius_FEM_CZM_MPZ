from __future__ import annotations

from argparse import Namespace
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from arrhenius_fracture import crack_backend
from arrhenius_fracture import mode_i_first_passage_v10_0_5_12_phase_c as old_entry
from arrhenius_fracture import mode_i_first_passage_v10_0_5_12_3_phase_c as entry
from arrhenius_fracture import mode_i_first_passage_v9_18_5 as v9185
from arrhenius_fracture.mpz_response_registry_v100512 import PARAMETER_SOURCE
import run_v10_0_5_12_3_phase_c_monotonic as campaign


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


def test_repair_command_selects_new_entry_and_correct_cluster_length(tmp_path: Path):
    command = campaign.build_command(
        "/example/python",
        _campaign_args(),
        "dbtt_primary",
        700,
        50.0,
        tmp_path / "case",
    )
    assert campaign.ENTRY_MODULE in command
    index = command.index("--rJ-cluster")
    assert float(command[index + 1]) == pytest.approx(30.0e-6)


def test_numpy_integration_compatibility_and_summary_metric(monkeypatch):
    monkeypatch.delattr(np, "trapz", raising=False)
    campaign.ensure_numpy_trapz_compat()
    assert hasattr(np, "trapz")
    rc = pd.DataFrame(
        {
            "crack_extension_um": [0.0, 50.0],
            "KJ_MPa_sqrt_m": [1.0, 3.0],
        }
    )
    result = campaign._base.metrics(rc, 50.0, 1.0)
    assert result["normalized_R_curve_area_MPa_sqrt_m"] == pytest.approx(2.0)


def test_entry_preserves_refinement_metadata_across_topology_rebuild(monkeypatch, tmp_path: Path):
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

    monkeypatch.setattr(old_entry._v10052, "main", fake_solver)
    out = tmp_path / "case"
    result = entry.main(
        [
            "--phase-c-option",
            "dbtt_primary",
            "--tip-refinement-radius-um",
            "330",
            "--selected-cluster-J-outer-um",
            "240",
            "--local-J-outer-um",
            "100",
            "--v10-material-source",
            PARAMETER_SOURCE,
            "--czm-opening-coupling",
            "clock_linear",
            "--mpz-length-um",
            "50",
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
    assert payload["point_release"] == "10.0.5.12.3"
    assert payload["run_completed_without_exception"] is True
    assert payload["mesh_refinement_runtime"]["actual_radius_verified"] is True
    assert payload["metadata_propagation_fix"]["physics_changed"] is False
