from __future__ import annotations

from argparse import Namespace
from types import SimpleNamespace

import numpy as np
import pytest

from arrhenius_fracture import mode_i_first_passage_v9_18_5 as v9185
from arrhenius_fracture import mode_i_first_passage_v9_18_5_2 as v91852
from arrhenius_fracture import mode_i_first_passage_v9_18_5_3 as v91853
from arrhenius_fracture import mode_i_first_passage_v10_0_5_13_2_barrier_only as entry
import run_v10_0_5_13_2_barrier_only_monotonic as campaign


def _install_fake_corridor(monkeypatch, *, qmin: float, hratio: float):
    mesh = SimpleNamespace(nn=41, ne=66)
    centers = np.array([[0.2, 0.0], [0.3, 0.0]], dtype=float)

    entry._quality_selected_corridor_mesh_v1005132._original = (
        lambda geom, mesh_cfg, seed=None, tip_center=None: object()
    )
    monkeypatch.setattr(v91853, "_candidate_counts", lambda *_: [2])
    monkeypatch.setattr(v91853, "_centers_for_count", lambda *_: centers.copy())
    monkeypatch.setattr(
        v91853,
        "_compact_without_quality_abort",
        lambda raw, tip_centers: (
            mesh,
            {
                "minimum_initial_triangle_quality": qmin,
                "minimum_initial_triangle_quality_required": 0.035,
            },
        ),
    )
    monkeypatch.setattr(
        v91853,
        "_corridor_resolution",
        lambda compact, start, stop, da_m: {
            "sample_x_m": [start, stop],
            "sample_hbar_tip_m": [hratio * da_m, hratio * da_m],
            "maximum_sampled_hbar_tip_m": hratio * da_m,
            "maximum_sampled_hbar_tip_over_da": hratio,
        },
    )
    return mesh


def _geom_and_mesh_cfg():
    return (
        SimpleNamespace(a0=0.2, Lx=1.0),
        SimpleNamespace(tip_h_fine=2.5e-6),
    )


def test_quality_valid_corridor_is_accepted_when_h_tip_over_da_warns(monkeypatch):
    mesh = _install_fake_corridor(monkeypatch, qmin=0.041, hratio=1.20)
    geom, mesh_cfg = _geom_and_mesh_cfg()
    monkeypatch.setenv("ARRHENIUS_COMMITTED_TARGET_EXTENSION_UM", "100")
    monkeypatch.setenv("ARRHENIUS_PHYSICAL_DA_UM", "5")
    monkeypatch.setenv("ARRHENIUS_MIN_INITIAL_TRIANGLE_QUALITY", "0.035")
    monkeypatch.setenv("ARRHENIUS_MAX_TIP_H_OVER_DA", "0.75")
    v91852._STARTUP_AUDIT.clear()

    selected = entry._quality_selected_corridor_mesh_v1005132(geom, mesh_cfg)

    assert selected is mesh
    audit = dict(v91852._STARTUP_AUDIT)
    assert audit["minimum_initial_triangle_quality"] == pytest.approx(0.041)
    assert audit["maximum_sampled_hbar_tip_over_da"] == pytest.approx(1.20)
    assert audit["tip_h_over_da_enforced_as_veto"] is False
    assert audit["tip_h_over_da_role"] == "audit_warning_only"
    assert audit["startup_resolution_warning"] is True
    assert audit["candidate_corridors"][0]["accepted"] is True
    assert audit["candidate_corridors"][0]["resolution_warning"] is True
    assert v9185._RUNTIME["mesh"] is mesh


def test_triangle_quality_floor_remains_a_fatal_startup_gate(monkeypatch):
    _install_fake_corridor(monkeypatch, qmin=0.030, hratio=0.40)
    geom, mesh_cfg = _geom_and_mesh_cfg()
    monkeypatch.setenv("ARRHENIUS_COMMITTED_TARGET_EXTENSION_UM", "100")
    monkeypatch.setenv("ARRHENIUS_PHYSICAL_DA_UM", "5")
    monkeypatch.setenv("ARRHENIUS_MIN_INITIAL_TRIANGLE_QUALITY", "0.035")
    monkeypatch.setenv("ARRHENIUS_MAX_TIP_H_OVER_DA", "0.75")
    v91852._STARTUP_AUDIT.clear()

    with pytest.raises(RuntimeError, match="initial triangle quality"):
        entry._quality_selected_corridor_mesh_v1005132(geom, mesh_cfg)

    audit = dict(v91852._STARTUP_AUDIT)
    assert audit["tip_h_over_da_enforced_as_veto"] is False
    assert audit["minimum_initial_triangle_quality_required"] == pytest.approx(0.035)


def test_campaign_routes_v1005132_entry(tmp_path):
    args = Namespace(
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
    cmd = campaign._build_command(
        "/example/python", args, "dbtt_primary", 700, 20.0, tmp_path / "case"
    )
    assert campaign.ENTRY_MODULE in cmd
    assert "arrhenius_fracture.mode_i_first_passage_v10_0_5_13_1_barrier_only" not in cmd
