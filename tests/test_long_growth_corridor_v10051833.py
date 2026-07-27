from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest

from arrhenius_fracture import mesh as mesh_module
from arrhenius_fracture.config import GeometryConfig, MeshConfig
from arrhenius_fracture import mode_i_first_passage_v9_18_5_2 as v91852
from arrhenius_fracture import mode_i_first_passage_v9_18_5_3 as v91853
from arrhenius_fracture import mode_i_first_passage_v10_0_5_18_3_3_four_class_long_growth_corridor as entry
from arrhenius_fracture.long_growth_corridor_v10051833 import (
    CORRIDOR_SCHEMA,
    target_aware_long_growth_corridor_mesh,
)


def _set_corridor_env(monkeypatch, target_um="1000"):
    monkeypatch.setenv("ARRHENIUS_COMMITTED_TARGET_EXTENSION_UM", target_um)
    monkeypatch.setenv("ARRHENIUS_PHYSICAL_DA_UM", "5")
    monkeypatch.setenv("ARRHENIUS_CORRIDOR_PROCESS_ZONE_UM", "50")
    monkeypatch.setenv("ARRHENIUS_CORRIDOR_GUARD_UM", "10")
    monkeypatch.setenv("ARRHENIUS_CORRIDOR_MAX_CENTER_GAP_UM", "100")
    monkeypatch.setenv("ARRHENIUS_MAX_CORRIDOR_H_OVER_LPZ", "0.25")
    monkeypatch.setenv("ARRHENIUS_MIN_INITIAL_TRIANGLE_QUALITY", "0.035")
    monkeypatch.setenv("ARRHENIUS_MAX_TIP_H_OVER_DA", "0.75")


def test_1000um_corridor_is_process_zone_resolved_and_compact(monkeypatch):
    _set_corridor_env(monkeypatch)
    geom = GeometryConfig()
    cfg = MeshConfig(nx=36, ny=72, tip_h_fine=2.5e-6, tip_ratio=1.15)
    target_aware_long_growth_corridor_mesh._original = mesh_module.make_tri_mesh
    mesh = target_aware_long_growth_corridor_mesh(geom, cfg, seed=1)
    audit = dict(v91852._STARTUP_AUDIT)

    assert audit["schema"] == CORRIDOR_SCHEMA
    assert audit["corridor_target_extension_um"] == pytest.approx(1000.0)
    assert audit["corridor_guard_um"] == pytest.approx(10.0)
    assert audit["full_requested_corridor_covered"] is True
    assert audit["target_propagated_before_mesh_construction"] is True
    assert audit["minimum_initial_triangle_quality"] >= 0.035
    assert audit["maximum_sampled_hbar_tip_over_L_pz"] <= 0.25
    assert audit["tip_h_over_da_enforced_as_veto"] is False
    assert 2 < audit["selected_center_count"] < 20
    assert mesh.nn < 6000
    assert mesh.ne < 12000
    assert np.all(np.isfinite(mesh.area_e))
    assert np.all(mesh.area_e > 0.0)


def test_entrypoint_propagates_target_before_base(monkeypatch, tmp_path):
    observed = {}

    def fake_base(argv):
        observed["target"] = os.environ.get("ARRHENIUS_COMMITTED_TARGET_EXTENSION_UM")
        observed["da"] = os.environ.get("ARRHENIUS_PHYSICAL_DA_UM")
        observed["lpz"] = os.environ.get("ARRHENIUS_CORRIDOR_PROCESS_ZONE_UM")
        observed["gap"] = os.environ.get("ARRHENIUS_CORRIDOR_MAX_CENTER_GAP_UM")
        observed["h_lpz"] = os.environ.get("ARRHENIUS_MAX_CORRIDOR_H_OVER_LPZ")
        observed["corridor"] = v91853._quality_selected_corridor_mesh
        return "ok"

    monkeypatch.setattr(entry._base, "main", fake_base)
    original = v91853._quality_selected_corridor_mesh
    result = entry.main([
        "--target-crack-extension-um", "1000",
        "--da-phys", "5e-6",
        "--mpz-length-um", "50",
        "--out", str(tmp_path),
    ])

    assert result == "ok"
    assert observed["target"] == "1000.0"
    assert float(observed["da"]) == pytest.approx(5.0e-6)
    assert observed["lpz"] == "50.0"
    assert observed["gap"] == "100"
    assert observed["h_lpz"] == "0.25"
    assert observed["corridor"] is target_aware_long_growth_corridor_mesh
    assert v91853._quality_selected_corridor_mesh is original


def test_entrypoint_rejects_missing_target():
    with pytest.raises(SystemExit, match="target-crack-extension"):
        entry.main(["--da-phys", "5e-6", "--mpz-length-um", "50"])


def test_runner_keeps_child_area_floor_and_uses_new_entrypoint():
    root = Path(__file__).resolve().parents[1]
    text = (root / "run_v10_0_5_18_3_3_four_class_focused_long_growth.sh").read_text()
    assert "mode_i_first_passage_v10_0_5_18_3_3_four_class_long_growth_corridor" in text
    assert "ARRHENIUS_CORRIDOR_MAX_CENTER_GAP_UM" in text
    assert "ARRHENIUS_MAX_CORRIDOR_H_OVER_LPZ" in text
    assert "child_area_ratio_floor_relaxed=false" in text
