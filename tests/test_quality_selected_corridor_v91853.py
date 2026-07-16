from __future__ import annotations

import numpy as np

from arrhenius_fracture.config import GeometryConfig, MeshConfig
from arrhenius_fracture import mesh as meshmod
from arrhenius_fracture import mode_i_first_passage_v9_18_5_2 as v91852
from arrhenius_fracture import mode_i_first_passage_v9_18_5_3 as v91853


def test_candidate_counts_include_three_center_70um_corridor():
    assert 3 in v91853._candidate_counts(70.0, 35.0)


def test_quality_selected_60um_corridor_passes_production_floors(monkeypatch):
    monkeypatch.setenv("ARRHENIUS_COMMITTED_TARGET_EXTENSION_UM", "60")
    monkeypatch.setenv("ARRHENIUS_CORRIDOR_GUARD_UM", "10")
    monkeypatch.setenv("ARRHENIUS_CORRIDOR_MAX_CENTER_GAP_UM", "35")
    monkeypatch.setenv("ARRHENIUS_PHYSICAL_DA_UM", "5")
    monkeypatch.setenv("ARRHENIUS_MIN_INITIAL_TRIANGLE_QUALITY", "0.035")
    monkeypatch.setenv("ARRHENIUS_MAX_TIP_H_OVER_DA", "0.75")

    cfg = MeshConfig(nx=36, ny=72, tip_h_fine=1.0e-6, tip_ratio=1.20)
    geom = GeometryConfig()
    v91853._quality_selected_corridor_mesh._original = meshmod.make_tri_mesh
    selected = v91853._quality_selected_corridor_mesh(geom, cfg, seed=42)

    audit = dict(v91852._STARTUP_AUDIT)
    assert selected.nn > 0 and selected.ne > 0
    assert audit["minimum_initial_triangle_quality"] >= 0.035
    assert audit["maximum_sampled_hbar_tip_over_da"] <= 0.75
    assert audit["selected_center_count"] >= 2
    incidence = np.bincount(selected.elems.ravel(), minlength=selected.nn)
    assert np.all(incidence > 0)
