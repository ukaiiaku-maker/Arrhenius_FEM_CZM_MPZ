from __future__ import annotations

import json

import numpy as np
import pytest

from arrhenius_fracture.config import GeometryConfig, MeshConfig
from arrhenius_fracture import mode_i_first_passage_v9_18_5_2 as v91852
from arrhenius_fracture import mode_i_first_passage_v9_18_5_3 as v91853
from arrhenius_fracture.mode_i_first_passage_v10_0_5_13_2_barrier_only import (
    _quality_selected_corridor_mesh_v1005132,
)
from arrhenius_fracture.mode_i_first_passage_v10_0_5_13_5_barrier_only import (
    expanded_candidate_counts_v1005135,
)
from arrhenius_fracture.physical_refinement_mesh_v100510 import (
    clear_physical_refinement_v100510,
    configure_physical_refinement_v100510,
    make_physical_refinement_mesh_v100510,
)


def test_expanded_counts_include_low_count_long_corridors():
    counts = expanded_candidate_counts_v1005135(110.0, 35.0)
    assert counts[0] == 2
    assert counts[-1] >= 8


def test_real_100um_physical_refinement_corridor_meets_quality_floor(monkeypatch):
    """Exercise the actual 330 um radial-cloud/Delaunay startup construction."""
    geom = GeometryConfig(Lx=2.0e-3, Ly=4.0e-3, a0=0.5e-3)
    mesh_cfg = MeshConfig(
        nx=36,
        ny=72,
        tip_h_fine=2.5e-6,
        tip_ratio=1.15,
    )
    monkeypatch.setenv("ARRHENIUS_COMMITTED_TARGET_EXTENSION_UM", "100")
    monkeypatch.setenv("ARRHENIUS_CORRIDOR_GUARD_UM", "10")
    monkeypatch.setenv("ARRHENIUS_CORRIDOR_MAX_CENTER_GAP_UM", "35")
    monkeypatch.setenv("ARRHENIUS_PHYSICAL_DA_UM", "5")
    monkeypatch.setenv("ARRHENIUS_MIN_ACCEPTED_TRIANGLE_QUALITY", "0.035")
    monkeypatch.setenv("ARRHENIUS_MAX_TIP_H_OVER_DA", "0.75")
    monkeypatch.setenv("ARRHENIUS_PREFINED_MODE_I_CORRIDOR", "1")

    fn = _quality_selected_corridor_mesh_v1005132
    had_original = hasattr(fn, "_original")
    saved_original = getattr(fn, "_original", None)
    saved_counts = v91853._candidate_counts
    fn._original = make_physical_refinement_mesh_v100510
    v91853._candidate_counts = expanded_candidate_counts_v1005135
    configure_physical_refinement_v100510(330.0e-6)
    v91852._STARTUP_AUDIT.clear()
    try:
        try:
            mesh = fn(geom, mesh_cfg, seed=42, tip_center=None)
        except Exception as exc:
            pytest.fail(
                f"100 um corridor construction failed: {exc}\n"
                + json.dumps(v91852._STARTUP_AUDIT, indent=2, default=str)
            )
    finally:
        clear_physical_refinement_v100510()
        v91853._candidate_counts = saved_counts
        if had_original:
            fn._original = saved_original
        else:
            delattr(fn, "_original")

    quality = v91852._triangle_quality(mesh.nodes, mesh.elems)
    assert np.all(np.isfinite(mesh.area_e))
    assert np.all(mesh.area_e > 0.0)
    assert float(np.min(quality)) >= 0.035
    audit = v91852._STARTUP_AUDIT
    assert audit["startup_resolution_warning"] in (True, False)
    assert audit["tip_h_over_da_enforced_as_veto"] is False
    assert audit["selected_center_count"] >= 2
    assert any(row.get("center_count") == 2 for row in audit["candidate_corridors"])
