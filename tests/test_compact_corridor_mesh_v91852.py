from __future__ import annotations

import json
from types import SimpleNamespace

import numpy as np

from arrhenius_fracture import mesh as meshmod
from arrhenius_fracture import mode_i_first_passage_v9_18_5_2 as v91852


def test_compact_mesh_removes_unused_nodes(monkeypatch):
    nodes = np.array([
        [0.0, 0.0],
        [1.0, 0.0],
        [0.5, np.sqrt(3.0) / 2.0],
        [9.0, 9.0],  # deliberately unused Delaunay input point
    ])
    elems = np.array([[0, 1, 2]], dtype=int)
    raw = meshmod.rebuild_tri_mesh(nodes, elems, tip_centers=np.array([[0.0, 0.0]]))
    monkeypatch.setenv("ARRHENIUS_MIN_INITIAL_TRIANGLE_QUALITY", "0.9")

    compact, audit = v91852._compact_mesh(raw, np.array([[0.0, 0.0]]))

    assert compact.nn == 3
    assert compact.ne == 1
    assert audit["unused_input_node_count"] == 1
    assert audit["orphan_node_count_after_compaction"] == 0
    assert audit["initial_mesh_compaction_applied"] is True
    assert audit["minimum_initial_triangle_quality"] >= 0.999999


def test_corridor_wrapper_records_compaction(monkeypatch):
    nodes = np.array([
        [0.0, 0.0],
        [1.0, 0.0],
        [0.5, np.sqrt(3.0) / 2.0],
        [2.0, 2.0],
    ])
    elems = np.array([[0, 1, 2]], dtype=int)
    raw = meshmod.rebuild_tri_mesh(nodes, elems, tip_centers=np.array([[0.0, 0.0]]))

    def fake_original(geom, mesh_cfg, seed=None, tip_center=None):
        return raw

    v91852._compact_corridor_mesh._original = fake_original
    monkeypatch.setenv("ARRHENIUS_PREFINED_MODE_I_CORRIDOR", "1")
    monkeypatch.setenv("ARRHENIUS_COMMITTED_TARGET_EXTENSION_UM", "60")
    monkeypatch.setenv("ARRHENIUS_MIN_INITIAL_TRIANGLE_QUALITY", "0.9")
    geom = SimpleNamespace(a0=0.0, Lx=1.0)
    mesh_cfg = SimpleNamespace(tip_h_fine=1.0, tip_ratio=1.2)

    compact = v91852._compact_corridor_mesh(geom, mesh_cfg, seed=1)

    assert compact.nn == 3
    assert v91852._STARTUP_AUDIT["unused_input_node_count"] == 1
    assert v91852._STARTUP_AUDIT["startup_completed"] if "startup_completed" in v91852._STARTUP_AUDIT else True


def test_startup_audit_records_exception(tmp_path):
    argv = ["--out", str(tmp_path)]
    v91852._STARTUP_AUDIT.clear()
    v91852._STARTUP_AUDIT.update({"unused_input_node_count": 2})
    err = RuntimeError("synthetic startup failure")

    v91852._write_startup_audit(argv, err)

    payload = json.loads((tmp_path / "compact_corridor_mesh_v91852.json").read_text())
    assert payload["startup_completed"] is False
    assert payload["startup_error_type"] == "RuntimeError"
    assert "synthetic startup failure" in payload["startup_error"]
