from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from arrhenius_fracture.field_snapshots_v913 import map_mpz_density_to_elements, render_field_snapshots_v913


def _snapshot():
    nodes = np.array([
        [0.0, -5e-6], [1e-5, -5e-6], [2e-5, -5e-6],
        [0.0, 5e-6], [1e-5, 5e-6], [2e-5, 5e-6],
    ])
    elems = np.array([[0, 1, 3], [1, 4, 3], [1, 2, 4], [2, 5, 4]], dtype=int)
    state = {
        "config": {"length_m": 2e-5, "n_bins": 2, "blunting_length_m": 5e-7},
        "mobile": [[1.0, 0.0], [0.0, 1.0]],
        "retained": [[2.0, 0.0], [0.0, 2.0]],
        "emitted_total": 8.0,
    }
    return {
        "step": 1, "KJ": 2e7, "a_tip": 0.0,
        "nodes": nodes, "elems": elems,
        "d": np.zeros(len(nodes)), "rho_gp": np.full(len(elems), 5e12),
        "s1_gp": np.linspace(1e8, 4e8, len(elems)),
        "epeq_gp": np.linspace(0.0, 1e-4, len(elems)),
        "front_paths": [(0, -1, np.array([[0.0, 0.0], [1e-5, 0.0]]))],
        "mpz_front_states": [{
            "front_id": 0, "xy_m": [0.0, 0.0], "direction": [1.0, 0.0], "state": state,
        }],
    }


def test_resolution_aware_mapping_is_nonzero_and_reports_coarse_graining():
    mapped, metadata = map_mpz_density_to_elements(_snapshot())
    assert mapped.shape == (4,)
    assert np.max(mapped) > 0.0
    assert metadata
    assert metadata[0]["display_coarse_grain_width_m"] >= metadata[0]["physical_width_m"]
    assert metadata[0]["active_line_count"] == 6.0
    assert metadata[0]["emitted_total"] == 8.0


def test_renderer_writes_full_and_tip_zoom_images(tmp_path):
    snap = _snapshot()
    mesh = SimpleNamespace(nodes=snap["nodes"], elems=snap["elems"])
    path = render_field_snapshots_v913(tmp_path, 700.0, mesh, [snap], max_cols=1)
    assert path is not None
    assert (tmp_path / "field_snapshots_700K.png").exists()
    assert (tmp_path / "field_snapshots_tip_zoom_700K.png").exists()
    assert (tmp_path / "field_snapshot_manifest_700K.json").exists()
    arrays = list((tmp_path / "field_snapshot_arrays_700K").glob("*.npz"))
    assert len(arrays) == 1
    payload = np.load(arrays[0])
    assert float(payload["emitted_total"]) == 8.0
    assert float(payload["retained_count"]) == 4.0
    assert float(payload["mobile_count"]) == 2.0
