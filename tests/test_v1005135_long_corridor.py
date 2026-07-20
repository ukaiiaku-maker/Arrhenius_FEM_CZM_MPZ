from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from arrhenius_fracture.config import GeometryConfig, MeshConfig
from arrhenius_fracture import mode_i_first_passage_v9_18_5_2 as v91852
from arrhenius_fracture import mode_i_first_passage_v9_18_5_3 as v91853
from arrhenius_fracture.mode_i_first_passage_v10_0_5_13_2_barrier_only import (
    _quality_selected_corridor_mesh_v1005132,
)
from arrhenius_fracture.mode_i_first_passage_v10_0_5_13_5_barrier_only import (
    _union_find_keep_indices_v1005135,
    expanded_candidate_counts_v1005135,
    make_physical_refinement_mesh_v1005135,
)
from arrhenius_fracture.physical_refinement_mesh_v100510 import (
    clear_physical_refinement_v100510,
    configure_physical_refinement_v100510,
)

AUDIT_PATH = Path("v1005135_candidate_audit.json")


def _persist_audit(extra=None):
    payload = dict(v91852._STARTUP_AUDIT)
    if extra:
        payload.update(extra)
    AUDIT_PATH.write_text(json.dumps(payload, indent=2, default=str))


def _candidate_geometry_diagnostics(geom, mesh_cfg):
    length_m = 110.0e-6
    out = []
    for count in expanded_candidate_counts_v1005135(110.0, 35.0):
        centers = v91853._centers_for_count(geom, length_m, count)
        raw = make_physical_refinement_mesh_v1005135(
            geom, mesh_cfg, seed=42, tip_center=centers
        )
        compact, audit = v91853._compact_without_quality_abort(raw, centers)
        quality = v91852._triangle_quality(compact.nodes, compact.elems)
        index = int(np.argmin(quality))
        tri = compact.nodes[compact.elems[index]]
        edge_lengths = [
            float(np.linalg.norm(tri[1] - tri[0])),
            float(np.linalg.norm(tri[2] - tri[1])),
            float(np.linalg.norm(tri[0] - tri[2])),
        ]
        centroid = np.mean(tri, axis=0)
        out.append(
            {
                "center_count": int(count),
                "centers_m": centers.tolist(),
                "minimum_triangle_quality": float(quality[index]),
                "worst_triangle_index": index,
                "worst_triangle_nodes_m": tri.tolist(),
                "worst_triangle_centroid_m": centroid.tolist(),
                "worst_triangle_edge_lengths_m": edge_lengths,
                "worst_triangle_nearest_center_distance_m": float(
                    np.min(np.linalg.norm(centers - centroid[None, :], axis=1))
                ),
                "robust_dedup_removed_nodes": int(
                    getattr(raw, "production_refinement_removed_near_duplicates", 0)
                ),
                "compaction_audit": audit,
            }
        )
    return out


def test_expanded_counts_include_low_count_long_corridors():
    counts = expanded_candidate_counts_v1005135(110.0, 35.0)
    assert counts[0] == 2
    assert counts[-1] >= 8


def test_radius_dedup_merges_points_across_rounding_bin_boundary():
    geom = GeometryConfig()
    points = np.array(
        [
            [601.261762741483e-6, 0.0],
            [601.31875e-6, 0.0],
            [602.1784640155221e-6, 3.766653210214303e-6],
        ],
        dtype=float,
    )
    keep = _union_find_keep_indices_v1005135(
        points, tolerance_m=0.125e-6, protected_start=len(points), geom=geom
    )
    assert len(keep) == 2


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
    fn._original = make_physical_refinement_mesh_v1005135
    v91853._candidate_counts = expanded_candidate_counts_v1005135
    configure_physical_refinement_v100510(330.0e-6)
    v91852._STARTUP_AUDIT.clear()
    try:
        try:
            mesh = fn(geom, mesh_cfg, seed=42, tip_center=None)
        except Exception as exc:
            _persist_audit(
                {
                    "test_exception_type": type(exc).__name__,
                    "test_exception": str(exc),
                    "candidate_geometry_diagnostics": _candidate_geometry_diagnostics(
                        geom, mesh_cfg
                    ),
                }
            )
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
    _persist_audit(
        {
            "test_mesh_node_count": int(mesh.nn),
            "test_mesh_triangle_count": int(mesh.ne),
            "test_minimum_triangle_quality": float(np.min(quality)),
            "test_dedup_policy": getattr(
                mesh, "production_refinement_policy", None
            ),
            "test_removed_near_duplicates": int(
                getattr(mesh, "production_refinement_removed_near_duplicates", 0)
            ),
        }
    )
    assert np.all(np.isfinite(mesh.area_e))
    assert np.all(mesh.area_e > 0.0)
    assert float(np.min(quality)) >= 0.035
    assert getattr(mesh, "production_refinement_removed_near_duplicates", 0) > 0
    audit = v91852._STARTUP_AUDIT
    assert audit["startup_resolution_warning"] in (True, False)
    assert audit["tip_h_over_da_enforced_as_veto"] is False
    assert audit["selected_center_count"] >= 2
    assert any(row.get("center_count") == 2 for row in audit["candidate_corridors"])
