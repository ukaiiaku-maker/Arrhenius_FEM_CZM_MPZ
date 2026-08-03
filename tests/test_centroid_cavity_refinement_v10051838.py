from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from arrhenius_fracture.mesh import rebuild_tri_mesh
from arrhenius_fracture import quality_aware_czm_patch_v10051835 as base
from arrhenius_fracture import quality_aware_czm_centroid_cavity_v10051838 as centroid
from arrhenius_fracture import quality_aware_czm_ray_handoff_v10051835 as handoff
from arrhenius_fracture import (
    mode_i_first_passage_v10_0_5_18_3_8_centroid_cavity as entry
)


class Backend:
    def __init__(self):
        self.geom = SimpleNamespace(
            Lx=2.0,
            Ly=3.0,
            a0=0.0,
            notch_half_thickness=1.0e-3,
        )
        self.min_triangle_quality = 0.035
        self.min_area_ratio = 0.08

    @staticmethod
    def _tip_geometric_node_ids(mesh, p0, front_id):
        distance = np.linalg.norm(
            np.asarray(mesh.nodes) - np.asarray(p0)[None, :], axis=1
        )
        minimum = float(np.min(distance))
        return np.where(distance <= minimum + 1.0e-12)[0].astype(int).tolist()

    @staticmethod
    def _incident_elements(elems, node_ids):
        return np.where(
            np.any(np.isin(np.asarray(elems), list(node_ids)), axis=1)
        )[0]

    @staticmethod
    def _signed_twice_area(nodes, elems):
        tri = np.asarray(nodes)[np.asarray(elems)]
        return (
            (tri[:, 1, 0] - tri[:, 0, 0])
            * (tri[:, 2, 1] - tri[:, 0, 1])
            - (tri[:, 1, 1] - tri[:, 0, 1])
            * (tri[:, 2, 0] - tri[:, 0, 0])
        )

    @classmethod
    def _triangle_quality(cls, nodes, elems):
        tri = np.asarray(nodes)[np.asarray(elems)]
        l2 = (
            np.sum((tri[:, 1] - tri[:, 0]) ** 2, axis=1)
            + np.sum((tri[:, 2] - tri[:, 1]) ** 2, axis=1)
            + np.sum((tri[:, 0] - tri[:, 2]) ** 2, axis=1)
        )
        a2 = np.abs(cls._signed_twice_area(nodes, elems))
        return 2.0 * np.sqrt(3.0) * a2 / np.maximum(l2, 1.0e-300)


def _state():
    # Triangle 0 defines the current tip fan. Triangle 1 is a healthy target
    # parent sharing the full edge (nodes 1,2) with that fan.
    nodes = np.asarray(
        [
            [0.5, -0.8],
            [0.0, 0.0],
            [1.0, 0.0],
            [0.5, np.sqrt(3.0) / 2.0],
        ],
        dtype=float,
    )
    elems = np.asarray([[0, 2, 1], [1, 2, 3]], dtype=int)
    mesh = rebuild_tri_mesh(nodes, elems, tip_centers=[[0.5, -0.8]])
    return base.State(
        mesh=mesh,
        boundary=object(),
        damage=np.zeros(mesh.nn),
        displacement=np.zeros(mesh.ndof),
        parent_to_root=np.arange(mesh.ne, dtype=int),
    )


def test_repeated_centroid_split_recovers_near_edge_exact_endpoint(monkeypatch):
    monkeypatch.setenv("ARRHENIUS_MIN_ACCEPTED_TRIANGLE_QUALITY", "0.035")
    monkeypatch.setenv("ARRHENIUS_MIN_ACCEPTED_CHILD_AREA_RATIO", "0.08")
    monkeypatch.setattr(centroid, "make_boundary_data", lambda mesh, geom: object())

    backend = Backend()
    state = _state()
    p0 = np.asarray([0.5, -0.8])
    target = np.asarray([0.5, 0.01])

    parent = handoff._target_triangle_any(state.mesh, target)
    initial = base._predicted_exact_target_metrics(
        backend, state.mesh, int(parent), target
    )
    assert initial["predicted_min_child_area_ratio"] < 0.08

    records = []
    for level in (1, 2, 3):
        state, record = centroid._centroid_split_candidate(
            backend, state, p0, target, 0, level
        )
        assert state is not None, record
        records.append(record)
        assert record["min_triangle_quality"] >= 0.035
        assert record["min_immediate_child_area_ratio"] >= 0.08
        assert record["target_fan_shared_node_count_after"] >= 2
        assert state.parent_to_root.tolist().count(1) >= 3
        if record["predicted_exact_target_pass"]:
            break

    assert len(records) == 2
    assert all(
        row["refinement_kind"] == "atomic_target_parent_centroid_steiner"
        for row in records
    )
    assert records[-1]["predicted_min_child_area_ratio"] >= 0.08
    assert records[-1]["predicted_min_triangle_quality"] >= 0.035


def test_centroid_fallback_activates_after_single_edge_family_exhaustion(monkeypatch):
    monkeypatch.setattr(centroid, "make_boundary_data", lambda mesh, geom: object())
    backend = Backend()
    state = _state()
    p0 = np.asarray([0.5, -0.8])
    target = np.asarray([0.5, 0.01])

    def exhausted(*args, **kwargs):
        return None, {
            "level": 1,
            "accepted": False,
            "reason": "no_target_tip_triangle",
            "target_cavity_failure": {
                "reason": "no_quality_safe_adjacent_target_cavity_bisection"
            },
        }

    monkeypatch.setattr(centroid, "_ORIGINAL_REFINE_ONCE", exhausted)
    refined, record = centroid.refine_once(
        backend, state, p0, target, 0, 1
    )
    assert refined is not None, record
    assert record["refinement_kind"] == "atomic_target_parent_centroid_steiner"
    assert record["min_immediate_child_area_ratio"] == np.testing.assert_allclose(
        record["min_immediate_child_area_ratio"], 1.0 / 3.0
    )


def test_v10051838_entry_renames_and_annotates_outputs(monkeypatch, tmp_path):
    observed = {}

    def fake_edge(argv):
        observed["install"] = entry._local37.install_ray_handoff
        for old_name in entry._RENAMES:
            payload = {"schema": "v10.0.5.18.3.7_test"}
            if old_name == entry._local37.PRODUCTION_MANIFEST:
                payload["physics_contract"] = {}
            (tmp_path / old_name).write_text(json.dumps(payload))
        return "ok"

    monkeypatch.setattr(entry._edge37, "main", fake_edge)
    result = entry.main(["--out", str(tmp_path)])
    assert result == "ok"
    assert observed["install"] is entry.install_centroid_cavity

    manifest = json.loads((tmp_path / entry.PRODUCTION_MANIFEST).read_text())
    physics = manifest["physics_contract"]
    assert manifest["point_release"] == "10.0.5.18.3.8"
    assert physics["centroid_steiner_cavity_refinement_active"] is True
    assert physics["centroid_split_nominal_child_parent_area_ratio"] == 1.0 / 3.0
    assert physics["equal_length_partition_retry"] is False
    assert physics["triangle_quality_floor_relaxed"] is False
    assert not (tmp_path / entry._local37.PRODUCTION_MANIFEST).exists()


def test_runner_and_package_contracts():
    root = Path(__file__).resolve().parents[1]
    runner = (root / "run_v10_0_5_18_3_8_four_class_centroid_cavity.sh").read_text()
    assert "mode_i_first_passage_v10_0_5_18_3_8_centroid_cavity" in runner
    assert "persistent_site_production_manifest_v10_0_5_18_3_8.json" in runner
    assert "ARRHENIUS_QUALITY_AWARE_PARTITIONS" not in runner
    pyproject = (root / "pyproject.toml").read_text()
    assert 'version = "10.0.5.18.3.8"' in pyproject
