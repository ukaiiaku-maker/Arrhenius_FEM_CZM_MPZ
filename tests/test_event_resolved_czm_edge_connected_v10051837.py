from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from arrhenius_fracture import event_resolved_czm_retry_v10051837 as retry
from arrhenius_fracture import (
    mode_i_first_passage_v10_0_5_18_3_7_event_resolved_edge_connected
    as entry
)


class Mesh:
    def __init__(self, nodes, elems):
        self.nodes = np.asarray(nodes, dtype=float)
        self.elems = np.asarray(elems, dtype=int)
        self.nn = len(self.nodes)
        self.ne = len(self.elems)
        self.hbar = 1.0
        self.hbar_tip = 1.0


class Backend:
    @staticmethod
    def _tip_geometric_node_ids(mesh, p0, front_id):
        distance = np.linalg.norm(
            mesh.nodes - np.asarray(p0, dtype=float)[None, :], axis=1
        )
        minimum = float(np.min(distance))
        return np.where(distance <= minimum + 1.0e-12)[0].astype(int).tolist()

    @staticmethod
    def _incident_elements(elems, node_ids):
        return np.where(
            np.any(np.isin(np.asarray(elems), list(node_ids)), axis=1)
        )[0]


def test_one_node_target_cavity_is_not_directly_reachable():
    mesh = Mesh(
        [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [2.0, 0.0], [2.0, 1.0]],
        [[0, 1, 2], [1, 3, 4]],
    )
    assert entry.edge_connected_local_endpoint_reachable(
        Backend(), mesh, np.array([0.0, 0.0]), np.array([1.7, 0.3]), 0
    ) is False


def test_shared_edge_target_cavity_is_directly_reachable():
    mesh = Mesh(
        [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]],
        [[0, 1, 2], [1, 3, 2]],
    )
    assert entry.edge_connected_local_endpoint_reachable(
        Backend(), mesh, np.array([0.0, 0.0]), np.array([0.7, 0.7]), 0
    ) is True


def test_target_in_tip_triangle_is_directly_reachable():
    mesh = Mesh(
        [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
        [[0, 1, 2]],
    )
    assert entry.edge_connected_local_endpoint_reachable(
        Backend(), mesh, np.array([0.0, 0.0]), np.array([0.2, 0.2]), 0
    ) is True


def test_entry_installs_and_restores_edge_connected_overlay(monkeypatch, tmp_path):
    original = retry._local_endpoint_reachable
    observed = {}

    def fake_base(argv):
        observed["installed"] = retry._local_endpoint_reachable
        for name in (
            entry._base.RETRY_AUDIT,
            entry._base.RESOLUTION_AUDIT,
            entry._base.PRODUCTION_MANIFEST,
        ):
            payload = {"physics_contract": {}} if name == entry._base.PRODUCTION_MANIFEST else {}
            (tmp_path / name).write_text(json.dumps(payload))
        return "ok"

    monkeypatch.setattr(entry._base, "main", fake_base)
    result = entry.main(["--out", str(tmp_path)])

    assert result == "ok"
    assert observed["installed"] is entry.edge_connected_local_endpoint_reachable
    assert retry._local_endpoint_reachable is original

    audit = json.loads((tmp_path / entry._base.RETRY_AUDIT).read_text())
    assert audit["direct_endpoint_refinement_requires_tip_fan_shared_edge"] is True
    assert audit["minimum_tip_fan_shared_node_count_for_direct_endpoint"] == 2
    assert audit["one_node_target_cavity_routes_to_exact_ray_crossing"] is True
    assert audit["equal_length_partition_retry"] is False

    manifest = json.loads((tmp_path / entry._base.PRODUCTION_MANIFEST).read_text())
    physics = manifest["physics_contract"]
    assert physics["direct_endpoint_refinement_requires_tip_fan_shared_edge"] is True
    assert physics["one_node_target_cavity_routes_to_exact_ray_crossing"] is True
    assert physics["hazard_or_event_length_changed"] is False
    assert physics["constitutive_physics_changed"] is False


def test_edge_connected_runner_contract():
    root = Path(__file__).resolve().parents[1]
    runner = (
        root / "run_v10_0_5_18_3_7_four_class_event_resolved_edge_connected.sh"
    ).read_text()
    assert "event_resolved_edge_connected" in runner
    assert "event_resolved_local_cavity" in runner
