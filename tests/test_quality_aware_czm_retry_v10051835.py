from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from arrhenius_fracture import mode_i_first_passage_v9_18_5_6 as v91856
from arrhenius_fracture import (
    mode_i_first_passage_v10_0_5_18_3_5_four_class_quality_aware_czm_retry
    as entry
)
from arrhenius_fracture import quality_aware_czm_retry_v10051835 as quality
from arrhenius_fracture import quality_aware_czm_patch_v10051835 as patch


class Mesh:
    def __init__(self, nodes, elems, hbar_tip=1.0):
        self.nodes = np.asarray(nodes, dtype=float)
        self.elems = np.asarray(elems, dtype=int)
        self.nn = len(self.nodes)
        self.ne = len(self.elems)
        tri = self.nodes[self.elems]
        self.area_e = 0.5 * np.abs(
            (tri[:, 1, 0] - tri[:, 0, 0])
            * (tri[:, 2, 1] - tri[:, 0, 1])
            - (tri[:, 1, 1] - tri[:, 0, 1])
            * (tri[:, 2, 0] - tri[:, 0, 0])
        )
        self.hbar = float(hbar_tip)
        self.hbar_tip = float(hbar_tip)


class RefineBackend:
    def __init__(self):
        self.geom = SimpleNamespace(Lx=10.0, Ly=10.0)
        self.min_triangle_quality = 0.035
        self.min_area_ratio = 0.08

    @staticmethod
    def _tip_geometric_node_ids(mesh, p0, front_id):
        distance = np.linalg.norm(mesh.nodes - np.asarray(p0)[None, :], axis=1)
        dmin = float(np.min(distance))
        return np.where(distance <= dmin + 1.0e-12)[0].astype(int).tolist()

    @staticmethod
    def _incident_elements(elems, node_ids):
        return np.where(
            np.any(np.isin(np.asarray(elems), list(node_ids)), axis=1)
        )[0]

    @staticmethod
    def _triangle_quality(nodes, elems):
        tri = np.asarray(nodes)[np.asarray(elems)]
        l2 = (
            np.sum((tri[:, 1] - tri[:, 0]) ** 2, axis=1)
            + np.sum((tri[:, 2] - tri[:, 1]) ** 2, axis=1)
            + np.sum((tri[:, 0] - tri[:, 2]) ** 2, axis=1)
        )
        a2 = np.abs(
            (tri[:, 1, 0] - tri[:, 0, 0])
            * (tri[:, 2, 1] - tri[:, 0, 1])
            - (tri[:, 1, 1] - tri[:, 0, 1])
            * (tri[:, 2, 0] - tri[:, 0, 0])
        )
        return 2.0 * np.sqrt(3.0) * a2 / np.maximum(l2, 1.0e-300)

    @staticmethod
    def _insert_point_on_edge(mesh, displacement, q, edge_i, edge_j):
        q = np.asarray(q, dtype=float)
        nodes = np.vstack([mesh.nodes, q[None, :]])
        new_id = len(mesh.nodes)
        parents = np.where(
            np.any(mesh.elems == int(edge_i), axis=1)
            & np.any(mesh.elems == int(edge_j), axis=1)
        )[0]
        elems = mesh.elems.copy()
        appended = []
        appended_parent = []
        for parent in parents:
            conn = [int(value) for value in mesh.elems[int(parent)]]
            third = [value for value in conn if value not in (edge_i, edge_j)][0]
            elems[int(parent)] = [third, int(edge_i), new_id]
            appended.append([third, new_id, int(edge_j)])
            appended_parent.append(int(parent))
        elems = np.vstack([elems, np.asarray(appended, dtype=int)])
        refined = Mesh(nodes, elems, hbar_tip=mesh.hbar_tip)
        parent_map = np.concatenate(
            [
                np.arange(mesh.ne, dtype=int),
                np.asarray(appended_parent, dtype=int),
            ]
        )
        u = np.asarray(displacement, dtype=float).reshape(-1, 2)
        uq = 0.5 * (u[int(edge_i)] + u[int(edge_j)])
        u = np.vstack([u, uq[None, :]]).reshape(-1)
        return (
            refined,
            u,
            "ok",
            {"n_new_bulk_elements": len(appended)},
            parent_map,
        )


def _state(mesh):
    return quality._State(
        mesh=mesh,
        boundary=object(),
        damage=np.zeros(mesh.nn),
        displacement=np.zeros(2 * mesh.nn),
        parent_to_root=np.arange(mesh.ne, dtype=int),
    )


def _predicted(backend, state, p0, target):
    elem = patch._target_tip_triangle(
        backend, state.mesh, p0, target, front_id=0
    )
    assert elem is not None
    return patch._predicted_exact_target_metrics(
        backend, state.mesh, elem, target
    )


def test_archived_step67_candidate_is_rejected_without_relaxing_floors():
    issues = quality._floor_issues(
        0.02108230,
        0.01727567,
        0.035,
        0.08,
    )
    assert issues == [
        "triangle_quality=2.108230e-02<3.500000e-02",
        "child_area_ratio=1.727567e-02<8.000000e-02",
    ]


def test_endpoint_centered_refinement_splits_a_tip_radial_edge(monkeypatch):
    monkeypatch.setattr(patch, "make_boundary_data", lambda mesh, geom: object())
    monkeypatch.setenv("ARRHENIUS_MIN_ACCEPTED_TRIANGLE_QUALITY", "0.035")
    monkeypatch.setenv("ARRHENIUS_MIN_ACCEPTED_CHILD_AREA_RATIO", "0.08")

    mesh = Mesh(
        [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
        [[0, 1, 2]],
    )
    refined, record = quality._refine_once(
        RefineBackend(),
        _state(mesh),
        p0=np.array([0.0, 0.0]),
        target=np.array([0.1, 0.1]),
        front_id=0,
        level=1,
    )

    assert refined is not None
    assert record["accepted"] is True
    assert record["refinement_kind"] == (
        "endpoint_centered_tip_radial_midpoint"
    )
    assert record["tip_node"] in (record["edge_i"], record["edge_j"])
    assert {record["edge_i"], record["edge_j"]} != {1, 2}
    assert record["min_immediate_child_area_ratio"] == pytest.approx(0.5)
    assert record["min_triangle_quality"] >= 0.035
    assert record["predicted_limiting_margin"] > 0.0
    assert refined.mesh.ne == 2
    assert refined.mesh.nn == 4
    assert refined.parent_to_root.tolist() == [0, 0]


def test_archived_step67_geometry_reaches_a_quality_safe_endpoint(monkeypatch):
    monkeypatch.setattr(patch, "make_boundary_data", lambda mesh, geom: object())
    monkeypatch.setenv("ARRHENIUS_MIN_ACCEPTED_TRIANGLE_QUALITY", "0.035")
    monkeypatch.setenv("ARRHENIUS_MIN_ACCEPTED_CHILD_AREA_RATIO", "0.08")

    # Reconstructed exactly from the archived step-67 p0/p1 diagnostic and
    # the first two opposite-edge midpoint audit records.
    p0 = np.array([0.0005861866667123705, 1.200896179212802e-05])
    target = np.array([0.0005887734388763731, 1.2527425261258408e-05])
    midpoint_1 = np.array([0.0005971874999999998, 1.0369721710082628e-05])
    midpoint_2 = np.array([0.0005987187499999999, 1.3021924509172447e-05])
    vertex_a = 2.0 * midpoint_2 - midpoint_1
    vertex_b = 2.0 * midpoint_1 - vertex_a

    backend = RefineBackend()
    state = _state(
        Mesh(
            [p0, vertex_a, vertex_b],
            [[0, 1, 2]],
            hbar_tip=1.1196217271646168e-05,
        )
    )

    initial = _predicted(backend, state, p0, target)
    assert initial["predicted_min_triangle_quality"] == pytest.approx(
        0.021082302341594004
    )
    assert initial["predicted_min_child_area_ratio"] == pytest.approx(
        0.017275666997289947
    )

    records = []
    final = initial
    for level in range(1, 9):
        state, record = quality._refine_once(
            backend,
            state,
            p0=p0,
            target=target,
            front_id=0,
            level=level,
        )
        assert state is not None, record
        records.append(record)
        assert record["min_triangle_quality"] >= 0.035
        assert record["min_immediate_child_area_ratio"] >= 0.08
        final = _predicted(backend, state, p0, target)
        if (
            final["predicted_min_triangle_quality"] >= 0.035
            and final["predicted_min_child_area_ratio"] >= 0.08
        ):
            break

    assert len(records) <= 4
    assert all(
        row["refinement_kind"]
        == "endpoint_centered_tip_radial_midpoint"
        for row in records
    )
    assert final["predicted_min_triangle_quality"] >= 0.035
    assert final["predicted_min_child_area_ratio"] >= 0.08
    assert state.parent_to_root.tolist() == [0] * state.mesh.ne


def test_quality_failure_reaches_patch_retry_before_final_veto(monkeypatch):
    quality.reset_audit()
    root_mesh = SimpleNamespace(ne=1)
    failed = SimpleNamespace(inserted=False, reason="v91856_quality_veto")
    accepted = SimpleNamespace(inserted=True)
    attempts = []

    def fake_call(
        original,
        self,
        state,
        root_kwargs,
        p0,
        p1,
        direction,
    ):
        attempts.append(len(attempts))
        return failed if len(attempts) == 1 else accepted

    monkeypatch.setattr(quality, "_call_on_state", fake_call)
    monkeypatch.setattr(
        quality,
        "_refine_once",
        lambda self, state, *args: (
            state,
            {"accepted": True, "level": args[-1]},
        ),
    )
    monkeypatch.setattr(
        quality,
        "_partition_retry",
        lambda *args, **kwargs: pytest.fail(
            "partition retry should not run after successful patch replay"
        ),
    )

    backend = SimpleNamespace(
        _transaction_snapshot=lambda: {},
        _transaction_rollback=lambda snap: None,
        advance_log=[],
    )
    quality.quality_aware_strict_advance_v10051835._original = lambda *a, **k: None
    result = quality.quality_aware_strict_advance_v10051835(
        backend,
        mesh=root_mesh,
        boundary=object(),
        damage=np.zeros(1),
        displacement=np.zeros(2),
        p0=np.array([0.0, 0.0]),
        p1=np.array([1.0, 0.0]),
        direction=np.array([1.0, 0.0]),
        front_id=0,
    )

    assert result is accepted
    assert attempts == [0, 1]


def test_entrypoint_installs_quality_wrapper_and_rewrites_v35_manifests(
    monkeypatch,
    tmp_path,
):
    original_wrapper = v91856._strict_quality_advance_v91856
    observed = {}

    def fake_base(argv):
        observed["wrapper"] = v91856._strict_quality_advance_v91856
        for name in (
            entry._base.PRODUCTION_MANIFEST,
            entry._base.SELECTION_MANIFEST,
            entry._base.TRANSFER_MANIFEST,
            entry._base.SEED_MANIFEST,
        ):
            payload = {
                "point_release": "10.0.5.18.3.4",
                "model": "old",
            }
            if name == entry._base.PRODUCTION_MANIFEST:
                payload["physics_contract"] = {}
            if name == entry._base.SELECTION_MANIFEST:
                payload["policy"] = {}
            (tmp_path / name).write_text(json.dumps(payload))
        (tmp_path / "compact_corridor_mesh_v91852.json").write_text("{}")
        return "ok"

    monkeypatch.setattr(entry._base, "main", fake_base)
    result = entry.main(["--out", str(tmp_path)])

    assert result == "ok"
    assert (
        observed["wrapper"]
        is quality.quality_aware_strict_advance_v10051835
    )
    assert v91856._strict_quality_advance_v91856 is original_wrapper

    production = json.loads(
        (tmp_path / entry.PRODUCTION_MANIFEST).read_text()
    )
    assert production["point_release"] == "10.0.5.18.3.5"
    physics = production["physics_contract"]
    assert physics["quality_gate_evaluated_inside_retry_transaction"] is True
    assert physics["triangle_quality_floor_relaxed"] is False
    assert physics["child_area_ratio_floor_relaxed"] is False
    assert not (tmp_path / entry._base.PRODUCTION_MANIFEST).exists()

    quality_audit = json.loads((tmp_path / entry.QUALITY_AUDIT).read_text())
    assert quality_audit["run_completed_without_exception"] is True
    assert quality_audit["quality_gate_inside_retry_transaction"] is True


def test_v35_runner_and_package_contracts():
    root = Path(__file__).resolve().parents[1]
    runner = (
        root / "run_v10_0_5_18_3_5_four_class_focused_long_growth.sh"
    ).read_text()
    assert (
        "mode_i_first_passage_v10_0_5_18_3_5_four_class_"
        "quality_aware_czm_retry"
    ) in runner
    assert "ARRHENIUS_QUALITY_AWARE_PATCH_LEVELS" in runner
    assert "ARRHENIUS_MIN_ACCEPTED_TRIANGLE_QUALITY" in runner
    assert "ARRHENIUS_MIN_ACCEPTED_CHILD_AREA_RATIO" in runner
    assert "shape_regular_local_tip_patch_bisection=true" in runner

    sweep = (
        root
        / "run_v10_0_5_18_3_5_four_class_14T_400um_quality_aware_sweep.sh"
    ).read_text()
    assert "quality_gate_inside_retry_transaction=true" in sweep
    assert "triangle_quality_floor_relaxed=false" in sweep
    assert "child_area_ratio_floor_relaxed=false" in sweep

    pyproject = (root / "pyproject.toml").read_text()
    assert 'version = "10.0.5.18.3.5"' in pyproject
