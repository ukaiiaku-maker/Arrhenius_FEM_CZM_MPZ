from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from arrhenius_fracture import adaptive_czm_quality_subdivision_v100518 as quality
from arrhenius_fracture.adaptive_czm_forward_support_prerefine_v100521 import (
    _immediate_metrics,
    installed_forward_support_prerefine_v100521,
)
from arrhenius_fracture.adaptive_czm_balanced_support_prerefine_v100522 import (
    installed_balanced_support_prerefine_v100522,
)
from arrhenius_fracture.crack_backend import AdaptiveCZMBackend
from arrhenius_fracture.mesh import rebuild_tri_mesh


def _backend_and_mesh():
    geom = SimpleNamespace(Lx=10.0, Ly=10.0)
    backend = AdaptiveCZMBackend(
        geom=geom,
        min_area_ratio=0.08,
        min_triangle_quality=0.035,
        max_node_move_factor=0.01,
    )
    nodes = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]], dtype=float)
    elems = np.array([[0, 1, 2]], dtype=int)
    mesh = rebuild_tri_mesh(nodes, elems, tip_centers=[np.array([0.0, 0.0])])
    return backend, mesh


def test_balanced_support_repairs_noncollinear_seven_percent_child():
    backend, mesh = _backend_and_mesh()
    p0 = np.array([0.0, 0.0])
    target = np.array([0.83, 0.07])
    displacement = np.zeros(2 * mesh.nn)

    direct = backend._insert_target_in_incident_triangle(
        mesh, displacement, p0, target, 0
    )
    direct_mesh, _, _, direct_reason, _, direct_map = direct
    assert direct_reason == "ok"
    _, direct_area = _immediate_metrics(backend, mesh, direct_mesh, direct_map)
    assert np.isclose(direct_area, 0.07)
    assert direct_area < backend.min_area_ratio

    with installed_forward_support_prerefine_v100521():
        forward_only = backend._insert_target_in_incident_triangle(
            mesh, displacement, p0, target, 0
        )
    assert forward_only[3] == "ok"

    with installed_forward_support_prerefine_v100521():
        with installed_balanced_support_prerefine_v100522():
            repaired = backend._insert_target_in_incident_triangle(
                mesh, displacement, p0, target, 0
            )

    new_mesh, _, endpoint, reason, meta, parent_map = repaired
    assert reason == "ok_balanced_support_prerefine"
    assert np.allclose(endpoint, target, rtol=0.0, atol=1e-14)
    assert meta["v100522_balanced_support_prerefine"] is True
    assert meta["exact_target_preserved"] is True
    assert meta["support_node_is_intact_off_ray"] is True
    assert meta["v100522_min_immediate_child_area_ratio"] >= 0.08
    assert meta["v100522_min_immediate_triangle_quality"] >= 0.035
    assert len(parent_map) == new_mesh.ne
    assert np.all(np.asarray(parent_map) == 0)


def test_quality_assessment_uses_balanced_support_subtransactions(monkeypatch):
    backend = SimpleNamespace(
        advance_log=[
            {
                "v100522_balanced_support_prerefine": True,
                "v100522_min_immediate_child_area_ratio": 0.11,
                "quality_prerefine_subtransaction_count": 2,
            }
        ]
    )

    def base_assess(self, old_mesh, result, kwargs, log_start):
        return (
            ["child_area_ratio=1.000000e-02<8.000000e-02"],
            {
                "min_child_area_ratio": 0.01,
                "child_area_ratio_floor": 0.08,
                "child_area_ratio_reference": "composed_origin",
                "multi_generation_atomic_event": False,
                "topology_subtransaction_count": 1,
                "issues": ["child_area_ratio=1.000000e-02<8.000000e-02"],
                "accepted": False,
            },
            False,
        )

    monkeypatch.setattr(quality, "_assess", base_assess)
    with installed_balanced_support_prerefine_v100522():
        issues, row, _ = quality._assess(backend, None, None, {}, 0)

    assert issues == []
    assert row["accepted"] is True
    assert row["min_child_area_ratio"] == 0.11
    assert row["multi_generation_atomic_event"] is True
    assert row["topology_subtransaction_count"] == 2
    assert row["child_area_ratio_reference"] == (
        "minimum_immediate_balanced_support_prerefine_subtransaction"
    )


def test_v100522_entry_preserves_endpoint_quality_and_physics_contract():
    text = open(
        "arrhenius_fracture/mode_i_first_passage_v10_0_5_22_four_class_standalone.py"
    ).read()
    assert '"exact_crack_endpoint_preserved": True' in text
    assert '"exact_stochastic_event_length_preserved": True' in text
    assert '"triangle_quality_floor_relaxed": False' in text
    assert '"child_area_ratio_floor_relaxed": False' in text
    assert '"constitutive_physics_changed_by_geometry_repair": False' in text
