from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from arrhenius_fracture import crack_backend as cb
from arrhenius_fracture import adaptive_czm_quality_subdivision_v100518 as uniform
from arrhenius_fracture import adaptive_czm_quality_partition_v100519 as adaptive
from arrhenius_fracture import adaptive_czm_partition_failure_audit_v100520 as audit


class _Mesh:
    def __init__(self):
        self.ne = 1
        self.nn = 3
        self.hbar_tip = 1.0
        self.nodes = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
        self.elems = np.array([[0, 1, 2]], dtype=int)
        self.area_e = np.array([0.5])


class _Backend:
    def __init__(self):
        self.advance_log = []
        self.cohesive_network = SimpleNamespace(elements=[])
        self._v91856_consecutive_geometry_vetoes = 0
        self._v91856_last_veto_reason = None

    def _transaction_snapshot(self):
        return {"n_log": len(self.advance_log)}

    def _transaction_rollback(self, snap):
        del self.advance_log[snap["n_log"] :]


def _always_reject(self, **kwargs):
    return cb.CrackAdvanceResult(
        mesh=kwargs["mesh"],
        boundary=kwargs["boundary"],
        damage=kwargs["damage"],
        displacement=kwargs["displacement"],
        moved=0.0,
        inserted=False,
        reason="no_quality_safe_exact_ray_move;hrefine:target_not_in_tip_incident_triangle;no_ray_exit",
    )


def test_terminal_raw_backend_reasons_survive_atomic_rollback(monkeypatch):
    backend = _Backend()
    mesh = _Mesh()

    uniform._v91856._AUDIT["quality_vetoes"] = []
    uniform._v9185._RUNTIME["quality_vetoes"] = []
    monkeypatch.setattr(uniform, "_subdivision_counts", lambda: (2,))
    monkeypatch.setenv("ARRHENIUS_QUALITY_PARTITION_MAX_DEPTH", "2")
    monkeypatch.setattr(
        uniform._v91856,
        "_record_or_raise",
        lambda self, kwargs, result: result,
    )
    uniform._quality_subdividing_advance_v100518._original = _always_reject

    with adaptive.installed_adaptive_quality_partition_v100519():
        with audit.installed_partition_failure_audit_v100520():
            result = uniform._quality_subdividing_advance_v100518(
                backend,
                mesh=mesh,
                boundary=object(),
                damage=np.zeros(mesh.nn),
                displacement=np.zeros(2 * mesh.nn),
                p0=np.array([0.0, 0.0]),
                p1=np.array([1.0, 0.0]),
                direction=np.array([1.0, 0.0]),
                front_id=0,
            )

    assert result.inserted is False
    assert backend.advance_log == []
    rows = uniform._v91856._AUDIT["quality_vetoes"]
    assert len(rows) == 1
    row = rows[0]
    assert row["schema"] == audit.MODEL_ID
    assert row["geometry_transaction_rolled_back"] is True
    assert row["physical_event_consumed"] is False
    assert row["raw_backend_rejection_count"] > 0
    assert row["dominant_raw_backend_reason"].endswith("no_ray_exit")
    assert row["raw_backend_reason_counts"][row["dominant_raw_backend_reason"]] > 0
    assert row["requested_da_m"] == 1.0
    assert row["minimum_rejected_segment_m"] > 0.0
    assert uniform._v9185._RUNTIME["quality_vetoes"][-1]["schema"] == audit.MODEL_ID


def test_v100520_wrapper_declares_diagnostic_only_contract():
    text = open(
        "arrhenius_fracture/mode_i_first_passage_v10_0_5_20_four_class_standalone.py"
    ).read()
    assert "geometry_behavior_changed_by_failure_audit\": False" in text
    assert "quality_floors_changed_by_failure_audit\": False" in text
    assert "constitutive_physics_changed_by_failure_audit\": False" in text
    assert "raw_backend_rejection_reasons_persisted_after_rollback\": True" in text
