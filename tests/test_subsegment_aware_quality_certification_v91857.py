from __future__ import annotations

import json
from pathlib import Path

import pytest

from audit_v91857_subsegment_quality import certify_case


def _write_case(tmp_path: Path) -> Path:
    case = tmp_path / "case"
    czm = case / "czm_0700K"
    czm.mkdir(parents=True)

    accepted = []
    advance = []
    specs = {1: [2.0, 3.0], 2: [1.0, 1.5, 2.5]}
    for event_id, lengths_um in specs.items():
        n = len(lengths_um)
        for index, length_um in enumerate(lengths_um):
            accepted.append({
                "accepted": True,
                "min_triangle_quality": 0.25,
                "triangle_quality_floor": 0.035,
                "min_child_area_ratio": 0.20,
                "child_area_ratio_floor": 0.08,
            })
            row = {
                "physical_event_id": event_id,
                "physical_subsegment_index": index,
                "physical_subsegment_count": n,
                "length_m": length_um * 1.0e-6,
                "v916_committed": True,
            }
            if index == n - 1:
                row.update({
                    "v91856_quality_gate_passed": True,
                    "v91856_min_triangle_quality": 0.25,
                    "v91856_min_child_area_ratio": 0.20,
                    "v91856_active_tip_h_over_da": 0.6,
                    "v91856_resolution_warning": False,
                })
            advance.append(row)

    (case / "explicit_quality_wrapper_chain_v91856.json").write_text(json.dumps({
        "run_completed_without_exception": True,
        "accepted_events": accepted,
        "resolution_warnings": [],
        "quality_vetoes": [],
        "consecutive_veto_abort": None,
    }))
    (czm / "czm_advance_log.json").write_text(json.dumps(advance))
    return case


def test_certifies_physical_events_and_recursive_subsegments(tmp_path):
    case = _write_case(tmp_path)
    result = certify_case(case, target_um=10.0, da_um=5.0)
    assert result["certified"] is True
    assert result["physical_event_count"] == 2
    assert result["committed_subsegment_count"] == 5
    assert result["accepted_quality_transaction_count"] == 5
    assert result["final_physical_event_quality_marker_count"] == 2
    assert (case / "subsegment_aware_quality_certification_v91857.json").exists()


def test_rejects_quality_transaction_count_mismatch(tmp_path):
    case = _write_case(tmp_path)
    path = case / "explicit_quality_wrapper_chain_v91856.json"
    audit = json.loads(path.read_text())
    audit["accepted_events"].pop()
    path.write_text(json.dumps(audit))
    with pytest.raises(RuntimeError, match="quality/subsegment count mismatch"):
        certify_case(case, target_um=10.0, da_um=5.0)


def test_rejects_bad_physical_event_length(tmp_path):
    case = _write_case(tmp_path)
    path = case / "czm_0700K" / "czm_advance_log.json"
    rows = json.loads(path.read_text())
    rows[0]["length_m"] *= 0.5
    path.write_text(json.dumps(rows))
    with pytest.raises(RuntimeError, match="does not equal physical da"):
        certify_case(case, target_um=10.0, da_um=5.0)


def test_rejects_missing_final_event_quality_marker(tmp_path):
    case = _write_case(tmp_path)
    path = case / "czm_0700K" / "czm_advance_log.json"
    rows = json.loads(path.read_text())
    rows[1].pop("v91856_quality_gate_passed")
    path.write_text(json.dumps(rows))
    with pytest.raises(RuntimeError, match="lacks v9.18.5.6 quality marker"):
        certify_case(case, target_um=10.0, da_um=5.0)
