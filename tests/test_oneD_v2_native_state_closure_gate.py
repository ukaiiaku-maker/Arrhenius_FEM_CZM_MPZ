from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"analysis_outputs/oneD_v2_native_state_closure"


def rows(name):
    with (OUT/name).open(newline="") as stream: return list(csv.DictReader(stream))


def decision(): return json.loads((OUT/"oneD_v2_level3_decision_v5.json").read_text())


def test_opening_is_explicit_not_inferred():
    data=rows("oneD_v2_state_domain_shadow.csv")
    assert data and {r["opening_source"] for r in data}=={"EXPLICIT_FIXTURE_INPUT"}


def test_exact_field_snapshots_are_candidate_independent():
    data=rows("oneD_v2_state_domain_shadow.csv")
    for backend in ("PF","FEMCZM"):
        subset=[r for r in data if r["backend"]==backend]
        assert len({r["field_snapshot_hash"] for r in subset})==1
        assert len({r["geometry_fingerprint"] for r in subset})==1


def test_state_changes_do_not_trigger_fixed_geometry_field_solve():
    for name in ("pf_native_state_shadow.json","fem_native_state_shadow.json"):
        payload=json.loads((OUT/name).read_text())
        assert payload["field_solves"]==1
        assert len(payload["records"])==36


def test_current_radius_is_serialized_with_each_probe_state():
    data=rows("oneD_v2_state_domain_shadow.csv")
    assert all(float(r["tip_radius_m"])>0 for r in data)
    assert all(float(r["front_width_m"])>0 for r in data)
    script=(ROOT/"scripts/qualify_oneD_v2_native_state_closure.py").read_text()
    assert "probe.evaluate_raw(snapshot, _pz_metrics(lanes.source), opening)" in script


def test_barriers_rates_use_current_local_drive_and_state():
    data=rows("oneD_v2_state_domain_shadow.csv")
    assert all(r["drive_reliable"]=="True" for r in data)
    assert all(r["barrier_rate_status"]=="EXACT_SOURCE_MATCH" for r in data)
    assert all(r["hazard_action_status"]=="EXACT_SOURCE_MATCH" for r in data)


def test_one_event_post_state_matches_source_for_four_cases():
    data=rows("oneD_v2_one_event_source_closure.csv")
    assert len(data)==4
    assert all(r["status"]=="ONE_EVENT_SOURCE_CLOSURE_PASS" for r in data)
    assert all(r["post_event_state_exact"]=="True" for r in data)


def test_renewal_occurs_exactly_once_per_one_event_fixture():
    data=rows("oneD_v2_state_domain_shadow.csv")
    events=[r for r in data if r["phase"]=="CLEAVAGE_FIRST_PASSAGE_AND_RENEWAL"]
    assert len(events)==24
    assert all(r["event_fired"]=="True" for r in events)
    assert all(int(r["event_count"])==1 for r in events)


def test_straight_path_scope_is_explicit():
    assert decision()["path_model"]=="STANDARDIZED_STRAIGHT_REDUCED_PATH"
    micro=rows("oneD_v2_level3_exact_microtrajectories_v2.csv")
    assert {r["backend"] for r in micro}=={"PF","FEMCZM"}
    assert {r["path_model"] for r in micro if r["backend"]=="FEMCZM"}=={"STANDARDIZED_STRAIGHT_REDUCED_PATH"}


def test_surrogate_remains_disabled(): assert decision()["surrogate_authorized"] is False


def test_no_forward_baseline_launches_early(): assert decision()["forward_baselines_authorized"] is False


def test_canonical_parameters_are_unchanged(): assert decision()["canonical_parameters_changed"] is False


def test_no_new_two_dimensional_run_occurred():
    d=decision();assert d["new_stochastic_2D_runs"]==d["new_PF_runs"]==d["new_FEMCZM_runs"]==0


def test_all_four_exact_microtrajectories_qualify():
    data=rows("oneD_v2_level3_exact_microtrajectories_v2_summary.csv")
    assert len(data)==4
    assert all(int(r["accepted_events"])==3 for r in data)
    assert all(r["status"]=="CONTROLLED_SOURCE_COMPOSITION_QUALIFIED" for r in data)
    assert all(r["trajectory_semantics"]=="CONTROLLED_SOURCE_COMPOSITION_NOT_FORWARD_BASELINE" for r in data)


def test_scientific_fingerprint_manifest_matches_files():
    manifest=json.loads((OUT/"oneD_v2_native_state_fingerprints.json").read_text())
    for name,expected in manifest.items():
        assert hashlib.sha256((OUT/name).read_bytes()).hexdigest()==expected


def test_level3_decision_is_qualified_but_baselines_remain_gated():
    d=decision();assert d["Level3"]=="CONTROLLED_SOURCE_COMPOSITION_QUALIFIED"
    assert d["microtrajectories_launched"]==4
    assert d["forward_baselines_authorized"] is False
