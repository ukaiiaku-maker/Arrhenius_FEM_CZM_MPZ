import json
from pathlib import Path


OUT = Path("analysis_outputs/oneD_v2_capability_v2")


def load(name):
    return json.loads((OUT / name).read_text())


def test_rejected_trials_restore_exact_snapshot():
    audit = load("oneD_v2_pf_subdivision_fixture_audit.json")
    assert audit["rejected_trials"] > 0
    assert audit["all_rejected_trials_restore_snapshot_fingerprint"] is True


def test_exactly_once_event_boundaries():
    audit = load("oneD_v2_pf_subdivision_fixture_audit.json")
    assert audit["first_passage_accept_count"] == 1
    assert audit["geometry_commit_count"] == 1
    assert audit["post_event_renewal_count"] == 1


def test_nonexistent_rng_operations_are_not_applicable():
    audit = load("oneD_v2_pf_subdivision_fixture_audit.json")
    assert audit["driver_rng_restore"].startswith("NOT_APPLICABLE")
    assert audit["event_length_draw"].startswith("NOT_APPLICABLE")
    assert audit["direction_draw"].startswith("NOT_APPLICABLE")


def test_observer_is_physically_neutral():
    audit = load("oneD_v2_pf_subdivision_fixture_audit.json")
    assert audit["diagnostics_off_on_physical_state_identical"] is True
    assert all(audit["physical_artifact_identity"].values())


def test_pf_normal_path_independent_of_unsupported_rollback():
    audit = load("oneD_v2_pf_subdivision_fixture_audit.json")
    assert audit["classification"] == "PF_NORMAL_PATH_QUALIFIED_FAIL_CLOSED_ON_VETO"
    assert audit["rollback_and_continue"] == "UNSUPPORTED_BY_PRODUCTION"


def test_corrected_fem_transaction_lineage_is_qualified():
    audit = load("oneD_v2_fem_transaction_lineage.json")
    assert audit["qualified_commit_is_ancestor"] is True
    assert audit["lineage_status"] == "QUALIFIED_FROM_EXISTING_TRANSACTIONAL_EVIDENCE"
    allowed = {"IDENTICAL", "NONSEMANTIC_CHANGE", "SEMANTIC_CHANGE_REQUALIFIED"}
    assert {row["change_classification"] for row in audit["functions"]} <= allowed


def test_campaign_stays_closed_while_maps_are_missing():
    gate = load("oneD_v2_forward_gate_v3.json")
    assert gate["four_1000K_baselines_authorized"] is False
    assert gate["campaign_outputs_populated"] is False
    assert len(gate["blocking_gates"]) == 3
