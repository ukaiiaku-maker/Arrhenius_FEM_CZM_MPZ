import json
from pathlib import Path

import pandas as pd


OUT = Path("analysis_outputs/oneD_v2_mechanics_maps_and_baselines")


def test_lifecycle_manifest_keeps_backend_capabilities_separate():
    manifest = json.loads((OUT / "oneD_v2_backend_lifecycle_manifest.json").read_text())
    assert manifest["PF"]["threshold_RNG_fields_not_present"] == "NOT_APPLICABLE"
    assert manifest["PF"]["late_veto_rollback_and_continue"] == "UNSUPPORTED_BY_PRODUCTION"
    assert manifest["FEMCZM"]["late_veto_rollback_and_continue"] == "QUALIFIED"


def test_numerical_event_is_not_physical_avalanche():
    tx = pd.read_csv(OUT / "pf_2d_event_transactions_v2.csv")
    ava = pd.read_csv(OUT / "pf_2d_physical_avalanches_v2.csv")
    assert len(tx) == 416
    assert len(ava) == 3
    assert (ava.event_transaction_count > 1).all()


def test_pf_target_is_right_censored_not_arrested():
    ava = pd.read_csv(OUT / "pf_2d_physical_avalanches_v2.csv")
    final = ava.groupby("material_class", sort=False).tail(1)
    assert final.target_right_censored.all()
    assert "arrest" not in " ".join(ava.columns).lower()


def test_onset_candidates_are_pre_event_states():
    onset = pd.read_csv(OUT / "pf_2d_onset_candidates_v2.csv")
    assert onset.onset_role.str.endswith("PRE_EVENT").all()
    assert (onset.pre_event_projected_extension_um <= onset.post_event_projected_extension_um).all()


def test_interior_states_are_not_resistance_points():
    interior = pd.read_csv(OUT / "pf_2d_in_avalanche_drive_v2.csv")
    assert interior.semantic_role.eq("INTERIOR_PRE_EVENT_NOT_RESISTANCE_POINT").all()


def test_archived_source_fingerprints_are_present():
    tx = pd.read_csv(OUT / "pf_2d_event_transactions_v2.csv")
    assert tx.source_steps_sha256.str.fullmatch(r"[0-9a-f]{64}").all()
    assert tx.post_event_process_zone_fingerprint.str.fullmatch(r"[0-9a-f]{64}").all()
