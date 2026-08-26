from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_outputs" / "oneD_v2_taylor_peierls_spatial_transfer"


def test_spatial_screen_and_selected_registry_are_complete_and_fixed_field() -> None:
    candidates = pd.read_parquet(OUT / "oneD_v2_spatial_1d_candidate_summary.parquet")
    vectors = pd.read_parquet(OUT / "oneD_v2_spatial_1d_response_feature_vectors.parquet")
    selected = pd.read_csv(OUT / "oneD_v2_spatial_transfer_selected_candidates.csv")
    registry = pd.read_csv(OUT / "oneD_v2_spatial_transfer_pf_registry.csv")
    assert len(candidates) == 513
    assert candidates.status.eq("complete").all()
    assert vectors.selection_eligible.all()
    assert candidates.achieved_projected_extension_um.ge(300.0).all()
    assert vectors.first_onset_relative_error_to_control.abs().le(0.2).all()
    assert len(selected) == len(registry) == 8
    assert selected.candidate_id.is_unique and selected.selection_role.is_unique
    for field in ("cleavage_barrier_sha256", "emission_barrier_sha256"):
        assert selected[field].nunique() == 1


def test_pf_v2_outputs_are_complete_and_non_toughening() -> None:
    summary = pd.read_csv(OUT / "pf_2d_spatial_transfer_summary.csv")
    transactions = pd.read_csv(OUT / "pf_2d_spatial_transfer_event_transactions.csv")
    onsets = pd.read_csv(OUT / "pf_2d_spatial_transfer_onsets.csv")
    avalanches = pd.read_csv(OUT / "pf_2d_spatial_transfer_physical_avalanches.csv")
    assert len(summary) == 8
    assert summary.event_transaction_count.eq(59).all()
    assert summary.physical_avalanche_count.eq(2).all()
    assert summary.N_reinit.eq(1).all()
    assert summary.deltaK_reinit_MPa_sqrt_m.lt(0).all()
    assert not summary.positive_reload_separated_resistance.any()
    assert summary.target_right_censored.all()
    assert len(onsets) == len(avalanches) == 16
    assert set(onsets.onset_role) == {
        "INITIAL_ONSET_PRE_EVENT", "REINITIATION_ONSET_PRE_EVENT"
    }
    expected = transactions.groupby(["candidate_id", "physical_avalanche_index"], sort=False).head(1)
    pd.testing.assert_series_equal(
        onsets.pre_event_step.reset_index(drop=True), expected.pre_event_step.reset_index(drop=True)
    )


def test_profiles_distances_and_rank_transfer_are_complete() -> None:
    states = pd.read_parquet(OUT / "pf_2d_spatial_transfer_state_features.parquet")
    profiles = pd.read_parquet(OUT / "pf_2d_spatial_transfer_state_profiles.parquet")
    distances = pd.read_csv(OUT / "pf_2d_spatial_transfer_pairwise_distances.csv")
    transfer = pd.read_csv(OUT / "oneD_to_pf_spatial_pairwise_rank_transfer.csv")
    assert states.candidate_id.nunique() == profiles.candidate_id.nunique() == 8
    assert {"INITIAL_ONSET", "REINITIATION_ONSET_01", "CHECKPOINT_300UM"}.issubset(states.state_label)
    assert set(profiles.region) == {"ACTIVE_AHEAD_OF_TIP", "WAKE_BEHIND_TIP"}
    expected_lab = profiles.laboratory_position_m - profiles.tip_relative_position_m
    assert np.isfinite(expected_lab).all()
    assert len(distances) == len(transfer) == 28
    assert distances.distance.gt(0).all()
    assert np.isfinite(transfer[["distance_1d", "distance_2d"]]).all().all()


def test_final_decision_and_provenance_close_conditional_gate() -> None:
    decision = json.loads((OUT / "oneD_v2_taylor_peierls_spatial_transfer_final_decision.json").read_text())
    provenance = json.loads((OUT / "oneD_v2_taylor_peierls_spatial_transfer_provenance.json").read_text())
    assert decision["positive_candidate_count"] == 0
    assert not decision["conditional_local_multiseed_confirmation_launched"]
    assert decision["final_decision"] == (
        "TAYLOR_PEIERLS_REDISTRIBUTION_ALONE_INSUFFICIENT_DESPITE_2D_STATE_DIVERSITY"
    )
    assert provenance["candidate_count"] == 8
    assert provenance["fresh_pf_cases"] == 7 and provenance["reused_pf_cases"] == 1
    assert provenance["maximum_pf_workers"] == 2
    assert provenance["new_femczm_runs"] == 0
    assert not provenance["memory_term_added"]
    assert not provenance["non_taylor_peierls_material_fields_varied"]
    assert not provenance["cleavage_or_emission_barriers_changed"]
