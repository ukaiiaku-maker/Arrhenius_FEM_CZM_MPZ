"""Fail-closed semantic tests for the DBTT Taylor/Peierls paper audit."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_outputs" / "taylor_peierls_spatial_coupling_paper_audit"
EXPECTED = {
    "v913_zeroD_sobol_0202500",
    "oneD_v2_dbtt_TP_f07fe99faea93479",
    "oneD_v2_dbtt_TP_4895f9e5b44deea5",
    "oneD_v2_dbtt_TP_7e668ee637fc3ac5",
    "oneD_v2_dbtt_TP_b60d78111740b058",
    "oneD_v2_dbtt_TP_bd5be1610f6e1bce",
    "oneD_v2_dbtt_TP_9555ff54d637c974",
    "oneD_v2_dbtt_TP_f2817e7998cb7be6",
}
TP_FIELDS = {
    "peierls_H0_eV", "peierls_activation_entropy_kB", "peierls_exp_a", "peierls_exp_n",
    "taylor_H0_eV", "taylor_activation_entropy_kB", "taylor_exp_a", "taylor_exp_n",
    "taylor_corr_rho_c_m2", "taylor_corr_scale",
}
DIAGNOSTIC = "COUNTERFACTUAL_DIAGNOSTIC_NOT_PRODUCTION_PHYSICS"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_01_exact_eight_ids_and_material_hashes_are_preserved():
    table = pd.read_csv(OUT / "PF_TP_CASE_PROVENANCE_TABLE.csv")
    assert set(table.candidate_id) == EXPECTED and len(table) == 8
    assert table.full_material_sha256.str.fullmatch(r"[0-9a-f]{64}").all()


def test_02_cleavage_and_emission_barriers_are_identical():
    table = pd.read_csv(OUT / "PF_TP_CASE_PROVENANCE_TABLE.csv")
    assert table.cleavage_barrier_sha256.nunique() == 1
    assert table.emission_barrier_sha256.nunique() == 1


def test_03_only_taylor_peierls_material_coordinates_vary():
    bank = pd.read_csv(OUT / "taylor_peierls_microstructure_option_bank.csv")
    vectors = [json.loads(value) for value in bank.full_material_vector_json]
    fields = set(vectors[0])
    assert all(set(vector) == fields for vector in vectors)
    varying = {field for field in fields if len({vector[field] for vector in vectors}) > 1}
    assert varying == TP_FIELDS


def test_04_original_and_pf_observed_roles_are_distinct_provenance_fields():
    roles = pd.read_csv(OUT / "pf_tp_observed_microstructure_roles.csv")
    assert {"original_1d_selection_role", "pf_observed_microstructure_role"} <= set(roles)
    assert roles.original_1d_selection_role.notna().all()
    assert roles.pf_observed_microstructure_role.notna().all()
    assert not (roles.original_1d_selection_role == roles.pf_observed_microstructure_role).any()


def test_05_lab_and_tip_relative_ledgers_use_the_same_counts_and_coordinates():
    tip = pd.read_parquet(OUT / "pf_tp_tip_relative_state_ledger.parquet")
    lab = pd.read_parquet(OUT / "pf_tp_lab_frame_state_ledger.parquet")
    assert len(tip) == len(lab)
    assert np.array_equal(tip.mobile_count.to_numpy(), lab.mobile_count.to_numpy())
    assert np.array_equal(tip.retained_count.to_numpy(), lab.retained_count.to_numpy())
    np.testing.assert_allclose(
        lab.laboratory_position_m - lab.current_tip_laboratory_m,
        tip.tip_relative_position_m,
        rtol=0.0, atol=2e-19,
    )
    assert (tip.loc[tip.region.eq("WAKE_BEHIND"), "tip_relative_position_m"] < 0).all()


def test_06_event_translation_conserves_or_explicitly_accounts_for_state():
    table = pd.read_csv(OUT / "pf_tp_event_state_translation.csv")
    assert len(table) == 472
    assert table.state_conservation_residual.abs().max() <= 1e-8
    assert table.operator_mobile_conservation_residual.abs().max() <= 1e-8
    assert table.operator_retained_conservation_residual.abs().max() <= 1e-8
    assert np.isfinite(table.inferred_wake_discard_or_unobserved_removal).all()


def test_07_frozen_swaps_leave_authoritative_production_files_unchanged():
    table = pd.read_csv(OUT / "PF_TP_CASE_PROVENANCE_TABLE.csv")
    for row in table.itertuples(index=False):
        root = Path(row.run_directory)
        assert sha(root / "steps_1100K.csv") == row.steps_sha256
        assert sha(root / "anisotropic_emission_audit_v10174.json") == row.observer_sha256
        assert sha(root / "stochastic_avalanche_geometry_events.json") == row.geometry_events_sha256


def test_08_every_counterfactual_is_labelled_nonproduction():
    for name in ("pf_tp_state_swap_matrix.csv", "pf_tp_geometry_swap_matrix.csv",
                 "pf_tp_component_ablation_matrix.csv", "pf_tp_resistance_component_decomposition.csv"):
        table = pd.read_csv(OUT / name)
        assert set(table.diagnostic_semantics) == {DIAGNOSTIC}


def test_09_native_and_local_resistance_are_not_conflated():
    table = pd.read_csv(OUT / "pf_tp_resistance_component_decomposition.csv")
    native = table[table.component.eq("NATIVE_APPARENT_DELTA_K_REINIT")].value
    local = table[table.component.eq("LOCAL_SOURCE_OPENING_EQUIVALENT_DELTA_AT_REFERENCE_1UM")].value
    assert len(native) == len(local) == 8
    assert (native < 0).all() and (local > 0).all()


def test_10_barrier_classification_uses_source_formulas_not_raw_h0():
    actual = pd.read_parquet(OUT / "pf_tp_actual_state_barriers_and_rates.parquet")
    classification = pd.read_csv(OUT / "pf_tp_candidate_plausibility_classification.csv")
    assert actual.source_formula.str.contains(r"ExpFloorBarrier|TransportBarrier|_transport_rates|gammainc").any()
    assert classification.classification_uses_source_evaluated_barriers_not_raw_H0.all()


def test_11_emission_relative_ratios_use_the_same_physical_state():
    ratios = pd.read_parquet(OUT / "pf_tp_source_relevant_barrier_rate_ratios.parquet")
    assert ratios.same_physical_state.all()
    assert set(ratios.temperature_K.unique()) == set(np.arange(300.0, 1200.0 + 1e-9, 50.0))
    required = {"deltaG_peierls_over_emission", "deltaG_taylor_over_emission",
                "log10_rate_peierls_over_emission", "log10_rate_encounter_over_emission",
                "log10_rate_taylor_completion_over_encounter"}
    assert required <= set(ratios) and np.isfinite(ratios[list(required)].to_numpy()).all()


def test_12_floor_saturation_inactivity_and_stress_domains_are_reported():
    table = pd.read_csv(OUT / "pf_tp_candidate_plausibility_classification.csv")
    for prefix in ("peierls", "taylor_completion_single_hit", "taylor_completion_multiplicity", "encounter_retention"):
        for suffix in ("floor_fraction", "saturated_fraction", "inactive_fraction",
                       "loaded_floor_fraction", "loaded_saturated_fraction", "loaded_inactive_fraction",
                       "observed_stress_min_Pa", "observed_stress_max_Pa", "observed_regime"):
            assert f"{prefix}_{suffix}" in table


def test_13_every_production_figure_has_three_formats_and_source_data():
    figure_names = {path.stem for path in (OUT / "figures").glob("*.pdf")}
    assert len(figure_names) == 11
    for name in figure_names:
        for suffix in (".pdf", ".svg", ".png"):
            assert (OUT / "figures" / f"{name}{suffix}").stat().st_size > 1000
        assert (OUT / "figure_source_data" / f"{name}_source.parquet").is_file()


def test_14_figure_labels_do_not_call_native_trajectory_an_rcurve():
    for svg in (OUT / "figures").glob("*.svg"):
        assert "r-curve" not in svg.read_text(errors="ignore").lower()


def test_15_every_paper_claim_resolves_to_scoped_evidence():
    claims = pd.read_csv(OUT / "paper_claims_and_evidence.csv")
    assert claims[["claim", "status", "supporting_artifact", "supporting_rows_checkpoints",
                   "scope", "limitation", "allowed_manuscript_wording", "wording_to_avoid"]].notna().all().all()
    for artifact in claims.supporting_artifact:
        assert (OUT / artifact).exists() or (ROOT / "analysis_outputs" / "oneD_v2_taylor_peierls_spatial_transfer" / artifact).exists()


def test_16_option_registry_contains_material_coordinates_only():
    bank = pd.read_csv(OUT / "taylor_peierls_microstructure_option_bank.csv")
    forbidden = {"L_pz_um_recommended", "n_bins_recommended", "mesh", "wake", "observer",
                 "numerical_tolerance", "cache", "avalanche_grouping"}
    for encoded in bank.full_material_vector_json:
        assert not (set(json.loads(encoded)) & forbidden)


def test_17_every_option_is_explicitly_not_fatigue_evaluated():
    bank = pd.read_csv(OUT / "taylor_peierls_microstructure_option_bank.csv")
    assert not bank.fatigue_evaluated.any()
    assert set(bank.fatigue_validation_status) == {"NOT_EVALUATED"}


def test_18_provenance_certifies_no_fatigue_repository_or_run_change():
    payload = json.loads((OUT / "PF_TP_PAPER_AUDIT_PROVENANCE.json").read_text())
    assert payload["fatigue_code_modified"] is False
    assert payload["fatigue_calculations_run"] is False
    assert payload["new_femczm_run_count"] == 0


def test_19_provenance_certifies_no_new_stochastic_pf_trajectory():
    payload = json.loads((OUT / "PF_TP_PAPER_AUDIT_PROVENANCE.json").read_text())
    assert payload["new_stochastic_pf_trajectory_count"] == 0
    assert payload["production_trajectory_changed"] is False
    assert payload["new_memory_law"] is False
    assert payload["cleavage_barrier_changed"] is False
    assert payload["emission_barrier_changed"] is False


def test_20_scientific_fingerprint_matches_generated_bundle():
    payload = json.loads((OUT / "PF_TP_PAPER_AUDIT_PROVENANCE.json").read_text())
    digest = hashlib.sha256()
    paths = [OUT / relative for relative in payload["scientific_fingerprint_artifacts"]]
    for path in sorted(paths, key=lambda value: value.name):
        digest.update(path.name.encode())
        digest.update(path.read_bytes())
    assert digest.hexdigest() == payload["scientific_fingerprint_sha256"]


def test_producer_commit_is_the_last_commit_that_changed_the_generator():
    payload = json.loads((OUT / "PF_TP_PAPER_AUDIT_PROVENANCE.json").read_text())
    expected = subprocess.check_output(
        ["git", "log", "-1", "--format=%H", "--", "scripts/run_pf_tp_spatial_coupling_paper_audit.py"],
        cwd=ROOT, text=True,
    ).strip()
    assert payload["producer_code_commit"] == expected


def test_complete_factorial_has_all_pair_interactions_and_higher_residual():
    ablations = pd.read_csv(OUT / "pf_tp_component_ablation_matrix.csv")
    assert len(ablations[ablations.factorial_member.eq(True)]) == 8 * 2 * 16
    decomposition = pd.read_csv(OUT / "pf_tp_resistance_component_decomposition.csv")
    assert decomposition.component.str.endswith("PAIR_INTERACTION").sum() == 8 * 6
    assert decomposition.component.eq("NONADDITIVE_HIGHER_ORDER_INTERACTION_RESIDUAL").sum() == 8


def test_fail_closed_option_tiers_match_manifest():
    bank = pd.read_csv(OUT / "taylor_peierls_microstructure_option_bank.csv")
    manifest = json.loads((OUT / "taylor_peierls_microstructure_option_manifest.json").read_text())
    assert bank.complete_mechanism_set.all() and len(bank) == 8
    assert bank.physically_interpretable_shortlist.sum() == len(manifest["physically_interpretable_shortlist_candidate_ids"])
    assert bank.minimal_future_test_set.sum() == len(manifest["minimal_future_test_set_candidate_ids"])
