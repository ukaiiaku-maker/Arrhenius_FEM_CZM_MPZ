from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess

import numpy as np
import pandas as pd

from reduced_fracture_v2.taylor_peierls_contract import (
    SEARCH_WHITELIST,
    barrier_hashes,
    full_material_hash,
)


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_outputs" / "oneD_v2_taylor_peierls_rcurve_search"
BASE = "3d64c918431874f1ee1c215e4441b86672a52cc6"


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_complete_population_and_material_contract() -> None:
    population = pd.read_parquet(OUT / "oneD_v2_taylor_peierls_complete_population.parquet")
    assert len(population) == 8194
    assert population.full_material_sha256.is_unique
    assert set(population.target_class) == {"Peak", "DBTT"}
    assert int(population.stage1_fully_evaluated.sum()) == 1026
    for material, local in population.groupby("target_class"):
        control = local[local.search_method.eq("HISTORICAL_CONTROL")].iloc[0]
        expected = barrier_hashes(control)
        assert local.cleavage_barrier_sha256.nunique() == 1
        assert local.emission_barrier_sha256.nunique() == 1
        assert local.cleavage_barrier_sha256.iloc[0] == expected["cleavage_barrier_sha256"]
        assert local.emission_barrier_sha256.iloc[0] == expected["emission_barrier_sha256"]
        for row in local.sample(16, random_state=7).itertuples():
            payload = json.loads(row.canonical_parameter_json)
            assert full_material_hash(payload) == row.full_material_sha256
            assert not any("lifecycle" in key or "memory" in key or "provider" in key for key in payload)


def test_pareto_bank_shortlist_and_registry_are_complete_and_distinct() -> None:
    pareto = pd.read_parquet(OUT / "oneD_v2_taylor_peierls_pareto_library.parquet")
    bank_csv = pd.read_csv(OUT / "oneD_v2_taylor_peierls_option_bank.csv")
    bank_parquet = pd.read_parquet(OUT / "oneD_v2_taylor_peierls_option_bank.parquet")
    shortlist = pd.read_csv(OUT / "oneD_v2_taylor_peierls_shortlist.csv")
    registry = pd.read_csv(OUT / "oneD_v2_taylor_peierls_material_registry.csv")
    assert len(pareto) == 1026 and pareto.pareto_nondominated.all()
    assert len(bank_csv) == len(bank_parquet) == len(registry) == 32
    assert len(shortlist) == 12
    assert set(bank_csv.candidate_id) == set(bank_parquet.candidate_id) == set(registry.candidate_id)
    assert set(shortlist.candidate_id) < set(bank_csv.candidate_id)
    assert "provider" not in registry.columns
    assert set(SEARCH_WHITELIST) <= set(registry.columns)
    for control in ("v913_zeroD_sobol_0242980", "v913_zeroD_sobol_0202500"):
        assert control in set(pareto.candidate_id)
        assert control in set(bank_csv.candidate_id)
        assert control in set(shortlist.candidate_id)


def test_long_multiseed_has_no_clipping_or_extrapolation() -> None:
    long = pd.read_parquet(OUT / "oneD_v2_taylor_peierls_long_validation.parquet")
    multi = pd.read_parquet(OUT / "oneD_v2_taylor_peierls_multiseed_results.parquet")
    assert len(long) == 858 and len(multi) == 468
    assert set(long.status) == {"TARGET_RIGHT_CENSORED"}
    assert (long.map_oracle_fallback_count == 0).all()
    assert set(multi.hazard_seed.groupby(multi.target_class).nunique()) == {3}
    for keys, local in multi.groupby(
        ["target_class", "provider", "temperature_K", "hazard_seed"]
    ):
        assert float(local.physical_avalanche_count.max() - local.physical_avalanche_count.min()) == 0.0


def test_exact_oracle_and_fem_read_only_policy() -> None:
    exact = pd.read_csv(OUT / "oneD_v2_taylor_peierls_exact_finalists.csv")
    assert exact.exact_probe_reliable.all()
    assert exact[exact.provider.eq("PF")].exact_check_status.eq(
        "PF_EXACT_FIELD_AND_TENSOR_PROBE_PASS"
    ).all()
    assert exact[exact.provider.eq("FEMCZM")].exact_check_status.eq(
        "PRECOMPUTED_EXACT_MAP_BOUNDED_NO_NEW_FEM_SOLVE"
    ).all()
    assert not exact.new_FEMCZM_solve.any()


def test_pf_onsets_are_reload_separated_and_in_avalanche_drive_is_not_resistance() -> None:
    transactions = pd.read_csv(OUT / "pf_2d_taylor_peierls_event_transactions.csv")
    onsets = pd.read_csv(OUT / "pf_2d_taylor_peierls_onset_candidates.csv")
    assert set(transactions.in_avalanche_drive_interpretation) == {
        "PF_MODEL_NATIVE_DRIVING_TRAJECTORY_NOT_RESISTANCE"
    }
    assert set(onsets.resistance_candidate_policy) == {"RELOAD_SEPARATED_PRE_EVENT_ONLY"}
    assert set(onsets.onset_role) <= {
        "INITIAL_ONSET_PRE_EVENT", "REINITIATION_ONSET_PRE_EVENT"
    }
    initial = onsets[onsets.onset_role.eq("INITIAL_ONSET_PRE_EVENT")]
    assert (initial.physical_avalanche_index == 0).all()
    reinit = onsets[onsets.onset_role.eq("REINITIATION_ONSET_PRE_EVENT")]
    assert (reinit.physical_avalanche_index > 0).all()
    assert onsets.groupby(
        ["candidate_id", "temperature_K", "hazard_seed", "target_um", "physical_avalanche_index"]
    ).size().eq(1).all()


def test_pf_profiles_are_actual_spatial_state_with_valid_provenance() -> None:
    profiles_path = OUT / "pf_2d_taylor_peierls_state_profiles.parquet"
    profiles = pd.read_parquet(profiles_path)
    assert len(profiles) == 11520
    assert set(profiles.region) == {"ACTIVE_AHEAD_OF_TIP", "WAKE_BEHIND_TIP"}
    assert set(profiles.profile_source) == {
        "EXISTING_PF_MPZ_STATE_DEFAULT_OFF_OBSERVER_NO_FEEDBACK"
    }
    active = profiles[profiles.region.eq("ACTIVE_AHEAD_OF_TIP")]
    assert np.isfinite(active.distance_from_tip_m).all()
    assert np.isfinite(active.peierls_rate_s_by_bin).all()
    assert np.isfinite(active.taylor_completion_rate_s_by_bin).all()
    assert np.isfinite(active.encounter_rate_s_by_bin).all()
    assert (active.mobile_count >= 0).all() and (active.retained_count >= 0).all()
    manifest = json.loads((OUT / "oneD_v2_taylor_peierls_provenance_manifest.json").read_text())
    assert manifest["artifacts"][profiles_path.name] == _sha(profiles_path)
    assert manifest["maximum_concurrent_heavy_PF_workers"] == 2
    assert manifest["new_FEMCZM_run_count"] == 0


def test_pf_transfer_rejects_both_variants_and_keeps_apparent_local_metrics_distinct() -> None:
    summary = pd.read_csv(OUT / "pf_2d_taylor_peierls_transfer_summary.csv")
    assert len(summary) == 16
    assert not summary.positive_reload_separated_resistance.any()
    assert not summary.direct_PF_success.any()
    assert set(summary.direct_PF_decision) == {"REJECT_NO_POSITIVE_RELOAD_SEPARATED_RESISTANCE"}
    onsets = pd.read_csv(OUT / "pf_2d_taylor_peierls_onset_candidates.csv")
    assert "pre_event_native_KJ_MPa_sqrt_m" in onsets
    assert "pre_event_K_effective_local_equivalent_MPa_sqrt_m" in onsets
    assert not np.allclose(
        onsets.pre_event_native_KJ_MPa_sqrt_m,
        onsets.pre_event_K_effective_local_equivalent_MPa_sqrt_m,
    )


def test_weak_and_ceramic_regression_and_final_decision() -> None:
    regression = pd.read_csv(OUT / "oneD_v2_weak_ceramic_compact_regression.csv")
    assert len(regression) == 12
    assert "REGRESSION_FAILURE" not in set(regression.regression_status)
    decision = json.loads((OUT / "oneD_v2_taylor_peierls_final_decision.json").read_text())
    assert decision["Peak"]["decision"] == "RETAIN_CONTROL"
    assert decision["DBTT"]["decision"] == "RETAIN_CONTROL"
    assert decision["barriers_changed"] is False
    assert decision["memory_terms_added"] is False
    assert decision["lifecycle_terms_changed"] is False
    assert decision["new_FEMCZM_runs"] == 0
    rows = decision["recommended_full_precision_taylor_peierls_rows"]
    population = pd.read_parquet(OUT / "oneD_v2_taylor_peierls_complete_population.parquet")
    for material, candidate_id in (
        ("Peak", "v913_zeroD_sobol_0242980"),
        ("DBTT", "v913_zeroD_sobol_0202500"),
    ):
        expected = population[population.candidate_id.eq(candidate_id)].iloc[0]
        assert rows[material]["candidate_id"] == candidate_id
        for field in SEARCH_WHITELIST:
            assert rows[material][field] == float(expected[field])


def test_required_figures_and_reports_exist() -> None:
    figures = (
        "PEAK_TAYLOR_PEIERLS_RESPONSE_MAP.png", "DBTT_TAYLOR_PEIERLS_RESPONSE_MAP.png",
        "TRANSPORT_RETENTION_TIMESCALE_MAP.png", "CHI_RET_VS_RCURVE_PROPENSITY.png",
        "APPARENT_VS_LOCAL_TOUGHENING.png", "RADIUS_BACKSTRESS_RETAINED_STATE_EVOLUTION.png",
        "CONTROL_VS_TAYLOR_PEIERLS_FINALISTS_REDUCED.png",
        "CONTROL_VS_TAYLOR_PEIERLS_FINALISTS_PF.png",
        "TAYLOR_PEIERLS_OPTION_BANK_MAP.png", "FINAL_TAYLOR_PEIERLS_DECISION.png",
    )
    for name in figures:
        path = OUT / name
        assert path.is_file() and path.stat().st_size > 10_000
    reports = (
        "ONE_D_V2_TAYLOR_PEIERLS_SEARCH_CONTRACT.md",
        "ONE_D_V2_TAYLOR_PEIERLS_TIMESCALE_ANALYSIS.md",
        "ONE_D_V2_PEAK_TAYLOR_PEIERLS_SEARCH.md",
        "ONE_D_V2_DBTT_TAYLOR_PEIERLS_SEARCH.md",
        "ONE_D_V2_TAYLOR_PEIERLS_MULTISEED_VALIDATION.md",
        "PF_2D_TAYLOR_PEIERLS_TRANSFER_VALIDATION.md",
        "ONE_D_V2_APPARENT_VS_LOCAL_TOUGHENING_DECOMPOSITION.md",
        "ONE_D_V2_TAYLOR_PEIERLS_OPTION_BANK.md",
        "ONE_D_V2_TAYLOR_PEIERLS_FINAL_DECISION.md",
    )
    for name in reports:
        assert (ROOT / name).is_file()


def test_completed_historical_and_v2_registries_are_unchanged() -> None:
    paths = (
        "analysis_outputs/oneD_v2_peak_dbtt_rcurve_search/oneD_v2_peak_dbtt_rcurve_registry.csv",
        "analysis_outputs/oneD_v2_peak_dbtt_rcurve_search/oneD_v2_fracture_option_bank_material_registry.csv",
        "analysis_outputs/oneD_v2_terminal_predictive_program/oneD_v2_new_four_class_registry.csv",
    )
    subprocess.run(["git", "diff", "--quiet", BASE, "--", *paths], cwd=ROOT, check=True)
