from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from reduced_fracture_v2.predictive import (
    ProviderMechanicsMap,
    SourceDriveMap,
    SourceDriveMapBoundsError,
)


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_outputs" / "oneD_v2_predictive_model"
IDS = {
    "Peak": "v913_zeroD_sobol_0242980",
    "DBTT": "v913_zeroD_sobol_0202500",
    "weak-T": "v913_zeroD_sobol_0129902",
    "ceramic-like": "v913_zeroD_sobol_0077080",
}


def test_source_drive_map_is_bounded_and_has_no_clipping():
    source = SourceDriveMap.from_csv(OUT / "oneD_v2_pf_source_drive_map.csv")
    assert np.array_equal(
        source.radius_um,
        np.asarray((1.0, 2.0, 4.0, 8.0, 12.0, 25.0, 50.0, 100.0)),
    )
    assert len(source.evaluate(500.0e-6, 25.0e-6)) == 2
    with pytest.raises(SourceDriveMapBoundsError):
        source.evaluate(500.0e-6, 100.0001e-6)
    with pytest.raises(SourceDriveMapBoundsError):
        source.evaluate(1000.0001e-6, 25.0e-6)


def test_mechanics_map_fails_closed_outside_qualified_extension():
    mechanics = ProviderMechanicsMap.from_csv(
        "PF",
        ROOT
        / "analysis_outputs/oneD_v2_mechanics_maps_and_baselines/oneD_v2_pf_mechanics_map.csv",
    )
    assert np.isfinite(mechanics.native_K_per_opening(1000.0e-6))
    with pytest.raises(ValueError, match="outside qualified map"):
        mechanics.native_K_per_opening(1000.1e-6)


def test_current_matrix_and_long_baselines_are_complete_artifacts():
    current = pd.read_csv(OUT / "oneD_v2_current_four_class_results.csv")
    long = pd.read_csv(OUT / "oneD_v2_1000um_baseline_results.csv")
    assert len(current) == 40
    assert set(current.material_class) == set(IDS)
    assert set(current.provider) == {"PF", "FEMCZM"}
    assert set(current.temperature_K) == {300.0, 600.0, 900.0, 1100.0, 1200.0}
    assert len(long) == 4
    dbtt = long[long.material_class == "DBTT"]
    assert set(dbtt.status) == {"RIGHT_CENSORED_DRIVE_MAP_BOUND"}
    assert (dbtt.terminal_extension_um < 1000.0).all()
    peak = long[long.material_class == "Peak"]
    assert set(peak.status) == {"TARGET_RIGHT_CENSORED"}
    assert np.allclose(peak.terminal_extension_um, 1000.0)


def test_parameter_search_is_fail_closed_and_multiobjective():
    manifest = json.loads((OUT / "oneD_v2_parameter_search_manifest.json").read_text())
    pareto = pd.read_csv(OUT / "oneD_v2_pareto_candidates.csv")
    population = pd.read_parquet(OUT / "oneD_v2_search_population.parquet")
    assert manifest["screen_candidate_count"] == 187
    assert manifest["validation_candidate_count"] == 16
    assert manifest["eligible_shared_dbtt_rows"] == []
    assert population[population.search_stage == "SCREEN"].candidate_id.nunique() == 187
    assert (
        population[population.search_stage == "REPEATED_VALIDATION"].candidate_id.nunique()
        == 16
    )
    objectives = {
        "class_topology_objective",
        "onset_envelope_objective",
        "precursor_count_objective",
        "avalanche_topology_objective",
        "event_size_objective",
        "provider_robustness_objective",
        "state_realism_objective",
    }
    assert objectives <= set(pareto)
    assert not pareto.full_dbtt_contract_pass.any()
    assert set(pareto.selection_status) == {"REJECTED_FULL_DBTT_CONTRACT"}


def test_final_registry_retains_exact_reference_ids_without_false_promotion():
    registry = pd.read_csv(OUT / "oneD_v2_new_four_class_registry.csv")
    assert dict(zip(registry.material_class, registry.candidate_id)) == IDS
    assert set(registry.parameter_decision) == {
        "CURRENT_ROW_RETAINED_NO_ADMISSIBLE_REPLACEMENT"
    }
    assert registry.loc[
        registry.material_class == "DBTT", "final_selection_status"
    ].item() == "RETAINED_PRODUCTION_REFERENCE_REDUCED_NOT_QUALIFIED"
    assert registry.loc[
        registry.material_class == "ceramic-like", "reduced_provider_robust"
    ].item()


def test_pf_transfer_uses_existing_authoritative_cases_only():
    transfer = pd.read_csv(OUT / "oneD_v2_pf_transfer_results.csv")
    assert len(transfer) == 12
    assert set(transfer.temperature_K) == {300.0, 1000.0, 1200.0}
    assert set(transfer.transfer_data_role) == {"EXISTING_AUTHORITATIVE_2D_PF"}
    assert not transfer.new_PF_run_launched.any()
    first = transfer.pivot(
        index="material_class",
        columns="temperature_K",
        values="pf_2D_initial_native_KJ_MPa_sqrt_m",
    )
    assert first.loc["Peak", 1000.0] > first.loc["Peak", 300.0]
    assert first.loc["Peak", 1200.0] < first.loc["Peak", 1000.0]
    assert first.loc["DBTT", 1000.0] > first.loc["DBTT", 300.0]
    assert np.ptp(first.loc["weak-T"].to_numpy(float)) < 2.0
    assert first.loc["ceramic-like", 1200.0] < first.loc["ceramic-like", 300.0]


def test_final_decision_matches_machine_outputs_and_run_policy():
    decision = json.loads((OUT / "oneD_v2_final_decision.json").read_text())
    assert decision["overall_status"] == "PARTIALLY_QUALIFIED_NO_NEW_ROW_PROMOTED"
    assert decision["controlled_common_kernel_parity"] == "PASS_EXACT"
    assert decision["new_parameter_rows_promoted"] == 0
    assert decision["new_2D_PF_runs"] == 0
    assert decision["new_2D_FEMCZM_runs"] == 0
    assert decision["production_trajectories_altered"] is False
    assert decision["physical_formulas_altered"] is False
    assert decision["final_candidate_ids"] == IDS
    assert decision["minimum_later_FEMCZM_matrix"] == [
        "Peak, 1000 K, theta=0, tip-only",
        "DBTT, 1000 K, theta=0, tip-only",
    ]


@pytest.mark.parametrize(
    "name",
    [
        "ONE_D_V2_PROVIDER_VS_OWN_2D_BASELINES.png",
        "ONE_D_V2_FOUR_CLASS_ONSET_ENVELOPES.png",
        "ONE_D_V2_AVALANCHE_TOPOLOGY_VS_TEMPERATURE.png",
        "ONE_D_V2_CURRENT_VS_FINAL_CLASSIFICATIONS.png",
        "ONE_D_V2_DIRECT_PF_TRANSFER_VALIDATION.png",
        "ONE_D_V2_FINAL_FOUR_CLASS_SUMMARY.png",
    ],
)
def test_required_final_figure_exists(name):
    path = OUT / name
    assert path.is_file()
    assert path.stat().st_size > 10_000
