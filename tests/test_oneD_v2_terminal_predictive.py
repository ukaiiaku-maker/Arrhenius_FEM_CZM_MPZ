from __future__ import annotations

from pathlib import Path
import sys
import hashlib
import json

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import run_oneD_v2_natural_forward as natural
from reduced_fracture_v2.predictive import (
    FEMCZM_PREDICTIVE_LIFECYCLE_V2,
    PF_PREDICTIVE_LIFECYCLE_V2,
    ProviderMechanicsMap,
    SourceDriveMapBoundsError,
    provider_loading_map,
)


class _PFEngine:
    def __init__(self):
        self.call = None

    def predict_clock_increment(self, K, T, dt):
        self.call = (K, T, dt)
        return 0.125


class _FEMEngine:
    def __init__(self):
        self.call = None

    def sigma_tip(self, K):
        self.call = ("sigma_tip", K)
        return 2.0 * K

    def lambda_cleave(self, sigma, T):
        self.call = ("lambda_cleave", sigma, T)
        return 0.5, None, None


def test_natural_forward_uses_backend_owned_predictors():
    pf = _PFEngine()
    fem = _FEMEngine()
    assert natural._predict(pf, "pf", 1.0, 2.0, 1000.0, 3.0) == 0.125
    assert pf.call == (1.0, 1000.0, 3.0)
    assert natural._predict(fem, "fem", 1.0, 2.0, 1000.0, 3.0) == 1.5
    assert fem.call == ("lambda_cleave", 2.0, 1000.0)


def test_fem_natural_forward_excludes_rejected_aggregate_closure():
    source = (ROOT / "scripts/run_oneD_v2_natural_forward.py").read_text()
    assert "aggregate_emission=True" not in source
    assert "constructor=_fem_constructor" in source


def test_natural_forward_has_no_threshold_priming():
    source = (ROOT / "scripts/run_oneD_v2_natural_forward.py").read_text()
    assert "_threshold_prime" not in source
    assert '"threshold_priming":False' in source


def test_provider_lifecycle_reductions_are_distinct_and_candidate_independent():
    assert PF_PREDICTIVE_LIFECYCLE_V2.backend == "PF"
    assert FEMCZM_PREDICTIVE_LIFECYCLE_V2.backend == "FEMCZM"
    assert PF_PREDICTIVE_LIFECYCLE_V2.threshold_mode == "DETERMINISTIC_UNIT_ACTION"
    assert FEMCZM_PREDICTIVE_LIFECYCLE_V2.threshold_mode == "EXPONENTIAL_SOURCE_THRESHOLD"
    mechanics = ProviderMechanicsMap(
        "PF", np.array([0.0, 100.0]), np.ones(2), np.ones(2)
    )
    first = provider_loading_map(mechanics, seed=1, target_extension_m=20.0e-6)
    second = provider_loading_map(mechanics, seed=999, target_extension_m=20.0e-6)
    assert first.threshold_actions == second.threshold_actions == (1.0,) * 4
    assert first.path_advances_m == second.path_advances_m == (5.0e-6,) * 4


def test_terminal_native_baselines_match_topology_without_silent_bounds():
    path = ROOT / "analysis_outputs/oneD_v2_terminal_predictive_program/oneD_v2_native_baseline_results.csv"
    frame = pd.read_csv(path).set_index(["provider", "material_class"])
    assert frame.status.eq("TARGET_RIGHT_CENSORED").all()
    expected = {("PF", "Peak"): 1, ("PF", "DBTT"): 2,
                ("FEMCZM", "Peak"): 1, ("FEMCZM", "DBTT"): 3}
    assert frame.physical_avalanche_count.astype(int).to_dict() == expected
    assert frame.topology_match.astype(bool).all()


def test_source_map_out_of_domain_fails_closed():
    from scripts.run_oneD_v2_predictive_campaign import inputs

    _, _, providers = inputs()
    drive = providers["PF"][1]
    try:
        drive.evaluate(0.0, 1000.1e-6)
    except SourceDriveMapBoundsError:
        pass
    else:
        raise AssertionError("source-drive map silently extrapolated")


def test_terminal_sensitivity_and_search_artifacts_are_deterministic_and_shared():
    out = ROOT / "analysis_outputs/oneD_v2_terminal_predictive_program"
    sensitivity = out / "oneD_v2_parameter_sensitivity_results.parquet"
    manifest = json.loads((out / "oneD_v2_parameter_sensitivity_manifest.json").read_text())
    assert hashlib.sha256(sensitivity.read_bytes()).hexdigest() == manifest["results_sha256"]
    frame = pd.read_parquet(sensitivity)
    assert len(frame) == 972
    assert frame.status.eq("TARGET_RIGHT_CENSORED").all()
    registry = pd.read_csv(out / "oneD_v2_new_four_class_registry.csv")
    assert len(registry) == 4
    assert registry.same_material_row_both_providers.astype(bool).all()
    assert registry.candidate_id.is_unique
    pareto = pd.read_csv(out / "oneD_v2_pareto_candidates.csv")
    eligible = pareto[pareto.full_class_contract_pass]
    assert eligible[["search_class", "candidate_id"]].to_dict("records") == [{
        "search_class": "weak-T",
        "candidate_id": "oneD_v2_focused_weak_T_0016",
    }]


def test_terminal_pf_transfer_and_final_matrix_are_bounded():
    out = ROOT / "analysis_outputs/oneD_v2_terminal_predictive_program"
    transfer = pd.read_csv(out / "oneD_v2_pf_transfer_results.csv")
    assert len(transfer) == 12
    assert int(transfer.new_PF_run_launched.sum()) == 6
    assert not transfer.new_FEMCZM_run_launched.any()
    assert transfer.pf_2D_target_right_censored.astype(bool).all()
    assert (transfer.pf_2D_event_count > transfer.pf_2D_physical_avalanche_count).all()
    assert int(transfer.topology_match.sum()) == 11
    final = pd.read_csv(out / "oneD_v2_final_four_class_results.csv")
    assert len(final) == 40
    assert final.status.eq("TARGET_RIGHT_CENSORED").all()


def test_terminal_reports_and_predecessor_audits_are_preserved():
    required = [
        "ONE_D_V2_FINAL_MODEL_ARCHITECTURE.md",
        "ONE_D_V2_NATIVE_BASELINE_VALIDATION.md",
        "ONE_D_V2_PARAMETER_SENSITIVITY_VALIDATION.md",
        "ONE_D_V2_DOMAIN_OF_USEFULNESS.md",
        "ONE_D_V2_CURRENT_PARAMETER_DIAGNOSTIC.md",
        "ONE_D_V2_PARAMETER_SEARCH.md",
        "ONE_D_V2_NEW_FOUR_CLASS_SELECTION.md",
        "ONE_D_V2_TO_PF_TRANSFER_VALIDATION.md",
        "ONE_D_V2_FINAL_PARAMETER_AND_2D_DECISION.md",
    ]
    assert all((ROOT / name).is_file() for name in required)
    predecessor = ROOT / "analysis_outputs/oneD_v2_predictive_model/earlier_audit_reports"
    assert predecessor.is_dir()
    assert len(list(predecessor.glob("*.md"))) >= 7


def test_terminal_decision_and_domain_are_fail_closed():
    out = ROOT / "analysis_outputs/oneD_v2_terminal_predictive_program"
    decision = json.loads((out / "oneD_v2_final_decision.json").read_text())
    domain = json.loads((out / "oneD_v2_domain_of_usefulness.json").read_text())
    assert decision["mission_status"] == "COMPLETE"
    assert decision["new_bounded_PF_case_count"] == 6
    assert decision["new_2D_FEMCZM_runs"] == 0
    assert decision["production_physical_formulas_changed"] is False
    assert decision["canonical_production_trajectories_changed"] is False
    assert domain["common"]["outside_domain_policy"].endswith(
        "fail closed; no clipping or extrapolation"
    )
    assert domain["common"]["target_termination"] == "RIGHT_CENSORED_NOT_PHYSICAL_ARREST"


def test_terminal_provenance_hashes_and_figure_contract():
    out = ROOT / "analysis_outputs/oneD_v2_terminal_predictive_program"
    manifest = json.loads((out / "oneD_v2_provenance_manifest.json").read_text())
    assert manifest["producer_branch"] == "codex/oneD-v2-terminal-predictive-program"
    assert manifest["maximum_concurrent_PF_workers"] <= 2
    assert manifest["new_2D_PF_runs"] == 6
    assert manifest["new_2D_FEMCZM_runs"] == 0
    assert manifest["canonical_PF_registry_modified"] is False
    assert manifest["canonical_production_trajectories_changed"] is False
    assert manifest["production_physical_formulas_changed"] is False
    for name, expected in manifest["required_output_sha256"].items():
        assert hashlib.sha256((out / name).read_bytes()).hexdigest() == expected
    png_signature = b"\x89PNG\r\n\x1a\n"
    for name, expected in manifest["required_figure_sha256"].items():
        data = (out / name).read_bytes()
        assert data.startswith(png_signature)
        assert len(data) > 10_000
        assert hashlib.sha256(data).hexdigest() == expected
