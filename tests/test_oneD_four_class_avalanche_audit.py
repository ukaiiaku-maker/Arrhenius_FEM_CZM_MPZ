"""Fail-closed semantic tests for the four-class 1-D V2 audit."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd

from scripts.audit_oneD_four_class_v2 import ACTIVE, CLASSES, REJECTED, canon, group_boundaries


ROOT = Path(__file__).parents[1]
OUT = ROOT / "analysis_outputs/oneD_four_class_audit/v2"
REGISTRY = Path("/private/tmp/oneD-pf-consistency-audit/arrhenius_fracture/data/materials/v10_2_27_v913_four_class_paper_registry.csv")


def read(name):
    return pd.read_csv(OUT / name)


def test_01_exact_canonical_option_keys():
    assert {v[1] for v in CLASSES.values()} == set(read("oneD_parameter_consistency.csv").option_key)


def test_02_exact_candidate_ids():
    assert set(CLASSES) == set(read("oneD_parameter_consistency.csv").candidate_id)


def test_03_full_precision_registry_equality():
    reg = pd.read_csv(REGISTRY, dtype=str)
    got = dict(zip(read("oneD_parameter_consistency.csv").candidate_id, read("oneD_parameter_consistency.csv").active_parameter_fingerprint))
    assert all(got[r.candidate_id] == canon({f: str(getattr(r, f)) for f in ACTIVE}) for r in reg.itertuples())


def test_04_rejected_candidates_excluded():
    assert not (set(read("oneD_parameter_consistency.csv").candidate_id) & REJECTED)


def test_05_backup_candidates_excluded():
    assert len(read("oneD_parameter_consistency.csv")) == 4


def test_06_exact_historical_reproduction_where_archive_exists():
    x = read("oneD_historical_reproduction_results.csv")
    matched = x[x.material_class.isin(["weak-T", "ceramic-like"])]
    assert len(matched) == 10 and matched.difference_MPa_sqrt_m.abs().max() < 1e-12


def test_07_one_physical_engine_not_represented_as_two():
    inv = json.loads((OUT / "oneD_fracture_model_inventory.json").read_text())
    ids = [x["id"] for x in inv["implementations"]]
    assert len(ids) == len(set(ids)) and ids.count("AUTONOMOUS_V913_ONE_D_DRIVER") == 1


def test_08_frozen_state_equality_is_only_required_where_intended():
    x = read("oneD_cross_implementation_frozen_states.csv")
    assert set(x.comparison_status) == {"INTENDED_MODEL_DIFFERENCE_NOT_COMPARABLE"}


def test_09_event_transactions_differ_from_avalanches():
    assert len(read("oneD_event_transactions_v2.csv")) > len(read("oneD_physical_avalanches_v2.csv"))


def test_10_pre_and_post_event_fields_are_distinct():
    cols = set(read("oneD_event_transactions_v2.csv").columns)
    assert {"pre_event_native_drive_MPa_sqrt_m", "post_event_native_drive_at_fixed_load_MPa_sqrt_m"} <= cols


def test_11_onset_candidates_are_pre_event():
    assert set(read("oneD_onset_resistance_candidates_v2.csv").onset_role) <= {"INITIAL_ONSET_PRE_EVENT", "REINITIATION_ONSET_PRE_EVENT"}


def test_12_in_avalanche_states_are_not_resistance_points():
    assert not read("oneD_in_avalanche_driving_force_v2.csv").is_resistance_point.any()


def test_13_target_termination_is_right_censored_not_arrested():
    x = read("oneD_physical_avalanches_v2.csv")
    assert x.right_censored_at_target.any() and not x.classification.str.contains("ARREST", case=False).any()


def test_14_grouping_is_deterministic():
    e = [{"applied_displacement_m": x} for x in (1.0, 1.000001, 1.000002, 1.01)]
    assert group_boundaries(e) == group_boundaries(e)


def test_15_missing_grouping_evidence_is_unclassified():
    e = [{"applied_displacement_m": 1.0}, {"applied_displacement_m": 1.0}]
    assert group_boundaries(e) == (None, None)


def test_16_grouping_and_threshold_analysis_advance_no_rng():
    source = (ROOT / "scripts/audit_oneD_four_class_v2.py").read_text()
    assert "np.random" not in source and "random." not in source


def test_17_loading_contract_is_explicit():
    dep = json.loads((OUT / "oneD_driving_force_dependency_graph.json").read_text())
    assert dep["autonomous_v913"]["external"] == "applied displacement U"


def test_18_unchanged_campaign_has_no_parameter_changes():
    x = read("oneD_parameter_consistency.csv")
    assert x.registry_full_precision_equal.all() and set(x.closure_status) == {"PASS"}


def test_19_legacy_rise_is_retained_but_deprecated():
    assert "legacy_eventwise_endpoint_K_MPa_sqrt_m" in read("oneD_four_class_temperature_summary_v2.csv")
    assert "deprecated" in (OUT / "ONE_D_DRIVING_FORCE_SEMANTICS_AUDIT.md").read_text().lower()


def test_20_no_geometry_mixed_comparison_enters_v2_onsets():
    x = read("oneD_onset_resistance_candidates_v2.csv")
    assert set(x.onset_native_drive_kind) == {"MODEL_NATIVE_EVENT_INDEXED_K"} and not x.is_structural_G.any()


def test_21_no_production_2d_artifact_is_an_output():
    assert not any(p.name.startswith("fem_") and p.suffix in {".csv", ".json"} for p in OUT.iterdir())


def test_22_scientific_fingerprints_reproduce_files():
    manifest = json.loads((OUT / "oneD_historical_reproduction_manifest.json").read_text())
    for name, expected in manifest["scientific_fingerprints"].items():
        assert hashlib.sha256((OUT / name).read_bytes()).hexdigest() == expected
