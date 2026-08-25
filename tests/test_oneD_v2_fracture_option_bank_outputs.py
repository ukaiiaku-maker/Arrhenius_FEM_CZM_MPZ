import hashlib
import json
from pathlib import Path

import pandas as pd

from arrhenius_fracture.emergent_gnd_contract_v913 import ACTIVE_CANDIDATE_PARAMETER_FIELDS
from scripts.finalize_oneD_v2_peak_dbtt_R_and_option_bank import (
    BANK_VERSION,
    CONTROLS,
    FOCUSED,
    FORBIDDEN_MATERIAL_FIELDS,
    ROOT,
    canonical_json,
    material_hash,
)


OUT = ROOT / "analysis_outputs" / "oneD_v2_peak_dbtt_rcurve_search"
BANK = OUT / "oneD_v2_fracture_option_bank.csv"
MATERIALS = OUT / "oneD_v2_fracture_option_bank_material_registry.csv"
RESPONSES = OUT / "oneD_v2_fracture_option_bank_response_features.csv"
MANIFEST = OUT / "oneD_v2_fracture_option_bank_manifest.json"
PROVENANCE = OUT / "oneD_v2_peak_dbtt_R_provenance_manifest.json"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, float_precision="round_trip")


def test_material_identity_is_canonical_full_precision_and_collapsed():
    bank = load_csv(BANK)
    materials = load_csv(MATERIALS)
    assert bank.candidate_version.eq(BANK_VERSION).all()
    assert bank.candidate_id.notna().all()
    assert not bank.duplicated(["target_response_class", "candidate_id"]).any()
    repeated_ids = bank[bank.candidate_id.duplicated(False)]
    assert repeated_ids.groupby("candidate_id").parameter_sha256.nunique().eq(1).all()
    assert materials.parameter_sha256.is_unique
    assert len(materials) == bank.parameter_sha256.nunique()
    for _, row in materials.iterrows():
        assert row.parameter_sha256 == material_hash(row)
        payload = json.loads(row.canonical_parameter_json)
        assert set(payload) == set(ACTIVE_CANDIDATE_PARAMETER_FIELDS)
        assert all(payload[field] == float(row[field]) for field in ACTIVE_CANDIDATE_PARAMETER_FIELDS)
        assert row.canonical_parameter_json == canonical_json(row)


def test_material_registry_excludes_backend_lifecycle_fields_and_records_domain():
    materials = load_csv(MATERIALS)
    assert not (set(materials.columns) & FORBIDDEN_MATERIAL_FIELDS)
    assert materials.active_parameter_count.eq(len(ACTIVE_CANDIDATE_PARAMETER_FIELDS)).all()
    assert set(materials.active_parameter_schema) == {"V913_ACTIVE_CANDIDATE_29"}
    assert materials.qualified_domain_sha256.nunique() == 1


def test_controls_fixed_rows_and_focused_variants_are_preserved():
    bank = load_csv(BANK)
    ids = set(bank.candidate_id)
    assert set(CONTROLS.values()) <= ids
    assert {item for values in FOCUSED.values() for item in values} <= ids
    controls = bank[bank.candidate_id.isin(CONTROLS.values())]
    assert controls.option_status.str.contains("HISTORICAL_CONTROL").all()
    assert controls.option_status.str.contains("CURRENT_V2_SELECTED").all()
    focused = bank[bank.candidate_id.isin({item for values in FOCUSED.values() for item in values})]
    assert focused.option_status.str.contains("R_ENRICHED_OPTION").any()


def test_curated_bank_and_shortlist_span_all_classes_and_responses():
    bank = load_csv(BANK)
    shortlist = load_csv(OUT / "oneD_v2_future_joint_search_shortlist.csv")
    expected_bank = {"Peak": 16, "DBTT": 16, "weak-T": 14, "ceramic-like": 12}
    expected_shortlist = {"Peak": 6, "DBTT": 6, "weak-T": 4, "ceramic-like": 4}
    assert bank.groupby("target_response_class").size().to_dict() == expected_bank
    assert shortlist.groupby("target_response_class").size().to_dict() == expected_shortlist
    for material, local in bank.groupby("target_response_class"):
        assert local.parameter_sha256.nunique() >= len(local) - 1
        assert local.option_role.nunique() >= 4
        assert any(local[field].nunique() > 3 for field in ACTIVE_CANDIDATE_PARAMETER_FIELDS)
    responses = load_csv(RESPONSES)
    assert set(bank.candidate_id) <= set(responses.candidate_id)
    keys = ["target_response_class", "candidate_id", "parameter_sha256"]
    joined = responses[keys].drop_duplicates()
    assert not joined.duplicated(["target_response_class", "candidate_id"]).any()
    expected = bank[keys].merge(
        joined, on=keys, validate="one_to_one"
    )
    assert len(expected) == len(bank)


def test_direct_pf_references_and_manifest_artifacts_are_hash_qualified():
    bank = load_csv(BANK)
    for payload in bank.direct_PF_artifacts_json:
        for artifact in json.loads(payload):
            path = Path(artifact["source_steps_file"])
            assert path.is_file()
            assert sha(path) == artifact["source_steps_sha256"]
    manifest = json.loads(MANIFEST.read_text())
    for name, expected in manifest["artifacts"].items():
        assert sha(OUT / name) == expected
    for relative, expected in manifest["immutable_input_artifacts"].items():
        assert sha(ROOT / relative) == expected


def test_domain_hash_and_femczm_decision_are_cross_output_consistent():
    bank = load_csv(BANK)
    materials = load_csv(MATERIALS)
    manifest = json.loads(MANIFEST.read_text())
    provenance = json.loads(PROVENANCE.read_text())
    decision = json.loads((OUT / "oneD_v2_peak_dbtt_R_final_decision.json").read_text())
    expected = manifest["qualified_domain_sha256"]
    assert expected == sha(ROOT / manifest["qualified_domain_path"])
    assert set(bank.qualified_domain_sha256) == {expected}
    assert set(materials.qualified_domain_sha256) == {expected}
    assert provenance["qualified_domain_sha256"] == expected
    sentence = "No FEM/CZM material row is changed; no new FEM/CZM simulation was run."
    assert decision["FEMCZM_decision"] == sentence
    assert sentence in (ROOT / "ONE_D_V2_PEAK_DBTT_R_FINAL_DECISION.md").read_text()
    assert manifest["new_FEMCZM_runs"] == provenance["new_FEMCZM_runs"] == 0


def test_fatigue_is_explicitly_not_evaluated_everywhere():
    bank = load_csv(BANK)
    complete = pd.read_parquet(OUT / "oneD_v2_complete_material_option_population.parquet")
    shortlist = load_csv(OUT / "oneD_v2_future_joint_search_shortlist.csv")
    for frame in (bank, complete, shortlist):
        assert not frame.fatigue_evaluated.astype(bool).any()
        assert set(frame.fatigue_validation_status) == {"NOT_EVALUATED"}
    manifest = json.loads(MANIFEST.read_text())
    assert manifest["fatigue_evaluated"] is False
    assert manifest["fatigue_validation_status"] == "NOT_EVALUATED"


def test_deterministic_membership_and_material_fingerprints_reproduce():
    bank = load_csv(BANK)
    materials = load_csv(MATERIALS)
    manifest = json.loads(MANIFEST.read_text())
    membership = bank[["target_response_class", "candidate_id", "parameter_sha256"]].sort_values(
        ["target_response_class", "candidate_id"]
    ).to_dict("records")
    identity = materials[[
        "candidate_id", "parameter_sha256", *ACTIVE_CANDIDATE_PARAMETER_FIELDS
    ]].sort_values("parameter_sha256").to_dict("records")
    assert hashlib.sha256(json.dumps(
        membership, sort_keys=True, separators=(",", ":")
    ).encode()).hexdigest() == manifest["bank_membership_fingerprint"]
    assert hashlib.sha256(json.dumps(
        identity, sort_keys=True, separators=(",", ":")
    ).encode()).hexdigest() == manifest["material_identity_fingerprint"]
