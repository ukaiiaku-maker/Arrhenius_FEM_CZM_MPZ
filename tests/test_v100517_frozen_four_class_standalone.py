from __future__ import annotations

from pathlib import Path

from arrhenius_fracture.frozen_four_class_registry_v100517 import (
    EXPECTED_ACTIVE_FINGERPRINT_SHA256,
    EXPECTED_OPTIONS,
    FROZEN_REGISTRY_SHA256,
    load_frozen_parameter_option,
    validate_frozen_registry,
)

RUNNER = Path("run_v10_0_5_17_frozen_four_class_300_1000K_20um_smoke.sh")


def test_frozen_registry_contains_exact_four_classes_without_pf_dependency():
    payload = validate_frozen_registry()
    assert payload["n_options"] == 4
    assert payload["options_in_production_order"] == list(EXPECTED_OPTIONS)
    assert payload["classes"] == ["peak", "DBTT", "weakT", "ceramic"]
    assert payload["registry_sha256"] == FROZEN_REGISTRY_SHA256
    assert (
        payload["active_parameter_fingerprint_sha256"]
        == EXPECTED_ACTIVE_FINGERPRINT_SHA256
    )
    assert payload["phase_field_solver_dependency"] is False
    assert payload["phase_field_result_comparison_enabled"] is False


def test_each_frozen_option_loads_exact_candidate_and_closed_sources():
    for option, (candidate_id, material_class) in EXPECTED_OPTIONS.items():
        candidate, audit = load_frozen_parameter_option(option)
        assert candidate.option_key == option
        assert candidate.candidate_id == candidate_id
        assert audit["material_class"] == material_class
        assert audit["selected_through_frozen_FEM_registry"] is True
        assert audit["phase_field_solver_dependency"] is False
        assert audit["phase_field_result_comparison_enabled"] is False
        assert audit["persistent_sites"] is True
        assert audit["finite_source_inventory"] is False
        assert audit["source_refresh"] is False
        assert audit["explicit_recovery"] is False


def test_standalone_smoke_has_no_live_pf_dependency():
    text = RUNNER.read_text()
    forbidden = (
        "PFROOT",
        "--pf-repo-root",
        "PF_BRANCH",
        "PF_COMMIT",
        "phase-field checkout",
    )
    for token in forbidden:
        assert token not in text
    assert 'TEMPERATURES=${TEMPERATURES:-"300 1000"}' in text
    assert "TARGET_EXT_UM=${TARGET_EXT_UM:-20}" in text
    assert "MAX_JOBS=${MAX_JOBS:-2}" in text
    assert "RUN_J_SCALE_DIAGNOSTIC=${RUN_J_SCALE_DIAGNOSTIC:-0}" in text
    assert "phase_field_solver_dependency=false" in text
    assert "phase_field_result_comparison_enabled=false" in text
    assert "mode_i_first_passage_v10_0_5_17_four_class_standalone" in text
    assert (
        "a876ea042bf291ca9607b586205d75ce2b5cb95c668fef53d713d611dfe59cde"
        in text
    )
