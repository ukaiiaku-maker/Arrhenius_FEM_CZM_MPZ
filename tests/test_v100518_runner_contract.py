from __future__ import annotations

from pathlib import Path

RUNNER = Path("run_v10_0_5_18_four_class_focused_100um.sh")
ENTRY = Path(
    "arrhenius_fracture/"
    "mode_i_first_passage_v10_0_5_18_four_class_stochastic_emission.py"
)


def test_theta_is_runtime_input_not_material_selector():
    runner = RUNNER.read_text()
    entry = ENTRY.read_text()
    assert "THETA=${THETA:-30}" in runner
    assert '--crystal-theta-deg "$THETA"' in runner
    assert "theta_is_runtime_input=true" in runner
    assert "--crystal-theta-deg" not in entry
    assert "EXPECTED_OPTIONS" in entry


def test_focused_matrix_contains_only_new_primaries_by_default():
    text = RUNNER.read_text()
    assert "v913_paper_weakT01_0129902_persistent_sites" in text
    assert "v913_paper_ceramic01_0077080_persistent_sites" in text
    assert "v913_paper_weakT02_0008816_persistent_sites" not in text
    assert "v913_paper_ceramic02_0008536_persistent_sites" not in text
    assert 'TEMPERATURES=${TEMPERATURES:-"300 600 900 1100 1200"}' in text
    assert "TARGET_EXT_UM=${TARGET_EXT_UM:-100}" in text


def test_runner_requires_stochastic_cleavage_length_and_emission():
    text = RUNNER.read_text()
    required = (
        "CLEAVAGE_HAZARD_MODE=exponential",
        "CLEAVAGE_EVENT_LENGTH_MODE=threshold_scaled",
        "EMISSION_HAZARD_SEED=\"$CASE_SEED\"",
        "EMISSION_MAX_ACTION_SUBSTEP",
        "EMISSION_MAX_EVENTS_PER_HALF_STEP",
        "continuous_fractional_MPZ_translation=true",
        "microstructure_advance_precedes_cohesive_checkpoint=true",
    )
    for token in required:
        assert token in text


def test_case_seed_changes_with_option_temperature_and_replicate():
    text = RUNNER.read_text()
    assert "option={option}" in text
    assert "T={temperature}" in text
    assert "replicate={replicate}" in text
    assert "stream=signed_emission" in text
    assert "channel={channel}" in text


def test_baseline_mechanics_and_moving_process_zone_are_retained():
    text = RUNNER.read_text()
    required = (
        "--crack-backend adaptive_czm",
        "--crystal-aniso",
        "--crystal-compete",
        "--max-fronts 1",
        "--j-decomposition cluster",
        "--mpz-length-um 50",
        "--mpz-n-bins 80",
        "--adaptive-events",
        "--da-phys 5e-6",
    )
    for token in required:
        assert token in text
