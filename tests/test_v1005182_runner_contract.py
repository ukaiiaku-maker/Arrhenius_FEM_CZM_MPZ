from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT / "run_v10_0_5_18_2_four_class_focused_100um.sh"
PYPROJECT = ROOT / "pyproject.toml"
ENTRY = (
    ROOT
    / "arrhenius_fracture"
    / "mode_i_first_passage_v10_0_5_18_2_four_class_state_coupled_stochastic_emission.py"
)
ENGINE = (
    ROOT
    / "arrhenius_fracture"
    / "persistent_site_state_coupled_stochastic_emission_v1005182.py"
)


def test_package_version_matches_point_release():
    text = PYPROJECT.read_text()
    assert 'version = "10.0.5.18.2"' in text


def test_state_coupled_entry_and_runner_exist():
    assert ENTRY.is_file()
    assert ENGINE.is_file()
    assert RUNNER.is_file()


def test_runner_uses_state_coupled_exact_event_entry():
    text = RUNNER.read_text()
    assert (
        "arrhenius_fracture.mode_i_first_passage_v10_0_5_18_2_"
        "four_class_state_coupled_stochastic_emission"
    ) in text
    assert "emission_accepted_update=exact_event_localized" in text
    assert "emission_post_onset_mean_field_burst=false" in text
    assert "emission_signed_channel_factors_applied_once=true" in text
    assert "EMISSION_ADAPT_MAX_LOG_HAZARD_CHANGE" in text
    assert "EMISSION_ADAPT_MAX_BACKSTRESS_FRACTION" in text
    assert "EMISSION_ADAPT_MAX_GEOMETRY_FRACTION" in text
    assert "persistent_site_production_manifest_v10_0_5_18_2.json" in text
    assert "mode_i_first_passage_v10_0_5_18_1" not in text


def test_engine_does_not_contain_mean_field_burst_path():
    text = ENGINE.read_text()
    assert "post_onset_mean_field_burst\": False" in text
    assert "exact_event_localized_one_activation_packets" in text
    assert "emit_persistent(" in text
    assert "return adaptive_increment" in text
    assert "EMISSION_MAX_EVENTS_PER_HALF_STEP" not in text


def test_entry_patches_absolute_opening_only_for_two_channel_path():
    text = ENTRY.read_text()
    assert "two_channel_absolute_opening_drives" in text
    assert "CalibratedTipEngineMixin._mm_drives" in text
    assert "scalar_emission_factor_applied_to_opening_K" in text
    assert "emission_backstress_is_physical_limiter" in text
