from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PYPROJECT = ROOT / "pyproject.toml"
RUNNER = ROOT / "run_v10_0_5_18_3_four_class_focused_100um.sh"
ENTRY = (
    ROOT
    / "arrhenius_fracture"
    / "mode_i_first_passage_v10_0_5_18_3_four_class_load_ramp_stochastic_emission.py"
)
ENGINE = (
    ROOT
    / "arrhenius_fracture"
    / "persistent_site_load_ramp_stochastic_emission_v1005183.py"
)
RESILIENCE = ROOT / "arrhenius_fracture" / "numerical_resilience_v1005183.py"


def test_package_metadata_matches_v1005183():
    text = PYPROJECT.read_text()
    assert 'version = "10.0.5.18.3"' in text
    assert "local linear-K load-ramp kinetics" in text


def test_v1005183_files_exist():
    assert RUNNER.is_file()
    assert ENTRY.is_file()
    assert ENGINE.is_file()
    assert RESILIENCE.is_file()


def test_runner_uses_load_ramp_entry_without_v1005182_state_predictor():
    text = RUNNER.read_text()
    assert (
        "arrhenius_fracture.mode_i_first_passage_v10_0_5_18_3_"
        "four_class_load_ramp_stochastic_emission"
    ) in text
    assert "global_FEM_deterministic_emission_state_predictor=false" in text
    assert "local_linear_K_ramp_kinetics=true" in text
    assert "emission_post_onset_mean_field_burst=false" in text
    assert "emission_signed_channel_factors_applied_once=true" in text
    assert "EMISSION_RAMP_MAX_LOG_RATE_CHANGE" in text
    assert "EMISSION_RAMP_ACTION_FLOOR" in text
    assert "TWO_CHANNEL_PROBE_FALLBACK_MAX_CALLS" in text
    assert "EMISSION_ADAPT_MAX_LOG_HAZARD_CHANGE" not in text
    assert "EMISSION_ADAPT_MAX_BACKSTRESS_FRACTION" not in text
    assert "EMISSION_ADAPT_MAX_GEOMETRY_FRACTION" not in text
    assert "mode_i_first_passage_v10_0_5_18_2" not in text
    assert "persistent_site_production_manifest_v10_0_5_18_3.json" in text


def test_load_ramp_engine_has_no_deterministic_expected_emission_gate():
    text = ENGINE.read_text()
    assert "emit_persistent(" not in text
    assert "global_FEM_state_change_predictor\": False" in text
    assert "local_linear_K_ramp\": True" in text
    assert "exact_event_localized_one_activation_packets" in text
    assert "post_onset_mean_field_burst\": False" in text


def test_entry_installs_and_restores_only_v1005183_runtime_hooks():
    text = ENTRY.read_text()
    assert "ProductionLoadRampStochasticEmissionFrontEngineV1005183" in text
    assert "two_channel_absolute_opening_drives_v1005183" in text
    assert "make_robust_process_zone_traction_probe" in text
    assert "RetryingAdaptiveCZMBackendV1005183" in text
    assert "finally:" in text
    assert "_mm.process_zone_traction_probe = saved[\"traction_probe\"]" in text
    assert "_crack_backend.AdaptiveCZMBackend = saved[\"adaptive_czm\"]" in text


def test_numerical_resilience_preserves_event_length_and_does_not_refit_physics():
    text = RESILIENCE.read_text()
    assert "physical_event_requested_length_m" in text
    assert "partition_retry_exhausted" in text
    assert "requested_length_m" in text
    assert "parameter" not in text.lower()
    assert "barrier" not in text.lower()
