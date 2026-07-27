from pathlib import Path


def test_v10051831_runner_uses_single_layer_entrypoint_and_manifest():
    root = Path(__file__).resolve().parents[1]
    text = (root / "run_v10_0_5_18_3_1_four_class_focused_100um.sh").read_text()
    assert "10.0.5.18.3.1" in text
    assert (
        "mode_i_first_passage_v10_0_5_18_3_1_four_class_"
        "single_layer_load_ramp_stochastic_emission"
    ) in text
    assert "persistent_site_production_manifest_v10_0_5_18_3_1.json" in text
    assert "outer_emission_action_limiter_active=false" in text
    assert "event_driven_emission_transport_active=true" in text
    assert "fixed_inner_action_substep_per_event=false" in text
    assert "emission_threshold_crossings_localized_individually=true" in text
    assert "emission_events_batched=false" in text
    assert "inner_exact_emission_action_control=true" in text
    assert "EMISSION_EVENT_HORIZON_FACTOR" in text
    assert "EMISSION_INNER_MAX_LOG_HAZARD_CHANGE" in text
    assert "EMISSION_ADAPT_" not in text
    assert (
        "mode_i_first_passage_v10_0_5_18_2_four_class_"
        "state_coupled_stochastic_emission"
    ) not in text


def test_package_version_is_v10051831():
    root = Path(__file__).resolve().parents[1]
    text = (root / "pyproject.toml").read_text()
    assert 'version = "10.0.5.18.3.1"' in text
