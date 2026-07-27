from pathlib import Path


def test_v10051832_runner_uses_joint_K_entry_and_manifest():
    root = Path(__file__).resolve().parents[1]
    text = (root / "run_v10_0_5_18_3_2_four_class_focused_100um.sh").read_text()
    assert "10.0.5.18.3.2" in text
    assert (
        "mode_i_first_passage_v10_0_5_18_3_2_four_class_"
        "joint_K_ramp_stochastic_emission"
    ) in text
    assert "persistent_site_production_manifest_v10_0_5_18_3_2.json" in text
    assert "outer_emission_action_limiter_active=false" in text
    assert "outer_emission_rate_variation_limiter_active=false" in text
    assert "K_ramp_evaluated_inside_event_horizon=true" in text
    assert "EMISSION_INNER_MAX_LOG_HAZARD_CHANGE" in text
    assert "EMISSION_ADAPT_" not in text
    assert "mode_i_first_passage_v10_0_5_18_2" not in text


def test_package_version_is_v10051832():
    root = Path(__file__).resolve().parents[1]
    text = (root / "pyproject.toml").read_text()
    assert 'version = "10.0.5.18.3.2"' in text
    assert "joint linear-K event-driven local kinetics" in text
