from pathlib import Path

from arrhenius_fracture.persistent_site_joint_K_single_trial_v10051832 import (
    PersistentSiteJointKSingleTrialFrontEngineV10051832,
    SINGLE_TRIAL_SCHEMA,
)


def test_single_trial_audit_preserves_exact_events_and_pf_operator():
    emission = PersistentSiteJointKSingleTrialFrontEngineV10051832.audit_payload()[
        "stochastic_emission"
    ]
    assert emission["single_trial_event_horizon_active"] is True
    assert emission["single_trial_event_horizon_schema"] == SINGLE_TRIAL_SCHEMA
    assert emission["state_dependent_log_rate_rejection_active"] is False
    assert emission["log_rate_change_is_diagnostic_only"] is True
    assert emission["PF_transport_operator_preserved"] is True
    assert emission["PF_transport_recursive_composition_error_control"] is False
    assert emission["threshold_crossings_localized_individually"] is True
    assert emission["events_batched"] is False
    assert emission["accepted_update"] == "exact_event_localized_one_activation_packets"


def test_single_trial_runner_uses_final_entry_and_explicit_policy():
    root = Path(__file__).resolve().parents[1]
    runner = (
        root
        / "run_v10_0_5_18_3_2_four_class_joint_K_single_trial_focused_100um.sh"
    ).read_text()
    assert (
        "mode_i_first_passage_v10_0_5_18_3_2_four_class_"
        "joint_K_single_trial_stochastic_emission"
    ) in runner
    assert "state_dependent_log_rate_rejection_active=false" in runner
    assert "log_rate_change_is_diagnostic_only=true" in runner
    assert "PF_transport_operator_preserved=true" in runner
    assert "events_batched=false" in runner
    assert "EMISSION_EVENT_HORIZON_FACTOR" in runner
    assert "EMISSION_MAX_EVENTS_PER_HALF_STEP" in runner
