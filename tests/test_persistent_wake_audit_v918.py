from __future__ import annotations

import json
from pathlib import Path

from arrhenius_fracture.coupled_event_audit_v918 import audit_case_v918


def _event(event_id: int, next_nucleation_wake: float = 0.0):
    rows = []
    for i, q in enumerate((0.25, 0.5, 0.75, 1.0), 1):
        rows.append({
            "damage_progress": q,
            "cleavage_hazard_increment": 0.25,
            "remote_displacement_increment_m": 0.0,
            "dt_s": 2.5,
            "lambda_c_current_s-1": 0.1,
            "lambda_c_raw_current_s-1": 1.0,
            "lambda_e_current_s-1": 1.0e6,
            "Kshield_Pa_sqrt_m": 2.0e4,
            "retained_count": 1.0,
        })
    return {
        "event_id": event_id,
        "status": "complete_committed",
        "progress_law": "dq_dt=lambda_c_absolute",
        "relaxation_completed": True,
        "final_damage": 1.0,
        "opening_hazard_integral": 1.0,
        "mpz_renewal_deferred": True,
        "mpz_renewal_committed": True,
        "committed_nfire": 1,
        "committed_distance_m": 5.0e-6,
        "source_sites_refreshed_during_opening": 1.0,
        "dN_emit_relaxation": 1.0,
        "persistent_wake_state_committed": True,
        "wake_retained_count_at_nucleation": next_nucleation_wake,
        "substeps": rows,
        "wake_on_commit": {
            "active_retained_count_precommit": 1.0,
            "wake_retained_count_precommit": next_nucleation_wake,
            "active_retained_count_postcommit": 0.0,
            "wake_retained_count_postcommit": 1.0 + next_nucleation_wake,
            "wake_retained": 1.0,
            "wake_retained_discarded": 0.0,
            "wake_K_shield_Pa_sqrt_m_postcommit": 2.0e4,
            "total_K_shield_Pa_sqrt_m_postcommit": 2.0e4,
        },
    }


def test_audit_distinguishes_precommit_from_persistent_carryover(tmp_path: Path):
    events = [_event(1, 0.0), _event(2, 1.0)]
    payload = {
        "schema": "persistent_plastic_wake_hazard_event_v918_v1",
        "events": events,
        "active_event_id_at_exit": None,
        "adaptive_dt_used_as_hold_cap": False,
        "nominal_loading_dt_s": 8.4,
        "committed_target_reached": True,
        "post_target_renewals_suppressed": 0,
    }
    (tmp_path / "persistent_wake_event_relaxation_v918.json").write_text(
        json.dumps(payload)
    )
    (tmp_path / "absolute_hazard_event_relaxation_v917.json").write_text(
        json.dumps(payload)
    )

    out = audit_case_v918(tmp_path, 700.0, 10.0)
    assert out["persistent_wake_commit_gate_passed"]
    assert out["retained_state_conservation_passed"]
    assert out["wake_shielding_observed"]
    assert out["wake_carryover_to_later_event_observed"]
    assert out["no_uncommitted_trial_at_exit"]
    assert not out["adaptive_dt_used_as_hold_cap"]
    assert out["interpretation"] == (
        "persistent_wake_conserved_and_shielding_carries_between_events"
    )
