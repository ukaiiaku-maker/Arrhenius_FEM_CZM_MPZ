from __future__ import annotations

from arrhenius_fracture import pf_theta0_first_passage_event_audit_v10051840 as event_audit
from arrhenius_fracture.pf_theta0_stochastic_fingerprint_v10051840 import NOT_MATERIALIZED


def test_empty_record_reports_not_materialized_for_every_field():
    record = event_audit.empty_first_passage_event_audit()
    payload = record.as_dict()
    assert payload["schema"] == event_audit.MODEL_ID
    fixed_fields = {k: v for k, v in payload.items() if k not in ("schema", "context")}
    assert all(v == NOT_MATERIALIZED for v in fixed_fields.values()), fixed_fields


def test_context_is_free_form_and_does_not_affect_the_fixed_contract():
    record = event_audit.empty_first_passage_event_audit(step=177, physical_time_s=1486.8)
    payload = record.as_dict()
    assert payload["context"] == {"step": 177, "physical_time_s": 1486.8}
    assert payload["threshold_action"] == NOT_MATERIALIZED


def test_record_accepts_real_values_for_a_subset_of_fields():
    record = event_audit.FirstPassageEventAudit(
        threshold_action=0.87,
        hazard_rng_digest="abc123",
        proposed_physical_extension_m=5.0e-6,
        commit_or_veto="committed",
    )
    payload = record.as_dict()
    assert payload["threshold_action"] == 0.87
    assert payload["commit_or_veto"] == "committed"
    # Untouched fields remain honestly NOT_MATERIALIZED, not defaulted to
    # zero or an empty string.
    assert payload["renewed_threshold_action"] == NOT_MATERIALIZED
