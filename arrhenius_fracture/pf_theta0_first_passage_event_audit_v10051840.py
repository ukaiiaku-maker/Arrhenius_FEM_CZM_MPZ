"""Gate 3/4 first-passage / cleavage-event audit schema (FEM_CZM_HANDOFF.md
section 8).

This module defines the record shape a real hazard-crossing event must be
audited with -- threshold and RNG provenance, event-length draw and
factor, proposed physical extension and endpoint, energy-gate inputs and
result, cohesive/geometry transaction identity, commit-or-veto result,
rollback identity, and the renewed threshold/RNG digest after a
successful commit.

Deliberately schema-only for now: this is defined ahead of the real
production run (which is still blocked on the frozen PF reference
artifacts, see CLAUDE_PROGRESS.md's "EXTERNAL BLOCKER") so the Gate 3/4
comparison tooling has a stable target to build against, but it is not
yet wired into `sharp_front.run_2d`'s actual event-commit path (the
`info['fired']` block, which several other modules AST/text-patch -- see
the "sharp_front.py architecture" notes in CLAUDE_PROGRESS.md). Wiring
this in is deferred until a real hazard crossing needs to be audited
end-to-end, per explicit instruction: do not populate this for every
prefracture step, only when an actual crossing occurs.

Every field defaults to NOT_MATERIALIZED. A field must only ever be
filled in with a value that was actually produced by a real event --
never fabricated or guessed to satisfy this schema.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from .pf_theta0_stochastic_fingerprint_v10051840 import NOT_MATERIALIZED

MODEL_ID = "pf_theta0_first_passage_event_audit_v10_0_5_18_4_0"


@dataclass(frozen=True)
class FirstPassageEventAudit:
    schema: str = MODEL_ID

    # Threshold and RNG provenance at the moment of the hazard crossing.
    threshold_action: float | str = NOT_MATERIALIZED
    hazard_rng_digest: str = NOT_MATERIALIZED
    emission_rng_digest: str = NOT_MATERIALIZED

    # Event-length draw and the resulting proposed physical extension.
    event_length_factor: float | str = NOT_MATERIALIZED
    event_length_draw_m: float | str = NOT_MATERIALIZED
    proposed_physical_extension_m: float | str = NOT_MATERIALIZED
    direction_xy: tuple[float, float] | str = NOT_MATERIALIZED
    proposed_endpoint_xy_m: tuple[float, float] | str = NOT_MATERIALIZED

    # Energy gate.
    available_energy_J_per_m: float | str = NOT_MATERIALIZED
    required_energy_J_per_m: float | str = NOT_MATERIALIZED
    energy_gate_result: str = NOT_MATERIALIZED  # "passed" | "failed" | NOT_MATERIALIZED

    # Cohesive/geometry transaction.
    cohesive_transaction_id: str = NOT_MATERIALIZED
    committed_endpoint_xy_m: tuple[float, float] | str = NOT_MATERIALIZED
    commit_or_veto: str = NOT_MATERIALIZED  # "committed" | "vetoed" | NOT_MATERIALIZED
    veto_reason: str = NOT_MATERIALIZED
    rollback_identity: str = NOT_MATERIALIZED

    # Renewal after a successful commit.
    renewed_threshold_action: float | str = NOT_MATERIALIZED
    renewed_hazard_rng_digest: str = NOT_MATERIALIZED
    B_after_renewal: float | str = NOT_MATERIALIZED
    N_em_after_renewal: float | str = NOT_MATERIALIZED

    # Free-form provenance (step index, physical time, git commit, etc.);
    # not part of the fixed contract above, so additions here never
    # require a schema version bump.
    context: dict[str, Any] = field(default_factory=dict)

    def as_dict(self) -> dict[str, Any]:
        payload = {
            "schema": self.schema,
            "threshold_action": self.threshold_action,
            "hazard_rng_digest": self.hazard_rng_digest,
            "emission_rng_digest": self.emission_rng_digest,
            "event_length_factor": self.event_length_factor,
            "event_length_draw_m": self.event_length_draw_m,
            "proposed_physical_extension_m": self.proposed_physical_extension_m,
            "direction_xy": self.direction_xy,
            "proposed_endpoint_xy_m": self.proposed_endpoint_xy_m,
            "available_energy_J_per_m": self.available_energy_J_per_m,
            "required_energy_J_per_m": self.required_energy_J_per_m,
            "energy_gate_result": self.energy_gate_result,
            "cohesive_transaction_id": self.cohesive_transaction_id,
            "committed_endpoint_xy_m": self.committed_endpoint_xy_m,
            "commit_or_veto": self.commit_or_veto,
            "veto_reason": self.veto_reason,
            "rollback_identity": self.rollback_identity,
            "renewed_threshold_action": self.renewed_threshold_action,
            "renewed_hazard_rng_digest": self.renewed_hazard_rng_digest,
            "B_after_renewal": self.B_after_renewal,
            "N_em_after_renewal": self.N_em_after_renewal,
        }
        if self.context:
            payload["context"] = dict(self.context)
        return payload


def empty_first_passage_event_audit(**context: Any) -> FirstPassageEventAudit:
    """A schema-conformant record with every field NOT_MATERIALIZED,
    optionally tagged with free-form context (step, physical time, ...).
    Fill in only the fields a real event actually produced."""
    return FirstPassageEventAudit(context=dict(context))


__all__ = [
    "MODEL_ID",
    "FirstPassageEventAudit",
    "empty_first_passage_event_audit",
]
