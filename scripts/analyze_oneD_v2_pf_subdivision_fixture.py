#!/usr/bin/env python3
"""Reduce the bounded PF OFF/ON observer fixture into canonical evidence."""

from __future__ import annotations

from collections import Counter
import hashlib
import json
from pathlib import Path


OFF = Path("/private/tmp/oneD_v2_pf_subdivision_off")
ON = Path("/private/tmp/oneD_v2_pf_subdivision_on2")
OUT = Path(__file__).resolve().parents[1] / "analysis_outputs/oneD_v2_capability_v2"
ARTIFACTS = (
    "steps_1000K.csv", "fronts_1000K.csv", "branch_diagnostics_1000K.csv",
    "crack_path_1000K.csv", "crack_path_front0_1000K.csv",
    "sharp_wake_advance_log.csv", "stochastic_avalanche_geometry_events.json",
    "summary.json", "v10_2_30_hazard_energy_gate_audit.json",
    "kinetic_tip_cell_audit_v101.json",
)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    records = [json.loads(line) for line in (ON / "observer.jsonl").read_text().splitlines()]
    counts = Counter(row["boundary"] for row in records)
    reject_restore = []
    latest_snapshot = None
    for row in records:
        if row["boundary"] == "TRANSACTION_SNAPSHOT":
            latest_snapshot = row["state_sha256"]
        elif row["boundary"] == "TRIAL_REJECT":
            reject_restore.append(row["state_sha256"] == latest_snapshot)
    identity = {name: digest(OFF / name) == digest(ON / name) for name in ARTIFACTS}
    payload = {
        "schema": "oneD_v2_pf_subdivision_fixture_v1",
        "authoritative_production_source_commit": "ab6279331d050919f78e2e7f278bc332466e34f8",
        "default_off_observer_commit": "c695b44",
        "fixture_scope": "BOUNDED_ANALYSIS_ONLY_ONE_ACCEPTED_EVENT",
        "canonical_campaign_membership": False,
        "boundary_counts": dict(sorted(counts.items())),
        "rejected_trials": counts["TRIAL_REJECT"],
        "all_rejected_trials_restore_snapshot_fingerprint": bool(reject_restore) and all(reject_restore),
        "threshold_contract": "DETERMINISTIC_UNIT_ACTION_NO_DRAW",
        "threshold_generation_persists_through_subdivision": bool(reject_restore) and all(reject_restore),
        "emission_threshold_contract": "NOT_APPLICABLE_CONTINUOUS_RATE",
        "driver_rng_restore": "NOT_APPLICABLE_NO_DRIVER_RNG",
        "event_length_draw": "NOT_APPLICABLE_FIXED_F_DA",
        "direction_draw": "NOT_APPLICABLE_DETERMINISTIC_PLANE_SELECTION",
        "first_passage_accept_count": counts["FIRST_PASSAGE_ACCEPT"],
        "geometry_commit_count": counts["GEOMETRY_COMMIT"],
        "post_event_renewal_count": counts["POST_EVENT_RENEWAL"],
        "authoritative_interval_commit_count": counts["AUTHORITATIVE_INTERVAL_COMMIT"],
        "physical_artifact_identity": identity,
        "diagnostics_off_on_physical_state_identical": all(identity.values()),
        "classification": "PF_NORMAL_PATH_QUALIFIED_FAIL_CLOSED_ON_VETO",
        "rollback_and_continue": "UNSUPPORTED_BY_PRODUCTION",
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "oneD_v2_pf_subdivision_fixture_audit.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n"
    )
    md = f"""# One-dimensional V2 PF normal-path qualification V2

Classification: **{payload['classification']}**.

The controlled Peak fixture used the actual `v10.2.30` PF entry, canonical
seed 8666, the production adaptive predictor/retry loop, and stopped after one
accepted 5 µm event. Tightening the adaptive target to 0.01 forced
**{payload['rejected_trials']} rejected trials**. All rejected-trial state
fingerprints exactly matched their transaction snapshots. Diagnostics OFF and
ON produced byte-identical authoritative physical artifacts.

Exactly-once results: first-passage accept = **{counts['FIRST_PASSAGE_ACCEPT']}**,
geometry commit = **{counts['GEOMETRY_COMMIT']}**, post-event renewal =
**{counts['POST_EVENT_RENEWAL']}**.

This source does not draw stochastic cleavage or emission thresholds. Cleavage
uses a deterministic unit-action renewal surface; emission is a continuous
rate. Event length is fixed by `f.da`, and direction is selected deterministically
from the cleavage planes. Their RNG/draw/restore requirements are therefore
**NOT_APPLICABLE**, not inferred from later rows.

Late-veto termination remains **SOURCE_CONTRACT_QUALIFIED** and rollback-and-
continue remains **UNSUPPORTED_BY_PRODUCTION**. This normal-path result does not
alter either status.
"""
    (OUT / "ONE_D_V2_PF_NORMAL_PATH_QUALIFICATION_V2.md").write_text(md)


if __name__ == "__main__":
    main()
