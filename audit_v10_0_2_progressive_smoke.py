#!/usr/bin/env python3
"""Fail-closed audit for the first 5 um v10.0.2 progressive event.

Unlike the v10.0 audit, transactional trial-damage rejections are allowed when
all rejected attempts are restored and subsequently accepted at reduced
physical time.  The audit still requires exactly one committed physical
checkpoint, no geometry-quality veto, no duplicate MPZ translation, and exact
target-stop accounting.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

SCHEMA = "kinetic_campaign_czm_progressive_smoke_v10_0_2"


def load(path: Path) -> Any:
    return json.loads(path.read_text())


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit(message)


def find_one(root: Path, name: str) -> Path:
    direct = root / name
    if direct.exists():
        return direct
    matches = sorted(root.rglob(name))
    require(len(matches) == 1, f"expected exactly one {name}; found {len(matches)}")
    return matches[0]


def committed_log_rows(root: Path) -> tuple[Path, list[dict[str, Any]]]:
    candidates = sorted(root.rglob("czm_advance_log.json"))
    require(len(candidates) == 1, f"expected one czm_advance_log.json; found {len(candidates)}")
    payload = load(candidates[0])
    rows = payload if isinstance(payload, list) else payload.get("rows", [])
    require(isinstance(rows, list) and rows, "CZM advance log is empty")
    return candidates[0], rows


def nondecreasing(values: list[float], tol: float = 1.0e-12) -> bool:
    return all(b + tol >= a for a, b in zip(values, values[1:]))


def audit(root: Path, target_um: float = 5.0) -> dict[str, Any]:
    root = Path(root)
    progressive_path = find_one(
        root, "kinetic_campaign_czm_progressive_2d_v10_0_2.json"
    )
    model_path = find_one(root, "kinetic_campaign_czm_v10_0_audit.json")
    quality_path = find_one(root, "explicit_quality_wrapper_chain_v91856.json")
    progressive = load(progressive_path)
    model = load(model_path)
    quality = load(quality_path)

    require(progressive.get("full_progressive_trial_loop_active") is True,
            "dedicated v10.0.2 progressive loop was not active")
    require(model.get("full_progressive_trial_loop_active") is True,
            "model audit did not certify progressive mode")
    require(model.get("point_release") == "10.0.2",
            f"wrong point release in model audit: {model.get('point_release')}")
    require(model.get("rejected_damage_retry_integrated") is True,
            "reduced-dt retry is not certified in model audit")
    require(model.get("unused_event_time_carry_integrated") is True,
            "unused-time carry is not certified in model audit")
    require(model.get("same_load_re_equilibration_after_commit") is True,
            "same-load re-equilibration is not certified")
    require(model.get("dot_ep_transactionally_restored") is True,
            "dot_ep rollback is not certified")

    require(progressive.get("rejected_step_retry_active") is True,
            "runtime retry controller was not active")
    require(progressive.get("unused_time_carry_active") is True,
            "runtime unused-time carry was not active")
    require(progressive.get("same_load_re_equilibration_after_commit") is True,
            "runtime same-load re-equilibration flag is absent")
    require(progressive.get("one_topology_event_per_equilibrium_state") is True,
            "topology/equilibrium lifecycle gate is not active")
    require(int(progressive.get("trial_insertions", 0)) == 1,
            f"expected one trial insertion; got {progressive.get('trial_insertions')}")
    require(int(progressive.get("committed_events", 0)) == 1,
            f"expected one committed event; got {progressive.get('committed_events')}")
    require(int(progressive.get("max_commits_in_outer_interval", 0)) <= 1,
            "more than one checkpoint was committed in one outer interval")
    require(int(progressive.get("full_rollbacks", 0)) == 0,
            "unexpected full topology rollback occurred")
    require(float(progressive.get("mpz_advance_on_commit_m", math.nan)) == 0.0,
            "MPZ was translated again at outer commitment")

    records = progressive.get("records", [])
    require(isinstance(records, list) and records, "progressive step records are empty")
    event_ids = {int(row.get("trial_event_id", -1)) for row in records}
    require(len(event_ids) == 1 and min(event_ids) >= 1,
            f"expected one positive trial event id; found {sorted(event_ids)}")

    damages = [float(row.get("trial_cohesive_damage", math.nan)) for row in records]
    require(all(math.isfinite(x) and -1e-12 <= x <= 1.0 + 1e-12 for x in damages),
            "trial damage contains invalid values")
    require(nondecreasing(damages), "accepted trial cohesive damage is not monotonic")
    require(math.isclose(damages[-1], 1.0, rel_tol=0.0, abs_tol=1.0e-10),
            f"final trial damage is not one: {damages[-1]}")
    require(bool(records[-1].get("fired", False)),
            "final progressive record did not commit")
    require(sum(bool(row.get("fired", False)) for row in records) == 1,
            "more than one committed checkpoint appears in accepted records")

    runtime_rejections = int(progressive.get("damage_rejections", 0))
    record_retries = sum(int(row.get("retry_count", 0)) for row in records)
    require(runtime_rejections == record_retries,
            f"retry accounting mismatch: runtime={runtime_rejections}, records={record_retries}")
    require(all(float(row.get("trial_requested_dt_s", 0.0)) > 0.0 for row in records),
            "accepted lifecycle record has non-positive requested dt")

    micro = float(records[-1].get("micro_advance_total_m", math.nan))
    checkpoint = float(records[-1].get("checkpoint_committed_total_m", math.nan))
    target_m = float(target_um) * 1.0e-6
    tol_m = max(1.0e-12, 1.0e-8 * target_m)
    require(math.isclose(micro, target_m, rel_tol=0.0, abs_tol=tol_m),
            f"micro advance {micro:.16g} m does not equal target {target_m:.16g} m")
    require(math.isclose(checkpoint, target_m, rel_tol=0.0, abs_tol=tol_m),
            f"checkpoint total {checkpoint:.16g} m does not equal target {target_m:.16g} m")
    require(float(records[-1].get("dt_unused_s", 0.0)) >= 0.0,
            "negative unused physical time")

    require(quality.get("run_completed_without_exception") is True,
            f"quality wrapper recorded an exception: {quality.get('runtime_error')}")
    require(not quality.get("quality_vetoes", []), "geometry quality vetoes were recorded")
    require(quality.get("consecutive_veto_abort") is None,
            "consecutive geometry-veto abort was recorded")
    accepted = quality.get("accepted_events", [])
    require(isinstance(accepted, list) and accepted,
            "quality wrapper has no accepted topology transactions")
    for row in accepted:
        require(bool(row.get("accepted", False)), "quality transaction was not accepted")
        require(float(row["min_triangle_quality"]) + 1e-15 >= float(row["triangle_quality_floor"]),
                "accepted topology falls below triangle-quality floor")
        require(float(row["min_child_area_ratio"]) + 1e-15 >= float(row["child_area_ratio_floor"]),
                "accepted topology falls below child-area floor")
        require(not row.get("issues", []),
                f"accepted quality row has issues: {row.get('issues')}")

    log_path, log_rows = committed_log_rows(root)
    committed = [row for row in log_rows if str(row.get("status", "committed")) == "committed"]
    require(committed, "no committed CZM log rows")
    log_event_ids = {
        int(row.get("physical_event_index", row.get("event_index", -1)))
        for row in committed
    }
    require(len(log_event_ids) == 1,
            f"CZM log contains multiple physical events: {sorted(log_event_ids)}")
    length = sum(float(row.get("length_m", 0.0)) for row in committed)
    require(math.isclose(length, target_m, rel_tol=0.0, abs_tol=tol_m),
            f"committed CZM subsegments sum to {length:.16g} m, expected {target_m:.16g} m")
    require(all(float(row.get("damage", 0.0)) >= 1.0 - 1.0e-10 for row in committed),
            "one or more committed CZM subsegments are not fully damaged")
    require(all(float(row.get("clock", 1.0)) >= 1.0 - 1.0e-10 for row in committed),
            "one or more committed CZM subsegments do not have clock=1")
    require(all(float(row.get("mpz_advance_on_commit_m", 0.0)) == 0.0 for row in committed),
            "CZM log reports a second MPZ translation at commitment")

    payload = {
        "schema": SCHEMA,
        "passed": True,
        "root": str(root.resolve()),
        "target_extension_um": float(target_um),
        "progressive_record_count": len(records),
        "trial_event_id": next(iter(event_ids)),
        "trial_insertions": 1,
        "committed_events": 1,
        "damage_rejections_retried": runtime_rejections,
        "micro_advance_total_m": micro,
        "checkpoint_committed_total_m": checkpoint,
        "dt_unused_s": float(records[-1].get("dt_unused_s", 0.0)),
        "carried_time_s": float(progressive.get("carried_time_s", 0.0)),
        "accepted_quality_transactions": len(accepted),
        "minimum_triangle_quality": min(float(row["min_triangle_quality"]) for row in accepted),
        "minimum_child_area_ratio": min(float(row["min_child_area_ratio"]) for row in accepted),
        "committed_czm_subsegments": len(committed),
        "committed_czm_length_m": length,
        "progressive_audit": str(progressive_path),
        "quality_audit": str(quality_path),
        "czm_advance_log": str(log_path),
        "no_duplicate_mpz_translation": True,
        "one_topology_event_only": True,
        "accepted_damage_monotonic": True,
        "transactional_retry_accounting": True,
    }
    out = root / "progressive_one_segment_smoke_v10_0_2.json"
    out.write_text(json.dumps(payload, indent=2))
    return payload


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("output")
    p.add_argument("--target-um", type=float, default=5.0)
    args = p.parse_args(argv)
    payload = audit(Path(args.output), args.target_um)
    print(
        "V10.0.2 PROGRESSIVE SMOKE CERTIFIED: "
        f"event={payload['trial_event_id']} "
        f"extension_um={payload['target_extension_um']:.6f} "
        f"retried_rejections={payload['damage_rejections_retried']} "
        f"quality_transactions={payload['accepted_quality_transactions']} "
        f"qmin={payload['minimum_triangle_quality']:.6f} "
        f"area_ratio_min={payload['minimum_child_area_ratio']:.6f}"
    )


if __name__ == "__main__":
    main()
