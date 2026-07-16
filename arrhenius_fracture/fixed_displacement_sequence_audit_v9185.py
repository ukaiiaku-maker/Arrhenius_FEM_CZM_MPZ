"""Time-aware classification of clustered fixed-displacement crack events."""
from __future__ import annotations

import csv
import json
import math
import os
from pathlib import Path
from typing import Any


def _float(value: Any, default: float = 0.0) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if math.isfinite(out) else default


def write_time_aware_sequence_audit(case_dir: str | Path) -> dict[str, Any] | None:
    case = Path(case_dir)
    clustered_path = case / "R_curve_load_events_clustered.csv"
    event_path = case / "kinetic_trial_event_relaxation_v916.json"
    if not clustered_path.exists() or not event_path.exists():
        return None

    try:
        clustered = list(csv.DictReader(clustered_path.open()))
        event_payload = json.loads(event_path.read_text())
    except Exception:
        return None
    events = event_payload.get("events", []) if isinstance(event_payload, dict) else []
    event_time = {
        int(e.get("event_id", -1)): _float(e.get("relaxation_time_s"), 0.0)
        for e in events if isinstance(e, dict)
    }
    nominal_dt = max(_float(os.environ.get("ARRHENIUS_NOMINAL_LOADING_DT_S", 8.4), 8.4), 0.0)
    rapid_limit = max(_float(
        os.environ.get("ARRHENIUS_RAPID_CASCADE_MAX_EVENT_TIME_S", nominal_dt),
        nominal_dt,
    ), 0.0)

    out_rows = []
    n_rapid = 0
    n_delayed = 0
    rapid_events = 0
    delayed_events = 0
    for row in clustered:
        out = dict(row)
        start = int(_float(row.get("raw_event_start"), -1))
        end = int(_float(row.get("raw_event_end"), start))
        times = [event_time.get(i, 0.0) for i in range(start, end + 1)] if start >= 0 else []
        count = int(_float(row.get("topology_event_count"), 0))
        out["event_opening_time_sum_s"] = sum(times)
        out["event_opening_time_max_s"] = max(times) if times else 0.0
        out["rapid_cascade_event_time_limit_s"] = rapid_limit
        if count > 1:
            if times and max(times) <= rapid_limit:
                out["time_aware_classification"] = "rapid_fixed_displacement_cascade"
                n_rapid += 1
                rapid_events += count
            else:
                out["time_aware_classification"] = "delayed_fixed_displacement_multievent_growth"
                n_delayed += 1
                delayed_events += count
        else:
            out["time_aware_classification"] = "single_topology_event"
        out_rows.append(out)

    out_csv = case / "R_curve_load_events_time_aware_v9185.csv"
    if out_rows:
        fields = list(out_rows[0])
        with out_csv.open("w", newline="") as fp:
            writer = csv.DictWriter(fp, fieldnames=fields)
            writer.writeheader()
            writer.writerows(out_rows)

    summary = {
        "schema": "time_aware_fixed_displacement_sequence_audit_v9185_v1",
        "source_cluster_file": clustered_path.name,
        "source_event_file": event_path.name,
        "rapid_cascade_max_event_time_s": rapid_limit,
        "n_rapid_fixed_displacement_cascades": n_rapid,
        "n_delayed_fixed_displacement_multievent_sequences": n_delayed,
        "topology_events_in_rapid_cascades": rapid_events,
        "topology_events_in_delayed_sequences": delayed_events,
        "classification_rule": (
            "multi-event fixed-displacement group is rapid only when every physical "
            "cohesive-opening time is no greater than the configured rapid limit"
        ),
    }
    (case / "fixed_displacement_sequence_audit_v9185.json").write_text(
        json.dumps(summary, indent=2, default=str)
    )
    return summary


def scan_campaign(root: str | Path) -> list[dict[str, Any]]:
    root = Path(root)
    results = []
    for source in root.rglob("R_curve_load_events_clustered.csv"):
        payload = write_time_aware_sequence_audit(source.parent)
        if payload is not None:
            results.append({"case_dir": str(source.parent), **payload})
    if results:
        (root / "fixed_displacement_sequence_campaign_audit_v9185.json").write_text(
            json.dumps(results, indent=2, default=str)
        )
    return results
