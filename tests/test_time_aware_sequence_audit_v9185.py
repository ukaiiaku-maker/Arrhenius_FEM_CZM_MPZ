from __future__ import annotations

import csv
import json

from arrhenius_fracture.fixed_displacement_sequence_audit_v9185 import (
    write_time_aware_sequence_audit,
)


def test_long_opening_event_is_delayed_not_unstable(tmp_path, monkeypatch):
    monkeypatch.setenv("ARRHENIUS_NOMINAL_LOADING_DT_S", "8.4")
    rows = [{
        "load_event_id": "1",
        "classification": "unstable_same_load_cascade",
        "raw_event_start": "3",
        "raw_event_end": "9",
        "topology_event_count": "7",
    }]
    with (tmp_path / "R_curve_load_events_clustered.csv").open("w", newline="") as fp:
        writer = csv.DictWriter(fp, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    events = {
        "events": [
            {"event_id": i, "relaxation_time_s": t}
            for i, t in zip(range(3, 10), [0.08, 0.02, 0.004, 0.0007, 0.0003, 0.00007, 73.6])
        ]
    }
    (tmp_path / "kinetic_trial_event_relaxation_v916.json").write_text(
        json.dumps(events)
    )
    summary = write_time_aware_sequence_audit(tmp_path)
    assert summary is not None
    assert summary["n_rapid_fixed_displacement_cascades"] == 0
    assert summary["n_delayed_fixed_displacement_multievent_sequences"] == 1
    with (tmp_path / "R_curve_load_events_time_aware_v9185.csv").open() as fp:
        result = next(csv.DictReader(fp))
    assert result["time_aware_classification"] == (
        "delayed_fixed_displacement_multievent_growth"
    )
    assert float(result["event_opening_time_max_s"]) == 73.6
