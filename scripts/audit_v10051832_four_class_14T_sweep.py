from __future__ import annotations

import argparse
import json
import pathlib
import re
import sys
from typing import Any

OPTIONS = {
    "v913_paper_peak01_0242980_persistent_sites": "peak",
    "v913_paper_dbtt01_0202500_persistent_sites": "DBTT",
    "v913_paper_weakT01_0129902_persistent_sites": "weakT",
    "v913_paper_ceramic01_0077080_persistent_sites": "ceramic",
}
TEMPERATURES = [300, 400, 500, 600, 700, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200]
MANIFEST_NAME = "persistent_site_production_manifest_v10_0_5_18_3_2.json"
EVENTS_NAME = "stochastic_geometry_events_v10_0_5_16.json"
FAILURE_SIGNATURES = [
    "Traceback (most recent call last)",
    "max_internal_steps",
    "proposal budget",
    "EMISSION_MAX_EVENTS_PER_HALF_STEP",
    "adaptive-CZM rejected",
    "unreliable signed two-channel",
    "Killed",
    "Abort trap",
]


def _load_json(path: pathlib.Path) -> Any:
    return json.loads(path.read_text())


def _case_record(root: pathlib.Path, option: str, temperature: int) -> dict[str, Any]:
    case = root / option / f"T{temperature:04d}K"
    manifest_path = case / MANIFEST_NAME
    events_path = case / EVENTS_NAME
    console_path = case / "console.log"
    console = console_path.read_text(errors="replace") if console_path.is_file() else ""
    signatures = [item for item in FAILURE_SIGNATURES if item in console]
    steps = [int(value) for value in re.findall(rf"\[T={temperature}K\]\s+step\s+(\d+)", console)]

    manifest: dict[str, Any] = {}
    if manifest_path.is_file():
        try:
            manifest = dict(_load_json(manifest_path) or {})
        except Exception as exc:
            signatures.append(f"manifest parse error: {exc}")

    accepted_events = 0
    rejected_events = 0
    extension_um = 0.0
    if events_path.is_file():
        try:
            events = list(_load_json(events_path) or [])
            accepted = [
                row
                for row in events
                if bool(row.get("inserted", False)) and float(row.get("moved_m", 0.0)) > 0.0
            ]
            rejected = [row for row in events if not bool(row.get("inserted", False))]
            accepted_events = len(accepted)
            rejected_events = len(rejected)
            extension_um = sum(float(row.get("moved_m", 0.0)) for row in accepted) * 1.0e6
        except Exception as exc:
            signatures.append(f"event parse error: {exc}")

    if manifest.get("run_completed_without_exception") is True:
        status = "SUCCESS"
    elif manifest_path.is_file() or signatures:
        status = "FAILED"
    elif console_path.is_file():
        status = "RUNNING_OR_INCOMPLETE"
    else:
        status = "NOT_STARTED"

    return {
        "parameter_option": option,
        "material_class": OPTIONS[option],
        "temperature_K": temperature,
        "status": status,
        "largest_printed_FEM_step": max(steps, default=0),
        "accepted_geometry_events": accepted_events,
        "rejected_geometry_events": rejected_events,
        "committed_extension_um": extension_um,
        "run_completed_without_exception": manifest.get("run_completed_without_exception"),
        "point_release": manifest.get("point_release"),
        "failure_signatures": signatures,
        "case_directory": str(case),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("campaign_root", type=pathlib.Path)
    parser.add_argument("--json", action="store_true")
    parser.add_argument("--strict-complete", action="store_true")
    args = parser.parse_args()

    root = args.campaign_root.resolve()
    records = [
        _case_record(root, option, temperature)
        for option in OPTIONS
        for temperature in TEMPERATURES
    ]
    counts = {
        status: sum(record["status"] == status for record in records)
        for status in ("SUCCESS", "FAILED", "RUNNING_OR_INCOMPLETE", "NOT_STARTED")
    }
    payload = {
        "campaign_root": str(root),
        "expected_cases": len(records),
        "counts": counts,
        "records": records,
    }

    if args.json:
        print(json.dumps(payload, indent=2))
    else:
        print(f"Campaign: {root}")
        print(
            "Completion: "
            f"success={counts['SUCCESS']}/56 "
            f"failed={counts['FAILED']} "
            f"running_or_incomplete={counts['RUNNING_OR_INCOMPLETE']} "
            f"not_started={counts['NOT_STARTED']}"
        )
        print()
        print("class     T_K   status                  step   events  reject  extension_um")
        for record in records:
            print(
                f"{record['material_class']:<9} "
                f"{record['temperature_K']:>4}  "
                f"{record['status']:<22} "
                f"{record['largest_printed_FEM_step']:>6} "
                f"{record['accepted_geometry_events']:>7} "
                f"{record['rejected_geometry_events']:>7} "
                f"{record['committed_extension_um']:>13.6f}"
            )
            if record["failure_signatures"]:
                print("  failures: " + "; ".join(record["failure_signatures"]))

    if args.strict_complete:
        complete = counts["SUCCESS"] == len(records) and counts["FAILED"] == 0
        return 0 if complete else 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
