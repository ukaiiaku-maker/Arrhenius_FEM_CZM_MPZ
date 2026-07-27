from __future__ import annotations

import argparse
import json
import math
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
TEMPERATURES = [
    300,
    400,
    500,
    600,
    700,
    800,
    850,
    900,
    950,
    1000,
    1050,
    1100,
    1150,
    1200,
]
MANIFEST_NAME = "persistent_site_production_manifest_v10_0_5_18_3_3.json"
EVENTS_NAME = "stochastic_geometry_events_v10_0_5_16.json"
CORRIDOR_NAME = "compact_corridor_mesh_v91852.json"
SWEPT_POLICY = "fixed_physical_radius_swept_capsule_same_size_law"
PHYSICAL_RADIUS_M = 330.0e-6
FAILURE_SIGNATURES = [
    "STARTUP_FAILURE",
    "found no corridor",
    "Traceback (most recent call last)",
    "max_internal_steps",
    "proposal budget",
    "EMISSION_MAX_EVENTS_PER_HALF_STEP",
    "adaptive-CZM rejected",
    "v91854_quality_veto",
    "v91856_quality_veto",
    "child_area_ratio<",
    "unreliable signed two-channel",
    "Killed",
    "Abort trap",
]


def _load_json(path: pathlib.Path) -> Any:
    return json.loads(path.read_text())


def _close(value: Any, expected: float, *, atol: float = 1.0e-15) -> bool:
    try:
        observed = float(value)
    except (TypeError, ValueError):
        return False
    return math.isclose(observed, expected, rel_tol=1.0e-12, abs_tol=atol)


def _case_record(
    root: pathlib.Path,
    option: str,
    temperature: int,
    target_extension_um: float,
) -> dict[str, Any]:
    case = root / option / f"T{temperature:04d}K"
    manifest_path = case / MANIFEST_NAME
    events_path = case / EVENTS_NAME
    corridor_path = case / CORRIDOR_NAME
    console_path = case / "console.log"

    console = (
        console_path.read_text(errors="replace")
        if console_path.is_file()
        else ""
    )
    signatures = [item for item in FAILURE_SIGNATURES if item in console]
    steps = [
        int(value)
        for value in re.findall(
            rf"\[T={temperature}K\]\s+step\s+(\d+)",
            console,
        )
    ]

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
                if bool(row.get("inserted", False))
                and float(row.get("moved_m", 0.0)) > 0.0
            ]
            rejected = [
                row
                for row in events
                if not bool(row.get("inserted", False))
            ]
            accepted_events = len(accepted)
            rejected_events = len(rejected)
            extension_um = (
                sum(float(row.get("moved_m", 0.0)) for row in accepted)
                * 1.0e6
            )
        except Exception as exc:
            signatures.append(f"event parse error: {exc}")

    corridor: dict[str, Any] = {}
    if corridor_path.is_file():
        try:
            corridor = dict(_load_json(corridor_path) or {})
        except Exception as exc:
            signatures.append(f"corridor parse error: {exc}")

    corridor_target_um = corridor.get("corridor_target_extension_um")
    center_gap_um = corridor.get("center_gap_um")
    qmin = corridor.get("minimum_initial_triangle_quality")
    h_over_lpz = corridor.get("maximum_sampled_hbar_tip_over_L_pz")
    physical_provider = corridor.get("production_physical_provider_detected")
    swept_active = corridor.get("swept_physical_refinement_active")
    swept_policy = corridor.get("swept_physical_refinement_policy")
    mesh_provider = corridor.get("mesh_provider")
    radius_m = corridor.get("production_refinement_radius_m")

    completion_errors: list[str] = []
    if manifest.get("run_completed_without_exception") is True:
        if manifest.get("point_release") != "10.0.5.18.3.3":
            completion_errors.append("wrong point release")
        if extension_um < target_extension_um - 1.0e-3:
            completion_errors.append(
                f"extension {extension_um:.6f} um below target {target_extension_um:.6f} um"
            )
        if rejected_events:
            completion_errors.append(
                f"{rejected_events} rejected geometry event(s)"
            )
        if corridor_target_um != target_extension_um:
            completion_errors.append(
                f"corridor target {corridor_target_um!r} does not equal {target_extension_um} um"
            )
        if corridor.get("full_requested_corridor_covered") is not True:
            completion_errors.append("requested corridor not fully covered")
        if corridor.get("center_gap_requirement_passed") is not True:
            completion_errors.append("center-gap requirement failed")
        if center_gap_um is None or float(center_gap_um) > 100.0 + 1.0e-12:
            completion_errors.append(f"center gap {center_gap_um!r} exceeds 100 um")
        if qmin is None or float(qmin) < 0.035:
            completion_errors.append(f"triangle quality {qmin!r} below 0.035")
        if h_over_lpz is None or float(h_over_lpz) > 0.25:
            completion_errors.append(
                f"h_tip/L_pz {h_over_lpz!r} exceeds 0.25"
            )
        if corridor.get("tip_h_over_da_enforced_as_veto") is not False:
            completion_errors.append("h_tip/da incorrectly enforced as veto")
        if physical_provider is not True:
            completion_errors.append("production physical-refinement provider not detected")
        if swept_active is not True:
            completion_errors.append("swept physical-refinement capsule not active")
        if swept_policy != SWEPT_POLICY:
            completion_errors.append(
                f"wrong swept physical-refinement policy {swept_policy!r}"
            )
        if mesh_provider != SWEPT_POLICY:
            completion_errors.append(f"wrong mesh provider {mesh_provider!r}")
        if not _close(radius_m, PHYSICAL_RADIUS_M):
            completion_errors.append(
                f"physical refinement radius {radius_m!r} is not 330 um"
            )
        if corridor.get("constitutive_physics_changed") is not False:
            completion_errors.append("corridor reports constitutive physics change")
        if signatures:
            completion_errors.append("failure signature present")

    if manifest.get("run_completed_without_exception") is True and not completion_errors:
        status = "SUCCESS"
    elif manifest_path.is_file() or signatures or completion_errors:
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
        "target_extension_um": target_extension_um,
        "run_completed_without_exception": manifest.get(
            "run_completed_without_exception"
        ),
        "point_release": manifest.get("point_release"),
        "corridor_target_extension_um": corridor_target_um,
        "corridor_center_gap_um": center_gap_um,
        "corridor_minimum_triangle_quality": qmin,
        "corridor_maximum_h_tip_over_L_pz": h_over_lpz,
        "production_physical_provider_detected": physical_provider,
        "swept_physical_refinement_active": swept_active,
        "swept_physical_refinement_policy": swept_policy,
        "production_refinement_radius_m": radius_m,
        "failure_signatures": signatures,
        "completion_errors": completion_errors,
        "case_directory": str(case),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("campaign_root", type=pathlib.Path)
    parser.add_argument("--target-extension-um", type=float, default=400.0)
    parser.add_argument("--json", action="store_true")
    parser.add_argument("--strict-complete", action="store_true")
    args = parser.parse_args()

    if args.target_extension_um <= 0.0:
        parser.error("--target-extension-um must be positive")

    root = args.campaign_root.resolve()
    records = [
        _case_record(
            root,
            option,
            temperature,
            args.target_extension_um,
        )
        for option in OPTIONS
        for temperature in TEMPERATURES
    ]
    counts = {
        status: sum(record["status"] == status for record in records)
        for status in (
            "SUCCESS",
            "FAILED",
            "RUNNING_OR_INCOMPLETE",
            "NOT_STARTED",
        )
    }
    payload = {
        "campaign_root": str(root),
        "target_extension_um": args.target_extension_um,
        "expected_cases": len(records),
        "counts": counts,
        "records": records,
    }

    if args.json:
        print(json.dumps(payload, indent=2))
    else:
        print(f"Campaign: {root}")
        print(f"Target extension: {args.target_extension_um:.6f} um")
        print(
            "Completion: "
            f"success={counts['SUCCESS']}/56 "
            f"failed={counts['FAILED']} "
            f"running_or_incomplete={counts['RUNNING_OR_INCOMPLETE']} "
            f"not_started={counts['NOT_STARTED']}"
        )
        print()
        print(
            "class     T_K   status                  step   events  reject  extension_um  gap_um  swept"
        )
        for record in records:
            gap = record["corridor_center_gap_um"]
            gap_text = "-" if gap is None else f"{float(gap):.3f}"
            swept_text = (
                "yes"
                if record["swept_physical_refinement_active"] is True
                else "no"
            )
            print(
                f"{record['material_class']:<9} "
                f"{record['temperature_K']:>4}  "
                f"{record['status']:<22} "
                f"{record['largest_printed_FEM_step']:>6} "
                f"{record['accepted_geometry_events']:>7} "
                f"{record['rejected_geometry_events']:>7} "
                f"{record['committed_extension_um']:>13.6f} "
                f"{gap_text:>8} "
                f"{swept_text:>6}"
            )
            issues = [
                *record["failure_signatures"],
                *record["completion_errors"],
            ]
            if issues:
                print("  issues: " + "; ".join(issues))

    if args.strict_complete:
        complete = (
            counts["SUCCESS"] == len(records)
            and counts["FAILED"] == 0
            and counts["RUNNING_OR_INCOMPLETE"] == 0
            and counts["NOT_STARTED"] == 0
        )
        return 0 if complete else 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
