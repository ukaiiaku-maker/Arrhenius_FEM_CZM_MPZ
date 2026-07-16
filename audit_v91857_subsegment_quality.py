#!/usr/bin/env python3
"""Certify v9.18.5.6 geometry at physical-event and recursive-subsegment levels."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any


def _fail(message: str) -> None:
    raise RuntimeError(message)


def _load_json(path: Path) -> Any:
    if not path.exists():
        _fail(f"required file missing: {path}")
    return json.loads(path.read_text())


def _close(a: float, b: float, *, atol: float = 1.0e-12, rtol: float = 1.0e-8) -> bool:
    return math.isclose(float(a), float(b), abs_tol=atol, rel_tol=rtol)


def certify_case(case_dir: Path, target_um: float, da_um: float) -> dict[str, Any]:
    case_dir = Path(case_dir)
    audit_path = case_dir / "explicit_quality_wrapper_chain_v91856.json"
    audit = _load_json(audit_path)

    logs = sorted(case_dir.glob("czm_*/czm_advance_log.json"))
    if len(logs) != 1:
        _fail(f"expected exactly one CZM advance log under {case_dir}, found {len(logs)}")
    advance_path = logs[0]
    advance = _load_json(advance_path)
    if not isinstance(advance, list):
        _fail(f"CZM advance log is not a list: {advance_path}")

    if not audit.get("run_completed_without_exception", False):
        _fail(f"quality wrapper reports runtime failure: {audit_path}")
    if audit.get("quality_vetoes"):
        _fail(f"production quality vetoes present: {len(audit['quality_vetoes'])}")
    if audit.get("consecutive_veto_abort") is not None:
        _fail("consecutive-veto abort was recorded")

    accepted = list(audit.get("accepted_events", []))
    if not accepted:
        _fail("quality audit contains no accepted transactions")
    if not all(bool(row.get("accepted", False)) for row in accepted):
        _fail("quality audit contains a nonaccepted transaction")

    committed = [
        row for row in advance
        if bool(row.get("v916_committed", False))
        and row.get("physical_event_id") is not None
    ]
    if not committed:
        _fail("CZM log contains no committed physical subsegments")
    if len(accepted) != len(committed):
        _fail(
            "quality/subsegment count mismatch: "
            f"accepted_quality_transactions={len(accepted)} committed_subsegments={len(committed)}"
        )

    expected_physical = int(round(float(target_um) / float(da_um)))
    if expected_physical <= 0:
        _fail("expected physical event count is not positive")

    event_ids = sorted({int(row["physical_event_id"]) for row in committed})
    expected_ids = list(range(1, expected_physical + 1))
    if event_ids != expected_ids:
        _fail(f"physical event IDs mismatch: got={event_ids} expected={expected_ids}")

    da_m = float(da_um) * 1.0e-6
    target_m = float(target_um) * 1.0e-6
    event_rows: list[dict[str, Any]] = []
    final_markers: list[dict[str, Any]] = []

    for event_id in event_ids:
        segs = [row for row in committed if int(row["physical_event_id"]) == event_id]
        declared = {int(row.get("physical_subsegment_count", -1)) for row in segs}
        if len(declared) != 1:
            _fail(f"event {event_id}: inconsistent declared subsegment counts {declared}")
        ndecl = declared.pop()
        if ndecl != len(segs):
            _fail(f"event {event_id}: actual subsegments={len(segs)} declared={ndecl}")

        indices = sorted(int(row.get("physical_subsegment_index", -1)) for row in segs)
        if indices != list(range(ndecl)):
            _fail(f"event {event_id}: invalid subsegment indices {indices}")

        length_m = sum(float(row.get("length_m", 0.0)) for row in segs)
        if not _close(length_m, da_m, atol=max(1.0e-13, 1.0e-8 * da_m)):
            _fail(
                f"event {event_id}: committed subsegment sum={length_m:.16e} "
                f"does not equal physical da={da_m:.16e}"
            )

        final = max(segs, key=lambda row: int(row["physical_subsegment_index"]))
        if not bool(final.get("v91856_quality_gate_passed", False)):
            _fail(f"event {event_id}: final subsegment lacks v9.18.5.6 quality marker")
        final_markers.append(final)

        event_rows.append({
            "physical_event_id": event_id,
            "subsegment_count": ndecl,
            "committed_length_um": length_m * 1.0e6,
            "min_triangle_quality": float(final["v91856_min_triangle_quality"]),
            "min_child_area_ratio": float(final["v91856_min_child_area_ratio"]),
            "active_tip_h_over_da": float(final["v91856_active_tip_h_over_da"]),
            "resolution_warning": bool(final.get("v91856_resolution_warning", False)),
        })

    total_length = sum(float(row.get("length_m", 0.0)) for row in committed)
    if not _close(total_length, target_m, atol=max(1.0e-13, 1.0e-8 * target_m)):
        _fail(
            f"total committed subsegment length={total_length:.16e} "
            f"does not equal target={target_m:.16e}"
        )

    qmin = min(float(row["min_triangle_quality"]) for row in accepted)
    amin = min(float(row["min_child_area_ratio"]) for row in accepted)
    qfloor = max(float(row["triangle_quality_floor"]) for row in accepted)
    afloor = max(float(row["child_area_ratio_floor"]) for row in accepted)
    if qmin < qfloor:
        _fail(f"accepted transaction qmin={qmin:.16e} below floor={qfloor:.16e}")
    if amin < afloor:
        _fail(f"accepted transaction area ratio={amin:.16e} below floor={afloor:.16e}")

    result = {
        "schema": "subsegment_aware_quality_certification_v91857_v1",
        "certified": True,
        "case_dir": str(case_dir),
        "physical_event_count": len(event_ids),
        "expected_physical_event_count": expected_physical,
        "committed_subsegment_count": len(committed),
        "accepted_quality_transaction_count": len(accepted),
        "final_physical_event_quality_marker_count": len(final_markers),
        "total_committed_extension_um": total_length * 1.0e6,
        "physical_da_um": float(da_um),
        "minimum_accepted_triangle_quality": qmin,
        "triangle_quality_floor": qfloor,
        "minimum_accepted_child_area_ratio": amin,
        "child_area_ratio_floor": afloor,
        "all_resolution_warning_count": len(audit.get("resolution_warnings", [])),
        "final_physical_event_resolution_warning_count": sum(
            bool(row.get("v91856_resolution_warning", False)) for row in final_markers
        ),
        "quality_veto_count": 0,
        "event_rows": event_rows,
        "source_quality_audit": str(audit_path),
        "source_czm_advance_log": str(advance_path),
        "constitutive_physics_changed": False,
    }
    out = case_dir / "subsegment_aware_quality_certification_v91857.json"
    out.write_text(json.dumps(result, indent=2))
    return result


def certify_campaign(case_root: Path, target_um: float, da_um: float) -> list[dict[str, Any]]:
    root = Path(case_root)
    summary_path = root / "v9_13_campaign_summary.json"
    rows = _load_json(summary_path)
    if not isinstance(rows, list):
        _fail(f"campaign summary is not a list: {summary_path}")

    failed = [
        row for row in rows
        if int(row.get("subprocess_returncode", row.get("returncode", 1)) or 0) != 0
    ]
    if failed:
        _fail(
            "inner solver failures: "
            + repr([(row.get("class"), row.get("subprocess_returncode"), row.get("log")) for row in failed])
        )

    results = []
    for row in rows:
        if not bool(row.get("target_completed", False)):
            _fail(f"campaign row is not target-complete: {row.get('case_dir')}")
        results.append(certify_case(Path(row["case_dir"]), target_um, da_um))
    return results


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("case_root", type=Path)
    parser.add_argument("target_um", type=float)
    parser.add_argument("da_um", type=float)
    args = parser.parse_args()

    for result in certify_campaign(args.case_root, args.target_um, args.da_um):
        print(
            "QUALITY CERTIFIED "
            f"{Path(result['case_dir']).name}: "
            f"physical_events={result['physical_event_count']} "
            f"committed_subsegments={result['committed_subsegment_count']} "
            f"accepted_quality_transactions={result['accepted_quality_transaction_count']} "
            f"final_markers={result['final_physical_event_quality_marker_count']} "
            f"qmin={result['minimum_accepted_triangle_quality']:.6g} "
            f"area_ratio_min={result['minimum_accepted_child_area_ratio']:.6g} "
            f"final_resolution_warnings={result['final_physical_event_resolution_warning_count']} "
            f"all_resolution_warnings={result['all_resolution_warning_count']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
