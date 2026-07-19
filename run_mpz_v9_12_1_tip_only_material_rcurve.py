#!/usr/bin/env python3
"""v9.12.1 full-field material R-curve runner.

This wrapper preserves the validated v9.12 FEM/CZM/MPZ solver command and fixes
only campaign bookkeeping and publication gating.  A shared temperature summary
is accepted only when it was created or changed by the successful case
invocation being recorded; stale output from another material class is rejected.
"""
from __future__ import annotations

import json
from pathlib import Path
import shutil
from typing import Any

import numpy as np
import pandas as pd

import run_mpz_v9_12_tip_only_material_rcurve as _legacy
from arrhenius_fracture.material_rcurve_audit_v9121 import (
    audit_campaign as _audit_campaign_v9121,
)

_legacy_run_case = _legacy.run_case
_ACTIVE_THETA_DEG = 45.0


def _finite_value(value) -> bool:
    if value is None:
        return False
    try:
        return bool(np.isfinite(float(value)))
    except (TypeError, ValueError):
        return bool(str(value).strip())


def _merge_observables(row: dict[str, Any], summary: dict[str, Any]) -> dict[str, Any]:
    merged = dict(row)
    for key, value in summary.items():
        if _finite_value(value):
            merged[key] = value
    return merged


def _file_signature(path: Path):
    if not path.is_file():
        return None
    stat = path.stat()
    return (int(stat.st_mtime_ns), int(stat.st_size))


def audit_campaign(campaign_root, seed, T_K, classes, bulk_mode="tip_only"):
    payload = _audit_campaign_v9121(
        campaign_root,
        seed,
        T_K,
        classes=classes,
        bulk_mode=bulk_mode,
        theta_deg=_ACTIVE_THETA_DEG,
    )
    # Compatibility fields consumed by the inherited console reporter.
    payload["missing_full_field_images"] = [
        row["material_class"]
        for row in payload["cases"]
        if not bool(row["full_field_image_present"])
    ]
    return payload


def run_case(args, base_seed: int, class_name: str, root: Path) -> dict[str, Any]:
    global _ACTIVE_THETA_DEG
    _ACTIVE_THETA_DEG = float(args.crystal_theta_deg)
    class_name = _legacy.normalize_class_name(class_name)
    run_root = root / f"seed_{base_seed}" / "tip_only"
    case_dir = (
        run_root
        / class_name
        / f"T{int(round(args.T_K))}_th{args.crystal_theta_deg:g}"
    )
    local_summary = case_dir / "rcurve_temperature_summary.csv"
    shared_summary = run_root / "rcurve_temperature_summary.csv"
    local_before = _file_signature(local_summary)
    shared_before = _file_signature(shared_summary)

    row = _legacy_run_case(args, base_seed, class_name, root)
    reused = bool(row.get("solver_output_reused", False))
    try:
        subprocess_ok = int(row.get("subprocess_returncode", -999)) == 0
    except (TypeError, ValueError):
        subprocess_ok = False
    local_after = _file_signature(local_summary)
    shared_after = _file_signature(shared_summary)
    local_changed = local_after is not None and local_after != local_before
    shared_changed = shared_after is not None and shared_after != shared_before
    summary_source = None
    summary_fresh = False

    if reused:
        if local_after is not None and local_after[1] > 0:
            summary_source = local_summary
            summary_fresh = True
    elif subprocess_ok:
        if local_changed and local_after[1] > 0:
            summary_source = local_summary
            summary_fresh = True
        elif shared_changed and shared_after[1] > 0:
            # The v9.11 runner writes this file at --outroot. Capture it before
            # the next class subprocess can overwrite it.
            case_dir.mkdir(parents=True, exist_ok=True)
            shutil.copy2(shared_summary, local_summary)
            summary_source = local_summary
            summary_fresh = True

    summary = _legacy.read_csv_row(summary_source) if summary_source else {}
    merged = _merge_observables(row, summary)
    fp_path = case_dir / "anisotropic_calibrated_tip_first_passage_summary.json"
    try:
        fp = json.loads(fp_path.read_text()) if fp_path.exists() else {}
    except Exception:
        fp = {}
    final_um = merged.get("final_extension_um")
    try:
        target_reached = float(final_um) >= float(args.target_extension_um) * (1.0 - 1.0e-3)
    except (TypeError, ValueError):
        target_reached = False
    contract = {
        "schema": "v9.12.1_case_completion_contract",
        "material_class": class_name,
        "base_seed": int(base_seed),
        "T_K": float(args.T_K),
        "crystal_theta_deg": float(args.crystal_theta_deg),
        "requested_target_extension_um": float(args.target_extension_um),
        "final_extension_um": final_um,
        "target_extension_reached": bool(target_reached),
        "subprocess_returncode": int(merged.get("subprocess_returncode", -999)),
        "solver_status": str(merged.get("status", "unknown")),
        "control_state": str(fp.get("control_state", "unknown")),
        "summary_source_path": str(summary_source) if summary_source else None,
        "summary_fresh_for_this_invocation": bool(summary_fresh),
        "local_summary_changed_during_invocation": bool(local_changed),
        "shared_summary_changed_during_invocation": bool(shared_changed),
        "shared_summary_captured_case_locally": bool(
            summary_source == local_summary and shared_changed and not reused
        ),
        "solver_output_reused": reused,
    }
    case_dir.mkdir(parents=True, exist_ok=True)
    (case_dir / "v9_12_1_case_contract.json").write_text(
        json.dumps(contract, indent=2, default=str)
    )
    merged.update(
        {
            "v9_12_1_summary_source_path": contract["summary_source_path"],
            "v9_12_1_summary_fresh": contract["summary_fresh_for_this_invocation"],
            "v9_12_1_shared_summary_captured": contract[
                "shared_summary_captured_case_locally"
            ],
            "v9_12_1_target_extension_reached": contract[
                "target_extension_reached"
            ],
        }
    )
    pd.DataFrame([merged]).to_csv(case_dir / "v9_12_case_summary.csv", index=False)
    (case_dir / "v9_12_case_summary.json").write_text(
        json.dumps(merged, indent=2, default=str)
    )
    return merged


def main() -> None:
    original_run_case = _legacy.run_case
    original_audit = _legacy.audit_campaign
    try:
        _legacy.run_case = run_case
        _legacy.audit_campaign = audit_campaign
        _legacy.main()
    finally:
        _legacy.run_case = original_run_case
        _legacy.audit_campaign = original_audit


if __name__ == "__main__":
    main()
