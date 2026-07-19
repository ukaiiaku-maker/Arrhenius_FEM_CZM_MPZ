#!/usr/bin/env python3
"""v10.0.5.7 repaired full-field material R-curve campaign.

The mechanics remain the validated v9.12 FEM/CZM implementation. This wrapper
repairs the v9.12 summary-path bug, persists the authoritative v9.11 solver row
inside each case directory, and replaces the vacuous publication audit with the
strict v10.0.5.7 gate.
"""
from __future__ import annotations

import json
from pathlib import Path
import shutil
from typing import Any

import pandas as pd

import run_mpz_v9_12_tip_only_material_rcurve as _base
from arrhenius_fracture.material_rcurve_audit_v10057 import audit_campaign

POINT_RELEASE = "10.0.5.7"
_ORIGINAL_RUN_CASE = _base.run_case


def _authoritative_temperature_summary(
    run_root: Path,
    class_name: str,
    T_K: float,
) -> tuple[dict[str, Any], Path]:
    """Read the actual v9.11 summary location and select the requested case."""
    source = run_root / "rcurve_temperature_summary.csv"
    if not source.exists() or source.stat().st_size == 0:
        return {}, source
    try:
        frame = pd.read_csv(source)
    except (pd.errors.EmptyDataError, OSError, ValueError):
        return {}, source
    if frame.empty:
        return {}, source

    selected = frame
    if "class" in selected.columns:
        selected = selected[selected["class"].astype(str) == str(class_name)]
    if "T_K" in selected.columns:
        T = pd.to_numeric(selected["T_K"], errors="coerce")
        selected = selected[(T - float(T_K)).abs() <= 0.5]
    if selected.empty:
        return {}, source
    return selected.iloc[0].to_dict(), source


def _persist_repaired_case_summary(
    *,
    legacy_row: dict[str, Any],
    solver_row: dict[str, Any],
    source: Path,
    case_dir: Path,
) -> dict[str, Any]:
    """Merge authoritative solver lifecycle fields with v9.12 metadata."""
    repaired = dict(legacy_row)
    repaired.update(solver_row)

    # v9.12 campaign metadata must not be overwritten by generic v9.11 columns.
    for key in (
        "class",
        "base_seed",
        "effective_stochastic_seed",
        "rng_coupling",
        "stochastic_emission",
        "propagation_control",
        "subprocess_returncode",
        "solver_output_reused",
        "field_snapshot_image",
        "field_snapshot_image_present",
        "case_dir",
        "log",
    ):
        if key in legacy_row:
            repaired[key] = legacy_row[key]

    copied = False
    copied_path = case_dir / "rcurve_temperature_summary_v9_11.csv"
    if source.exists() and source.stat().st_size > 0:
        shutil.copy2(source, copied_path)
        copied = True

    repaired.update(
        {
            "point_release": POINT_RELEASE,
            "solver_summary_source": str(source),
            "solver_summary_source_present": bool(source.exists()),
            "solver_summary_copied_to_case": copied,
            "solver_summary_copy": str(copied_path) if copied else None,
        }
    )
    if not solver_row:
        repaired.update(
            {
                "status": "missing_solver_summary",
                "target_completed": False,
                "summary_path_error": (
                    "v9.11 did not produce the authoritative root-level "
                    "rcurve_temperature_summary.csv"
                ),
            }
        )

    pd.DataFrame([repaired]).to_csv(case_dir / "v9_12_case_summary.csv", index=False)
    (case_dir / "v9_12_case_summary.json").write_text(
        json.dumps(repaired, indent=2, default=str)
    )
    return repaired


def run_case(args, base_seed: int, class_name: str, root: Path) -> dict[str, Any]:
    """Run the validated v9.12 case and repair its bookkeeping immediately."""
    legacy_row = _ORIGINAL_RUN_CASE(args, base_seed, class_name, root)
    normalized = _base.normalize_class_name(class_name)
    run_root = root / f"seed_{base_seed}" / "tip_only"
    case_dir = run_root / normalized / f"T{int(round(args.T_K))}_th{args.crystal_theta_deg:g}"
    solver_row, source = _authoritative_temperature_summary(
        run_root, normalized, args.T_K
    )
    repaired = _persist_repaired_case_summary(
        legacy_row=legacy_row,
        solver_row=solver_row,
        source=source,
        case_dir=case_dir,
    )
    print(
        f"REPAIRED {normalized:7s} status={repaired.get('status')} "
        f"target_completed={repaired.get('target_completed')} "
        f"control_state={repaired.get('control_state')}"
    )
    return repaired


def main() -> None:
    saved_run_case = _base.run_case
    saved_audit = _base.audit_campaign
    _base.run_case = run_case
    _base.audit_campaign = audit_campaign
    try:
        _base.main()
    finally:
        _base.run_case = saved_run_case
        _base.audit_campaign = saved_audit


if __name__ == "__main__":
    main()


__all__ = [
    "POINT_RELEASE",
    "_authoritative_temperature_summary",
    "_persist_repaired_case_summary",
    "run_case",
    "main",
]
