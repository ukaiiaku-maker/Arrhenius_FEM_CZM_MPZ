#!/usr/bin/env python3
"""Audited launcher for the v10.0.5.6 KJ audit and bracket workflow."""
from __future__ import annotations

import json
from pathlib import Path

import run_v10_0_5_6_kj_audit_bracket as _base


AUDITED_CAMPAIGN = (
    Path(__file__).resolve().parent
    / "run_v10_0_5_6_stochastic_delta_sigma_audited.py"
)
KJ_LEFM_RATIO_MIN = 0.50
KJ_LEFM_RATIO_MAX = 1.50


def _read_csv_with_booleans(path):
    rows = _ORIGINAL_READ_CSV(path)
    for row in rows:
        for key in (
            "first_passage_observed",
            "right_censored",
            "reached_cycle_horizon",
            "reached_target_extension",
        ):
            value = row.get(key)
            if isinstance(value, str):
                token = value.strip().lower()
                if token in {"true", "1", "yes"}:
                    row[key] = True
                elif token in {"false", "0", "no", ""}:
                    row[key] = False
    return rows


def _classify_without_numerical_censoring(rows):
    materialized = [dict(row) for row in rows]
    censored = [row for row in materialized if bool(row.get("right_censored"))]
    if censored:
        stresses = [float(row["delta_sigma_requested_MPa"]) for row in censored]
        raise RuntimeError(
            "first-passage bracket contains numerically censored cases at "
            f"Delta-sigma={stresses} MPa; increase MAX_BLOCKS or revise block controls"
        )
    return _ORIGINAL_CLASSIFY(materialized)


def _run_audit_fail_closed(args):
    payload = _ORIGINAL_RUN_AUDIT(args)
    out = Path(args.out).resolve()
    selected_path = out / _base.SELECTED_JSON
    selected = json.loads(selected_path.read_text())
    if selected.get("status") == "plateau_selected":
        ratio_min = float(selected.get("plateau_KJ_over_K_LEFM_min", float("nan")))
        ratio_max = float(selected.get("plateau_KJ_over_K_LEFM_max", float("nan")))
        accepted = (
            ratio_min >= KJ_LEFM_RATIO_MIN
            and ratio_max <= KJ_LEFM_RATIO_MAX
        )
        selected["KJ_LEFM_ratio_acceptance"] = [
            KJ_LEFM_RATIO_MIN,
            KJ_LEFM_RATIO_MAX,
        ]
        selected["KJ_LEFM_ratio_audit_passed"] = bool(accepted)
        if not accepted:
            selected["status"] = "plateau_selected_but_KJ_LEFM_mismatch"
        selected_path.write_text(json.dumps(selected, indent=2, default=str))
        audit_path = out / _base.AUDIT_JSON
        if audit_path.exists():
            audit = json.loads(audit_path.read_text())
            audit["selection"] = selected
            audit["KJ_LEFM_ratio_acceptance"] = [
                KJ_LEFM_RATIO_MIN,
                KJ_LEFM_RATIO_MAX,
            ]
            audit_path.write_text(json.dumps(audit, indent=2, default=str))
        payload["selection"] = selected
    return payload


_ORIGINAL_READ_CSV = _base._read_csv
_ORIGINAL_RUN_AUDIT = _base.run_audit
_ORIGINAL_CLASSIFY = _base.classify_first_passage_rows


def main(argv=None) -> int:
    saved_campaign = _base.CAMPAIGN
    saved_reader = _base._read_csv
    saved_audit = _base.run_audit
    saved_classify = _base.classify_first_passage_rows
    _base.CAMPAIGN = AUDITED_CAMPAIGN
    _base._read_csv = _read_csv_with_booleans
    _base.run_audit = _run_audit_fail_closed
    _base.classify_first_passage_rows = _classify_without_numerical_censoring
    try:
        return int(_base.main(argv) or 0)
    finally:
        _base.CAMPAIGN = saved_campaign
        _base._read_csv = saved_reader
        _base.run_audit = saved_audit
        _base.classify_first_passage_rows = saved_classify


if __name__ == "__main__":
    raise SystemExit(main())
