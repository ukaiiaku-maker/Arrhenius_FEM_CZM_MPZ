#!/usr/bin/env python3
"""Audited launcher for the v10.0.5.6 KJ audit and bracket workflow."""
from __future__ import annotations

import json
import math
import os
from pathlib import Path

import run_v10_0_5_6_kj_audit_bracket as _base


AUDITED_CAMPAIGN = (
    Path(__file__).resolve().parent
    / "run_v10_0_5_6_stochastic_delta_sigma_audited.py"
)
KJ_LEFM_RATIO_MIN = 0.50
KJ_LEFM_RATIO_MAX = 1.50
STRESS_RELATIVE_TOLERANCE = 1.0e-6
KJ_LINEARITY_RELATIVE_TOLERANCE = 0.02
LINEARITY_JSON = "remote_stress_KJ_linearity_v10_0_5_6.json"


_BOOLEAN_COLUMNS = (
    "first_passage_observed",
    "right_censored",
    "reached_cycle_horizon",
    "reached_target_extension",
)


def _normalize_boolean_value(value, *, column: str, path: Path):
    """Normalize audited CSV booleans after the generic reader converts numbers.

    The base reader converts the text tokens ``0`` and ``1`` to floating-point
    values before this wrapper sees them.  Accept both textual and numeric boolean
    encodings, preserve missing values, and fail closed on any other token.
    """
    if value is None or isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        numeric = float(value)
        if numeric == 0.0:
            return False
        if numeric == 1.0:
            return True
        if math.isnan(numeric):
            return value
        raise ValueError(
            f"invalid numeric boolean {value!r} for column {column!r} in {path}"
        )
    if isinstance(value, str):
        token = value.strip().lower()
        if token in {"true", "1", "yes"}:
            return True
        if token in {"false", "0", "no", ""}:
            return False
        raise ValueError(
            f"invalid textual boolean {value!r} for column {column!r} in {path}"
        )
    raise TypeError(
        f"unsupported boolean value {value!r} for column {column!r} in {path}"
    )


def _read_csv_with_booleans(path):
    rows = _ORIGINAL_READ_CSV(path)
    csv_path = Path(path)
    for row in rows:
        for key in _BOOLEAN_COLUMNS:
            if key in row:
                row[key] = _normalize_boolean_value(
                    row.get(key), column=key, path=csv_path
                )
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


def _single_step_options():
    return [
        "--cycles-max",
        "1e-9",
        "--block-cycles",
        "1e-9",
        "--max-block-cycles",
        "1e-9",
        "--max-blocks",
        "1",
        "--target-extension-um",
        "1e-12",
        "--target-dN-store",
        "0.05",
        "--target-dN-emit",
        "inf",
        "--target-dN-mobile",
        "inf",
        "--target-dN-escape",
        "0.25",
    ]


def _linearity_verification(args, selected):
    out = Path(args.out).resolve()
    first = dict(selected["selected_row"])
    contour_um = 1.0e6 * float(selected["selected_outer_radius_m"])
    second_delta = 2.0 * float(args.audit_delta_sigma_MPa)
    second_root = out / "linearity_verification"
    command = _base._common_campaign_command(
        args,
        out=second_root,
        delta_sigma_MPa=second_delta,
        contour_um=contour_um,
    )
    command.extend(_single_step_options())
    env = os.environ.copy()
    env.update(
        {
            "ARRHENIUS_EVENT_STATISTICS": "deterministic",
            "ARRHENIUS_STOCHASTIC_EMISSION": "0",
            "ARRHENIUS_STOCHASTIC_BLOCKS": "0",
            "ARRHENIUS_VHCF_FEM_CACHE": "0",
        }
    )
    result_path = second_root / "remote_stress_KJ_audit_v10_0_5_6.csv"
    if not result_path.exists():
        _base._run(command, env=env, log=second_root / "launcher.log")
    rows = _base._read_csv(result_path)
    if len(rows) != 1:
        raise RuntimeError(f"expected one linearity row in {result_path}")
    second = dict(rows[0])

    sigma1 = float(first["sigma_gross_MPa"])
    sigma2 = float(second["sigma_gross_MPa"])
    K1 = float(first["KJ_MPa_sqrt_m"])
    K2 = float(second["KJ_MPa_sqrt_m"])
    requested1 = float(args.audit_delta_sigma_MPa) / (1.0 - float(args.R))
    requested2 = second_delta / (1.0 - float(args.R))
    sigma_ratio = sigma2 / sigma1
    K_ratio = K2 / K1
    stress_error1 = abs(sigma1 / requested1 - 1.0)
    stress_error2 = abs(sigma2 / requested2 - 1.0)
    linearity_error = abs(K_ratio / sigma_ratio - 1.0)
    passed = (
        stress_error1 <= STRESS_RELATIVE_TOLERANCE
        and stress_error2 <= STRESS_RELATIVE_TOLERANCE
        and linearity_error <= KJ_LINEARITY_RELATIVE_TOLERANCE
    )
    audit = {
        "schema": "remote_stress_KJ_linearity_v10_0_5_6",
        "first_delta_sigma_MPa": float(args.audit_delta_sigma_MPa),
        "second_delta_sigma_MPa": second_delta,
        "first_sigma_max_actual_MPa": sigma1,
        "second_sigma_max_actual_MPa": sigma2,
        "first_sigma_max_requested_MPa": requested1,
        "second_sigma_max_requested_MPa": requested2,
        "first_KJ_MPa_sqrt_m": K1,
        "second_KJ_MPa_sqrt_m": K2,
        "sigma_ratio": sigma_ratio,
        "KJ_ratio": K_ratio,
        "first_stress_relative_error": stress_error1,
        "second_stress_relative_error": stress_error2,
        "KJ_linearity_relative_error": linearity_error,
        "stress_relative_tolerance": STRESS_RELATIVE_TOLERANCE,
        "KJ_linearity_relative_tolerance": KJ_LINEARITY_RELATIVE_TOLERANCE,
        "passed": bool(passed),
        "constitutive_physics_changed": False,
    }
    (out / LINEARITY_JSON).write_text(json.dumps(audit, indent=2))
    return audit


def _run_audit_fail_closed(args):
    payload = _ORIGINAL_RUN_AUDIT(args)
    out = Path(args.out).resolve()
    selected_path = out / _base.SELECTED_JSON
    selected = json.loads(selected_path.read_text())
    if selected.get("status") == "plateau_selected":
        ratio_min = float(selected.get("plateau_KJ_over_K_LEFM_min", float("nan")))
        ratio_max = float(selected.get("plateau_KJ_over_K_LEFM_max", float("nan")))
        ratio_accepted = (
            math.isfinite(ratio_min)
            and math.isfinite(ratio_max)
            and ratio_min >= KJ_LEFM_RATIO_MIN
            and ratio_max <= KJ_LEFM_RATIO_MAX
        )
        selected["KJ_LEFM_ratio_acceptance"] = [
            KJ_LEFM_RATIO_MIN,
            KJ_LEFM_RATIO_MAX,
        ]
        selected["KJ_LEFM_ratio_audit_passed"] = bool(ratio_accepted)
        if ratio_accepted:
            linearity = _linearity_verification(args, selected)
            selected["remote_stress_KJ_linearity"] = linearity
            if not bool(linearity["passed"]):
                selected["status"] = "plateau_selected_but_linearity_failed"
        else:
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
