#!/usr/bin/env python3
"""Normalize legacy v9.11 summary files after a certified v10.0.3 run.

This utility changes reporting only. It does not recompute mechanics, kinetics,
cohesive damage, crack extension, or any fitted/material parameter.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

REPORTING_POINT_RELEASE = "10.0.3.1"
SCHEMA = "kinetic_campaign_czm_reporting_normalization_v10_0_3_1"


def _load(path: Path) -> Any:
    if not path.exists():
        raise FileNotFoundError(path)
    return json.loads(path.read_text())


def _finite_or_none(value: Any) -> float | None:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if math.isfinite(out) else None


def _temperature(row: dict[str, Any]) -> float | None:
    return _finite_or_none(row.get("T_K", row.get("T")))


def _match_result(results: list[dict[str, Any]], legacy: dict[str, Any]) -> dict[str, Any]:
    target = _temperature(legacy)
    if target is None and len(results) == 1:
        return results[0]
    matches = [row for row in results if _temperature(row) == target]
    if len(matches) != 1:
        raise RuntimeError(
            f"could not uniquely match legacy summary temperature {target!r}; "
            f"found {len(matches)} matches"
        )
    return matches[0]


def _write_csv(path: Path, payload: dict[str, Any]) -> None:
    with path.open("w", newline="") as fp:
        writer = csv.DictWriter(fp, fieldnames=list(payload))
        writer.writeheader()
        writer.writerow(payload)


def normalize(root: str | Path) -> dict[str, Any]:
    root = Path(root)
    results_path = root / "mode_i_v10_0_3_results.json"
    model_audit_path = root / "kinetic_campaign_czm_v10_0_3_audit.json"
    legacy_json_path = root / "anisotropic_calibrated_tip_first_passage_summary.json"
    legacy_csv_path = root / "anisotropic_calibrated_tip_first_passage_summary.csv"

    results = _load(results_path)
    model_audit = _load(model_audit_path)
    legacy = _load(legacy_json_path)
    if not isinstance(results, list) or not results or not all(isinstance(x, dict) for x in results):
        raise RuntimeError("mode_i_v10_0_3_results.json must contain a nonempty list of objects")
    if not isinstance(legacy, dict):
        raise RuntimeError("legacy anisotropic summary must be one JSON object")

    result = _match_result(results, legacy)
    runtime = model_audit.get("runtime", {}) if isinstance(model_audit, dict) else {}
    result_checks = model_audit.get("result_checks", []) if isinstance(model_audit, dict) else []
    check = next(
        (
            row for row in result_checks
            if isinstance(row, dict) and _temperature(row) == _temperature(result)
        ),
        {},
    )

    B_final = _finite_or_none(result.get("B_final"))
    if B_final is None:
        B_final = _finite_or_none(check.get("B_final"))
    if B_final is None:
        raise RuntimeError("certified v10.0.3 output does not contain a finite B_final")

    before = {
        "model": legacy.get("model"),
        "B_final": legacy.get("B_final"),
        "front_state_model": legacy.get("front_state_model"),
        "front_state_model_detail": legacy.get("front_state_model_detail"),
    }

    normalized = dict(legacy)
    normalized.update({
        "legacy_wrapper_model": legacy.get("model"),
        "model": result.get("model"),
        "integration_point_release": result.get("point_release", "10.0.3"),
        "reporting_point_release": REPORTING_POINT_RELEASE,
        "front_state_model": result.get("front_state_model", "kinetic_campaign_czm"),
        "front_state_model_detail": result.get(
            "front_state_model_detail",
            "pf_v10_1_7_1_campaign_calibrated_continuous_tip_reset_safe_v1003",
        ),
        "B_final": B_final,
        "crack_extension_final_m": result.get(
            "crack_extension_final_m", check.get("crack_extension_m")
        ),
        "max_N_em": result.get("max_N_em", check.get("max_N_em")),
        "source_budget_total": result.get(
            "source_budget_total", runtime.get("source_budget_total")
        ),
        "source_population_bound": runtime.get("source_population_bound"),
        "progressive_runtime_audit": result.get("progressive_runtime_audit"),
        "full_progressive_trial_loop_active": bool(
            model_audit.get("full_progressive_trial_loop_active", False)
        ),
        "live_binding_capture_verified": bool(
            model_audit.get("live_binding_capture_verified", False)
        ),
        "reporting_normalized": True,
        "reporting_normalization_physics_changed": False,
    })

    if normalized["front_state_model"] != "kinetic_campaign_czm":
        raise RuntimeError("normalization did not produce the campaign state label")
    if not normalized["full_progressive_trial_loop_active"]:
        raise RuntimeError("refusing to normalize an uncertified progressive run")
    if not normalized["live_binding_capture_verified"]:
        raise RuntimeError("refusing to normalize a run without verified live bindings")

    legacy_json_path.write_text(json.dumps(normalized, indent=2, default=str))
    _write_csv(legacy_csv_path, normalized)

    normalized_results = []
    for row in results:
        out = dict(row)
        out["reporting_point_release"] = REPORTING_POINT_RELEASE
        out["reporting_normalized"] = True
        out["reporting_normalization_physics_changed"] = False
        normalized_results.append(out)
    normalized_results_path = root / "mode_i_v10_0_3_1_results.json"
    normalized_results_path.write_text(json.dumps(normalized_results, indent=2, default=str))

    audit = {
        "schema": SCHEMA,
        "reporting_point_release": REPORTING_POINT_RELEASE,
        "integration_point_release": result.get("point_release", "10.0.3"),
        "physics_recomputed": False,
        "mechanics_changed": False,
        "kinetics_changed": False,
        "cohesive_state_changed": False,
        "material_parameters_changed": False,
        "normalized_files": [
            legacy_json_path.name,
            legacy_csv_path.name,
            normalized_results_path.name,
        ],
        "before": before,
        "after": {
            "model": normalized.get("model"),
            "B_final": normalized.get("B_final"),
            "front_state_model": normalized.get("front_state_model"),
            "front_state_model_detail": normalized.get("front_state_model_detail"),
        },
    }
    audit_path = root / "reporting_normalization_v10_0_3_1.json"
    audit_path.write_text(json.dumps(audit, indent=2, default=str))
    return audit


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    args = parser.parse_args()
    audit = normalize(args.root)
    print("V10.0.3.1 REPORTING NORMALIZATION PASSED")
    print(json.dumps(audit["after"], indent=2, default=str))


if __name__ == "__main__":
    main()
