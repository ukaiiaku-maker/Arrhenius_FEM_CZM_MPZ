#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path

from analyze_v10_0_4_penalty_convergence import load_case


@dataclass(frozen=True)
class MatrixCase:
    case_id: str
    material: str
    temperature_K: float
    Kc_MPa_sqrt_m: float
    event_Uapp_m: float
    event_Ftop_N: float
    max_N_em: float
    source_population_bound: float
    parameter_fingerprint_sha256: str
    mode_classification: str
    certified: bool


def relative_endpoint_change(low: float, high: float) -> float:
    return abs(float(high) - float(low)) / max(0.5 * (abs(float(low)) + abs(float(high))), 1e-30)


def load_matrix_case(root: Path) -> MatrixCase:
    base = load_case(root)
    metadata = json.loads((root / "v10_0_4_case_metadata.json").read_text())
    summary = json.loads(
        (root / "anisotropic_calibrated_tip_first_passage_summary.json").read_text()
    )
    fingerprint = metadata.get("parameter_fingerprint_sha256")
    if not fingerprint:
        raise RuntimeError(f"{root}: missing parameter fingerprint")
    return MatrixCase(
        case_id=root.name,
        material=base.material,
        temperature_K=base.temperature_K,
        Kc_MPa_sqrt_m=base.Kc_MPa_sqrt_m,
        event_Uapp_m=base.event_Uapp_m,
        event_Ftop_N=base.event_Ftop_N,
        max_N_em=base.max_N_em,
        source_population_bound=base.source_population_bound,
        parameter_fingerprint_sha256=str(fingerprint),
        mode_classification=str(summary.get("mode_classification", summary.get("mode", "unknown"))),
        certified=bool(base.certified and base.quality_veto_count == 0 and base.full_rollbacks == 0),
    )


def certify(
    cases: list[MatrixCase],
    expected_materials: set[str],
    expected_temperatures: list[float],
    ceramic_endpoint_tol: float,
    weakT_endpoint_tol: float,
    dbtt_ratio_min: float,
    dbtt_delta_min: float,
    dbtt_contrast_margin: float,
) -> dict:
    failures: list[str] = []
    if not all(row.certified for row in cases):
        failures.append("one_or_more_cases_uncertified")
    grouped: dict[str, list[MatrixCase]] = {}
    for row in cases:
        grouped.setdefault(row.material, []).append(row)
    if set(grouped) != expected_materials:
        failures.append("material_set_mismatch")

    expected_T = sorted(float(x) for x in expected_temperatures)
    class_metrics = {}
    fingerprints = {}
    for material in sorted(expected_materials):
        rows = sorted(grouped.get(material, []), key=lambda row: row.temperature_K)
        temperatures = [row.temperature_K for row in rows]
        if len(rows) != len(expected_T) or any(
            not math.isclose(a, b, rel_tol=0.0, abs_tol=1e-9)
            for a, b in zip(temperatures, expected_T)
        ):
            failures.append(f"{material}_temperature_schedule_mismatch")
            continue
        fps = {row.parameter_fingerprint_sha256 for row in rows}
        if len(fps) != 1:
            failures.append(f"{material}_fingerprint_changes_with_temperature")
        else:
            fingerprints[material] = next(iter(fps))
        low, high = rows[0], rows[-1]
        endpoint_change = relative_endpoint_change(low.Kc_MPa_sqrt_m, high.Kc_MPa_sqrt_m)
        endpoint_ratio = high.Kc_MPa_sqrt_m / max(low.Kc_MPa_sqrt_m, 1e-30)
        values = [row.Kc_MPa_sqrt_m for row in rows]
        class_metrics[material] = {
            "temperatures_K": temperatures,
            "Kc_MPa_sqrt_m": values,
            "low_T_Kc_MPa_sqrt_m": low.Kc_MPa_sqrt_m,
            "high_T_Kc_MPa_sqrt_m": high.Kc_MPa_sqrt_m,
            "endpoint_relative_change": endpoint_change,
            "endpoint_ratio_high_over_low": endpoint_ratio,
            "Kc_span_MPa_sqrt_m": max(values) - min(values),
            "mode_classifications": [row.mode_classification for row in rows],
        }

    if len(set(fingerprints.values())) != len(expected_materials):
        failures.append("material_parameter_fingerprints_not_distinct")

    ceramic_change = class_metrics.get("ceramic", {}).get("endpoint_relative_change", math.inf)
    weak_change = class_metrics.get("weakT", {}).get("endpoint_relative_change", math.inf)
    dbtt_ratio = class_metrics.get("DBTT", {}).get("endpoint_ratio_high_over_low", -math.inf)
    dbtt_delta = (
        class_metrics.get("DBTT", {}).get("high_T_Kc_MPa_sqrt_m", -math.inf)
        - class_metrics.get("DBTT", {}).get("low_T_Kc_MPa_sqrt_m", math.inf)
    )
    dbtt_change = class_metrics.get("DBTT", {}).get("endpoint_relative_change", -math.inf)

    if ceramic_change > ceramic_endpoint_tol:
        failures.append("ceramic_temperature_dependence_too_large")
    if weak_change > weakT_endpoint_tol:
        failures.append("weakT_temperature_dependence_too_large")
    if dbtt_ratio < dbtt_ratio_min:
        failures.append("DBTT_high_low_ratio_too_small")
    if dbtt_delta < dbtt_delta_min:
        failures.append("DBTT_high_low_toughness_difference_too_small")
    if dbtt_change < max(ceramic_change, weak_change) + dbtt_contrast_margin:
        failures.append("DBTT_response_not_distinct_from_weak_temperature_classes")

    limits = {
        "ceramic_endpoint_relative_change_max": ceramic_endpoint_tol,
        "weakT_endpoint_relative_change_max": weakT_endpoint_tol,
        "DBTT_high_low_ratio_min": dbtt_ratio_min,
        "DBTT_high_low_delta_MPa_sqrt_m_min": dbtt_delta_min,
        "DBTT_contrast_margin_min": dbtt_contrast_margin,
    }
    certified = not failures
    return {
        "schema": "v10_0_4_parameterization_matrix_certification",
        "certified": certified,
        "integration_certified": all(row.certified for row in cases),
        "class_behavior_certified": certified,
        "expected_materials": sorted(expected_materials),
        "expected_temperatures_K": expected_T,
        "limits": limits,
        "class_metrics": class_metrics,
        "parameter_fingerprints": fingerprints,
        "failures": failures,
        "cases": [asdict(row) for row in sorted(
            cases, key=lambda row: (row.material, row.temperature_K)
        )],
        "short_growth_authorized": certified,
        "long_growth_authorized": False,
        "temperature_sweep_authorized": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("--ceramic-endpoint-tol", type=float, default=0.35)
    parser.add_argument("--weakT-endpoint-tol", type=float, default=0.35)
    parser.add_argument("--dbtt-ratio-min", type=float, default=1.5)
    parser.add_argument("--dbtt-delta-min", type=float, default=2.0)
    parser.add_argument("--dbtt-contrast-margin", type=float, default=0.15)
    args = parser.parse_args()

    plan = json.loads((args.root / "v10_0_4_parameter_matrix_plan.json").read_text())
    expected_materials = set(plan["materials"])
    expected_temperatures = [float(x) for x in plan["temperatures_K"]]
    roots = sorted(path.parent for path in args.root.glob("*/v10_0_4_case_metadata.json"))
    cases = [load_matrix_case(path) for path in roots]
    if len(cases) != int(plan["n_cases"]):
        raise RuntimeError(f"expected {plan['n_cases']} completed cases; found {len(cases)}")

    report = certify(
        cases,
        expected_materials,
        expected_temperatures,
        args.ceramic_endpoint_tol,
        args.weakT_endpoint_tol,
        args.dbtt_ratio_min,
        args.dbtt_delta_min,
        args.dbtt_contrast_margin,
    )
    (args.root / "v10_0_4_parameterization_matrix_certification.json").write_text(
        json.dumps(report, indent=2, default=str)
    )
    rows = report["cases"]
    with (args.root / "v10_0_4_parameterization_matrix_cases.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps({
        "certified": report["certified"],
        "class_metrics": report["class_metrics"],
        "failures": report["failures"],
        "short_growth_authorized": report["short_growth_authorized"],
    }, indent=2))
    if not report["certified"]:
        raise SystemExit("V10.0.4 PARAMETERIZATION MATRIX FAILED")
    print("V10.0.4 THREE-PARAMETERIZATION MATRIX CERTIFIED")


if __name__ == "__main__":
    main()
