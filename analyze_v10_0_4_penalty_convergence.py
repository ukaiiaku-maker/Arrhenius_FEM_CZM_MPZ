#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
from typing import Iterable


@dataclass(frozen=True)
class PenaltyCase:
    case_id: str
    normal_penalty: float
    tangent_penalty: float
    material: str
    temperature_K: float
    Kc_MPa_sqrt_m: float
    event_Uapp_m: float
    event_Ftop_N: float
    event_sigma_tip_Pa: float
    max_N_em: float
    source_population_bound: float
    certified: bool
    quality_veto_count: int
    full_rollbacks: int
    committed_events: int


def finite(value, name: str) -> float:
    out = float(value)
    if not math.isfinite(out):
        raise RuntimeError(f"non-finite {name}: {value!r}")
    return out


def rel_span(values: Iterable[float]) -> float:
    vals = [finite(x, "convergence value") for x in values]
    if len(vals) < 2:
        return 0.0
    center = sorted(vals)[len(vals) // 2]
    return (max(vals) - min(vals)) / max(abs(center), 1.0e-30)


def close(a: float, b: float) -> bool:
    return math.isclose(float(a), float(b), rel_tol=1.0e-12, abs_tol=0.0)


def first_existing(root: Path, names: list[str]) -> Path:
    for name in names:
        path = root / name
        if path.exists():
            return path
    raise FileNotFoundError(f"none of {names!r} exists in {root}")


def result_Kc(summary: dict) -> float:
    for key in (
        "KJ_reference_first_MPa_sqrt_m",
        "Kcleave_calibrated_first_MPa_sqrt_m",
        "Kc_first_existing_MPa_sqrt_m",
        "Kc_first_MPa_sqrt_m",
    ):
        value = summary.get(key)
        if value is not None:
            return finite(value, key)
    raise RuntimeError("normalized summary contains no first-event toughness")


def load_case(root: Path) -> PenaltyCase:
    metadata = json.loads((root / "v10_0_4_case_metadata.json").read_text())
    certification = json.loads(
        (root / "v10_0_3_progressive_integration_certification.json").read_text()
    )
    summary = json.loads(
        (root / "anisotropic_calibrated_tip_first_passage_summary.json").read_text()
    )
    runtime = json.loads(
        (root / "kinetic_campaign_czm_progressive_2d_v10_0_3.json").read_text()
    )
    quality = json.loads((root / "explicit_quality_wrapper_chain_v91856.json").read_text())
    steps_path = first_existing(root, sorted(p.name for p in root.glob("steps_*K.csv")))
    with steps_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    fired = [row for row in rows if float(row.get("n_fire", 0.0)) > 0.0]
    if len(fired) != 1:
        raise RuntimeError(f"{root}: expected one fired accepted row; found {len(fired)}")
    event = fired[0]
    return PenaltyCase(
        case_id=root.name,
        normal_penalty=finite(metadata["czm_penalty_normal_Pa_per_m"], "normal penalty"),
        tangent_penalty=finite(metadata["czm_penalty_tangent_Pa_per_m"], "tangent penalty"),
        material=str(metadata["material"]),
        temperature_K=finite(metadata["temperature_K"], "temperature"),
        Kc_MPa_sqrt_m=result_Kc(summary),
        event_Uapp_m=finite(event["Uapp_m"], "event Uapp"),
        event_Ftop_N=finite(event["Ftop_N"], "event Ftop"),
        event_sigma_tip_Pa=finite(event["sigma_tip_Pa"], "event sigma tip"),
        max_N_em=finite(certification["max_N_em"], "max N_em"),
        source_population_bound=finite(
            certification["source_population_bound"], "source population bound"
        ),
        certified=bool(metadata.get("certified") and certification.get("certified")),
        quality_veto_count=len(quality.get("quality_vetoes") or []),
        full_rollbacks=int(runtime.get("full_rollbacks", -1)),
        committed_events=int(runtime.get("committed_events", -1)),
    )


def certify(
    cases: list[PenaltyCase],
    reference: float,
    normal_Kc_tol: float,
    normal_U_tol: float,
    tangent_Kc_tol: float,
    tangent_U_tol: float,
) -> dict:
    if len(cases) < 5:
        raise RuntimeError(f"expected at least five unique penalty cases; found {len(cases)}")
    if not all(row.certified for row in cases):
        raise RuntimeError("one or more penalty cases did not pass the v10.0.3 certification")
    if any(row.quality_veto_count for row in cases):
        raise RuntimeError("one or more penalty cases recorded a geometry quality veto")
    if any(row.full_rollbacks != 0 for row in cases):
        raise RuntimeError("one or more penalty cases required a full topology rollback")
    if any(row.committed_events != 1 for row in cases):
        raise RuntimeError("one or more penalty cases did not commit exactly one event")
    if any(row.max_N_em > row.source_population_bound + 1e-8 for row in cases):
        raise RuntimeError("one or more penalty cases violated the source-population bound")

    materials = {row.material for row in cases}
    temperatures = {row.temperature_K for row in cases}
    if len(materials) != 1 or len(temperatures) != 1:
        raise RuntimeError("penalty convergence cases must use one material and temperature")

    normal_axis = [row for row in cases if close(row.tangent_penalty, reference)]
    tangent_axis = [row for row in cases if close(row.normal_penalty, reference)]
    if len(normal_axis) < 3:
        raise RuntimeError("normal-penalty axis must contain at least three points")
    if len(tangent_axis) < 3:
        raise RuntimeError("tangential-penalty axis must contain at least three points")
    baseline = [
        row for row in cases
        if close(row.normal_penalty, reference) and close(row.tangent_penalty, reference)
    ]
    if len(baseline) != 1:
        raise RuntimeError(f"expected one reference-penalty case; found {len(baseline)}")

    metrics = {
        "normal_Kc_relative_span": rel_span(row.Kc_MPa_sqrt_m for row in normal_axis),
        "normal_event_displacement_relative_span": rel_span(row.event_Uapp_m for row in normal_axis),
        "normal_event_force_relative_span": rel_span(row.event_Ftop_N for row in normal_axis),
        "tangent_Kc_relative_span": rel_span(row.Kc_MPa_sqrt_m for row in tangent_axis),
        "tangent_event_displacement_relative_span": rel_span(row.event_Uapp_m for row in tangent_axis),
        "tangent_event_force_relative_span": rel_span(row.event_Ftop_N for row in tangent_axis),
    }
    limits = {
        "normal_Kc_relative_span": normal_Kc_tol,
        "normal_event_displacement_relative_span": normal_U_tol,
        "tangent_Kc_relative_span": tangent_Kc_tol,
        "tangent_event_displacement_relative_span": tangent_U_tol,
    }
    failures = [key for key, limit in limits.items() if metrics[key] > limit]
    certified = not failures
    return {
        "schema": "v10_0_4_cohesive_penalty_convergence",
        "certified": certified,
        "material": next(iter(materials)),
        "temperature_K": next(iter(temperatures)),
        "reference_normal_penalty_Pa_per_m": reference,
        "reference_tangent_penalty_Pa_per_m": reference,
        "recommended_normal_penalty_Pa_per_m": reference if certified else None,
        "recommended_tangent_penalty_Pa_per_m": reference if certified else None,
        "metrics": metrics,
        "limits": limits,
        "failures": failures,
        "cases": [asdict(row) for row in sorted(
            cases, key=lambda x: (x.normal_penalty, x.tangent_penalty)
        )],
        "long_growth_authorized": False,
        "parameterization_matrix_authorized": certified,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("--reference-penalty", type=float, default=1e18)
    parser.add_argument("--normal-Kc-tol", type=float, default=0.01)
    parser.add_argument("--normal-U-tol", type=float, default=0.01)
    parser.add_argument("--tangent-Kc-tol", type=float, default=0.0025)
    parser.add_argument("--tangent-U-tol", type=float, default=0.005)
    args = parser.parse_args()

    roots = sorted(
        path.parent for path in args.root.glob("*/v10_0_4_case_metadata.json")
    )
    cases = [load_case(path) for path in roots]
    report = certify(
        cases,
        args.reference_penalty,
        args.normal_Kc_tol,
        args.normal_U_tol,
        args.tangent_Kc_tol,
        args.tangent_U_tol,
    )
    json_path = args.root / "v10_0_4_penalty_convergence_certification.json"
    json_path.write_text(json.dumps(report, indent=2, default=str))
    csv_path = args.root / "v10_0_4_penalty_convergence_cases.csv"
    with csv_path.open("w", newline="") as handle:
        rows = report["cases"]
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps({
        "certified": report["certified"],
        "metrics": report["metrics"],
        "failures": report["failures"],
        "recommended_normal_penalty_Pa_per_m": report["recommended_normal_penalty_Pa_per_m"],
        "recommended_tangent_penalty_Pa_per_m": report["recommended_tangent_penalty_Pa_per_m"],
    }, indent=2))
    if not report["certified"]:
        raise SystemExit("V10.0.4 COHESIVE PENALTY CONVERGENCE FAILED")
    print("V10.0.4 COHESIVE PENALTY CONVERGENCE CERTIFIED")


if __name__ == "__main__":
    main()
