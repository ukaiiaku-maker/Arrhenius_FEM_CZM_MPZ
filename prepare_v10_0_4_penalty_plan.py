#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def parse_values(text: str) -> list[float]:
    values = [float(x) for x in str(text).replace(",", " ").split() if x]
    if not values:
        raise ValueError("at least one penalty value is required")
    if any(x <= 0.0 for x in values):
        raise ValueError("cohesive penalties must be positive")
    return values


def label(value: float) -> str:
    text = f"{value:.6g}".replace("+", "")
    return text.replace(".", "p")


def build_plan(normal_values: list[float], tangent_values: list[float], reference: float):
    cases: dict[tuple[float, float], dict] = {}
    for value in normal_values:
        key = (float(value), float(reference))
        cases[key] = {
            "case_id": f"normal_{label(value)}__tangent_{label(reference)}",
            "axis": "normal",
            "normal_penalty": float(value),
            "tangent_penalty": float(reference),
        }
    for value in tangent_values:
        key = (float(reference), float(value))
        row = {
            "case_id": f"normal_{label(reference)}__tangent_{label(value)}",
            "axis": "tangent",
            "normal_penalty": float(reference),
            "tangent_penalty": float(value),
        }
        if key in cases:
            cases[key]["axis"] = "baseline"
        else:
            cases[key] = row
    ordered = sorted(cases.values(), key=lambda row: (
        row["normal_penalty"], row["tangent_penalty"], row["case_id"]
    ))
    return ordered


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--normal-penalties", default="5e17 1e18 2e18")
    parser.add_argument("--tangent-penalties", default="5e17 1e18 2e18")
    parser.add_argument("--reference-penalty", type=float, default=1e18)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    normal = parse_values(args.normal_penalties)
    tangent = parse_values(args.tangent_penalties)
    reference = float(args.reference_penalty)
    if reference <= 0.0:
        raise ValueError("reference penalty must be positive")
    if not any(abs(x - reference) <= 1e-12 * reference for x in normal):
        raise ValueError("normal penalties must include the reference penalty")
    if not any(abs(x - reference) <= 1e-12 * reference for x in tangent):
        raise ValueError("tangent penalties must include the reference penalty")

    cases = build_plan(normal, tangent, reference)
    args.out.mkdir(parents=True, exist_ok=False)
    csv_path = args.out / "v10_0_4_penalty_plan.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=[
            "case_id", "axis", "normal_penalty", "tangent_penalty"
        ])
        writer.writeheader()
        writer.writerows(cases)
    payload = {
        "schema": "v10_0_4_penalty_plan",
        "reference_penalty_Pa_per_m": reference,
        "normal_penalties_Pa_per_m": normal,
        "tangent_penalties_Pa_per_m": tangent,
        "n_unique_cases": len(cases),
        "cases": cases,
    }
    (args.out / "v10_0_4_penalty_plan.json").write_text(json.dumps(payload, indent=2))
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
