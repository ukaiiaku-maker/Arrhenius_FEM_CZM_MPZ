#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from arrhenius_fracture.pf_equivalent_material_manifest import normalize_material_class


def parse_tokens(text: str) -> list[str]:
    return [x for x in str(text).replace(",", " ").split() if x]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--penalty-certification", type=Path, required=True)
    parser.add_argument("--materials", default="ceramic weakT DBTT")
    parser.add_argument("--temperatures", default="300 700 1100")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    penalty = json.loads(args.penalty_certification.read_text())
    if not penalty.get("certified") or not penalty.get("parameterization_matrix_authorized"):
        raise RuntimeError("cohesive-penalty convergence is not certified")
    normal = float(penalty["recommended_normal_penalty_Pa_per_m"])
    tangent = float(penalty["recommended_tangent_penalty_Pa_per_m"])

    materials = []
    for value in parse_tokens(args.materials):
        normalized = normalize_material_class(value)
        if normalized not in materials:
            materials.append(normalized)
    required = {"ceramic", "weakT", "DBTT"}
    if set(materials) != required:
        raise ValueError(f"completion matrix requires exactly {sorted(required)}; got {materials}")

    temperatures = sorted({float(x) for x in parse_tokens(args.temperatures)})
    if len(temperatures) < 3:
        raise ValueError("completion matrix requires at least three temperatures")
    if temperatures[0] > 300.0 or temperatures[-1] < 1100.0:
        raise ValueError("temperature matrix must span at least 300--1100 K")

    cases = []
    for material in materials:
        for temperature in temperatures:
            tag = f"{temperature:g}".replace(".", "p")
            cases.append({
                "case_id": f"{material}_T{tag}K",
                "material": material,
                "temperature_K": temperature,
                "normal_penalty": normal,
                "tangent_penalty": tangent,
            })

    args.out.mkdir(parents=True, exist_ok=False)
    with (args.out / "v10_0_4_parameter_matrix_plan.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(cases[0]))
        writer.writeheader()
        writer.writerows(cases)
    payload = {
        "schema": "v10_0_4_parameterization_matrix_plan",
        "penalty_certification": str(args.penalty_certification.resolve()),
        "normal_penalty_Pa_per_m": normal,
        "tangent_penalty_Pa_per_m": tangent,
        "materials": materials,
        "temperatures_K": temperatures,
        "n_cases": len(cases),
        "cases": cases,
    }
    (args.out / "v10_0_4_parameter_matrix_plan.json").write_text(json.dumps(payload, indent=2))
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
