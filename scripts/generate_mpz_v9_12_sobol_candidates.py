#!/usr/bin/env python3
"""Generate v9.12 candidate rows from a base row and explicit search bounds."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

from scipy.stats import qmc


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-registry", required=True)
    parser.add_argument("--base-candidate-id", required=True)
    parser.add_argument("--bounds-json", required=True)
    parser.add_argument("--n", type=int, required=True)
    parser.add_argument("--seed", type=int, default=912)
    parser.add_argument("--out", required=True)
    parser.add_argument("--prefix", default="v912_emergent_gnd")
    return parser.parse_args()


def load_base(path: str, candidate_id: str) -> tuple[list[str], dict[str, str]]:
    with Path(path).open(newline="") as fp:
        reader = csv.DictReader(fp)
        rows = list(reader)
        fields = list(reader.fieldnames or [])
    matches = [row for row in rows if row.get("candidate_id") == candidate_id]
    if len(matches) != 1:
        raise RuntimeError(
            f"expected one base row for {candidate_id!r}; found {len(matches)}"
        )
    return fields, dict(matches[0])


def transform(u: float, spec: dict[str, object]) -> float:
    low = float(spec["low"])
    high = float(spec["high"])
    scale = str(spec.get("scale", "linear")).lower()
    if high < low:
        raise ValueError(f"invalid bounds: low={low}, high={high}")
    if scale == "linear":
        return low + float(u) * (high - low)
    if scale == "log10":
        if low <= 0.0 or high <= 0.0:
            raise ValueError("log10 bounds must be positive")
        return 10.0 ** (
            math.log10(low) + float(u) * (math.log10(high) - math.log10(low))
        )
    raise ValueError(f"unknown bound scale: {scale}")


def main() -> int:
    args = parse_args()
    fields, base = load_base(args.base_registry, args.base_candidate_id)
    bounds_payload = json.loads(Path(args.bounds_json).read_text())
    bounds = bounds_payload.get("search_bounds", bounds_payload)
    names = list(bounds)
    if not names:
        raise RuntimeError("search bounds are empty")
    for name in names:
        if name not in fields:
            fields.append(name)
    for extra in ("campaign_parent_id", "campaign_generator", "campaign_seed"):
        if extra not in fields:
            fields.append(extra)

    m = int(math.ceil(math.log2(max(args.n, 1))))
    sampler = qmc.Sobol(d=len(names), scramble=True, seed=args.seed)
    points = sampler.random_base2(m=m)[: args.n]

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fp:
        writer = csv.DictWriter(fp, fieldnames=fields)
        writer.writeheader()
        for index, point in enumerate(points):
            row = dict(base)
            row["candidate_id"] = f"{args.prefix}_{index:06d}"
            row["campaign_parent_id"] = args.base_candidate_id
            row["campaign_generator"] = "scipy.stats.qmc.Sobol"
            row["campaign_seed"] = str(args.seed)
            for name, u in zip(names, point):
                row[name] = f"{transform(float(u), bounds[name]):.17g}"
            writer.writerow(row)
    print(
        f"GENERATED candidates={args.n} dimensions={len(names)} out={out}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
