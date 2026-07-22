#!/usr/bin/env python3
"""Generate local Sobol families around selected v9.12 parent candidates."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

from scipy.stats import qmc


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--base-registry", required=True)
    p.add_argument("--policy-json", required=True)
    p.add_argument("--per-seed", type=int, default=512)
    p.add_argument("--seed", type=int, default=19120)
    p.add_argument("--prefix", default="v912_targeted_local")
    p.add_argument("--out", required=True)
    return p.parse_args()


def read_rows(path: str | Path) -> tuple[list[str], dict[str, dict[str, str]]]:
    with Path(path).open(newline="") as fp:
        reader = csv.DictReader(fp)
        rows = list(reader)
        fields = list(reader.fieldnames or [])
    by_id = {str(row["candidate_id"]): dict(row) for row in rows}
    if len(by_id) != len(rows):
        raise RuntimeError("candidate_id must be unique in base registry")
    return fields, by_id


def local_interval(center: float, spec: dict[str, Any]) -> tuple[float, float, str]:
    mode = str(spec.get("mode", "linear_delta")).lower()
    global_low = float(spec["low"])
    global_high = float(spec["high"])
    if global_high <= global_low:
        raise ValueError(f"invalid global interval {global_low}, {global_high}")

    if mode == "linear_delta":
        half = float(spec["half_width"])
        low = max(global_low, center - half)
        high = min(global_high, center + half)
        scale = "linear"
    elif mode == "log10_delta":
        if center <= 0.0 or global_low <= 0.0 or global_high <= 0.0:
            raise ValueError("log10 local intervals require positive values")
        half = float(spec["half_width_decades"])
        c = math.log10(center)
        low = max(global_low, 10.0 ** (c - half))
        high = min(global_high, 10.0 ** (c + half))
        scale = "log10"
    else:
        raise ValueError(f"unknown local interval mode: {mode}")

    if high < low:
        raise RuntimeError(
            f"empty local interval around {center}: low={low}, high={high}"
        )
    return low, high, scale


def transform(u: float, low: float, high: float, scale: str) -> float:
    if high == low:
        return low
    if scale == "linear":
        return low + float(u) * (high - low)
    if scale == "log10":
        return 10.0 ** (
            math.log10(low) + float(u) * (math.log10(high) - math.log10(low))
        )
    raise ValueError(scale)


def main() -> int:
    a = parse_args()
    if a.per_seed < 1:
        raise ValueError("--per-seed must be positive")

    fields, registry = read_rows(a.base_registry)
    policy = json.loads(Path(a.policy_json).read_text())
    dimensions = policy.get("search_dimensions", {})
    seed_families = policy.get("seed_families", {})
    if not dimensions or not seed_families:
        raise RuntimeError("policy must define search_dimensions and seed_families")

    dimension_names = list(dimensions)
    seed_items: list[tuple[str, str]] = []
    for family, seed_ids in seed_families.items():
        for candidate_id in seed_ids:
            seed_items.append((str(family), str(candidate_id)))

    missing = [candidate_id for _, candidate_id in seed_items if candidate_id not in registry]
    if missing:
        raise RuntimeError(f"missing parent candidates: {missing}")

    extras = [
        "campaign_parent_id",
        "campaign_parent_family",
        "campaign_generator",
        "campaign_seed",
        "campaign_local_policy",
    ]
    for name in dimension_names + extras:
        if name not in fields:
            fields.append(name)

    output = Path(a.out)
    output.parent.mkdir(parents=True, exist_ok=True)
    m = int(math.ceil(math.log2(a.per_seed)))
    total = 0

    with output.open("w", newline="") as fp:
        writer = csv.DictWriter(fp, fieldnames=fields)
        writer.writeheader()

        for seed_index, (family, parent_id) in enumerate(seed_items):
            parent = registry[parent_id]
            sampler = qmc.Sobol(
                d=len(dimension_names),
                scramble=True,
                seed=int(a.seed) + seed_index,
            )
            points = sampler.random_base2(m=m)[: a.per_seed]
            intervals: dict[str, tuple[float, float, str]] = {}
            for name in dimension_names:
                if name not in parent or parent[name] in (None, ""):
                    raise KeyError(f"parent {parent_id} lacks {name}")
                intervals[name] = local_interval(float(parent[name]), dimensions[name])

            short_parent = parent_id.rsplit("_", 1)[-1]
            for local_index, point in enumerate(points):
                row = dict(parent)
                row["candidate_id"] = (
                    f"{a.prefix}_{family}_{short_parent}_{local_index:04d}"
                )
                row["campaign_parent_id"] = parent_id
                row["campaign_parent_family"] = family
                row["campaign_generator"] = "local_sobol_v1"
                row["campaign_seed"] = str(int(a.seed) + seed_index)
                row["campaign_local_policy"] = str(Path(a.policy_json))
                for name, u in zip(dimension_names, point):
                    low, high, scale = intervals[name]
                    row[name] = f"{transform(float(u), low, high, scale):.17g}"
                writer.writerow(row)
                total += 1

    print(
        "TARGETED_LOCAL_POOL "
        f"parents={len(seed_items)} per_seed={a.per_seed} "
        f"candidates={total} dimensions={len(dimension_names)} out={output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
