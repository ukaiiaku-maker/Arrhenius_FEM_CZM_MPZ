#!/usr/bin/env python3
"""Extract the promoted peak candidates and validate a long loading map."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from arrhenius_fracture.dbtt_long_alignment_v913 import (
    loading_map_coverage_um,
    truthy,
)
from arrhenius_fracture.emergent_gnd_contract_v913 import (
    ACTIVE_CANDIDATE_PARAMETER_FIELDS,
    effective_candidate_parameters,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-ranking", type=Path, required=True)
    parser.add_argument("--loading-map", type=Path, required=True)
    parser.add_argument("--out-registry", type=Path, required=True)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--selection-column", default="y__peak_like_1d")
    parser.add_argument("--selected-count", type=int, default=25)
    parser.add_argument("--target-extension-um", type=float, default=100.0)
    parser.add_argument("--coverage-margin-um", type=float, default=0.0)
    return parser.parse_args()


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise RuntimeError(f"empty ranking CSV: {path}")
    return rows


def canonical_row(row: dict[str, str]) -> dict[str, Any]:
    normalized: dict[str, Any] = dict(row)
    for field in ACTIVE_CANDIDATE_PARAMETER_FIELDS:
        raw_feature = normalized.get(f"x_raw__{field}")
        if normalized.get(field) in (None, "") and raw_feature not in (None, ""):
            normalized[field] = raw_feature
    active = effective_candidate_parameters(normalized)
    output: dict[str, Any] = {
        "candidate_id": str(row["candidate_id"]),
        "campaign_parent_id": row.get("campaign_parent_id", ""),
        "campaign_parent_family": row.get("campaign_parent_family", ""),
        "source_search_rank": row.get("search_rank", ""),
        "source_peak_temperature_K": row.get("y__peak_temperature_K", ""),
        "source_peak_prominence_MPa_sqrt_m": row.get("y__peak_prominence", ""),
    }
    output.update({field: active[field] for field in ACTIVE_CANDIDATE_PARAMETER_FIELDS})
    return output


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    if args.selected_count < 1:
        raise ValueError("selected count must be positive")
    if args.target_extension_um <= 0.0 or args.coverage_margin_um < 0.0:
        raise ValueError("extension must be positive and margin nonnegative")

    ranking = read_rows(args.source_ranking)
    if args.selection_column not in ranking[0]:
        raise KeyError(f"selection column is absent: {args.selection_column}")
    selected = [row for row in ranking if truthy(row.get(args.selection_column))]
    if "search_rank" in ranking[0]:
        selected.sort(key=lambda row: float(row.get("search_rank") or float("inf")))
    if len(selected) < args.selected_count:
        raise RuntimeError(
            f"ranking contains only {len(selected)} selected rows; "
            f"requested {args.selected_count}"
        )
    selected = selected[: args.selected_count]
    ids = [str(row["candidate_id"]) for row in selected]
    if len(set(ids)) != len(ids):
        raise RuntimeError("selected candidate IDs are not unique")

    loading_payload = json.loads(args.loading_map.read_text())
    coverage_um = loading_map_coverage_um(loading_payload)
    required_um = float(args.target_extension_um) + float(args.coverage_margin_um)
    if coverage_um + 1.0e-9 < required_um:
        raise RuntimeError(
            f"loading-map coverage {coverage_um:.9g} um is below the required "
            f"{required_um:.9g} um; generate a longer mechanically calibrated map"
        )

    output_rows = [canonical_row(row) for row in selected]
    write_csv(args.out_registry, output_rows)
    ids_path = args.out_registry.with_suffix(".ids.txt")
    ids_path.write_text("\n".join(ids) + "\n")

    manifest_path = args.manifest or args.out_registry.with_suffix(".manifest.json")
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest = {
        "schema": "v9.13_long_peak25_preparation_v1",
        "source_ranking": str(args.source_ranking.resolve()),
        "source_ranking_sha256": sha256_path(args.source_ranking),
        "selection_column": args.selection_column,
        "selected_count": len(output_rows),
        "candidate_ids": ids,
        "loading_map": str(args.loading_map.resolve()),
        "loading_map_sha256": sha256_path(args.loading_map),
        "loading_map_coverage_um": coverage_um,
        "target_extension_um": args.target_extension_um,
        "coverage_margin_um": args.coverage_margin_um,
        "out_registry": str(args.out_registry.resolve()),
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n")
    print(
        "V913_LONG_PEAK25_PREPARED "
        f"selected={len(output_rows)} coverage_um={coverage_um:.9g} "
        f"target_um={args.target_extension_um:g} out={args.out_registry}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
