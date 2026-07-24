#!/usr/bin/env python3
"""Prepare a quality-diverse v9.13 1-D screen from corrected zero-D promotions."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import numpy as np
import pandas as pd

from arrhenius_fracture.emergent_gnd_contract_v913 import (
    ACTIVE_CANDIDATE_PARAMETER_FIELDS,
    effective_candidate_parameters,
)


PROMOTABLE_TIERS = ("strict_gate", "relaxed_desired_peak")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-registry", type=Path, required=True)
    parser.add_argument("--out-registry", type=Path, required=True)
    parser.add_argument("--selected-count", type=int, default=384)
    return parser.parse_args()


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _numeric(frame: pd.DataFrame, name: str, default: float) -> pd.Series:
    if name not in frame:
        return pd.Series(np.full(len(frame), default), index=frame.index, dtype=float)
    return pd.to_numeric(frame[name], errors="coerce").fillna(default)


def _completion_mask(source: pd.DataFrame) -> pd.Series:
    """Return completion status for full metrics or corrected registry schemas.

    The corrected downstream registry deliberately omits ``zeroD_complete`` because
    its repair step already excludes incomplete rows before writing promotion tiers.
    When the column is present, retain the stronger row-level check so legacy or
    synthetic tables cannot promote an incomplete case.
    """
    if "zeroD_complete" in source:
        return _numeric(source, "zeroD_complete", 0.0).astype(bool)
    return pd.Series(np.ones(len(source), dtype=bool), index=source.index)


def select_rows(source: pd.DataFrame, selected_count: int) -> pd.DataFrame:
    if selected_count < 1:
        raise ValueError("selected count must be positive")
    required = {
        "candidate_id",
        "promotion_tier",
        *ACTIVE_CANDIDATE_PARAMETER_FIELDS,
    }
    missing = sorted(required - set(source.columns))
    if missing:
        raise RuntimeError(f"source registry is missing required columns: {missing}")
    if source["candidate_id"].astype(str).duplicated().any():
        raise RuntimeError("candidate_id must be unique")

    source = source.copy()
    tiers = source["promotion_tier"].astype(str)
    unknown_tiers = sorted(set(tiers) - set(PROMOTABLE_TIERS))
    if unknown_tiers:
        raise RuntimeError(
            "source registry contains non-promotable tiers: "
            f"{unknown_tiers}; expected only {list(PROMOTABLE_TIERS)}"
        )
    for row in source.to_dict(orient="records"):
        effective_candidate_parameters(row)

    complete = _completion_mask(source)
    strict = source[(tiers == "strict_gate") & complete].copy()
    relaxed = source[(tiers == "relaxed_desired_peak") & complete].copy()

    strict["_objective"] = _numeric(strict, "zeroD_objective", float("inf"))
    strict["_rank"] = _numeric(strict, "zeroD_rank", float("inf"))
    strict = strict.sort_values(
        ["_objective", "_rank", "candidate_id"],
        kind="stable",
    )

    relaxed["_diversity"] = _numeric(relaxed, "diversity_rank", float("inf"))
    relaxed["_objective"] = _numeric(relaxed, "zeroD_objective", float("inf"))
    relaxed = relaxed.sort_values(
        ["_diversity", "_objective", "candidate_id"],
        kind="stable",
    )

    if len(strict) > selected_count:
        raise RuntimeError(
            f"selected count {selected_count} is smaller than strict-pass count {len(strict)}"
        )
    relaxed_needed = selected_count - len(strict)
    if len(relaxed) < relaxed_needed:
        raise RuntimeError(
            "corrected registry does not contain enough complete relaxed desired peaks: "
            f"needed={relaxed_needed}, available={len(relaxed)}"
        )

    strict = strict.drop(columns=["_objective", "_rank"])
    relaxed = relaxed.head(relaxed_needed).drop(columns=["_diversity", "_objective"])
    selected = pd.concat([strict, relaxed], ignore_index=True, sort=False)
    selected.insert(0, "oneD_selection_rank", np.arange(1, len(selected) + 1))
    selected.insert(
        1,
        "oneD_selection_tier",
        np.where(
            selected["promotion_tier"].astype(str) == "strict_gate",
            "strict_zeroD",
            "relaxed_diverse_zeroD",
        ),
    )
    return selected


def main() -> int:
    args = parse_args()
    if not args.source_registry.is_file():
        raise FileNotFoundError(args.source_registry)
    source = pd.read_csv(args.source_registry)
    selected = select_rows(source, int(args.selected_count))

    args.out_registry.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.out_registry.with_suffix(args.out_registry.suffix + ".tmp")
    selected.to_csv(temporary, index=False)
    temporary.replace(args.out_registry)

    strict_count = int(
        (selected["oneD_selection_tier"].astype(str) == "strict_zeroD").sum()
    )
    relaxed_count = len(selected) - strict_count
    manifest: dict[str, Any] = {
        "schema": "v9.13_zero_d_to_one_d_screen_selection_v2",
        "source_registry": str(args.source_registry.resolve()),
        "source_registry_sha256": sha256_path(args.source_registry),
        "out_registry": str(args.out_registry.resolve()),
        "out_registry_sha256": sha256_path(args.out_registry),
        "source_count": int(len(source)),
        "selected_count": int(len(selected)),
        "strict_zeroD_count": strict_count,
        "relaxed_diverse_zeroD_count": relaxed_count,
        "source_completion_column_present": "zeroD_complete" in source.columns,
        "completion_validation": (
            "explicit_zeroD_complete"
            if "zeroD_complete" in source.columns
            else "corrected_promotion_tier_contract"
        ),
        "all_source_strict_passes_selected": bool(
            set(
                source.loc[
                    source["promotion_tier"].astype(str) == "strict_gate",
                    "candidate_id",
                ].astype(str)
            ).issubset(set(selected["candidate_id"].astype(str)))
        ),
        "selection_order": [
            "all complete strict gates ordered by zeroD objective",
            "complete relaxed desired peaks ordered by corrected diversity rank",
        ],
    }
    manifest_path = args.out_registry.with_suffix(".manifest.json")
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")

    if not manifest["all_source_strict_passes_selected"]:
        raise RuntimeError("not all strict zero-D passes were selected")
    print(
        "V913_ZERO_D_TO_1D_PREPARED "
        f"selected={len(selected)} strict={strict_count} relaxed={relaxed_count} "
        f"completion={manifest['completion_validation']} "
        f"registry={args.out_registry}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
