#!/usr/bin/env python3
"""Select a bracket-balanced candidate set for mechanism ablations.

The completed v9.10.4.8 four-temperature table is used directly.  Within each
transition bracket, strict endpoint-gate passes are ranked first, followed by
the lowest corrected moving-interface objective among the remaining valid
candidates.  The default produces one candidate from each occupied bracket.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def select_candidates(
    results: pd.DataFrame,
    *,
    per_bracket: int = 1,
    expected_brackets: int | None = 6,
) -> pd.DataFrame:
    required = {
        "candidate_id",
        "transition_bracket",
        "moving_1d_valid",
        "moving_1d_accept",
        "moving_1d_objective",
    }
    missing = sorted(required.difference(results.columns))
    if missing:
        raise ValueError(f"input is missing required columns: {missing}")

    valid = results[results.moving_1d_valid.fillna(False).astype(bool)].copy()
    if valid.empty:
        raise ValueError("no valid moving-interface candidates were found")

    selected: list[pd.DataFrame] = []
    for bracket, group in valid.groupby("transition_bracket", sort=True):
        passed = group[group.moving_1d_accept.fillna(False).astype(bool)].sort_values(
            ["moving_1d_objective", "candidate_id"]
        )
        remaining = group[~group.candidate_id.astype(str).isin(passed.candidate_id.astype(str))].sort_values(
            ["moving_1d_objective", "candidate_id"]
        )
        ordered = pd.concat([passed, remaining], ignore_index=False)
        keep = ordered.head(int(per_bracket)).copy()
        passed_ids = set(passed.candidate_id.astype(str))
        keep["ablation_selection_basis"] = [
            "strict_1d_pass" if str(cid) in passed_ids else "best_valid_near_miss"
            for cid in keep.candidate_id
        ]
        keep["ablation_rank_within_bracket"] = np.arange(1, len(keep) + 1)
        selected.append(keep)

    out = pd.concat(selected, ignore_index=True).sort_values(
        ["coarse_transition_low_T_K", "ablation_rank_within_bracket"]
    )
    n_brackets = int(out.transition_bracket.nunique())
    if expected_brackets is not None and n_brackets != int(expected_brackets):
        raise ValueError(
            f"expected {expected_brackets} occupied transition brackets; found {n_brackets}"
        )
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--per-bracket", type=int, default=1)
    ap.add_argument("--expected-brackets", type=int, default=6)
    args = ap.parse_args()

    results = pd.read_csv(args.input)
    selected = select_candidates(
        results,
        per_bracket=args.per_bracket,
        expected_brackets=args.expected_brackets,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    selected.to_csv(args.out, index=False)

    columns = [
        "candidate_id",
        "transition_bracket",
        "ablation_selection_basis",
        "moving_1d_accept",
        "moving_1d_objective",
        "moving_1d_reason",
        "moving_1d_low_endpoint_K",
        "moving_1d_high_endpoint_K",
        "moving_1d_edge_ratio",
        "refinement_transition_temperatures_K",
    ]
    print(selected[[c for c in columns if c in selected.columns]].to_string(index=False), flush=True)
    print(f"wrote: {args.out}", flush=True)


if __name__ == "__main__":
    main()
