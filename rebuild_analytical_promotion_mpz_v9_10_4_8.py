#!/usr/bin/env python3
"""Rebuild an analytical promotion manifest without rerunning the Sobol screen.

Within each occupied transition bracket, analytical screen passes are retained
first and the remaining positions are filled by the lowest-objective valid
candidates.  This fixes the v9.10.4.7 behavior that returned fewer than
``per_bracket_keep`` candidates whenever a bracket contained at least one pass.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def rebuild(results: pd.DataFrame, per_bracket_keep: int) -> pd.DataFrame:
    valid = results[results.analysis_valid.fillna(False).astype(bool)].copy()
    selected: list[pd.DataFrame] = []
    for bracket, group in valid.groupby("transition_bracket", sort=True):
        passed = group[group.screen_pass.fillna(False).astype(bool)].sort_values("objective")
        remaining = group[~group.candidate_id.isin(passed.candidate_id)].sort_values("objective")
        ordered = pd.concat([passed, remaining], ignore_index=False)
        keep = ordered.head(int(per_bracket_keep)).copy()
        pass_ids = set(passed.candidate_id.astype(str))
        keep["selection_basis"] = [
            "analytical_screen_pass" if str(cid) in pass_ids else "best_valid_fill"
            for cid in keep.candidate_id
        ]
        keep["rank_within_bracket"] = np.arange(1, len(keep) + 1)
        selected.append(keep)
    if not selected:
        return valid.head(0).copy()
    return pd.concat(selected, ignore_index=True).sort_values(
        ["coarse_transition_low_T_K", "rank_within_bracket"]
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--per-bracket-keep", type=int, default=8)
    args = ap.parse_args()

    results = pd.read_csv(args.input)
    selected = rebuild(results, args.per_bracket_keep)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    selected.to_csv(args.out, index=False)

    queue_rows: list[dict[str, object]] = []
    for _, row in selected.iterrows():
        schedule = [float(x) for x in json.loads(str(row.refinement_transition_temperatures_K))]
        for order, T in enumerate(schedule, start=1):
            queue_rows.append(
                {
                    "candidate_id": row.candidate_id,
                    "transition_bracket": row.transition_bracket,
                    "temperature_order": order,
                    "T_K": T,
                    "selection_basis": row.selection_basis,
                }
            )
    queue_path = args.out.with_name(args.out.stem + "_queue.csv")
    pd.DataFrame(queue_rows).to_csv(queue_path, index=False)

    summary = (
        selected.groupby("transition_bracket", as_index=False)
        .agg(
            n_promoted=("candidate_id", "count"),
            n_screen_pass=("screen_pass", "sum"),
            best_objective=("objective", "min"),
        )
    )
    summary_path = args.out.with_name(args.out.stem + "_bracket_summary.csv")
    summary.to_csv(summary_path, index=False)
    print(summary.to_string(index=False), flush=True)
    print(f"wrote manifest: {args.out}", flush=True)
    print(f"wrote queue: {queue_path}", flush=True)


if __name__ == "__main__":
    main()
