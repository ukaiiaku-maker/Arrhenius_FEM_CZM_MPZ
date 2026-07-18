#!/usr/bin/env python3
"""Rescore completed v9.10.4.7 four-temperature results with the v9.10.4.8 gate."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from evaluate_dynamic_1d_mpz_v9_10_4_8 import corrected_transition_metrics


LEGACY_COLUMNS = (
    "moving_1d_objective",
    "moving_1d_accept",
    "moving_1d_reason",
    "moving_1d_robust_ratio",
)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--per-bracket-keep", type=int, default=2)
    args = ap.parse_args()

    source = pd.read_csv(args.input)
    rescored = source.copy()
    for name in LEGACY_COLUMNS:
        if name in rescored.columns:
            rescored[f"legacy_{name}"] = rescored[name]

    records: list[dict[str, object]] = []
    for _, row in rescored.iterrows():
        T = np.asarray(json.loads(str(row.moving_1d_temperatures_K)), dtype=float)
        full = np.asarray(json.loads(str(row.moving_1d_full_K_json)), dtype=float)
        off = np.asarray(json.loads(str(row.moving_1d_off_K_json)), dtype=float)
        records.append(corrected_transition_metrics(T, full, off))
    metric_df = pd.DataFrame(records, index=rescored.index)
    for column in metric_df.columns:
        rescored[column] = metric_df[column]

    args.out.mkdir(parents=True, exist_ok=True)
    rescored = rescored.sort_values(
        ["coarse_transition_low_T_K", "moving_1d_objective"]
    )
    accepted = rescored[rescored.moving_1d_accept.fillna(False).astype(bool)].copy()

    promoted: list[pd.DataFrame] = []
    for bracket, group in accepted.groupby("transition_bracket", sort=True):
        keep = group.sort_values("moving_1d_objective").head(args.per_bracket_keep).copy()
        keep["moving_1d_selection_basis"] = "corrected_endpoint_gate"
        keep["moving_1d_rank_within_bracket"] = np.arange(1, len(keep) + 1)
        promoted.append(keep)
    promotion = pd.concat(promoted, ignore_index=True) if promoted else accepted.head(0).copy()

    rescored.to_csv(args.out / "dynamic_1d_all_candidates_rescored.csv", index=False)
    accepted.to_csv(args.out / "dynamic_1d_accepted_rescored.csv", index=False)
    promotion.to_csv(args.out / "short_growth_promotion_manifest.csv", index=False)

    review = rescored.sort_values("moving_1d_objective").copy()
    review.to_csv(args.out / "dynamic_1d_candidate_review.csv", index=False)

    bracket_rows: list[dict[str, object]] = []
    for bracket, group in rescored.groupby("transition_bracket", sort=True):
        group_accepted = group[group.moving_1d_accept.fillna(False).astype(bool)]
        bracket_rows.append(
            {
                "transition_bracket": bracket,
                "n_candidates": int(len(group)),
                "n_accepted": int(len(group_accepted)),
                "best_available_objective": float(group.moving_1d_objective.min()),
                "best_accepted_objective": (
                    float(group_accepted.moving_1d_objective.min())
                    if not group_accepted.empty
                    else float("nan")
                ),
            }
        )
    pd.DataFrame(bracket_rows).to_csv(
        args.out / "dynamic_1d_bracket_summary_rescored.csv", index=False
    )

    report = {
        "status": "V9_10_4_8_DYNAMIC_1D_RESCORE_COMPLETE",
        "source": str(args.input),
        "n_candidates": int(len(rescored)),
        "n_accepted": int(len(accepted)),
        "n_promoted": int(len(promotion)),
        "n_accepted_brackets": int(accepted.transition_bracket.nunique()) if not accepted.empty else 0,
        "gate": "ENDPOINT_FACTOR_TWO_OVER_CANDIDATE_SPECIFIC_100K_BRACKET",
        "next_stage_manifest": str(args.out / "short_growth_promotion_manifest.csv"),
    }
    (args.out / "dynamic_1d_rescore_summary.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2), flush=True)
    if not accepted.empty:
        columns = [
            "candidate_id",
            "transition_bracket",
            "moving_1d_low_endpoint_K",
            "moving_1d_high_endpoint_K",
            "moving_1d_edge_ratio",
            "moving_1d_transition_width_K",
            "moving_1d_plasticity_off_ratio",
            "moving_1d_mechanistic_fraction",
            "moving_1d_reason",
        ]
        print(accepted[columns].to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
