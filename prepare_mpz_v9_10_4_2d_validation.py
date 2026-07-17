#!/usr/bin/env python3
"""Prepare a compact 2-D validation matrix from the long-growth shortlist."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-manifest", type=Path, required=True)
    ap.add_argument("--candidate-count", type=int, default=3)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    candidates = pd.read_csv(args.input_manifest).head(args.candidate_count).copy()
    rows = []
    for _, row in candidates.iterrows():
        low = int(round(float(row["init_transition_low_K"])))
        high = int(round(float(row["init_transition_high_K"])))
        temps = sorted({max(300, low - 100), low, high, min(1100, high + 100)})
        for T in temps:
            rows.append(
                {
                    "candidate_id": row["candidate_id"],
                    "T_K": T,
                    "validation_role": (
                        "below_transition" if T < low else
                        "lower_bracket" if T == low else
                        "upper_bracket" if T == high else
                        "above_transition"
                    ),
                    **{key: row[key] for key in row.index if key not in {"T_K", "validation_role"}},
                }
            )
    out = args.out.resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, index=False)
    report = {
        "status": "V9_10_4_2D_VALIDATION_MANIFEST_READY",
        "n_candidates": int(len(candidates)),
        "n_cases": int(len(rows)),
        "output": str(out),
        "required_modes": ["full", "plasticity_off", "backstress_off", "shielding_off", "blunting_off"],
    }
    out.with_suffix(".json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
