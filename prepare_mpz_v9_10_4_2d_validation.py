#!/usr/bin/env python3
"""Prepare a dynamic 2-D validation matrix from the long-growth shortlist."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def _transition_points(row: pd.Series) -> list[float]:
    stored = str(row.get("refinement_transition_temperatures_K", "")).strip()
    if stored and stored.lower() != "nan":
        return [float(x) for x in json.loads(stored)]
    low = float(row.get("coarse_transition_low_T_K", row["init_transition_low_K"]))
    high = float(row.get("coarse_transition_high_T_K", row["init_transition_high_K"]))
    return [float(x) for x in np.linspace(low, high, 4)]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-manifest", type=Path, required=True)
    ap.add_argument("--candidate-count", type=int, default=3)
    ap.add_argument("--temperature-min", type=float, default=300.0)
    ap.add_argument("--temperature-max", type=float, default=1100.0)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    candidates = pd.read_csv(args.input_manifest).head(args.candidate_count).copy()
    rows = []
    for _, row in candidates.iterrows():
        transition = _transition_points(row)
        low = min(transition)
        high = max(transition)
        below = max(float(args.temperature_min), low - 100.0)
        above = min(float(args.temperature_max), high + 100.0)
        schedule = [(below, "below_transition")]
        schedule.extend(
            (T, f"transition_point_{index + 1}_of_{len(transition)}")
            for index, T in enumerate(transition)
        )
        schedule.append((above, "above_transition"))
        seen: set[float] = set()
        for T, role in schedule:
            key = round(float(T), 10)
            if key in seen:
                continue
            seen.add(key)
            rows.append(
                {
                    "candidate_id": row["candidate_id"],
                    "T_K": float(T),
                    "validation_role": role,
                    **{key: row[key] for key in row.index if key not in {"T_K", "validation_role"}},
                }
            )
    out = args.out.resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, index=False)
    report = {
        "status": "V9_10_4_3_DYNAMIC_2D_VALIDATION_MANIFEST_READY",
        "n_candidates": int(len(candidates)),
        "n_cases": int(len(rows)),
        "output": str(out),
        "temperature_selection": "one shelf point below + four candidate-specific transition points + one shelf point above",
        "required_modes": ["full", "plasticity_off", "backstress_off", "shielding_off", "blunting_off"],
    }
    out.with_suffix(".json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
