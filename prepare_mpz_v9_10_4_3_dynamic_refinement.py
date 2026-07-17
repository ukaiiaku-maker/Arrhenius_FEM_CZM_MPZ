#!/usr/bin/env python3
"""Attach candidate-specific four-point DBTT refinement schedules.

Input is the coarse first-passage shortlist evaluated on the complete 300--1100
K grid.  Each row already contains the best adjacent coarse transition bracket.
This script adds four temperatures across that bracket and broad low/high shelf
anchors for the refined first-passage and growth stages.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from arrhenius_fracture.dbtt_temperature_schedule_v91043 import schedule_from_bracket


def parse_floats(text: str) -> list[float]:
    return [float(x) for x in str(text).replace(",", " ").split() if x]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-manifest", type=Path, required=True)
    ap.add_argument("--coarse-temperatures", default="300 400 500 600 700 800 900 1000 1100")
    ap.add_argument("--refinement-points", type=int, default=4)
    ap.add_argument("--shelf-anchor-count", type=int, default=2)
    ap.add_argument("--max-candidates", type=int, default=24)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    coarse = parse_floats(args.coarse_temperatures)
    frame = pd.read_csv(args.input_manifest).head(args.max_candidates).copy()
    rows = []
    rejected = []
    for _, row in frame.iterrows():
        cid = str(row.get("candidate_id", "candidate"))
        low = row.get("coarse_transition_low_T_K", row.get("transition_transition_low_K"))
        high = row.get("coarse_transition_high_T_K", row.get("transition_transition_high_K"))
        try:
            schedule = schedule_from_bracket(
                coarse,
                float(low),
                float(high),
                refinement_points=args.refinement_points,
                shelf_anchor_count=args.shelf_anchor_count,
            )
        except Exception as exc:
            rejected.append({"candidate_id": cid, "reason": str(exc)})
            print(f"[schedule-rejected] candidate={cid} reason={exc}", flush=True)
            continue
        data = row.to_dict()
        data.update(schedule.to_columns())
        data["temperature_schedule_mode"] = "COARSE_100K_PLUS_DYNAMIC_4POINT_BRACKET"
        rows.append(data)
        print(
            f"[schedule] candidate={cid} bracket={schedule.transition_low_K:g}-"
            f"{schedule.transition_high_K:g}K transition={list(schedule.transition_temperatures_K)} "
            f"evaluation={list(schedule.evaluation_temperatures_K)}",
            flush=True,
        )

    out = args.out.resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, index=False)
    rejected_path = out.with_name(out.stem + "_rejected.csv")
    pd.DataFrame(rejected).to_csv(rejected_path, index=False)
    summary = {
        "status": "V9_10_4_3_DYNAMIC_REFINEMENT_MANIFEST_COMPLETE",
        "input_manifest": str(args.input_manifest),
        "output_manifest": str(out),
        "n_input": int(len(frame)),
        "n_scheduled": int(len(rows)),
        "n_rejected": int(len(rejected)),
        "coarse_temperatures_K": coarse,
        "refinement_points": int(args.refinement_points),
        "shelf_anchor_count": int(args.shelf_anchor_count),
    }
    out.with_suffix(".summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
