#!/usr/bin/env python3
"""Candidate-specific four-point first-passage refinement for narrow DBTT.

Each input candidate was first evaluated on the complete 100 K coarse grid.
The candidate's selected adjacent coarse bracket is held fixed during local
refinement. Four points resolve that bracket, while broad shelf anchors retain
the factor-of-two and shelf-flatness constraints.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.optimize import minimize

import optimize_mpz_v9_10_4_narrow_dbtt as first
from arrhenius_fracture.dbtt_temperature_schedule_v91043 import (
    DynamicTemperatureSchedule,
    schedule_from_candidate_row,
)
from arrhenius_fracture.reduced_campaign_front_v9104 import (
    ReducedFrontSettings,
    TransitionRequirements,
)

PARAMETER_NAMES = first.PARAMETER_NAMES


class WindowedFirstPassageObjective:
    def __init__(
        self,
        schedule: DynamicTemperatureSchedule,
        settings: ReducedFrontSettings,
        *,
        cleavage_slope_mode: str,
    ) -> None:
        self.schedule = schedule
        self.base = first.NarrowDBTTObjective(
            np.asarray(schedule.evaluation_temperatures_K, dtype=float),
            settings,
            cleavage_slope_mode=cleavage_slope_mode,
            transition_requirements=TransitionRequirements(),
        )
        self.bounds = [self.base.bounds_dict[name] for name in PARAMETER_NAMES]

    def evaluate(self, x: np.ndarray, *, details: bool = False) -> dict[str, Any]:
        result = dict(self.base.evaluate(np.asarray(x, dtype=float), details=details))
        selected_low = float(result.get("transition_transition_low_K", np.nan))
        selected_high = float(result.get("transition_transition_high_K", np.nan))
        low = float(self.schedule.transition_low_K)
        high = float(self.schedule.transition_high_K)
        span = max(high - low, 1.0e-12)
        finite = np.isfinite(selected_low) and np.isfinite(selected_high)
        in_window = bool(
            finite
            and selected_low >= low - 1.0e-8
            and selected_high <= high + 1.0e-8
        )
        if finite:
            distance = (
                max(low - selected_low, 0.0)
                + max(selected_high - high, 0.0)
            ) / span
        else:
            distance = 10.0
        window_loss = 0.0 if in_window else 100.0 + 100.0 * distance * distance
        result["coarse_transition_low_T_K"] = low
        result["coarse_transition_high_T_K"] = high
        result["refined_transition_in_coarse_bracket"] = in_window
        result["refined_transition_window_loss"] = float(window_loss)
        result.update(self.schedule.to_columns())
        result["objective"] = float(result.get("objective", 1.0e12) + window_loss)
        result["objective_mode"] = "DYNAMIC_4POINT_DBTT_FIRST_PASSAGE_REFINEMENT"
        return result

    def __call__(self, x: np.ndarray) -> float:
        return float(self.evaluate(x, details=False)["objective"])


def accepted(detail: dict[str, Any]) -> tuple[bool, str]:
    base_passed, reason = first.accepted(detail)
    if not base_passed:
        return False, reason
    if not bool(detail.get("refined_transition_in_coarse_bracket", False)):
        return False, "refined_transition_left_coarse_bracket"
    return True, "dynamic_four_point_first_passage_gate_passed"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-manifest", type=Path, required=True)
    ap.add_argument("--coarse-temperatures", default="300 400 500 600 700 800 900 1000 1100")
    ap.add_argument("--cleavage-slope-mode", choices=["fixed_zero", "narrow"], default="fixed_zero")
    ap.add_argument("--max-candidates", type=int, default=24)
    ap.add_argument("--local-maxiter", type=int, default=150)
    ap.add_argument("--Kdot", type=float, default=0.005)
    ap.add_argument("--Kmax", type=float, default=80.0)
    ap.add_argument("--target-extension-um", type=float, default=5.0)
    ap.add_argument("--max-dK-substep", type=float, default=0.05)
    ap.add_argument("--max-K-shield", type=float, default=1.0)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    coarse = first.parse_floats(args.coarse_temperatures)
    settings = ReducedFrontSettings(
        Kdot_MPa_sqrt_m_s=float(args.Kdot),
        Kmax_MPa_sqrt_m=float(args.Kmax),
        max_dK_substep_MPa_sqrt_m=float(args.max_dK_substep),
        target_extension_um=float(args.target_extension_um),
        max_K_shield_MPa_sqrt_m=float(args.max_K_shield),
    )
    manifest = pd.read_csv(args.input_manifest).head(args.max_candidates)
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, Any]] = []
    temperature_rows: list[dict[str, Any]] = []
    event_rows: list[dict[str, Any]] = []

    for rank, row in manifest.iterrows():
        parent_id = str(row.get("candidate_id", rank))
        schedule = schedule_from_candidate_row(row, coarse)
        objective = WindowedFirstPassageObjective(
            schedule,
            settings,
            cleavage_slope_mode=args.cleavage_slope_mode,
        )
        x0 = np.asarray([float(row[name]) for name in PARAMETER_NAMES], dtype=float)
        print(
            f"[candidate-start] parent={parent_id} bracket="
            f"{schedule.transition_low_K:g}-{schedule.transition_high_K:g}K "
            f"evaluation={list(schedule.evaluation_temperatures_K)}",
            flush=True,
        )
        local = minimize(
            objective,
            x0,
            method="Powell",
            bounds=objective.bounds,
            options={"maxiter": args.local_maxiter, "xtol": 1.0e-4, "ftol": 1.0e-5},
        )
        candidates = [
            (objective(x0), x0, "coarse_input"),
            (float(local.fun), np.asarray(local.x, dtype=float), "dynamic_local_refined"),
        ]
        for source_index, (_, x, source) in enumerate(sorted(candidates, key=lambda item: item[0])):
            detail = objective.evaluate(x, details=True)
            p = detail.pop("parameters", None) or first.decode(x)
            tdetail = detail.pop("temperature_detail", None) or []
            edetail = detail.pop("event_detail", None) or []
            candidate_id = f"DBTT_v91043_refined_{rank:02d}_{source_index:02d}"
            is_accepted, reason = accepted(detail)
            output_row = {
                "candidate_id": candidate_id,
                "parent_candidate_id": parent_id,
                "candidate_source": source,
                "stage": "refined_first_passage",
                "accepted_for_short_growth": bool(is_accepted),
                "acceptance_reason": reason,
                "temperature_schedule_mode": "COARSE_100K_PLUS_DYNAMIC_4POINT_BRACKET",
                **{name: float(x[i]) for i, name in enumerate(PARAMETER_NAMES)},
                **{key: float(value) for key, value in p.items() if key not in PARAMETER_NAMES},
                **detail,
            }
            rows.append(output_row)
            temperature_rows.extend({"candidate_id": candidate_id, **record} for record in tdetail)
            event_rows.extend({"candidate_id": candidate_id, **record} for record in edetail)
            print(
                f"[candidate-result] candidate={candidate_id} objective={output_row['objective']:.6g} "
                f"accepted={is_accepted} reason={reason}",
                flush=True,
            )

    results = pd.DataFrame(rows).sort_values("objective").drop_duplicates(PARAMETER_NAMES)
    accepted_df = results[results.accepted_for_short_growth.astype(bool)].copy()
    promotion = (accepted_df if not accepted_df.empty else results).head(12).copy()
    results.to_csv(out / "narrow_dbtt_refined_first_passage_all_candidates.csv", index=False)
    accepted_df.to_csv(out / "narrow_dbtt_refined_first_passage_accepted.csv", index=False)
    promotion.to_csv(out / "short_growth_promotion_manifest.csv", index=False)
    pd.DataFrame(temperature_rows).to_csv(out / "narrow_dbtt_refined_first_passage_temperature_detail.csv", index=False)
    pd.DataFrame(event_rows).to_csv(out / "narrow_dbtt_refined_first_passage_event_detail.csv", index=False)
    summary = {
        "status": "V9_10_4_3_DYNAMIC_FIRST_PASSAGE_REFINEMENT_COMPLETE",
        "n_candidates": int(len(results)),
        "n_accepted": int(len(accepted_df)),
        "target_extension_um": float(args.target_extension_um),
        "next_stage_manifest": str(out / "short_growth_promotion_manifest.csv"),
    }
    (out / "narrow_dbtt_refined_first_passage_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
