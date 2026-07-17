#!/usr/bin/env python3
"""Dynamic-schedule short- and long-growth DBTT refinement.

Each promoted candidate retains the 100 K coarse bracket selected during the
broad first-passage sweep.  Four temperatures resolve that bracket, and broad
shelf anchors remain in the objective.  Initiation and plateau transitions must
both remain inside the candidate's bracket.
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
import refine_mpz_v9_10_4_growth as legacy_growth
from arrhenius_fracture.dbtt_temperature_schedule_v91043 import (
    DynamicTemperatureSchedule,
    schedule_from_candidate_row,
)
from arrhenius_fracture.reduced_campaign_front_v9104 import ReducedFrontSettings

PARAMETER_NAMES = first.PARAMETER_NAMES


class WindowedGrowthObjective:
    def __init__(
        self,
        schedule: DynamicTemperatureSchedule,
        settings: ReducedFrontSettings,
        *,
        cleavage_slope_mode: str,
    ) -> None:
        self.schedule = schedule
        self.base = legacy_growth.GrowthObjective(
            np.asarray(schedule.evaluation_temperatures_K, dtype=float),
            settings,
            cleavage_slope_mode=cleavage_slope_mode,
        )
        self.bounds = self.base.bounds

    def _inside(self, low_value: Any, high_value: Any) -> bool:
        try:
            selected_low = float(low_value)
            selected_high = float(high_value)
        except (TypeError, ValueError):
            return False
        return bool(
            np.isfinite(selected_low)
            and np.isfinite(selected_high)
            and selected_low >= self.schedule.transition_low_K - 1.0e-8
            and selected_high <= self.schedule.transition_high_K + 1.0e-8
        )

    def evaluate(self, x: np.ndarray, *, details: bool = False) -> dict[str, Any]:
        result = dict(self.base.evaluate(np.asarray(x, dtype=float), details=details))
        init_inside = self._inside(
            result.get("init_transition_low_K"),
            result.get("init_transition_high_K"),
        )
        plateau_inside = self._inside(
            result.get("plateau_transition_low_K"),
            result.get("plateau_transition_high_K"),
        )
        window_loss = 0.0 if init_inside and plateau_inside else 200.0
        result["init_transition_in_coarse_bracket"] = init_inside
        result["plateau_transition_in_coarse_bracket"] = plateau_inside
        result["dynamic_transition_window_loss"] = float(window_loss)
        result.update(self.schedule.to_columns())
        result["objective"] = float(result.get("objective", 1.0e12) + window_loss)
        result["objective_mode"] = "DYNAMIC_4POINT_DBTT_GROWTH_REFINEMENT"
        return result

    def __call__(self, x: np.ndarray) -> float:
        return float(self.evaluate(x, details=False)["objective"])


def growth_acceptance(detail: dict[str, Any]) -> tuple[bool, str]:
    passed, reason = legacy_growth.growth_acceptance(detail)
    if not passed:
        return False, reason
    if not bool(detail.get("init_transition_in_coarse_bracket", False)):
        return False, "initiation_transition_left_coarse_bracket"
    if not bool(detail.get("plateau_transition_in_coarse_bracket", False)):
        return False, "plateau_transition_left_coarse_bracket"
    return True, "dynamic_four_point_growth_gate_passed"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-manifest", type=Path, required=True)
    ap.add_argument("--stage", choices=["short", "long"], required=True)
    ap.add_argument("--coarse-temperatures", default="300 400 500 600 700 800 900 1000 1100")
    ap.add_argument("--target-extension-um", type=float, default=None)
    ap.add_argument("--cleavage-slope-mode", choices=["fixed_zero", "narrow"], default="fixed_zero")
    ap.add_argument("--max-candidates", type=int, default=12)
    ap.add_argument("--local-maxiter", type=int, default=100)
    ap.add_argument("--Kdot", type=float, default=0.005)
    ap.add_argument("--Kmax", type=float, default=80.0)
    ap.add_argument("--max-dK-substep", type=float, default=0.05)
    ap.add_argument("--max-K-shield", type=float, default=1.0)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    extension = (
        float(args.target_extension_um)
        if args.target_extension_um is not None
        else (100.0 if args.stage == "short" else 500.0)
    )
    coarse = first.parse_floats(args.coarse_temperatures)
    settings = ReducedFrontSettings(
        Kdot_MPa_sqrt_m_s=float(args.Kdot),
        Kmax_MPa_sqrt_m=float(args.Kmax),
        max_dK_substep_MPa_sqrt_m=float(args.max_dK_substep),
        target_extension_um=extension,
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
        objective = WindowedGrowthObjective(
            schedule,
            settings,
            cleavage_slope_mode=args.cleavage_slope_mode,
        )
        x0 = np.asarray([float(row[name]) for name in PARAMETER_NAMES], dtype=float)
        print(
            f"[candidate-start] stage={args.stage} parent={parent_id} bracket="
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
            (objective(x0), x0, "input"),
            (float(local.fun), np.asarray(local.x), "dynamic_local_refined"),
        ]
        for source_index, (_, x, source) in enumerate(sorted(candidates, key=lambda item: item[0])):
            detail = objective.evaluate(x, details=True)
            p = detail.pop("parameters", first.decode(x))
            tdetail = detail.pop("temperature_detail", [])
            edetail = detail.pop("event_detail", [])
            candidate_id = f"DBTT_v91043_{args.stage}_{rank:02d}_{source_index:02d}"
            accepted, reason = growth_acceptance(detail)
            output_row = {
                "candidate_id": candidate_id,
                "parent_candidate_id": parent_id,
                "candidate_source": source,
                "stage": args.stage,
                "target_extension_um": extension,
                "accepted_for_next_stage": bool(accepted),
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
                f"[candidate-result] stage={args.stage} candidate={candidate_id} "
                f"objective={output_row['objective']:.6g} accepted={accepted} reason={reason}",
                flush=True,
            )

    results = pd.DataFrame(rows).sort_values("objective").drop_duplicates(PARAMETER_NAMES)
    accepted_df = results[results.accepted_for_next_stage.astype(bool)].copy()
    promotion = (accepted_df if not accepted_df.empty else results).head(8).copy()
    prefix = f"narrow_dbtt_v91043_{args.stage}_growth"
    results.to_csv(out / f"{prefix}_all_candidates.csv", index=False)
    accepted_df.to_csv(out / f"{prefix}_accepted.csv", index=False)
    promotion.to_csv(out / f"{prefix}_promotion_manifest.csv", index=False)
    pd.DataFrame(temperature_rows).to_csv(out / f"{prefix}_temperature_detail.csv", index=False)
    pd.DataFrame(event_rows).to_csv(out / f"{prefix}_event_detail.csv", index=False)
    summary = {
        "status": f"V9_10_4_3_DYNAMIC_{args.stage.upper()}_GROWTH_COMPLETE",
        "target_extension_um": extension,
        "n_candidates": int(len(results)),
        "n_accepted": int(len(accepted_df)),
        "promotion_manifest": str(out / f"{prefix}_promotion_manifest.csv"),
    }
    (out / f"{prefix}_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
