#!/usr/bin/env python3
"""Dynamic-schedule short- and long-growth DBTT refinement.

Each promoted candidate retains the 100 K coarse bracket selected during the
broad first-passage sweep. Four temperatures resolve the complete bracket, and
broad shelf anchors remain in the objective. Initiation and plateau are scored
against that fixed bracket rather than one ~33 K refined subinterval.
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
    fixed_bracket_transition_metrics,
    schedule_from_candidate_row,
)
from arrhenius_fracture.reduced_campaign_front_v9104 import ReducedFrontSettings

PARAMETER_NAMES = first.PARAMETER_NAMES


def _values_at(frame: pd.DataFrame, temperatures: tuple[float, ...], column: str) -> np.ndarray:
    values = []
    for temperature in temperatures:
        matches = frame.loc[np.isclose(frame.T_K, temperature, rtol=0.0, atol=1.0e-7), column]
        if len(matches) != 1:
            raise ValueError(f"expected one {column} value at {temperature} K; found {len(matches)}")
        values.append(float(matches.iloc[0]))
    return np.asarray(values, dtype=float)


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

    def evaluate(self, x: np.ndarray, *, details: bool = False) -> dict[str, Any]:
        x = np.asarray(x, dtype=float)
        base_result = dict(self.base.evaluate(x, details=True))
        schedule_columns = self.schedule.to_columns()
        temperature_detail = base_result.get("temperature_detail") or []
        parameters = base_result.get("parameters")

        if not temperature_detail or parameters is None:
            result = dict(base_result)
            result.update(schedule_columns)
            result["init_transition_in_coarse_bracket"] = False
            result["plateau_transition_in_coarse_bracket"] = False
            result["objective_mode"] = "DYNAMIC_4POINT_DBTT_GROWTH_REFINEMENT"
            if not details:
                result.pop("temperature_detail", None)
                result.pop("event_detail", None)
                result.pop("parameters", None)
            return result

        frame = pd.DataFrame(temperature_detail).sort_values("T_K")
        required = {
            "T_K", "completed", "K_init_proxy", "K_plateau_proxy",
            "K_init_plasticity_off", "K_plateau_plasticity_off",
        }
        if not required.issubset(frame.columns) or not frame.completed.astype(bool).all():
            result = {
                "objective": 1.0e6,
                "completion_loss": 1.0e5,
                "parameters": parameters,
                "temperature_detail": temperature_detail,
                "event_detail": base_result.get("event_detail") or [],
                "init_transition_in_coarse_bracket": False,
                "plateau_transition_in_coarse_bracket": False,
                "objective_mode": "DYNAMIC_4POINT_DBTT_GROWTH_REFINEMENT",
                **schedule_columns,
            }
            if not details:
                result.pop("temperature_detail", None)
                result.pop("event_detail", None)
                result.pop("parameters", None)
            return result

        init_transition = fixed_bracket_transition_metrics(
            frame.T_K,
            frame.K_init_proxy,
            self.schedule,
            plasticity_off_toughness=frame.K_init_plasticity_off,
        )
        plateau_transition = fixed_bracket_transition_metrics(
            frame.T_K,
            frame.K_plateau_proxy,
            self.schedule,
            plasticity_off_toughness=frame.K_plateau_plasticity_off,
        )
        if not bool(init_transition.get("valid", False)) or not bool(plateau_transition.get("valid", False)):
            result = {
                "objective": 1.0e12,
                "completion_loss": 0.0,
                "init_transition_reason": init_transition.get("reason", "invalid_transition"),
                "plateau_transition_reason": plateau_transition.get("reason", "invalid_transition"),
                "parameters": parameters,
                "temperature_detail": temperature_detail,
                "event_detail": base_result.get("event_detail") or [],
                "init_transition_in_coarse_bracket": False,
                "plateau_transition_in_coarse_bracket": False,
                "objective_mode": "DYNAMIC_4POINT_DBTT_GROWTH_REFINEMENT",
                **schedule_columns,
            }
            if not details:
                result.pop("temperature_detail", None)
                result.pop("event_detail", None)
                result.pop("parameters", None)
            return result

        frame = frame.assign(growth_increment=frame.K_plateau_proxy - frame.K_init_proxy)
        high_growth = float(np.median(_values_at(frame, self.schedule.high_anchor_temperatures_K, "growth_increment")))
        low_growth = float(np.median(_values_at(frame, self.schedule.low_anchor_temperatures_K, "growth_increment")))
        no_collapse_loss = max(-high_growth, 0.0) / 2.0
        low_rcurve_loss = max(low_growth - 3.0, 0.0) / 1.5

        objective = float(
            15.0 * float(init_transition["loss"])
            + 15.0 * float(plateau_transition["loss"])
            + 5.0 * no_collapse_loss**2
            + 2.0 * low_rcurve_loss**2
        )
        result: dict[str, Any] = {
            "objective": objective,
            "completion_loss": 0.0,
            "split_mismatch": 0,
            "high_growth_increment": high_growth,
            "low_growth_increment": low_growth,
            "init_transition_in_coarse_bracket": True,
            "plateau_transition_in_coarse_bracket": True,
            "objective_mode": "DYNAMIC_4POINT_DBTT_GROWTH_REFINEMENT",
            **schedule_columns,
            **{f"init_{key}": value for key, value in init_transition.items() if key != "penalties"},
            **{f"plateau_{key}": value for key, value in plateau_transition.items() if key != "penalties"},
        }
        for key, value in init_transition.get("penalties", {}).items():
            result[f"init_penalty_{key}"] = float(value)
        for key, value in plateau_transition.get("penalties", {}).items():
            result[f"plateau_penalty_{key}"] = float(value)
        if details:
            result["parameters"] = parameters
            result["temperature_detail"] = frame.to_dict(orient="records")
            result["event_detail"] = base_result.get("event_detail") or []
        return result

    def __call__(self, x: np.ndarray) -> float:
        return float(self.evaluate(x, details=False)["objective"])


def growth_acceptance(detail: dict[str, Any]) -> tuple[bool, str]:
    checks = [
        (float(detail.get("init_shelf_ratio", 0.0)) >= 2.0, "initiation_ratio_below_two"),
        (float(detail.get("plateau_shelf_ratio", 0.0)) >= 1.8, "plateau_ratio_too_small"),
        (float(detail.get("init_jump_concentration", 0.0)) >= 0.75, "initiation_rise_not_in_selected_100K_bracket"),
        (float(detail.get("plateau_jump_concentration", 0.0)) >= 0.65, "plateau_rise_not_in_selected_100K_bracket"),
        (float(detail.get("init_transition_width_K", np.inf)) <= 100.0, "initiation_T10_T90_width_exceeds_100K"),
        (float(detail.get("plateau_transition_width_K", np.inf)) <= 100.0, "plateau_T10_T90_width_exceeds_100K"),
        (float(detail.get("init_transition_monotonic_fraction", 0.0)) >= 0.90, "initiation_transition_not_monotonic"),
        (float(detail.get("plateau_transition_monotonic_fraction", 0.0)) >= 0.85, "plateau_transition_not_monotonic"),
        (float(detail.get("high_growth_increment", -99.0)) >= 0.0, "high_temperature_branch_collapses"),
        (float(detail.get("low_growth_increment", 99.0)) <= 3.0, "low_temperature_Rcurve_too_large"),
    ]
    for passed, reason in checks:
        if not passed:
            return False, reason
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
            p = detail.pop("parameters", None) or first.decode(x)
            tdetail = detail.pop("temperature_detail", None) or []
            edetail = detail.pop("event_detail", None) or []
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
