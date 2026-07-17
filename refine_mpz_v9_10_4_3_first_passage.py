#!/usr/bin/env python3
"""Candidate-specific four-point first-passage refinement for narrow DBTT.

Each input candidate was first evaluated on the complete 100 K coarse grid.
The selected coarse bracket is held fixed. Four points resolve the full bracket,
while broad shelf anchors retain the factor-of-two and shelf-flatness tests.
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
    fixed_bracket_transition_metrics,
    schedule_from_candidate_row,
)
from arrhenius_fracture.reduced_campaign_front_v9104 import (
    ReducedFrontSettings,
    TransitionRequirements,
)

PARAMETER_NAMES = first.PARAMETER_NAMES


def _values_at(frame: pd.DataFrame, temperatures: tuple[float, ...], column: str) -> np.ndarray:
    values = []
    for temperature in temperatures:
        matches = frame.loc[np.isclose(frame.T_K, temperature, rtol=0.0, atol=1.0e-7), column]
        if len(matches) != 1:
            raise ValueError(f"expected one {column} value at {temperature} K; found {len(matches)}")
        values.append(float(matches.iloc[0]))
    return np.asarray(values, dtype=float)


class WindowedFirstPassageObjective:
    def __init__(
        self,
        schedule: DynamicTemperatureSchedule,
        settings: ReducedFrontSettings,
        *,
        cleavage_slope_mode: str,
    ) -> None:
        self.schedule = schedule
        self.cleavage_slope_mode = cleavage_slope_mode
        self.base = first.NarrowDBTTObjective(
            np.asarray(schedule.evaluation_temperatures_K, dtype=float),
            settings,
            cleavage_slope_mode=cleavage_slope_mode,
            transition_requirements=TransitionRequirements(),
        )
        self.bounds = [self.base.bounds_dict[name] for name in PARAMETER_NAMES]

    def evaluate(self, x: np.ndarray, *, details: bool = False) -> dict[str, Any]:
        x = np.asarray(x, dtype=float)
        base_result = dict(self.base.evaluate(x, details=True))
        schedule_columns = self.schedule.to_columns()
        temperature_detail = base_result.get("temperature_detail") or []
        parameters = base_result.get("parameters")

        if not temperature_detail or parameters is None:
            result = dict(base_result)
            result.update(schedule_columns)
            result["coarse_transition_low_T_K"] = self.schedule.transition_low_K
            result["coarse_transition_high_T_K"] = self.schedule.transition_high_K
            result["refined_transition_in_coarse_bracket"] = False
            result["objective_mode"] = "DYNAMIC_4POINT_DBTT_FIRST_PASSAGE_REFINEMENT"
            if not details:
                result.pop("temperature_detail", None)
                result.pop("event_detail", None)
                result.pop("parameters", None)
            return result

        frame = pd.DataFrame(temperature_detail).sort_values("T_K")
        required = {"T_K", "K_init_proxy", "K_init_plasticity_off", "completed"}
        if not required.issubset(frame.columns) or not frame.completed.astype(bool).all():
            result = {
                "objective": 1.0e6 + 1.0e5 * int((~frame.get("completed", pd.Series(False, index=frame.index)).astype(bool)).sum()),
                "completion_loss": 1.0e5,
                "parameters": parameters,
                "temperature_detail": temperature_detail,
                "event_detail": base_result.get("event_detail") or [],
                "refined_transition_in_coarse_bracket": False,
                "objective_mode": "DYNAMIC_4POINT_DBTT_FIRST_PASSAGE_REFINEMENT",
                **schedule_columns,
            }
            if not details:
                result.pop("temperature_detail", None)
                result.pop("event_detail", None)
                result.pop("parameters", None)
            return result

        transition = fixed_bracket_transition_metrics(
            frame.T_K,
            frame.K_init_proxy,
            self.schedule,
            plasticity_off_toughness=frame.K_init_plasticity_off,
        )
        if not bool(transition.get("valid", False)):
            result = {
                "objective": 1.0e12,
                "transition_valid": False,
                "transition_reason": transition.get("reason", "invalid_transition"),
                "parameters": parameters,
                "temperature_detail": temperature_detail,
                "event_detail": base_result.get("event_detail") or [],
                "refined_transition_in_coarse_bracket": False,
                "objective_mode": "DYNAMIC_4POINT_DBTT_FIRST_PASSAGE_REFINEMENT",
                **schedule_columns,
            }
            if not details:
                result.pop("temperature_detail", None)
                result.pop("event_detail", None)
                result.pop("parameters", None)
            return result

        plastic_increment = frame.K_init_proxy.to_numpy() - frame.K_init_plasticity_off.to_numpy()
        frame = frame.assign(K_init_plastic_increment=plastic_increment)
        low_plastic = float(np.median(_values_at(frame, self.schedule.low_anchor_temperatures_K, "K_init_plastic_increment")))
        high_plastic = float(np.median(_values_at(frame, self.schedule.high_anchor_temperatures_K, "K_init_plastic_increment")))
        total_jump = max(float(transition["high_shelf"] - transition["low_shelf"]), 1.0e-12)
        mechanistic_fraction = (high_plastic - low_plastic) / total_jump
        mechanism_loss = max(0.60 - mechanistic_fraction, 0.0) / 0.15

        slope_regularization = 0.0
        if self.cleavage_slope_mode == "narrow":
            slope_regularization = (
                float(parameters["cleave_gT_eV_per_K"]) / 5.0e-4
            ) ** 2 + (
                float(parameters["cleave_sT_GPa_per_K"]) / 4.0e-4
            ) ** 2

        objective = float(
            20.0 * float(transition["loss"])
            + 10.0 * mechanism_loss**2
            + 2.0 * slope_regularization
        )
        result: dict[str, Any] = {
            "objective": objective,
            "completion_loss": 0.0,
            "transition_loss": float(transition["loss"]),
            "mechanism_loss": float(mechanism_loss**2),
            "cleavage_slope_regularization": float(slope_regularization),
            "mechanistic_fraction": float(mechanistic_fraction),
            "low_plastic_increment": low_plastic,
            "high_plastic_increment": high_plastic,
            "barrier_order_margin_eV": float(base_result.get("barrier_order_margin_eV", np.nan)),
            "min_raw_barrier_eV": float(base_result.get("min_raw_barrier_eV", np.nan)),
            "objective_mode": "DYNAMIC_4POINT_DBTT_FIRST_PASSAGE_REFINEMENT",
            "cleavage_slope_mode": self.cleavage_slope_mode,
            "refined_transition_in_coarse_bracket": True,
            **schedule_columns,
            **{f"transition_{key}": value for key, value in transition.items() if key != "penalties"},
        }
        for key, value in transition.get("penalties", {}).items():
            result[f"transition_penalty_{key}"] = float(value)
        if details:
            result["parameters"] = parameters
            result["temperature_detail"] = frame.to_dict(orient="records")
            result["event_detail"] = base_result.get("event_detail") or []
        return result

    def __call__(self, x: np.ndarray) -> float:
        return float(self.evaluate(x, details=False)["objective"])


def accepted(detail: dict[str, Any]) -> tuple[bool, str]:
    checks = [
        (float(detail.get("transition_shelf_ratio", 0.0)) >= 2.0, "shelf_ratio_below_two"),
        (float(detail.get("transition_robust_shelf_ratio", 0.0)) >= 1.8, "robust_ratio_too_small"),
        (float(detail.get("transition_jump_concentration", 0.0)) >= 0.75, "rise_not_concentrated_in_selected_100K_bracket"),
        (float(detail.get("transition_transition_width_K", np.inf)) <= 100.0, "T10_T90_width_exceeds_100K"),
        (float(detail.get("transition_transition_monotonic_fraction", 0.0)) >= 0.90, "transition_not_monotonic"),
        (float(detail.get("transition_low_span_fraction", 1.0)) <= 0.15, "low_shelf_not_flat"),
        (float(detail.get("transition_high_span_fraction", 1.0)) <= 0.20, "high_shelf_not_flat"),
        (float(detail.get("transition_plasticity_off_ratio", 99.0)) <= 1.25, "cleavage_only_temperature_cheat"),
        (float(detail.get("mechanistic_fraction", -99.0)) >= 0.60, "plastic_mechanism_fraction_too_small"),
    ]
    for passed, reason in checks:
        if not passed:
            return False, reason
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
