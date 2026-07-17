#!/usr/bin/env python3
"""Short- and long-growth refinement for the v9.10.4 narrow-DBTT campaign."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.optimize import minimize

import optimize_mpz_v9_10_4_narrow_dbtt as first
from arrhenius_fracture.reduced_campaign_front_v9104 import (
    ReducedFrontSettings,
    TransitionRequirements,
    best_adjacent_transition,
    simulate_reduced_response,
)

PARAMETER_NAMES = first.PARAMETER_NAMES


class GrowthObjective:
    def __init__(
        self,
        temperatures: np.ndarray,
        settings: ReducedFrontSettings,
        *,
        cleavage_slope_mode: str,
    ) -> None:
        self.temperatures = np.asarray(temperatures, dtype=float)
        self.settings = settings
        self.cleavage_slope_mode = cleavage_slope_mode
        self.bounds_dict = first.search_bounds(cleavage_slope_mode)
        self.bounds = [self.bounds_dict[name] for name in PARAMETER_NAMES]
        self.requirements = TransitionRequirements()

    def evaluate(self, x: np.ndarray, *, details: bool = False) -> dict[str, Any]:
        x = np.asarray(x, dtype=float)
        p = first.decode(x)
        full_rows: list[dict[str, Any]] = []
        off_rows: list[dict[str, Any]] = []
        event_rows: list[dict[str, Any]] = []
        for T in self.temperatures:
            full = simulate_reduced_response(p, float(T), self.settings, mode="full")
            off = simulate_reduced_response(p, float(T), self.settings, mode="plasticity_off")
            events = full.pop("events", [])
            off.pop("events", None)
            full_rows.append({"T_K": float(T), **full})
            off_rows.append({"T_K": float(T), **off})
            if details:
                event_rows.extend({"T_K": float(T), "mode": "full", **row} for row in events)
        full_df = pd.DataFrame(full_rows).sort_values("T_K")
        off_df = pd.DataFrame(off_rows).sort_values("T_K")
        incomplete = int((~full_df.completed.astype(bool)).sum() + (~off_df.completed.astype(bool)).sum())
        if incomplete:
            return {"objective": 1.0e6 + 1.0e5 * incomplete, "completion_loss": 1.0e5 * incomplete}

        init_transition = best_adjacent_transition(
            full_df.T_K,
            full_df.K_init_proxy,
            plasticity_off_toughness=off_df.K_init_proxy,
            requirements=self.requirements,
        )
        plateau_transition = best_adjacent_transition(
            full_df.T_K,
            full_df.K_plateau_proxy,
            plasticity_off_toughness=off_df.K_plateau_proxy,
            requirements=self.requirements,
        )
        split_mismatch = abs(
            int(init_transition["split_index"]) - int(plateau_transition["split_index"])
        )
        j = int(init_transition["split_index"])
        high = full_df.iloc[j + 1 :]
        low = full_df.iloc[: j + 1]
        high_growth = float(np.median(high.K_plateau_proxy - high.K_init_proxy))
        low_growth = float(np.median(low.K_plateau_proxy - low.K_init_proxy))
        no_collapse_loss = max(-high_growth, 0.0) / 2.0
        low_rcurve_loss = max(low_growth - 3.0, 0.0) / 1.5

        objective = float(
            15.0 * init_transition["loss"]
            + 15.0 * plateau_transition["loss"]
            + 5.0 * split_mismatch**2
            + 5.0 * no_collapse_loss**2
            + 2.0 * low_rcurve_loss**2
        )
        result: dict[str, Any] = {
            "objective": objective,
            "completion_loss": 0.0,
            "split_mismatch": int(split_mismatch),
            "high_growth_increment": high_growth,
            "low_growth_increment": low_growth,
            **{f"init_{key}": value for key, value in init_transition.items() if key != "penalties"},
            **{f"plateau_{key}": value for key, value in plateau_transition.items() if key != "penalties"},
        }
        if details:
            merged = full_df.merge(
                off_df[["T_K", "K_init_proxy", "K_plateau_proxy"]].rename(
                    columns={
                        "K_init_proxy": "K_init_plasticity_off",
                        "K_plateau_proxy": "K_plateau_plasticity_off",
                    }
                ),
                on="T_K",
                how="left",
            )
            result["parameters"] = p
            result["temperature_detail"] = merged.to_dict(orient="records")
            result["event_detail"] = event_rows
        return result

    def __call__(self, x: np.ndarray) -> float:
        return float(self.evaluate(x, details=False)["objective"])


def growth_acceptance(detail: dict[str, Any]) -> tuple[bool, str]:
    checks = [
        (float(detail.get("init_shelf_ratio", 0.0)) >= 2.0, "initiation_ratio_below_two"),
        (float(detail.get("plateau_shelf_ratio", 0.0)) >= 1.8, "plateau_ratio_too_small"),
        (float(detail.get("init_jump_concentration", 0.0)) >= 0.70, "initiation_transition_too_broad"),
        (float(detail.get("plateau_jump_concentration", 0.0)) >= 0.65, "plateau_transition_too_broad"),
        (int(detail.get("split_mismatch", 99)) == 0, "initiation_and_plateau_splits_differ"),
        (float(detail.get("high_growth_increment", -99.0)) >= 0.0, "high_temperature_branch_collapses"),
        (float(detail.get("low_growth_increment", 99.0)) <= 3.0, "low_temperature_Rcurve_too_large"),
    ]
    for passed, reason in checks:
        if not passed:
            return False, reason
    return True, "narrow_DBTT_growth_gate_passed"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-manifest", type=Path, required=True)
    ap.add_argument("--stage", choices=["short", "long"], required=True)
    ap.add_argument("--temperatures", default="300 400 500 600 700 800 900 1000 1100")
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
    temperatures = np.asarray(first.parse_floats(args.temperatures), dtype=float)
    settings = ReducedFrontSettings(
        Kdot_MPa_sqrt_m_s=float(args.Kdot),
        Kmax_MPa_sqrt_m=float(args.Kmax),
        max_dK_substep_MPa_sqrt_m=float(args.max_dK_substep),
        target_extension_um=extension,
        max_K_shield_MPa_sqrt_m=float(args.max_K_shield),
    )
    objective = GrowthObjective(
        temperatures,
        settings,
        cleavage_slope_mode=args.cleavage_slope_mode,
    )
    manifest = pd.read_csv(args.input_manifest).head(args.max_candidates)
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    temperature_rows: list[dict[str, Any]] = []
    event_rows: list[dict[str, Any]] = []

    for rank, row in manifest.iterrows():
        x0 = np.asarray([float(row[name]) for name in PARAMETER_NAMES], dtype=float)
        local = minimize(
            objective,
            x0,
            method="Powell",
            bounds=objective.bounds,
            options={"maxiter": args.local_maxiter, "xtol": 1.0e-4, "ftol": 1.0e-5},
        )
        candidates = [(objective(x0), x0, "input"), (float(local.fun), np.asarray(local.x), "local_refined")]
        for source_index, (_, x, source) in enumerate(sorted(candidates, key=lambda item: item[0])):
            detail = objective.evaluate(x, details=True)
            p = detail.pop("parameters")
            tdetail = detail.pop("temperature_detail")
            edetail = detail.pop("event_detail")
            candidate_id = f"DBTT_v9104_{args.stage}_{rank:02d}_{source_index:02d}"
            accepted, reason = growth_acceptance(detail)
            output_row = {
                "candidate_id": candidate_id,
                "parent_candidate_id": str(row.get("candidate_id", rank)),
                "candidate_source": source,
                "stage": args.stage,
                "target_extension_um": extension,
                "accepted_for_next_stage": bool(accepted),
                "acceptance_reason": reason,
                **{name: float(x[i]) for i, name in enumerate(PARAMETER_NAMES)},
                **{key: float(value) for key, value in p.items() if key not in PARAMETER_NAMES},
                **detail,
            }
            rows.append(output_row)
            temperature_rows.extend({"candidate_id": candidate_id, **record} for record in tdetail)
            event_rows.extend({"candidate_id": candidate_id, **record} for record in edetail)
            print(
                f"stage={args.stage} candidate={candidate_id} objective={output_row['objective']:.6g} "
                f"accepted={accepted} reason={reason}",
                flush=True,
            )

    results = pd.DataFrame(rows).sort_values("objective").drop_duplicates(PARAMETER_NAMES)
    accepted_df = results[results.accepted_for_next_stage.astype(bool)].copy()
    promotion = (accepted_df if not accepted_df.empty else results).head(8).copy()
    prefix = f"narrow_dbtt_{args.stage}_growth"
    results.to_csv(out / f"{prefix}_all_candidates.csv", index=False)
    accepted_df.to_csv(out / f"{prefix}_accepted.csv", index=False)
    promotion.to_csv(out / f"{prefix}_promotion_manifest.csv", index=False)
    pd.DataFrame(temperature_rows).to_csv(out / f"{prefix}_temperature_detail.csv", index=False)
    pd.DataFrame(event_rows).to_csv(out / f"{prefix}_event_detail.csv", index=False)
    summary = {
        "status": f"V9_10_4_NARROW_DBTT_{args.stage.upper()}_GROWTH_COMPLETE",
        "target_extension_um": extension,
        "n_candidates": int(len(results)),
        "n_accepted": int(len(accepted_df)),
        "promotion_manifest": str(out / f"{prefix}_promotion_manifest.csv"),
    }
    (out / f"{prefix}_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
