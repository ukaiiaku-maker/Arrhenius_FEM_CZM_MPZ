#!/usr/bin/env python3
"""Mechanism-first narrow-DBTT search using the PF-equivalent reduced front.

Stage 1 searches first passage on the complete 300--1100 K grid. The DBTT
location is free: each candidate is scored at every admissible adjacent split,
and the best split is used. The default ``fixed_zero`` cleavage-slope mode
prevents the historical large activation-entropy shortcut. ``narrow`` may be
used only after the fixed-zero feasibility search is assessed.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.optimize import differential_evolution, minimize

import optimize_mpz_v9_10_unified_global as base
import optimize_mpz_v9_10_2_independent_shape_global as v102
from arrhenius_fracture.reduced_campaign_front_v9104 import (
    ReducedFrontSettings,
    TransitionRequirements,
    best_adjacent_transition,
    simulate_reduced_response,
)

PARAMETER_NAMES = tuple(v102.PARAMETER_NAMES)


def parse_floats(text: str) -> list[float]:
    return [float(value) for value in str(text).replace(",", " ").split() if value]


def search_bounds(cleavage_slope_mode: str) -> dict[str, tuple[float, float]]:
    bounds = dict(v102.BOUNDS)
    if cleavage_slope_mode == "fixed_zero":
        bounds["cleave_gT_eV_per_K"] = (0.0, 0.0)
        bounds["cleave_sT_GPa_per_K"] = (0.0, 0.0)
    elif cleavage_slope_mode == "narrow":
        bounds["cleave_gT_eV_per_K"] = (-0.0015, 0.0015)
        bounds["cleave_sT_GPa_per_K"] = (-0.0010, 0.0010)
    else:
        raise ValueError(f"unknown cleavage slope mode: {cleavage_slope_mode}")
    return bounds


def decode(x: np.ndarray) -> dict[str, float]:
    return v102.decode(np.asarray(x, dtype=float))


class NarrowDBTTObjective:
    def __init__(
        self,
        temperatures: np.ndarray,
        settings: ReducedFrontSettings,
        *,
        cleavage_slope_mode: str,
        transition_requirements: TransitionRequirements,
    ) -> None:
        self.temperatures = np.asarray(temperatures, dtype=float)
        self.settings = settings
        self.cleavage_slope_mode = cleavage_slope_mode
        self.requirements = transition_requirements
        self.bounds_dict = search_bounds(cleavage_slope_mode)
        self.bounds = np.asarray([self.bounds_dict[name] for name in PARAMETER_NAMES], dtype=float)

    def evaluate(self, x: np.ndarray, *, details: bool = False) -> dict[str, Any]:
        x = np.asarray(x, dtype=float)
        if x.shape != (len(PARAMETER_NAMES),) or not np.all(np.isfinite(x)):
            return {"objective": 1.0e12, "invalid_parameter_vector": True}
        outside = np.maximum(self.bounds[:, 0] - x, 0.0) + np.maximum(x - self.bounds[:, 1], 0.0)
        if np.any(outside > 0.0):
            return {"objective": 1.0e10 + 1.0e7 * float(np.sum(outside**2))}

        p = decode(x)
        model = v102.build_model(p, self.settings.Tref_K)
        stress_grid = np.linspace(0.0, 30.0e9, 31)
        order_margin = base.core.barrier_order_margin_eV(model, self.temperatures, stress_grid)
        raw_barriers = [
            model.raw_zero_stress_barrier_eV(mechanism, T)
            for mechanism in ("peierls", "taylor")
            for T in self.temperatures
        ]
        min_raw = float(np.min(raw_barriers))
        if order_margin < -1.0e-9 or min_raw <= 0.0:
            return {
                "objective": 1.0e8
                + 1.0e6 * max(-order_margin, 0.0)
                + 1.0e6 * max(-min_raw, 0.0),
                "barrier_order_margin_eV": float(order_margin),
                "min_raw_barrier_eV": min_raw,
            }

        full_rows: list[dict[str, Any]] = []
        off_rows: list[dict[str, Any]] = []
        event_rows: list[dict[str, Any]] = []
        for T in self.temperatures:
            full = simulate_reduced_response(p, float(T), self.settings, mode="full")
            off = simulate_reduced_response(p, float(T), self.settings, mode="plasticity_off")
            full_events = full.pop("events", [])
            off.pop("events", None)
            full_rows.append({"T_K": float(T), **full})
            off_rows.append({"T_K": float(T), **off})
            if details:
                event_rows.extend({"T_K": float(T), "mode": "full", **row} for row in full_events)

        full_df = pd.DataFrame(full_rows).sort_values("T_K")
        off_df = pd.DataFrame(off_rows).sort_values("T_K")
        incomplete = int((~full_df.completed.astype(bool)).sum() + (~off_df.completed.astype(bool)).sum())
        if incomplete:
            return {
                "objective": 1.0e6 + 1.0e5 * incomplete,
                "completion_loss": 1.0e5 * incomplete,
                "parameters": p if details else None,
                "temperature_detail": full_rows if details else None,
            }

        transition = best_adjacent_transition(
            full_df.T_K,
            full_df.K_init_proxy,
            plasticity_off_toughness=off_df.K_init_proxy,
            requirements=self.requirements,
        )
        j = int(transition["split_index"])
        plastic_increment = full_df.K_init_proxy.to_numpy() - off_df.K_init_proxy.to_numpy()
        low_plastic = float(np.median(plastic_increment[: j + 1]))
        high_plastic = float(np.median(plastic_increment[j + 1 :]))
        total_jump = max(float(transition["high_shelf"] - transition["low_shelf"]), 1.0e-12)
        mechanistic_fraction = (high_plastic - low_plastic) / total_jump
        mechanism_loss = max(0.60 - mechanistic_fraction, 0.0) / 0.15

        slope_regularization = 0.0
        if self.cleavage_slope_mode == "narrow":
            slope_regularization = (
                float(p["cleave_gT_eV_per_K"]) / 5.0e-4
            ) ** 2 + (
                float(p["cleave_sT_GPa_per_K"]) / 4.0e-4
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
            "barrier_order_margin_eV": float(order_margin),
            "min_raw_barrier_eV": min_raw,
            "objective_mode": "NARROW_DBTT_FREE_SPLIT_PF_EQUIVALENT_FIRST_PASSAGE",
            "cleavage_slope_mode": self.cleavage_slope_mode,
            **{f"transition_{key}": value for key, value in transition.items() if key != "penalties"},
        }
        for key, value in transition.get("penalties", {}).items():
            result[f"transition_penalty_{key}"] = float(value)
        if details:
            merged = full_df.merge(
                off_df[["T_K", "K_init_proxy"]].rename(
                    columns={"K_init_proxy": "K_init_plasticity_off"}
                ),
                on="T_K",
                how="left",
            )
            merged["K_init_plastic_increment"] = (
                merged.K_init_proxy - merged.K_init_plasticity_off
            )
            result["parameters"] = p
            result["temperature_detail"] = merged.to_dict(orient="records")
            result["event_detail"] = event_rows
        return result

    def __call__(self, x: np.ndarray) -> float:
        return float(self.evaluate(x, details=False)["objective"])


def accepted(detail: dict[str, Any]) -> tuple[bool, str]:
    checks = [
        (float(detail.get("transition_shelf_ratio", 0.0)) >= 2.0, "shelf_ratio_below_two"),
        (float(detail.get("transition_robust_shelf_ratio", 0.0)) >= 1.8, "robust_ratio_too_small"),
        (float(detail.get("transition_jump_concentration", 0.0)) >= 0.75, "transition_too_broad"),
        (float(detail.get("transition_low_span_fraction", 1.0)) <= 0.15, "low_shelf_not_flat"),
        (float(detail.get("transition_high_span_fraction", 1.0)) <= 0.20, "high_shelf_not_flat"),
        (float(detail.get("transition_plasticity_off_ratio", 99.0)) <= 1.25, "cleavage_only_temperature_cheat"),
        (float(detail.get("mechanistic_fraction", -99.0)) >= 0.60, "plastic_mechanism_fraction_too_small"),
    ]
    for passed, reason in checks:
        if not passed:
            return False, reason
    return True, "narrow_DBTT_first_passage_gate_passed"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--temperatures", default="300 400 500 600 700 800 900 1000 1100")
    ap.add_argument("--cleavage-slope-mode", choices=["fixed_zero", "narrow"], default="fixed_zero")
    ap.add_argument("--restarts", type=int, default=4)
    ap.add_argument("--de-maxiter", type=int, default=80)
    ap.add_argument("--de-popsize", type=int, default=8)
    ap.add_argument("--local-maxiter", type=int, default=250)
    ap.add_argument("--seed", type=int, default=9104001)
    ap.add_argument("--Kdot", type=float, default=0.005)
    ap.add_argument("--Kmax", type=float, default=80.0)
    ap.add_argument("--target-extension-um", type=float, default=5.0)
    ap.add_argument("--max-dK-substep", type=float, default=0.05)
    ap.add_argument("--max-K-shield", type=float, default=1.0)
    ap.add_argument("--population-keep", type=int, default=16)
    ap.add_argument("--shortlist-count", type=int, default=24)
    ap.add_argument("--out", type=Path, default=Path("runs/mpz_v9_10_4_narrow_dbtt_first_passage_v1"))
    args = ap.parse_args()

    temperatures = np.asarray(parse_floats(args.temperatures), dtype=float)
    settings = ReducedFrontSettings(
        Kdot_MPa_sqrt_m_s=float(args.Kdot),
        Kmax_MPa_sqrt_m=float(args.Kmax),
        max_dK_substep_MPa_sqrt_m=float(args.max_dK_substep),
        target_extension_um=float(args.target_extension_um),
        max_K_shield_MPa_sqrt_m=float(args.max_K_shield),
    )
    objective = NarrowDBTTObjective(
        temperatures,
        settings,
        cleavage_slope_mode=args.cleavage_slope_mode,
        transition_requirements=TransitionRequirements(),
    )
    out = args.out.resolve()
    checkpoints = out / "checkpoints"
    checkpoints.mkdir(parents=True, exist_ok=True)
    all_rows: list[dict[str, Any]] = []
    all_temperature: list[dict[str, Any]] = []
    all_events: list[dict[str, Any]] = []
    all_history: list[dict[str, Any]] = []
    bounds = [objective.bounds_dict[name] for name in PARAMETER_NAMES]

    for restart in range(args.restarts):
        checkpoint = checkpoints / f"restart_{restart:03d}.json"
        if checkpoint.exists():
            payload = json.loads(checkpoint.read_text())
            if payload.get("status") == "COMPLETE":
                all_rows.extend(payload.get("candidates", []))
                all_temperature.extend(payload.get("temperature_detail", []))
                all_events.extend(payload.get("event_detail", []))
                all_history.extend(payload.get("history", []))
                print(f"resumed completed restart {restart}", flush=True)
                continue

        history: list[dict[str, Any]] = []

        def callback(xk: np.ndarray, convergence: float) -> bool:
            rec = {
                "restart": restart,
                "generation": len(history),
                "objective": objective(xk),
                "convergence": float(convergence),
            }
            history.append(rec)
            print(
                f"restart={restart} generation={len(history)} "
                f"objective={rec['objective']:.6g}",
                flush=True,
            )
            return False

        de = differential_evolution(
            objective,
            bounds,
            maxiter=args.de_maxiter,
            popsize=args.de_popsize,
            seed=args.seed + 1009 * restart,
            init="sobol",
            polish=False,
            updating="immediate",
            workers=1,
            callback=callback,
            tol=1.0e-3,
        )
        local = minimize(
            objective,
            de.x,
            method="Powell",
            bounds=bounds,
            options={"maxiter": args.local_maxiter, "xtol": 1.0e-4, "ftol": 1.0e-5},
        )
        candidates: list[tuple[float, np.ndarray, str]] = [
            (float(de.fun), np.asarray(de.x, dtype=float), "de_best")
        ]
        if np.isfinite(local.fun):
            candidates.append((float(local.fun), np.asarray(local.x, dtype=float), "local_best"))
        order = np.argsort(np.asarray(de.population_energies, dtype=float))[: args.population_keep]
        for rank in order:
            candidates.append(
                (
                    float(de.population_energies[rank]),
                    np.asarray(de.population[rank], dtype=float),
                    f"population_{rank}",
                )
            )

        restart_rows: list[dict[str, Any]] = []
        restart_temperature: list[dict[str, Any]] = []
        restart_events: list[dict[str, Any]] = []
        seen: set[tuple[float, ...]] = set()
        for rank, (_, x, source) in enumerate(sorted(candidates, key=lambda row: row[0])):
            key = tuple(np.round(x, 10))
            if key in seen:
                continue
            seen.add(key)
            detail = objective.evaluate(x, details=True)
            p = detail.pop("parameters")
            temperature_detail = detail.pop("temperature_detail")
            event_detail = detail.pop("event_detail")
            candidate_id = f"DBTT_v9104_restart{restart:02d}_candidate{rank:02d}"
            is_accepted, reason = accepted(detail)
            row = {
                "candidate_id": candidate_id,
                "target_class": "DBTT",
                "restart": restart,
                "candidate_source": source,
                "accepted_for_short_growth": bool(is_accepted),
                "acceptance_reason": reason,
                "search_initialization": "FULL_SOBOL_NO_PRIOR_SHORTLIST",
                "reduced_model": "PF_EQUIVALENT_STRESS_SEPARATED_MOVING_TIP_V9104",
                **{name: float(x[i]) for i, name in enumerate(PARAMETER_NAMES)},
                **{key: float(value) for key, value in p.items() if key not in PARAMETER_NAMES},
                **detail,
            }
            restart_rows.append(row)
            restart_temperature.extend(
                {"candidate_id": candidate_id, **record} for record in temperature_detail
            )
            restart_events.extend(
                {"candidate_id": candidate_id, **record} for record in event_detail
            )

        payload = {
            "status": "COMPLETE",
            "restart": restart,
            "candidates": restart_rows,
            "temperature_detail": restart_temperature,
            "event_detail": restart_events,
            "history": history,
        }
        checkpoint.write_text(json.dumps(payload, indent=2, allow_nan=True))
        all_rows.extend(restart_rows)
        all_temperature.extend(restart_temperature)
        all_events.extend(restart_events)
        all_history.extend(history)

    results = pd.DataFrame(all_rows).sort_values("objective").drop_duplicates(PARAMETER_NAMES)
    accepted_df = results[results.accepted_for_short_growth.astype(bool)].copy()
    shortlist = (accepted_df if not accepted_df.empty else results).head(args.shortlist_count).copy()
    results.to_csv(out / "narrow_dbtt_first_passage_all_candidates.csv", index=False)
    accepted_df.to_csv(out / "narrow_dbtt_first_passage_accepted.csv", index=False)
    shortlist.to_csv(out / "narrow_dbtt_first_passage_shortlist.csv", index=False)
    shortlist.to_csv(out / "short_growth_promotion_manifest.csv", index=False)
    pd.DataFrame(all_temperature).to_csv(out / "narrow_dbtt_first_passage_temperature_detail.csv", index=False)
    pd.DataFrame(all_events).to_csv(out / "narrow_dbtt_first_passage_event_detail.csv", index=False)
    pd.DataFrame(all_history).to_csv(out / "narrow_dbtt_first_passage_generation_history.csv", index=False)
    summary = {
        "status": "V9_10_4_NARROW_DBTT_FIRST_PASSAGE_COMPLETE",
        "n_candidates": int(len(results)),
        "n_accepted": int(len(accepted_df)),
        "best_objective": float(results.iloc[0].objective),
        "cleavage_slope_mode": args.cleavage_slope_mode,
        "temperatures_K": temperatures.tolist(),
        "target_extension_um": float(args.target_extension_um),
        "next_stage_manifest": str(out / "short_growth_promotion_manifest.csv"),
    }
    (out / "narrow_dbtt_first_passage_summary.json").write_text(json.dumps(summary, indent=2))
    config = vars(args).copy()
    config["out"] = str(config["out"])
    config.update({"parameter_names": PARAMETER_NAMES, "bounds": objective.bounds_dict})
    (out / "narrow_dbtt_first_passage_config.json").write_text(json.dumps(config, indent=2))
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
