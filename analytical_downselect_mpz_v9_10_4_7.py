#!/usr/bin/env python3
"""Analytical/zero-D DBTT down-selection before moving-interface evaluation.

This stage samples the complete v9.10.2 parameter domain with a deterministic
Sobol design, evaluates a cheap first-event zero-D closure on the broad
300--1100 K grid, assigns every viable candidate to its best adjacent 100 K
transition bracket, and retains the best candidates within each bracket.

The zero-D closure preserves the mechanism ordering used by the moving front:
finite source depletion, emission back stress, Peierls transport, Taylor
retention/release, recovery, escape, blunting, shielding of cleavage only, and
an Arrhenius cleavage clock.  It does not translate the crack-tip process zone
or refresh sources because it stops at the first cleavage event.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, asdict
import json
import math
from pathlib import Path
import time
from typing import Any, Iterable

import numpy as np
import pandas as pd
from scipy.stats import qmc

import optimize_mpz_v9_10_unified_global as zero
import optimize_mpz_v9_10_2_independent_shape_global as v102
from arrhenius_fracture.config import EV_TO_J, KB
from arrhenius_fracture.moving_process_zone_v910 import MovingProcessZoneState
from arrhenius_fracture.reduced_campaign_front_v9104 import (
    TransitionRequirements,
    best_adjacent_transition,
)

PARAMETER_NAMES = tuple(v102.PARAMETER_NAMES)


@dataclass(frozen=True)
class AnalyticalSettings:
    temperatures_K: tuple[float, ...] = tuple(float(x) for x in range(300, 1101, 100))
    Kdot_MPa_sqrt_m_s: float = 0.005
    dK_MPa_sqrt_m: float = 0.5
    Kmax_MPa_sqrt_m: float = 80.0
    r0_m: float = 1.0e-6
    L_pz_m: float = 50.0e-6
    b_m: float = 2.74e-10
    G_Pa: float = 160.0e9
    poisson: float = 0.28
    rho0_m2: float = 5.0e12
    nu0_c_s: float = 1.0e12
    nu0_e_s: float = 1.0e11
    cleavage_hits: float = 3.0
    cleavage_tau_s: float = 1.0e-6
    Tref_K: float = 481.33
    mobile_shield_fraction: float = 0.0
    resolved_stress_fraction: float = 1.0
    max_K_shield_MPa_sqrt_m: float = 1.0


def parse_floats(text: str) -> tuple[float, ...]:
    return tuple(float(x) for x in str(text).replace(",", " ").split() if x)


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


def _cleavage_rate_s(G_eV: float, T_K: float, s: AnalyticalSettings) -> float:
    return float(
        zero.core.cleavage_effective_rate(
            np.asarray([G_eV]),
            np.asarray([T_K]),
            s.nu0_c_s,
            s.cleavage_hits,
            s.cleavage_tau_s,
        )[0]
    )


def _shield_coefficient(s: AnalyticalSettings) -> float:
    return (
        s.G_Pa
        * s.b_m
        / max(1.0 - s.poisson, 1.0e-12)
        / math.sqrt(2.0 * math.pi * max(0.5 * s.r0_m, s.b_m))
        / 1.0e6
    )


def analytical_first_passage(
    p: dict[str, float],
    T_K: float,
    s: AnalyticalSettings,
    *,
    plasticity_active: bool,
) -> dict[str, Any]:
    """Return the first cleavage event from the zero-D kinetic closure."""
    model = v102.build_model(p, s.Tref_K)
    capacity = max(2.0 * float(p["source_sites_per_system"]), 0.0)
    available = capacity
    mobile = 0.0
    retained = 0.0
    slip = 0.0
    B = 0.0
    shield_coefficient = _shield_coefficient(s)
    area = math.pi * max(s.r0_m, s.b_m) ** 2
    last: dict[str, float] = {}

    n_steps = int(math.ceil(s.Kmax_MPa_sqrt_m / s.dK_MPa_sqrt_m))
    for step in range(1, n_steps + 1):
        K = min(step * s.dK_MPa_sqrt_m, s.Kmax_MPa_sqrt_m)
        dt = s.dK_MPa_sqrt_m / max(s.Kdot_MPa_sqrt_m_s, 1.0e-30)
        r_eff = s.r0_m + max(float(p.get("c_blunt", 0.0)), 0.0) * s.b_m * max(slip, 0.0)
        denominator = math.sqrt(2.0 * math.pi * max(r_eff, 1.0e-30))
        sigma_open = K * 1.0e6 / denominator

        p_rate = 0.0
        t_rate = 0.0
        encounter = 0.0
        velocity = 0.0
        sigma_back = 0.0
        if plasticity_active:
            rho_active = max(s.rho0_m2 + (mobile + retained) / max(area, 1.0e-30), 1.0)
            sigma_back = (
                s.G_Pa
                * s.b_m
                * math.sqrt(rho_active)
                / max(s.resolved_stress_fraction, 1.0e-12)
            )
            sigma_emit = max(sigma_open - sigma_back, 0.0)
            Ge = zero.barrier_eV(p, "emit", sigma_emit, T_K, s.Tref_K)
            emit_rate = s.nu0_e_s * math.exp(
                float(np.clip(-Ge * EV_TO_J / (KB * T_K), -700.0, 0.0))
            )
            emitted = available * (1.0 - math.exp(-min(emit_rate * dt, 700.0)))
            available = max(available - emitted, 0.0)
            mobile += emitted
            slip += emitted

            rho_active = max(s.rho0_m2 + (mobile + retained) / max(area, 1.0e-30), 1.0)
            rates = model.rates(sigma_open, rho_active, T_K, s.b_m)
            p_rate = float(np.asarray(rates["peierls_rate_s"]))
            t_rate = float(np.asarray(rates["taylor_completion_rate_s"]))
            jump = float(np.asarray(rates["jump_length_m"]))
            velocity = max(jump * p_rate, 0.0)
            encounter = float(
                MovingProcessZoneState.encounter_rate_s(
                    p_rate, jump, rho_active, p["encounter_efficiency"]
                )
            )
            total = max(mobile, 0.0) + max(retained, 0.0)
            exchange = max(encounter + t_rate, 0.0)
            if exchange > 0.0 and total > 0.0:
                retained_eq = encounter / exchange * total
                retained = retained_eq + (retained - retained_eq) * math.exp(
                    -min(exchange * dt, 700.0)
                )
                retained = min(max(retained, 0.0), total)
                mobile = total - retained
            recovery_rate = max(float(p["retained_recovery_rate_s"]), 0.0)
            retained *= math.exp(-min(recovery_rate * dt, 700.0))
            escape_rate = velocity / max(s.L_pz_m, 1.0e-30)
            mobile *= math.exp(-min(escape_rate * dt, 700.0))

        K_shield = 0.0
        if plasticity_active:
            K_shield = shield_coefficient * (
                max(retained, 0.0)
                + max(s.mobile_shield_fraction, 0.0) * max(mobile, 0.0)
            )
            K_shield = float(np.clip(K_shield, 0.0, s.max_K_shield_MPa_sqrt_m))
        sigma_cleave = max(K - K_shield, 0.0) * 1.0e6 / denominator
        Gc = zero.barrier_eV(p, "cleave", sigma_cleave, T_K, s.Tref_K)
        lambda_c = _cleavage_rate_s(Gc, T_K, s)
        B += lambda_c * dt
        last = {
            "K_MPa_sqrt_m": K,
            "B": B,
            "available_sites": available,
            "mobile_count": mobile,
            "retained_count": retained,
            "slip_count": slip,
            "K_shield_MPa_sqrt_m": K_shield,
            "sigma_back_Pa": sigma_back,
            "peierls_rate_s": p_rate,
            "taylor_completion_rate_s": t_rate,
            "encounter_rate_s": encounter,
            "velocity_m_s": velocity,
            "lambda_c_s": lambda_c,
        }
        if B >= 1.0:
            return {"completed": True, "internal_K_steps": step, **last}

    return {
        "completed": False,
        "internal_K_steps": n_steps,
        "K_MPa_sqrt_m": float("nan"),
        **last,
    }


def _screen_pass(transition: dict[str, Any], mechanistic_fraction: float) -> tuple[bool, str]:
    checks = [
        (float(transition.get("shelf_ratio", 0.0)) >= 1.50, "analytical_ratio_too_small"),
        (float(transition.get("robust_shelf_ratio", 0.0)) >= 1.30, "analytical_robust_ratio_too_small"),
        (float(transition.get("jump_concentration", 0.0)) >= 0.35, "analytical_transition_too_broad"),
        (float(transition.get("low_span_fraction", 99.0)) <= 0.35, "analytical_low_shelf_not_flat"),
        (float(transition.get("high_span_fraction", 99.0)) <= 0.35, "analytical_high_shelf_not_flat"),
        (float(transition.get("plasticity_off_ratio", 99.0)) <= 1.50, "analytical_cleavage_only_T_dependence"),
        (mechanistic_fraction >= 0.20, "analytical_plastic_increment_too_small"),
    ]
    for passed, reason in checks:
        if not passed:
            return False, reason
    return True, "analytical_screen_passed"


def evaluate_candidate(payload: tuple[int, np.ndarray, dict[str, Any]]) -> dict[str, Any]:
    candidate_index, x, config = payload
    s = AnalyticalSettings(**config["settings"])
    temperatures = np.asarray(s.temperatures_K, dtype=float)
    p = v102.decode(np.asarray(x, dtype=float))
    candidate_id = f"DBTT_A{candidate_index:07d}"

    model = v102.build_model(p, s.Tref_K)
    stress_grid = np.linspace(0.0, 30.0e9, 31)
    order_margin = zero.core.barrier_order_margin_eV(model, temperatures, stress_grid)
    raw_barriers = [
        model.raw_zero_stress_barrier_eV(mechanism, T)
        for mechanism in ("peierls", "taylor")
        for T in temperatures
    ]
    min_raw = float(np.min(raw_barriers))
    base_row: dict[str, Any] = {
        "candidate_id": candidate_id,
        "candidate_index": int(candidate_index),
        "barrier_order_margin_eV": float(order_margin),
        "min_raw_barrier_eV": min_raw,
        **{name: float(x[i]) for i, name in enumerate(PARAMETER_NAMES)},
    }
    if order_margin < -1.0e-9 or min_raw <= 0.0:
        return {
            **base_row,
            "analysis_valid": False,
            "screen_pass": False,
            "screen_reason": "barrier_hierarchy_invalid",
            "objective": 1.0e8 + 1.0e6 * max(-order_margin, 0.0) + 1.0e6 * max(-min_raw, 0.0),
        }

    full_runs = [analytical_first_passage(p, float(T), s, plasticity_active=True) for T in temperatures]
    off_runs = [analytical_first_passage(p, float(T), s, plasticity_active=False) for T in temperatures]
    full_K = np.asarray([r.get("K_MPa_sqrt_m", np.nan) for r in full_runs], dtype=float)
    off_K = np.asarray([r.get("K_MPa_sqrt_m", np.nan) for r in off_runs], dtype=float)
    complete = bool(np.all(np.isfinite(full_K)) and np.all(np.isfinite(off_K)))
    if not complete:
        return {
            **base_row,
            "analysis_valid": False,
            "screen_pass": False,
            "screen_reason": "incomplete_first_passage",
            "objective": 1.0e6 + 1.0e5 * int(np.sum(~np.isfinite(full_K)) + np.sum(~np.isfinite(off_K))),
            "analytical_full_K_json": json.dumps(full_K.tolist()),
            "analytical_off_K_json": json.dumps(off_K.tolist()),
        }

    transition = best_adjacent_transition(
        temperatures,
        full_K,
        plasticity_off_toughness=off_K,
        requirements=TransitionRequirements(),
    )
    if not bool(transition.get("valid", True)):
        return {
            **base_row,
            "analysis_valid": False,
            "screen_pass": False,
            "screen_reason": str(transition.get("reason", "invalid_transition")),
            "objective": float(transition.get("loss", 1.0e12)),
            "analytical_full_K_json": json.dumps(full_K.tolist()),
            "analytical_off_K_json": json.dumps(off_K.tolist()),
        }

    split = int(transition["split_index"])
    plastic_increment = full_K - off_K
    low_plastic = float(np.median(plastic_increment[: split + 1]))
    high_plastic = float(np.median(plastic_increment[split + 1 :]))
    shelf_jump = max(float(transition["high_shelf"] - transition["low_shelf"]), 1.0e-12)
    mechanistic_fraction = (high_plastic - low_plastic) / shelf_jump
    mechanism_penalty = max(0.60 - mechanistic_fraction, 0.0) / 0.15
    objective = 20.0 * float(transition["loss"]) + 10.0 * mechanism_penalty**2
    passed, reason = _screen_pass(transition, mechanistic_fraction)
    low_T = float(transition["transition_low_K"])
    high_T = float(transition["transition_high_K"])
    refinement = np.linspace(low_T, high_T, 4)
    bracket_label = f"T{int(round(low_T)):04d}_{int(round(high_T)):04d}K"

    row = {
        **base_row,
        "analysis_valid": True,
        "screen_pass": bool(passed),
        "screen_reason": reason,
        "objective": float(objective),
        "mechanistic_fraction": float(mechanistic_fraction),
        "low_plastic_increment": low_plastic,
        "high_plastic_increment": high_plastic,
        "coarse_transition_low_T_K": low_T,
        "coarse_transition_high_T_K": high_T,
        "transition_bracket": bracket_label,
        "refinement_transition_temperatures_K": json.dumps([float(v) for v in refinement]),
        "analytical_full_K_json": json.dumps(full_K.tolist()),
        "analytical_off_K_json": json.dumps(off_K.tolist()),
        **{f"transition_{key}": value for key, value in transition.items() if key != "penalties"},
    }
    for key, value in transition.get("penalties", {}).items():
        row[f"transition_penalty_{key}"] = float(value)
    for name, value in p.items():
        if name not in row and np.isscalar(value):
            row[name] = float(value)
    return row


def sobol_parameter_vectors(
    n_samples: int,
    bounds: dict[str, tuple[float, float]],
    seed: int,
) -> np.ndarray:
    free = [i for i, name in enumerate(PARAMETER_NAMES) if bounds[name][1] > bounds[name][0]]
    fixed = [i for i, name in enumerate(PARAMETER_NAMES) if bounds[name][1] <= bounds[name][0]]
    if not free:
        vector = np.asarray([bounds[name][0] for name in PARAMETER_NAMES], dtype=float)
        return np.repeat(vector[None, :], n_samples, axis=0)
    m = int(math.ceil(math.log2(max(n_samples, 1))))
    unit = qmc.Sobol(d=len(free), scramble=True, seed=seed).random_base2(m)[:n_samples]
    vectors = np.zeros((n_samples, len(PARAMETER_NAMES)), dtype=float)
    for j, i in enumerate(free):
        low, high = bounds[PARAMETER_NAMES[i]]
        vectors[:, i] = low + unit[:, j] * (high - low)
    for i in fixed:
        vectors[:, i] = bounds[PARAMETER_NAMES[i]][0]
    return vectors


def _expand_selected_temperature_detail(selected: pd.DataFrame, temperatures: Iterable[float]) -> pd.DataFrame:
    records: list[dict[str, Any]] = []
    temps = list(float(T) for T in temperatures)
    for _, row in selected.iterrows():
        full = json.loads(row["analytical_full_K_json"])
        off = json.loads(row["analytical_off_K_json"])
        for T, Kf, Ko in zip(temps, full, off):
            records.append(
                {
                    "candidate_id": row.candidate_id,
                    "transition_bracket": row.transition_bracket,
                    "T_K": T,
                    "analytical_full_K_init": float(Kf),
                    "analytical_plasticity_off_K_init": float(Ko),
                    "analytical_plastic_increment": float(Kf - Ko),
                }
            )
    return pd.DataFrame(records)


def select_by_bracket(results: pd.DataFrame, per_bracket_keep: int) -> pd.DataFrame:
    valid = results[results.analysis_valid.astype(bool)].copy()
    if valid.empty:
        return valid
    selected: list[pd.DataFrame] = []
    for bracket, group in valid.groupby("transition_bracket", sort=True):
        passed = group[group.screen_pass.astype(bool)].sort_values("objective")
        source = passed if not passed.empty else group.sort_values("objective")
        keep = source.head(per_bracket_keep).copy()
        keep["selection_basis"] = "screen_pass" if not passed.empty else "best_available_in_bracket"
        keep["rank_within_bracket"] = np.arange(1, len(keep) + 1)
        selected.append(keep)
    return pd.concat(selected, ignore_index=True).sort_values(
        ["coarse_transition_low_T_K", "objective"]
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=8192)
    ap.add_argument("--seed", type=int, default=9104701)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--batch-size", type=int, default=128)
    ap.add_argument("--per-bracket-keep", type=int, default=8)
    ap.add_argument("--temperatures", default="300 400 500 600 700 800 900 1000 1100")
    ap.add_argument("--cleavage-slope-mode", choices=["fixed_zero", "narrow"], default="fixed_zero")
    ap.add_argument("--Kdot", type=float, default=0.005)
    ap.add_argument("--dK", type=float, default=0.5)
    ap.add_argument("--Kmax", type=float, default=80.0)
    ap.add_argument("--out", type=Path, default=Path("runs/mpz_v9_10_4_7_analytical_downselect_v1"))
    args = ap.parse_args()

    out = args.out.resolve()
    checkpoints = out / "checkpoints"
    checkpoints.mkdir(parents=True, exist_ok=True)
    temperatures = parse_floats(args.temperatures)
    settings = AnalyticalSettings(
        temperatures_K=temperatures,
        Kdot_MPa_sqrt_m_s=float(args.Kdot),
        dK_MPa_sqrt_m=float(args.dK),
        Kmax_MPa_sqrt_m=float(args.Kmax),
    )
    bounds = search_bounds(args.cleavage_slope_mode)
    vectors = sobol_parameter_vectors(args.samples, bounds, args.seed)
    config = {
        "settings": asdict(settings),
        "cleavage_slope_mode": args.cleavage_slope_mode,
    }
    (out / "analytical_downselect_config.json").write_text(
        json.dumps(
            {
                **vars(args),
                "out": str(args.out),
                "settings": asdict(settings),
                "bounds": bounds,
                "parameter_names": PARAMETER_NAMES,
            },
            indent=2,
        )
    )

    print("=" * 80, flush=True)
    print("v9.10.4.7 analytical DBTT down-selection", flush=True)
    print(f"samples={args.samples} workers={args.workers} batch_size={args.batch_size}", flush=True)
    print(f"temperatures={temperatures}", flush=True)
    print(f"cleavage_slope_mode={args.cleavage_slope_mode}", flush=True)
    print(f"out={out}", flush=True)
    print("=" * 80, flush=True)

    all_batches: list[pd.DataFrame] = []
    started = time.perf_counter()
    for start in range(0, args.samples, args.batch_size):
        end = min(start + args.batch_size, args.samples)
        path = checkpoints / f"batch_{start:07d}_{end:07d}.csv"
        if path.exists():
            batch = pd.read_csv(path)
            print(f"[resume-batch] {start}:{end} rows={len(batch)}", flush=True)
        else:
            payloads = [(i, vectors[i], config) for i in range(start, end)]
            batch_started = time.perf_counter()
            if args.workers > 1:
                with ProcessPoolExecutor(max_workers=args.workers) as pool:
                    rows = list(pool.map(evaluate_candidate, payloads, chunksize=4))
            else:
                rows = [evaluate_candidate(payload) for payload in payloads]
            batch = pd.DataFrame(rows)
            batch.to_csv(path, index=False)
            valid = int(batch.get("analysis_valid", pd.Series(dtype=bool)).fillna(False).astype(bool).sum())
            passed = int(batch.get("screen_pass", pd.Series(dtype=bool)).fillna(False).astype(bool).sum())
            print(
                f"[batch] {start}:{end} valid={valid} passed={passed} "
                f"elapsed={time.perf_counter() - batch_started:.1f}s",
                flush=True,
            )
        all_batches.append(batch)

    results = pd.concat(all_batches, ignore_index=True)
    results = results.sort_values(["analysis_valid", "objective"], ascending=[False, True])
    selected = select_by_bracket(results, args.per_bracket_keep)
    results.to_csv(out / "analytical_all_candidates.csv", index=False)
    selected.to_csv(out / "analytical_promotion_manifest.csv", index=False)

    bracket_rows: list[dict[str, Any]] = []
    queue_rows: list[dict[str, Any]] = []
    bracket_root = out / "promoted_by_bracket"
    bracket_root.mkdir(exist_ok=True)
    if not selected.empty:
        detail = _expand_selected_temperature_detail(selected, temperatures)
        detail.to_csv(out / "analytical_promoted_temperature_detail.csv", index=False)
        for bracket, group in selected.groupby("transition_bracket", sort=True):
            folder = bracket_root / str(bracket)
            folder.mkdir(exist_ok=True)
            group.to_csv(folder / "manifest.csv", index=False)
            schedule_records: list[dict[str, Any]] = []
            for _, row in group.iterrows():
                schedule = [float(x) for x in json.loads(row.refinement_transition_temperatures_K)]
                for order, T in enumerate(schedule, start=1):
                    record = {
                        "candidate_id": row.candidate_id,
                        "transition_bracket": bracket,
                        "coarse_transition_low_T_K": float(row.coarse_transition_low_T_K),
                        "coarse_transition_high_T_K": float(row.coarse_transition_high_T_K),
                        "temperature_order": order,
                        "T_K": T,
                    }
                    schedule_records.append(record)
                    queue_rows.append(record)
            pd.DataFrame(schedule_records).to_csv(folder / "moving_interface_temperature_schedule.csv", index=False)
            source = results[results.transition_bracket == bracket]
            bracket_rows.append(
                {
                    "transition_bracket": bracket,
                    "coarse_transition_low_T_K": float(group.coarse_transition_low_T_K.iloc[0]),
                    "coarse_transition_high_T_K": float(group.coarse_transition_high_T_K.iloc[0]),
                    "n_valid_candidates": int(len(source)),
                    "n_screen_pass": int(source.screen_pass.astype(bool).sum()),
                    "n_promoted": int(len(group)),
                    "best_objective": float(group.objective.min()),
                    "four_detailed_temperatures_K": group.refinement_transition_temperatures_K.iloc[0],
                }
            )
    pd.DataFrame(bracket_rows).to_csv(out / "analytical_bracket_summary.csv", index=False)
    pd.DataFrame(queue_rows).to_csv(out / "moving_interface_queue.csv", index=False)

    summary = {
        "status": "V9_10_4_7_ANALYTICAL_DOWNSELECT_COMPLETE",
        "n_samples": int(args.samples),
        "n_valid": int(results.analysis_valid.fillna(False).astype(bool).sum()),
        "n_screen_pass": int(results.screen_pass.fillna(False).astype(bool).sum()),
        "n_promoted": int(len(selected)),
        "n_transition_brackets": int(selected.transition_bracket.nunique()) if not selected.empty else 0,
        "per_bracket_keep": int(args.per_bracket_keep),
        "moving_interface_temperatures_per_candidate": 4,
        "wall_time_s": time.perf_counter() - started,
        "next_stage_manifest": str(out / "analytical_promotion_manifest.csv"),
        "next_stage_queue": str(out / "moving_interface_queue.csv"),
    }
    (out / "analytical_downselect_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
