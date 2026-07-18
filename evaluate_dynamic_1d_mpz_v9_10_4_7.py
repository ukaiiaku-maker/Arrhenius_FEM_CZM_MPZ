#!/usr/bin/env python3
"""Evaluate analytically promoted DBTT candidates in the 1-D moving front.

Each candidate carries a four-temperature schedule spanning its own analytically
predicted 100 K transition bracket.  No optimization is performed here: the
stage simply evaluates the expensive full and plasticity-off moving-interface
responses at those four temperatures, ranks candidates within each bracket,
and writes the next promotion manifest.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import time
from typing import Any

import numpy as np
import pandas as pd

import optimize_mpz_v9_10_2_independent_shape_global as v102
from arrhenius_fracture.reduced_campaign_front_v9104 import (
    ReducedFrontSettings,
    simulate_reduced_response,
)

PARAMETER_NAMES = tuple(v102.PARAMETER_NAMES)


def _parameters_from_row(row: pd.Series) -> dict[str, float]:
    x = np.asarray([float(row[name]) for name in PARAMETER_NAMES], dtype=float)
    return v102.decode(x)


def _schedule_from_row(row: pd.Series) -> list[float]:
    value = row.get("refinement_transition_temperatures_K", "")
    schedule = [float(x) for x in json.loads(str(value))]
    if len(schedule) != 4:
        raise ValueError(
            f"candidate {row.get('candidate_id')} must have exactly four detailed temperatures; "
            f"found {len(schedule)}"
        )
    if not np.all(np.diff(schedule) > 0.0):
        raise ValueError(f"candidate {row.get('candidate_id')} has a non-monotone temperature schedule")
    return schedule


def _transition_metrics(
    temperatures: np.ndarray,
    full_K: np.ndarray,
    off_K: np.ndarray,
) -> dict[str, Any]:
    if not np.all(np.isfinite(full_K)) or not np.all(np.isfinite(off_K)):
        return {
            "moving_1d_valid": False,
            "moving_1d_objective": 1.0e12,
            "moving_1d_accept": False,
            "moving_1d_reason": "incomplete_first_passage",
        }
    rise = float(full_K[-1] - full_K[0])
    ratio = float(full_K[-1] / max(full_K[0], 1.0e-12))
    robust_ratio = float(np.min(full_K[2:]) / max(np.max(full_K[:2]), 1.0e-12))
    increments = np.diff(full_K)
    total_variation = float(np.sum(np.abs(increments)))
    positive_variation = float(np.sum(np.maximum(increments, 0.0)))
    monotonic_fraction = 1.0 if total_variation <= 1.0e-12 else positive_variation / total_variation
    off_ratio = float(off_K[-1] / max(off_K[0], 1.0e-12))
    plastic = full_K - off_K
    mechanistic_fraction = float((plastic[-1] - plastic[0]) / max(rise, 1.0e-12))

    normalized = (full_K - full_K[0]) / max(rise, 1.0e-12)

    def crossing(level: float) -> float:
        if normalized[0] >= level:
            return float(temperatures[0])
        for i in range(len(temperatures) - 1):
            y0 = float(normalized[i])
            y1 = float(normalized[i + 1])
            if y0 < level <= y1 and y1 > y0:
                fraction = (level - y0) / (y1 - y0)
                return float(temperatures[i] + fraction * (temperatures[i + 1] - temperatures[i]))
        return float("nan")

    T10 = crossing(0.10)
    T90 = crossing(0.90)
    width = float(T90 - T10) if np.isfinite(T10) and np.isfinite(T90) else float("inf")
    penalties = {
        "ratio": max(1.80 - ratio, 0.0) / 0.20,
        "robust_ratio": max(1.50 - robust_ratio, 0.0) / 0.20,
        "monotonicity": max(0.90 - monotonic_fraction, 0.0) / 0.10,
        "plasticity_off_ratio": max(off_ratio - 1.25, 0.0) / 0.10,
        "mechanistic_fraction": max(0.50 - mechanistic_fraction, 0.0) / 0.15,
        "transition_width": max(width - 100.0, 0.0) / 25.0 if np.isfinite(width) else 20.0,
    }
    objective = float(sum(value * value for value in penalties.values()))
    checks = [
        (ratio >= 1.80, "moving_1d_ratio_too_small"),
        (robust_ratio >= 1.50, "moving_1d_robust_ratio_too_small"),
        (monotonic_fraction >= 0.90, "moving_1d_transition_nonmonotone"),
        (off_ratio <= 1.25, "moving_1d_cleavage_only_T_dependence"),
        (mechanistic_fraction >= 0.50, "moving_1d_plastic_increment_too_small"),
        (width <= 100.0, "moving_1d_transition_wider_than_bracket"),
    ]
    accepted = True
    reason = "moving_1d_transition_passed"
    for passed, failure in checks:
        if not passed:
            accepted = False
            reason = failure
            break
    return {
        "moving_1d_valid": True,
        "moving_1d_objective": objective,
        "moving_1d_accept": accepted,
        "moving_1d_reason": reason,
        "moving_1d_edge_ratio": ratio,
        "moving_1d_robust_ratio": robust_ratio,
        "moving_1d_rise_MPa_sqrt_m": rise,
        "moving_1d_monotonic_fraction": monotonic_fraction,
        "moving_1d_plasticity_off_ratio": off_ratio,
        "moving_1d_mechanistic_fraction": mechanistic_fraction,
        "moving_1d_T10_K": T10,
        "moving_1d_T90_K": T90,
        "moving_1d_transition_width_K": width,
        **{f"moving_1d_penalty_{key}": value for key, value in penalties.items()},
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=Path("runs/mpz_v9_10_4_7_dynamic_1d_v1"))
    ap.add_argument("--target-extension-um", type=float, default=5.0)
    ap.add_argument("--per-bracket-keep", type=int, default=2)
    ap.add_argument("--heartbeat-candidates", type=int, default=1)
    args = ap.parse_args()

    manifest = pd.read_csv(args.manifest)
    out = args.out.resolve()
    checkpoints = out / "checkpoints"
    checkpoints.mkdir(parents=True, exist_ok=True)
    settings = ReducedFrontSettings(target_extension_um=float(args.target_extension_um))
    started = time.perf_counter()
    summary_rows: list[dict[str, Any]] = []
    temperature_rows: list[dict[str, Any]] = []
    event_rows: list[dict[str, Any]] = []

    print("=" * 80, flush=True)
    print("v9.10.4.7 dynamic four-temperature 1-D evaluation", flush=True)
    print(f"manifest={args.manifest} candidates={len(manifest)}", flush=True)
    print(f"target_extension_um={args.target_extension_um}", flush=True)
    print(f"out={out}", flush=True)
    print("=" * 80, flush=True)

    for position, (_, row) in enumerate(manifest.iterrows(), start=1):
        cid = str(row.candidate_id)
        checkpoint = checkpoints / f"{cid}.json"
        if checkpoint.exists():
            payload = json.loads(checkpoint.read_text())
            if payload.get("status") == "COMPLETE":
                summary_rows.append(payload["summary"])
                temperature_rows.extend(payload.get("temperature_detail", []))
                event_rows.extend(payload.get("event_detail", []))
                print(f"[resume] {cid} ({position}/{len(manifest)})", flush=True)
                continue

        p = _parameters_from_row(row)
        schedule = np.asarray(_schedule_from_row(row), dtype=float)
        print(
            f"[candidate] {cid} ({position}/{len(manifest)}) "
            f"bracket={row.transition_bracket} temperatures={schedule.tolist()}",
            flush=True,
        )
        full_K: list[float] = []
        off_K: list[float] = []
        local_temperature: list[dict[str, Any]] = []
        local_events: list[dict[str, Any]] = []
        for temperature_index, T in enumerate(schedule, start=1):
            t0 = time.perf_counter()
            full = simulate_reduced_response(p, float(T), settings, mode="full")
            off = simulate_reduced_response(p, float(T), settings, mode="plasticity_off")
            full_events = full.pop("events", [])
            off.pop("events", None)
            Kf = float(full.get("K_init_proxy", np.nan))
            Ko = float(off.get("K_init_proxy", np.nan))
            full_K.append(Kf)
            off_K.append(Ko)
            record = {
                "candidate_id": cid,
                "transition_bracket": row.transition_bracket,
                "temperature_order": temperature_index,
                "T_K": float(T),
                "full_K_init": Kf,
                "plasticity_off_K_init": Ko,
                "plastic_increment": Kf - Ko,
                "full_completed": bool(full.get("completed", False)),
                "off_completed": bool(off.get("completed", False)),
                "full_internal_steps": int(full.get("internal_steps", -1)),
                "off_internal_steps": int(off.get("internal_steps", -1)),
                "elapsed_s": time.perf_counter() - t0,
            }
            local_temperature.append(record)
            local_events.extend(
                {"candidate_id": cid, "transition_bracket": row.transition_bracket, "T_K": float(T), **event}
                for event in full_events
            )
            print(
                f"[temperature] {cid} T={T:.3f}K full={Kf:.6g} off={Ko:.6g} "
                f"steps={record['full_internal_steps']} elapsed={record['elapsed_s']:.2f}s",
                flush=True,
            )

        metrics = _transition_metrics(schedule, np.asarray(full_K), np.asarray(off_K))
        summary = {
            **row.to_dict(),
            **metrics,
            "moving_1d_temperatures_K": json.dumps(schedule.tolist()),
            "moving_1d_full_K_json": json.dumps(full_K),
            "moving_1d_off_K_json": json.dumps(off_K),
        }
        payload = {
            "status": "COMPLETE",
            "summary": summary,
            "temperature_detail": local_temperature,
            "event_detail": local_events,
        }
        checkpoint.write_text(json.dumps(payload, indent=2, allow_nan=True))
        summary_rows.append(summary)
        temperature_rows.extend(local_temperature)
        event_rows.extend(local_events)
        print(
            f"[candidate-complete] {cid} objective={metrics['moving_1d_objective']:.6g} "
            f"accepted={metrics['moving_1d_accept']} reason={metrics['moving_1d_reason']}",
            flush=True,
        )

    results = pd.DataFrame(summary_rows).sort_values(
        ["coarse_transition_low_T_K", "moving_1d_objective"]
    )
    promoted: list[pd.DataFrame] = []
    for bracket, group in results.groupby("transition_bracket", sort=True):
        accepted = group[group.moving_1d_accept.astype(bool)].sort_values("moving_1d_objective")
        source = accepted if not accepted.empty else group.sort_values("moving_1d_objective")
        keep = source.head(args.per_bracket_keep).copy()
        keep["moving_1d_selection_basis"] = (
            "moving_1d_accept" if not accepted.empty else "best_available_in_bracket"
        )
        keep["moving_1d_rank_within_bracket"] = np.arange(1, len(keep) + 1)
        promoted.append(keep)
    promotion = pd.concat(promoted, ignore_index=True) if promoted else results.head(0).copy()

    results.to_csv(out / "dynamic_1d_all_candidates.csv", index=False)
    results[results.moving_1d_accept.astype(bool)].to_csv(
        out / "dynamic_1d_accepted.csv", index=False
    )
    promotion.to_csv(out / "short_growth_promotion_manifest.csv", index=False)
    pd.DataFrame(temperature_rows).to_csv(out / "dynamic_1d_temperature_detail.csv", index=False)
    pd.DataFrame(event_rows).to_csv(out / "dynamic_1d_event_detail.csv", index=False)

    bracket_summary = (
        results.groupby("transition_bracket", as_index=False)
        .agg(
            n_candidates=("candidate_id", "count"),
            n_accepted=("moving_1d_accept", "sum"),
            best_objective=("moving_1d_objective", "min"),
            transition_low_T_K=("coarse_transition_low_T_K", "first"),
            transition_high_T_K=("coarse_transition_high_T_K", "first"),
        )
    )
    bracket_summary.to_csv(out / "dynamic_1d_bracket_summary.csv", index=False)
    report = {
        "status": "V9_10_4_7_DYNAMIC_1D_COMPLETE",
        "n_candidates": int(len(results)),
        "n_accepted": int(results.moving_1d_accept.astype(bool).sum()),
        "n_promoted": int(len(promotion)),
        "n_transition_brackets": int(results.transition_bracket.nunique()),
        "temperatures_per_candidate": 4,
        "wall_time_s": time.perf_counter() - started,
        "next_stage_manifest": str(out / "short_growth_promotion_manifest.csv"),
    }
    (out / "dynamic_1d_summary.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
