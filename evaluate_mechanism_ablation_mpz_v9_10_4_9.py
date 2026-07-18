#!/usr/bin/env python3
"""Bracket-balanced reduced-front mechanism ablation for 2-D promotion.

For every candidate and each of its four candidate-specific transition
temperatures, evaluate the requested constitutive modes:

    full, plasticity_off, blunting_off, backstress_off, shielding_off.

An optional ``background_field_off`` diagnostic uses the full constitutive mode
with ``rho0_m2=0``.  Emitted defects still undergo Peierls transport and Taylor
retention, so this isolates sensitivity to the imposed pre-existing density
field without suppressing emission-generated plasticity.

The reported sensitivity fractions are finite-difference ablations, not an
additive energy decomposition; interactions between mechanisms are retained.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
import json
import math
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
REQUESTED_MODES = (
    "full",
    "plasticity_off",
    "blunting_off",
    "backstress_off",
    "shielding_off",
)
OPTIONAL_MODE = "background_field_off"


def parse_modes(text: str) -> tuple[str, ...]:
    modes = tuple(x for x in str(text).replace(",", " ").split() if x)
    allowed = set(REQUESTED_MODES) | {OPTIONAL_MODE}
    unknown = sorted(set(modes).difference(allowed))
    if unknown:
        raise ValueError(f"unknown ablation modes: {unknown}")
    if len(set(modes)) != len(modes):
        raise ValueError("ablation modes must be unique")
    return modes


def _parameters_from_row(row: pd.Series) -> dict[str, float]:
    x = np.asarray([float(row[name]) for name in PARAMETER_NAMES], dtype=float)
    return v102.decode(x)


def _schedule_from_row(row: pd.Series) -> np.ndarray:
    value = row.get("refinement_transition_temperatures_K", "")
    schedule = np.asarray([float(x) for x in json.loads(str(value))], dtype=float)
    if schedule.size != 4:
        raise ValueError(
            f"candidate {row.get('candidate_id')} must have four transition temperatures; "
            f"found {schedule.size}"
        )
    if not np.all(np.diff(schedule) > 0.0):
        raise ValueError(f"candidate {row.get('candidate_id')} has a nonmonotone schedule")
    return schedule


def curve_metrics(temperatures: np.ndarray, toughness: np.ndarray) -> dict[str, Any]:
    T = np.asarray(temperatures, dtype=float)
    K = np.asarray(toughness, dtype=float)
    if T.size != 4 or K.size != 4 or not np.all(np.isfinite(K)):
        return {
            "curve_valid": False,
            "low_endpoint_K": float("nan"),
            "high_endpoint_K": float("nan"),
            "rise_MPa_sqrt_m": float("nan"),
            "endpoint_ratio": float("nan"),
            "monotonic_fraction": float("nan"),
            "T10_K": float("nan"),
            "T90_K": float("nan"),
            "transition_width_K": float("nan"),
        }
    low = float(K[0])
    high = float(K[-1])
    rise = high - low
    ratio = high / max(low, 1.0e-12)
    increments = np.diff(K)
    variation = float(np.sum(np.abs(increments)))
    positive = float(np.sum(np.maximum(increments, 0.0)))
    monotonic = 1.0 if variation <= 1.0e-12 else positive / variation

    T10 = float("nan")
    T90 = float("nan")
    width = float("nan")
    if rise > 0.0:
        normalized = (K - low) / rise

        def crossing(level: float) -> float:
            if normalized[0] >= level:
                return float(T[0])
            for i in range(len(T) - 1):
                y0 = float(normalized[i])
                y1 = float(normalized[i + 1])
                if y0 < level <= y1 and y1 > y0:
                    f = (level - y0) / (y1 - y0)
                    return float(T[i] + f * (T[i + 1] - T[i]))
            return float("nan")

        T10 = crossing(0.10)
        T90 = crossing(0.90)
        if np.isfinite(T10) and np.isfinite(T90) and T90 >= T10:
            width = T90 - T10
    return {
        "curve_valid": True,
        "low_endpoint_K": low,
        "high_endpoint_K": high,
        "rise_MPa_sqrt_m": rise,
        "endpoint_ratio": ratio,
        "monotonic_fraction": monotonic,
        "T10_K": T10,
        "T90_K": T90,
        "transition_width_K": width,
    }


def _solver_configuration(
    mode: str,
    settings: ReducedFrontSettings,
) -> tuple[str, ReducedFrontSettings, str]:
    if mode == OPTIONAL_MODE:
        return "full", replace(settings, rho0_m2=0.0), "full_with_rho0_zero"
    return mode, settings, mode


def _first_event(result: dict[str, Any]) -> dict[str, Any]:
    events = result.get("events", [])
    if not events:
        return {}
    return dict(events[0])


def _safe_float(value: Any, default: float = float("nan")) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out


def _safe_div(numerator: float, denominator: float) -> float:
    if not np.isfinite(numerator) or not np.isfinite(denominator):
        return float("nan")
    return float(numerator / max(abs(denominator), 1.0e-12))


def _mode_summary(
    candidate_id: str,
    bracket: str,
    mode: str,
    temperatures: np.ndarray,
    records: list[dict[str, Any]],
) -> dict[str, Any]:
    ordered = sorted(records, key=lambda row: float(row["T_K"]))
    K = np.asarray([_safe_float(row.get("K_init_proxy")) for row in ordered], dtype=float)
    metrics = curve_metrics(temperatures, K)
    low = ordered[0]
    high = ordered[-1]
    return {
        "candidate_id": candidate_id,
        "transition_bracket": bracket,
        "mode": mode,
        **metrics,
        "completed_all_temperatures": bool(all(bool(row.get("completed", False)) for row in ordered)),
        "total_internal_steps": int(sum(int(row.get("internal_steps", 0)) for row in ordered)),
        "max_K_shield_MPa_sqrt_m": float(
            np.nanmax([_safe_float(row.get("max_K_shield_MPa_sqrt_m"), 0.0) for row in ordered])
        ),
        "max_sigma_back_Pa": float(
            np.nanmax([_safe_float(row.get("max_sigma_back_Pa"), 0.0) for row in ordered])
        ),
        "low_cumulative_emitted": _safe_float(low.get("cumulative_emitted"), 0.0),
        "high_cumulative_emitted": _safe_float(high.get("cumulative_emitted"), 0.0),
        "low_mobile_count": _safe_float(low.get("mobile_count"), 0.0),
        "high_mobile_count": _safe_float(high.get("mobile_count"), 0.0),
        "low_retained_count": _safe_float(low.get("retained_count"), 0.0),
        "high_retained_count": _safe_float(high.get("retained_count"), 0.0),
        "low_r_eff_m": _safe_float(low.get("r_eff_m")),
        "high_r_eff_m": _safe_float(high.get("r_eff_m")),
        "high_to_low_r_eff_ratio": _safe_div(
            _safe_float(high.get("r_eff_m")),
            _safe_float(low.get("r_eff_m")),
        ),
        "K_json": json.dumps(K.tolist()),
        "temperatures_K_json": json.dumps([float(x) for x in temperatures]),
    }


def candidate_sensitivity(
    row: pd.Series,
    mode_rows: dict[str, dict[str, Any]],
) -> dict[str, Any]:
    cid = str(row.candidate_id)
    bracket = str(row.transition_bracket)
    full = mode_rows.get("full", {})
    off = mode_rows.get("plasticity_off", {})
    blunt = mode_rows.get("blunting_off", {})
    back = mode_rows.get("backstress_off", {})
    shield = mode_rows.get("shielding_off", {})
    background = mode_rows.get(OPTIONAL_MODE, {})

    full_rise = _safe_float(full.get("rise_MPa_sqrt_m"))
    off_rise = _safe_float(off.get("rise_MPa_sqrt_m"))
    blunt_rise = _safe_float(blunt.get("rise_MPa_sqrt_m"))
    back_rise = _safe_float(back.get("rise_MPa_sqrt_m"))
    shield_rise = _safe_float(shield.get("rise_MPa_sqrt_m"))
    background_rise = _safe_float(background.get("rise_MPa_sqrt_m"))

    emission_coupled_rise = full_rise - off_rise
    emission_fraction = _safe_div(emission_coupled_rise, full_rise)
    blunting_fraction = _safe_div(full_rise - blunt_rise, full_rise)
    shielding_fraction = _safe_div(full_rise - shield_rise, full_rise)
    backstress_fraction = _safe_div(full_rise - back_rise, full_rise)
    background_fraction = _safe_div(full_rise - background_rise, full_rise)
    background_retained = _safe_div(background_rise, full_rise)

    full_ratio = _safe_float(full.get("endpoint_ratio"))
    off_ratio = _safe_float(off.get("endpoint_ratio"))
    low_K = _safe_float(full.get("low_endpoint_K"))
    high_K = _safe_float(full.get("high_endpoint_K"))
    monotonic = _safe_float(full.get("monotonic_fraction"))
    width = _safe_float(full.get("transition_width_K"))
    max_shield = _safe_float(full.get("max_K_shield_MPa_sqrt_m"), 0.0)
    emitted_growth = _safe_float(full.get("high_cumulative_emitted"), 0.0) - _safe_float(
        full.get("low_cumulative_emitted"), 0.0
    )

    penalties = {
        "factor_two": max(2.0 - full_ratio, 0.0) / 0.25 if np.isfinite(full_ratio) else 20.0,
        "low_floor": max(8.0 - low_K, 0.0) / 2.0 if np.isfinite(low_K) else 20.0,
        "low_ceiling": max(low_K - 25.0, 0.0) / 3.0 if np.isfinite(low_K) else 20.0,
        "high_ceiling": max(high_K - 70.0, 0.0) / 5.0 if np.isfinite(high_K) else 20.0,
        "opening_baseline": max(off_ratio - 1.25, 0.0) / 0.10 if np.isfinite(off_ratio) else 20.0,
        "emission_fraction": max(0.60 - emission_fraction, 0.0) / 0.15 if np.isfinite(emission_fraction) else 20.0,
        "blunting_control": max(0.50 - blunting_fraction, 0.0) / 0.15 if np.isfinite(blunting_fraction) else 20.0,
        "shielding_dependence": max(abs(shielding_fraction) - 0.20, 0.0) / 0.10 if np.isfinite(shielding_fraction) else 20.0,
        "background_dependence": max(0.75 - background_retained, 0.0) / 0.10 if np.isfinite(background_retained) else 20.0,
        "monotonicity": max(0.90 - monotonic, 0.0) / 0.10 if np.isfinite(monotonic) else 20.0,
        "transition_width": max(width - 100.0, 0.0) / 25.0 if np.isfinite(width) else 20.0,
        "absolute_shield": max(max_shield - 0.5, 0.0) / 0.25,
    }
    score = float(sum(value * value for value in penalties.values()))
    priority_checks = [
        np.isfinite(full_ratio) and full_ratio >= 2.0,
        np.isfinite(low_K) and 8.0 <= low_K <= 25.0,
        np.isfinite(high_K) and high_K <= 70.0,
        np.isfinite(off_ratio) and off_ratio <= 1.25,
        np.isfinite(emission_fraction) and emission_fraction >= 0.60,
        np.isfinite(blunting_fraction) and blunting_fraction >= 0.50,
        np.isfinite(shielding_fraction) and abs(shielding_fraction) <= 0.20,
        np.isfinite(background_retained) and background_retained >= 0.75,
        np.isfinite(monotonic) and monotonic >= 0.90,
        np.isfinite(width) and width <= 100.0,
        max_shield <= 0.5,
        emitted_growth > 0.0,
    ]
    priority = bool(all(priority_checks))
    reason = (
        "emission_opening_2d_priority"
        if priority
        else "review_ablation_sensitivities_before_2d"
    )

    return {
        **row.to_dict(),
        "ablation_full_rise_MPa_sqrt_m": full_rise,
        "ablation_opening_only_rise_MPa_sqrt_m": off_rise,
        "ablation_emission_coupled_rise_MPa_sqrt_m": emission_coupled_rise,
        "ablation_emission_fraction_of_full_rise": emission_fraction,
        "ablation_blunting_sensitivity_fraction": blunting_fraction,
        "ablation_shielding_sensitivity_fraction": shielding_fraction,
        "ablation_backstress_sensitivity_fraction": backstress_fraction,
        "ablation_background_field_sensitivity_fraction": background_fraction,
        "ablation_background_off_retained_rise_fraction": background_retained,
        "ablation_full_endpoint_ratio": full_ratio,
        "ablation_plasticity_off_endpoint_ratio": off_ratio,
        "ablation_full_transition_width_K": width,
        "ablation_full_max_K_shield_MPa_sqrt_m": max_shield,
        "ablation_high_minus_low_cumulative_emitted": emitted_growth,
        "two_d_emission_opening_score": score,
        "two_d_emission_opening_priority": priority,
        "two_d_promotion_reason": reason,
        "ablation_nonadditive_interactions_retained": True,
        **{f"two_d_penalty_{key}": value for key, value in penalties.items()},
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--target-extension-um", type=float, default=5.0)
    ap.add_argument(
        "--modes",
        default=" ".join(REQUESTED_MODES),
        help="space- or comma-separated requested modes",
    )
    ap.add_argument("--include-background-field-off", action="store_true")
    ap.add_argument("--expected-candidates", type=int, default=6)
    args = ap.parse_args()

    modes = list(parse_modes(args.modes))
    if args.include_background_field_off and OPTIONAL_MODE not in modes:
        modes.append(OPTIONAL_MODE)
    if "full" not in modes or "plasticity_off" not in modes:
        raise ValueError("full and plasticity_off are required for the sensitivity summary")

    manifest = pd.read_csv(args.manifest)
    if len(manifest) != int(args.expected_candidates):
        raise ValueError(
            f"expected {args.expected_candidates} candidates; manifest contains {len(manifest)}"
        )

    out = args.out.resolve()
    checkpoints = out / "checkpoints"
    checkpoints.mkdir(parents=True, exist_ok=True)
    base_settings = ReducedFrontSettings(target_extension_um=float(args.target_extension_um))
    started = time.perf_counter()

    temperature_rows: list[dict[str, Any]] = []
    event_rows: list[dict[str, Any]] = []
    mode_summary_rows: list[dict[str, Any]] = []
    candidate_summary_rows: list[dict[str, Any]] = []

    print("=" * 80, flush=True)
    print("v9.10.4.9 six-candidate mechanism ablation", flush=True)
    print(f"manifest={args.manifest} candidates={len(manifest)}", flush=True)
    print(f"modes={modes}", flush=True)
    print(f"target_extension_um={args.target_extension_um}", flush=True)
    print(f"out={out}", flush=True)
    print("=" * 80, flush=True)

    for candidate_position, (_, candidate) in enumerate(manifest.iterrows(), start=1):
        cid = str(candidate.candidate_id)
        bracket = str(candidate.transition_bracket)
        p = _parameters_from_row(candidate)
        temperatures = _schedule_from_row(candidate)
        print(
            f"[candidate] {cid} ({candidate_position}/{len(manifest)}) "
            f"bracket={bracket} temperatures={temperatures.tolist()}",
            flush=True,
        )
        candidate_mode_rows: dict[str, dict[str, Any]] = {}

        for mode in modes:
            solver_mode, settings, mode_description = _solver_configuration(mode, base_settings)
            local_records: list[dict[str, Any]] = []
            for temperature_position, T in enumerate(temperatures, start=1):
                mode_folder = checkpoints / cid / mode
                mode_folder.mkdir(parents=True, exist_ok=True)
                checkpoint = mode_folder / f"T{float(T):09.3f}.json"
                if checkpoint.exists():
                    payload = json.loads(checkpoint.read_text())
                    result = payload["result"]
                    print(f"[resume] {cid} mode={mode} T={T:.3f}K", flush=True)
                else:
                    t0 = time.perf_counter()
                    result = simulate_reduced_response(p, float(T), settings, mode=solver_mode)
                    payload = {
                        "status": "COMPLETE",
                        "candidate_id": cid,
                        "transition_bracket": bracket,
                        "mode": mode,
                        "solver_mode": solver_mode,
                        "mode_description": mode_description,
                        "T_K": float(T),
                        "rho0_m2": float(settings.rho0_m2),
                        "elapsed_s": time.perf_counter() - t0,
                        "result": result,
                    }
                    checkpoint.write_text(json.dumps(payload, indent=2, allow_nan=True))
                    print(
                        f"[result] {cid} mode={mode} T={T:.3f}K "
                        f"K={_safe_float(result.get('K_init_proxy')):.6g} "
                        f"steps={int(result.get('internal_steps', -1))} "
                        f"elapsed={payload['elapsed_s']:.2f}s",
                        flush=True,
                    )

                event = _first_event(result)
                record = {
                    "candidate_id": cid,
                    "transition_bracket": bracket,
                    "ablation_selection_basis": candidate.get("ablation_selection_basis", ""),
                    "mode": mode,
                    "solver_mode": solver_mode,
                    "mode_description": mode_description,
                    "temperature_order": temperature_position,
                    "T_K": float(T),
                    "rho0_m2": float(settings.rho0_m2),
                    "completed": bool(result.get("completed", False)),
                    "K_init_proxy": _safe_float(result.get("K_init_proxy")),
                    "internal_steps": int(result.get("internal_steps", -1)),
                    "n_events": int(result.get("n_events", 0)),
                    "max_K_shield_MPa_sqrt_m": _safe_float(
                        result.get("max_K_shield_MPa_sqrt_m"),
                        _safe_float(event.get("K_shield_MPa_sqrt_m"), 0.0),
                    ),
                    "max_sigma_back_Pa": _safe_float(
                        result.get("max_sigma_back_Pa"),
                        _safe_float(event.get("sigma_back_max_Pa"), 0.0),
                    ),
                    "cumulative_emitted": _safe_float(
                        result.get("cumulative_emitted"),
                        _safe_float(event.get("cumulative_emitted"), 0.0),
                    ),
                    "final_mobile_count": _safe_float(
                        result.get("final_mobile_count"),
                        _safe_float(event.get("mobile_count"), 0.0),
                    ),
                    "final_retained_count": _safe_float(
                        result.get("final_retained_count"),
                        _safe_float(event.get("retained_count"), 0.0),
                    ),
                    "mobile_count": _safe_float(event.get("mobile_count"), 0.0),
                    "retained_count": _safe_float(event.get("retained_count"), 0.0),
                    "available_sites": _safe_float(event.get("available_sites")),
                    "r_eff_m": _safe_float(event.get("r_eff_m")),
                    "sigma_open_Pa": _safe_float(event.get("sigma_open_Pa")),
                    "sigma_cleave_Pa": _safe_float(event.get("sigma_cleave_Pa")),
                    "sigma_emit_max_Pa": _safe_float(event.get("sigma_emit_max_Pa")),
                }
                local_records.append(record)
                temperature_rows.append(record)
                for event_row in result.get("events", []):
                    event_rows.append(
                        {
                            "candidate_id": cid,
                            "transition_bracket": bracket,
                            "mode": mode,
                            "solver_mode": solver_mode,
                            "T_K": float(T),
                            "rho0_m2": float(settings.rho0_m2),
                            **event_row,
                        }
                    )

            summary = _mode_summary(cid, bracket, mode, temperatures, local_records)
            mode_summary_rows.append(summary)
            candidate_mode_rows[mode] = summary
            print(
                f"[mode-complete] {cid} mode={mode} ratio={summary['endpoint_ratio']:.6g} "
                f"rise={summary['rise_MPa_sqrt_m']:.6g}",
                flush=True,
            )

        candidate_summary = candidate_sensitivity(candidate, candidate_mode_rows)
        candidate_summary_rows.append(candidate_summary)
        print(
            f"[candidate-complete] {cid} score={candidate_summary['two_d_emission_opening_score']:.6g} "
            f"priority={candidate_summary['two_d_emission_opening_priority']}",
            flush=True,
        )

        pd.DataFrame(temperature_rows).to_csv(
            out / "mechanism_ablation_temperature_detail.partial.csv", index=False
        )
        pd.DataFrame(mode_summary_rows).to_csv(
            out / "mechanism_ablation_mode_summary.partial.csv", index=False
        )
        pd.DataFrame(candidate_summary_rows).to_csv(
            out / "mechanism_ablation_candidate_summary.partial.csv", index=False
        )

    temperature_df = pd.DataFrame(temperature_rows)
    event_df = pd.DataFrame(event_rows)
    mode_df = pd.DataFrame(mode_summary_rows)
    candidate_df = pd.DataFrame(candidate_summary_rows).sort_values(
        ["two_d_emission_opening_score", "coarse_transition_low_T_K"]
    )
    priority_df = candidate_df[
        candidate_df.two_d_emission_opening_priority.fillna(False).astype(bool)
    ].copy()

    temperature_df.to_csv(out / "mechanism_ablation_temperature_detail.csv", index=False)
    event_df.to_csv(out / "mechanism_ablation_event_detail.csv", index=False)
    mode_df.to_csv(out / "mechanism_ablation_mode_summary.csv", index=False)
    candidate_df.to_csv(out / "mechanism_ablation_candidate_summary.csv", index=False)
    candidate_df.to_csv(out / "two_d_candidate_ranking.csv", index=False)
    priority_df.to_csv(out / "two_d_emission_opening_priority.csv", index=False)

    report = {
        "status": "V9_10_4_9_MECHANISM_ABLATION_COMPLETE",
        "n_candidates": int(len(manifest)),
        "n_modes": int(len(modes)),
        "modes": modes,
        "temperatures_per_candidate": 4,
        "solver_cases": int(len(manifest) * len(modes) * 4),
        "n_two_d_emission_opening_priority": int(len(priority_df)),
        "wall_time_s": time.perf_counter() - started,
        "ranking": str(out / "two_d_candidate_ranking.csv"),
        "priority_manifest": str(out / "two_d_emission_opening_priority.csv"),
        "interpretation": (
            "Ablation sensitivities are signed finite differences and are not additive; "
            "mechanism interactions are retained."
        ),
    }
    (out / "mechanism_ablation_summary.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
