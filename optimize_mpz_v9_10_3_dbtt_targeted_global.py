#!/usr/bin/env python3
"""Target-aware DBTT search over the full v9.10.2 independent-shape space.

Version 9.10.2 established that independent EXP-floor shapes are necessary, but
its DBTT objective only rewarded a plateau rise, an initiation rise, and a high-
temperature R-curve threshold.  That objective admitted an essentially zero
300 K toughness shelf and therefore assigned objective zero to a large,
physically unsuitable region.

Version 9.10.3 retains the complete 29-dimensional v9.10.2 search domain and the
same unified mobile/retained kinetics.  It changes only the DBTT response
objective and acceptance gate:

* use the temperature-resolved DBTT targets already stored in
  ``mpz_three_class_design_targets.csv``;
* add a smooth guard against a vanishing or excessive low-temperature shelf;
* require a finite brittle shelf and a developed high-temperature branch for
  promotion;
* impose no shielding, blunting, source, Peierls, or Taylor dominance rule.

Every restart still begins from a fresh full Sobol population.  No previous
shortlist is used as an initial population.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
import sys
from typing import Any

import numpy as np
import pandas as pd

import optimize_mpz_v9_10_unified_global as base
import optimize_mpz_v9_10_2_independent_shape_global as v102


LOW_SHELF_MIN = 8.0
LOW_SHELF_MAX = 25.0
LOW_PLATEAU_MAX = 28.0
LOW_RCURVE_MAX = 3.0
HIGH_INIT_MIN = 25.0
HIGH_PLATEAU_MIN = 35.0
PLATEAU_RISE_MIN = 15.0
HIGH_RCURVE_MIN = 5.0


def _finite(value: Any) -> bool:
    try:
        return bool(np.isfinite(float(value)))
    except (TypeError, ValueError):
        return False


def dbtt_target_components(merged: pd.DataFrame) -> dict[str, float]:
    """Return continuous target-aware DBTT response losses.

    Rows with incomplete/nonfinite predictions are skipped here; the completion
    loss is handled by the objective.  The low-shelf guard is intentionally a
    response constraint, not a mechanism constraint.
    """
    components: dict[str, float] = {
        "K_init_loss": 0.0,
        "K_plateau_loss": 0.0,
        "early_rise_loss": 0.0,
        "plateau_rise_loss": 0.0,
        "delta_window_loss": 0.0,
        "DBTT_low_shelf_guard_loss": 0.0,
        "DBTT_high_branch_guard_loss": 0.0,
        "DBTT_transition_guard_loss": 0.0,
    }

    valid = merged.copy()
    required = [
        "K_init_proxy",
        "K_plateau_proxy",
        "delta_KR_proxy",
        "early_rise_per_100um_proxy",
        "plateau_rise_per_100um_proxy",
    ]
    valid = valid[
        np.logical_and.reduce([
            pd.to_numeric(valid[name], errors="coerce").notna().to_numpy()
            for name in required
        ])
    ].copy()
    if valid.empty:
        return components

    for _, row in valid.iterrows():
        w = float(row.weight)
        components["K_init_loss"] += w * base.core.huber(
            (float(row.K_init_proxy) - float(row.K_init_target))
            / max(float(row.K_init_scale), 1.0e-9)
        )
        components["K_plateau_loss"] += w * base.core.huber(
            (float(row.K_plateau_proxy) - float(row.K_plateau_target))
            / max(float(row.K_plateau_scale), 1.0e-9)
        )
        components["early_rise_loss"] += w * base.core.huber(
            (
                float(row.early_rise_per_100um_proxy)
                - float(row.early_rise_per_100um_target)
            )
            / max(float(row.early_rise_scale), 1.0e-9)
        )
        components["plateau_rise_loss"] += w * base.core.huber(
            (
                float(row.plateau_rise_per_100um_proxy)
                - float(row.plateau_rise_per_100um_target)
            )
            / max(float(row.plateau_rise_scale), 1.0e-9)
        )
        if float(row.delta_KR_proxy) < float(row.delta_KR_min):
            components["delta_window_loss"] += w * base.core.huber(
                (float(row.delta_KR_min) - float(row.delta_KR_proxy)) / 1.5
            )
        elif float(row.delta_KR_proxy) > float(row.delta_KR_max):
            components["delta_window_loss"] += w * base.core.huber(
                (float(row.delta_KR_proxy) - float(row.delta_KR_max)) / 1.5
            )

    ordered = valid.sort_values("T_K")
    low = ordered.iloc[0]
    high = ordered[ordered.T_K >= 900.0]
    if high.empty:
        high = ordered.tail(max(1, len(ordered) // 2))

    low_init = float(low.K_init_proxy)
    low_plateau = float(low.K_plateau_proxy)
    high_init = float(high.K_init_proxy.median())
    high_plateau = float(high.K_plateau_proxy.median())
    high_dkr = float(high.delta_KR_proxy.median())
    plateau_rise = high_plateau - low_plateau

    # Strongly remove the zero-shelf shortcut while leaving a broad admissible
    # brittle shelf.  The target data continue to rank candidates inside it.
    components["DBTT_low_shelf_guard_loss"] = (
        25.0
        * base.core.huber(max(LOW_SHELF_MIN - low_init, 0.0) / 2.0)
        + 10.0
        * base.core.huber(max(low_init - LOW_SHELF_MAX, 0.0) / 3.0)
        + 10.0
        * base.core.huber(max(low_plateau - LOW_PLATEAU_MAX, 0.0) / 3.0)
    )
    components["DBTT_high_branch_guard_loss"] = (
        5.0 * base.core.huber(max(HIGH_INIT_MIN - high_init, 0.0) / 3.0)
        + 5.0 * base.core.huber(max(HIGH_PLATEAU_MIN - high_plateau, 0.0) / 4.0)
        + 5.0 * base.core.huber(max(HIGH_RCURVE_MIN - high_dkr, 0.0) / 2.0)
    )
    components["DBTT_transition_guard_loss"] = 5.0 * base.core.huber(
        max(PLATEAU_RISE_MIN - plateau_rise, 0.0) / 3.0
    )
    return components


class DBTTTargetObjective:
    """Full v9.10.2 parameter search with a target-aware DBTT objective."""

    def __init__(self, settings: base.ZeroDSettings):
        if settings.target_class != "DBTT":
            raise ValueError("v9.10.3 target-aware optimizer is DBTT-only")
        self.s = settings
        self.bounds = np.asarray(v102.bounds_list(), dtype=float)

    def evaluate(self, x: np.ndarray, details: bool = False) -> dict[str, Any]:
        x = np.asarray(x, dtype=float)
        if not np.all(np.isfinite(x)):
            return {"objective": 1.0e12, "nonfinite_parameter_vector": True}
        outside = np.maximum(self.bounds[:, 0] - x, 0.0) + np.maximum(
            x - self.bounds[:, 1], 0.0
        )
        if np.any(outside > 0.0):
            return {
                "objective": 1.0e10 + 1.0e7 * float(np.sum(outside**2))
            }

        p = v102.decode(x)
        model = v102.build_model(p, self.s.Tref_K)
        stress_grid = np.linspace(0.0, 30.0e9, 31)
        order_margin = base.core.barrier_order_margin_eV(
            model, self.s.temperatures, stress_grid
        )
        raw = [
            model.raw_zero_stress_barrier_eV(mechanism, T)
            for mechanism in ("peierls", "taylor")
            for T in self.s.temperatures
        ]
        min_raw = float(np.min(raw))
        if order_margin < -1.0e-9 or min_raw <= 0.0:
            return {
                "objective": 1.0e8
                + 1.0e6 * max(-order_margin, 0.0)
                + 1.0e6 * max(-min_raw, 0.0),
                "barrier_order_margin_eV": order_margin,
                "min_raw_barrier_eV": min_raw,
            }

        rows: list[dict[str, Any]] = []
        event_detail: list[dict[str, Any]] = []
        for T in self.s.temperatures:
            response = base.simulate_zero_d_rcurve(p, float(T), self.s)
            events = response.pop("events", [])
            rows.append({"T_K": float(T), **response})
            if details:
                event_detail.extend(
                    {"T_K": float(T), **event} for event in events
                )

        pred = pd.DataFrame(rows)
        merged = self.s.targets.merge(pred, on="T_K", how="left")
        incomplete = int(np.sum(~pred.completed.astype(bool)))
        components: dict[str, float] = {
            "completion_loss": 500.0 * incomplete
        }
        components.update(dbtt_target_components(merged))
        objective = float(sum(components.values()))
        if not np.isfinite(objective):
            return {
                "objective": 1.0e12,
                "nonfinite_objective_replaced": True,
            }

        result: dict[str, Any] = {
            "objective": objective,
            **components,
            "barrier_order_margin_eV": order_margin,
            "min_raw_barrier_eV": min_raw,
            "min_peierls_traverse_number": float(
                pred.min_peierls_traverse_number.min()
            ),
            "plateau_temperature_rise": float(
                pred.sort_values("T_K").K_plateau_proxy.iloc[-1]
                - pred.sort_values("T_K").K_plateau_proxy.iloc[0]
            ),
            "max_K_shield_MPa_sqrt_m": float(
                pred.max_K_shield_MPa_sqrt_m.max()
            ),
            "parameters": p,
            "objective_mode": "DBTT_TARGET_AWARE_WITH_FINITE_LOW_SHELF",
        }
        if details:
            result["temperature_detail"] = merged.to_dict(orient="records")
            result["event_detail"] = event_detail
        return result

    def __call__(self, x: np.ndarray) -> float:
        return float(self.evaluate(x, details=False)["objective"])


def dbtt_acceptance(
    target_class: str,
    detail: pd.DataFrame,
    summary: dict[str, Any],
) -> tuple[bool, str]:
    """Broad response-only promotion gate for the corrected DBTT search."""
    if target_class != "DBTT":
        return base.acceptance(target_class, detail, summary)
    if not detail.completed.astype(bool).all():
        return False, "incomplete_zero_d_growth"

    ordered = detail.sort_values("T_K")
    low = ordered.iloc[0]
    high = ordered[ordered.T_K >= 900.0]
    if high.empty:
        high = ordered.tail(max(1, len(ordered) // 2))

    low_init = float(low.K_init_proxy)
    low_plateau = float(low.K_plateau_proxy)
    low_dkr = float(low.delta_KR_proxy)
    high_init = float(high.K_init_proxy.median())
    high_plateau = float(high.K_plateau_proxy.median())
    high_dkr = float(high.delta_KR_proxy.median())
    plateau_rise = high_plateau - low_plateau

    checks = [
        (
            LOW_SHELF_MIN <= low_init <= LOW_SHELF_MAX,
            "DBTT_low_shelf_outside_window",
        ),
        (
            LOW_SHELF_MIN <= low_plateau <= LOW_PLATEAU_MAX,
            "DBTT_low_plateau_outside_window",
        ),
        (low_dkr <= LOW_RCURVE_MAX, "DBTT_low_Rcurve_too_large"),
        (high_init >= HIGH_INIT_MIN, "DBTT_high_initiation_too_low"),
        (high_plateau >= HIGH_PLATEAU_MIN, "DBTT_high_plateau_too_low"),
        (plateau_rise >= PLATEAU_RISE_MIN, "DBTT_transition_too_small"),
        (high_dkr >= HIGH_RCURVE_MIN, "DBTT_high_Rcurve_too_small"),
    ]
    for passed, reason in checks:
        if not passed:
            return False, reason
    return True, "DBTT_target_aware_zeroD_gate_passed"


def _argument_value(flag: str, default: str) -> str:
    try:
        index = sys.argv.index(flag)
    except ValueError:
        return default
    if index + 1 >= len(sys.argv):
        return default
    return sys.argv[index + 1]


def _mark_outputs() -> None:
    outroot = Path(
        _argument_value(
            "--out", "runs/mpz_v9_10_3_dbtt_targeted_global_search_v1"
        )
    )
    class_dir = outroot.resolve() / "DBTT"
    summary_path = class_dir / "unified_global_summary.json"
    config_path = class_dir / "unified_global_config.json"
    if summary_path.exists():
        summary = json.loads(summary_path.read_text())
        summary.update(
            {
                "status": "V9_10_3_DBTT_TARGET_AWARE_GLOBAL_SEARCH_COMPLETE",
                "objective_mode": "DBTT_TARGET_AWARE_WITH_FINITE_LOW_SHELF",
                "parameter_count": len(v102.PARAMETER_NAMES),
                "full_search_space": True,
                "prior_shortlist_used": False,
            }
        )
        summary_path.write_text(json.dumps(summary, indent=2))
    if config_path.exists():
        config = json.loads(config_path.read_text())
        config.update(
            {
                "objective_mode": "DBTT_TARGET_AWARE_WITH_FINITE_LOW_SHELF",
                "low_shelf_window_MPa_sqrt_m": [
                    LOW_SHELF_MIN,
                    LOW_SHELF_MAX,
                ],
                "high_branch_minima": {
                    "K_init_MPa_sqrt_m": HIGH_INIT_MIN,
                    "K_plateau_MPa_sqrt_m": HIGH_PLATEAU_MIN,
                    "delta_KR_MPa_sqrt_m": HIGH_RCURVE_MIN,
                    "plateau_rise_MPa_sqrt_m": PLATEAU_RISE_MIN,
                },
                "mechanism_dominance_constraint": False,
            }
        )
        config_path.write_text(json.dumps(config, indent=2))


def main() -> None:
    # Retain the complete v9.10.2 parameterization and physics while replacing
    # only the DBTT objective and acceptance gate.
    base.PARAMETER_NAMES = v102.PARAMETER_NAMES
    base.BOUNDS = v102.BOUNDS
    base.bounds_list = v102.bounds_list
    base.decode = v102.decode
    base.build_model = v102.build_model
    base.UnifiedObjective = DBTTTargetObjective
    base.acceptance = dbtt_acceptance
    base.main()
    _mark_outputs()


if __name__ == "__main__":
    main()
