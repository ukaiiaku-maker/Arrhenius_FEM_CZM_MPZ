#!/usr/bin/env python3
"""Analyze staged v9.13 weak-T and ceramic autonomous 1-D R-curves."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import pandas as pd

from arrhenius_fracture.emergent_gnd_contract_v913 import ACTIVE_CANDIDATE_PARAMETER_FIELDS
from arrhenius_fracture.weakT_ceramic_objectives_v913 import (
    add_class_scores,
    diverse_select,
)
from arrhenius_fracture.zero_d_search_v913 import _load_policy


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-root", type=Path, required=True)
    parser.add_argument("--candidate-registry", type=Path, required=True)
    parser.add_argument(
        "--policy-json",
        type=Path,
        default=Path("mpz_v9_13_zero_d_weakT_ceramic_search_policy.json"),
    )
    parser.add_argument("--target-extension-um", type=float, required=True)
    parser.add_argument("--promote-per-class", type=int, default=10)
    parser.add_argument(
        "--temperatures-K",
        nargs="+",
        type=float,
        default=(
            300,
            400,
            500,
            600,
            700,
            800,
            900,
            950,
            1000,
            1050,
            1100,
            1150,
            1200,
            1250,
            1300,
        ),
    )
    parser.add_argument("--out", type=Path, required=True)
    return parser.parse_args()


def _tag(value: float) -> str:
    return f"{float(value):g}".replace(".", "p")


def _events(payload: dict[str, Any]) -> list[dict[str, Any]]:
    events = [dict(row) for row in payload.get("events", [])]
    events.sort(key=lambda row: float(row["cumulative_projected_extension_m"]))
    return events


def checkpoint_K(events: Sequence[dict[str, Any]], extension_um: float) -> float:
    target_m = max(float(extension_um), 0.0) * 1.0e-6
    if not events:
        return float("nan")
    if target_m <= 0.0:
        return float(events[0]["K_MPa_sqrt_m"])
    for event in events:
        if float(event["cumulative_projected_extension_m"]) + 1.0e-15 >= target_m:
            return float(event["K_MPa_sqrt_m"])
    return float("nan")


def extension_average_K(events: Sequence[dict[str, Any]], extension_um: float) -> float:
    target_m = float(extension_um) * 1.0e-6
    if not events or target_m <= 0.0:
        return float("nan")
    x = np.asarray(
        [float(row["cumulative_projected_extension_m"]) for row in events],
        dtype=float,
    )
    y = np.asarray([float(row["K_MPa_sqrt_m"]) for row in events], dtype=float)
    if not np.isfinite(x).all() or not np.isfinite(y).all() or x[-1] + 1.0e-15 < target_m:
        return float("nan")
    keep = np.concatenate(([True], np.diff(x) > 1.0e-15))
    x = x[keep]
    y = y[keep]
    x_aug = np.concatenate(([0.0], x))
    y_aug = np.concatenate(([y[0]], y))
    below = x_aug < target_m
    x_used = x_aug[below]
    y_used = y_aug[below]
    y_target = float(np.interp(target_m, x_aug, y_aug))
    x_used = np.concatenate((x_used, [target_m]))
    y_used = np.concatenate((y_used, [y_target]))
    return float(np.trapezoid(y_used, x_used) / target_m)


def load_cases(
    case_root: Path,
    registry: pd.DataFrame,
    temperatures: Sequence[float],
    target_extension_um: float,
) -> pd.DataFrame:
    records: list[dict[str, Any]] = []
    for candidate_id in registry["candidate_id"].astype(str):
        for temperature in temperatures:
            path = case_root / candidate_id / f"T{_tag(temperature)}K.json"
            if not path.is_file():
                records.append(
                    {
                        "candidate_id": candidate_id,
                        "temperature_K": float(temperature),
                        "status": "missing",
                        "complete": False,
                    }
                )
                continue
            payload = json.loads(path.read_text())
            events = _events(payload)
            achieved = float(payload.get("achieved_projected_extension_um", 0.0))
            complete = bool(
                str(payload.get("status", "")) == "complete"
                and achieved + 1.0e-9 >= target_extension_um
            )
            records.append(
                {
                    "candidate_id": candidate_id,
                    "temperature_K": float(temperature),
                    "status": str(payload.get("status", "")),
                    "complete": complete,
                    "achieved_extension_um": achieved,
                    "K_initial_MPa_sqrt_m": checkpoint_K(events, 0.0),
                    "K_50um_MPa_sqrt_m": checkpoint_K(events, 50.0),
                    "K_250um_MPa_sqrt_m": checkpoint_K(events, 250.0),
                    "K_500um_MPa_sqrt_m": checkpoint_K(events, 500.0),
                    "K_target_MPa_sqrt_m": checkpoint_K(events, target_extension_um),
                    "K_extension_average_MPa_sqrt_m": extension_average_K(
                        events, target_extension_um
                    ),
                    "n_events": len(events),
                    "max_backstress_GPa": float(payload.get("max_backstress_GPa", float("nan"))),
                    "max_tip_radius_um": float(payload.get("max_tip_radius_um", float("nan"))),
                    "min_front_width_um": float(payload.get("min_front_width_um", float("nan"))),
                    "source_case_json": str(path.resolve()),
                }
            )
    return pd.DataFrame(records).sort_values(["candidate_id", "temperature_K"])


def candidate_table(
    cases: pd.DataFrame,
    registry: pd.DataFrame,
    temperatures: Sequence[float],
) -> pd.DataFrame:
    source = registry.set_index(registry["candidate_id"].astype(str), drop=False)
    records: list[dict[str, Any]] = []
    expected = [float(value) for value in temperatures]
    for candidate_id, local in cases.groupby("candidate_id", sort=True):
        local = local.sort_values("temperature_K")
        source_row = source.loc[str(candidate_id)].to_dict()
        record: dict[str, Any] = {
            "candidate_id": str(candidate_id),
            "target_class": str(source_row.get("target_class", "")),
            **{field: source_row[field] for field in ACTIVE_CANDIDATE_PARAMETER_FIELDS},
        }
        complete = bool(
            len(local) == len(expected)
            and local["temperature_K"].astype(float).tolist() == expected
            and local["complete"].fillna(False).astype(bool).all()
            and np.isfinite(local["K_initial_MPa_sqrt_m"].to_numpy(float)).all()
            and np.isfinite(local["K_target_MPa_sqrt_m"].to_numpy(float)).all()
        )
        record["oneD_complete"] = complete
        for row in local.to_dict(orient="records"):
            tag = _tag(float(row["temperature_K"]))
            record[f"oneD_K_initial_T{tag}"] = row.get("K_initial_MPa_sqrt_m", float("nan"))
            record[f"oneD_K_developed_T{tag}"] = row.get("K_target_MPa_sqrt_m", float("nan"))
            record[f"oneD_K_average_T{tag}"] = row.get(
                "K_extension_average_MPa_sqrt_m", float("nan")
            )
            for checkpoint in (50, 250, 500):
                record[f"oneD_K{checkpoint}_T{tag}"] = row.get(
                    f"K_{checkpoint}um_MPa_sqrt_m", float("nan")
                )
        records.append(record)
    frame = pd.DataFrame(records)
    frame = add_class_scores(
        frame,
        temperatures,
        initial_prefix="oneD_K_initial_T",
        developed_prefix="oneD_K_developed_T",
        metric_prefix="oneD_",
    )
    frame["oneD_weakT_gate"] = frame["oneD_weakT_gate"].astype(bool) & frame[
        "oneD_complete"
    ].astype(bool)
    frame["oneD_ceramic_gate"] = frame["oneD_ceramic_gate"].astype(bool) & frame[
        "oneD_complete"
    ].astype(bool)
    return frame


def selected_registry(selection: pd.DataFrame, material_class: str) -> pd.DataFrame:
    fields = ["candidate_id", *ACTIVE_CANDIDATE_PARAMETER_FIELDS]
    result = selection[fields].copy()
    result.insert(1, "target_class", material_class)
    stem = "weakT" if material_class == "weakT_FCC_like" else "ceramic"
    result.insert(2, "oneD_strict_gate_passed", selection[f"oneD_{stem}_gate"].astype(bool).to_numpy())
    result.insert(3, "oneD_class_score", selection[f"oneD_{stem}_score"].to_numpy())
    return result


def main() -> int:
    args = parse_args()
    if args.target_extension_um <= 0.0 or args.promote_per_class < 1:
        raise ValueError("target extension and promotion count must be positive")
    for path in (args.case_root, args.candidate_registry, args.policy_json):
        if not path.exists():
            raise FileNotFoundError(path)
    temperatures = tuple(sorted(set(float(value) for value in args.temperatures_K)))
    registry = pd.read_csv(args.candidate_registry)
    required = {"candidate_id", "target_class", *ACTIVE_CANDIDATE_PARAMETER_FIELDS}
    missing = sorted(required - set(registry.columns))
    if missing:
        raise RuntimeError(f"candidate registry is missing fields: {missing}")
    if registry["candidate_id"].astype(str).duplicated().any():
        raise RuntimeError("candidate registry contains duplicate candidate IDs")
    policy = _load_policy(args.policy_json)
    cases = load_cases(
        args.case_root.resolve(),
        registry,
        temperatures,
        float(args.target_extension_um),
    )
    candidates = candidate_table(cases, registry, temperatures)
    out = args.out.expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)
    cases.to_csv(out / "case_temperature_metrics.csv", index=False)
    candidates.to_csv(out / "candidate_metrics.csv", index=False)

    weak_ranked = candidates[candidates["target_class"] == "weakT_FCC_like"].sort_values(
        ["oneD_weakT_gate", "oneD_weakT_score", "candidate_id"],
        ascending=[False, True, True],
        kind="stable",
    )
    ceramic_ranked = candidates[candidates["target_class"] == "ceramic_like"].sort_values(
        ["oneD_ceramic_gate", "oneD_ceramic_score", "candidate_id"],
        ascending=[False, True, True],
        kind="stable",
    )
    if weak_ranked.empty or ceramic_ranked.empty:
        raise RuntimeError("both target classes must be represented in the candidate registry")
    weak_ranked.to_csv(out / "weakT_ranked.csv", index=False)
    ceramic_ranked.to_csv(out / "ceramic_ranked.csv", index=False)
    weak_selected = diverse_select(
        weak_ranked,
        policy,
        count=min(args.promote_per_class, len(weak_ranked)),
        score_column="oneD_weakT_score",
        gate_column="oneD_weakT_gate",
    )
    ceramic_selected = diverse_select(
        ceramic_ranked,
        policy,
        count=min(args.promote_per_class, len(ceramic_ranked)),
        score_column="oneD_ceramic_score",
        gate_column="oneD_ceramic_gate",
    )
    weak_registry = selected_registry(weak_selected, "weakT_FCC_like")
    ceramic_registry = selected_registry(ceramic_selected, "ceramic_like")
    weak_registry.to_csv(out / "weakT_promoted_registry.csv", index=False)
    ceramic_registry.to_csv(out / "ceramic_promoted_registry.csv", index=False)
    combined = pd.concat([weak_registry, ceramic_registry], ignore_index=True, sort=False)
    combined.to_csv(out / "promoted_registry.csv", index=False)
    summary = {
        "schema": "v9.13_weakT_ceramic_1d_analysis_v1",
        "target_extension_um": float(args.target_extension_um),
        "temperature_grid_K": list(temperatures),
        "candidate_count": int(len(candidates)),
        "complete_candidate_count": int(candidates["oneD_complete"].astype(bool).sum()),
        "weakT_strict_count": int(weak_ranked["oneD_weakT_gate"].astype(bool).sum()),
        "ceramic_strict_count": int(ceramic_ranked["oneD_ceramic_gate"].astype(bool).sum()),
        "weakT_promoted_count": int(len(weak_registry)),
        "ceramic_promoted_count": int(len(ceramic_registry)),
        "promoted_registry": str((out / "promoted_registry.csv").resolve()),
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(
        "V913_WEAKT_CERAMIC_1D_ANALYSIS_COMPLETE "
        f"target_um={args.target_extension_um:g} candidates={len(candidates)} "
        f"weak_strict={summary['weakT_strict_count']} "
        f"ceramic_strict={summary['ceramic_strict_count']} out={out}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
