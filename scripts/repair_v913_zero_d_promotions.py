#!/usr/bin/env python3
"""Repair zero-D ranking and promotion outputs without rerunning exact cases."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil

import numpy as np
import pandas as pd

import scripts.run_v913_zero_d_large_search as base
from scripts.run_v913_zero_d_large_search_safe import (
    SAFETY_SCHEMA,
    install_safety_patch,
    json_safe,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument(
        "--policy-json",
        type=Path,
        default=Path("mpz_v9_13_zero_d_large_search_policy.json"),
    )
    parser.add_argument("--promote-count", type=int, default=512)
    parser.add_argument(
        "--replace-primary",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Also replace the legacy primary ranking/promotion files after backups.",
    )
    return parser.parse_args()


def atomic_csv(frame: pd.DataFrame, path: Path) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(temporary, index=False)
    temporary.replace(path)


def atomic_json(payload: dict, path: Path) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(json_safe(payload), indent=2, sort_keys=True, allow_nan=False)
        + "\n"
    )
    temporary.replace(path)


def main() -> int:
    args = parse_args()
    run_dir = args.run_dir.resolve()
    ranked_path = run_dir / "zero_d_ranked_candidates.csv"
    contract_path = run_dir / "run_contract.json"
    summary_path = run_dir / "summary.json"
    for required in (ranked_path, contract_path, args.policy_json):
        if not required.is_file():
            raise FileNotFoundError(required)

    install_safety_patch()
    policy = base._load_policy(args.policy_json)
    contract = json.loads(contract_path.read_text())
    objective = contract["contract"]["objective"]
    ranked = pd.read_csv(ranked_path)

    original_gate = ranked["zeroD_gate_pass"].astype(bool).to_numpy()
    original_incomplete_gate = int(
        np.sum(
            original_gate
            & ~pd.to_numeric(ranked["zeroD_complete"], errors="coerce")
            .fillna(0)
            .astype(bool)
            .to_numpy()
        )
    )

    corrected = base._score_frame(
        ranked,
        "zeroD",
        minimum_prominence=float(objective["minimum_prominence"]),
        minimum_drop=float(objective["minimum_post_peak_drop"]),
        maximum_rebound=float(objective["maximum_high_temperature_rebound"]),
        peak_min=float(objective["peak_temperature_min_K"]),
        peak_max=float(objective["peak_temperature_max_K"]),
    )
    corrected = corrected.sort_values(
        ["zeroD_gate_pass", "zeroD_objective", "proxy_objective"],
        ascending=[False, True, True],
        kind="stable",
    ).reset_index(drop=True)
    corrected["zeroD_rank"] = np.arange(1, len(corrected) + 1)

    promoted = base._diverse_selection(corrected, policy, args.promote_count)
    registry_columns = ["candidate_id", *base.ACTIVE_CANDIDATE_PARAMETER_FIELDS]
    registry = promoted[registry_columns].copy()
    registry.insert(1, "zeroD_rank", promoted["zeroD_rank"].to_numpy())
    registry.insert(2, "zeroD_objective", promoted["zeroD_objective"].to_numpy())
    registry.insert(
        3,
        "zeroD_peak_temperature_K",
        promoted["zeroD_peak_temperature_K"].to_numpy(),
    )
    registry.insert(
        4,
        "zeroD_prominence_MPa_sqrt_m",
        promoted["zeroD_two_sided_prominence_MPa_sqrt_m"].to_numpy(),
    )
    registry.insert(
        5,
        "zeroD_high_temperature_rebound_MPa_sqrt_m",
        promoted["zeroD_high_temperature_rebound_MPa_sqrt_m"].to_numpy(),
    )
    registry.insert(6, "promotion_tier", promoted["promotion_tier"].to_numpy())

    corrected_ranked_path = run_dir / "zero_d_ranked_candidates_corrected.csv"
    corrected_metrics_path = run_dir / "promoted_metrics_corrected.csv"
    corrected_registry_path = run_dir / "promoted_registry_corrected.csv"
    atomic_csv(corrected, corrected_ranked_path)
    atomic_csv(promoted, corrected_metrics_path)
    atomic_csv(registry, corrected_registry_path)

    strict_mask = promoted["promotion_tier"].astype(str) == "strict_gate"
    repair_summary = {
        "schema": "v9.13_zero_d_promotion_repair_v1",
        "safety_schema": SAFETY_SCHEMA,
        "source_ranked_table": str(ranked_path),
        "corrected_ranked_table": str(corrected_ranked_path),
        "corrected_promoted_metrics": str(corrected_metrics_path),
        "corrected_promoted_registry": str(corrected_registry_path),
        "exact_evaluated": int(len(corrected)),
        "original_reported_gate_pass_count": int(np.sum(original_gate)),
        "original_incomplete_gate_pass_count": original_incomplete_gate,
        "corrected_complete_gate_pass_count": int(
            corrected["zeroD_gate_pass"].astype(bool).sum()
        ),
        "promoted_count": int(len(promoted)),
        "promoted_complete_count": int(
            promoted["zeroD_complete"].astype(bool).sum()
        ),
        "promoted_strict_gate_count": int(np.sum(strict_mask)),
        "promotion_tier_counts": {
            str(key): int(value)
            for key, value in promoted["promotion_tier"].value_counts().items()
        },
        "all_corrected_strict_passes_promoted": bool(
            set(
                corrected.loc[
                    corrected["zeroD_gate_pass"].astype(bool), "candidate_id"
                ].astype(str)
            ).issubset(set(promoted["candidate_id"].astype(str)))
        ),
        "primary_files_replaced": bool(args.replace_primary),
    }
    atomic_json(repair_summary, run_dir / "promotion_repair_summary.json")

    if args.replace_primary:
        for name in (
            "zero_d_ranked_candidates.csv",
            "promoted_metrics.csv",
            "promoted_registry.csv",
        ):
            source = run_dir / name
            backup = run_dir / f"{name}.pre_completion_promotion_fix"
            if source.exists() and not backup.exists():
                shutil.copy2(source, backup)
        shutil.copy2(corrected_ranked_path, run_dir / "zero_d_ranked_candidates.csv")
        shutil.copy2(corrected_metrics_path, run_dir / "promoted_metrics.csv")
        shutil.copy2(corrected_registry_path, run_dir / "promoted_registry.csv")

        if summary_path.is_file():
            summary = json.loads(summary_path.read_text())
            summary["exact_gate_pass_count_pre_completion_fix"] = summary.get(
                "exact_gate_pass_count"
            )
            summary["exact_gate_pass_count"] = repair_summary[
                "corrected_complete_gate_pass_count"
            ]
            summary["promotion_repair"] = repair_summary
            atomic_json(summary, summary_path)

    print(
        "V913_ZERO_D_PROMOTION_REPAIRED "
        f"exact={len(corrected)} "
        f"strict_complete={repair_summary['corrected_complete_gate_pass_count']} "
        f"promoted={len(promoted)} "
        f"registry={corrected_registry_path}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
