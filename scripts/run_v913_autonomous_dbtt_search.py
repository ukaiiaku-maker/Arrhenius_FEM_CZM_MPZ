#!/usr/bin/env python3
"""Search v9.12 candidates with the calibrated autonomous v9.13 R-curve model.

The response objective is not redefined here.  This driver applies the
directional-DBTT and genuine-peak trajectory metrics already used by the
v9.12 active-learning workflow to the autonomous K-versus-extension response.
Candidate parameter rows are read-only.
"""

from __future__ import annotations

import argparse
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
import csv
from datetime import datetime, timezone
import hashlib
import json
import multiprocessing
import os
from pathlib import Path
import platform
import time
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

from arrhenius_fracture.emergent_gnd_campaign_v913 import (
    candidate_from_registry_row,
)
from arrhenius_fracture.emergent_gnd_contract_v913 import (
    ACTIVE_CANDIDATE_PARAMETER_FIELDS,
    PERSISTENT_INACTIVE_REGISTRY_FIELDS,
    candidate_feature_record,
    candidate_parameter_fingerprint,
    effective_candidate_parameters,
)
from arrhenius_fracture.emergent_gnd_rcurve_v913 import (
    RCurveLoadingMap,
    run_autonomous_rcurve,
)
from scripts.augment_mpz_v9_12_directional_peak_targets import (
    add_directional_peak_classifications,
    add_trajectory_metrics,
)
from scripts.run_mpz_v9_13_persistent_top5 import load_physics


_WORKER_PHYSICS: Any = None
_WORKER_LOADING_MAP: RCurveLoadingMap | None = None
_WORKER_SETTINGS: dict[str, float] = {}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-registry", type=Path, required=True)
    parser.add_argument(
        "--base-physics-json",
        type=Path,
        default=Path("mpz_v9_13_v10222_transfer_common_physics.json"),
    )
    parser.add_argument(
        "--loading-map",
        type=Path,
        default=Path(
            "runs/v9_13_v10222_rcurve_targets_v1/v10_2_22_rcurve_loading_map.json"
        ),
    )
    parser.add_argument(
        "--policy-json",
        type=Path,
        default=Path("mpz_v9_12_targeted_local_search_policy.json"),
    )
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--families", nargs="*", default=("peak",))
    parser.add_argument("--candidate-ids", nargs="*", default=())
    parser.add_argument(
        "--per-parent",
        type=int,
        default=128,
        help="Nested Sobol prefix per parent; zero selects every matching row.",
    )
    parser.add_argument(
        "--parent-offset",
        type=int,
        default=0,
        help="Start index within each parent family for a later acquisition wave.",
    )
    parser.add_argument(
        "--temperatures",
        nargs="+",
        type=float,
        default=(700.0, 800.0, 900.0, 1000.0, 1100.0, 1200.0),
    )
    parser.add_argument("--checkpoint-um", type=float, default=25.0)
    parser.add_argument("--target-extension-um", type=float, default=25.0)
    parser.add_argument("--translation-action-exponent", type=float, default=0.95)
    parser.add_argument("--max-hazard-increment", type=float, default=0.05)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument(
        "--progress-interval-s",
        type=float,
        default=60.0,
        help=(
            "Maximum time between aggregate progress heartbeats while cases "
            "are running."
        ),
    )
    parser.add_argument("--promote-count", type=int, default=48)
    parser.add_argument("--low-max-K", type=float, default=700.0)
    parser.add_argument("--high-min-K", type=float, default=1000.0)
    parser.add_argument("--peak-min-K", type=float, default=800.0)
    parser.add_argument("--peak-max-K", type=float, default=1000.0)
    parser.add_argument("--direction-threshold", type=float, default=5.0)
    parser.add_argument("--peak-threshold", type=float, default=1.0)
    parser.add_argument(
        "--resume",
        action=argparse.BooleanOptionalAction,
        default=True,
    )
    return parser.parse_args()


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise RuntimeError(f"empty CSV: {path}")
    return rows


def _write_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    if not rows:
        return
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=fields,
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _write_json_atomic(path: Path, payload: Mapping[str, Any]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    temporary.replace(path)


def _progress_payload(
    *,
    state: str,
    phase: str,
    started_at_utc: str,
    elapsed_s: float,
    completed_cases: int,
    resumed_cases: int,
    total_cases: int,
    jobs: int,
    contract_sha256: str,
    last_case: Mapping[str, Any] | None = None,
    error: str | None = None,
) -> dict[str, Any]:
    new_cases = max(int(completed_cases) - int(resumed_cases), 0)
    rate_per_s = new_cases / elapsed_s if elapsed_s > 0.0 else 0.0
    remaining = max(int(total_cases) - int(completed_cases), 0)
    eta_s = remaining / rate_per_s if rate_per_s > 0.0 else None
    active_workers = min(max(int(jobs), 0), remaining)
    payload: dict[str, Any] = {
        "schema": "v9.13_autonomous_dbtt_progress_v1",
        "state": str(state),
        "phase": str(phase),
        "started_at_utc": str(started_at_utc),
        "updated_at_utc": _utc_now(),
        "run_contract_sha256": str(contract_sha256),
        "completed_cases": int(completed_cases),
        "resumed_cases": int(resumed_cases),
        "newly_completed_cases": new_cases,
        "remaining_cases": remaining,
        "total_cases": int(total_cases),
        "progress_fraction": (
            float(completed_cases) / float(total_cases)
            if total_cases
            else 1.0
        ),
        "active_workers_upper_bound": active_workers,
        "elapsed_s": float(elapsed_s),
        "new_cases_per_hour": 3600.0 * rate_per_s,
        "eta_s": eta_s,
        "last_case": dict(last_case) if last_case is not None else None,
        "error": error,
    }
    return payload


def _progress_line(payload: Mapping[str, Any]) -> str:
    eta = payload.get("eta_s")
    eta_text = "unknown" if eta is None else f"{float(eta):.1f}"
    return (
        "V913_DBTT_PROGRESS "
        f"state={payload['state']} phase={payload['phase']} "
        f"completed={payload['completed_cases']}/{payload['total_cases']} "
        f"resumed={payload['resumed_cases']} "
        f"remaining={payload['remaining_cases']} "
        f"active_lte={payload['active_workers_upper_bound']} "
        f"elapsed_s={float(payload['elapsed_s']):.1f} "
        f"new_cases_per_hour={float(payload['new_cases_per_hour']):.3f} "
        f"eta_s={eta_text}"
    )


def _temperature_tag(temperature_K: float) -> str:
    rounded = round(float(temperature_K))
    if np.isclose(temperature_K, rounded):
        return str(int(rounded))
    return f"{float(temperature_K):g}".replace(".", "p")


def _case_path(root: Path, candidate_id: str, temperature_K: float) -> Path:
    return root / candidate_id / f"T{_temperature_tag(temperature_K)}K.json"


def _candidate_fingerprint(rows: Sequence[Mapping[str, str]]) -> str:
    return candidate_parameter_fingerprint(rows)


def _sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _canonical_sha256(payload: Mapping[str, Any]) -> str:
    text = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _run_contract(
    args: argparse.Namespace,
    *,
    all_rows: Sequence[Mapping[str, str]],
    selected_rows: Sequence[Mapping[str, str]],
    temperatures: Sequence[float],
) -> dict[str, Any]:
    repo_root = Path(__file__).resolve().parents[1]
    source_paths = tuple(
        sorted((repo_root / "arrhenius_fracture").glob("*.py"))
    ) + (
        Path(__file__).resolve(),
        repo_root / "scripts/augment_mpz_v9_12_directional_peak_targets.py",
        repo_root / "scripts/run_mpz_v9_13_persistent_top5.py",
    )
    payload = {
        "schema": "v9.13_autonomous_dbtt_run_contract_v2",
        "model_source_sha256": {
            str(path.relative_to(repo_root)): _sha256_path(path)
            for path in source_paths
        },
        "runtime": {
            "python": platform.python_version(),
            "numpy": np.__version__,
        },
        "inputs": {
            "candidate_registry_sha256": _sha256_path(args.candidate_registry),
            "candidate_parameter_sha256": _candidate_fingerprint(all_rows),
            "selected_candidate_parameter_sha256": _candidate_fingerprint(
                selected_rows
            ),
            "base_physics_sha256": _sha256_path(args.base_physics_json),
            "loading_map_sha256": _sha256_path(args.loading_map),
            "policy_sha256": _sha256_path(args.policy_json),
        },
        "candidate_contract": {
            "active_fields": list(ACTIVE_CANDIDATE_PARAMETER_FIELDS),
            "inactive_legacy_fields": list(PERSISTENT_INACTIVE_REGISTRY_FIELDS),
            "candidate_parameters_refit": False,
            "source_refresh_active": False,
            "explicit_recovery_active": False,
        },
        "selection": {
            "families": list(args.families),
            "candidate_ids": list(args.candidate_ids),
            "per_parent": int(args.per_parent),
            "parent_offset": int(args.parent_offset),
            "selected_candidates": len(selected_rows),
        },
        "simulation": {
            "temperatures_K": list(temperatures),
            "checkpoint_um": float(args.checkpoint_um),
            "target_extension_um": float(args.target_extension_um),
            "translation_mode": "hazard_coupled",
            "translation_action_exponent": float(
                args.translation_action_exponent
            ),
            "max_hazard_increment": float(args.max_hazard_increment),
        },
        "objective": {
            "response": f"K({args.checkpoint_um:g} um, T)",
            "low_max_K": float(args.low_max_K),
            "high_min_K": float(args.high_min_K),
            "peak_min_K": float(args.peak_min_K),
            "peak_max_K": float(args.peak_max_K),
            "direction_threshold": float(args.direction_threshold),
            "peak_threshold": float(args.peak_threshold),
        },
    }
    return {
        "sha256": _canonical_sha256(payload),
        "contract": payload,
    }


def _establish_run_contract(
    path: Path,
    current: Mapping[str, Any],
    *,
    case_root: Path,
) -> str:
    expected = str(current["sha256"])
    if path.exists():
        previous = json.loads(path.read_text())
        actual = str(previous.get("sha256", ""))
        if actual != expected or previous.get("contract") != current.get("contract"):
            raise RuntimeError(
                "output directory belongs to a different autonomous-search "
                "contract; choose a new --out directory instead of mixing cases "
                f"(existing={actual or '<missing>'}, requested={expected})"
            )
    else:
        existing_cases = list(case_root.glob("*/T*K.json"))
        if existing_cases:
            raise RuntimeError(
                "case JSON files exist without run_contract.json; refusing "
                "unsafe legacy resume. Choose a new --out directory."
            )
        path.write_text(json.dumps(current, indent=2) + "\n")
    return expected


def _validate_resumed_payload(
    payload: Mapping[str, Any],
    *,
    candidate_id: str,
    temperature_K: float,
    contract_sha256: str,
    loading_map_seed: int,
) -> None:
    errors: list[str] = []
    if str(payload.get("run_contract_sha256", "")) != contract_sha256:
        errors.append("run_contract_sha256")
    if str(payload.get("candidate_id", "")) != candidate_id:
        errors.append("candidate_id")
    try:
        stored_temperature = float(payload.get("temperature_K"))
    except (TypeError, ValueError):
        stored_temperature = float("nan")
    if not np.isclose(stored_temperature, temperature_K, rtol=0.0, atol=1.0e-12):
        errors.append("temperature_K")
    try:
        stored_seed = int(payload.get("seed"))
    except (TypeError, ValueError):
        stored_seed = -1
    if stored_seed != int(loading_map_seed):
        errors.append("seed")
    if errors:
        raise RuntimeError(
            f"unsafe resume payload for {candidate_id} at {temperature_K:g} K; "
            f"invalid fields: {errors}"
        )


def _select_rows(
    rows: Sequence[dict[str, str]],
    *,
    families: Sequence[str],
    candidate_ids: Sequence[str],
    per_parent: int,
    parent_offset: int,
) -> list[dict[str, str]]:
    by_id = {row["candidate_id"]: row for row in rows}
    if len(by_id) != len(rows):
        raise RuntimeError("candidate_id must be unique")
    if candidate_ids:
        missing = sorted(set(candidate_ids) - set(by_id))
        if missing:
            raise KeyError(f"candidate IDs are absent from registry: {missing}")
        return [by_id[candidate_id] for candidate_id in candidate_ids]

    family_filter = {str(value).lower() for value in families}
    filtered = [
        row
        for row in rows
        if not family_filter
        or str(row.get("campaign_parent_family", "")).lower() in family_filter
    ]
    if not filtered:
        raise RuntimeError("family selection produced no candidates")
    if per_parent == 0:
        return filtered
    if per_parent < 0 or parent_offset < 0:
        raise ValueError("--per-parent and --parent-offset must be nonnegative")

    grouped: dict[str, list[dict[str, str]]] = {}
    for row in filtered:
        parent = str(row.get("campaign_parent_id") or "all")
        grouped.setdefault(parent, []).append(row)
    selected: list[dict[str, str]] = []
    for parent in sorted(grouped):
        ordered = sorted(grouped[parent], key=lambda row: row["candidate_id"])
        selected.extend(ordered[parent_offset : parent_offset + per_parent])
    if not selected:
        raise RuntimeError("parent offset and count selected no candidates")
    return selected


def _initialize_worker(
    physics_path: str,
    loading_map_path: str,
    settings: Mapping[str, float],
) -> None:
    global _WORKER_PHYSICS, _WORKER_LOADING_MAP, _WORKER_SETTINGS
    _WORKER_PHYSICS, _ = load_physics(Path(physics_path))
    _WORKER_LOADING_MAP = RCurveLoadingMap.from_dict(
        json.loads(Path(loading_map_path).read_text())
    )
    _WORKER_SETTINGS = dict(settings)


def _checkpoint_from_payload(payload: Mapping[str, Any], extension_um: float) -> float:
    events = list(payload.get("events", []))
    if not events:
        return float("nan")
    target_m = float(extension_um) * 1.0e-6
    for event in events:
        if float(event["cumulative_projected_extension_m"]) >= target_m:
            return float(event["K_MPa_sqrt_m"])
    return float(events[-1]["K_MPa_sqrt_m"])


def _case_record(
    payload: Mapping[str, Any],
    *,
    checkpoint_um: float,
) -> dict[str, Any]:
    return {
        "candidate_id": str(payload["candidate_id"]),
        "temperature_K": float(payload["temperature_K"]),
        "status": str(payload["status"]),
        "seed": int(payload["seed"]),
        "checkpoint_um": float(checkpoint_um),
        "K_checkpoint_MPa_sqrt_m": _checkpoint_from_payload(payload, checkpoint_um),
        "K_first_MPa_sqrt_m": float(payload["K_first_MPa_sqrt_m"]),
        "K_10um_MPa_sqrt_m": float(payload["K_10um_MPa_sqrt_m"]),
        "K_25um_MPa_sqrt_m": float(payload["K_25um_MPa_sqrt_m"]),
        "K_50um_MPa_sqrt_m": float(payload["K_50um_MPa_sqrt_m"]),
        "achieved_projected_extension_um": float(
            payload["achieved_projected_extension_um"]
        ),
        "max_backstress_GPa": float(payload["max_backstress_GPa"]),
        "min_front_width_um": float(payload["min_front_width_um"]),
        "max_tip_radius_um": float(payload["max_tip_radius_um"]),
        "max_source_multiplicity": float(payload["max_source_multiplicity"]),
        "n_events": len(payload.get("events", [])),
    }


def _run_case(
    row: Mapping[str, str],
    temperature_K: float,
) -> tuple[dict[str, Any], float]:
    if _WORKER_LOADING_MAP is None or _WORKER_PHYSICS is None:
        raise RuntimeError("worker was not initialized")
    started = time.perf_counter()
    print(
        "V913_DBTT_CASE_START "
        f"candidate={row['candidate_id']} T={float(temperature_K):g} "
        f"worker_pid={os.getpid()}",
        flush=True,
    )
    candidate = candidate_from_registry_row(row)
    result = run_autonomous_rcurve(
        candidate,
        _WORKER_PHYSICS,
        _WORKER_LOADING_MAP,
        float(temperature_K),
        target_projected_extension_m=(_WORKER_SETTINGS["target_extension_um"] * 1.0e-6),
        max_hazard_increment=_WORKER_SETTINGS["max_hazard_increment"],
        translation_mode="hazard_coupled",
        translation_action_exponent=(_WORKER_SETTINGS["translation_action_exponent"]),
    )
    return result.as_dict(), time.perf_counter() - started


def _objective_tables(
    selected_rows: Sequence[dict[str, str]],
    case_records: Sequence[Mapping[str, Any]],
    *,
    temperatures: Sequence[float],
    checkpoint_um: float,
    policy_path: Path,
    low_max_K: float,
    high_min_K: float,
    peak_min_K: float,
    peak_max_K: float,
    direction_threshold: float,
    peak_threshold: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    # Read the policy here because it remains the immutable source of the
    # established objective/generator provenance.  Surrogate inputs, however,
    # are every active v9.13 candidate constant, not the legacy generator's
    # inactive recovery/source-refresh coordinates.
    json.loads(policy_path.read_text())
    by_case = {
        (str(row["candidate_id"]), float(row["temperature_K"])): row
        for row in case_records
    }
    response_prefix = f"y__K{_temperature_tag(checkpoint_um)}_1d_"
    records: list[dict[str, Any]] = []
    for candidate in selected_rows:
        candidate_id = candidate["candidate_id"]
        record: dict[str, Any] = {
            "candidate_id": candidate_id,
            "campaign_parent_id": candidate.get("campaign_parent_id", ""),
            "campaign_parent_family": candidate.get("campaign_parent_family", ""),
        }
        record.update(candidate_feature_record(candidate))
        complete = True
        for temperature in temperatures:
            case = by_case.get((candidate_id, float(temperature)))
            value = float("nan")
            if case is None:
                complete = False
            else:
                value = float(case["K_checkpoint_MPa_sqrt_m"])
                complete &= str(case["status"]) == "complete" and np.isfinite(value)
            record[f"{response_prefix}T{_temperature_tag(temperature)}K"] = value
        record["y__complete_1d"] = bool(complete)
        records.append(record)

    table = pd.DataFrame(records)
    add_trajectory_metrics(
        table,
        prefix=response_prefix,
        out_prefix="y__",
        low_max_K=low_max_K,
        high_min_K=high_min_K,
        peak_min_K=peak_min_K,
        peak_max_K=peak_max_K,
    )
    add_directional_peak_classifications(
        table,
        out_prefix="y__",
        peak_min_K=peak_min_K,
        direction_threshold=direction_threshold,
        peak_threshold=peak_threshold,
    )
    table["objective_response"] = f"K_at_{checkpoint_um:g}_um"
    table["objective_provenance"] = (
        "v9_12_directional_peak_targets_on_v9_13_autonomous_R_curve"
    )

    ranking = table.copy()
    ranking["_peak_like"] = ranking["y__peak_like_1d"].fillna(False).astype(bool)
    ranking = ranking.sort_values(
        by=[
            "_peak_like",
            "y__peak_prominence",
            "y__peak_drop",
            "y__peak_rise",
            "candidate_id",
        ],
        ascending=[False, False, False, False, True],
        kind="stable",
    ).drop(columns=["_peak_like"])
    ranking.insert(0, "search_rank", np.arange(1, len(ranking) + 1))
    return table, ranking


def _event_rows(payloads: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for payload in payloads:
        for event in payload.get("events", []):
            rows.append(
                {
                    "candidate_id": payload["candidate_id"],
                    "temperature_K": payload["temperature_K"],
                    "status": payload["status"],
                    **event,
                }
            )
    return rows


def main() -> int:
    args = parse_args()
    if args.checkpoint_um <= 0.0 or args.target_extension_um <= 0.0:
        raise ValueError("checkpoint and target extension must be positive")
    if args.checkpoint_um > args.target_extension_um:
        raise ValueError("checkpoint extension must not exceed target extension")
    if int(args.jobs) < 1:
        raise ValueError("--jobs must be at least one")
    if not np.isfinite(args.progress_interval_s) or args.progress_interval_s <= 0.0:
        raise ValueError("--progress-interval-s must be positive and finite")
    args.out.mkdir(parents=True, exist_ok=True)
    case_root = args.out / "cases"
    case_root.mkdir(parents=True, exist_ok=True)
    all_rows = _read_csv(args.candidate_registry)
    for row in all_rows:
        effective_candidate_parameters(row)
    selected_rows = _select_rows(
        all_rows,
        families=args.families,
        candidate_ids=args.candidate_ids,
        per_parent=args.per_parent,
        parent_offset=args.parent_offset,
    )
    temperatures = sorted({float(value) for value in args.temperatures})
    if not temperatures or not all(np.isfinite(temperatures)):
        raise ValueError("temperatures must be nonempty and finite")
    required_windows = {
        "low": any(value <= args.low_max_K for value in temperatures),
        "peak": any(
            args.peak_min_K <= value <= args.peak_max_K
            for value in temperatures
        ),
        "high": any(value >= args.high_min_K for value in temperatures),
    }
    missing_windows = [
        name for name, present in required_windows.items() if not present
    ]
    if missing_windows:
        raise ValueError(
            "temperature grid does not cover objective windows: "
            f"{missing_windows}"
        )
    loading_map = RCurveLoadingMap.from_dict(
        json.loads(args.loading_map.read_text())
    )
    contract = _run_contract(
        args,
        all_rows=all_rows,
        selected_rows=selected_rows,
        temperatures=temperatures,
    )
    contract_sha256 = _establish_run_contract(
        args.out / "run_contract.json",
        contract,
        case_root=case_root,
    )
    settings = {
        "target_extension_um": float(args.target_extension_um),
        "max_hazard_increment": float(args.max_hazard_increment),
        "translation_action_exponent": float(args.translation_action_exponent),
    }

    payloads: dict[tuple[str, float], dict[str, Any]] = {}
    pending: list[tuple[dict[str, str], float]] = []
    for row in selected_rows:
        for temperature in temperatures:
            path = _case_path(case_root, row["candidate_id"], temperature)
            if args.resume and path.exists():
                payload = json.loads(path.read_text())
                _validate_resumed_payload(
                    payload,
                    candidate_id=row["candidate_id"],
                    temperature_K=temperature,
                    contract_sha256=contract_sha256,
                    loading_map_seed=loading_map.seed,
                )
                payloads[(row["candidate_id"], temperature)] = payload
            else:
                pending.append((row, temperature))

    progress_path = args.out / "autonomous_dbtt_progress.json"
    total_cases = len(selected_rows) * len(temperatures)
    resumed_cases = len(payloads)
    search_started_at_utc = _utc_now()
    search_started = time.perf_counter()
    last_case: dict[str, Any] | None = None

    def report_progress(
        *,
        state: str = "running",
        phase: str = "cases",
        emit: bool = True,
        error: str | None = None,
    ) -> dict[str, Any]:
        progress = _progress_payload(
            state=state,
            phase=phase,
            started_at_utc=search_started_at_utc,
            elapsed_s=time.perf_counter() - search_started,
            completed_cases=len(payloads),
            resumed_cases=resumed_cases,
            total_cases=total_cases,
            jobs=int(args.jobs),
            contract_sha256=contract_sha256,
            last_case=last_case,
            error=error,
        )
        _write_json_atomic(progress_path, progress)
        if emit:
            print(_progress_line(progress), flush=True)
        return progress

    print(
        "V913_DBTT_SEARCH_START "
        f"selected={len(selected_rows)} cases={total_cases} "
        f"resumed={len(payloads)} pending={len(pending)} jobs={args.jobs}",
        flush=True,
    )
    report_progress()
    initializer_args = (
        str(args.base_physics_json),
        str(args.loading_map),
        settings,
    )

    def accept(payload: dict[str, Any], wall_s: float) -> None:
        nonlocal last_case
        key = (str(payload["candidate_id"]), float(payload["temperature_K"]))
        payload["run_contract_sha256"] = contract_sha256
        _validate_resumed_payload(
            payload,
            candidate_id=key[0],
            temperature_K=key[1],
            contract_sha256=contract_sha256,
            loading_map_seed=loading_map.seed,
        )
        payloads[key] = payload
        path = _case_path(case_root, *key)
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_suffix(path.suffix + ".tmp")
        temporary.write_text(json.dumps(payload, indent=2) + "\n")
        temporary.replace(path)
        last_case = {
            "candidate_id": key[0],
            "temperature_K": key[1],
            "status": str(payload["status"]),
            "K_checkpoint_MPa_sqrt_m": _checkpoint_from_payload(
                payload,
                args.checkpoint_um,
            ),
            "case_wall_s": float(wall_s),
        }
        records = [
            _case_record(item, checkpoint_um=args.checkpoint_um)
            for item in payloads.values()
        ]
        _write_csv(
            args.out / "case_results_checkpoint.csv",
            sorted(
                records,
                key=lambda row: (row["candidate_id"], row["temperature_K"]),
            ),
        )
        print(
            "V913_DBTT_CASE_COMPLETE "
            f"candidate={key[0]} T={key[1]:g} status={payload['status']} "
            f"K={_checkpoint_from_payload(payload, args.checkpoint_um):.8g} "
            f"wall_s={wall_s:.3f} complete={len(payloads)}/"
            f"{total_cases}",
            flush=True,
        )
        report_progress(emit=False)

    if pending:
        context = multiprocessing.get_context("spawn")
        with ProcessPoolExecutor(
            max_workers=max(int(args.jobs), 1),
            mp_context=context,
            initializer=_initialize_worker,
            initargs=initializer_args,
        ) as executor:
            futures = {
                executor.submit(_run_case, row, temperature): (
                    row["candidate_id"],
                    temperature,
                )
                for row, temperature in pending
            }
            remaining = set(futures)
            heartbeat_deadline = (
                time.perf_counter() + float(args.progress_interval_s)
            )
            while remaining:
                timeout_s = max(heartbeat_deadline - time.perf_counter(), 0.0)
                done, remaining = wait(
                    remaining,
                    timeout=timeout_s,
                    return_when=FIRST_COMPLETED,
                )
                if not done:
                    report_progress()
                    heartbeat_deadline = (
                        time.perf_counter() + float(args.progress_interval_s)
                    )
                    continue
                for future in done:
                    candidate_id, temperature = futures[future]
                    try:
                        payload, wall_s = future.result()
                    except Exception as exc:
                        message = (
                            f"{type(exc).__name__}: {exc}; "
                            f"candidate={candidate_id}; T={temperature:g}"
                        )
                        report_progress(
                            state="failed",
                            phase="cases",
                            error=message,
                        )
                        raise
                    accept(payload, wall_s)
                if time.perf_counter() >= heartbeat_deadline:
                    report_progress()
                    heartbeat_deadline = (
                        time.perf_counter() + float(args.progress_interval_s)
                    )

    report_progress(phase="postprocess")
    case_records = [
        _case_record(payload, checkpoint_um=args.checkpoint_um)
        for payload in payloads.values()
    ]
    table, ranking = _objective_tables(
        selected_rows,
        case_records,
        temperatures=temperatures,
        checkpoint_um=args.checkpoint_um,
        policy_path=args.policy_json,
        low_max_K=args.low_max_K,
        high_min_K=args.high_min_K,
        peak_min_K=args.peak_min_K,
        peak_max_K=args.peak_max_K,
        direction_threshold=args.direction_threshold,
        peak_threshold=args.peak_threshold,
    )
    table.to_csv(args.out / "autonomous_dbtt_training_table.csv", index=False)
    ranking.to_csv(args.out / "ranked_candidates.csv", index=False)
    pool_records = []
    for candidate in all_rows:
        pool_records.append(
            {
                "candidate_id": candidate["candidate_id"],
                "campaign_parent_id": candidate.get("campaign_parent_id", ""),
                "campaign_parent_family": candidate.get("campaign_parent_family", ""),
                **candidate_feature_record(candidate),
            }
        )
    pd.DataFrame(pool_records).to_csv(
        args.out / "candidate_pool_features.csv",
        index=False,
    )
    _write_csv(
        args.out / "R_curve_events.csv",
        _event_rows(list(payloads.values())),
    )

    promote_count = min(max(int(args.promote_count), 0), len(ranking))
    promoted_ids = ranking.head(promote_count)["candidate_id"].astype(str)
    rank_columns = [
        "candidate_id",
        "search_rank",
        "y__peak_like_1d",
        "y__peak_prominence",
        "y__peak_temperature_K",
        "y__peak_rise",
        "y__peak_drop",
    ]
    promoted = pd.DataFrame(all_rows)
    promoted = promoted[promoted["candidate_id"].isin(promoted_ids)].merge(
        ranking[rank_columns],
        on="candidate_id",
        how="left",
        validate="one_to_one",
    )
    promoted = promoted.sort_values("search_rank")
    promoted.to_csv(args.out / "promoted_registry.csv", index=False)

    expected_cases = total_cases
    complete_grid = (
        len(payloads) == expected_cases
        and all(str(payload.get("status")) == "complete" for payload in payloads.values())
    )
    status_counts: dict[str, int] = {}
    for payload in payloads.values():
        status = str(payload.get("status", "missing"))
        status_counts[status] = status_counts.get(status, 0) + 1
    manifest = {
        "schema": "v9.13_autonomous_dbtt_candidate_search",
        "run_contract_sha256": contract_sha256,
        "candidate_parameters_refit": False,
        "candidate_registry": str(args.candidate_registry.resolve()),
        "candidate_registry_file_sha256": _sha256_path(args.candidate_registry),
        "candidate_parameter_sha256": _candidate_fingerprint(all_rows),
        "selected_candidate_parameter_sha256": _candidate_fingerprint(selected_rows),
        "active_candidate_parameter_fields": list(
            ACTIVE_CANDIDATE_PARAMETER_FIELDS
        ),
        "inactive_legacy_registry_fields": list(
            PERSISTENT_INACTIVE_REGISTRY_FIELDS
        ),
        "selected_candidates": len(selected_rows),
        "temperatures_K": temperatures,
        "complete_grid": complete_grid,
        "completed_cases": len(payloads),
        "case_status_counts": status_counts,
        "families": list(args.families),
        "per_parent": int(args.per_parent),
        "parent_offset": int(args.parent_offset),
        "checkpoint_um": float(args.checkpoint_um),
        "target_extension_um": float(args.target_extension_um),
        "translation_action_exponent": float(args.translation_action_exponent),
        "max_hazard_increment": float(args.max_hazard_increment),
        "objective": {
            "implementation": ("scripts/augment_mpz_v9_12_directional_peak_targets.py"),
            "response": f"K({args.checkpoint_um:g} um, T)",
            "low_max_K": float(args.low_max_K),
            "high_min_K": float(args.high_min_K),
            "peak_min_K": float(args.peak_min_K),
            "peak_max_K": float(args.peak_max_K),
            "direction_threshold": float(args.direction_threshold),
            "peak_threshold": float(args.peak_threshold),
        },
        "peak_like_count": int(
            ranking["y__peak_like_1d"].fillna(False).astype(bool).sum()
        ),
        "promoted_candidates": promote_count,
        "host_cpu_count": os.cpu_count(),
    }
    (args.out / "search_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )
    print(
        "V913_DBTT_SEARCH_COMPLETE "
        f"complete_grid={str(complete_grid).lower()} "
        f"peak_like={manifest['peak_like_count']} promoted={promote_count} "
        f"out={args.out}",
        flush=True,
    )
    report_progress(
        state="complete" if complete_grid else "failed",
        phase="complete",
    )
    return 0 if complete_grid else 2


if __name__ == "__main__":
    raise SystemExit(main())
