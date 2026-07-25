#!/usr/bin/env python3
"""Large persistent-site zero-D search for v9.13 DBTT candidates."""
from __future__ import annotations

import argparse
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import platform
import time
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

from arrhenius_fracture.emergent_gnd_campaign_v913 import candidate_from_registry_row
from arrhenius_fracture.emergent_gnd_contract_v913 import (
    ACTIVE_CANDIDATE_PARAMETER_FIELDS,
    PERSISTENT_INACTIVE_REGISTRY_FIELDS,
    candidate_parameter_fingerprint,
    effective_candidate_parameters,
)
from arrhenius_fracture.emergent_gnd_rcurve_v913 import RCurveLoadingMap
from arrhenius_fracture.zero_d_persistent_v913 import (
    ZeroDRunSettings,
    local_peak_metrics,
    run_zero_d_rcurve,
)
from arrhenius_fracture.zero_d_search_v913 import (
    FIXED_ACTIVE_FIELDS,
    VARIABLE_FIELDS,
    _curve_metrics_matrix,
    _load_policy,
    _proxy_response_batch,
    _sample_rows,
    _score_frame,
)
from scripts.run_mpz_v9_13_persistent_top5 import load_physics

_WORKER_PHYSICS: Any = None
_WORKER_LOADING_MAP: RCurveLoadingMap | None = None
_WORKER_TEMPERATURES: tuple[float, ...] = ()
_WORKER_SETTINGS: ZeroDRunSettings | None = None
_WORKER_CHECKPOINT_M: float = 0.0
_WORKER_METRIC_SETTINGS: dict[str, float] = {}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--anchor-registry", type=Path, required=True)
    parser.add_argument("--base-physics-json", type=Path, required=True)
    parser.add_argument("--loading-map", type=Path, required=True)
    parser.add_argument(
        "--policy-json",
        type=Path,
        default=Path("mpz_v9_13_zero_d_large_search_policy.json"),
    )
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--samples", type=int, default=262_144)
    parser.add_argument("--proxy-batch-size", type=int, default=4096)
    parser.add_argument("--exact-count", type=int, default=4096)
    parser.add_argument("--promote-count", type=int, default=512)
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--seed", type=int, default=913_100)
    parser.add_argument(
        "--temperatures-K",
        nargs="+",
        type=float,
        default=(700, 800, 900, 950, 1000, 1050, 1100, 1200, 1300, 1400),
    )
    parser.add_argument("--proxy-extension-um", type=float, default=50.0)
    parser.add_argument("--exact-extension-um", type=float, default=50.0)
    parser.add_argument("--checkpoint-um", type=float, default=50.0)
    parser.add_argument("--load-increment-factor", type=float, default=2.0)
    parser.add_argument("--minimum-prominence", type=float, default=5.0)
    parser.add_argument("--minimum-post-peak-drop", type=float, default=5.0)
    parser.add_argument("--maximum-high-temperature-rebound", type=float, default=3.0)
    parser.add_argument("--peak-temperature-min-K", type=float, default=850.0)
    parser.add_argument("--peak-temperature-max-K", type=float, default=1100.0)
    parser.add_argument("--proxy-target-cleavage-rate-s", type=float, default=1.0e-3)
    parser.add_argument("--proxy-history-events", type=float, default=4.0)
    parser.add_argument("--progress-interval-s", type=float, default=60.0)
    parser.add_argument(
        "--resume", action=argparse.BooleanOptionalAction, default=True
    )
    return parser.parse_args()


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _canonical_sha256(payload: Mapping[str, Any]) -> str:
    text = json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _json_safe(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    return value


def _write_json_atomic(path: Path, payload: Mapping[str, Any]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(
            _json_safe(payload),
            indent=2,
            sort_keys=True,
            allow_nan=False,
        )
        + "\n"
    )
    temporary.replace(path)


def _read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise RuntimeError(f"empty CSV: {path}")
    return rows


def _proxy_batch_path(root: Path, batch_index: int) -> Path:
    return root / f"proxy_batch_{batch_index:05d}.csv.gz"


def _load_anchors(
    rows: Sequence[Mapping[str, str]], policy: Mapping[str, Any]
) -> pd.DataFrame:
    by_id = {str(row["candidate_id"]): row for row in rows}
    selected = []
    missing = []
    for candidate_id in policy["anchor_candidate_ids"]:
        if candidate_id not in by_id:
            missing.append(candidate_id)
            continue
        record = {
            "candidate_id": candidate_id,
            **effective_candidate_parameters(by_id[candidate_id]),
        }
        selected.append(record)
    if missing:
        raise RuntimeError(f"anchor registry is missing policy anchors: {missing}")
    return pd.DataFrame(selected)


def _run_contract(
    args: argparse.Namespace,
    policy: Mapping[str, Any],
    anchor_rows: Sequence[Mapping[str, str]],
) -> dict[str, Any]:
    payload = {
        "schema": "v9.13_persistent_zero_d_large_search_contract_v1",
        "created_at_utc": _utc_now(),
        "runtime": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
        },
        "inputs": {
            "anchor_registry_sha256": _sha256_path(args.anchor_registry),
            "anchor_active_parameter_sha256": candidate_parameter_fingerprint(
                anchor_rows
            ),
            "base_physics_sha256": _sha256_path(args.base_physics_json),
            "loading_map_sha256": _sha256_path(args.loading_map),
            "policy_sha256": _sha256_path(args.policy_json),
        },
        "candidate_contract": {
            "active_fields": list(ACTIVE_CANDIDATE_PARAMETER_FIELDS),
            "variable_fields": list(VARIABLE_FIELDS),
            "fixed_fields": dict(FIXED_ACTIVE_FIELDS),
            "inactive_legacy_fields": list(PERSISTENT_INACTIVE_REGISTRY_FIELDS),
            "finite_source_inventory": False,
            "source_refresh_active": False,
            "explicit_recovery_active": False,
        },
        "search": {
            "samples": args.samples,
            "proxy_batch_size": args.proxy_batch_size,
            "exact_count": args.exact_count,
            "promote_count": args.promote_count,
            "seed": args.seed,
            "temperatures_K": list(args.temperatures_K),
            "proxy_extension_um": args.proxy_extension_um,
            "exact_extension_um": args.exact_extension_um,
            "checkpoint_um": args.checkpoint_um,
            "load_increment_factor": args.load_increment_factor,
            "proxy_target_cleavage_rate_s": args.proxy_target_cleavage_rate_s,
            "proxy_history_events": args.proxy_history_events,
        },
        "objective": {
            "minimum_prominence": args.minimum_prominence,
            "minimum_post_peak_drop": args.minimum_post_peak_drop,
            "maximum_high_temperature_rebound": (
                args.maximum_high_temperature_rebound
            ),
            "peak_temperature_min_K": args.peak_temperature_min_K,
            "peak_temperature_max_K": args.peak_temperature_max_K,
        },
        "policy": policy,
    }
    stable = dict(payload)
    stable.pop("created_at_utc", None)
    return {"sha256": _canonical_sha256(stable), "contract": payload}


def _establish_contract(path: Path, current: Mapping[str, Any]) -> str:
    expected = str(current["sha256"])
    if path.exists():
        previous = json.loads(path.read_text())
        if str(previous.get("sha256", "")) != expected:
            raise RuntimeError(
                "output directory belongs to a different zero-D search contract; "
                "choose a new --out directory"
            )
    else:
        _write_json_atomic(path, current)
    return expected


def _write_progress(path: Path, **updates: Any) -> None:
    payload = {
        "schema": "v9.13_persistent_zero_d_large_search_progress_v1",
        "updated_at_utc": _utc_now(),
        **updates,
    }
    _write_json_atomic(path, payload)


def _worker_initialize(
    physics_json: str,
    loading_map_json: str,
    temperatures: tuple[float, ...],
    exact_extension_m: float,
    load_increment_factor: float,
    checkpoint_m: float,
    metric_settings: dict[str, float],
) -> None:
    global _WORKER_PHYSICS, _WORKER_LOADING_MAP, _WORKER_TEMPERATURES
    global _WORKER_SETTINGS, _WORKER_CHECKPOINT_M, _WORKER_METRIC_SETTINGS
    _WORKER_PHYSICS, _physics_metadata = load_physics(Path(physics_json))
    _WORKER_LOADING_MAP = RCurveLoadingMap.from_dict(
        json.loads(Path(loading_map_json).read_text())
    )
    _WORKER_TEMPERATURES = temperatures
    _WORKER_SETTINGS = ZeroDRunSettings(
        target_projected_extension_m=exact_extension_m,
        load_increment_factor=load_increment_factor,
    )
    _WORKER_CHECKPOINT_M = checkpoint_m
    _WORKER_METRIC_SETTINGS = metric_settings


def _exact_candidate_worker(row: dict[str, Any]) -> dict[str, Any]:
    assert _WORKER_LOADING_MAP is not None
    assert _WORKER_SETTINGS is not None
    candidate = candidate_from_registry_row(row)
    curve = []
    statuses = []
    maxima = {
        "max_backstress_GPa": 0.0,
        "max_tip_radius_um": 0.0,
        "min_front_width_um": float("inf"),
    }
    temperature_rows = []
    for temperature in _WORKER_TEMPERATURES:
        result = run_zero_d_rcurve(
            candidate,
            _WORKER_PHYSICS,
            _WORKER_LOADING_MAP,
            temperature,
            settings=_WORKER_SETTINGS,
        )
        K = result.checkpoint_K(_WORKER_CHECKPOINT_M)
        curve.append(K)
        statuses.append(result.status)
        maxima["max_backstress_GPa"] = max(
            maxima["max_backstress_GPa"], result.max_backstress_Pa * 1.0e-9
        )
        maxima["max_tip_radius_um"] = max(
            maxima["max_tip_radius_um"], result.max_tip_radius_m * 1.0e6
        )
        maxima["min_front_width_um"] = min(
            maxima["min_front_width_um"], result.min_front_width_m * 1.0e6
        )
        temperature_rows.append(
            {
                "T_K": float(temperature),
                "status": result.status,
                "K_checkpoint_MPa_sqrt_m": K,
                "n_events": len(result.events),
                "achieved_extension_um": (
                    result.achieved_projected_extension_m * 1.0e6
                ),
                "max_backstress_GPa": result.max_backstress_Pa * 1.0e-9,
                "max_tip_radius_um": result.max_tip_radius_m * 1.0e6,
                "min_front_width_um": result.min_front_width_m * 1.0e6,
            }
        )
    metrics = local_peak_metrics(
        _WORKER_TEMPERATURES,
        curve,
        desired_min_K=_WORKER_METRIC_SETTINGS["peak_min"],
        desired_max_K=_WORKER_METRIC_SETTINGS["peak_max"],
    )
    record = dict(row)
    for temperature, value in zip(_WORKER_TEMPERATURES, curve):
        tag = f"{float(temperature):g}".replace(".", "p")
        record[f"zeroD_K_T{tag}"] = value
    record.update(
        {
            "zeroD_peak_temperature_K": metrics["peak_temperature_K"],
            "zeroD_peak_value_MPa_sqrt_m": metrics["peak_value"],
            "zeroD_two_sided_prominence_MPa_sqrt_m": (
                metrics["two_sided_prominence"]
            ),
            "zeroD_post_peak_drop_MPa_sqrt_m": metrics["post_peak_drop"],
            "zeroD_high_temperature_rebound_MPa_sqrt_m": (
                metrics["high_temperature_rebound"]
            ),
            "zeroD_peak_internal": int(bool(metrics["peak_internal"])),
            "zeroD_peak_in_desired_window": int(
                bool(metrics["peak_in_desired_window"])
            ),
            "zeroD_complete": int(
                all(status == "complete" for status in statuses)
            ),
            **maxima,
            "temperature_detail": temperature_rows,
        }
    )
    return record


def _exact_case_path(root: Path, candidate_id: str) -> Path:
    return root / f"{candidate_id}.json"


def _run_exact_stage(
    pool: pd.DataFrame,
    args: argparse.Namespace,
    exact_root: Path,
    progress_path: Path,
) -> pd.DataFrame:
    exact_root.mkdir(parents=True, exist_ok=True)
    completed: list[dict[str, Any]] = []
    pending_rows: list[dict[str, Any]] = []
    for row in pool.to_dict(orient="records"):
        path = _exact_case_path(exact_root, str(row["candidate_id"]))
        if args.resume and path.exists():
            completed.append(json.loads(path.read_text()))
        else:
            pending_rows.append(row)
    total = len(pool)
    started = time.monotonic()
    metric_settings = {
        "peak_min": float(args.peak_temperature_min_K),
        "peak_max": float(args.peak_temperature_max_K),
    }
    _write_progress(
        progress_path,
        state="running",
        phase="exact_zero_d",
        completed=len(completed),
        total=total,
        pending=len(pending_rows),
    )
    if pending_rows:
        context = None
        try:
            import multiprocessing

            context = multiprocessing.get_context("fork")
        except (ValueError, RuntimeError):
            context = None
        with ProcessPoolExecutor(
            max_workers=max(int(args.jobs), 1),
            mp_context=context,
            initializer=_worker_initialize,
            initargs=(
                str(args.base_physics_json.resolve()),
                str(args.loading_map.resolve()),
                tuple(float(v) for v in args.temperatures_K),
                float(args.exact_extension_um) * 1.0e-6,
                float(args.load_increment_factor),
                float(args.checkpoint_um) * 1.0e-6,
                metric_settings,
            ),
        ) as executor:
            iterator = iter(pending_rows)
            active: dict[Any, dict[str, Any]] = {}
            for _ in range(min(max(int(args.jobs), 1), len(pending_rows))):
                row = next(iterator, None)
                if row is not None:
                    active[executor.submit(_exact_candidate_worker, row)] = row
            last_heartbeat = time.monotonic()
            while active:
                done, _ = wait(
                    active,
                    timeout=max(float(args.progress_interval_s), 1.0),
                    return_when=FIRST_COMPLETED,
                )
                if not done:
                    now = time.monotonic()
                    if now - last_heartbeat >= args.progress_interval_s:
                        print(
                            "V913_ZERO_D_PROGRESS "
                            f"phase=exact completed={len(completed)}/{total} "
                            f"active={len(active)} elapsed_s={now-started:.1f}",
                            flush=True,
                        )
                        _write_progress(
                            progress_path,
                            state="running",
                            phase="exact_zero_d",
                            completed=len(completed),
                            total=total,
                            active=len(active),
                            elapsed_s=now - started,
                        )
                        last_heartbeat = now
                    continue
                for future in done:
                    source = active.pop(future)
                    try:
                        record = future.result()
                    except Exception as exc:
                        raise RuntimeError(
                            "exact zero-D candidate failed: "
                            f"{source['candidate_id']}"
                        ) from exc
                    path = _exact_case_path(
                        exact_root, str(record["candidate_id"])
                    )
                    _write_json_atomic(path, record)
                    completed.append(record)
                    row = next(iterator, None)
                    if row is not None:
                        active[executor.submit(_exact_candidate_worker, row)] = row
                now = time.monotonic()
                print(
                    "V913_ZERO_D_PROGRESS "
                    f"phase=exact completed={len(completed)}/{total} "
                    f"active={len(active)} elapsed_s={now-started:.1f}",
                    flush=True,
                )
                _write_progress(
                    progress_path,
                    state="running",
                    phase="exact_zero_d",
                    completed=len(completed),
                    total=total,
                    active=len(active),
                    elapsed_s=now - started,
                )
    frame = pd.DataFrame(completed)
    if "temperature_detail" in frame:
        frame = frame.drop(columns=["temperature_detail"])
    return frame


def _normalize_features(
    frame: pd.DataFrame, policy: Mapping[str, Any]
) -> np.ndarray:
    columns = []
    for name in VARIABLE_FIELDS:
        values = pd.to_numeric(frame[name], errors="coerce").to_numpy(float)
        spec = policy["search_dimensions"][name]
        low = float(spec["low"])
        high = float(spec["high"])
        if str(spec["mode"]) == "log10_delta":
            values = np.log10(np.maximum(values, 1.0e-300))
            low = math.log10(low)
            high = math.log10(high)
        columns.append(
            np.clip((values - low) / max(high - low, 1.0e-30), 0.0, 1.0)
        )
    return np.column_stack(columns)


def _diverse_selection(
    ranked: pd.DataFrame,
    policy: Mapping[str, Any],
    count: int,
) -> pd.DataFrame:
    if len(ranked) <= count:
        return ranked.copy()
    pool_count = min(len(ranked), max(count * 12, count))
    pool = ranked.head(pool_count).reset_index(drop=True)
    features = _normalize_features(pool, policy)
    selected = [0]
    min_distance = np.sum(np.square(features - features[0]), axis=1)
    min_distance[0] = -np.inf
    objective = pd.to_numeric(
        pool["zeroD_objective"], errors="coerce"
    ).to_numpy(float)
    scale = (
        np.nanpercentile(objective[np.isfinite(objective)], 90)
        if np.any(np.isfinite(objective))
        else 1.0
    )
    scale = max(float(scale), 1.0e-12)
    for _ in range(1, count):
        quality = np.exp(-np.clip(objective / scale, 0.0, 50.0))
        acquisition = min_distance * (0.25 + 0.75 * quality)
        acquisition[selected] = -np.inf
        index = int(np.argmax(acquisition))
        selected.append(index)
        distance = np.sum(np.square(features - features[index]), axis=1)
        min_distance = np.minimum(min_distance, distance)
    result = pool.iloc[selected].copy()
    result["diversity_rank"] = np.arange(1, len(result) + 1)
    return result


def main() -> int:
    args = parse_args()
    if args.samples < 1 or args.proxy_batch_size < 1:
        raise ValueError("samples and proxy batch size must be positive")
    if args.exact_count < args.promote_count:
        raise ValueError("exact-count must be at least promote-count")
    for path in (
        args.anchor_registry,
        args.base_physics_json,
        args.loading_map,
        args.policy_json,
    ):
        if not path.is_file():
            raise FileNotFoundError(path)
    temperatures = tuple(sorted(set(float(v) for v in args.temperatures_K)))
    if len(temperatures) < 5:
        raise ValueError("at least five temperatures are required")
    args.temperatures_K = temperatures
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    proxy_root = out / "proxy_batches"
    exact_root = out / "exact_cases"
    proxy_root.mkdir(parents=True, exist_ok=True)
    policy = _load_policy(args.policy_json)
    anchor_rows = _read_csv_rows(args.anchor_registry)
    anchors = _load_anchors(anchor_rows, policy)
    physics, _physics_metadata = load_physics(args.base_physics_json)
    loading_map = RCurveLoadingMap.from_dict(
        json.loads(args.loading_map.read_text())
    )
    loading_map.validate()
    if (
        sum(loading_map.projected_advances_m) + 1.0e-15
        < args.exact_extension_um * 1.0e-6
    ):
        raise RuntimeError(
            "loading map does not cover requested exact zero-D extension"
        )
    contract = _run_contract(args, policy, anchor_rows)
    contract_sha = _establish_contract(out / "run_contract.json", contract)
    progress_path = out / "progress.json"

    n_batches = int(math.ceil(args.samples / args.proxy_batch_size))
    proxy_started = time.monotonic()
    for batch_index in range(n_batches):
        start = batch_index * args.proxy_batch_size
        count = min(args.proxy_batch_size, args.samples - start)
        path = _proxy_batch_path(proxy_root, batch_index)
        if args.resume and path.exists():
            continue
        sampled = _sample_rows(
            start=start,
            count=count,
            total_samples=args.samples,
            seed=args.seed,
            anchors=anchors,
            policy=policy,
        )
        evaluated = _proxy_response_batch(
            sampled,
            temperatures,
            physics=physics,
            loading_map=loading_map,
            target_rate_s=args.proxy_target_cleavage_rate_s,
            history_events=args.proxy_history_events,
            target_extension_m=args.proxy_extension_um * 1.0e-6,
        )
        evaluated = _score_frame(
            evaluated,
            "proxy",
            minimum_prominence=args.minimum_prominence,
            minimum_drop=args.minimum_post_peak_drop,
            maximum_rebound=args.maximum_high_temperature_rebound,
            peak_min=args.peak_temperature_min_K,
            peak_max=args.peak_temperature_max_K,
        )
        evaluated["run_contract_sha256"] = contract_sha
        temporary = path.with_suffix(path.suffix + ".tmp")
        evaluated.to_csv(temporary, index=False, compression="gzip")
        temporary.replace(path)
        elapsed = time.monotonic() - proxy_started
        completed_samples = start + count
        print(
            "V913_ZERO_D_PROGRESS "
            f"phase=proxy completed={completed_samples}/{args.samples} "
            f"batch={batch_index+1}/{n_batches} elapsed_s={elapsed:.1f}",
            flush=True,
        )
        _write_progress(
            progress_path,
            state="running",
            phase="proxy",
            completed_samples=completed_samples,
            total_samples=args.samples,
            completed_batches=batch_index + 1,
            total_batches=n_batches,
            elapsed_s=elapsed,
        )

    proxy_pool: pd.DataFrame | None = None
    keep_count = min(
        args.samples, max(args.exact_count * 8, args.exact_count)
    )
    proxy_pass_count = 0
    for batch_index in range(n_batches):
        batch = pd.read_csv(_proxy_batch_path(proxy_root, batch_index))
        proxy_pass_count += int(batch["proxy_gate_pass"].astype(bool).sum())
        candidate = batch.nsmallest(
            min(keep_count, len(batch)), "proxy_objective"
        )
        proxy_pool = (
            candidate
            if proxy_pool is None
            else pd.concat([proxy_pool, candidate], ignore_index=True)
        )
        proxy_pool = proxy_pool.nsmallest(
            min(keep_count, len(proxy_pool)), "proxy_objective"
        )
    assert proxy_pool is not None
    exact_input = proxy_pool.nsmallest(
        args.exact_count, "proxy_objective"
    ).copy()
    exact_input.to_csv(out / "proxy_exact_input.csv", index=False)
    _write_json_atomic(
        out / "proxy_summary.json",
        {
            "schema": "v9.13_persistent_zero_d_proxy_summary_v1",
            "samples": args.samples,
            "proxy_gate_pass_count": proxy_pass_count,
            "exact_count": len(exact_input),
            "contract_sha256": contract_sha,
        },
    )

    exact = _run_exact_stage(exact_input, args, exact_root, progress_path)
    exact = _score_frame(
        exact,
        "zeroD",
        minimum_prominence=args.minimum_prominence,
        minimum_drop=args.minimum_post_peak_drop,
        maximum_rebound=args.maximum_high_temperature_rebound,
        peak_min=args.peak_temperature_min_K,
        peak_max=args.peak_temperature_max_K,
    )
    exact = exact.sort_values(
        ["zeroD_gate_pass", "zeroD_objective", "proxy_objective"],
        ascending=[False, True, True],
        kind="stable",
    ).reset_index(drop=True)
    exact["zeroD_rank"] = np.arange(1, len(exact) + 1)
    exact.to_csv(out / "zero_d_ranked_candidates.csv", index=False)

    promoted = _diverse_selection(exact, policy, args.promote_count)
    registry_columns = ["candidate_id", *ACTIVE_CANDIDATE_PARAMETER_FIELDS]
    registry = promoted[registry_columns].copy()
    registry.insert(1, "zeroD_rank", promoted["zeroD_rank"].to_numpy())
    registry.insert(
        2, "zeroD_objective", promoted["zeroD_objective"].to_numpy()
    )
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
    registry.to_csv(out / "promoted_registry.csv", index=False)
    promoted.to_csv(out / "promoted_metrics.csv", index=False)

    summary = {
        "schema": "v9.13_persistent_zero_d_large_search_summary_v1",
        "status": "complete",
        "contract_sha256": contract_sha,
        "samples": args.samples,
        "proxy_gate_pass_count": proxy_pass_count,
        "exact_evaluated": len(exact),
        "exact_gate_pass_count": int(
            exact["zeroD_gate_pass"].astype(bool).sum()
        ),
        "promoted_count": len(promoted),
        "promoted_registry": str((out / "promoted_registry.csv").resolve()),
        "full_one_dimensional_validation_required": True,
        "completed_at_utc": _utc_now(),
    }
    _write_json_atomic(out / "summary.json", summary)
    _write_progress(
        progress_path,
        state="complete",
        phase="complete",
        **summary,
    )
    print(
        "V913_ZERO_D_SEARCH_COMPLETE "
        f"samples={args.samples} exact={len(exact)} "
        f"promoted={len(promoted)} out={out}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
