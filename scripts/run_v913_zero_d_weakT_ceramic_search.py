#!/usr/bin/env python3
"""Run a class-specific persistent-site zero-D search for weak-T and ceramic rows.

This is a new Sobol population, not a reranking of the DBTT/peak exact subset.  It
uses the current v9.13 persistent-site physics and permits temperature coefficients
and source densities excluded by the peak-oriented search policy.
"""
from __future__ import annotations

import argparse
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
import multiprocessing
from pathlib import Path
import time
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

from arrhenius_fracture.emergent_gnd_campaign_v913 import candidate_from_registry_row
from arrhenius_fracture.emergent_gnd_contract_v913 import (
    ACTIVE_CANDIDATE_PARAMETER_FIELDS,
    candidate_parameter_fingerprint,
    effective_candidate_parameters,
)
from arrhenius_fracture.emergent_gnd_rcurve_v913 import RCurveLoadingMap
from arrhenius_fracture.weakT_ceramic_objectives_v913 import (
    add_class_scores,
    diverse_select,
)
from arrhenius_fracture.zero_d_persistent_v913 import ZeroDRunSettings, run_zero_d_rcurve
from arrhenius_fracture.zero_d_search_v913 import (
    _load_policy,
    _neutral_cleavage_K,
    _proxy_response_batch,
    _sample_rows,
)
from scripts.run_mpz_v9_13_persistent_top5 import load_physics

_WORKER_PHYSICS: Any = None
_WORKER_LOADING_MAP: RCurveLoadingMap | None = None
_WORKER_TEMPERATURES: tuple[float, ...] = ()
_WORKER_SETTINGS: ZeroDRunSettings | None = None
_WORKER_CHECKPOINT_M = 0.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--anchor-registry", type=Path, required=True)
    parser.add_argument("--base-physics-json", type=Path, required=True)
    parser.add_argument("--loading-map", type=Path, required=True)
    parser.add_argument(
        "--policy-json",
        type=Path,
        default=Path("mpz_v9_13_zero_d_weakT_ceramic_search_policy.json"),
    )
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--samples", type=int, default=262_144)
    parser.add_argument("--proxy-batch-size", type=int, default=4096)
    parser.add_argument("--exact-per-class", type=int, default=512)
    parser.add_argument("--promote-per-class", type=int, default=48)
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--seed", type=int, default=913_260)
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
    parser.add_argument("--proxy-extension-um", type=float, default=100.0)
    parser.add_argument("--exact-extension-um", type=float, default=100.0)
    parser.add_argument("--checkpoint-um", type=float, default=100.0)
    parser.add_argument("--load-increment-factor", type=float, default=2.0)
    parser.add_argument("--proxy-target-cleavage-rate-s", type=float, default=1.0e-3)
    parser.add_argument("--proxy-history-events", type=float, default=8.0)
    parser.add_argument("--progress-interval-s", type=float, default=60.0)
    parser.add_argument(
        "--resume", action=argparse.BooleanOptionalAction, default=True
    )
    return parser.parse_args()


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _tag(value: float) -> str:
    return f"{float(value):g}".replace(".", "p")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _canonical_sha(payload: Mapping[str, Any]) -> str:
    text = json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _json_safe(value: Any) -> Any:
    if isinstance(value, np.generic):
        return _json_safe(value.item())
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    return value


def _write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(_json_safe(payload), indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    temporary.replace(path)


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise RuntimeError(f"empty CSV: {path}")
    return rows


def _physical_surface_gate(frame: pd.DataFrame, temperatures: Sequence[float]) -> np.ndarray:
    gate = np.ones(len(frame), dtype=bool)
    Tref = pd.to_numeric(frame["Tref_K"], errors="coerce").to_numpy(float)
    for prefix in ("cleave", "emit"):
        G00 = pd.to_numeric(frame[f"{prefix}_G00_eV"], errors="coerce").to_numpy(float)
        gT = pd.to_numeric(frame[f"{prefix}_gT_eV_per_K"], errors="coerce").to_numpy(float)
        sig0 = pd.to_numeric(frame[f"{prefix}_sigc0_GPa"], errors="coerce").to_numpy(float)
        sT = pd.to_numeric(frame[f"{prefix}_sT_GPa_per_K"], errors="coerce").to_numpy(float)
        for temperature in (min(temperatures), max(temperatures)):
            gate &= G00 + gT * (float(temperature) - Tref) >= 0.05
            gate &= sig0 + sT * (float(temperature) - Tref) >= 0.10
    return gate


def _score_proxy(
    sampled: pd.DataFrame,
    temperatures: Sequence[float],
    *,
    physics: Any,
    loading_map: RCurveLoadingMap,
    args: argparse.Namespace,
) -> pd.DataFrame:
    evaluated = _proxy_response_batch(
        sampled,
        temperatures,
        physics=physics,
        loading_map=loading_map,
        target_rate_s=float(args.proxy_target_cleavage_rate_s),
        history_events=float(args.proxy_history_events),
        target_extension_m=float(args.proxy_extension_um) * 1.0e-6,
    )
    for temperature in temperatures:
        evaluated[f"proxy_K_initial_T{_tag(temperature)}"] = _neutral_cleavage_K(
            evaluated,
            float(temperature),
            physics=physics,
            target_rate_s=float(args.proxy_target_cleavage_rate_s),
        )
    evaluated = add_class_scores(
        evaluated,
        temperatures,
        initial_prefix="proxy_K_initial_T",
        developed_prefix="proxy_K_T",
        metric_prefix="proxy_",
    )
    physical = _physical_surface_gate(evaluated, temperatures)
    evaluated["proxy_physical_surface_gate"] = physical
    for material_class in ("weakT", "ceramic"):
        gate = f"proxy_{material_class}_gate"
        score = f"proxy_{material_class}_score"
        evaluated[gate] = evaluated[gate].astype(bool) & physical
        evaluated[score] = pd.to_numeric(evaluated[score], errors="coerce") + np.where(
            physical, 0.0, 1000.0
        )
    return evaluated


def _worker_initialize(
    physics_path: str,
    loading_map_path: str,
    temperatures: tuple[float, ...],
    exact_extension_m: float,
    checkpoint_m: float,
    load_increment_factor: float,
) -> None:
    global _WORKER_PHYSICS, _WORKER_LOADING_MAP, _WORKER_TEMPERATURES
    global _WORKER_SETTINGS, _WORKER_CHECKPOINT_M
    _WORKER_PHYSICS, _ = load_physics(Path(physics_path))
    _WORKER_LOADING_MAP = RCurveLoadingMap.from_dict(
        json.loads(Path(loading_map_path).read_text())
    )
    _WORKER_TEMPERATURES = temperatures
    _WORKER_SETTINGS = ZeroDRunSettings(
        target_projected_extension_m=exact_extension_m,
        load_increment_factor=load_increment_factor,
    )
    _WORKER_CHECKPOINT_M = checkpoint_m


def _exact_worker(row: dict[str, Any]) -> dict[str, Any]:
    if _WORKER_LOADING_MAP is None or _WORKER_SETTINGS is None:
        raise RuntimeError("zero-D worker was not initialized")
    candidate = candidate_from_registry_row(row)
    record = dict(row)
    complete = True
    maxima = {
        "zeroD_max_backstress_GPa": 0.0,
        "zeroD_max_tip_radius_um": 0.0,
        "zeroD_min_front_width_um": float("inf"),
    }
    temperature_detail: list[dict[str, Any]] = []
    for temperature in _WORKER_TEMPERATURES:
        result = run_zero_d_rcurve(
            candidate,
            _WORKER_PHYSICS,
            _WORKER_LOADING_MAP,
            float(temperature),
            settings=_WORKER_SETTINGS,
        )
        initial = result.events[0].K_MPa_sqrt_m if result.events else float("nan")
        developed = result.checkpoint_K(_WORKER_CHECKPOINT_M)
        tag = _tag(temperature)
        record[f"zeroD_K_initial_T{tag}"] = initial
        record[f"zeroD_K_developed_T{tag}"] = developed
        record[f"zeroD_status_T{tag}"] = result.status
        complete &= result.status == "complete" and math.isfinite(initial) and math.isfinite(developed)
        maxima["zeroD_max_backstress_GPa"] = max(
            maxima["zeroD_max_backstress_GPa"], result.max_backstress_Pa * 1.0e-9
        )
        maxima["zeroD_max_tip_radius_um"] = max(
            maxima["zeroD_max_tip_radius_um"], result.max_tip_radius_m * 1.0e6
        )
        maxima["zeroD_min_front_width_um"] = min(
            maxima["zeroD_min_front_width_um"], result.min_front_width_m * 1.0e6
        )
        temperature_detail.append(
            {
                "temperature_K": float(temperature),
                "status": result.status,
                "K_initial_MPa_sqrt_m": initial,
                "K_developed_MPa_sqrt_m": developed,
                "achieved_extension_um": result.achieved_projected_extension_m * 1.0e6,
                "n_events": len(result.events),
            }
        )
    record["zeroD_complete"] = bool(complete)
    record.update(maxima)
    record["temperature_detail"] = temperature_detail
    return record


def _run_exact(
    pool: pd.DataFrame,
    args: argparse.Namespace,
    exact_root: Path,
) -> pd.DataFrame:
    exact_root.mkdir(parents=True, exist_ok=True)
    completed: list[dict[str, Any]] = []
    pending: list[dict[str, Any]] = []
    for row in pool.to_dict(orient="records"):
        path = exact_root / f"{row['candidate_id']}.json"
        if args.resume and path.is_file():
            completed.append(json.loads(path.read_text()))
        else:
            pending.append(row)
    started = time.monotonic()
    context = None
    try:
        context = multiprocessing.get_context("fork")
    except (ValueError, RuntimeError):
        context = None
    if pending:
        with ProcessPoolExecutor(
            max_workers=max(int(args.jobs), 1),
            mp_context=context,
            initializer=_worker_initialize,
            initargs=(
                str(args.base_physics_json.resolve()),
                str(args.loading_map.resolve()),
                tuple(float(value) for value in args.temperatures_K),
                float(args.exact_extension_um) * 1.0e-6,
                float(args.checkpoint_um) * 1.0e-6,
                float(args.load_increment_factor),
            ),
        ) as executor:
            iterator = iter(pending)
            active: dict[Any, dict[str, Any]] = {}
            for _ in range(min(max(int(args.jobs), 1), len(pending))):
                row = next(iterator, None)
                if row is not None:
                    active[executor.submit(_exact_worker, row)] = row
            last_heartbeat = time.monotonic()
            while active:
                done, _ = wait(
                    active,
                    timeout=max(float(args.progress_interval_s), 1.0),
                    return_when=FIRST_COMPLETED,
                )
                if not done:
                    now = time.monotonic()
                    if now - last_heartbeat >= float(args.progress_interval_s):
                        print(
                            "V913_WEAKT_CERAMIC_ZERO_D_PROGRESS "
                            f"completed={len(completed)}/{len(pool)} active={len(active)} "
                            f"elapsed_s={now-started:.1f}",
                            flush=True,
                        )
                        last_heartbeat = now
                    continue
                for future in done:
                    source = active.pop(future)
                    try:
                        record = future.result()
                    except Exception as exc:
                        raise RuntimeError(
                            f"exact zero-D candidate failed: {source['candidate_id']}"
                        ) from exc
                    _write_json(exact_root / f"{record['candidate_id']}.json", record)
                    completed.append(record)
                    row = next(iterator, None)
                    if row is not None:
                        active[executor.submit(_exact_worker, row)] = row
                now = time.monotonic()
                print(
                    "V913_WEAKT_CERAMIC_ZERO_D_PROGRESS "
                    f"completed={len(completed)}/{len(pool)} active={len(active)} "
                    f"elapsed_s={now-started:.1f}",
                    flush=True,
                )
    frame = pd.DataFrame(completed)
    if "temperature_detail" in frame:
        frame = frame.drop(columns=["temperature_detail"])
    return frame


def _registry(selection: pd.DataFrame, material_class: str) -> pd.DataFrame:
    fields = ["candidate_id", *ACTIVE_CANDIDATE_PARAMETER_FIELDS]
    result = selection[fields].copy()
    result.insert(1, "target_class", material_class)
    gate = f"zeroD_{'weakT' if material_class == 'weakT_FCC_like' else 'ceramic'}_gate"
    score = f"zeroD_{'weakT' if material_class == 'weakT_FCC_like' else 'ceramic'}_score"
    result.insert(2, "zeroD_strict_gate_passed", selection[gate].astype(bool).to_numpy())
    result.insert(3, "zeroD_class_score", pd.to_numeric(selection[score], errors="coerce").to_numpy())
    return result


def main() -> int:
    args = parse_args()
    if args.samples < 1 or args.proxy_batch_size < 1:
        raise ValueError("sample and batch counts must be positive")
    if args.exact_per_class < args.promote_per_class:
        raise ValueError("exact-per-class must be at least promote-per-class")
    for path in (
        args.anchor_registry,
        args.base_physics_json,
        args.loading_map,
        args.policy_json,
    ):
        if not path.is_file():
            raise FileNotFoundError(path)
    temperatures = tuple(sorted(set(float(value) for value in args.temperatures_K)))
    if min(temperatures) > 300.0 or max(temperatures) < 1200.0:
        raise ValueError("temperature grid must include 300 K and at least 1200 K")
    args.temperatures_K = temperatures
    out = args.out.expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)
    proxy_root = out / "proxy_batches"
    exact_root = out / "exact_cases"
    proxy_root.mkdir(parents=True, exist_ok=True)

    policy = _load_policy(args.policy_json)
    anchor_rows = _read_csv(args.anchor_registry)
    by_id = {str(row["candidate_id"]): row for row in anchor_rows}
    anchors = []
    for candidate_id in policy["anchor_candidate_ids"]:
        if candidate_id not in by_id:
            raise RuntimeError(f"anchor registry lacks policy candidate {candidate_id}")
        anchors.append({"candidate_id": candidate_id, **effective_candidate_parameters(by_id[candidate_id])})
    anchors_frame = pd.DataFrame(anchors)
    physics, _ = load_physics(args.base_physics_json)
    loading_map = RCurveLoadingMap.from_dict(json.loads(args.loading_map.read_text()))
    loading_map.validate()
    coverage_um = sum(loading_map.projected_advances_m) * 1.0e6
    if coverage_um + 1.0e-9 < max(args.proxy_extension_um, args.exact_extension_um):
        raise RuntimeError(
            f"loading map coverage {coverage_um:.6g} um is below requested zero-D extension"
        )

    contract_payload = {
        "schema": "v9.13_weakT_ceramic_zero_d_search_contract_v1",
        "anchor_registry_sha256": _sha256(args.anchor_registry),
        "anchor_parameter_fingerprint_sha256": candidate_parameter_fingerprint(anchor_rows),
        "physics_sha256": _sha256(args.base_physics_json),
        "loading_map_sha256": _sha256(args.loading_map),
        "policy_sha256": _sha256(args.policy_json),
        "samples": int(args.samples),
        "proxy_batch_size": int(args.proxy_batch_size),
        "exact_per_class": int(args.exact_per_class),
        "promote_per_class": int(args.promote_per_class),
        "seed": int(args.seed),
        "temperatures_K": list(temperatures),
        "proxy_extension_um": float(args.proxy_extension_um),
        "exact_extension_um": float(args.exact_extension_um),
        "checkpoint_um": float(args.checkpoint_um),
    }
    contract = {"sha256": _canonical_sha(contract_payload), "contract": contract_payload}
    contract_path = out / "run_contract.json"
    if contract_path.is_file():
        previous = json.loads(contract_path.read_text())
        if previous != contract:
            raise RuntimeError("output directory belongs to a different search contract")
    else:
        _write_json(contract_path, contract)

    n_batches = int(math.ceil(args.samples / args.proxy_batch_size))
    proxy_started = time.monotonic()
    for batch_index in range(n_batches):
        start = batch_index * args.proxy_batch_size
        count = min(args.proxy_batch_size, args.samples - start)
        path = proxy_root / f"proxy_batch_{batch_index:05d}.csv.gz"
        if args.resume and path.is_file():
            continue
        sampled = _sample_rows(
            start=start,
            count=count,
            total_samples=args.samples,
            seed=args.seed,
            anchors=anchors_frame,
            policy=policy,
        )
        evaluated = _score_proxy(
            sampled,
            temperatures,
            physics=physics,
            loading_map=loading_map,
            args=args,
        )
        evaluated["run_contract_sha256"] = contract["sha256"]
        temporary = path.with_suffix(path.suffix + ".tmp")
        evaluated.to_csv(temporary, index=False, compression="gzip")
        temporary.replace(path)
        print(
            "V913_WEAKT_CERAMIC_PROXY_PROGRESS "
            f"completed={start+count}/{args.samples} batch={batch_index+1}/{n_batches} "
            f"elapsed_s={time.monotonic()-proxy_started:.1f}",
            flush=True,
        )

    keep = max(args.exact_per_class * 24, args.exact_per_class)
    weak_pool: pd.DataFrame | None = None
    ceramic_pool: pd.DataFrame | None = None
    for batch_index in range(n_batches):
        batch = pd.read_csv(proxy_root / f"proxy_batch_{batch_index:05d}.csv.gz")
        weak = batch.sort_values(
            ["proxy_weakT_gate", "proxy_weakT_score"],
            ascending=[False, True],
            kind="stable",
        ).head(min(keep, len(batch)))
        ceramic = batch.sort_values(
            ["proxy_ceramic_gate", "proxy_ceramic_score"],
            ascending=[False, True],
            kind="stable",
        ).head(min(keep, len(batch)))
        weak_pool = weak if weak_pool is None else pd.concat([weak_pool, weak], ignore_index=True)
        ceramic_pool = ceramic if ceramic_pool is None else pd.concat([ceramic_pool, ceramic], ignore_index=True)
        weak_pool = weak_pool.sort_values(
            ["proxy_weakT_gate", "proxy_weakT_score"],
            ascending=[False, True],
            kind="stable",
        ).head(keep)
        ceramic_pool = ceramic_pool.sort_values(
            ["proxy_ceramic_gate", "proxy_ceramic_score"],
            ascending=[False, True],
            kind="stable",
        ).head(keep)
    assert weak_pool is not None and ceramic_pool is not None
    weak_exact = diverse_select(
        weak_pool,
        policy,
        count=args.exact_per_class,
        score_column="proxy_weakT_score",
        gate_column="proxy_weakT_gate",
    )
    ceramic_candidates = ceramic_pool[
        ~ceramic_pool["candidate_id"].astype(str).isin(set(weak_exact["candidate_id"].astype(str)))
    ]
    ceramic_exact = diverse_select(
        ceramic_candidates,
        policy,
        count=args.exact_per_class,
        score_column="proxy_ceramic_score",
        gate_column="proxy_ceramic_gate",
    )
    weak_exact.to_csv(out / "proxy_weakT_exact_input.csv", index=False)
    ceramic_exact.to_csv(out / "proxy_ceramic_exact_input.csv", index=False)
    exact_input = pd.concat([weak_exact, ceramic_exact], ignore_index=True, sort=False)
    exact_input = exact_input.drop_duplicates("candidate_id", keep="first")
    exact_input.to_csv(out / "proxy_combined_exact_input.csv", index=False)

    exact = _run_exact(exact_input, args, exact_root)
    exact = add_class_scores(
        exact,
        temperatures,
        initial_prefix="zeroD_K_initial_T",
        developed_prefix="zeroD_K_developed_T",
        metric_prefix="zeroD_",
    )
    exact["zeroD_weakT_gate"] = exact["zeroD_weakT_gate"].astype(bool) & exact["zeroD_complete"].astype(bool)
    exact["zeroD_ceramic_gate"] = exact["zeroD_ceramic_gate"].astype(bool) & exact["zeroD_complete"].astype(bool)
    weak_ranked = exact.sort_values(
        ["zeroD_weakT_gate", "zeroD_weakT_score", "candidate_id"],
        ascending=[False, True, True],
        kind="stable",
    ).reset_index(drop=True)
    ceramic_ranked = exact.sort_values(
        ["zeroD_ceramic_gate", "zeroD_ceramic_score", "candidate_id"],
        ascending=[False, True, True],
        kind="stable",
    ).reset_index(drop=True)
    weak_ranked.to_csv(out / "zeroD_weakT_ranked.csv", index=False)
    ceramic_ranked.to_csv(out / "zeroD_ceramic_ranked.csv", index=False)

    weak_promoted = diverse_select(
        weak_ranked,
        policy,
        count=args.promote_per_class,
        score_column="zeroD_weakT_score",
        gate_column="zeroD_weakT_gate",
    )
    ceramic_available = ceramic_ranked[
        ~ceramic_ranked["candidate_id"].astype(str).isin(set(weak_promoted["candidate_id"].astype(str)))
    ]
    ceramic_promoted = diverse_select(
        ceramic_available,
        policy,
        count=args.promote_per_class,
        score_column="zeroD_ceramic_score",
        gate_column="zeroD_ceramic_gate",
    )
    weak_registry = _registry(weak_promoted, "weakT_FCC_like")
    ceramic_registry = _registry(ceramic_promoted, "ceramic_like")
    weak_registry.to_csv(out / "weakT_promoted_registry.csv", index=False)
    ceramic_registry.to_csv(out / "ceramic_promoted_registry.csv", index=False)
    combined = pd.concat([weak_registry, ceramic_registry], ignore_index=True, sort=False)
    combined.to_csv(out / "combined_promoted_registry.csv", index=False)

    summary = {
        "schema": "v9.13_weakT_ceramic_zero_d_search_summary_v1",
        "status": "complete",
        "completed_at_utc": _utc_now(),
        "contract_sha256": contract["sha256"],
        "samples": int(args.samples),
        "exact_evaluated": int(len(exact)),
        "weakT_strict_exact_count": int(exact["zeroD_weakT_gate"].astype(bool).sum()),
        "ceramic_strict_exact_count": int(exact["zeroD_ceramic_gate"].astype(bool).sum()),
        "weakT_promoted_count": int(len(weak_registry)),
        "ceramic_promoted_count": int(len(ceramic_registry)),
        "combined_promoted_registry": str((out / "combined_promoted_registry.csv").resolve()),
        "one_dimensional_validation_required": True,
    }
    _write_json(out / "summary.json", summary)
    print(
        "V913_WEAKT_CERAMIC_ZERO_D_SEARCH_COMPLETE "
        f"samples={args.samples} exact={len(exact)} promoted={len(combined)} out={out}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
