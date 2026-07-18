#!/usr/bin/env python3
"""Crash-safe first-passage search wrapper for v9.10.4.5.

The v9.10.4 search objective intentionally returns a compact result when a
candidate is incomplete or rejected before the transition calculation. The
result-export loop, however, assumes that every detailed evaluation contains
``parameters``, ``temperature_detail``, and ``event_detail``. This wrapper
stabilizes that schema before delegating to the existing search driver.

It also checkpoints the completed differential-evolution population before
Powell refinement and detailed export. If a restart later crashes, the same
command reloads that DE population and resumes after the expensive global
search rather than repeating all generations.

No objective, constitutive rate, timestep, transition gate, or optimization
algorithm is changed.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
import sys
from typing import Any

import numpy as np
from scipy.optimize import OptimizeResult

import optimize_mpz_v9_10_4_narrow_dbtt as _base


_ORIGINAL_EVALUATE = _base.NarrowDBTTObjective.evaluate
_ORIGINAL_DIFFERENTIAL_EVOLUTION = _base.differential_evolution
_PENDING_RESTARTS: list[int] | None = None
_PENDING_POSITION = 0


def stabilize_detailed_result(
    result: dict[str, Any],
    x: np.ndarray,
    *,
    details: bool,
) -> dict[str, Any]:
    """Return a stable export schema for every detailed objective result."""
    stable = dict(result)
    if not details:
        return stable

    parameters = stable.get("parameters")
    if not isinstance(parameters, dict):
        try:
            parameters = _base.decode(np.asarray(x, dtype=float))
        except Exception:
            parameters = {}
    stable["parameters"] = parameters

    temperature_detail = stable.get("temperature_detail")
    stable["temperature_detail"] = (
        list(temperature_detail) if isinstance(temperature_detail, list) else []
    )

    event_detail = stable.get("event_detail")
    stable["event_detail"] = list(event_detail) if isinstance(event_detail, list) else []

    if float(stable.get("completion_loss", 0.0)) > 0.0:
        stable.setdefault("evaluation_status", "INCOMPLETE_CANDIDATE")
    elif bool(stable.get("invalid_parameter_vector", False)):
        stable.setdefault("evaluation_status", "INVALID_PARAMETER_VECTOR")
    elif "transition_shelf_ratio" not in stable:
        stable.setdefault("evaluation_status", "EARLY_REJECTED_CANDIDATE")
    else:
        stable.setdefault("evaluation_status", "COMPLETE_CANDIDATE")
    return stable


def _stable_evaluate(self, x: np.ndarray, *, details: bool = False) -> dict[str, Any]:
    result = _ORIGINAL_EVALUATE(self, x, details=details)
    return stabilize_detailed_result(result, np.asarray(x, dtype=float), details=details)


def _argument_value(flag: str, default: str) -> str:
    try:
        index = sys.argv.index(flag)
    except ValueError:
        return default
    if index + 1 >= len(sys.argv):
        return default
    return str(sys.argv[index + 1])


def _output_directory() -> Path:
    return Path(
        _argument_value(
            "--out",
            "runs/mpz_v9_10_4_narrow_dbtt_first_passage_v1",
        )
    ).resolve()


def _restart_count() -> int:
    return int(_argument_value("--restarts", "4"))


def _pending_restart_indices() -> list[int]:
    out = _output_directory()
    checkpoints = out / "checkpoints"
    pending: list[int] = []
    for restart in range(_restart_count()):
        completed = checkpoints / f"restart_{restart:03d}.json"
        if completed.exists():
            try:
                payload = json.loads(completed.read_text())
            except Exception:
                payload = {}
            if payload.get("status") == "COMPLETE":
                continue
        pending.append(restart)
    return pending


def _next_restart_index() -> int:
    global _PENDING_RESTARTS, _PENDING_POSITION
    if _PENDING_RESTARTS is None:
        _PENDING_RESTARTS = _pending_restart_indices()
    if _PENDING_POSITION >= len(_PENDING_RESTARTS):
        raise RuntimeError("differential-evolution call count exceeds pending restarts")
    restart = int(_PENDING_RESTARTS[_PENDING_POSITION])
    _PENDING_POSITION += 1
    return restart


def _de_state_path(restart: int) -> Path:
    path = _output_directory() / "checkpoints" / f"restart_{restart:03d}_de_state.npz"
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def _save_de_result(path: Path, result: OptimizeResult) -> None:
    """Atomically save the DE state needed by the existing export loop."""
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("wb") as handle:
        np.savez_compressed(
            handle,
            x=np.asarray(result.x, dtype=float),
            fun=np.asarray(float(result.fun), dtype=float),
            population=np.asarray(result.population, dtype=float),
            population_energies=np.asarray(result.population_energies, dtype=float),
            nit=np.asarray(int(getattr(result, "nit", -1)), dtype=int),
            nfev=np.asarray(int(getattr(result, "nfev", -1)), dtype=int),
        )
    os.replace(temporary, path)


def _load_de_result(path: Path) -> OptimizeResult:
    with np.load(path, allow_pickle=False) as data:
        return OptimizeResult(
            x=np.asarray(data["x"], dtype=float),
            fun=float(np.asarray(data["fun"]).item()),
            population=np.asarray(data["population"], dtype=float),
            population_energies=np.asarray(data["population_energies"], dtype=float),
            nit=int(np.asarray(data["nit"]).item()),
            nfev=int(np.asarray(data["nfev"]).item()),
            success=True,
            message="resumed completed differential-evolution state",
        )


def _checkpointing_differential_evolution(*args: Any, **kwargs: Any) -> OptimizeResult:
    restart = _next_restart_index()
    state_path = _de_state_path(restart)
    if state_path.exists():
        result = _load_de_result(state_path)
        print(
            f"[resume-de] restart={restart} objective={float(result.fun):.6g} "
            f"state={state_path}",
            flush=True,
        )
        return result

    print(f"[de-start] restart={restart}", flush=True)
    result = _ORIGINAL_DIFFERENTIAL_EVOLUTION(*args, **kwargs)
    _save_de_result(state_path, result)
    print(
        f"[de-checkpoint] restart={restart} objective={float(result.fun):.6g} "
        f"state={state_path}",
        flush=True,
    )
    return result


_base.NarrowDBTTObjective.evaluate = _stable_evaluate
_base.differential_evolution = _checkpointing_differential_evolution


if __name__ == "__main__":
    _base.main()
