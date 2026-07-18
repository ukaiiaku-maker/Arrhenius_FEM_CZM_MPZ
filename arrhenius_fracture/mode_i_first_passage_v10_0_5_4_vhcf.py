"""v10.0.5.4 VHCF first-passage wrapper.

This point release retains the audited v10.0.5.3 FEM/CZM/MPZ implementation
and adds fail-closed verification that adaptive cycle blocks are predicted by
the authoritative tensor-resolved kinetic engine. It also records the actual
physical termination condition instead of treating exhaustion of the numerical
outer-block budget as a completed fatigue result.
"""
from __future__ import annotations

import csv
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import sys
from typing import Any, Callable

from . import mode_i_first_passage_v10_0_5_3_fatigue as _v10053_original
from . import mode_i_first_passage_v10_0_5_3_fatigue_audited as _v10053_audited

POINT_RELEASE = "10.0.5.4"
MODEL_ID = "FEM_CZM_Mode_I_VHCF_first_passage_v10_0_5_4"
COMPLETION_MANIFEST = "run_completion_v10_0_5_4_vhcf.json"
VHCF_AUDIT = "progressive_fatigue_v10_0_5_4_vhcf.json"


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _option_value(args: list[str], option: str, default: str | None = None) -> str | None:
    try:
        index = args.index(option)
    except ValueError:
        return default
    return args[index + 1] if index + 1 < len(args) else default


def _float_option(
    args: list[str], option: str, default: float = math.inf
) -> float:
    value = _option_value(args, option)
    return default if value is None else float(value)


def _read_step_rows(out: Path) -> list[dict[str, str]]:
    matches = sorted(out.glob("steps_*K.csv"))
    if len(matches) != 1:
        raise RuntimeError(
            f"v10.0.5.4 expected exactly one steps CSV in {out}; found {matches}"
        )
    with matches[0].open(newline="") as handle:
        return list(csv.DictReader(handle))


def _as_float(row: dict[str, str], key: str, default: float = 0.0) -> float:
    try:
        value = float(row.get(key, default))
    except (TypeError, ValueError):
        return default
    return value if math.isfinite(value) else default


def classify_termination_v10054(
    rows: list[dict[str, str]],
    *,
    cycles_max: float,
    max_blocks: int,
    target_extension_um: float,
) -> dict[str, Any]:
    """Classify physical completion separately from numerical censoring."""

    cycles_total = sum(max(_as_float(row, "fatigue_cycles"), 0.0) for row in rows)
    last = rows[-1] if rows else {}
    extension_m = max(_as_float(last, "crack_extension_m"), 0.0)
    extension_target_m = max(float(target_extension_um), 0.0) * 1.0e-6
    n_fire_total = sum(max(int(round(_as_float(row, "n_fire"))), 0) for row in rows)

    cycle_tol = max(1.0e-6, 1.0e-10 * max(abs(cycles_max), 1.0))
    extension_tol = max(1.0e-15, 1.0e-10 * max(extension_target_m, 1.0e-12))
    reached_cycle_horizon = (
        math.isfinite(cycles_max) and cycles_total >= cycles_max - cycle_tol
    )
    reached_target_extension = (
        extension_target_m > 0.0
        and extension_m >= extension_target_m - extension_tol
    )
    exhausted_outer_blocks = len(rows) >= max(int(max_blocks), 1)

    if reached_target_extension:
        termination = "target_extension"
        status = "complete"
    elif reached_cycle_horizon:
        termination = "cycle_horizon"
        status = "complete"
    elif exhausted_outer_blocks:
        termination = "right_censored_max_blocks"
        status = "right_censored"
    else:
        termination = "terminated_before_physical_horizon"
        status = "incomplete"

    return {
        "status": status,
        "termination": termination,
        "right_censored": status == "right_censored",
        "reached_cycle_horizon": reached_cycle_horizon,
        "reached_target_extension": reached_target_extension,
        "exhausted_outer_blocks": exhausted_outer_blocks,
        "cycles_total": cycles_total,
        "cycle_horizon": cycles_max,
        "cycle_horizon_fraction": (
            cycles_total / cycles_max
            if math.isfinite(cycles_max) and cycles_max > 0.0
            else None
        ),
        "n_outer_blocks": len(rows),
        "max_outer_blocks": int(max_blocks),
        "crack_extension_m": extension_m,
        "target_extension_m": extension_target_m,
        "first_passage_observed": n_fire_total > 0,
        "n_fire_total": n_fire_total,
        "final_cleavage_clock_B": _as_float(last, "B"),
        "final_mobile_count": _as_float(last, "mpz_mobile_count"),
        "final_retained_count": _as_float(last, "mpz_retained_count"),
        "final_available_site_fraction": _as_float(
            last, "mpz_available_site_fraction", 1.0
        ),
    }


def _monitored_predictor_factory(
    original_factory: Callable[..., Callable[..., Any]],
    counters: dict[str, int],
) -> Callable[..., Callable[..., Any]]:
    """Count authoritative versus legacy predictor dispatches at runtime."""

    def factory(original: Callable[..., Any]) -> Callable[..., Any]:
        wrapped = original_factory(original)

        def integrate_one_cycle(controller, front, waveform, T_K):
            if hasattr(front, "predict_fatigue_cycle"):
                counters["authoritative_engine_calls"] += 1
            else:
                counters["legacy_fallback_calls"] += 1
            return wrapped(controller, front, waveform, T_K)

        integrate_one_cycle._v10054_predictor_monitor = True
        return integrate_one_cycle

    return factory


def main(argv: list[str] | None = None):
    args = list(sys.argv[1:] if argv is None else argv)
    out_value = _option_value(args, "--out")
    if out_value is None:
        raise SystemExit("v10.0.5.4 VHCF requires --out")
    out = Path(out_value)
    out.mkdir(parents=True, exist_ok=True)

    if _option_value(args, "--cycle-block-mode", "") != "hazard_limited":
        raise SystemExit("v10.0.5.4 VHCF requires --cycle-block-mode hazard_limited")

    cycles_max = _float_option(args, "--cycles-max")
    max_block_cycles = _float_option(args, "--max-block-cycles")
    max_blocks = int(_float_option(args, "--steps", 0.0))
    target_extension_um = _float_option(
        args, "--target-crack-extension-um", 0.0
    )
    if (
        math.isfinite(cycles_max)
        and math.isfinite(max_block_cycles)
        and max_block_cycles < cycles_max
    ):
        raise SystemExit(
            "v10.0.5.4 VHCF refuses an artificial cycle-jump ceiling below the "
            "physical cycle horizon; use --max-block-cycles inf or a value >= "
            "--cycles-max"
        )

    status_path = out / COMPLETION_MANIFEST
    status: dict[str, Any] = {
        "schema": "authoritative_run_completion_v10_0_5_4_vhcf",
        "point_release": POINT_RELEASE,
        "model": MODEL_ID,
        "started_utc": _utc_now(),
        "completed_utc": None,
        "status": "running",
        "run_completed_without_exception": False,
        "constitutive_physics_changed_in_v10054": False,
    }
    status_path.write_text(json.dumps(status, indent=2, default=str))

    counters = {
        "authoritative_engine_calls": 0,
        "legacy_fallback_calls": 0,
    }
    original_dispatch_factory = _v10053_original._fatigue_predictor_dispatch
    _v10053_original._fatigue_predictor_dispatch = _monitored_predictor_factory(
        original_dispatch_factory, counters
    )

    try:
        results = _v10053_audited.main(args)
        rows = _read_step_rows(out)
        termination = classify_termination_v10054(
            rows,
            cycles_max=cycles_max,
            max_blocks=max_blocks,
            target_extension_um=target_extension_um,
        )
        if counters["authoritative_engine_calls"] < 1:
            raise RuntimeError(
                "v10.0.5.4 completed without an authoritative engine cycle prediction"
            )
        if counters["legacy_fallback_calls"] != 0:
            raise RuntimeError(
                "v10.0.5.4 invoked the legacy scalar cycle predictor "
                f"{counters['legacy_fallback_calls']} time(s)"
            )

        source_audit = out / _v10053_original.FATIGUE_AUDIT
        audit: dict[str, Any] = {
            "schema": "progressive_fatigue_v10_0_5_4_vhcf",
            "point_release": POINT_RELEASE,
            "model": MODEL_ID,
            "source_v10053_fatigue_audit": source_audit.name,
            "source_v10053_audit_exists": source_audit.exists(),
            "authoritative_tensor_mpz_predictor_calls": counters[
                "authoritative_engine_calls"
            ],
            "legacy_scalar_predictor_calls": counters["legacy_fallback_calls"],
            "authoritative_predictor_verified": True,
            "hazard_limited_cycle_blocks": True,
            "max_block_cycles": max_block_cycles,
            "max_block_cycles_is_unbounded": not math.isfinite(max_block_cycles),
            "physical_cycle_horizon": cycles_max,
            "first_passage_clock": "cumulative cleavage action B",
            "first_passage_threshold": 1.0,
            "cycle_jump_changes_constitutive_physics": False,
            "constitutive_physics_changed": False,
            **termination,
        }
        (out / VHCF_AUDIT).write_text(
            json.dumps(audit, indent=2, default=str)
        )

        status.update(
            {
                "completed_utc": _utc_now(),
                "status": termination["status"],
                "termination": termination["termination"],
                "run_completed_without_exception": True,
                "right_censored": termination["right_censored"],
                "vhcf_audit": VHCF_AUDIT,
                "authoritative_predictor_verified": True,
            }
        )
        status_path.write_text(json.dumps(status, indent=2, default=str))
        return results
    except BaseException as exc:
        status.update(
            {
                "completed_utc": _utc_now(),
                "status": "failed",
                "run_completed_without_exception": False,
                "runtime_error_type": type(exc).__name__,
                "runtime_error": str(exc),
                "authoritative_predictor_calls_before_failure": counters[
                    "authoritative_engine_calls"
                ],
                "legacy_predictor_calls_before_failure": counters[
                    "legacy_fallback_calls"
                ],
            }
        )
        status_path.write_text(json.dumps(status, indent=2, default=str))
        raise
    finally:
        _v10053_original._fatigue_predictor_dispatch = original_dispatch_factory


if __name__ == "__main__":
    main()


__all__ = [
    "POINT_RELEASE",
    "MODEL_ID",
    "COMPLETION_MANIFEST",
    "VHCF_AUDIT",
    "classify_termination_v10054",
    "main",
]
