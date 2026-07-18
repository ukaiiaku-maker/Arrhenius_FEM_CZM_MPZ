"""v10.0.5.3 progressive fatigue entry point.

The point release reuses the certified v10.0.5.2 mechanics and constitutive
implementation. Runtime source adapters route the existing fatigue load/block
selection into the progressive trial-CZM lifecycle rather than the retired
abrupt fatigue renewal path.
"""
from __future__ import annotations

import inspect
import json
from pathlib import Path
import sys
import textwrap
from typing import Any, Callable

import numpy as np

from . import fatigue_v1 as _fatigue_v1
from . import kinetic_progressive_2d_v1002 as _v1002
from . import mode_i_first_passage_v10_0_5_2_parallel as _v10052_entry
from .kinetic_campaign_czm_v10053 import engine_factory_v10053
from .tensor_resolved_coupling_v1005 import latest_tensor_drive

POINT_RELEASE = "10.0.5.3"
MODEL_ID = "FEM_CZM_Mode_I_progressive_cycle_block_fatigue_v10_0_5_3"
COMPLETION_MANIFEST = "run_completion_v10_0_5_3_fatigue.json"
FATIGUE_AUDIT = "progressive_fatigue_v10_0_5_3.json"


def _utc_now() -> str:
    from datetime import datetime, timezone

    return datetime.now(timezone.utc).isoformat()


def _replace_unique(source: str, old: str, new: str, name: str) -> str:
    count = source.count(old)
    if count != 1:
        raise RuntimeError(
            f"v10.0.5.3 expected exactly one {name} anchor; found {count}"
        )
    return source.replace(old, new)


def _patch_run_2d_source(source: str) -> str:
    """Route only the v10 progressive engine away from legacy fatigue commit."""
    old = (
        "            else:\n"
        "                if fatigue_mode:\n"
        "                    wave_acc = FatigueWaveform(\n"
    )
    new = (
        "            else:\n"
        "                if fatigue_mode and not kinetic_progressive:\n"
        "                    wave_acc = FatigueWaveform(\n"
    )
    return _replace_unique(source, old, new, "single-front fatigue dispatch")


def _build_progressive_run_2d_v10053(original_run_2d):
    """Compile the v10.0.2 lifecycle with guarded cyclic-loading support."""
    original_getsource = _v1002.inspect.getsource

    def patched_getsource(obj):
        source = original_getsource(obj)
        if obj is original_run_2d:
            source = _patch_run_2d_source(source)
        return source

    builder_source = textwrap.dedent(inspect.getsource(_ORIGINAL_V1002_BUILDER))
    builder_source = _replace_unique(
        builder_source,
        "if deflect or fatigue_mode or max_fronts != 1:",
        "if deflect or max_fronts != 1:",
        "progressive fatigue guard",
    )
    builder_source = _replace_unique(
        builder_source,
        "single-front monotonic Mode I with branching disabled",
        "single-front Mode I (monotonic or v10.0.5.3 fatigue) with branching disabled",
        "progressive mode diagnostic",
    )
    builder_source = _replace_unique(
        builder_source,
        "total_dt_s=dt_cur,",
        "total_dt_s=(fatigue_cycles_accepted / max(float(getattr(args, 'frequency_Hz', 1.0e3) or 1.0e3), 1.0e-300) if fatigue_mode else dt_cur),",
        "cycle-block physical time",
    )
    target_block = """                        def target_reached_v1002():
                            return (
                                np.isfinite(target_m_v1002)
                                and eng.checkpoint_advance_total_m
                                >= target_m_v1002 - max(1.0e-12 * da_phys, 1.0e-15)
                            )
"""
    target_replacement = """                        fatigue_commits_v10053 = 0

                        def target_reached_v1002():
                            return (
                                (fatigue_mode and fatigue_commits_v10053 >= 1)
                                or (
                                    np.isfinite(target_m_v1002)
                                    and eng.checkpoint_advance_total_m
                                    >= target_m_v1002 - max(1.0e-12 * da_phys, 1.0e-15)
                                )
                            )
"""
    builder_source = _replace_unique(
        builder_source,
        target_block,
        target_replacement,
        "one-fatigue-commit-per-outer-block target",
    )
    commit_nonlocal = """                        def on_commit_v1002(context_v1002, result_v1002):
                            nonlocal a_tip, Kc_first, Kc_first_step
                            nonlocal x, y, cx_e, cy_e, adj
"""
    commit_nonlocal_replacement = """                        def on_commit_v1002(context_v1002, result_v1002):
                            nonlocal a_tip, Kc_first, Kc_first_step
                            nonlocal x, y, cx_e, cy_e, adj
                            nonlocal fatigue_commits_v10053
                            fatigue_commits_v10053 += 1
"""
    builder_source = _replace_unique(
        builder_source,
        commit_nonlocal,
        commit_nonlocal_replacement,
        "fatigue commit counter",
    )
    accepted_anchor = (
        "                        _V1002_RUNTIME['damage_rejections'] += int(\n"
        "                            lifecycle_v1002.rejected_attempts)\n"
    )
    accepted_replacement = (
        "                        if fatigue_mode:\n"
        "                            fatigue_cycles_accepted = (\n"
        "                                float(lifecycle_v1002.consumed_dt_s)\n"
        "                                * max(float(getattr(args, 'frequency_Hz', 1.0e3) or 1.0e3), 0.0)\n"
        "                            )\n"
        + accepted_anchor
    )
    builder_source = _replace_unique(
        builder_source,
        accepted_anchor,
        accepted_replacement,
        "consumed-cycle accounting",
    )

    namespace = dict(_ORIGINAL_V1002_BUILDER.__globals__)
    exec(compile(builder_source, "<v10_0_5_3_progressive_builder>", "exec"), namespace)
    builder = namespace[_ORIGINAL_V1002_BUILDER.__name__]

    _v1002.inspect.getsource = patched_getsource
    try:
        transformed = builder(original_run_2d)
    finally:
        _v1002.inspect.getsource = original_getsource
    transformed._v10053_progressive_fatigue = True
    transformed._v10053_legacy_fatigue_commit_bypassed = True
    transformed._v10053_constitutive_physics_changed = False
    return transformed


_ORIGINAL_V1002_BUILDER = _v1002.build_progressive_run_2d_v1002


def _fatigue_predictor_dispatch(original: Callable[..., Any]) -> Callable[..., Any]:
    def integrate_one_cycle(controller, front, waveform, T_K):
        if hasattr(front, "predict_fatigue_cycle"):
            drive = latest_tensor_drive()
            system_weights = np.asarray(
                drive["slip_system_drive_factors"], dtype=float
            )
            pz = front.predict_fatigue_cycle(
                waveform,
                T_K,
                controller.cfg.n_phase,
                system_weights=system_weights,
            )
            return _fatigue_v1.CycleHazardResult(
                mu_emit=float(pz.get("dN_emit_per_cycle", 0.0)),
                mu_peierls=float(pz.get("dN_peierls_per_cycle", 0.0)),
                mu_taylor=float(pz.get("dN_taylor_per_cycle", 0.0)),
                mu_escape=float(pz.get("dN_escape_per_cycle", 0.0)),
                mu_cleave=float(pz.get("mu_cleave_per_cycle", 0.0)),
                store_per_cycle=float(
                    pz.get("dN_store_per_cycle", pz.get("dN_emit_per_cycle", 0.0))
                ),
                mobile_per_cycle=float(
                    pz.get("dN_mobile_per_cycle", pz.get("dN_emit_per_cycle", 0.0))
                ),
                escape_per_cycle=float(pz.get("dN_escape_per_cycle", 0.0)),
                peierls_per_cycle=float(pz.get("dN_peierls_per_cycle", 0.0)),
                taylor_per_cycle=float(pz.get("dN_taylor_per_cycle", 0.0)),
                avg_sigma_tip=float(pz.get("avg_sigma_tip", 0.0)),
                max_sigma_tip=float(pz.get("max_sigma_tip", 0.0)),
                avg_sigma_emit_eff=float(pz.get("avg_sigma_emit_eff", 0.0)),
                storage_fraction=0.0,
            )
        return original(controller, front, waveform, T_K)

    integrate_one_cycle._v10053_fatigue_predictor_dispatch = True
    return integrate_one_cycle


def _option_value(args: list[str], option: str) -> str | None:
    try:
        index = args.index(option)
    except ValueError:
        return None
    return args[index + 1] if index + 1 < len(args) else None


def _write_status(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, default=str))


def main(argv: list[str] | None = None):
    args = list(sys.argv[1:] if argv is None else argv)
    if "--fatigue-cycles" not in args:
        raise SystemExit("v10.0.5.3 fatigue requires --fatigue-cycles")
    if "--fatigue-hold-load" not in args:
        raise SystemExit(
            "v10.0.5.3 fatigue requires --fatigue-hold-load after the amplitude ramp"
        )
    out_value = _option_value(args, "--out")
    if out_value is None:
        raise SystemExit("v10.0.5.3 fatigue requires --out")
    out = Path(out_value)
    status_path = out / COMPLETION_MANIFEST
    status = {
        "schema": "authoritative_run_completion_v10_0_5_3_fatigue",
        "point_release": POINT_RELEASE,
        "model": MODEL_ID,
        "started_utc": _utc_now(),
        "completed_utc": None,
        "status": "running",
        "run_completed_without_exception": False,
        "constitutive_physics_changed_in_v10053": False,
    }
    _write_status(status_path, status)

    original_factory = _v10052_entry.engine_factory_v10052
    original_builder = _v1002.build_progressive_run_2d_v1002
    original_predictor = _fatigue_v1.FatigueCycleHazardController.integrate_one_cycle
    _v10052_entry.engine_factory_v10052 = engine_factory_v10053
    _v1002.build_progressive_run_2d_v1002 = _build_progressive_run_2d_v10053
    _fatigue_v1.FatigueCycleHazardController.integrate_one_cycle = (
        _fatigue_predictor_dispatch(original_predictor)
    )

    try:
        results = _v10052_entry.main(args)
    except BaseException as exc:
        status.update(
            {
                "completed_utc": _utc_now(),
                "status": "failed",
                "runtime_error_type": type(exc).__name__,
                "runtime_error": str(exc),
            }
        )
        _write_status(status_path, status)
        raise
    finally:
        _v10052_entry.engine_factory_v10052 = original_factory
        _v1002.build_progressive_run_2d_v1002 = original_builder
        _fatigue_v1.FatigueCycleHazardController.integrate_one_cycle = original_predictor

    progressive_path = out / "kinetic_campaign_czm_progressive_2d_v10_0_3.json"
    progressive = json.loads(progressive_path.read_text())
    records = progressive.get("records")
    if records is None and isinstance(progressive.get("v1002_runtime_after_run"), dict):
        records = progressive["v1002_runtime_after_run"].get("records")
    records = records if isinstance(records, list) else []
    accepted_time = sum(float(row.get("dt_consumed_s", 0.0)) for row in records)
    frequency = float(_option_value(args, "--frequency-Hz") or 1.0e3)
    audit = {
        "schema": "progressive_fatigue_v10_0_5_3",
        "point_release": POINT_RELEASE,
        "model": MODEL_ID,
        "source_progressive_runtime": progressive_path.name,
        "accepted_substep_record_count": len(records),
        "accepted_physical_time_s_from_records": accepted_time,
        "accepted_cycles_from_records": accepted_time * frequency,
        "frequency_Hz": frequency,
        "legacy_fatigue_state_commit_bypassed": True,
        "progressive_trial_czm_lifecycle_active": bool(
            progressive.get("full_progressive_trial_loop_active", False)
        ),
        "v10_0_5_2_constitutive_physics_reused": True,
        "fatigue_loading_changes_constitutive_physics": False,
        "cycle_jump_quadrature": "midpoint phase exposure",
        "maximum_commits_per_accepted_outer_block": 1,
        "unused_cycle_time_after_commit_is_not_accepted": True,
    }
    if not audit["progressive_trial_czm_lifecycle_active"]:
        raise RuntimeError("v10.0.5.3 fatigue completed without progressive lifecycle")
    (out / FATIGUE_AUDIT).write_text(json.dumps(audit, indent=2))

    status.update(
        {
            "completed_utc": _utc_now(),
            "status": "complete",
            "run_completed_without_exception": True,
            "fatigue_audit": FATIGUE_AUDIT,
        }
    )
    _write_status(status_path, status)
    return results


if __name__ == "__main__":
    main()
