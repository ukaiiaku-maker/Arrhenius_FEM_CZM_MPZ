"""v10.0.5.2 Mode-I entry point with complete per-channel diagnostics.

The certified v10.0.5 mechanics, tensor projection, barriers, exact depletion,
transport, cohesive lifecycle, and geometry gates are unchanged.  This wrapper
adds lossless channel aggregation and an authoritative invocation completion
manifest for long-growth gates.
"""
from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import sys
from typing import Any, Callable

import numpy as np

from . import j_integral as _j_integral
from . import kinetic_progressive_2d_v1002 as _v1002
from . import mode_i_first_passage_v10_0_3_progressive as _v1003
from .kinetic_campaign_czm_v10052 import engine_factory_v10052
from .mode_i_first_passage_v10_0_5_parallel import (
    _option_value,
    _write_v1005_results,
    make_progressive_formatter_v1005,
)
from .tensor_resolved_coupling_v1005 import (
    TensorResolvedDriveConfig,
    TensorResolvedKineticCohesiveStepper,
    make_tensor_resolved_J_wrapper,
    reset_tensor_drive_runtime,
    tensor_drive_runtime_payload,
)

MODEL_ID = "FEM_CZM_Mode_I_parallel_opening_tensor_emission_v10_0_5_2"
POINT_RELEASE = "10.0.5.2"
COMPLETION_MANIFEST = "run_completion_v10_0_5_2.json"
CHANNEL_AUDIT = "parallel_channel_diagnostics_v10_0_5_2.json"
RESULTS_NAME = "mode_i_v10_0_5_2_results.json"
PROGRESSIVE_NAME = "kinetic_campaign_czm_progressive_2d_v10_0_3.json"
QUALITY_AUDIT_NAME = "explicit_quality_wrapper_chain_v91856.json"


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _records(payload: dict[str, Any]) -> list[dict[str, Any]]:
    rows = payload.get("records")
    if rows is None and isinstance(payload.get("v1002_runtime_after_run"), dict):
        rows = payload["v1002_runtime_after_run"].get("records")
    if not isinstance(rows, list) or any(not isinstance(row, dict) for row in rows):
        raise RuntimeError("progressive runtime lacks accepted-substep records")
    return [dict(row) for row in rows]


def _channel_indices(row: dict[str, Any]) -> list[int]:
    indices = []
    prefix = "slip_drive_factor_"
    for key in row:
        if str(key).startswith(prefix):
            suffix = str(key)[len(prefix):]
            if suffix.isdigit():
                indices.append(int(suffix))
    indices = sorted(set(indices))
    if indices != list(range(len(indices))):
        raise RuntimeError(f"non-contiguous slip-trace channel indices: {indices}")
    return indices


def make_progressive_formatter_v10052(
    original_formatter: Callable[..., dict[str, Any]],
) -> Callable[..., dict[str, Any]]:
    """Add fail-closed channel partition checks to the v10.0.5 formatter."""

    base_formatter = make_progressive_formatter_v1005(original_formatter)

    def formatter(engine, step_result, KJ, N_em_pre):
        out = dict(base_formatter(engine, step_result, KJ, N_em_pre))
        indices = _channel_indices(out)
        if not indices:
            raise RuntimeError("v10.0.5.2 formatter found no slip-trace channels")
        rates = []
        emitted = []
        for index in indices:
            rate_key = f"lambda_emit_s-1_{index}"
            emit_key = f"dN_emit_{index}"
            if rate_key not in out or emit_key not in out:
                raise RuntimeError(
                    "v10.0.5.2 incomplete per-channel diagnostics: "
                    f"missing {rate_key if rate_key not in out else emit_key}"
                )
            rate = float(out[rate_key])
            count = float(out[emit_key])
            if not np.isfinite(rate) or rate < 0.0:
                raise RuntimeError(f"invalid per-channel emission rate {rate_key}={rate}")
            if not np.isfinite(count) or count < -1.0e-14:
                raise RuntimeError(f"invalid per-channel emitted increment {emit_key}={count}")
            rates.append(rate)
            emitted.append(max(count, 0.0))

        scalar_total = float(out.get("dN_emit_block", sum(emitted)))
        vector_total = float(sum(emitted))
        tolerance = max(1.0e-12, 1.0e-10 * max(abs(scalar_total), 1.0))
        residual = vector_total - scalar_total
        if abs(residual) > tolerance:
            raise RuntimeError(
                "v10.0.5.2 channel emission partition mismatch: "
                f"sum={vector_total:.16g}, scalar={scalar_total:.16g}, "
                f"residual={residual:.16g}"
            )
        out.update({
            "per_channel_strang_diagnostics_complete": True,
            "per_channel_emitted_increment_semantics": (
                "sum_over_all_strang_half_steps"
            ),
            "per_channel_hazard_semantics": (
                "last_rate_after_final_strang_half_step"
            ),
            "per_channel_emission_partition_residual": residual,
            "per_channel_emission_partition_verified": True,
            "per_channel_hazard_count": len(rates),
            "per_channel_emitted_increment_count": len(emitted),
        })
        return out

    formatter.__name__ = "_v10_format_progressive_info_v10052"
    formatter._v10052_complete_channel_reporting = True
    formatter._v10052_original = original_formatter
    return formatter


def _write_status(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, default=str))


def _write_results(out: Path, results: Any, channel_audit: str) -> Path:
    rows = []
    for source in results or []:
        row = dict(source)
        row.update({
            "model": MODEL_ID,
            "point_release": POINT_RELEASE,
            "front_state_model": "kinetic_campaign_czm",
            "front_state_model_detail": (
                "pf_v10_1_7_1_parallel_opening_tensor_emission_reset_safe_v10052"
            ),
            "tensor_resolved_parallel_coupling": True,
            "per_channel_strang_diagnostics_complete": True,
            "per_channel_emitted_increment_semantics": (
                "sum_over_all_strang_half_steps"
            ),
            "per_channel_hazard_semantics": (
                "last_rate_after_final_strang_half_step"
            ),
            "constitutive_physics_changed_in_v10052": False,
            "parallel_channel_diagnostics_audit": channel_audit,
            "authoritative_completion_manifest": COMPLETION_MANIFEST,
        })
        rows.append(row)
    path = out / RESULTS_NAME
    path.write_text(json.dumps(rows, indent=2, default=str))
    return path


def main(argv: list[str] | None = None):
    args = list(sys.argv[1:] if argv is None else argv)
    out_value = _option_value(args, "--out")
    if out_value is None:
        raise RuntimeError("v10.0.5.2 requires --out")
    out = Path(out_value)
    status_path = out / COMPLETION_MANIFEST
    status = {
        "schema": "authoritative_run_completion_v10_0_5_2",
        "point_release": POINT_RELEASE,
        "started_utc": _utc_now(),
        "completed_utc": None,
        "status": "running",
        "run_completed_without_exception": False,
        "runtime_error_type": None,
        "runtime_error": None,
        "authoritative_for_this_invocation": True,
        "constitutive_physics_changed_in_v10052": False,
    }
    _write_status(status_path, status)

    theta = float(_option_value(args, "--crystal-theta-deg") or 45.0)
    reset_tensor_drive_runtime(
        TensorResolvedDriveConfig(
            crystal_theta_deg=theta,
            probe_radius_m=float(
                os.environ.get("ARRHENIUS_TENSOR_DRIVE_PROBE_RADIUS_M", "1e-5")
            ),
            sector_half_angle_deg=float(
                os.environ.get(
                    "ARRHENIUS_TENSOR_DRIVE_SECTOR_HALF_ANGLE_DEG", "25"
                )
            ),
            damage_cutoff=float(
                os.environ.get("ARRHENIUS_TENSOR_DRIVE_DAMAGE_CUTOFF", "0.85")
            ),
            min_elements=int(
                os.environ.get("ARRHENIUS_TENSOR_DRIVE_MIN_ELEMENTS", "3")
            ),
            schmid_reference=0.5,
        )
    )

    original_J = _j_integral.compute_J_integral
    original_stepper = _v1002.KineticCohesiveStepper
    original_formatter = _v1002._v10_format_progressive_info
    original_factory = _v1003.engine_factory_v1003
    _j_integral.compute_J_integral = make_tensor_resolved_J_wrapper(original_J)
    _v1002.KineticCohesiveStepper = TensorResolvedKineticCohesiveStepper
    _v1002._v10_format_progressive_info = make_progressive_formatter_v10052(
        original_formatter
    )
    _v1003.engine_factory_v1003 = engine_factory_v10052

    try:
        results = _v1003.main(args)
    except BaseException as exc:
        status.update({
            "completed_utc": _utc_now(),
            "status": "failed",
            "run_completed_without_exception": False,
            "runtime_error_type": type(exc).__name__,
            "runtime_error": str(exc),
        })
        _write_status(status_path, status)
        raise
    finally:
        _j_integral.compute_J_integral = original_J
        _v1002.KineticCohesiveStepper = original_stepper
        _v1002._v10_format_progressive_info = original_formatter
        _v1003.engine_factory_v1003 = original_factory

    try:
        tensor_audit = tensor_drive_runtime_payload()
        if not tensor_audit.get("implementation_certified", False):
            raise RuntimeError(
                "v10.0.5.2 tensor-resolved coupling audit failed: "
                + json.dumps(tensor_audit, default=str)
            )
        if int(tensor_audit.get("nonzero_emission_drive_capture_count", 0)) < 1:
            raise RuntimeError("v10.0.5.2 captured no nonzero slip-trace drive")

        tensor_path = out / "parallel_opening_emission_v10_0_5_audit.json"
        tensor_path.write_text(json.dumps(tensor_audit, indent=2, default=str))
        legacy_results_path = _write_v1005_results(out, results, tensor_path.name)

        progressive_path = out / PROGRESSIVE_NAME
        progressive = json.loads(progressive_path.read_text())
        records = _records(progressive)
        if not records:
            raise RuntimeError("v10.0.5.2 progressive runtime contains no records")

        partition_residual_max = 0.0
        channel_count = None
        for record_index, row in enumerate(records):
            indices = _channel_indices(row)
            if not indices:
                raise RuntimeError(
                    f"progressive record {record_index} has no channel drives"
                )
            if channel_count is None:
                channel_count = len(indices)
            elif len(indices) != channel_count:
                raise RuntimeError("slip-trace channel count changed during the run")
            for index in indices:
                for prefix in ("lambda_emit_s-1_", "dN_emit_"):
                    key = f"{prefix}{index}"
                    if key not in row or not np.isfinite(float(row[key])):
                        raise RuntimeError(
                            f"progressive record {record_index} lacks finite {key}"
                        )
            residual = abs(float(row.get("per_channel_emission_partition_residual", 0.0)))
            partition_residual_max = max(partition_residual_max, residual)
            if row.get("per_channel_emission_partition_verified") is not True:
                raise RuntimeError(
                    f"progressive record {record_index} lacks partition verification"
                )

        quality_path = out / QUALITY_AUDIT_NAME
        quality = json.loads(quality_path.read_text()) if quality_path.is_file() else {}
        quality_consistent = quality.get("run_completed_without_exception") is True
        if not quality_consistent:
            raise RuntimeError(
                "v10.0.5.2 completed solver return conflicts with the invocation "
                "quality audit: " + json.dumps({
                    "runtime_error_type": quality.get("runtime_error_type"),
                    "runtime_error": quality.get("runtime_error"),
                }, default=str)
            )

        channel_audit = {
            "schema": "parallel_channel_diagnostics_v10_0_5_2",
            "point_release": POINT_RELEASE,
            "source_progressive_runtime": PROGRESSIVE_NAME,
            "accepted_substep_record_count": len(records),
            "slip_trace_channel_count": int(channel_count or 0),
            "finite_hazard_rows": len(records) * int(channel_count or 0),
            "finite_emitted_increment_rows": len(records) * int(channel_count or 0),
            "per_channel_strang_diagnostics_complete": True,
            "per_channel_emitted_increment_semantics": (
                "sum_over_all_strang_half_steps"
            ),
            "per_channel_hazard_semantics": (
                "last_rate_after_final_strang_half_step"
            ),
            "maximum_emission_partition_residual": partition_residual_max,
            "quality_audit_consistent_with_completed_invocation": True,
            "constitutive_physics_changed_in_v10052": False,
            "implementation_certified": True,
        }
        channel_path = out / CHANNEL_AUDIT
        channel_path.write_text(json.dumps(channel_audit, indent=2, default=str))
        results_path = _write_results(out, results, channel_path.name)

        source_files = [
            tensor_path,
            legacy_results_path,
            progressive_path,
            quality_path,
            channel_path,
            results_path,
        ]
        hashes = {path.name: _sha256(path) for path in source_files if path.is_file()}
        status.update({
            "completed_utc": _utc_now(),
            "status": "complete",
            "run_completed_without_exception": True,
            "runtime_error_type": None,
            "runtime_error": None,
            "quality_audit_consistent_with_completed_invocation": True,
            "per_channel_strang_diagnostics_complete": True,
            "accepted_substep_record_count": len(records),
            "slip_trace_channel_count": int(channel_count or 0),
            "output_sha256": hashes,
            "parallel_channel_diagnostics_audit": CHANNEL_AUDIT,
            "results": RESULTS_NAME,
        })
        _write_status(status_path, status)
    except BaseException as exc:
        status.update({
            "completed_utc": _utc_now(),
            "status": "failed_post_run_audit",
            "run_completed_without_exception": False,
            "runtime_error_type": type(exc).__name__,
            "runtime_error": str(exc),
        })
        _write_status(status_path, status)
        raise

    print("V10.0.5.2 COMPLETE PER-CHANNEL DIAGNOSTICS GATE PASSED")
    print(json.dumps({
        "accepted_substep_records": len(records),
        "slip_trace_channel_count": int(channel_count or 0),
        "maximum_emission_partition_residual": partition_residual_max,
        "quality_audit_consistent": True,
        "constitutive_physics_changed": False,
        "completion_manifest": str(status_path),
        "channel_audit": str(channel_path),
        "results": str(results_path),
    }, indent=2))
    return results


if __name__ == "__main__":
    main()
