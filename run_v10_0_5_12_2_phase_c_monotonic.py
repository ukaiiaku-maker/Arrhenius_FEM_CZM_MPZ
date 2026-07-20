#!/usr/bin/env python3
"""Phase-C v10.0.5.12.2 runner with live, prefixed case reporting.

The underlying constitutive/mechanics entry remains v10.0.5.12.1.  This point
release changes campaign observability only: solver progress is still written to
each case ``run.log`` and selected records are mirrored to the campaign console.
A per-case ``phase_c_live_status.json`` is refreshed during execution.
"""
from __future__ import annotations

from datetime import datetime, timezone
import json
import math
from pathlib import Path
import re
import shlex
import subprocess
import threading
import time

import run_v10_0_5_12_phase_c_monotonic as _base
import run_v10_0_5_12_1_phase_c_monotonic as _radius

POINT_RELEASE = "10.0.5.12.2"
LIVE_STATUS = "phase_c_live_status.json"
_PRINT_LOCK = threading.Lock()

_PROGRESS_RE = re.compile(
    r"\[T=(?P<T>[-+0-9.eE]+)K\]\s+step\s+(?P<step>\d+)"
    r".*?KJ=\s*(?P<KJ>[-+0-9.eE]+)"
    r".*?B=\s*(?P<B>[-+0-9.eE]+)"
    r".*?N_em=\s*(?P<N>[-+0-9.eE]+)"
    r".*?a=(?P<a_mm>[-+0-9.eE]+)mm"
)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def emit(message: str, *, error: bool = False) -> None:
    """Print one complete, immediately flushed campaign-console record."""
    import sys

    with _PRINT_LOCK:
        print(message, file=sys.stderr if error else sys.stdout, flush=True)


def parse_progress_line(line: str) -> dict[str, float | int] | None:
    """Parse one standard sharp-front progress line."""
    match = _PROGRESS_RE.search(str(line))
    if match is None:
        return None
    try:
        return {
            "T_K": float(match.group("T")),
            "step": int(match.group("step")),
            "KJ_MPa_sqrt_m": float(match.group("KJ")),
            "B": float(match.group("B")),
            "N_em": float(match.group("N")),
            "a_mm": float(match.group("a_mm")),
        }
    except (TypeError, ValueError):
        return None


def selected_diagnostic_line(line: str) -> bool:
    text = str(line)
    needles = (
        "mesh GRADED",
        "da_phys=",
        "cluster J outer radius=",
        "process zone RESOLVED",
        "sources-only plasticity:",
        "crack backend:",
        "TARGET",
        "target crack extension",
    )
    return any(needle in text for needle in needles)


def write_live_status(case_dir: Path, payload: dict) -> None:
    path = case_dir / LIVE_STATUS
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, default=str))
    tmp.replace(path)


def _case_prefix(option_key: str, T_K: int) -> str:
    return f"[{option_key} T={int(T_K):04d}K]"


def run_case(py, args, root, option_key, T_K, target_um):
    """Run one case while mirroring useful progress records to stdout."""
    case_dir = root / option_key / f"T{int(T_K):04d}"
    case_dir.mkdir(parents=True, exist_ok=True)
    prefix = _case_prefix(option_key, T_K)

    if args.skip_existing and (case_dir / ".phase_c_complete").is_file():
        emit(f"{prefix} REUSE complete case: {case_dir}")
        row = _base.summarize(case_dir, option_key, T_K, target_um, 0, True, args.registry)
        emit(
            f"{prefix} COMPLETE reused ext={row.get('final_extension_um')} um "
            f"target={target_um:g} um"
        )
        return row

    option = _base.load_option(option_key, args.registry)
    cmd = _radius.build_command(py, args, option_key, T_K, target_um, case_dir)
    (case_dir / "command.txt").write_text(shlex.join(cmd) + "\n")
    (case_dir / "phase_c_case_input.json").write_text(
        json.dumps(
            {
                "schema": "phase_c_case_input_v10_0_5_12_2",
                "point_release": POINT_RELEASE,
                "option": option.audit_payload(),
                "T_K": int(T_K),
                "target_extension_um": target_um,
                "command": cmd,
            },
            indent=2,
            default=str,
        )
    )

    emit(
        f"{prefix} START candidate={option.candidate_id} target={target_um:g} um "
        f"MPZ={option.mpz_length_um:g} um/{option.mpz_n_bins} bins"
    )
    emit(f"{prefix} LOG {case_dir / 'run.log'}")
    if args.dry_run:
        emit(f"{prefix} DRY_RUN command prepared")
        return {
            "option_key": option_key,
            "candidate_id": option.candidate_id,
            "T_K": int(T_K),
            "status": "dry_run",
            "target_extension_um": target_um,
            "case_dir": str(case_dir),
        }

    started_wall = utc_now()
    started_clock = time.monotonic()
    process = None
    initial_a_mm = None
    last_progress = None
    runtime = {
        "schema": "phase_c_live_status_v10_0_5_12_2",
        "point_release": POINT_RELEASE,
        "status": "starting",
        "option_key": option.option_key,
        "candidate_id": option.candidate_id,
        "T_K": int(T_K),
        "target_extension_um": float(target_um),
        "case_dir": str(case_dir),
        "started_utc": started_wall,
        "updated_utc": started_wall,
        "completed_utc": None,
        "pid": None,
        "returncode": None,
        "elapsed_s": 0.0,
        "last_progress": None,
        "command": cmd,
    }
    write_live_status(case_dir, runtime)

    try:
        with (case_dir / "run.log").open("w", buffering=1) as log:
            process = subprocess.Popen(
                cmd,
                env=_base.case_environment(args, target_um),
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
            )
            runtime.update(status="running", pid=process.pid, updated_utc=utc_now())
            write_live_status(case_dir, runtime)
            emit(f"{prefix} RUNNING pid={process.pid}")

            if process.stdout is None:
                raise RuntimeError("Phase-C child process has no stdout pipe")

            for raw in process.stdout:
                log.write(raw)
                line = raw.rstrip("\n")
                progress = parse_progress_line(line)
                if progress is not None:
                    a_mm = float(progress["a_mm"])
                    if initial_a_mm is None:
                        initial_a_mm = a_mm
                    extension_um = max((a_mm - initial_a_mm) * 1000.0, 0.0)
                    fraction = extension_um / float(target_um) if target_um > 0 else math.nan
                    progress.update(
                        extension_um=extension_um,
                        target_extension_um=float(target_um),
                        fraction_complete=fraction,
                        advance_event="<< ADVANCE" in line,
                        raw=line,
                    )
                    last_progress = progress
                    percent = 100.0 * min(max(fraction, 0.0), 1.0)
                    marker = " ADVANCE" if progress["advance_event"] else ""
                    emit(
                        f"{prefix} PROGRESS step={progress['step']} "
                        f"KJ={progress['KJ_MPa_sqrt_m']:.3f} MPa√m "
                        f"ext={extension_um:.1f}/{target_um:g} um ({percent:.1f}%) "
                        f"B={progress['B']:.3f} N_em={progress['N_em']:.2f}{marker}"
                    )
                    runtime.update(
                        status="running",
                        updated_utc=utc_now(),
                        elapsed_s=time.monotonic() - started_clock,
                        last_progress=last_progress,
                    )
                    write_live_status(case_dir, runtime)
                elif selected_diagnostic_line(line):
                    emit(f"{prefix} CONFIG {line.strip()}")

            returncode = process.wait()

        row = _base.summarize(
            case_dir, option_key, T_K, target_um, returncode, False, args.registry
        )
        elapsed = time.monotonic() - started_clock
        final_status = str(row.get("status", "unknown"))
        runtime.update(
            status=final_status,
            updated_utc=utc_now(),
            completed_utc=utc_now(),
            returncode=int(returncode),
            elapsed_s=elapsed,
            last_progress=last_progress,
            final_extension_um=row.get("final_extension_um"),
            completion_manifest_passed=row.get("completion_manifest_passed"),
            production_manifest_passed=row.get("production_manifest_passed"),
        )
        write_live_status(case_dir, runtime)

        label = final_status.upper()
        message = (
            f"{prefix} {label} returncode={returncode} "
            f"ext={row.get('final_extension_um')}/{target_um:g} um "
            f"elapsed={elapsed / 60.0:.1f} min "
            f"completion_manifest={row.get('completion_manifest_passed')} "
            f"production_manifest={row.get('production_manifest_passed')}"
        )
        emit(message, error=final_status != "complete")
        return row

    except BaseException as exc:
        if process is not None and process.poll() is None:
            process.terminate()
            try:
                process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait()
        elapsed = time.monotonic() - started_clock
        runtime.update(
            status="runner_exception",
            updated_utc=utc_now(),
            completed_utc=utc_now(),
            returncode=None if process is None else process.returncode,
            elapsed_s=elapsed,
            last_progress=last_progress,
            runtime_error_type=type(exc).__name__,
            runtime_error=str(exc),
        )
        write_live_status(case_dir, runtime)
        emit(
            f"{prefix} FAILED {type(exc).__name__}: {exc} "
            f"elapsed={elapsed / 60.0:.1f} min",
            error=True,
        )
        raise


def main() -> None:
    saved_build = _base.build_command
    saved_run_case = _base.run_case
    saved_point_release = _base.POINT_RELEASE
    _base.build_command = _radius.build_command
    _base.run_case = run_case
    _base.POINT_RELEASE = POINT_RELEASE
    try:
        _base.main()
    finally:
        _base.build_command = saved_build
        _base.run_case = saved_run_case
        _base.POINT_RELEASE = saved_point_release


if __name__ == "__main__":
    main()


__all__ = [
    "POINT_RELEASE",
    "LIVE_STATUS",
    "parse_progress_line",
    "selected_diagnostic_line",
    "run_case",
    "main",
]
