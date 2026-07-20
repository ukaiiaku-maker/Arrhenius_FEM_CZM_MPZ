#!/usr/bin/env python3
"""Run the v10.0.5.13 four-option barrier-only FEM/CZM matrix."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import threading

import pandas as pd

from arrhenius_fracture.barrier_only_response_registry_v100513 import (
    PRIMARY_OPTION_KEYS,
    TWO_D_STATE_POLICY,
    load_barrier_option,
)
from arrhenius_fracture.mode_i_first_passage_v10_0_5_13_barrier_only import (
    PRODUCTION_MANIFEST,
)
from run_four_class_exp_floor_czm_500um_sweep import (
    completion_status,
    extract_r_curve,
)

POINT_RELEASE = "10.0.5.13"
TEMPERATURES_FULL = tuple(range(300, 1201, 100))
STATUS_FILE = "barrier_only_case_status_v10_0_5_13.json"
CAMPAIGN_STATUS = "barrier_only_campaign_status_v10_0_5_13.json"
_PRINT_LOCK = threading.Lock()


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _values(text: str, cast=str):
    return [cast(x) for x in str(text).replace(",", " ").split() if x]


def _fstr(value: float) -> str:
    value = float(value)
    return "inf" if math.isinf(value) else f"{value:.16g}"


def _read_json(path: Path) -> dict:
    if not path.is_file():
        return {}
    try:
        data = json.loads(path.read_text())
    except Exception:
        return {}
    if isinstance(data, list):
        return data[0] if data and isinstance(data[0], dict) else {}
    return data if isinstance(data, dict) else {}


def _write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, default=str))
    tmp.replace(path)


def _resolve_python(args) -> str:
    if args.python_bin:
        path = Path(args.python_bin).expanduser().resolve()
        if not path.is_file():
            raise SystemExit(f"PYTHON_BIN does not exist: {path}")
        return str(path)
    if os.environ.get("CONDA_DEFAULT_ENV") == args.conda_env:
        return sys.executable
    cp = subprocess.run(
        ["conda", "run", "-n", args.conda_env, "python", "-c", "import sys; print(sys.executable)"],
        text=True,
        capture_output=True,
    )
    if cp.returncode != 0:
        raise SystemExit(f"could not resolve Python in {args.conda_env!r}:\n{cp.stderr}")
    lines = [line.strip() for line in cp.stdout.splitlines() if line.strip()]
    if not lines:
        raise SystemExit(f"Conda environment {args.conda_env!r} returned no Python path")
    return lines[-1]


def _selected_matrix(args):
    if args.mode == "smoke":
        defaults = (["dbtt_primary"], [700], 20.0)
    else:
        defaults = (list(PRIMARY_OPTION_KEYS), list(TEMPERATURES_FULL), 100.0)
    options = _values(args.options) if args.options else defaults[0]
    temperatures = (
        [int(round(x)) for x in _values(args.temperatures, float)]
        if args.temperatures
        else defaults[1]
    )
    target = float(args.target_extension_um if args.target_extension_um is not None else defaults[2])
    options = list(dict.fromkeys(options))
    temperatures = list(dict.fromkeys(temperatures))
    for option in options:
        load_barrier_option(option, args.registry)
    if not options or not temperatures or target <= 0.0 or any(T <= 0 for T in temperatures):
        raise SystemExit("matrix requires positive temperatures and target extension")
    return options, temperatures, target


def _build_command(py: str, args, option_key: str, T_K: int, target_um: float, case_dir: Path):
    cmd = [
        py,
        "-m",
        "arrhenius_fracture.mode_i_first_passage_v10_0_5_13_barrier_only",
        "--barrier-option",
        option_key,
        "--tip-refinement-radius-um",
        _fstr(args.tip_refinement_radius_um),
        "--selected-cluster-J-outer-um",
        _fstr(args.cluster_J_outer_um),
        "--local-J-outer-um",
        _fstr(args.local_J_outer_um),
        "--mode",
        "2d",
        "--bulk-plasticity-mode",
        "bulk_same_pt_km",
        "--temperatures",
        str(int(T_K)),
        "--steps",
        str(args.steps),
        "--nx",
        str(args.nx),
        "--ny",
        str(args.ny),
        "--tip-h-fine",
        _fstr(args.tip_h_fine),
        "--tip-ratio",
        _fstr(args.tip_ratio),
        "--dU",
        _fstr(args.dU),
        "--dt",
        _fstr(args.dt),
        "--n-stagger",
        str(args.n_stagger),
        "--print-every",
        str(args.print_every),
        "--adaptive-events",
        "--adaptive-event-target",
        _fstr(args.adaptive_event_target),
        "--adaptive-min-frac",
        "1e-8",
        "--adaptive-grow",
        "4",
        "--da-phys",
        _fstr(args.da_um * 1.0e-6),
        "--target-crack-extension-um",
        _fstr(target_um),
        "--crystal-aniso",
        "--crystal-compete",
        "--crystal-theta-deg",
        _fstr(args.theta_deg),
        "--crystal-C11",
        "523e9",
        "--crystal-C12",
        "203e9",
        "--crystal-C44",
        "160e9",
        "--cleave-gamma-aniso",
        "0.3",
        "--crystal-material",
        "w",
        "--max-fronts",
        "1",
        "--crack-backend",
        "adaptive_czm",
        "--czm-max-angle-error-deg",
        "35",
        "--j-decomposition",
        "cluster",
        "--mpz-length-um",
        _fstr(TWO_D_STATE_POLICY["mpz_length_um"]),
        "--mpz-n-bins",
        str(TWO_D_STATE_POLICY["mpz_n_bins"]),
        "--save-snapshots",
        str(args.save_snapshots),
        "--snapshot-cols",
        str(args.snapshot_cols),
        "--snapshot-by-crack-extension-um",
        _fstr(args.snapshot_interval_um),
        "--no-plots",
        "--out",
        str(case_dir),
    ]
    if args.registry:
        cmd[3:3] = ["--barrier-registry", str(args.registry.resolve())]
    return cmd


def _case_environment(args, target_um: float):
    env = os.environ.copy()
    env.update(
        {
            "ARRHENIUS_COMMITTED_TARGET_EXTENSION_UM": _fstr(target_um),
            "ARRHENIUS_PREFINED_MODE_I_CORRIDOR": "1",
            "ARRHENIUS_MIN_ACCEPTED_TRIANGLE_QUALITY": _fstr(args.min_triangle_quality),
            "ARRHENIUS_MIN_ACCEPTED_CHILD_AREA_RATIO": _fstr(args.min_child_area_ratio),
            "ARRHENIUS_MAX_TIP_H_OVER_DA": _fstr(args.max_tip_h_over_da),
            "ARRHENIUS_MAX_IDENTICAL_GEOMETRY_VETOES": str(args.max_identical_geometry_vetoes),
            "ARRHENIUS_EVENT_STATISTICS": "deterministic",
            "ARRHENIUS_STOCHASTIC_EMISSION": "0",
            "ARRHENIUS_PROPAGATION_CONTROL": "raw",
            "PYTHONUNBUFFERED": "1",
        }
    )
    return env


def _case_is_complete(case_dir: Path, option_key: str, target_um: float) -> bool:
    status = _read_json(case_dir / STATUS_FILE)
    production = _read_json(case_dir / PRODUCTION_MANIFEST)
    bulk = _read_json(case_dir / "bulk_state_v9_11_summary.json")
    complete, extension_um = completion_status(case_dir, target_um)
    return bool(
        status.get("status") == "complete"
        and status.get("option_key") == option_key
        and math.isclose(float(status.get("target_extension_um", -1.0)), target_um, abs_tol=1.0e-9)
        and complete
        and extension_um is not None
        and production.get("run_completed_without_exception") is True
        and production.get("candidate_state_fields_applied") is False
        and production.get("mesh_refinement_runtime", {}).get("actual_radius_verified") is True
        and production.get("barrier_option", {}).get("option_key") == option_key
        and bulk.get("bulk_explicit_mobile_retained_state") is True
        and int(bulk.get("bulk_state_update_calls", 0)) > 0
    )


def _archive_and_clean(case_dir: Path, root: Path, option_key: str, T_K: int) -> None:
    if not case_dir.exists() or not any(case_dir.iterdir()):
        return
    archive = root / "interrupted_case_logs"
    archive.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%dT%H%M%S")
    diagnostic = archive / f"{option_key}_T{T_K}K_{stamp}.txt"
    parts = [f"case_dir={case_dir}", f"archived_utc={_utc_now()}"]
    for name in ("command.txt", STATUS_FILE, PRODUCTION_MANIFEST):
        path = case_dir / name
        if path.is_file():
            parts.extend([f"\n===== {name} =====", path.read_text(errors="replace")])
    log = case_dir / "run.log"
    if log.is_file():
        parts.append("\n===== run.log tail =====")
        parts.extend(log.read_text(errors="replace").splitlines()[-200:])
    diagnostic.write_text("\n".join(parts) + "\n")
    shutil.rmtree(case_dir)
    with _PRINT_LOCK:
        print(f"RESTART clean: {option_key}/T{T_K}K; diagnostic={diagnostic}", flush=True)


def _first_passage(case_dir: Path):
    data = _read_json(case_dir / "anisotropic_calibrated_tip_first_passage_summary.json")
    for key in ("KJ_reference_first_MPa_sqrt_m", "Kc_first_existing_MPa_sqrt_m"):
        try:
            value = float(data[key])
        except (KeyError, TypeError, ValueError):
            continue
        if math.isfinite(value):
            return value
    return None


def _summarize(case_dir: Path, option_key: str, T_K: int, target_um: float, returncode: int, reused: bool):
    option = load_barrier_option(option_key)
    production = _read_json(case_dir / PRODUCTION_MANIFEST)
    bulk = _read_json(case_dir / "bulk_state_v9_11_summary.json")
    complete, extension_um = completion_status(case_dir, target_um)
    if returncode != 0 or production.get("run_completed_without_exception") is not True:
        status_name = "failed"
    elif complete:
        status_name = "complete"
    else:
        status_name = "right_censored"
    rc = extract_r_curve(case_dir, T_K)
    row = {
        "schema": "barrier_only_case_status_v10_0_5_13",
        "point_release": POINT_RELEASE,
        "option_key": option.option_key,
        "candidate_id": option.candidate_id,
        "barrier_fingerprint_sha256": option.barrier_fingerprint_sha256,
        "T_K": int(T_K),
        "status": status_name,
        "returncode": int(returncode),
        "reused": bool(reused),
        "target_completed": bool(complete),
        "final_extension_um": extension_um,
        "target_extension_um": float(target_um),
        "K_FP_MPa_sqrt_m": _first_passage(case_dir),
        "n_growth_events": int(len(rc)),
        "candidate_state_fields_applied": production.get("candidate_state_fields_applied"),
        "mesh_refinement_verified": production.get("mesh_refinement_runtime", {}).get("actual_radius_verified"),
        "bulk_explicit_mobile_retained_state": bulk.get("bulk_explicit_mobile_retained_state"),
        "bulk_state_update_calls": bulk.get("bulk_state_update_calls"),
        "case_dir": str(case_dir),
        "completed_utc": _utc_now(),
    }
    _write_json(case_dir / STATUS_FILE, row)
    return row


def _run_case(py: str, args, root: Path, option_key: str, T_K: int, target_um: float):
    case_dir = root / option_key / f"T{int(T_K):04d}"
    if args.skip_finished and _case_is_complete(case_dir, option_key, target_um):
        with _PRINT_LOCK:
            print(f"SKIP complete: {option_key}/T{T_K}K target={target_um:g}um", flush=True)
        return _summarize(case_dir, option_key, T_K, target_um, 0, True)
    if case_dir.exists():
        _archive_and_clean(case_dir, root, option_key, T_K)
    case_dir.mkdir(parents=True, exist_ok=True)
    cmd = _build_command(py, args, option_key, T_K, target_um, case_dir)
    (case_dir / "command.txt").write_text(shlex.join(cmd) + "\n")
    _write_json(
        case_dir / STATUS_FILE,
        {
            "schema": "barrier_only_case_status_v10_0_5_13",
            "point_release": POINT_RELEASE,
            "option_key": option_key,
            "T_K": int(T_K),
            "status": "running",
            "target_extension_um": float(target_um),
            "started_utc": _utc_now(),
            "command": cmd,
        },
    )
    with _PRINT_LOCK:
        print(f"START: {option_key}/T{T_K}K -> {case_dir}", flush=True)
    log_path = case_dir / "run.log"
    with log_path.open("w", buffering=1) as log:
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            env=_case_environment(args, target_um),
        )
        assert process.stdout is not None
        for line in process.stdout:
            log.write(line)
            with _PRINT_LOCK:
                print(f"[{option_key}/T{T_K}K] {line}", end="", flush=True)
        returncode = process.wait()
    row = _summarize(case_dir, option_key, T_K, target_um, returncode, False)
    with _PRINT_LOCK:
        print(
            f"FINISHED: {option_key}/T{T_K}K status={row['status']} "
            f"extension={row['final_extension_um']}um rc={returncode}",
            flush=True,
        )
    return row


def _preflight(py: str, run_tests: bool):
    compile_cmd = [
        py,
        "-m",
        "py_compile",
        "arrhenius_fracture/barrier_only_response_registry_v100513.py",
        "arrhenius_fracture/mode_i_first_passage_v10_0_5_13_barrier_only.py",
        "run_v10_0_5_13_barrier_only_monotonic.py",
    ]
    cp = subprocess.run(compile_cmd, text=True)
    if cp.returncode != 0:
        raise SystemExit(f"compile failed: {shlex.join(compile_cmd)}")
    if run_tests:
        test_cmd = [
            py,
            "-m",
            "pytest",
            "-q",
            "tests/test_v100513_barrier_only.py",
            "tests/test_v1005123_phase_c_repairs.py",
        ]
        cp = subprocess.run(test_cmd, text=True)
        if cp.returncode != 0:
            raise SystemExit(f"tests failed: {shlex.join(test_cmd)}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("smoke", "full"), default="full")
    parser.add_argument("--python-bin", default="")
    parser.add_argument("--conda-env", default="arrhenius-fem-czm")
    parser.add_argument("--registry", type=Path, default=None)
    parser.add_argument("--outroot", type=Path, required=True)
    parser.add_argument("--options", default="")
    parser.add_argument("--temperatures", default="")
    parser.add_argument("--target-extension-um", type=float, default=None)
    parser.add_argument("--max-jobs", type=int, default=2)
    parser.add_argument("--skip-finished", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--run-tests", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--steps", type=int, default=50000)
    parser.add_argument("--nx", type=int, default=36)
    parser.add_argument("--ny", type=int, default=72)
    parser.add_argument("--tip-h-fine", type=float, default=2.5e-6)
    parser.add_argument("--tip-ratio", type=float, default=1.15)
    parser.add_argument("--dU", type=float, default=2.0e-7)
    parser.add_argument("--dt", type=float, default=8.4)
    parser.add_argument("--n-stagger", type=int, default=2)
    parser.add_argument("--print-every", type=int, default=25)
    parser.add_argument("--adaptive-event-target", type=float, default=0.15)
    parser.add_argument("--da-um", type=float, default=5.0)
    parser.add_argument("--theta-deg", type=float, default=45.0)
    parser.add_argument("--tip-refinement-radius-um", type=float, default=330.0)
    parser.add_argument("--cluster-J-outer-um", type=float, default=240.0)
    parser.add_argument("--local-J-outer-um", type=float, default=100.0)
    parser.add_argument("--save-snapshots", type=int, default=3)
    parser.add_argument("--snapshot-cols", type=int, default=3)
    parser.add_argument("--snapshot-interval-um", type=float, default=50.0)
    parser.add_argument("--min-triangle-quality", type=float, default=0.035)
    parser.add_argument("--min-child-area-ratio", type=float, default=0.08)
    parser.add_argument("--max-tip-h-over-da", type=float, default=0.75)
    parser.add_argument("--max-identical-geometry-vetoes", type=int, default=12)
    args = parser.parse_args()

    py = _resolve_python(args)
    options, temperatures, target_um = _selected_matrix(args)
    root = args.outroot.expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    _preflight(py, args.run_tests)

    plan = {
        "schema": "barrier_only_campaign_plan_v10_0_5_13",
        "point_release": POINT_RELEASE,
        "created_utc": _utc_now(),
        "python": py,
        "options": options,
        "temperatures_K": temperatures,
        "target_extension_um": target_um,
        "max_jobs": args.max_jobs,
        "two_d_state_policy": TWO_D_STATE_POLICY,
        "case_count": len(options) * len(temperatures),
    }
    _write_json(root / "barrier_only_campaign_plan_v10_0_5_13.json", plan)
    print(json.dumps(plan, indent=2, default=str), flush=True)

    rows = []
    jobs = [(option, T) for option in options for T in temperatures]
    with ThreadPoolExecutor(max_workers=max(int(args.max_jobs), 1)) as pool:
        future_map = {
            pool.submit(_run_case, py, args, root, option, T, target_um): (option, T)
            for option, T in jobs
        }
        for future in as_completed(future_map):
            option, T = future_map[future]
            try:
                row = future.result()
            except BaseException as exc:
                row = {
                    "option_key": option,
                    "T_K": T,
                    "status": "controller_exception",
                    "runtime_error_type": type(exc).__name__,
                    "runtime_error": str(exc),
                    "target_extension_um": target_um,
                }
                with _PRINT_LOCK:
                    print(f"FAILED: {option}/T{T}K: {type(exc).__name__}: {exc}", flush=True)
            rows.append(row)
            _write_json(
                root / CAMPAIGN_STATUS,
                {
                    **plan,
                    "updated_utc": _utc_now(),
                    "completed_cases": len(rows),
                    "remaining_cases": len(jobs) - len(rows),
                    "rows": sorted(rows, key=lambda x: (x.get("option_key", ""), x.get("T_K", 0))),
                },
            )

    frame = pd.DataFrame(rows).sort_values(["option_key", "T_K"])
    frame.to_csv(root / "barrier_only_campaign_summary_v10_0_5_13.csv", index=False)
    _write_json(
        root / "barrier_only_campaign_summary_v10_0_5_13.json",
        {**plan, "completed_utc": _utc_now(), "rows": frame.to_dict(orient="records")},
    )
    failures = int((~frame["status"].isin(["complete"])).sum()) if not frame.empty else len(jobs)
    print(frame.to_string(index=False), flush=True)
    print(f"v10.0.5.13 barrier-only runner finished: failures={failures}", flush=True)
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
