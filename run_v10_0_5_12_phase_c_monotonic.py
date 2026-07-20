#!/usr/bin/env python3
"""Run the deterministic four-option Phase-C monotonic FEM/CZM matrix."""
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

import numpy as np
import pandas as pd

from arrhenius_fracture.mpz_response_registry_v100512 import (
    PARAMETER_SOURCE,
    PRIMARY_OPTION_KEYS,
    default_registry_path,
    load_option,
    normalize_option_key,
)
from run_four_class_exp_floor_czm_500um_sweep import (
    completion_status,
    extract_r_curve,
    plot_case_r_curve,
)

POINT_RELEASE = "10.0.5.12"
TEMPERATURES_FULL = tuple(range(300, 1201, 100))
TEMPERATURES_ANCHOR = (300, 700, 900, 1200)
PRODUCTION_MANIFEST = "phase_c_production_manifest_v10_0_5_12.json"
COMPLETION_MANIFEST = "run_completion_v10_0_5_2.json"


def values(text, cast=str):
    return [cast(x) for x in str(text).replace(",", " ").split() if x]


def fstr(value):
    value = float(value)
    return "inf" if math.isinf(value) else f"{value:.16g}"


def read_json(path: Path):
    if not path.is_file():
        return {}
    try:
        data = json.loads(path.read_text())
    except Exception:
        return {}
    if isinstance(data, list):
        return data[0] if data and isinstance(data[0], dict) else {}
    return data if isinstance(data, dict) else {}


def resolve_python(args):
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
    lines = [x.strip() for x in cp.stdout.splitlines() if x.strip()]
    if not lines:
        raise SystemExit(f"Conda environment {args.conda_env!r} returned no Python path")
    return lines[-1]


def default_matrix(mode):
    if mode == "smoke":
        return ("dbtt_primary",), (700,), 50.0
    if mode == "anchors":
        return PRIMARY_OPTION_KEYS, TEMPERATURES_ANCHOR, 500.0
    if mode == "full":
        return PRIMARY_OPTION_KEYS, TEMPERATURES_FULL, 500.0
    raise ValueError(mode)


def selected_matrix(args):
    options0, temperatures0, target0 = default_matrix(args.mode)
    options = [normalize_option_key(x) for x in values(args.options)] if args.options else list(options0)
    temperatures = [int(round(x)) for x in values(args.temperatures, float)] if args.temperatures else list(temperatures0)
    target = float(args.target_extension_um if args.target_extension_um is not None else target0)
    options = list(dict.fromkeys(options))
    temperatures = list(dict.fromkeys(temperatures))
    if not options or not temperatures or target <= 0 or any(T <= 0 for T in temperatures):
        raise SystemExit("Phase-C matrix requires positive temperatures and target extension")
    return options, temperatures, target


def build_command(py, args, option_key, T_K, target_um, case_dir):
    option = load_option(option_key, args.registry)
    cmd = [
        py, "-m", "arrhenius_fracture.mode_i_first_passage_v10_0_5_12_phase_c",
        "--phase-c-option", option.option_key,
        "--tip-refinement-radius-um", fstr(args.tip_refinement_radius_um),
        "--selected-cluster-J-outer-um", fstr(args.cluster_J_outer_um),
        "--local-J-outer-um", fstr(args.local_J_outer_um),
        "--v10-material-source", PARAMETER_SOURCE,
        "--czm-opening-coupling", "clock_linear",
        "--mode", "2d", "--bulk-plasticity-mode", "tip_only",
        "--temperatures", str(int(T_K)), "--steps", str(args.steps),
        "--nx", str(args.nx), "--ny", str(args.ny),
        "--tip-h-fine", fstr(args.tip_h_fine), "--tip-ratio", fstr(args.tip_ratio),
        "--dU", fstr(args.dU), "--dt", fstr(args.dt),
        "--n-stagger", str(args.n_stagger), "--print-every", str(args.print_every),
        "--adaptive-events", "--adaptive-event-target", fstr(args.adaptive_event_target),
        "--adaptive-min-frac", "1e-8", "--adaptive-grow", "4",
        "--da-phys", fstr(args.da_um * 1e-6),
        "--target-crack-extension-um", fstr(target_um),
        "--crystal-aniso", "--crystal-compete", "--crystal-theta-deg", fstr(args.theta_deg),
        "--crystal-C11", "523e9", "--crystal-C12", "203e9", "--crystal-C44", "160e9",
        "--cleave-gamma-aniso", "0.3", "--crystal-material", "w",
        "--max-fronts", "1", "--crack-backend", "adaptive_czm",
        "--czm-max-angle-error-deg", "35", "--j-decomposition", "cluster",
        "--rJ-cluster", fstr(args.cluster_J_outer_um * 1e-6),
        "--rJ-outer", fstr(args.local_J_outer_um * 1e-6),
        "--mpz-length-um", fstr(option.mpz_length_um),
        "--mpz-n-bins", str(option.mpz_n_bins),
        "--save-snapshots", str(args.save_snapshots),
        "--snapshot-cols", str(args.snapshot_cols),
        "--snapshot-by-crack-extension-um", fstr(args.snapshot_interval_um),
        "--out", str(case_dir),
    ]
    if args.registry:
        cmd[3:3] = ["--phase-c-registry", str(args.registry.resolve())]
    return cmd


def case_environment(args, target_um):
    env = os.environ.copy()
    env.update({
        "ARRHENIUS_COMMITTED_TARGET_EXTENSION_UM": fstr(target_um),
        "ARRHENIUS_PREFINED_MODE_I_CORRIDOR": "1",
        "ARRHENIUS_MIN_ACCEPTED_TRIANGLE_QUALITY": fstr(args.min_triangle_quality),
        "ARRHENIUS_MIN_ACCEPTED_CHILD_AREA_RATIO": fstr(args.min_child_area_ratio),
        "ARRHENIUS_MAX_TIP_H_OVER_DA": fstr(args.max_tip_h_over_da),
        "ARRHENIUS_MAX_TRIAL_DAMAGE_CHANGE": fstr(args.max_trial_damage_change),
        "ARRHENIUS_MIN_TRIAL_RETRY_DT_S": fstr(args.min_trial_retry_dt_s),
        "ARRHENIUS_MAX_TRIAL_RETRIES": str(args.max_trial_retries),
        "ARRHENIUS_MAX_ACCEPTED_SUBSTEPS_PER_INTERVAL": str(args.max_accepted_substeps_per_interval),
        "ARRHENIUS_TENSOR_DRIVE_PROBE_RADIUS_M": fstr(args.tensor_probe_radius_m),
        "ARRHENIUS_TENSOR_DRIVE_SECTOR_HALF_ANGLE_DEG": fstr(args.tensor_sector_half_angle_deg),
        "ARRHENIUS_TENSOR_DRIVE_MIN_ELEMENTS": str(args.tensor_min_elements),
        "ARRHENIUS_EVENT_STATISTICS": "mean_field",
        "ARRHENIUS_STOCHASTIC_EMISSION": "0",
        "ARRHENIUS_VHCF_FEM_CACHE": "0",
        "PYTHONUNBUFFERED": "1",
    })
    return env


def first_passage(case_dir):
    data = read_json(case_dir / "anisotropic_calibrated_tip_first_passage_summary.json")
    for key in ("KJ_reference_first_MPa_sqrt_m", "Kc_first_existing_MPa_sqrt_m"):
        try:
            value = float(data[key])
        except (KeyError, TypeError, ValueError):
            continue
        if math.isfinite(value):
            return value
    return None


def median_window(rc, lo, hi):
    if rc.empty:
        return math.nan
    values_ = rc.loc[(rc.crack_extension_um >= lo) & (rc.crack_extension_um <= hi), "KJ_MPa_sqrt_m"].to_numpy(float)
    return float(np.nanmedian(values_)) if len(values_) else math.nan


def metrics(rc, target_um, Kfp):
    if rc.empty:
        return {"K_FP_MPa_sqrt_m": Kfp, "K0_0_50um_median": math.nan,
                "Kmid_200_300um_median": math.nan, "Klate_350_target_median": math.nan,
                "delta_KR_late_minus_K0": math.nan, "normalized_R_curve_area_MPa_sqrt_m": math.nan,
                "n_growth_events": 0}
    K0 = median_window(rc, 0, min(50.5, target_um + 0.5))
    Kmid = median_window(rc, min(200, 0.4 * target_um), min(300.5, 0.6 * target_um + 0.5))
    Klate = median_window(rc, min(350, max(0.7 * target_um, 50)), target_um + 0.5)
    ordered = rc.sort_values("crack_extension_um")
    x = ordered.crack_extension_um.to_numpy(float)
    y = ordered.KJ_MPa_sqrt_m.to_numpy(float)
    good = np.isfinite(x) & np.isfinite(y) & (x >= 0) & (x <= target_um + 0.5)
    x, y = x[good], y[good]
    area = math.nan
    if len(x) >= 2 and np.nanmax(x) >= target_um - 0.5:
        x, idx = np.unique(x, return_index=True)
        y = y[idx]
        if x[0] > 0:
            x = np.r_[0.0, x]
            y = np.r_[float(Kfp) if Kfp is not None else y[0], y]
        yt = float(np.interp(target_um, x, y))
        keep = x < target_um
        area = float(np.trapz(np.r_[y[keep], yt], np.r_[x[keep], target_um]) / target_um)
    return {
        "K_FP_MPa_sqrt_m": Kfp,
        "K0_0_50um_median": K0,
        "Kmid_200_300um_median": Kmid,
        "Klate_350_target_median": Klate,
        "delta_KR_late_minus_K0": Klate - K0 if math.isfinite(Klate) and math.isfinite(K0) else math.nan,
        "normalized_R_curve_area_MPa_sqrt_m": area,
        "n_growth_events": int(len(rc)),
    }


def summarize(case_dir, option_key, T_K, target_um, returncode, reused, registry=None):
    option = load_option(option_key, registry)
    complete, extension_um = completion_status(case_dir, target_um)
    lifecycle = read_json(case_dir / COMPLETION_MANIFEST)
    production = read_json(case_dir / PRODUCTION_MANIFEST)
    lifecycle_ok = lifecycle.get("run_completed_without_exception") is True
    production_ok = (
        production.get("run_completed_without_exception") is True
        and production.get("option", {}).get("candidate_id") == option.candidate_id
        and production.get("option", {}).get("fingerprint_sha256") == option.fingerprint_sha256
        and production.get("mesh_refinement_runtime", {}).get("actual_radius_verified") is True
    )
    if returncode != 0 or not lifecycle_ok or not production_ok:
        status = "failed"
    elif complete:
        status = "complete"
    else:
        status = "right_censored"
    rc = extract_r_curve(case_dir, T_K)
    row = {
        "option_key": option.option_key, "candidate_id": option.candidate_id,
        "parameter_fingerprint_sha256": option.fingerprint_sha256,
        "canonical_class": option.canonical_class, "T_K": int(T_K), "status": status,
        "returncode": int(returncode), "reused": bool(reused), "target_completed": bool(complete),
        "final_extension_um": extension_um, "target_extension_um": float(target_um),
        "completion_manifest_passed": lifecycle_ok, "production_manifest_passed": production_ok,
        "mpz_length_um": option.mpz_length_um, "mpz_n_bins": option.mpz_n_bins,
        "case_dir": str(case_dir), **metrics(rc, target_um, first_passage(case_dir)),
    }
    pd.DataFrame([row]).to_csv(case_dir / "phase_c_case_summary.csv", index=False)
    (case_dir / "phase_c_case_summary.json").write_text(json.dumps(row, indent=2, default=str))
    plot_case_r_curve(rc, case_dir, option.option_key, T_K, target_um)
    if status == "complete":
        (case_dir / ".phase_c_complete").touch()
    return row


def run_case(py, args, root, option_key, T_K, target_um):
    case_dir = root / option_key / f"T{int(T_K):04d}"
    case_dir.mkdir(parents=True, exist_ok=True)
    if args.skip_existing and (case_dir / ".phase_c_complete").is_file():
        print(f"REUSE {option_key:16s} T={T_K:4d} K")
        return summarize(case_dir, option_key, T_K, target_um, 0, True, args.registry)
    option = load_option(option_key, args.registry)
    cmd = build_command(py, args, option_key, T_K, target_um, case_dir)
    (case_dir / "command.txt").write_text(shlex.join(cmd) + "\n")
    (case_dir / "phase_c_case_input.json").write_text(json.dumps({
        "schema": "phase_c_case_input_v10_0_5_12", "point_release": POINT_RELEASE,
        "option": option.audit_payload(), "T_K": int(T_K), "target_extension_um": target_um,
        "command": cmd}, indent=2, default=str))
    print(f"START {option_key:16s} T={T_K:4d} K -> {case_dir}")
    if args.dry_run:
        return {"option_key": option_key, "candidate_id": option.candidate_id, "T_K": int(T_K),
                "status": "dry_run", "target_extension_um": target_um, "case_dir": str(case_dir)}
    with (case_dir / "run.log").open("w") as log:
        cp = subprocess.run(cmd, env=case_environment(args, target_um), stdout=log,
                            stderr=subprocess.STDOUT, text=True)
    row = summarize(case_dir, option_key, T_K, target_um, cp.returncode, False, args.registry)
    print(f"{row['status'].upper():14s} {option_key:16s} T={T_K:4d} K ext={row.get('final_extension_um')}")
    return row


def preflight(py, run_tests):
    commands = [
        [py, "-m", "py_compile", "arrhenius_fracture/mpz_response_registry_v100512.py",
         "arrhenius_fracture/mode_i_first_passage_v10_0_5_12_phase_c.py",
         "run_v10_0_5_12_phase_c_monotonic.py"],
        [py, "-c", "from arrhenius_fracture.mpz_response_registry_v100512 import load_registry; "
         "r=load_registry(); assert len(r)==4; print({k:v.candidate_id for k,v in r.items()})"],
    ]
    if run_tests:
        commands.append([py, "-m", "pytest", "-q", "tests/test_v100512_phase_c.py"])
    for command in commands:
        cp = subprocess.run(command, text=True)
        if cp.returncode != 0:
            raise SystemExit(f"Phase-C preflight failed: {shlex.join(command)}")


def provenance(py):
    def run(command):
        cp = subprocess.run(command, text=True, capture_output=True)
        return cp.stdout.strip() if cp.returncode == 0 else None
    return {
        "recorded_utc": datetime.now(timezone.utc).isoformat(),
        "python_executable": str(Path(py).resolve()),
        "conda_default_env": os.environ.get("CONDA_DEFAULT_ENV"),
        "git_branch": run(["git", "branch", "--show-current"]),
        "git_commit": run(["git", "rev-parse", "HEAD"]),
        "git_status_porcelain": run(["git", "status", "--porcelain"]),
        "editable_package_path": run([py, "-c", "import pathlib, arrhenius_fracture; "
                                      "print(pathlib.Path(arrhenius_fracture.__file__).resolve())"]),
    }


def parser():
    p = argparse.ArgumentParser()
    p.add_argument("--mode", choices=("smoke", "anchors", "full"), default="anchors")
    p.add_argument("--options"); p.add_argument("--temperatures")
    p.add_argument("--target-extension-um", type=float); p.add_argument("--registry", type=Path)
    p.add_argument("--outroot", type=Path, default=Path("runs/v10_0_5_12_phase_c_500um_theta45_v1"))
    p.add_argument("--conda-env", default="arrhenius-fem-czm"); p.add_argument("--python-bin")
    p.add_argument("--max-jobs", type=int, default=2); p.add_argument("--steps", type=int, default=50000)
    p.add_argument("--nx", type=int, default=36); p.add_argument("--ny", type=int, default=72)
    p.add_argument("--tip-h-fine", type=float, default=2.5e-6); p.add_argument("--tip-ratio", type=float, default=1.15)
    p.add_argument("--dU", type=float, default=2e-7); p.add_argument("--dt", type=float, default=8.4)
    p.add_argument("--n-stagger", type=int, default=2); p.add_argument("--print-every", type=int, default=25)
    p.add_argument("--adaptive-event-target", type=float, default=0.15); p.add_argument("--da-um", type=float, default=5)
    p.add_argument("--theta-deg", type=float, default=45); p.add_argument("--tip-refinement-radius-um", type=float, default=330)
    p.add_argument("--cluster-J-outer-um", type=float, default=240); p.add_argument("--local-J-outer-um", type=float, default=100)
    p.add_argument("--save-snapshots", type=int, default=11); p.add_argument("--snapshot-cols", type=int, default=6)
    p.add_argument("--snapshot-interval-um", type=float, default=50); p.add_argument("--min-triangle-quality", type=float, default=0.035)
    p.add_argument("--min-child-area-ratio", type=float, default=0.08); p.add_argument("--max-tip-h-over-da", type=float, default=0.75)
    p.add_argument("--max-trial-damage-change", type=float, default=0.02); p.add_argument("--min-trial-retry-dt-s", type=float, default=1e-18)
    p.add_argument("--max-trial-retries", type=int, default=64); p.add_argument("--max-accepted-substeps-per-interval", type=int, default=10000)
    p.add_argument("--tensor-probe-radius-m", type=float, default=1e-5); p.add_argument("--tensor-sector-half-angle-deg", type=float, default=25)
    p.add_argument("--tensor-min-elements", type=int, default=3)
    p.add_argument("--skip-existing", action=argparse.BooleanOptionalAction, default=True)
    p.add_argument("--run-tests", action=argparse.BooleanOptionalAction, default=True)
    p.add_argument("--dry-run", action="store_true"); p.add_argument("--fail-fast", action="store_true")
    return p


def main():
    args = parser().parse_args()
    if args.max_jobs < 1:
        raise SystemExit("--max-jobs must be positive")
    options, temperatures, target = selected_matrix(args)
    py = resolve_python(args)
    preflight(py, args.run_tests)
    root = args.outroot.resolve(); root.mkdir(parents=True, exist_ok=True)
    source_registry = args.registry.resolve() if args.registry else default_registry_path()
    registry_snapshot = root / "phase_c_registry_snapshot.json"; shutil.copy2(source_registry, registry_snapshot)
    prov = provenance(py); (root / "phase_c_provenance.json").write_text(json.dumps(prov, indent=2))
    matrix = []
    for option_key in options:
        option = load_option(option_key, args.registry)
        for T_K in temperatures:
            matrix.append({"option_key": option.option_key, "candidate_id": option.candidate_id,
                           "parameter_fingerprint_sha256": option.fingerprint_sha256,
                           "canonical_class": option.canonical_class, "T_K": int(T_K),
                           "target_extension_um": target, "mpz_length_um": option.mpz_length_um,
                           "mpz_n_bins": option.mpz_n_bins,
                           "case_dir": str(root / option.option_key / f"T{int(T_K):04d}")})
    pd.DataFrame(matrix).to_csv(root / "phase_c_matrix.csv", index=False)
    (root / "phase_c_matrix.json").write_text(json.dumps(matrix, indent=2))
    campaign = {"schema": "phase_c_campaign_v10_0_5_12", "point_release": POINT_RELEASE,
                "created_utc": datetime.now(timezone.utc).isoformat(), "mode": args.mode,
                "options": options, "temperatures_K": temperatures, "target_extension_um": target,
                "case_count": len(matrix), "max_jobs": args.max_jobs, "python": py,
                "registry_snapshot": str(registry_snapshot), "provenance": prov,
                "event_statistics": "mean_field", "stochastic_emission": False,
                "branching_enabled": False, "max_fronts": 1,
                "tip_refinement_radius_um": args.tip_refinement_radius_um,
                "selected_cluster_J_outer_um": args.cluster_J_outer_um,
                "local_J_outer_um": args.local_J_outer_um, "argv": sys.argv}
    (root / "phase_c_campaign.json").write_text(json.dumps(campaign, indent=2, default=str))
    rows, failures = [], []
    with ThreadPoolExecutor(max_workers=args.max_jobs) as pool:
        futures = {pool.submit(run_case, py, args, root, item["option_key"], item["T_K"], target): item for item in matrix}
        for future in as_completed(futures):
            item = futures[future]
            try:
                row = future.result()
            except BaseException as exc:
                failures.append(exc); row = {**item, "status": "runner_exception", "runtime_error": str(exc)}
                print(f"RUNNER_EXCEPTION {item['option_key']} T={item['T_K']} K: {exc}", file=sys.stderr)
                if args.fail_fast:
                    for pending in futures: pending.cancel()
            rows.append(row)
            pd.DataFrame(rows).sort_values(["option_key", "T_K"]).to_csv(root / "phase_c_summary.partial.csv", index=False)
            (root / "phase_c_summary.partial.json").write_text(json.dumps(rows, indent=2, default=str))
    frame = pd.DataFrame(rows).sort_values(["option_key", "T_K"])
    frame.to_csv(root / "phase_c_summary.csv", index=False)
    (root / "phase_c_summary.json").write_text(json.dumps(rows, indent=2, default=str))
    final = {**campaign, "completed_utc": datetime.now(timezone.utc).isoformat(),
             "status_counts": frame["status"].value_counts(dropna=False).to_dict()}
    (root / "phase_c_completion.json").write_text(json.dumps(final, indent=2, default=str))
    print(frame.to_string(index=False)); print("wrote", root / "phase_c_summary.csv")
    if failures or (not args.dry_run and any(frame.status != "complete")):
        raise SystemExit("one or more Phase-C cases failed or were incomplete; inspect phase_c_summary.csv")


if __name__ == "__main__":
    main()
