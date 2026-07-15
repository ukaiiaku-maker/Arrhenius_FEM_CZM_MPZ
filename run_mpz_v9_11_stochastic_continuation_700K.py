#!/usr/bin/env python3
"""Run the three v9.11 classes in both bulk modes with stochastic event clocks.

Each seed is a reproducible realization of the same class parameterization. The
solver uses exponential integrated-hazard thresholds for cleavage, binomial
finite-site emission, and event-to-event reload continuation. Raw cohesive-edge
events and cascade-aware load events are both retained.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import traceback

import numpy as np
import pandas as pd

from arrhenius_fracture.bulk_state_v911 import VALID_BULK_MODES
from arrhenius_fracture.mpz_parameterization_v911 import normalize_class_name
from arrhenius_fracture.rcurve_postprocess_v911 import write_cascade_aware_outputs

CLASSES = ("ceramic", "weakT", "DBTT")
_RAW_PROPAGATION_KEYS = (
    "Kprop_200_500um_median",
    "Kprop_200_500um_mean",
    "Kprop_200_500um_p10",
    "Kprop_200_500um_p90",
    "delta_KR_median_minus_init",
    "n_growth_events",
)


def values(text: str, cast=int):
    return [cast(x) for x in str(text).replace(",", " ").split() if x]


def read_json(path: Path) -> dict:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text())
    except Exception:
        return {}


def read_first_csv_row(path: Path) -> dict:
    if not path.exists() or path.stat().st_size == 0:
        return {}
    try:
        frame = pd.read_csv(path)
    except (pd.errors.EmptyDataError, OSError, ValueError):
        return {}
    return frame.iloc[0].to_dict() if not frame.empty else {}


def empty_continuation_metrics() -> dict:
    return {
        "Kload_200_500um_median": np.nan,
        "Kload_200_500um_mean": np.nan,
        "Kload_200_500um_p10": np.nan,
        "Kload_200_500um_p90": np.nan,
        "delta_Kload_median_minus_init": np.nan,
    }


def continuation_metrics(case_dir: Path, K_init) -> dict:
    path = case_dir / "R_curve_load_events_clustered.csv"
    if not path.exists() or path.stat().st_size == 0:
        return empty_continuation_metrics()
    try:
        frame = pd.read_csv(path)
    except (pd.errors.EmptyDataError, OSError, ValueError):
        return empty_continuation_metrics()
    required = {"crack_extension_um", "KJ_onset_MPa_sqrt_m"}
    if frame.empty or not required.issubset(frame.columns):
        vals = np.asarray([], dtype=float)
    else:
        x = pd.to_numeric(frame["crack_extension_um"], errors="coerce")
        k = pd.to_numeric(frame["KJ_onset_MPa_sqrt_m"], errors="coerce")
        vals = k[(x >= 200.0) & (x <= 500.5)].dropna().to_numpy(float)
    if vals.size:
        med = float(np.median(vals))
        mean = float(np.mean(vals))
        p10 = float(np.percentile(vals, 10))
        p90 = float(np.percentile(vals, 90))
    else:
        med = mean = p10 = p90 = np.nan
    try:
        init = float(K_init)
    except (TypeError, ValueError):
        init = np.nan
    return {
        "Kload_200_500um_median": med,
        "Kload_200_500um_mean": mean,
        "Kload_200_500um_p10": p10,
        "Kload_200_500um_p90": p90,
        "delta_Kload_median_minus_init": (
            med - init if np.isfinite(med) and np.isfinite(init) else np.nan
        ),
    }


def summary_matches_case(row: dict, class_name: str, mode: str, T_K: float) -> bool:
    if not row:
        return False
    try:
        temp_match = int(round(float(row.get("T_K")))) == int(round(float(T_K)))
    except (TypeError, ValueError):
        temp_match = False
    return (
        str(row.get("class", "")) == class_name
        and str(row.get("bulk_plasticity_mode", "")) == mode
        and temp_match
    )


def solver_outputs_complete_enough_to_reuse(case_dir: Path, T_K: float) -> bool:
    tag = int(round(float(T_K)))
    return all((case_dir / name).exists() for name in (
        f"steps_{tag:04d}K.csv",
        "anisotropic_calibrated_tip_first_passage_summary.json",
        "bulk_state_v9_11_summary.json",
        "run.log",
    ))


def run_case(args, cls: str, mode: str, seed: int, root: Path) -> dict:
    class_name = normalize_class_name(cls)
    run_root = root / f"seed_{seed}" / mode
    case_dir = run_root / class_name / f"T{int(round(args.T_K))}_th{args.crystal_theta_deg:g}"
    case_dir.mkdir(parents=True, exist_ok=True)
    log_dir = root / "matrix_logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log = log_dir / f"{class_name}_{mode}_seed{seed}_{int(args.T_K)}K.log"
    temp_summary_path = run_root / "rcurve_temperature_summary.csv"

    cmd = [
        sys.executable,
        "run_mpz_v9_11_mode_i_rcurve_3T.py",
        "--parameter-root", str(args.parameter_root),
        "--material-class", class_name,
        "--bulk-plasticity-mode", mode,
        "--temperatures", str(int(round(args.T_K))),
        "--outroot", str(run_root),
        "--target-extension-um", str(args.target_extension_um),
        "--steps", str(args.steps),
        "--nx", str(args.nx),
        "--ny", str(args.ny),
        "--tip-h-fine", str(args.tip_h_fine),
        "--tip-ratio", str(args.tip_ratio),
        "--dU", str(args.dU),
        "--dt", str(args.dt),
        "--n-stagger", str(args.n_stagger),
        "--print-every", str(args.print_every),
        "--adaptive-event-target", str(args.adaptive_event_target),
        "--da-phys-um", str(args.da_phys_um),
        "--mpz-length-um", str(args.mpz_length_um),
        "--mpz-n-bins", str(args.mpz_n_bins),
        "--crystal-theta-deg", str(args.crystal_theta_deg),
        "--save-snapshots", str(args.save_snapshots),
        "--snapshot-cols", str(args.snapshot_cols),
        "--snapshot-by-extension-um", str(args.snapshot_by_extension_um),
        "--no-make-solver-plots",
    ]
    cmd.append("--skip-existing" if args.skip_existing else "--no-skip-existing")

    env = os.environ.copy()
    env.update({
        "ARRHENIUS_EVENT_STATISTICS": "stochastic",
        "ARRHENIUS_STOCHASTIC_EMISSION": "1" if args.stochastic_emission else "0",
        "ARRHENIUS_STOCHASTIC_SEED": str(seed),
        "ARRHENIUS_PROPAGATION_CONTROL": "event_reload",
        "ARRHENIUS_RELOAD_RELATIVE_U": str(args.reload_relative_U),
        "ARRHENIUS_RELOAD_ABSOLUTE_U_M": str(args.reload_absolute_U_m),
        "ARRHENIUS_RELOAD_RELATIVE_K": str(args.reload_relative_K),
        "ARRHENIUS_RELOAD_ABSOLUTE_K_PA_SQRT_M": str(args.reload_absolute_K),
    })
    (case_dir / "stochastic_continuation_run.json").write_text(json.dumps({
        "class": class_name,
        "bulk_plasticity_mode": mode,
        "seed": seed,
        "event_statistics": "stochastic",
        "stochastic_emission": bool(args.stochastic_emission),
        "propagation_control": "event_reload",
        "reload_relative_U": args.reload_relative_U,
        "reload_absolute_U_m": args.reload_absolute_U_m,
        "reload_relative_K": args.reload_relative_K,
        "reload_absolute_K_Pa_sqrt_m": args.reload_absolute_K,
        "command": cmd,
    }, indent=2))

    existing_summary = read_first_csv_row(temp_summary_path)
    reuse_solver = bool(
        args.skip_existing
        and summary_matches_case(existing_summary, class_name, mode, args.T_K)
        and solver_outputs_complete_enough_to_reuse(case_dir, args.T_K)
    )

    print(f"START {class_name:7s} {mode:18s} seed={seed}")
    if reuse_solver:
        subprocess_returncode = int(existing_summary.get("returncode", 0) or 0)
        solver_reused = True
        print(f"REUSE {class_name:7s} {mode:18s} seed={seed} existing solver outputs")
    else:
        with log.open("w") as fp:
            cp = subprocess.run(cmd, env=env, stdout=fp, stderr=subprocess.STDOUT, text=True)
        subprocess_returncode = int(cp.returncode)
        solver_reused = False

    cascade = write_cascade_aware_outputs(
        case_dir,
        args.T_K,
        relative_load_tolerance=args.cluster_relative_load_tolerance,
        absolute_load_tolerance_m=args.cluster_absolute_load_tolerance_m,
    )
    fp_summary = read_json(case_dir / "anisotropic_calibrated_tip_first_passage_summary.json")
    bulk = read_json(case_dir / "bulk_state_v9_11_summary.json")
    run_summary = read_first_csv_row(temp_summary_path)

    # Preserve legacy metrics explicitly as serialized-topology diagnostics. They
    # are not interpreted as resistance after same-load cascades were discovered.
    for key in _RAW_PROPAGATION_KEYS:
        if key in run_summary:
            run_summary[f"serialized_topology_{key}"] = run_summary.pop(key)
    load_metrics = continuation_metrics(
        case_dir, run_summary.get("K_init_MPa_sqrt_m")
    )

    row = {
        **run_summary,
        "class": class_name,
        "bulk_plasticity_mode": mode,
        "seed": seed,
        "T_K": float(args.T_K),
        "event_statistics": "stochastic",
        "stochastic_emission": bool(args.stochastic_emission),
        "propagation_control": "event_reload",
        "subprocess_returncode": subprocess_returncode,
        "solver_output_reused": solver_reused,
        "case_dir": str(case_dir),
        "matrix_log": str(log),
        "B_target_final": fp_summary.get("B_target_final", fp_summary.get("B_target")),
        "stochastic_event_index_final": fp_summary.get("stochastic_event_index_final"),
        **bulk,
        **cascade,
        **load_metrics,
    }
    (case_dir / "stochastic_continuation_case_summary.json").write_text(
        json.dumps(row, indent=2, default=str)
    )
    pd.DataFrame([row]).to_csv(
        case_dir / "stochastic_continuation_case_summary.csv", index=False
    )
    print(
        f"DONE  {class_name:7s} {mode:18s} seed={seed} rc={subprocess_returncode} "
        f"status={row.get('status')} ext={row.get('final_extension_um')} "
        f"load_events={row.get('n_independent_load_events')}"
    )
    return row


def write_campaign_summary(root: Path, rows: list[dict], final: bool = False) -> None:
    stem = "stochastic_continuation_700K_summary" if final else "stochastic_continuation_700K_summary.partial"
    frame = pd.DataFrame(rows)
    frame.to_csv(root / f"{stem}.csv", index=False)
    (root / f"{stem}.json").write_text(json.dumps(rows, indent=2, default=str))


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--parameter-root", type=Path, default=Path("mpz_v9_11_parameters"))
    p.add_argument("--outroot", type=Path, default=Path("runs/mpz_v9_11_stochastic_continuation_700K_v1"))
    p.add_argument("--seeds", default="1 2 3")
    p.add_argument("--classes", default="ceramic weakT DBTT")
    p.add_argument("--bulk-modes", default="tip_only bulk_same_pt_km")
    p.add_argument("--T-K", type=float, default=700.0)
    p.add_argument("--target-extension-um", type=float, default=500.0)
    p.add_argument("--steps", type=int, default=12000)
    p.add_argument("--nx", type=int, default=36)
    p.add_argument("--ny", type=int, default=72)
    p.add_argument("--tip-h-fine", type=float, default=1.0e-6)
    p.add_argument("--tip-ratio", type=float, default=1.20)
    p.add_argument("--dU", type=float, default=2.0e-7)
    p.add_argument("--dt", type=float, default=8.4)
    p.add_argument("--n-stagger", type=int, default=2)
    p.add_argument("--print-every", type=int, default=25)
    p.add_argument("--adaptive-event-target", type=float, default=0.15)
    p.add_argument("--da-phys-um", type=float, default=5.0)
    p.add_argument("--mpz-length-um", type=float, default=100.0)
    p.add_argument("--mpz-n-bins", type=int, default=200)
    p.add_argument("--crystal-theta-deg", type=float, default=45.0)
    p.add_argument("--save-snapshots", type=int, default=12)
    p.add_argument("--snapshot-cols", type=int, default=4)
    p.add_argument("--snapshot-by-extension-um", type=float, default=50.0)
    p.add_argument("--stochastic-emission", action=argparse.BooleanOptionalAction, default=True)
    p.add_argument("--reload-relative-U", type=float, default=1.0e-4)
    p.add_argument("--reload-absolute-U-m", type=float, default=1.0e-12)
    p.add_argument("--reload-relative-K", type=float, default=1.0e-4)
    p.add_argument("--reload-absolute-K", type=float, default=1.0e3)
    p.add_argument("--cluster-relative-load-tolerance", type=float, default=1.0e-4)
    p.add_argument("--cluster-absolute-load-tolerance-m", type=float, default=1.0e-12)
    p.add_argument("--skip-existing", action=argparse.BooleanOptionalAction, default=True)
    args = p.parse_args()

    seeds = values(args.seeds, int)
    classes = [normalize_class_name(x) for x in values(args.classes, str)]
    modes = values(args.bulk_modes, str)
    bad_modes = sorted(set(modes).difference(VALID_BULK_MODES))
    if bad_modes:
        raise SystemExit(f"unknown bulk modes: {bad_modes}")
    root = args.outroot.resolve()
    root.mkdir(parents=True, exist_ok=True)

    rows: list[dict] = []
    for seed in seeds:
        for cls in classes:
            for mode in modes:
                try:
                    row = run_case(args, cls, mode, seed, root)
                except Exception as exc:
                    error_path = root / "matrix_logs" / f"{normalize_class_name(cls)}_{mode}_seed{seed}_driver_error.log"
                    error_path.parent.mkdir(parents=True, exist_ok=True)
                    error_path.write_text(traceback.format_exc())
                    row = {
                        "class": normalize_class_name(cls),
                        "bulk_plasticity_mode": mode,
                        "seed": seed,
                        "T_K": float(args.T_K),
                        "status": "campaign_driver_error",
                        "campaign_driver_error": f"{type(exc).__name__}: {exc}",
                        "campaign_driver_error_log": str(error_path),
                    }
                    print(f"ERROR {row['class']:7s} {mode:18s} seed={seed}: {exc}", file=sys.stderr)
                rows.append(row)
                write_campaign_summary(root, rows, final=False)

    frame = pd.DataFrame(rows)
    write_campaign_summary(root, rows, final=True)
    print(frame.to_string(index=False))
    print("wrote", root / "stochastic_continuation_700K_summary.csv")


if __name__ == "__main__":
    main()
