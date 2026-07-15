#!/usr/bin/env python3
"""Run long-growth Mode-I v9.11 FEM/CZM R-curves at three temperatures.

This driver uses the already selected moving-reference-frame ceramic, weakT, or
DBTT material manifest.  It does not fit or recalibrate any barrier.  Each case
runs the full anisotropic elastic-plastic FEM, domain-integral K driver,
adaptive-CZM geometry, and v9.11 moving process zone to a requested projected
crack extension.  Branching is disabled for this constitutive transfer gate.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import shlex
import subprocess
import sys

import numpy as np
import pandas as pd

from arrhenius_fracture.mpz_parameterization_v911 import normalize_class_name
from run_four_class_exp_floor_czm_500um_sweep import (
    completion_status,
    make_class_r_curve_overlays,
    process_case_r_curve,
    propagation_metrics,
)


def values(text: str, cast=float):
    return [cast(x) for x in str(text).replace(",", " ").split() if x]


def fs(value: float) -> str:
    value = float(value)
    return "inf" if math.isinf(value) else f"{value:.16g}"


def manifest_for(root: Path, class_name: str) -> Path:
    path = root / class_name / "spatial_promotion_manifest.csv"
    if not path.exists():
        raise FileNotFoundError(path)
    return path


def build_command(py: str, args, class_name: str, manifest: Path, T_K: int, case_dir: Path) -> list[str]:
    cmd = [
        py, "-m", "arrhenius_fracture.mode_i_first_passage_v9_11",
        "--mpz-material-manifest", str(manifest),
        "--mpz-material-class", class_name,
        "--mpz-length-um", fs(args.mpz_length_um),
        "--mpz-n-bins", str(args.mpz_n_bins),
        "--mpz-profile-sector-half-angle-deg", fs(args.mpz_profile_sector_half_angle_deg),
        "--mpz-profile-damage-cutoff", fs(args.mpz_profile_damage_cutoff),
        "--mode", "2d",
        "--nx", str(args.nx), "--ny", str(args.ny),
        "--tip-h-fine", fs(args.tip_h_fine), "--tip-ratio", fs(args.tip_ratio),
        "--dU", fs(args.dU), "--dt", fs(args.dt), "--steps", str(args.steps),
        "--n-stagger", str(args.n_stagger), "--print-every", str(args.print_every),
        "--target-crack-extension-um", fs(args.target_extension_um),
        "--max-fronts", "1",
        "--adaptive-events", "--adaptive-event-target", fs(args.adaptive_event_target),
        "--adaptive-min-frac", "1e-8", "--adaptive-grow", "4",
        "--da-phys", fs(args.da_phys_um * 1.0e-6),
        "--j-decomposition", "cluster",
        "--rJ-cluster", fs(args.rJ_cluster_um * 1.0e-6),
        "--rJ-outer", fs(args.rJ_outer_um * 1.0e-6),
        "--temperatures", str(int(T_K)),
        "--crack-backend", "adaptive_czm",
        "--czm-max-angle-error-deg", "35",
        "--crystal-aniso", "--crystal-compete",
        "--crystal-theta-deg", fs(args.crystal_theta_deg),
        "--crystal-C11", "523e9", "--crystal-C12", "203e9", "--crystal-C44", "160e9",
        "--cleave-gamma-aniso", "0.3", "--crystal-material", "w",
        "--multihit-m", "3", "--multihit-tau", "1e-6",
        "--sigma-cap-GPa", "0",
        "--save-snapshots", str(args.save_snapshots),
        "--snapshot-cols", str(args.snapshot_cols),
        "--out", str(case_dir),
    ]
    if args.snapshot_by_extension_um > 0.0:
        cmd += ["--snapshot-by-crack-extension-um", fs(args.snapshot_by_extension_um)]
    if not args.make_solver_plots:
        cmd.append("--no-plots")
    return cmd


def read_first_passage(case_dir: Path) -> dict:
    path = case_dir / "anisotropic_calibrated_tip_first_passage_summary.json"
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text())
    except Exception:
        return {}


def run_case(py: str, args, class_name: str, manifest: Path, T_K: int, root: Path) -> dict:
    case_dir = root / class_name / f"T{int(T_K)}_th{args.crystal_theta_deg:g}"
    case_dir.mkdir(parents=True, exist_ok=True)
    complete, extension_um = completion_status(case_dir, args.target_extension_um)
    if complete and args.skip_existing:
        rc = process_case_r_curve(case_dir, class_name, int(T_K), args.target_extension_um)
        fp = read_first_passage(case_dir)
        row = {
            "class": class_name,
            "T_K": int(T_K),
            "status": "skipped_complete",
            "final_extension_um": extension_um,
            "target_extension_um": args.target_extension_um,
            "case_dir": str(case_dir),
            "K_init_MPa_sqrt_m": fp.get("KJ_reference_first_MPa_sqrt_m"),
            **propagation_metrics(rc, fp.get("KJ_reference_first_MPa_sqrt_m")),
        }
        print(f"SKIP {class_name} T={T_K} K: extension={extension_um:.1f} um, events={len(rc)}")
        return row

    cmd = build_command(py, args, class_name, manifest, T_K, case_dir)
    (case_dir / "command.txt").write_text(shlex.join(cmd) + "\n")
    run_audit = {
        "driver": "run_mpz_v9_11_mode_i_rcurve_3T",
        "class": class_name,
        "T_K": int(T_K),
        "target_extension_um": float(args.target_extension_um),
        "branching_enabled": False,
        "max_fronts": 1,
        "material_manifest": str(manifest.resolve()),
        "command": cmd,
    }
    (case_dir / "rcurve_run_audit.json").write_text(json.dumps(run_audit, indent=2))

    print(f"START {class_name} T={T_K} K -> {case_dir}")
    with (case_dir / "run.log").open("w") as log:
        cp = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, text=True)

    complete, extension_um = completion_status(case_dir, args.target_extension_um)
    rc = process_case_r_curve(case_dir, class_name, int(T_K), args.target_extension_um)
    fp = read_first_passage(case_dir)
    if cp.returncode == 0 and complete:
        status = "complete"
        (case_dir / ".long_growth_complete").touch()
    elif cp.returncode == 0:
        status = "incomplete"
    else:
        status = "failed"
        tail = "\n".join((case_dir / "run.log").read_text(errors="replace").splitlines()[-30:])
        print(tail, file=sys.stderr)

    Kinit = fp.get("KJ_reference_first_MPa_sqrt_m")
    row = {
        "class": class_name,
        "T_K": int(T_K),
        "status": status,
        "returncode": int(cp.returncode),
        "final_extension_um": extension_um,
        "target_extension_um": float(args.target_extension_um),
        "target_completed": bool(complete),
        "case_dir": str(case_dir),
        "K_init_MPa_sqrt_m": Kinit,
        "B_final": fp.get("B_final"),
        "B_first_fire_residual": fp.get("B_first_fire_residual"),
        "active_mpz_length_um": fp.get("active_mpz_length_um"),
        **propagation_metrics(rc, Kinit),
    }
    print(
        f"{status.upper():10s} {class_name} T={T_K} K "
        f"extension={extension_um} um, events={len(rc)}, "
        f"Kinit={Kinit}"
    )
    return row


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--parameter-root", type=Path, default=Path("mpz_v9_11_parameters"))
    p.add_argument("--material-class", default="DBTT")
    p.add_argument("--temperatures", default="300 700 1100")
    p.add_argument("--outroot", type=Path, default=Path("runs/mpz_v9_11_DBTT_Rcurve_3T_v1"))
    p.add_argument("--target-extension-um", type=float, default=500.0)
    p.add_argument("--steps", type=int, default=6000)
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
    p.add_argument("--mpz-profile-sector-half-angle-deg", type=float, default=45.0)
    p.add_argument("--mpz-profile-damage-cutoff", type=float, default=0.85)
    p.add_argument("--rJ-cluster-um", type=float, default=20.0)
    p.add_argument("--rJ-outer-um", type=float, default=25.0)
    p.add_argument("--crystal-theta-deg", type=float, default=45.0)
    p.add_argument("--save-snapshots", type=int, default=12)
    p.add_argument("--snapshot-cols", type=int, default=4)
    p.add_argument("--snapshot-by-extension-um", type=float, default=50.0)
    p.add_argument("--make-solver-plots", action=argparse.BooleanOptionalAction, default=False)
    p.add_argument("--skip-existing", action=argparse.BooleanOptionalAction, default=True)
    args = p.parse_args()

    class_name = normalize_class_name(args.material_class)
    temps = [int(round(x)) for x in values(args.temperatures, float)]
    if len(temps) != 3:
        raise SystemExit("--temperatures must contain exactly three temperatures")
    if args.target_extension_um <= 0.0:
        raise SystemExit("--target-extension-um must be positive")

    parameter_root = args.parameter_root.resolve()
    manifest = manifest_for(parameter_root, class_name)
    root = args.outroot.resolve()
    root.mkdir(parents=True, exist_ok=True)

    rows = [run_case(sys.executable, args, class_name, manifest, T, root) for T in temps]
    frame = pd.DataFrame(rows).sort_values("T_K")
    frame.to_csv(root / "rcurve_3T_summary.csv", index=False)
    (root / "rcurve_3T_summary.json").write_text(
        json.dumps(rows, indent=2, default=str)
    )
    make_class_r_curve_overlays(root, [class_name], temps, args.target_extension_um)

    incomplete = frame[frame["status"] != "complete"]
    print(frame.to_string(index=False))
    print("wrote", root / "rcurve_3T_summary.csv")
    if not incomplete.empty:
        raise SystemExit(
            "one or more R-curve cases did not reach the requested extension; "
            "inspect the per-case run.log and rerun with more --steps if needed"
        )


if __name__ == "__main__":
    main()
