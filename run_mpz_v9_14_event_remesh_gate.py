#!/usr/bin/env python3
"""Run the v9.14 event-localized, conservatively remeshed 700 K gate."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Any

import pandas as pd

from arrhenius_fracture.event_remesh_audit_v914 import audit_campaign
from arrhenius_fracture.mpz_parameterization_v911 import normalize_class_name
from arrhenius_fracture.rcurve_postprocess_v911 import write_cascade_aware_outputs

CLASSES = ("ceramic", "weakT", "DBTT")
ADAPTIVE_EVENT_COORDINATE = "absolute_integrated_hazard_action"


def values(text: str, cast=str):
    return [cast(x) for x in str(text).replace(",", " ").split() if x]


def read_csv_row(path: Path, material_class: str | None = None, T_K: float | None = None) -> dict[str, Any]:
    if not path.exists() or path.stat().st_size == 0:
        return {}
    try:
        frame = pd.read_csv(path)
    except (pd.errors.EmptyDataError, OSError, ValueError):
        return {}
    if material_class is not None and "class" in frame.columns:
        frame = frame[
            frame["class"].astype(str).map(normalize_class_name)
            == normalize_class_name(material_class)
        ]
    if T_K is not None:
        for name in ("T_K", "T"):
            if name in frame.columns:
                q = pd.to_numeric(frame[name], errors="coerce")
                frame = frame[(q - float(T_K)).abs() < 0.5]
                break
    return frame.iloc[-1].to_dict() if not frame.empty else {}


def build_command(args, class_name: str, run_root: Path, force_rerun: bool) -> list[str]:
    cmd = [
        sys.executable,
        "run_mpz_v9_14_mode_i_rcurve.py",
        "--parameter-root", str(args.parameter_root),
        "--material-class", class_name,
        "--bulk-plasticity-mode", "tip_only",
        "--temperatures", str(int(round(args.T_K))),
        "--outroot", str(run_root),
        "--target-extension-um", str(args.target_extension_um),
        "--steps", str(args.steps),
        "--nx", str(args.nx), "--ny", str(args.ny),
        "--tip-h-fine", str(args.tip_h_fine),
        "--tip-ratio", str(args.tip_ratio),
        "--dU", str(args.dU), "--dt", str(args.dt),
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
        "--make-solver-plots",
    ]
    cmd.append("--no-skip-existing" if force_rerun else "--skip-existing")
    return cmd


def _copy_temperature_summary(
    run_root: Path, case_dir: Path, class_name: str, T_K: float
) -> tuple[dict[str, Any], str | None]:
    source = run_root / "rcurve_temperature_summary.csv"
    row = read_csv_row(source, class_name, T_K)
    if not source.exists():
        return row, None
    destination = case_dir / "rcurve_temperature_summary_v914.csv"
    shutil.copy2(source, destination)
    return row, str(destination)


def run_case(args, seed: int, class_name: str, root: Path) -> dict[str, Any]:
    class_name = normalize_class_name(class_name)
    run_root = root / f"seed_{seed}" / "tip_only"
    case_dir = run_root / class_name / f"T{int(round(args.T_K))}_th{args.crystal_theta_deg:g}"
    case_dir.mkdir(parents=True, exist_ok=True)
    marker = case_dir / ".long_growth_complete"
    required = [
        case_dir / f"field_snapshots_{int(round(args.T_K))}K.png",
        case_dir / f"field_snapshots_tip_zoom_{int(round(args.T_K))}K.png",
        case_dir / f"field_snapshot_manifest_{int(round(args.T_K))}K.json",
        case_dir / f"czm_{int(round(args.T_K)):04d}K" / "event_remesh_audit_v914.json",
        case_dir / "event_equilibrium_audit_v914.json",
    ]
    force_rerun = not (
        args.skip_existing and marker.exists() and all(path.exists() for path in required)
    )
    cmd = build_command(args, class_name, run_root, force_rerun)
    env = os.environ.copy()
    env.update({
        "ARRHENIUS_EVENT_STATISTICS": args.event_statistics,
        "ARRHENIUS_STOCHASTIC_EMISSION": "1" if args.stochastic_emission else "0",
        "ARRHENIUS_STOCHASTIC_SEED": str(seed),
        "ARRHENIUS_PROPAGATION_CONTROL": "raw",
        "ARRHENIUS_EVENT_REMESH_V914": "1",
        "ARRHENIUS_EVENT_REMESH_TARGET_H_M": str(args.event_remesh_target_h_m),
        "ARRHENIUS_EVENT_REMESH_PATCH_RADIUS_M": str(args.event_remesh_patch_radius_um * 1.0e-6),
        "ARRHENIUS_EVENT_REMESH_MAX_EDGE_SPLITS": str(args.event_remesh_max_edge_splits),
        "ARRHENIUS_EVENT_REMESH_TARGET_EDGE_FACTOR": str(args.event_remesh_target_edge_factor),
        "ARRHENIUS_EVENT_REMESH_MIN_QUALITY": str(args.event_remesh_min_quality),
        "ARRHENIUS_EVENT_REMESH_REQUIRE_EQUILIBRIUM": "1",
    })
    config = {
        "schema": "v9.14_event_remesh_material_transfer",
        "material_class": class_name,
        "seed": int(seed),
        "event_statistics": args.event_statistics,
        "stochastic_emission": bool(args.stochastic_emission),
        "propagation_control": "raw",
        "bulk_plasticity_mode": "tip_only",
        "effective_crack_backend": "event_remesh_czm",
        "adaptive_event_coordinate": ADAPTIVE_EVENT_COORDINATE,
        "adaptive_event_action_tolerance": float(args.adaptive_event_target),
        "one_physical_event_per_renewal": True,
        "same_time_same_load_post_event_equilibrium": True,
        "event_remesh_target_h_m": float(args.event_remesh_target_h_m),
        "event_remesh_patch_radius_um": float(args.event_remesh_patch_radius_um),
        "event_remesh_max_edge_splits": int(args.event_remesh_max_edge_splits),
        "target_extension_um": float(args.target_extension_um),
        "command": cmd,
    }
    (case_dir / "v9_14_run_config.json").write_text(json.dumps(config, indent=2))
    log_dir = root / "matrix_logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log = log_dir / f"{class_name}_tip_only_seed{seed}_{int(args.T_K)}K.log"

    print(
        f"START {class_name:7s} seed={seed} event-remesh "
        f"target_h={args.event_remesh_target_h_m:g} m "
        f"absolute_dB_tol={args.adaptive_event_target:g}"
    )
    if not force_rerun:
        returncode, reused = 0, True
        print(f"REUSE {class_name:7s}: completion marker and v9.14 audits exist")
    else:
        with log.open("w") as fp:
            cp = subprocess.run(cmd, env=env, stdout=fp, stderr=subprocess.STDOUT, text=True)
        returncode, reused = int(cp.returncode), False

    cascade = write_cascade_aware_outputs(
        case_dir,
        args.T_K,
        relative_load_tolerance=args.cluster_relative_load_tolerance,
        absolute_load_tolerance_m=args.cluster_absolute_load_tolerance_m,
    )
    temp_summary, copied_summary = _copy_temperature_summary(
        run_root, case_dir, class_name, args.T_K
    )
    row = {
        **temp_summary,
        **cascade,
        "class": class_name,
        "seed": int(seed),
        "T_K": float(args.T_K),
        "event_statistics": args.event_statistics,
        "stochastic_emission": bool(args.stochastic_emission),
        "propagation_control": "raw",
        "effective_crack_backend": "event_remesh_czm",
        "adaptive_event_coordinate": ADAPTIVE_EVENT_COORDINATE,
        "adaptive_event_action_tolerance": float(args.adaptive_event_target),
        "subprocess_returncode": returncode,
        "solver_output_reused": reused,
        "completion_marker_present": marker.exists(),
        "event_remesh_audit_present": required[3].exists(),
        "event_equilibrium_audit_present": required[4].exists(),
        "copied_temperature_summary": copied_summary,
        "case_dir": str(case_dir),
        "log": str(log),
    }
    pd.DataFrame([row]).to_csv(case_dir / "v9_14_case_summary.csv", index=False)
    (case_dir / "v9_14_case_summary.json").write_text(
        json.dumps(row, indent=2, default=str)
    )
    print(
        f"DONE  {class_name:7s} rc={returncode} status={row.get('status')} "
        f"ext={row.get('final_extension_um')} remesh={required[3].exists()} "
        f"equilibrium={required[4].exists()}"
    )
    return row


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--parameter-root", type=Path, default=Path("mpz_v9_11_parameters"))
    p.add_argument("--outroot", type=Path, default=Path("runs/mpz_v9_14_event_remesh_700K_v1"))
    p.add_argument("--seeds", default="1")
    p.add_argument("--classes", default="ceramic weakT DBTT")
    p.add_argument("--T-K", type=float, default=700.0)
    p.add_argument("--target-extension-um", type=float, default=50.0)
    p.add_argument("--steps", type=int, default=12000)
    p.add_argument("--nx", type=int, default=36)
    p.add_argument("--ny", type=int, default=72)
    p.add_argument("--tip-h-fine", type=float, default=1.0e-6)
    p.add_argument("--tip-ratio", type=float, default=1.20)
    p.add_argument("--dU", type=float, default=2.0e-7)
    p.add_argument("--dt", type=float, default=8.4)
    p.add_argument("--n-stagger", type=int, default=2)
    p.add_argument("--print-every", type=int, default=25)
    p.add_argument(
        "--adaptive-event-target",
        type=float,
        default=0.01,
        help="maximum accepted absolute integrated-hazard increment dB",
    )
    p.add_argument("--da-phys-um", type=float, default=5.0)
    p.add_argument("--mpz-length-um", type=float, default=100.0)
    p.add_argument("--mpz-n-bins", type=int, default=200)
    p.add_argument("--crystal-theta-deg", type=float, default=45.0)
    p.add_argument("--save-snapshots", type=int, default=5)
    p.add_argument("--snapshot-cols", type=int, default=5)
    p.add_argument("--snapshot-by-extension-um", type=float, default=10.0)
    p.add_argument("--event-statistics", choices=("deterministic", "stochastic"), default="deterministic")
    p.add_argument("--stochastic-emission", action=argparse.BooleanOptionalAction, default=False)
    p.add_argument("--event-remesh-target-h-m", type=float, default=1.0e-6)
    p.add_argument("--event-remesh-patch-radius-um", type=float, default=25.0)
    p.add_argument("--event-remesh-max-edge-splits", type=int, default=256)
    p.add_argument("--event-remesh-target-edge-factor", type=float, default=1.25)
    p.add_argument("--event-remesh-min-quality", type=float, default=0.02)
    p.add_argument("--cluster-relative-load-tolerance", type=float, default=1.0e-4)
    p.add_argument("--cluster-absolute-load-tolerance-m", type=float, default=1.0e-12)
    p.add_argument("--skip-existing", action=argparse.BooleanOptionalAction, default=True)
    args = p.parse_args()

    root = args.outroot.resolve()
    root.mkdir(parents=True, exist_ok=True)
    classes = [normalize_class_name(x) for x in values(args.classes, str)]
    rows: list[dict[str, Any]] = []
    for seed in values(args.seeds, int):
        for cls in classes:
            rows.append(run_case(args, seed, cls, root))
            pd.DataFrame(rows).to_csv(root / "v9_14_campaign_summary.partial.csv", index=False)
        audit = audit_campaign(root, seed, args.T_K, classes=classes, bulk_mode="tip_only")
        print(
            f"AUDIT seed={seed}: numerical={audit['numerical_event_remesh_gate_passed']} "
            f"material={audit['material_transfer_gate_passed_v914']} "
            f"interpretation={audit['interpretation']} "
            f"failed_numerical={audit['failed_numerical_remesh_cases']}"
        )
    frame = pd.DataFrame(rows)
    frame.to_csv(root / "v9_14_campaign_summary.csv", index=False)
    (root / "v9_14_campaign_summary.json").write_text(
        json.dumps(rows, indent=2, default=str)
    )
    print(frame.to_string(index=False))
    print("wrote", root / "v9_14_campaign_summary.csv")


if __name__ == "__main__":
    main()
