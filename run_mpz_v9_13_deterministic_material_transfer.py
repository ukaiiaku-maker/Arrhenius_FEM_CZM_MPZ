#!/usr/bin/env python3
"""Run the full 2-D v9.13 material-transfer gate.

The default gate is deterministic in integrated-hazard action and uses the
calibrated expected finite-site emission response.  This isolates constitutive
material differences from single-realization Poisson noise.  Stochastic mode is
available for subsequent ensemble work but cannot pass the deterministic gate.
"""
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

from arrhenius_fracture.material_rcurve_audit_v913 import audit_campaign
from arrhenius_fracture.mpz_parameterization_v911 import normalize_class_name
from arrhenius_fracture.rcurve_postprocess_v911 import write_cascade_aware_outputs

CLASSES = ("ceramic", "weakT", "DBTT")
CLASS_SEED_OFFSETS = {"ceramic": 1_000_003, "weakT": 2_000_003, "DBTT": 3_000_003}


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
        frame = frame[frame["class"].astype(str).map(normalize_class_name) == normalize_class_name(material_class)]
    if T_K is not None:
        for name in ("T_K", "T"):
            if name in frame.columns:
                q = pd.to_numeric(frame[name], errors="coerce")
                frame = frame[(q - float(T_K)).abs() < 0.5]
                break
    return frame.iloc[-1].to_dict() if not frame.empty else {}


def effective_seed(base_seed: int, material_class: str, coupling: str) -> int:
    if coupling == "common":
        return int(base_seed)
    return int(base_seed) + CLASS_SEED_OFFSETS[normalize_class_name(material_class)]


def build_command(args, class_name: str, run_root: Path, force_rerun: bool) -> list[str]:
    cmd = [
        sys.executable, "run_mpz_v9_13_mode_i_rcurve.py",
        "--parameter-root", str(args.parameter_root),
        "--material-class", class_name,
        "--bulk-plasticity-mode", "tip_only",
        "--temperatures", str(int(round(args.T_K))),
        "--outroot", str(run_root),
        "--target-extension-um", str(args.target_extension_um),
        "--steps", str(args.steps),
        "--nx", str(args.nx), "--ny", str(args.ny),
        "--tip-h-fine", str(args.tip_h_fine), "--tip-ratio", str(args.tip_ratio),
        "--dU", str(args.dU), "--dt", str(args.dt),
        "--n-stagger", str(args.n_stagger), "--print-every", str(args.print_every),
        "--adaptive-event-target", str(args.adaptive_event_target),
        "--da-phys-um", str(args.da_phys_um),
        "--mpz-length-um", str(args.mpz_length_um), "--mpz-n-bins", str(args.mpz_n_bins),
        "--crystal-theta-deg", str(args.crystal_theta_deg),
        "--save-snapshots", str(args.save_snapshots),
        "--snapshot-cols", str(args.snapshot_cols),
        "--snapshot-by-extension-um", str(args.snapshot_by_extension_um),
        "--make-solver-plots",
    ]
    cmd.append("--no-skip-existing" if force_rerun else "--skip-existing")
    return cmd


def _copy_temperature_summary(run_root: Path, case_dir: Path, class_name: str, T_K: float) -> tuple[dict[str, Any], str | None]:
    # The v9.11 base runner writes this at its --outroot, not under class_name.
    source = run_root / "rcurve_temperature_summary.csv"
    row = read_csv_row(source, class_name, T_K)
    if not source.exists():
        return row, None
    destination = case_dir / "rcurve_temperature_summary_v913.csv"
    shutil.copy2(source, destination)
    return row, str(destination)


def run_case(args, base_seed: int, class_name: str, root: Path) -> dict[str, Any]:
    class_name = normalize_class_name(class_name)
    run_root = root / f"seed_{base_seed}" / "tip_only"
    case_dir = run_root / class_name / f"T{int(round(args.T_K))}_th{args.crystal_theta_deg:g}"
    case_dir.mkdir(parents=True, exist_ok=True)
    field_image = case_dir / f"field_snapshots_{int(round(args.T_K))}K.png"
    zoom_image = case_dir / f"field_snapshots_tip_zoom_{int(round(args.T_K))}K.png"
    manifest = case_dir / f"field_snapshot_manifest_{int(round(args.T_K))}K.json"
    marker = case_dir / ".long_growth_complete"
    force_rerun = not (
        args.skip_existing and marker.exists() and field_image.exists()
        and zoom_image.exists() and manifest.exists()
    )
    seed_used = effective_seed(base_seed, class_name, args.rng_coupling)
    cmd = build_command(args, class_name, run_root, force_rerun)
    env = os.environ.copy()
    env.update({
        "ARRHENIUS_EVENT_STATISTICS": args.event_statistics,
        "ARRHENIUS_STOCHASTIC_EMISSION": "1" if args.stochastic_emission else "0",
        "ARRHENIUS_STOCHASTIC_SEED": str(seed_used),
        "ARRHENIUS_PROPAGATION_CONTROL": args.propagation_control,
    })
    config = {
        "schema": "v9.13_deterministic_material_transfer",
        "material_class": class_name, "base_seed": int(base_seed),
        "effective_stochastic_seed": int(seed_used), "rng_coupling": args.rng_coupling,
        "event_statistics": args.event_statistics,
        "stochastic_emission": bool(args.stochastic_emission),
        "propagation_control": args.propagation_control,
        "bulk_plasticity_mode": "tip_only", "target_extension_um": float(args.target_extension_um),
        "protocol_role": (
            "deterministic_mean_material_transfer_gate"
            if args.event_statistics == "deterministic" and not args.stochastic_emission
            else "stochastic_ensemble_realization"
        ),
        "full_field_output_required": True,
        "tip_zoom_output_required": True,
        "command": cmd,
    }
    (case_dir / "v9_13_run_config.json").write_text(json.dumps(config, indent=2))
    log_dir = root / "matrix_logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log = log_dir / f"{class_name}_tip_only_seed{base_seed}_{int(args.T_K)}K.log"

    print(
        f"START {class_name:7s} seed={base_seed} event_statistics={args.event_statistics} "
        f"stochastic_emission={args.stochastic_emission} propagation={args.propagation_control}"
    )
    if not force_rerun:
        returncode, reused = 0, True
        print(f"REUSE {class_name:7s}: completion marker and required field outputs exist")
    else:
        with log.open("w") as fp:
            cp = subprocess.run(cmd, env=env, stdout=fp, stderr=subprocess.STDOUT, text=True)
        returncode, reused = int(cp.returncode), False

    cascade = write_cascade_aware_outputs(
        case_dir, args.T_K,
        relative_load_tolerance=args.cluster_relative_load_tolerance,
        absolute_load_tolerance_m=args.cluster_absolute_load_tolerance_m,
    )
    temp_summary, copied_summary = _copy_temperature_summary(run_root, case_dir, class_name, args.T_K)
    row = {
        **temp_summary, **cascade,
        "class": class_name, "base_seed": int(base_seed),
        "effective_stochastic_seed": int(seed_used), "rng_coupling": args.rng_coupling,
        "event_statistics": args.event_statistics,
        "stochastic_emission": bool(args.stochastic_emission),
        "propagation_control": args.propagation_control,
        "subprocess_returncode": returncode, "solver_output_reused": reused,
        "completion_marker_present": marker.exists(),
        "field_snapshot_image": str(field_image) if field_image.exists() else None,
        "field_snapshot_image_present": field_image.exists(),
        "field_snapshot_tip_zoom_image": str(zoom_image) if zoom_image.exists() else None,
        "field_snapshot_tip_zoom_present": zoom_image.exists(),
        "field_snapshot_manifest_present": manifest.exists(),
        "copied_temperature_summary": copied_summary,
        "case_dir": str(case_dir), "log": str(log),
    }
    pd.DataFrame([row]).to_csv(case_dir / "v9_13_case_summary.csv", index=False)
    (case_dir / "v9_13_case_summary.json").write_text(json.dumps(row, indent=2, default=str))
    print(
        f"DONE  {class_name:7s} rc={returncode} status={row.get('status')} "
        f"ext={row.get('final_extension_um')} fields={field_image.exists()} zoom={zoom_image.exists()}"
    )
    return row


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--parameter-root", type=Path, default=Path("mpz_v9_11_parameters"))
    p.add_argument("--outroot", type=Path, default=Path("runs/mpz_v9_13_deterministic_material_transfer_700K_v1"))
    p.add_argument("--seeds", default="1")
    p.add_argument("--classes", default="ceramic weakT DBTT")
    p.add_argument("--T-K", type=float, default=700.0)
    p.add_argument("--target-extension-um", type=float, default=100.0)
    p.add_argument("--steps", type=int, default=6000)
    p.add_argument("--nx", type=int, default=36); p.add_argument("--ny", type=int, default=72)
    p.add_argument("--tip-h-fine", type=float, default=1.0e-6); p.add_argument("--tip-ratio", type=float, default=1.20)
    p.add_argument("--dU", type=float, default=2.0e-7); p.add_argument("--dt", type=float, default=8.4)
    p.add_argument("--n-stagger", type=int, default=2); p.add_argument("--print-every", type=int, default=25)
    p.add_argument("--adaptive-event-target", type=float, default=0.15)
    p.add_argument("--da-phys-um", type=float, default=5.0)
    p.add_argument("--mpz-length-um", type=float, default=100.0); p.add_argument("--mpz-n-bins", type=int, default=200)
    p.add_argument("--crystal-theta-deg", type=float, default=45.0)
    p.add_argument("--save-snapshots", type=int, default=5); p.add_argument("--snapshot-cols", type=int, default=5)
    p.add_argument("--snapshot-by-extension-um", type=float, default=25.0)
    p.add_argument("--event-statistics", choices=("deterministic", "stochastic"), default="deterministic")
    p.add_argument("--stochastic-emission", action=argparse.BooleanOptionalAction, default=False)
    p.add_argument("--propagation-control", choices=("raw", "event_reload"), default="raw")
    p.add_argument("--rng-coupling", choices=("common", "independent"), default="common")
    p.add_argument("--cluster-relative-load-tolerance", type=float, default=1.0e-4)
    p.add_argument("--cluster-absolute-load-tolerance-m", type=float, default=1.0e-12)
    p.add_argument("--skip-existing", action=argparse.BooleanOptionalAction, default=True)
    args = p.parse_args()

    root = args.outroot.resolve(); root.mkdir(parents=True, exist_ok=True)
    seeds = values(args.seeds, int)
    classes = [normalize_class_name(x) for x in values(args.classes, str)]
    rows: list[dict[str, Any]] = []
    for seed in seeds:
        for cls in classes:
            rows.append(run_case(args, seed, cls, root))
        audit = audit_campaign(root, seed, args.T_K, classes=classes, bulk_mode="tip_only")
        print(
            f"AUDIT seed={seed}: transfer_gate={audit['material_transfer_gate_passed']} "
            f"interpretation={audit['interpretation']} failed={audit['failed_solver_cases']} "
            f"incomplete={audit['incomplete_cases']} similar={audit['strongly_similar_pairs']}"
        )
    frame = pd.DataFrame(rows)
    frame.to_csv(root / "v9_13_campaign_summary.csv", index=False)
    (root / "v9_13_campaign_summary.json").write_text(json.dumps(rows, indent=2, default=str))
    print(frame.to_string(index=False))
    print("wrote", root / "v9_13_campaign_summary.csv")


if __name__ == "__main__":
    main()
