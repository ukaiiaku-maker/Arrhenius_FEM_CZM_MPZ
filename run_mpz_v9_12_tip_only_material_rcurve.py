#!/usr/bin/env python3
"""Run the full 2-D v9.12 tip-source-only material transfer gate.

The v9.12 protocol deliberately separates three choices that were previously
conflated:

* cleavage first passage is stochastic through exponential integrated-hazard
  thresholds;
* crack-tip source emission uses the calibrated deterministic expectation by
  default, preserving the mean plastic response used during parameterization;
* fixed-displacement propagation uses the physical raw renewal clock by default.
  Same-load cohesive insertions are clustered and classified as an unstable jump,
  not converted into an artificial smooth R-curve by a prescribed reload gate.

Every solver case enables the existing full FEM field renderer.  The resulting
``field_snapshots_<T>K.png`` contains damage/crack path, dislocation density,
maximum principal FEM stress, and equivalent plastic strain at several accepted
states.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
from typing import Any

import pandas as pd

from arrhenius_fracture.material_rcurve_audit_v912 import audit_campaign
from arrhenius_fracture.mpz_parameterization_v911 import normalize_class_name
from arrhenius_fracture.rcurve_postprocess_v911 import write_cascade_aware_outputs

CLASSES = ("ceramic", "weakT", "DBTT")
CLASS_SEED_OFFSETS = {
    "ceramic": 1_000_003,
    "weakT": 2_000_003,
    "DBTT": 3_000_003,
}


def values(text: str, cast=str):
    return [cast(x) for x in str(text).replace(",", " ").split() if x]


def read_csv_row(path: Path) -> dict[str, Any]:
    if not path.exists() or path.stat().st_size == 0:
        return {}
    try:
        frame = pd.read_csv(path)
    except (pd.errors.EmptyDataError, OSError, ValueError):
        return {}
    return frame.iloc[0].to_dict() if not frame.empty else {}


def effective_seed(base_seed: int, material_class: str, coupling: str) -> int:
    if coupling == "common":
        return int(base_seed)
    return int(base_seed) + CLASS_SEED_OFFSETS[normalize_class_name(material_class)]


def build_command(args, class_name: str, run_root: Path, force_rerun: bool) -> list[str]:
    cmd = [
        sys.executable,
        "run_mpz_v9_11_mode_i_rcurve_3T.py",
        "--parameter-root", str(args.parameter_root),
        "--material-class", class_name,
        "--bulk-plasticity-mode", "tip_only",
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
        "--make-solver-plots",
    ]
    cmd.append("--no-skip-existing" if force_rerun else "--skip-existing")
    return cmd


def run_case(args, base_seed: int, class_name: str, root: Path) -> dict[str, Any]:
    class_name = normalize_class_name(class_name)
    run_root = root / f"seed_{base_seed}" / "tip_only"
    case_dir = run_root / class_name / f"T{int(round(args.T_K))}_th{args.crystal_theta_deg:g}"
    case_dir.mkdir(parents=True, exist_ok=True)
    field_image = case_dir / f"field_snapshots_{int(round(args.T_K))}K.png"
    complete_marker = case_dir / ".long_growth_complete"
    force_rerun = not (
        args.skip_existing and complete_marker.exists() and field_image.exists()
    )
    seed_used = effective_seed(base_seed, class_name, args.rng_coupling)
    cmd = build_command(args, class_name, run_root, force_rerun=force_rerun)
    env = os.environ.copy()
    env.update({
        "ARRHENIUS_EVENT_STATISTICS": "stochastic",
        "ARRHENIUS_STOCHASTIC_EMISSION": "1" if args.stochastic_emission else "0",
        "ARRHENIUS_STOCHASTIC_SEED": str(seed_used),
        "ARRHENIUS_PROPAGATION_CONTROL": args.propagation_control,
    })
    config = {
        "schema": "v9.12_full_field_material_transfer",
        "material_class": class_name,
        "base_seed": int(base_seed),
        "effective_stochastic_seed": int(seed_used),
        "rng_coupling": args.rng_coupling,
        "event_statistics": "stochastic",
        "stochastic_emission": bool(args.stochastic_emission),
        "propagation_control": args.propagation_control,
        "bulk_plasticity_mode": "tip_only",
        "full_field_output_required": True,
        "full_field_rows": [
            "damage_and_crack_path",
            "log10_dislocation_density",
            "maximum_principal_FEM_stress",
            "equivalent_plastic_strain",
        ],
        "command": cmd,
    }
    (case_dir / "v9_12_run_config.json").write_text(json.dumps(config, indent=2))
    log_dir = root / "matrix_logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log = log_dir / f"{class_name}_tip_only_seed{base_seed}_{int(args.T_K)}K.log"

    print(
        f"START {class_name:7s} tip_only seed={base_seed} effective_seed={seed_used} "
        f"propagation={args.propagation_control} stochastic_emission={args.stochastic_emission}"
    )
    if not force_rerun:
        returncode = 0
        reused = True
        print(f"REUSE {class_name:7s}: complete solver result and full-field image exist")
    else:
        with log.open("w") as fp:
            cp = subprocess.run(cmd, env=env, stdout=fp, stderr=subprocess.STDOUT, text=True)
        returncode = int(cp.returncode)
        reused = False

    cascade = write_cascade_aware_outputs(
        case_dir,
        args.T_K,
        relative_load_tolerance=args.cluster_relative_load_tolerance,
        absolute_load_tolerance_m=args.cluster_absolute_load_tolerance_m,
    )
    temp_summary = read_csv_row(run_root / class_name / "rcurve_temperature_summary.csv")
    row = {
        **temp_summary,
        **cascade,
        "class": class_name,
        "base_seed": int(base_seed),
        "effective_stochastic_seed": int(seed_used),
        "rng_coupling": args.rng_coupling,
        "stochastic_emission": bool(args.stochastic_emission),
        "propagation_control": args.propagation_control,
        "subprocess_returncode": returncode,
        "solver_output_reused": reused,
        "field_snapshot_image": str(field_image) if field_image.exists() else None,
        "field_snapshot_image_present": bool(field_image.exists()),
        "case_dir": str(case_dir),
        "log": str(log),
    }
    pd.DataFrame([row]).to_csv(case_dir / "v9_12_case_summary.csv", index=False)
    (case_dir / "v9_12_case_summary.json").write_text(
        json.dumps(row, indent=2, default=str)
    )
    print(
        f"DONE  {class_name:7s} rc={returncode} status={row.get('status')} "
        f"ext={row.get('final_extension_um')} load_events={row.get('n_independent_load_events')} "
        f"fields={field_image.exists()}"
    )
    return row


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--parameter-root", type=Path, default=Path("mpz_v9_11_parameters"))
    p.add_argument("--outroot", type=Path, default=Path("runs/mpz_v9_12_tip_only_material_rcurve_700K_v1"))
    p.add_argument("--seeds", default="1")
    p.add_argument("--classes", default="ceramic weakT DBTT")
    p.add_argument("--T-K", type=float, default=700.0)
    p.add_argument("--target-extension-um", type=float, default=100.0)
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
    p.add_argument("--crystal-theta-deg", type=float, default=45.0)
    p.add_argument("--save-snapshots", type=int, default=5)
    p.add_argument("--snapshot-cols", type=int, default=5)
    p.add_argument("--snapshot-by-extension-um", type=float, default=25.0)
    p.add_argument(
        "--stochastic-emission",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Opt-in discrete finite-site emission. Default preserves the calibrated mean emission response.",
    )
    p.add_argument(
        "--propagation-control",
        choices=("raw", "event_reload"),
        default="raw",
        help="raw is the physical fixed-displacement cascade; event_reload is diagnostic only.",
    )
    p.add_argument(
        "--rng-coupling",
        choices=("independent", "common"),
        default="independent",
        help="independent gives each material a distinct reproducible threshold stream; common is a variance-reduction diagnostic.",
    )
    p.add_argument("--cluster-relative-load-tolerance", type=float, default=1.0e-4)
    p.add_argument("--cluster-absolute-load-tolerance-m", type=float, default=1.0e-12)
    p.add_argument("--skip-existing", action=argparse.BooleanOptionalAction, default=True)
    args = p.parse_args()

    root = args.outroot.resolve()
    root.mkdir(parents=True, exist_ok=True)
    seeds = values(args.seeds, int)
    classes = [normalize_class_name(x) for x in values(args.classes, str)]
    rows: list[dict[str, Any]] = []
    for seed in seeds:
        for cls in classes:
            rows.append(run_case(args, seed, cls, root))
        audit = audit_campaign(root, seed, args.T_K, classes=classes, bulk_mode="tip_only")
        print(
            f"AUDIT seed={seed}: gate={audit['material_rcurve_gate_passed']} "
            f"interpretation={audit['interpretation']} missing_fields={audit['missing_full_field_images']}"
        )

    frame = pd.DataFrame(rows)
    frame.to_csv(root / "v9_12_campaign_summary.csv", index=False)
    (root / "v9_12_campaign_summary.json").write_text(
        json.dumps(rows, indent=2, default=str)
    )
    print(frame.to_string(index=False))
    print("wrote", root / "v9_12_campaign_summary.csv")


if __name__ == "__main__":
    main()
