#!/usr/bin/env python3
"""Run the v10.0.5.9 production-initialization J parity audit.

Two or more elastic production probes are launched at different grip openings.
Each probe traverses the actual audited v10.0.5.5 production initialization path,
records the post-equilibrium state, and evaluates production/full/no-exclusion J
metrics. The campaign then compares J/sigma^2 with the converged v10.0.5.8
fixed-grip reference and verifies quadratic elastic load scaling.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import subprocess
import sys
from typing import Any, Iterable

from arrhenius_fracture.production_j_parity_v10059 import (
    CONTOUR_CSV,
    PROBE_CSV,
    PROBE_JSON,
    SUMMARY_JSON,
    analyze_production_j_parity_v10059,
)

ENTRY_MODULE = (
    "arrhenius_fracture."
    "mode_i_first_passage_v10_0_5_9_production_j_probe"
)


def _values(text: str, scale: float = 1.0) -> list[float]:
    return [
        float(token) * scale
        for token in str(text).replace(",", " ").split()
        if token
    ]


def _token(value: float) -> str:
    return f"{value:g}".replace("-", "m").replace(".", "p")


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-json", type=Path, required=True)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("runs/v10_0_5_9_production_j_parity_v1"),
    )
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--openings-um", default="1 2")
    parser.add_argument("--contour-outer-um", default="180 240 300")
    parser.add_argument("--selected-outer-um", type=float, default=240.0)
    parser.add_argument("--parity-rel-tol", type=float, default=0.10)
    parser.add_argument("--elastic-scaling-rel-tol", type=float, default=0.02)
    parser.add_argument("--energy-closure-rel-tol", type=float, default=1.0e-6)
    parser.add_argument("--temperature-K", type=float, default=700.0)
    parser.add_argument("--material-class", default="ceramic")
    parser.add_argument("--nx", type=int, default=36)
    parser.add_argument("--ny", type=int, default=72)
    parser.add_argument("--tip-h-um", type=float, default=2.5)
    parser.add_argument("--tip-ratio", type=float, default=1.15)
    parser.add_argument("--cluster-J-outer-um", type=float, default=240.0)
    parser.add_argument("--local-J-outer-um", type=float, default=100.0)
    parser.add_argument("--L-pz-um", type=float, default=20.0)
    parser.add_argument("--mpz-length-um", type=float, default=100.0)
    parser.add_argument("--mpz-n-bins", type=int, default=80)
    parser.add_argument("--da-phys-um", type=float, default=5.0)
    parser.add_argument(
        "--crack-backend",
        choices=["sharp_wake", "edge_split_czm", "adaptive_czm"],
        default="adaptive_czm",
    )
    parser.add_argument(
        "--crystal-aniso",
        action=argparse.BooleanOptionalAction,
        default=True,
    )
    parser.add_argument("--crystal-theta-deg", type=float, default=45.0)
    return parser


def _probe_command(args, case_dir: Path, opening_m: float) -> list[str]:
    """Build one elastic production-probe command.

    v10.0.5.4 and every audited VHCF descendant require
    ``--cycle-block-mode hazard_limited``. The probe still consumes exactly one
    cycle because ``cycles-max``, ``block-cycles`` and ``max-block-cycles`` are
    all one; hazard-limited is therefore a wrapper-contract requirement, not a
    change to the elastic mechanics state being recorded.
    """
    command = [
        args.python,
        "-m",
        ENTRY_MODULE,
        "--v10-material-class",
        str(args.material_class),
        "--czm-opening-coupling",
        "clock_linear",
        "--mode",
        "2d",
        "--temperatures",
        f"{args.temperature_K:.17g}",
        "--out",
        str(case_dir),
        "--steps",
        "1",
        "--dU",
        f"{opening_m:.17g}",
        "--dt",
        "1e-12",
        "--n-stagger",
        "2",
        "--fatigue-cycles",
        "--fatigue-hold-load",
        "--R",
        "0.1",
        "--frequency-Hz",
        "1000",
        "--cycles-max",
        "1",
        "--block-cycles",
        "1",
        "--max-block-cycles",
        "1",
        "--min-block-cycles",
        "1e-9",
        "--cycle-block-mode",
        "hazard_limited",
        "--target-dB",
        "0.25",
        "--target-dN-store",
        "0.25",
        "--n-phase",
        "8",
        "--no-cyclic-mechanics",
        "--nx",
        str(args.nx),
        "--ny",
        str(args.ny),
        "--tip-h-fine",
        f"{args.tip_h_um * 1.0e-6:.17g}",
        "--tip-ratio",
        f"{args.tip_ratio:.17g}",
        "--da-phys",
        f"{args.da_phys_um * 1.0e-6:.17g}",
        "--rJ-outer",
        f"{args.local_J_outer_um * 1.0e-6:.17g}",
        "--rJ-cluster",
        f"{args.cluster_J_outer_um * 1.0e-6 / 8.0:.17g}",
        "--L-pz",
        f"{args.L_pz_um * 1.0e-6:.17g}",
        "--mpz-length-um",
        f"{args.mpz_length_um:.17g}",
        "--mpz-n-bins",
        str(args.mpz_n_bins),
        "--crack-backend",
        str(args.crack_backend),
        "--max-fronts",
        "1",
        "--target-crack-extension-um",
        "100",
        "--save-snapshots",
        "0",
        "--print-every",
        "1",
        "--no-plots",
    ]
    if args.crystal_aniso:
        command.extend(
            ["--crystal-aniso", "--crystal-theta-deg", f"{args.crystal_theta_deg:.17g}"]
        )
    return command


def main(argv: Iterable[str] | None = None) -> int:
    args = build_parser().parse_args(list(argv) if argv is not None else None)
    reference_path = args.reference_json.resolve()
    if not reference_path.exists():
        raise SystemExit(f"missing v10.0.5.8 reference JSON: {reference_path}")
    reference = json.loads(reference_path.read_text())
    if reference.get("schema") != "fixed_grip_elastic_convergence_v10_0_5_8":
        raise SystemExit(
            "--reference-json must be a v10.0.5.8 fixed-grip convergence summary"
        )
    if not bool(reference.get("passed", False)):
        raise SystemExit("the supplied v10.0.5.8 reference did not pass its mechanics gate")

    openings = sorted(set(_values(args.openings_um, 1.0e-6)))
    contours = sorted(set(_values(args.contour_outer_um, 1.0e-6)))
    if len(openings) < 2 or any(value <= 0.0 for value in openings):
        raise SystemExit("at least two positive --openings-um values are required")
    if len(contours) < 3 or any(value <= 0.0 for value in contours):
        raise SystemExit("at least three positive --contour-outer-um values are required")

    root = args.out.resolve()
    root.mkdir(parents=True, exist_ok=True)
    probes: list[dict[str, Any]] = []
    process_rows: list[dict[str, Any]] = []
    env = os.environ.copy()
    env["ARRHENIUS_V10059_CONTOURS_UM"] = " ".join(
        f"{value * 1.0e6:.17g}" for value in contours
    )
    env["ARRHENIUS_EVENT_STATISTICS"] = "mean_field"
    env["ARRHENIUS_STOCHASTIC_EMISSION"] = "0"
    env["ARRHENIUS_VHCF_FEM_CACHE"] = "0"

    for opening in openings:
        case_dir = root / f"opening_{_token(opening * 1.0e6)}um"
        case_dir.mkdir(parents=True, exist_ok=True)
        command = _probe_command(args, case_dir, opening)
        log_path = case_dir / "production_j_probe.log"
        with log_path.open("w") as log:
            process = subprocess.run(
                command,
                env=env,
                stdout=log,
                stderr=subprocess.STDOUT,
                text=True,
                check=False,
            )
        probe_path = case_dir / PROBE_JSON
        process_rows.append(
            {
                "requested_opening_um": opening * 1.0e6,
                "subprocess_returncode": int(process.returncode),
                "probe_exists": probe_path.exists(),
                "case_directory": str(case_dir),
                "log": str(log_path),
                "command": " ".join(command),
            }
        )
        if not probe_path.exists():
            raise RuntimeError(
                f"production probe did not write {probe_path}; see {log_path}"
            )
        probe = json.loads(probe_path.read_text())
        probe["subprocess_returncode"] = int(process.returncode)
        probe["case_directory"] = str(case_dir)
        probes.append(probe)

    summary = analyze_production_j_parity_v10059(
        reference=reference,
        probes=probes,
        selected_outer_radius_m=args.selected_outer_um * 1.0e-6,
        parity_relative_tolerance=args.parity_rel_tol,
        elastic_scaling_relative_tolerance=args.elastic_scaling_rel_tol,
        energy_closure_relative_tolerance=args.energy_closure_rel_tol,
    )
    summary.update(
        {
            "reference_json": str(reference_path),
            "probe_processes": process_rows,
            "production_configuration": {
                "material_class": args.material_class,
                "temperature_K": args.temperature_K,
                "nx": args.nx,
                "ny": args.ny,
                "tip_h_um": args.tip_h_um,
                "tip_ratio": args.tip_ratio,
                "cluster_J_outer_um": args.cluster_J_outer_um,
                "local_J_outer_um": args.local_J_outer_um,
                "L_pz_um": args.L_pz_um,
                "mpz_length_um": args.mpz_length_um,
                "mpz_n_bins": args.mpz_n_bins,
                "da_phys_um": args.da_phys_um,
                "crack_backend": args.crack_backend,
                "crystal_anisotropic": args.crystal_aniso,
                "crystal_theta_deg": args.crystal_theta_deg,
            },
            "subprocess_nonzero_count": sum(
                int(row["subprocess_returncode"] != 0) for row in process_rows
            ),
            "mechanics_probe_complete": all(row["probe_exists"] for row in process_rows),
        }
    )

    all_contours: list[dict[str, Any]] = []
    for probe in probes:
        for row in probe.get("contours", []):
            combined = {
                "Uapp_um": float(probe["Uapp_um"]),
                "sigma_gross_MPa": float(probe["sigma_gross_MPa"]),
                "case_directory": probe["case_directory"],
                **dict(row),
            }
            all_contours.append(combined)
    _write_csv(root / PROBE_CSV, list(summary.get("cases", [])))
    _write_csv(root / CONTOUR_CSV, all_contours)
    (root / SUMMARY_JSON).write_text(json.dumps(summary, indent=2, default=str))

    print(f"PRODUCTION J PARITY STATUS: {summary['status']}")
    print(root / SUMMARY_JSON)
    return 0 if bool(summary.get("passed", False)) else 2


if __name__ == "__main__":
    raise SystemExit(main())
