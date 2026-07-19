#!/usr/bin/env python3
"""Run production J parity with a fixed physical refinement radius."""
from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
from typing import Any, Iterable

from arrhenius_fracture.production_j_refinement_support_v100510 import (
    CONTOUR_CSV,
    PROBE_CSV,
    PROBE_JSON,
    RADIAL_CSV,
    SUMMARY_JSON,
    analyze_refinement_support_v100510,
)
from run_v10_0_5_9_production_j_parity import (
    _authoritative_mpz_length_um,
    _log_tail,
    _probe_command as _base_probe_command,
    _token,
    _values,
    _write_csv,
    build_parser as _base_parser,
)

ENTRY_MODULE = "arrhenius_fracture.mode_i_first_passage_v10_0_5_10_refinement_probe"
LAUNCH_FAILURE_JSON = "production_j_refinement_launch_failure_v10_0_5_10.json"


def build_parser():
    parser = _base_parser()
    parser.description = __doc__
    parser.set_defaults(
        out=Path("runs/v10_0_5_10_refinement_support_v1"),
        contour_outer_um="100 140 180 240 300",
    )
    parser.add_argument("--tip-refinement-radius-um", type=float, default=330.0)
    parser.add_argument("--accepted-contour-um", default="180 240 300")
    parser.add_argument("--contour-stability-rel-tol", type=float, default=0.10)
    parser.add_argument("--radial-edges-um", default="0 60 100 140 180 240 300 330")
    return parser


def _probe_command(args, case_dir: Path, opening_m: float) -> list[str]:
    command = _base_probe_command(args, case_dir, opening_m)
    if command[:2] != [args.python, "-m"]:
        raise RuntimeError("unexpected v10.0.5.9 probe command prefix")
    command[2] = ENTRY_MODULE
    return command


def main(argv: Iterable[str] | None = None) -> int:
    args = build_parser().parse_args(list(argv) if argv is not None else None)
    reference_path = args.reference_json.resolve()
    if not reference_path.exists():
        raise SystemExit(f"missing v10.0.5.8 reference JSON: {reference_path}")
    reference = json.loads(reference_path.read_text())
    if reference.get("schema") != "fixed_grip_elastic_convergence_v10_0_5_8":
        raise SystemExit("--reference-json must be the v10.0.5.8 fixed-grip summary")
    if not bool(reference.get("passed", False)):
        raise SystemExit("the supplied v10.0.5.8 reference did not pass")

    openings = sorted(set(_values(args.openings_um, 1.0e-6)))
    contours = sorted(set(_values(args.contour_outer_um, 1.0e-6)))
    accepted = sorted(set(_values(args.accepted_contour_um, 1.0e-6)))
    radial_edges = sorted(set(_values(args.radial_edges_um, 1.0e-6)))
    radius_m = float(args.tip_refinement_radius_um) * 1.0e-6
    if len(openings) < 2 or any(value <= 0.0 for value in openings):
        raise SystemExit("at least two positive --openings-um values are required")
    if len(contours) < 3 or any(value <= 0.0 for value in contours):
        raise SystemExit("at least three positive contour radii are required")
    if len(accepted) < 2 or any(value <= 0.0 for value in accepted):
        raise SystemExit("at least two positive accepted contours are required")
    if radius_m <= max(accepted):
        raise SystemExit("--tip-refinement-radius-um must exceed every accepted contour")
    if len(radial_edges) < 2 or radial_edges[0] < 0.0:
        raise SystemExit("--radial-edges-um must be increasing and nonnegative")

    mpz_length_um = _authoritative_mpz_length_um(args)
    requested_L_pz_um = None if args.L_pz_um is None else float(args.L_pz_um)
    root = args.out.resolve()
    root.mkdir(parents=True, exist_ok=True)
    probes: list[dict[str, Any]] = []
    processes: list[dict[str, Any]] = []

    env = os.environ.copy()
    env["ARRHENIUS_V10059_CONTOURS_UM"] = " ".join(
        f"{value * 1.0e6:.17g}" for value in contours
    )
    env["ARRHENIUS_V100510_REFINEMENT_RADIUS_M"] = f"{radius_m:.17g}"
    env["ARRHENIUS_V100510_RADIAL_EDGES_UM"] = " ".join(
        f"{value * 1.0e6:.17g}" for value in radial_edges
    )
    env["ARRHENIUS_EVENT_STATISTICS"] = "mean_field"
    env["ARRHENIUS_STOCHASTIC_EMISSION"] = "0"
    env["ARRHENIUS_VHCF_FEM_CACHE"] = "0"

    for opening in openings:
        case_dir = root / f"opening_{_token(opening * 1.0e6)}um"
        case_dir.mkdir(parents=True, exist_ok=True)
        command = _probe_command(args, case_dir, opening)
        log_path = case_dir / "production_j_refinement_probe.log"
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
        log_tail = _log_tail(log_path)
        row = {
            "requested_opening_um": opening * 1.0e6,
            "subprocess_returncode": int(process.returncode),
            "probe_exists": probe_path.exists(),
            "case_directory": str(case_dir),
            "log": str(log_path),
            "command": " ".join(command),
            "log_tail": log_tail,
        }
        processes.append(row)
        if not probe_path.exists():
            failure = {
                "schema": "production_j_refinement_launch_failure_v10_0_5_10",
                "point_release": "10.0.5.10",
                "mechanics_probe_complete": False,
                "tip_refinement_radius_um": args.tip_refinement_radius_um,
                "probe_processes": processes,
            }
            (root / LAUNCH_FAILURE_JSON).write_text(json.dumps(failure, indent=2))
            raise RuntimeError(
                f"refinement probe did not write {probe_path}; return code "
                f"{process.returncode}. Last lines of {log_path}:\n{log_tail}"
            )
        probe = json.loads(probe_path.read_text())
        probe["subprocess_returncode"] = int(process.returncode)
        probe["case_directory"] = str(case_dir)
        probes.append(probe)

    summary = analyze_refinement_support_v100510(
        reference=reference,
        probes=probes,
        selected_outer_radius_m=args.selected_outer_um * 1.0e-6,
        accepted_contours_m=accepted,
        parity_relative_tolerance=args.parity_rel_tol,
        elastic_scaling_relative_tolerance=args.elastic_scaling_rel_tol,
        energy_closure_relative_tolerance=args.energy_closure_rel_tol,
        contour_stability_relative_tolerance=args.contour_stability_rel_tol,
    )
    summary.update(
        {
            "reference_json": str(reference_path),
            "probe_processes": processes,
            "production_configuration": {
                "material_class": args.material_class,
                "temperature_K": args.temperature_K,
                "nx": args.nx,
                "ny": args.ny,
                "tip_h_um": args.tip_h_um,
                "tip_ratio": args.tip_ratio,
                "tip_refinement_radius_um": args.tip_refinement_radius_um,
                "cluster_J_outer_um": args.cluster_J_outer_um,
                "local_J_outer_um": args.local_J_outer_um,
                "authoritative_mpz_length_um": mpz_length_um,
                "legacy_L_pz_um_requested": requested_L_pz_um,
                "mpz_n_bins": args.mpz_n_bins,
                "da_phys_um": args.da_phys_um,
                "crack_backend": args.crack_backend,
                "crystal_anisotropic": args.crystal_aniso,
                "crystal_theta_deg": args.crystal_theta_deg,
            },
            "subprocess_nonzero_count": sum(
                int(row["subprocess_returncode"] != 0) for row in processes
            ),
            "mechanics_probe_complete": all(row["probe_exists"] for row in processes),
        }
    )

    contour_rows = []
    radial_rows = []
    for probe in probes:
        common = {
            "Uapp_um": float(probe["Uapp_um"]),
            "sigma_gross_MPa": float(probe["sigma_gross_MPa"]),
            "case_directory": probe["case_directory"],
        }
        contour_rows.extend({**common, **dict(row)} for row in probe.get("contours", []))
        radial_rows.extend({**common, **dict(row)} for row in probe.get("radial_mesh_support", []))

    _write_csv(root / PROBE_CSV, list(summary["base_v10_0_5_9_analysis"].get("cases", [])))
    _write_csv(root / CONTOUR_CSV, contour_rows)
    _write_csv(root / RADIAL_CSV, radial_rows)
    (root / SUMMARY_JSON).write_text(json.dumps(summary, indent=2, default=str))
    print(f"PRODUCTION REFINEMENT SUPPORT STATUS: {summary['status']}")
    print(root / SUMMARY_JSON)
    return 0 if bool(summary.get("passed", False)) else 2


if __name__ == "__main__":
    raise SystemExit(main())
