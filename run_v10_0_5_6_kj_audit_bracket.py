#!/usr/bin/env python3
"""Audit the primary J contour and find a stochastic first-passage load bracket.

Modes
-----
audit
    Run one-step elastic calibrations over candidate actual J outer radii, reject
    contours that intersect a specimen boundary, compare KJ/sigma with a finite-
    width single-edge-crack reference, and select a stable contour plateau.
bracket
    Read the selected contour, convert target KJ values to remote stress ranges,
    run ascending stochastic cases, expand as needed, and bisect the interval
    between no first passage and first passage at the selected cycle horizon.
all
    Perform audit and then bracket.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
import shlex
import subprocess
import sys
from typing import Any

from arrhenius_fracture.kj_audit_v10056 import (
    POINT_RELEASE,
    SpecimenGeometryV10056,
    classify_first_passage_rows,
    contour_geometry_audit,
    select_contour_plateau,
)


ROOT = Path(__file__).resolve().parent
CAMPAIGN = ROOT / "run_v10_0_5_6_stochastic_delta_sigma.py"
AUDIT_CSV = "KJ_contour_sweep_v10_0_5_6.csv"
AUDIT_JSON = "KJ_contour_audit_v10_0_5_6.json"
SELECTED_JSON = "selected_KJ_contour_v10_0_5_6.json"
BRACKET_CSV = "first_passage_bracket_cases_v10_0_5_6.csv"
BRACKET_JSON = "first_passage_bracket_v10_0_5_6.json"


def _float_token(value: float) -> str:
    return f"{value:.10g}".replace("-", "m").replace(".", "p")


def _read_csv(path: Path) -> list[dict[str, Any]]:
    with path.open(newline="") as handle:
        rows = []
        for raw in csv.DictReader(handle):
            row: dict[str, Any] = {}
            for key, value in raw.items():
                try:
                    row[key] = float(value)
                except (TypeError, ValueError):
                    row[key] = value
            rows.append(row)
        return rows


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


def _run(command: list[str], *, env: dict[str, str], log: Path) -> None:
    log.parent.mkdir(parents=True, exist_ok=True)
    print(f"\nRUN: {shlex.join(command)}", flush=True)
    with log.open("w") as handle:
        process = subprocess.Popen(
            command,
            cwd=ROOT,
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        if process.stdout is None:
            raise RuntimeError("child stdout pipe unavailable")
        for line in process.stdout:
            handle.write(line)
            handle.flush()
            print(line, end="", flush=True)
        code = process.wait()
    if code != 0:
        tail = "\n".join(log.read_text(errors="replace").splitlines()[-100:])
        raise RuntimeError(
            f"command failed with exit code {code}: {shlex.join(command)}\n{tail}"
        )


def _common_campaign_command(args, *, out: Path, delta_sigma_MPa: float, contour_um: float):
    return [
        args.python,
        str(CAMPAIGN),
        "--out",
        str(out),
        "--material-class",
        args.material_class,
        "--temperatures",
        f"{args.temperature_K:g}",
        "--delta-sigma-MPa",
        f"{delta_sigma_MPa:.17g}",
        "--R",
        f"{args.R:.17g}",
        "--frequency-Hz",
        f"{args.frequency_Hz:.17g}",
        "--cluster-J-outer-um",
        f"{contour_um:.17g}",
        "--specimen-width-m",
        f"{args.specimen_width_m:.17g}",
        "--specimen-height-m",
        f"{args.specimen_height_m:.17g}",
        "--initial-crack-m",
        f"{args.initial_crack_m:.17g}",
        "--nx",
        str(args.nx),
        "--ny",
        str(args.ny),
        "--tip-h-fine-m",
        f"{args.tip_h_fine_m:.17g}",
        "--tip-ratio",
        f"{args.tip_ratio:.17g}",
        "--n-phase",
        str(args.n_phase),
        "--save-snapshots",
        "0",
        "--print-every",
        str(args.print_every),
    ]


def _front_active_elements(case_directory: Path, temperature_K: float) -> int | None:
    exact = case_directory / f"fronts_{int(round(temperature_K)):04d}K.csv"
    matches = [exact] if exact.exists() else sorted(case_directory.glob("fronts_*K.csv"))
    if not matches:
        return None
    rows = _read_csv(matches[0])
    if not rows:
        return None
    return int(float(rows[0].get("J_active_elems", 0)))


def run_audit(args) -> dict[str, Any]:
    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    geometry = SpecimenGeometryV10056(
        width_m=args.specimen_width_m,
        height_m=args.specimen_height_m,
        initial_crack_m=args.initial_crack_m,
    ).validate()
    rows: list[dict[str, Any]] = []
    env = os.environ.copy()
    env.update(
        {
            "ARRHENIUS_EVENT_STATISTICS": "deterministic",
            "ARRHENIUS_STOCHASTIC_EMISSION": "0",
            "ARRHENIUS_STOCHASTIC_BLOCKS": "0",
            "ARRHENIUS_VHCF_FEM_CACHE": "0",
        }
    )

    for outer_um in sorted(set(float(v) for v in args.contour_outer_um)):
        geom = contour_geometry_audit(
            outer_radius_m=outer_um * 1.0e-6,
            geometry=geometry,
            safety_fraction=args.J_contour_safety_fraction,
        )
        if not geom["contour_closes_inside_body"]:
            rows.append(
                {
                    **geom,
                    "outer_radius_um": outer_um,
                    "audit_run_status": "rejected_before_run_boundary_intersection",
                }
            )
            continue
        case_root = out / "contours" / f"outer_{_float_token(outer_um)}um"
        command = _common_campaign_command(
            args,
            out=case_root,
            delta_sigma_MPa=args.audit_delta_sigma_MPa,
            contour_um=outer_um,
        )
        command.extend(
            [
                "--cycles-max",
                "1e-9",
                "--block-cycles",
                "1e-9",
                "--max-block-cycles",
                "1e-9",
                "--max-blocks",
                "1",
                "--target-extension-um",
                "1e-12",
                "--target-dN-store",
                "0.05",
                "--target-dN-emit",
                "inf",
                "--target-dN-mobile",
                "inf",
                "--target-dN-escape",
                "0.25",
            ]
        )
        if not (case_root / "K_vs_delta_sigma_v10_0_5_6.csv").exists():
            _run(command, env=env, log=case_root / "launcher.log")
        Krows = _read_csv(case_root / "remote_stress_KJ_audit_v10_0_5_6.csv")
        if len(Krows) != 1:
            raise RuntimeError(f"expected one K audit row in {case_root}")
        row = dict(Krows[0])
        case_directory = Path(
            _read_csv(case_root / "K_vs_delta_sigma_v10_0_5_6.csv")[0][
                "run_directory"
            ]
        )
        active = _front_active_elements(case_directory, args.temperature_K)
        row.update(
            {
                "outer_radius_um": outer_um,
                "J_active_elements": active if active is not None else 0,
                "audit_run_status": "complete",
                "case_root": str(case_root),
            }
        )
        rows.append(row)

    selection = select_contour_plateau(
        rows,
        relative_tolerance=args.plateau_relative_tolerance,
        minimum_points=args.plateau_minimum_points,
        minimum_active_elements=args.minimum_J_active_elements,
    )
    legacy = contour_geometry_audit(
        outer_radius_m=8.0e-3,
        geometry=geometry,
        safety_fraction=args.J_contour_safety_fraction,
    )
    payload = {
        "schema": "KJ_contour_audit_v10_0_5_6",
        "point_release": POINT_RELEASE,
        "temperature_K": args.temperature_K,
        "audit_delta_sigma_MPa": args.audit_delta_sigma_MPa,
        "geometry": {
            "width_m": geometry.width_m,
            "height_m": geometry.height_m,
            "initial_crack_m": geometry.initial_crack_m,
            "a_over_W": geometry.a_over_W,
            "nearest_tip_boundary_m": geometry.nearest_tip_boundary_m,
        },
        "legacy_default_cluster_outer_radius_m": 8.0e-3,
        "legacy_default_geometry_audit": legacy,
        "selection": selection,
        "constitutive_physics_changed": False,
    }
    _write_csv(out / AUDIT_CSV, rows)
    (out / AUDIT_JSON).write_text(json.dumps(payload, indent=2, default=str))
    (out / SELECTED_JSON).write_text(json.dumps(selection, indent=2, default=str))
    print(f"\nKJ AUDIT STATUS: {selection['status']}")
    if selection.get("selected_outer_radius_m") is not None:
        print(
            "selected cluster J outer radius = "
            f"{1e6 * float(selection['selected_outer_radius_m']):.6g} um"
        )
    print(out / AUDIT_CSV)
    return payload


def _target_K_to_delta_sigma_MPa(
    target_K_MPa_sqrt_m: float,
    KJ_per_sigma_sqrt_m: float,
    R: float,
) -> float:
    sigma_max_Pa = target_K_MPa_sqrt_m * 1.0e6 / KJ_per_sigma_sqrt_m
    return sigma_max_Pa * (1.0 - R) / 1.0e6


def run_bracket(args) -> dict[str, Any]:
    out = Path(args.out).resolve()
    selected_path = Path(args.selected_contour_json).resolve()
    selected = json.loads(selected_path.read_text())
    if selected.get("status") != "plateau_selected":
        raise SystemExit(
            f"cannot bracket first passage without a selected KJ plateau: {selected_path}"
        )
    contour_um = 1.0e6 * float(selected["selected_outer_radius_m"])
    slope = float(selected["selected_row"]["KJ_per_sigma_gross_sqrt_m"])
    if not math.isfinite(slope) or slope <= 0.0:
        raise SystemExit("selected KJ/sigma slope is invalid")

    env = os.environ.copy()
    env.update(
        {
            "ARRHENIUS_EVENT_STATISTICS": "stochastic",
            "ARRHENIUS_STOCHASTIC_EMISSION": "1",
            "ARRHENIUS_STOCHASTIC_BLOCKS": "1",
            "ARRHENIUS_STOCHASTIC_SEED": str(args.stochastic_seed),
            "ARRHENIUS_RARE_EVENT_TARGET": str(args.rare_event_target),
            "ARRHENIUS_TAU_LEAP_TARGET": str(args.tau_leap_target),
            "ARRHENIUS_TAU_SWITCH_EXPECTED_EVENTS": str(
                args.tau_switch_expected_events
            ),
            "ARRHENIUS_VHCF_FEM_CACHE": "0",
        }
    )
    cases: dict[str, dict[str, Any]] = {}

    def evaluate(target_K: float) -> dict[str, Any]:
        delta_sigma = _target_K_to_delta_sigma_MPa(target_K, slope, args.R)
        token = _float_token(delta_sigma)
        if token in cases:
            return cases[token]
        case_root = out / "bracket_cases" / f"DeltaSigma_{token}MPa"
        command = _common_campaign_command(
            args,
            out=case_root,
            delta_sigma_MPa=delta_sigma,
            contour_um=contour_um,
        )
        command.extend(
            [
                "--cycles-max",
                f"{args.cycles_max:.17g}",
                "--block-cycles",
                "1e4",
                "--max-block-cycles",
                "inf",
                "--max-blocks",
                str(args.max_blocks),
                "--target-extension-um",
                f"{args.target_extension_um:.17g}",
                "--target-dB",
                f"{args.target_dB:.17g}",
                "--target-dN-store",
                "0.05",
                "--target-dN-emit",
                "inf",
                "--target-dN-mobile",
                "inf",
                "--target-dN-escape",
                "0.25",
            ]
        )
        result_file = case_root / "K_vs_delta_sigma_v10_0_5_6.csv"
        if not result_file.exists():
            _run(command, env=env, log=case_root / "launcher.log")
        rows = _read_csv(result_file)
        if len(rows) != 1:
            raise RuntimeError(f"expected one bracket case row in {result_file}")
        row = dict(rows[0])
        row["target_KJmax_MPa_sqrt_m"] = target_K
        row["selected_cluster_J_outer_um"] = contour_um
        row["stochastic_seed"] = args.stochastic_seed
        cases[token] = row
        return row

    targets = sorted(set(float(value) for value in args.target_KJmax_MPa_sqrt_m))
    for target in targets:
        evaluate(target)
        state = classify_first_passage_rows(cases.values())
        if state["status"] == "bracketed":
            break

    state = classify_first_passage_rows(cases.values())
    while state["status"] != "bracketed" and len(cases) < args.maximum_cases:
        if state["n_first_passage"] == 0:
            next_K = max(float(row["target_KJmax_MPa_sqrt_m"]) for row in cases.values()) * args.expansion_factor
            if next_K > args.maximum_target_KJmax_MPa_sqrt_m:
                break
        else:
            next_K = min(float(row["target_KJmax_MPa_sqrt_m"]) for row in cases.values()) / args.expansion_factor
            if next_K < args.minimum_target_KJmax_MPa_sqrt_m:
                break
        evaluate(next_K)
        state = classify_first_passage_rows(cases.values())

    for _ in range(args.bisection_refinements):
        state = classify_first_passage_rows(cases.values())
        if state["status"] != "bracketed":
            break
        low = state["lower_no_first_passage"]
        high = state["upper_first_passage"]
        midpoint_K = 0.5 * (
            float(low["target_KJmax_MPa_sqrt_m"])
            + float(high["target_KJmax_MPa_sqrt_m"])
        )
        evaluate(midpoint_K)

    ordered = sorted(
        cases.values(), key=lambda row: float(row["delta_sigma_requested_MPa"])
    )
    state = classify_first_passage_rows(ordered)
    payload = {
        "schema": "first_passage_bracket_v10_0_5_6",
        "point_release": POINT_RELEASE,
        "temperature_K": args.temperature_K,
        "material_class": args.material_class,
        "cycles_max": args.cycles_max,
        "stochastic_seed": args.stochastic_seed,
        "selected_contour_json": str(selected_path),
        "selected_cluster_J_outer_um": contour_um,
        "selected_KJ_per_sigma_sqrt_m": slope,
        "bracket": state,
        "n_cases": len(ordered),
        "constitutive_physics_changed": False,
    }
    _write_csv(out / BRACKET_CSV, ordered)
    (out / BRACKET_JSON).write_text(json.dumps(payload, indent=2, default=str))
    print(f"\nFIRST-PASSAGE BRACKET STATUS: {state['status']}")
    if state.get("stress_interval_MPa") is not None:
        lo, hi = state["stress_interval_MPa"]
        print(f"Delta-sigma bracket = [{lo:.6g}, {hi:.6g}] MPa")
    print(out / BRACKET_CSV)
    return payload


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=["audit", "bracket", "all"])
    parser.add_argument("--out", required=True)
    parser.add_argument("--material-class", choices=["DBTT", "weakT", "ceramic"], default="DBTT")
    parser.add_argument("--temperature-K", type=float, default=700.0)
    parser.add_argument("--R", type=float, default=0.1)
    parser.add_argument("--frequency-Hz", type=float, default=1000.0, dest="frequency_Hz")
    parser.add_argument("--specimen-width-m", type=float, default=2.0e-3)
    parser.add_argument("--specimen-height-m", type=float, default=4.0e-3)
    parser.add_argument("--initial-crack-m", type=float, default=0.5e-3)
    parser.add_argument("--nx", type=int, default=40)
    parser.add_argument("--ny", type=int, default=80)
    parser.add_argument("--tip-h-fine-m", type=float, default=2.5e-6)
    parser.add_argument("--tip-ratio", type=float, default=1.2)
    parser.add_argument("--n-phase", type=int, default=96)
    parser.add_argument("--print-every", type=int, default=10)
    parser.add_argument("--python", default=sys.executable)

    parser.add_argument(
        "--contour-outer-um",
        type=float,
        nargs="+",
        default=[60, 80, 100, 140, 180, 240, 300, 360, 400],
    )
    parser.add_argument("--audit-delta-sigma-MPa", type=float, default=100.0)
    parser.add_argument("--J-contour-safety-fraction", type=float, default=0.80, dest="J_contour_safety_fraction")
    parser.add_argument("--plateau-relative-tolerance", type=float, default=0.10)
    parser.add_argument("--plateau-minimum-points", type=int, default=3)
    parser.add_argument("--minimum-J-active-elements", type=int, default=12)

    parser.add_argument("--selected-contour-json")
    parser.add_argument(
        "--target-KJmax-MPa-sqrt-m",
        type=float,
        nargs="+",
        default=[2, 4, 6, 8, 10, 12, 16, 20, 24],
    )
    parser.add_argument("--cycles-max", type=float, default=1.0e7)
    parser.add_argument("--max-blocks", type=int, default=10000)
    parser.add_argument("--target-extension-um", type=float, default=5.0)
    parser.add_argument("--target-dB", type=float, default=0.01)
    parser.add_argument("--stochastic-seed", type=int, default=1)
    parser.add_argument("--rare-event-target", type=float, default=0.25)
    parser.add_argument("--tau-leap-target", type=float, default=3.0)
    parser.add_argument("--tau-switch-expected-events", type=float, default=10.0)
    parser.add_argument("--expansion-factor", type=float, default=1.5)
    parser.add_argument("--minimum-target-KJmax-MPa-sqrt-m", type=float, default=0.25)
    parser.add_argument("--maximum-target-KJmax-MPa-sqrt-m", type=float, default=60.0)
    parser.add_argument("--maximum-cases", type=int, default=20)
    parser.add_argument("--bisection-refinements", type=int, default=3)
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    out = Path(args.out).resolve()
    if args.mode in {"audit", "all"}:
        run_audit(args)
    if args.mode in {"bracket", "all"}:
        if args.selected_contour_json is None:
            args.selected_contour_json = str(out / SELECTED_JSON)
        run_bracket(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
