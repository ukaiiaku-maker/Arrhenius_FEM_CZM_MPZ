#!/usr/bin/env python3
"""v10.0.5.6 stochastic delta-sigma runner with an audited KJ contour.

This wrapper keeps the v10.0.5.5 stochastic source/cleavage implementation but
requires an actual cluster-J outer radius that closes inside the specimen.  It
also replaces ambiguous compatibility diagnostics with authoritative stochastic
scheduler and moving-process-zone labels.
"""
from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
from typing import Iterable

import run_v10_0_5_4_vhcf_delta_sigma as _base

from arrhenius_fracture.kj_audit_v10056 import (
    POINT_RELEASE,
    SpecimenGeometryV10056,
    build_kj_audit_row,
    contour_geometry_audit,
    enrich_stochastic_block_rows,
)


ENTRY_MODULE = (
    "arrhenius_fracture."
    "mode_i_first_passage_v10_0_5_5_stochastic_vhcf_audited"
)
COMPLETION_MANIFEST = "run_completion_v10_0_5_5_stochastic_vhcf.json"
OUTPUT_BLOCKS = "fatigue_block_diagnostics_v10_0_5_6.csv"
OUTPUT_CASES = "K_vs_delta_sigma_v10_0_5_6.csv"
OUTPUT_K_AUDIT = "remote_stress_KJ_audit_v10_0_5_6.csv"
OUTPUT_MANIFEST = "campaign_manifest_v10_0_5_6.json"

_ORIGINAL_BUILD_PARSER = _base.build_parser
_ORIGINAL_BASE_COMMAND = _base._base_command


def build_parser() -> argparse.ArgumentParser:
    parser = _ORIGINAL_BUILD_PARSER()
    parser.description = __doc__
    parser.set_defaults(cycles_max=1.0e8, max_blocks=10000)
    parser.add_argument(
        "--cluster-J-outer-um",
        type=float,
        required=True,
        dest="cluster_J_outer_um",
        help=(
            "Actual outer radius of the primary/cluster J domain in micrometers. "
            "Use the selected value from the v10.0.5.6 contour audit."
        ),
    )
    parser.add_argument("--specimen-height-m", type=float, default=4.0e-3)
    parser.add_argument("--initial-crack-m", type=float, default=0.5e-3)
    parser.add_argument(
        "--J-contour-safety-fraction",
        type=float,
        default=0.80,
        dest="J_contour_safety_fraction",
    )
    parser.add_argument(
        "--KJ-LEFM-ratio-min",
        type=float,
        default=0.70,
        dest="KJ_LEFM_ratio_min",
    )
    parser.add_argument(
        "--KJ-LEFM-ratio-max",
        type=float,
        default=1.30,
        dest="KJ_LEFM_ratio_max",
    )
    parser.add_argument(
        "--fail-on-KJ-mismatch",
        action="store_true",
        dest="fail_on_KJ_mismatch",
    )
    return parser


def _base_command_v10056(args, outdir: Path, temperature: float, dU_m: float):
    command = list(_ORIGINAL_BASE_COMMAND(args, outdir, temperature, dU_m))
    command[command.index(_base.ENTRY_MODULE)] = ENTRY_MODULE
    cluster_outer_m = float(args.cluster_J_outer_um) * 1.0e-6
    command.extend(["--rJ-cluster", f"{cluster_outer_m / 8.0:.17g}"])
    return command


def _case_audit_path(case_row: dict) -> Path:
    return Path(str(case_row["run_directory"])) / "stochastic_vhcf_v10_0_5_5.json"


def _postprocess(root: Path, args) -> tuple[list[dict], list[dict], list[dict]]:
    cases = _base._read_rows(root / "K_vs_delta_sigma.csv")
    blocks = _base._read_rows(root / "fatigue_block_diagnostics_v10_0_5_4.csv")
    enriched_cases: list[dict] = []
    enriched_blocks: list[dict] = []
    K_audit_rows: list[dict] = []
    geometry = SpecimenGeometryV10056(
        width_m=float(args.specimen_width_m),
        height_m=float(args.specimen_height_m),
        initial_crack_m=float(args.initial_crack_m),
    ).validate()
    cluster_outer_m = float(args.cluster_J_outer_um) * 1.0e-6

    for case in cases:
        temperature = float(case["temperature_K"])
        delta_sigma = float(case["delta_sigma_requested_MPa"])
        audit_path = _case_audit_path(case)
        if not audit_path.exists():
            raise RuntimeError(f"missing stochastic audit: {audit_path}")
        audit = json.loads(audit_path.read_text())
        scheduler_records = list(audit.get("scheduler", {}).get("records", []))
        engines = list(audit.get("engines", []))
        if len(engines) != 1:
            raise RuntimeError(
                f"expected one single-front engine in {audit_path}; found {len(engines)}"
            )
        engine = dict(engines[0])
        selected = [
            row
            for row in blocks
            if math.isclose(float(row["temperature_K"]), temperature)
            and math.isclose(
                float(row["delta_sigma_requested_MPa"]), delta_sigma
            )
        ]
        enriched_blocks.extend(
            enrich_stochastic_block_rows(selected, scheduler_records)
        )

        updated = dict(case)
        updated.update(
            {
                "point_release": POINT_RELEASE,
                "cluster_J_outer_um": float(args.cluster_J_outer_um),
                "cleavage_threshold": float(engine.get("cleavage_threshold", math.nan)),
                "cleavage_event_index": int(engine.get("cleavage_event_index", 0)),
                "source_budget_total": float(engine.get("source_budget_total", math.nan)),
                "source_consumed_final": float(engine.get("source_budget_consumed", math.nan)),
                "source_remaining_final": float(engine.get("source_budget_remaining", math.nan)),
                "mobile_count_final": float(engine.get("mobile_count", math.nan)),
                "retained_count_final": float(engine.get("retained_count", math.nan)),
                "active_count_final": float(engine.get("active_count", math.nan)),
                "cumulative_emitted": float(engine.get("cumulative_emitted", math.nan)),
                # This is a count of nonzero system/channel realizations, not emitted sites.
                "stochastic_emission_channel_events": float(
                    engine.get("stochastic_emission_events", math.nan)
                ),
                "predictor_mean_field_calls": int(
                    engine.get("predictor_mean_field_calls", 0)
                ),
            }
        )
        enriched_cases.append(updated)

        sigma_actual_Pa = float(case["sigma_max_actual_MPa_first"]) * 1.0e6
        Ftop = sigma_actual_Pa * geometry.width_m
        Krow = build_kj_audit_row(
            Ftop_N_per_thickness=Ftop,
            KJ_Pa_sqrt_m=float(case["KJmax_first_MPa_sqrt_m"]) * 1.0e6,
            outer_radius_m=cluster_outer_m,
            geometry=geometry,
            safety_fraction=float(args.J_contour_safety_fraction),
        )
        Krow.update(
            {
                "temperature_K": temperature,
                "delta_sigma_requested_MPa": delta_sigma,
                "run_directory": case["run_directory"],
            }
        )
        K_audit_rows.append(Krow)

    _base._write_csv(root / OUTPUT_BLOCKS, enriched_blocks)
    _base._write_csv(root / OUTPUT_CASES, enriched_cases)
    _base._write_csv(root / OUTPUT_K_AUDIT, K_audit_rows)
    return enriched_cases, enriched_blocks, K_audit_rows


def main(argv: Iterable[str] | None = None) -> int:
    arg_list = list(argv) if argv is not None else None
    args = build_parser().parse_args(arg_list)
    geometry = SpecimenGeometryV10056(
        width_m=float(args.specimen_width_m),
        height_m=float(args.specimen_height_m),
        initial_crack_m=float(args.initial_crack_m),
    ).validate()
    contour = contour_geometry_audit(
        outer_radius_m=float(args.cluster_J_outer_um) * 1.0e-6,
        geometry=geometry,
        safety_fraction=float(args.J_contour_safety_fraction),
    )
    if not contour["contour_within_safety_limit"]:
        raise SystemExit(
            "cluster J contour is not geometrically valid: "
            f"outer={contour['outer_radius_m']:.6g} m, "
            f"safe_limit={contour['safe_outer_radius_limit_m']:.6g} m"
        )

    os.environ.setdefault("ARRHENIUS_EVENT_STATISTICS", "stochastic")
    os.environ.setdefault("ARRHENIUS_STOCHASTIC_EMISSION", "1")
    os.environ.setdefault("ARRHENIUS_STOCHASTIC_BLOCKS", "1")
    os.environ.setdefault("ARRHENIUS_RARE_EVENT_TARGET", "0.25")
    os.environ.setdefault("ARRHENIUS_TAU_LEAP_TARGET", "3.0")
    os.environ.setdefault("ARRHENIUS_TAU_SWITCH_EXPECTED_EVENTS", "10.0")
    os.environ.setdefault("ARRHENIUS_VHCF_FEM_CACHE", "0")

    saved_build = _base.build_parser
    saved_command = _base._base_command
    saved_entry = _base.ENTRY_MODULE
    saved_release = _base.POINT_RELEASE
    saved_completion = _base.COMPLETION_MANIFEST
    _base.build_parser = build_parser
    _base._base_command = _base_command_v10056
    _base.ENTRY_MODULE = ENTRY_MODULE
    _base.POINT_RELEASE = POINT_RELEASE
    _base.COMPLETION_MANIFEST = COMPLETION_MANIFEST
    try:
        status = int(_base.main(arg_list) or 0)
    finally:
        _base.build_parser = saved_build
        _base._base_command = saved_command
        _base.ENTRY_MODULE = saved_entry
        _base.POINT_RELEASE = saved_release
        _base.COMPLETION_MANIFEST = saved_completion

    root = Path(args.out).resolve()
    cases, blocks, Krows = _postprocess(root, args)
    mismatches = [
        row
        for row in Krows
        if not (
            float(args.KJ_LEFM_ratio_min)
            <= float(row["KJ_over_K_LEFM_gross"])
            <= float(args.KJ_LEFM_ratio_max)
        )
    ]
    manifest = {
        "schema": "v10_0_5_6_stochastic_delta_sigma_campaign",
        "point_release": POINT_RELEASE,
        "cluster_J_outer_um": float(args.cluster_J_outer_um),
        "contour_geometry": contour,
        "KJ_LEFM_ratio_acceptance": [
            float(args.KJ_LEFM_ratio_min),
            float(args.KJ_LEFM_ratio_max),
        ],
        "KJ_mismatch_count": len(mismatches),
        "n_cases": len(cases),
        "n_blocks": len(blocks),
        "outputs": {
            "cases": OUTPUT_CASES,
            "blocks": OUTPUT_BLOCKS,
            "KJ_audit": OUTPUT_K_AUDIT,
        },
        "constitutive_physics_changed": False,
    }
    (root / OUTPUT_MANIFEST).write_text(json.dumps(manifest, indent=2))
    if mismatches:
        print(
            "WARNING: KJ/LEFM mismatch outside accepted range for "
            f"{len(mismatches)} case(s); inspect {root / OUTPUT_K_AUDIT}"
        )
        if args.fail_on_KJ_mismatch:
            return 3
    print(root / OUTPUT_CASES)
    return status


if __name__ == "__main__":
    raise SystemExit(main())
