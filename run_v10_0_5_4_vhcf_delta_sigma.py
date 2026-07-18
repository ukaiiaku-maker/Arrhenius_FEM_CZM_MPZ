#!/usr/bin/env python3
"""v10.0.5.4 VHCF remote-stress fatigue campaign.

The campaign preserves the v10.0.5.3 FEM/CZM/MPZ constitutive implementation
but removes artificial cycle-jump ceilings, uses the authoritative engine
predictor, distinguishes physical completion from numerical censoring, and
writes block increments with unambiguous per-block/per-cycle semantics.
"""
from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import shutil
import sys
from typing import Iterable

import run_v10_0_5_3_delta_sigma_fatigue as _base

from arrhenius_fracture.mode_i_first_passage_v10_0_5_4_vhcf import (
    COMPLETION_MANIFEST,
    classify_termination_v10054,
)

POINT_RELEASE = "10.0.5.4"
ENTRY_MODULE = "arrhenius_fracture.mode_i_first_passage_v10_0_5_4_vhcf"

_read_rows = _base._read_rows
_find_steps = _base._find_steps
_write_csv = _base._write_csv
_make_plots = _base._make_plots
_run = _base._run


def _remove_option_pair(command: list[str], option: str) -> None:
    while option in command:
        index = command.index(option)
        del command[index:index + 2]


def _base_command(args, outdir: Path, temperature: float, dU_m: float) -> list[str]:
    command = list(_base._base_command(args, outdir, temperature, dU_m))
    command[command.index("arrhenius_fracture.mode_i_first_passage_v10_0_5_3_fatigue")] = (
        ENTRY_MODULE
    )

    # The inherited parser has no separate target-da option. Temporal accuracy is
    # controlled by the finite hazard/state increments below.
    _remove_option_pair(command, "--target-da-per-block-um")

    if "--crystal-compete" not in command:
        index = command.index("--crystal-aniso") + 1
        command.insert(index, "--crystal-compete")

    command.extend(
        [
            "--target-dN-escape",
            f"{args.target_dN_escape:.17g}",
            "--target-dN-peierls",
            f"{args.target_dN_peierls:.17g}",
            "--target-dN-taylor",
            f"{args.target_dN_taylor:.17g}",
        ]
    )
    if not args.resolve_cyclic_mechanics:
        command.append("--no-cyclic-mechanics")
    return command


def build_parser() -> argparse.ArgumentParser:
    parser = _base.build_parser()
    parser.description = __doc__
    parser.set_defaults(
        max_block_cycles=math.inf,
        max_blocks=500,
        cycles_max=1.0e14,
    )
    parser.add_argument(
        "--target-dN-escape",
        type=float,
        default=0.10,
        dest="target_dN_escape",
    )
    parser.add_argument(
        "--target-dN-peierls",
        type=float,
        default=math.inf,
        dest="target_dN_peierls",
    )
    parser.add_argument(
        "--target-dN-taylor",
        type=float,
        default=math.inf,
        dest="target_dN_taylor",
    )
    parser.add_argument(
        "--resolve-cyclic-mechanics",
        action="store_true",
        help=(
            "Resolve phase-by-phase FEM cyclic mechanics. The default tip-only "
            "VHCF path reuses the maximum-load FEM/tensor field and integrates "
            "the local hazards over the waveform."
        ),
    )
    parser.add_argument(
        "--fail-on-censor",
        action="store_true",
        help="Return a nonzero campaign status if any case exhausts max blocks.",
    )
    return parser


def _corrected_block_rows(
    raw_rows: list[dict[str, float]],
    *,
    temperature: float,
    delta_sigma_MPa: float,
) -> list[dict]:
    output: list[dict] = []
    cycles_cumulative = 0.0
    for raw in raw_rows:
        cycles = max(float(raw.get("fatigue_cycles", 0.0)), 0.0)
        cycles_cumulative += cycles

        # The inherited step table labels these three quantities "per_cycle",
        # but the progressive formatter supplies accepted block increments.
        dN_emit_block = max(float(raw.get("mu_emit_per_cycle", 0.0)), 0.0)
        dB_predictor_block = max(
            float(raw.get("mu_cleave_pred_per_cycle", 0.0)), 0.0
        )
        dN_escape_predictor_block = max(
            float(raw.get("mu_escape_per_cycle", 0.0)), 0.0
        )
        output.append(
            {
                "temperature_K": temperature,
                "delta_sigma_requested_MPa": delta_sigma_MPa,
                "step": int(float(raw.get("step", len(output) + 1))),
                "cycles_block": cycles,
                "cycles_cumulative": cycles_cumulative,
                "cycle_limiter_code": int(
                    float(raw.get("cycle_limiter_code", -1))
                ),
                "cycle_unlimited": float(raw.get("cycle_unlimited", cycles)),
                "dB_committed_block": float(raw.get("dB_block", 0.0)),
                "dB_predictor_block": dB_predictor_block,
                "dB_predictor_per_cycle": (
                    dB_predictor_block / cycles if cycles > 0.0 else 0.0
                ),
                "dN_emit_block": dN_emit_block,
                "dN_emit_per_cycle": (
                    dN_emit_block / cycles if cycles > 0.0 else 0.0
                ),
                "dN_store_block": float(raw.get("dN_store_block", 0.0)),
                "dN_mobile_block": float(raw.get("dN_mobile_block", 0.0)),
                "dN_escape_block": float(raw.get("dN_escape_block", 0.0)),
                "dN_escape_predictor_block": dN_escape_predictor_block,
                "dN_escape_predictor_per_cycle": (
                    dN_escape_predictor_block / cycles if cycles > 0.0 else 0.0
                ),
                "B": float(raw.get("B", 0.0)),
                "n_fire": int(float(raw.get("n_fire", 0.0))),
                "KJmax_MPa_sqrt_m": float(raw.get("KJ_Pa_sqrtm", 0.0))
                / 1.0e6,
                "crack_extension_um": float(
                    raw.get("crack_extension_m", 0.0)
                )
                * 1.0e6,
                "mobile_count": float(raw.get("mpz_mobile_count", 0.0)),
                "retained_count": float(raw.get("mpz_retained_count", 0.0)),
                "available_site_fraction": float(
                    raw.get("mpz_available_site_fraction", 1.0)
                ),
            }
        )
    return output


def _summarize_case_v10054(
    outdir: Path,
    temperature: float,
    delta_sigma_MPa: float,
    args,
) -> tuple[dict, list[dict], list[dict]]:
    case, local = _base._summarize_case(
        outdir,
        temperature,
        delta_sigma_MPa,
        args.R,
        args.specimen_width_m,
    )
    raw_rows = _read_rows(_find_steps(outdir, temperature))
    termination = classify_termination_v10054(
        raw_rows,
        cycles_max=float(args.cycles_max),
        max_blocks=int(args.max_blocks),
        target_extension_um=float(args.target_extension_um),
    )
    case.update(
        {
            "point_release": POINT_RELEASE,
            "physical_status": termination["status"],
            "termination": termination["termination"],
            "right_censored": termination["right_censored"],
            "reached_cycle_horizon": termination["reached_cycle_horizon"],
            "reached_target_extension": termination[
                "reached_target_extension"
            ],
            "first_passage_observed": termination["first_passage_observed"],
            "cycle_horizon_fraction": termination["cycle_horizon_fraction"],
            "max_block_cycles": args.max_block_cycles,
            "cyclic_mechanics_resolved": bool(args.resolve_cyclic_mechanics),
        }
    )
    corrected = _corrected_block_rows(
        raw_rows,
        temperature=temperature,
        delta_sigma_MPa=delta_sigma_MPa,
    )
    return case, local, corrected


def main(argv: Iterable[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.R >= 1.0:
        raise SystemExit("R must be less than 1 for a positive stress range")
    if (
        math.isfinite(args.cycles_max)
        and math.isfinite(args.max_block_cycles)
        and args.max_block_cycles < args.cycles_max
    ):
        raise SystemExit(
            "MAX_BLOCK_CYCLES must be inf or at least CYCLES_MAX for the "
            "v10.0.5.4 VHCF first-passage runner"
        )

    root = Path(args.out).resolve()
    root.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env.update(
        {
            "ARRHENIUS_CZM_OPENING_COUPLING": "clock_linear",
            "ARRHENIUS_MAX_TRIAL_DAMAGE_CHANGE": str(min(args.target_dB, 0.05)),
            "ARRHENIUS_COMMITTED_TARGET_EXTENSION_UM": str(
                args.target_extension_um
            ),
        }
    )

    case_rows: list[dict] = []
    local_rows: list[dict] = []
    corrected_block_rows: list[dict] = []
    calibration_rows: list[dict] = []

    for temperature in args.temperatures:
        calibration = root / "calibration" / f"T{int(round(temperature)):04d}K"
        if calibration.exists() and not args.keep_existing:
            shutil.rmtree(calibration)
        calibration_args = argparse.Namespace(**vars(args))
        calibration_args.max_blocks = 1
        calibration_args.cycles_max = 1.0e-9
        calibration_args.block_cycles = 1.0e-9
        calibration_args.max_block_cycles = 1.0e-9
        calibration_args.target_extension_um = 1.0e-12
        calibration_args.resolve_cyclic_mechanics = False
        command = _base_command(
            calibration_args, calibration, temperature, args.calibration_dU_m
        )
        if not (calibration / COMPLETION_MANIFEST).exists():
            _run(command, env, calibration / "run.log")

        crow = _read_rows(_find_steps(calibration, temperature))[0]
        sigma_trial = abs(float(crow["Ftop_N"])) / args.specimen_width_m
        if not math.isfinite(sigma_trial) or sigma_trial <= 0.0:
            raise RuntimeError(f"invalid calibration nominal stress {sigma_trial}")
        calibration_rows.append(
            {
                "temperature_K": temperature,
                "trial_dU_m": args.calibration_dU_m,
                "trial_sigma_max_MPa": sigma_trial / 1.0e6,
                "trial_KJmax_MPa_sqrt_m": float(crow["KJ_Pa_sqrtm"])
                / 1.0e6,
                "specimen_width_m": args.specimen_width_m,
            }
        )

        for delta_sigma in args.delta_sigma_MPa:
            sigma_max_target = delta_sigma * 1.0e6 / (1.0 - args.R)
            dU = args.calibration_dU_m * sigma_max_target / sigma_trial
            outdir = (
                root
                / f"T{int(round(temperature)):04d}K"
                / f"DeltaSigma_{_base._float_token(delta_sigma)}MPa"
            )
            completion = outdir / COMPLETION_MANIFEST
            if outdir.exists() and not args.keep_existing:
                shutil.rmtree(outdir)
            if not completion.exists():
                _run(
                    _base_command(args, outdir, temperature, dU),
                    env,
                    outdir / "run.log",
                )

            case, local, corrected = _summarize_case_v10054(
                outdir, temperature, delta_sigma, args
            )
            case["calibrated_dU_m"] = dU
            case_rows.append(case)
            local_rows.extend(local)
            corrected_block_rows.extend(corrected)

    _write_csv(root / "remote_stress_calibration.csv", calibration_rows)
    _write_csv(root / "K_vs_delta_sigma.csv", case_rows)
    _write_csv(root / "fatigue_growth_points.csv", local_rows)
    _write_csv(
        root / "fatigue_block_diagnostics_v10_0_5_4.csv",
        corrected_block_rows,
    )
    _make_plots(root, case_rows)

    censored = [row for row in case_rows if bool(row.get("right_censored"))]
    manifest = {
        "schema": "v10_0_5_4_vhcf_delta_sigma_campaign",
        "point_release": POINT_RELEASE,
        "material_class": args.material_class,
        "R": args.R,
        "frequency_Hz": args.frequency_Hz,
        "delta_sigma_MPa": args.delta_sigma_MPa,
        "temperatures_K": args.temperatures,
        "cycles_max": args.cycles_max,
        "max_block_cycles": args.max_block_cycles,
        "max_block_cycles_is_unbounded": not math.isfinite(
            args.max_block_cycles
        ),
        "cyclic_mechanics_resolved": bool(args.resolve_cyclic_mechanics),
        "right_censored_case_count": len(censored),
        "case_count": len(case_rows),
        "outputs": {
            "K_vs_delta_sigma": "K_vs_delta_sigma.csv",
            "growth_points": "fatigue_growth_points.csv",
            "block_diagnostics": "fatigue_block_diagnostics_v10_0_5_4.csv",
            "calibration": "remote_stress_calibration.csv",
        },
        "constitutive_physics_changed": False,
        "loading_change": (
            "cycle-integrated remote cyclic stress with hazard-limited "
            "first-passage jumps"
        ),
    }
    (root / "campaign_manifest.json").write_text(
        json.dumps(manifest, indent=2, default=str)
    )
    print(root / "K_vs_delta_sigma.csv", flush=True)

    if censored:
        print(
            f"WARNING: {len(censored)} case(s) are right-censored by MAX_BLOCKS",
            flush=True,
        )
        if args.fail_on_censor:
            return 3
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
