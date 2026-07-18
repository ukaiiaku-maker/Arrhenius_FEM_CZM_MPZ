#!/usr/bin/env python3
"""Calibrate and run v10.0.5.3 fatigue as a remote stress-range sweep.

For each requested Delta-sigma, the script converts to sigma_max=Delta-sigma/(1-R),
calibrates the imposed displacement from the model's own top reaction, runs the
progressive kinetic-CZM fatigue entry point, and writes K-versus-Delta-sigma plus
local/integrated crack-growth summaries.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Iterable


def _float_token(value: float) -> str:
    return f"{value:g}".replace("-", "m").replace(".", "p")


def _read_rows(path: Path) -> list[dict[str, float]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        rows = []
        for row in reader:
            converted = {}
            for key, value in row.items():
                try:
                    converted[key] = float(value)
                except (TypeError, ValueError):
                    converted[key] = value
            rows.append(converted)
        return rows


def _find_steps(outdir: Path, temperature: float) -> Path:
    exact = outdir / f"steps_{int(round(temperature)):04d}K.csv"
    if exact.exists():
        return exact
    matches = sorted(outdir.glob("steps_*K.csv"))
    if len(matches) != 1:
        raise RuntimeError(f"expected one steps CSV in {outdir}; found {matches}")
    return matches[0]


def _run(command: list[str], env: dict[str, str], log_path: Path) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w") as log:
        process = subprocess.run(
            command,
            env=env,
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
    if process.returncode != 0:
        raise RuntimeError(
            f"command failed with exit code {process.returncode}; see {log_path}"
        )


def _base_command(args, outdir: Path, temperature: float, dU_m: float) -> list[str]:
    return [
        args.python,
        "-m",
        "arrhenius_fracture.mode_i_first_passage_v10_0_5_3_fatigue",
        "--v10-material-class",
        args.material_class,
        "--czm-opening-coupling",
        "clock_linear",
        "--mode",
        "2d",
        "--temperatures",
        f"{temperature:g}",
        "--out",
        str(outdir),
        "--steps",
        str(args.max_blocks),
        "--dU",
        f"{dU_m:.17g}",
        "--dt",
        f"{args.outer_dt_s:.17g}",
        "--fatigue-cycles",
        "--fatigue-hold-load",
        "--R",
        f"{args.R:.17g}",
        "--frequency-Hz",
        f"{args.frequency_Hz:.17g}",
        "--cycles-max",
        f"{args.cycles_max:.17g}",
        "--block-cycles",
        f"{args.block_cycles:.17g}",
        "--max-block-cycles",
        f"{args.max_block_cycles:.17g}",
        "--min-block-cycles",
        f"{args.min_block_cycles:.17g}",
        "--cycle-block-mode",
        "hazard_limited",
        "--target-dB",
        f"{args.target_dB:.17g}",
        "--target-dN-store",
        f"{args.target_dN_store:.17g}",
        "--target-dN-emit",
        f"{args.target_dN_emit:.17g}",
        "--target-dN-mobile",
        f"{args.target_dN_mobile:.17g}",
        "--n-phase",
        str(args.n_phase),
        "--cyclic-mechanics-phases",
        str(args.cyclic_mechanics_phases),
        "--nx",
        str(args.nx),
        "--ny",
        str(args.ny),
        "--tip-h-fine",
        f"{args.tip_h_fine_m:.17g}",
        "--tip-ratio",
        f"{args.tip_ratio:.17g}",
        "--da-phys",
        f"{args.da_phys_m:.17g}",
        "--rJ-outer",
        f"{args.rJ_outer_m:.17g}",
        "--L-pz",
        f"{args.L_pz_m:.17g}",
        "--mpz-length-um",
        f"{args.mpz_length_um:.17g}",
        "--mpz-n-bins",
        str(args.mpz_n_bins),
        "--crack-backend",
        "adaptive_czm",
        "--crystal-aniso",
        "--crystal-theta-deg",
        f"{args.theta_deg:.17g}",
        "--max-fronts",
        "1",
        "--adaptive-events",
        "--adaptive-event-target",
        f"{args.adaptive_event_target:.17g}",
        "--adaptive-min-frac",
        f"{args.adaptive_min_frac:.17g}",
        "--target-crack-extension-um",
        f"{args.target_extension_um:.17g}",
        "--target-da-per-block-um",
        f"{args.target_da_per_block_um:.17g}",
        "--max-da-per-block-um",
        f"{args.max_da_per_block_um:.17g}",
        "--save-snapshots",
        str(args.save_snapshots),
        "--snapshot-cols",
        str(args.snapshot_cols),
        "--snapshot-by-crack-extension-um",
        f"{args.snapshot_by_extension_um:.17g}",
        "--print-every",
        str(args.print_every),
    ]


def _write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    keys = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def _summarize_case(
    outdir: Path,
    temperature: float,
    delta_sigma_MPa: float,
    R: float,
    specimen_width_m: float,
) -> tuple[dict, list[dict]]:
    rows = _read_rows(_find_steps(outdir, temperature))
    if not rows:
        raise RuntimeError(f"empty steps CSV in {outdir}")
    first = rows[0]
    last = rows[-1]
    sigma_max_actual_Pa = abs(float(first["Ftop_N"])) / specimen_width_m
    Kmax_first = float(first["KJ_Pa_sqrtm"])
    Kmax_last = float(last["KJ_Pa_sqrtm"])
    cycles_total = sum(max(float(row.get("fatigue_cycles", 0.0)), 0.0) for row in rows)
    crack_extension = max(float(last.get("crack_extension_m", 0.0)), 0.0)
    integrated_da_dN = crack_extension / cycles_total if cycles_total > 0.0 else 0.0
    case = {
        "temperature_K": temperature,
        "delta_sigma_requested_MPa": delta_sigma_MPa,
        "R": R,
        "sigma_max_requested_MPa": delta_sigma_MPa / max(1.0 - R, 1.0e-300),
        "sigma_max_actual_MPa_first": sigma_max_actual_Pa / 1.0e6,
        "delta_sigma_actual_MPa_first": sigma_max_actual_Pa * (1.0 - R) / 1.0e6,
        "KJmax_first_MPa_sqrt_m": Kmax_first / 1.0e6,
        "KJmin_first_MPa_sqrt_m": max(R * Kmax_first, 0.0) / 1.0e6,
        "DeltaKJ_first_MPa_sqrt_m": Kmax_first * (1.0 - R) / 1.0e6,
        "KJmax_last_MPa_sqrt_m": Kmax_last / 1.0e6,
        "DeltaKJ_last_MPa_sqrt_m": Kmax_last * (1.0 - R) / 1.0e6,
        "cycles_total": cycles_total,
        "crack_extension_um": crack_extension * 1.0e6,
        "integrated_da_dN_m_per_cycle": integrated_da_dN,
        "n_blocks": len(rows),
        "run_directory": str(outdir),
    }
    local = []
    cycles_cumulative = 0.0
    for row in rows:
        cycles = max(float(row.get("fatigue_cycles", 0.0)), 0.0)
        cycles_cumulative += cycles
        da = max(float(row.get("da_block_m", 0.0)), 0.0)
        K = float(row.get("KJ_Pa_sqrtm", 0.0))
        local.append(
            {
                "temperature_K": temperature,
                "delta_sigma_requested_MPa": delta_sigma_MPa,
                "cycles_block": cycles,
                "cycles_cumulative": cycles_cumulative,
                "crack_extension_um": float(row.get("crack_extension_m", 0.0)) * 1.0e6,
                "da_block_m": da,
                "local_da_dN_m_per_cycle": da / cycles if cycles > 0.0 else 0.0,
                "KJmax_MPa_sqrt_m": K / 1.0e6,
                "DeltaKJ_MPa_sqrt_m": K * (1.0 - R) / 1.0e6,
                "B": float(row.get("B", 0.0)),
                "dB_block": float(row.get("dB_block", 0.0)),
                "n_fire": int(float(row.get("n_fire", 0.0))),
            }
        )
    return case, local


def _make_plots(root: Path, case_rows: list[dict]) -> None:
    if not case_rows:
        return
    try:
        import matplotlib.pyplot as plt
    except Exception as exc:
        print(f"WARNING: plotting unavailable: {exc}")
        return
    temperatures = sorted({float(row["temperature_K"]) for row in case_rows})

    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    for temperature in temperatures:
        rows = sorted(
            (row for row in case_rows if float(row["temperature_K"]) == temperature),
            key=lambda row: float(row["delta_sigma_requested_MPa"]),
        )
        x = [float(row["delta_sigma_requested_MPa"]) for row in rows]
        y = [float(row["KJmax_first_MPa_sqrt_m"]) for row in rows]
        ax.plot(x, y, marker="o", label=f"Kmax, {temperature:g} K")
        yd = [float(row["DeltaKJ_first_MPa_sqrt_m"]) for row in rows]
        ax.plot(x, yd, marker="s", linestyle="--", label=f"DeltaK, {temperature:g} K")
    ax.set_xlabel(r"Remote stress range, $\Delta\sigma$ (MPa)")
    ax.set_ylabel(r"Stress-intensity range / maximum (MPa $\sqrt{m}$)")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(root / "K_vs_delta_sigma.png", dpi=300)
    fig.savefig(root / "K_vs_delta_sigma.pdf")
    plt.close(fig)

    positive = [row for row in case_rows if float(row["integrated_da_dN_m_per_cycle"]) > 0.0]
    if positive:
        fig, ax = plt.subplots(figsize=(7.0, 5.0))
        for temperature in temperatures:
            rows = sorted(
                (row for row in positive if float(row["temperature_K"]) == temperature),
                key=lambda row: float(row["DeltaKJ_first_MPa_sqrt_m"]),
            )
            if not rows:
                continue
            ax.plot(
                [float(row["DeltaKJ_first_MPa_sqrt_m"]) for row in rows],
                [float(row["integrated_da_dN_m_per_cycle"]) for row in rows],
                marker="o",
                label=f"{temperature:g} K",
            )
        ax.set_yscale("log")
        ax.set_xlabel(r"$\Delta K_J$ (MPa $\sqrt{m}$)")
        ax.set_ylabel(r"Integrated $da/dN$ (m cycle$^{-1}$)")
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(root / "da_dN_vs_delta_K.png", dpi=300)
        fig.savefig(root / "da_dN_vs_delta_K.pdf")
        plt.close(fig)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", required=True)
    p.add_argument("--delta-sigma-MPa", type=float, nargs="+", required=True)
    p.add_argument("--temperatures", type=float, nargs="+", default=[300.0])
    p.add_argument("--material-class", choices=["DBTT", "weakT", "ceramic"], default="DBTT")
    p.add_argument("--R", type=float, default=0.1)
    p.add_argument("--frequency-Hz", type=float, default=1.0e3, dest="frequency_Hz")
    p.add_argument("--cycles-max", type=float, default=1.0e12, dest="cycles_max")
    p.add_argument("--block-cycles", type=float, default=1.0e4, dest="block_cycles")
    p.add_argument("--max-block-cycles", type=float, default=1.0e8, dest="max_block_cycles")
    p.add_argument("--min-block-cycles", type=float, default=1.0e-6, dest="min_block_cycles")
    p.add_argument("--max-blocks", type=int, default=20000)
    p.add_argument("--target-dB", type=float, default=0.01, dest="target_dB")
    p.add_argument("--target-dN-store", type=float, default=0.01, dest="target_dN_store")
    p.add_argument("--target-dN-emit", type=float, default=0.10, dest="target_dN_emit")
    p.add_argument("--target-dN-mobile", type=float, default=0.10, dest="target_dN_mobile")
    p.add_argument("--n-phase", type=int, default=96, dest="n_phase")
    p.add_argument("--cyclic-mechanics-phases", type=int, default=16, dest="cyclic_mechanics_phases")
    p.add_argument("--theta-deg", type=float, default=45.0)
    p.add_argument("--nx", type=int, default=60)
    p.add_argument("--ny", type=int, default=120)
    p.add_argument("--tip-h-fine-m", type=float, default=2.5e-6)
    p.add_argument("--tip-ratio", type=float, default=1.2)
    p.add_argument("--da-phys-m", type=float, default=5.0e-6)
    p.add_argument("--rJ-outer-m", type=float, default=60.0e-6)
    p.add_argument("--L-pz-m", type=float, default=100.0e-6)
    p.add_argument("--mpz-length-um", type=float, default=100.0)
    p.add_argument("--mpz-n-bins", type=int, default=200)
    p.add_argument("--target-extension-um", type=float, default=250.0)
    p.add_argument("--target-da-per-block-um", type=float, default=5.0)
    p.add_argument("--max-da-per-block-um", type=float, default=10.0)
    p.add_argument("--save-snapshots", type=int, default=0)
    p.add_argument("--snapshot-cols", type=int, default=6)
    p.add_argument("--snapshot-by-extension-um", type=float, default=25.0)
    p.add_argument("--print-every", type=int, default=10)
    p.add_argument("--specimen-width-m", type=float, default=2.0e-3)
    p.add_argument("--calibration-dU-m", type=float, default=1.0e-7)
    p.add_argument("--outer-dt-s", type=float, default=1.0e-9)
    p.add_argument("--adaptive-event-target", type=float, default=0.05)
    p.add_argument("--adaptive-min-frac", type=float, default=1.0e-10)
    p.add_argument("--python", default=sys.executable)
    p.add_argument("--keep-existing", action="store_true")
    return p


def main(argv: Iterable[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.R >= 1.0:
        raise SystemExit("R must be less than 1 for a positive stress range")
    root = Path(args.out).resolve()
    root.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env.update(
        {
            "ARRHENIUS_CZM_OPENING_COUPLING": "clock_linear",
            "ARRHENIUS_MAX_TRIAL_DAMAGE_CHANGE": str(min(args.target_dB, 0.05)),
            "ARRHENIUS_COMMITTED_TARGET_EXTENSION_UM": str(args.target_extension_um),
        }
    )
    case_rows = []
    local_rows = []
    calibration_rows = []

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
        command = _base_command(
            calibration_args, calibration, temperature, args.calibration_dU_m
        )
        if not (calibration / "run_completion_v10_0_5_3_fatigue.json").exists():
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
                "trial_KJmax_MPa_sqrt_m": float(crow["KJ_Pa_sqrtm"]) / 1.0e6,
                "specimen_width_m": args.specimen_width_m,
            }
        )

        for delta_sigma in args.delta_sigma_MPa:
            sigma_max_target = delta_sigma * 1.0e6 / (1.0 - args.R)
            dU = args.calibration_dU_m * sigma_max_target / sigma_trial
            outdir = (
                root
                / f"T{int(round(temperature)):04d}K"
                / f"DeltaSigma_{_float_token(delta_sigma)}MPa"
            )
            completion = outdir / "run_completion_v10_0_5_3_fatigue.json"
            if outdir.exists() and not args.keep_existing:
                shutil.rmtree(outdir)
            if not completion.exists():
                _run(_base_command(args, outdir, temperature, dU), env, outdir / "run.log")
            case, local = _summarize_case(
                outdir, temperature, delta_sigma, args.R, args.specimen_width_m
            )
            case["calibrated_dU_m"] = dU
            case_rows.append(case)
            local_rows.extend(local)

    _write_csv(root / "remote_stress_calibration.csv", calibration_rows)
    _write_csv(root / "K_vs_delta_sigma.csv", case_rows)
    _write_csv(root / "fatigue_growth_points.csv", local_rows)
    _make_plots(root, case_rows)
    manifest = {
        "schema": "v10_0_5_3_delta_sigma_fatigue_campaign",
        "material_class": args.material_class,
        "R": args.R,
        "frequency_Hz": args.frequency_Hz,
        "delta_sigma_MPa": args.delta_sigma_MPa,
        "temperatures_K": args.temperatures,
        "outputs": {
            "K_vs_delta_sigma": "K_vs_delta_sigma.csv",
            "growth_points": "fatigue_growth_points.csv",
            "calibration": "remote_stress_calibration.csv",
        },
        "constitutive_physics_changed": False,
        "loading_change": "remote cyclic stress range mapped to displacement amplitude",
    }
    (root / "campaign_manifest.json").write_text(json.dumps(manifest, indent=2))
    print(root / "K_vs_delta_sigma.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
