#!/usr/bin/env python3
"""Screen selected MPZ classes at matched own-initiation stress before FEM use.

The audit reports absolute cleavage and emission rates, the self-consistent
one-renewal opening time, source-exhaustion time, retention lifetime, and
shielding persistence.  It is diagnostic only; it does not refit parameters.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
from scipy.special import gammainc

from arrhenius_fracture.config import KB, EV_TO_J
from arrhenius_fracture.moving_process_zone import MovingProcessZoneState
from arrhenius_fracture.mpz_parameterization_v911 import (
    TREF_K,
    build_mpz_config,
    load_selected_row,
    normalize_class_name,
)

DEFAULT_KINIT = {"ceramic": 11.820868, "weakT": 16.949365, "DBTT": 29.197177}


def _exp_floor_barrier_eV(row: dict, prefix: str, sigma_Pa: float, T_K: float) -> float:
    dT = float(T_K) - float(TREF_K)
    G0 = max(
        float(row[f"{prefix}_G00_eV"])
        + float(row[f"{prefix}_gT_eV_per_K"]) * dT,
        1.0e-12,
    )
    sigc = max(
        (
            float(row[f"{prefix}_sigc0_GPa"])
            + float(row[f"{prefix}_sT_GPa_per_K"]) * dT
        ) * 1.0e9,
        1.0,
    )
    floor = min(
        0.95 * G0,
        max(1.0e-4, float(row[f"{prefix}_floor_frac"]) * G0),
    )
    x = max(float(sigma_Pa), 0.0) / sigc
    return float(
        floor
        + (G0 - floor)
        * math.exp(
            -max(float(row[f"{prefix}_exp_a"]), 0.0)
            * x ** max(float(row[f"{prefix}_exp_n"]), 1.0e-9)
        )
    )


def _rate(nu0_s: float, barrier_eV: float, T_K: float) -> float:
    exponent = -float(barrier_eV) * EV_TO_J / max(KB * float(T_K), 1.0e-30)
    return float(max(float(nu0_s), 0.0) * math.exp(float(np.clip(exponent, -700.0, 0.0))))


def _multihit_rate(raw_rate_s: float, m: float = 3.0, tau_s: float = 1.0e-6) -> float:
    return float(gammainc(max(float(m), 1.0), min(max(raw_rate_s, 0.0) * tau_s, 1.0e12)) / tau_s)


def _time_grid(t_end: float) -> np.ndarray:
    early = np.geomspace(1.0e-12, 1.0, 120)
    middle = np.linspace(1.0, min(100.0, t_end), 120)
    late = np.linspace(100.0, t_end, 391) if t_end > 100.0 else np.array([t_end])
    return np.unique(np.concatenate(([0.0], early, middle, late, [t_end])))


def _simulate_hold(state: MovingProcessZoneState, T_K: float, sigma_Pa: float,
                   emission_rate_per_site_s: float, b_m: float,
                   checkpoints=(100.0, 4000.0)) -> dict[str, float]:
    end = max(float(x) for x in checkpoints)
    outputs: dict[float, dict[str, float]] = {}
    previous = 0.0
    for current in _time_grid(end)[1:]:
        dt = float(current - previous)
        state.evolve(
            dt,
            T_K,
            sigma_Pa,
            b_m,
            emission_hazard_integral=emission_rate_per_site_s * dt,
        )
        previous = float(current)
        for checkpoint in checkpoints:
            if checkpoint not in outputs and current >= checkpoint - 1.0e-12:
                outputs[checkpoint] = {
                    "emitted": float(state.emitted_total),
                    "retained": float(state.retained_count),
                    "mobile": float(state.mobile_count),
                    "Kshield": float(state.shielding_K(160.0e9, 0.28, b_m)),
                    "available_fraction": float(state.available_site_fraction),
                }
    result: dict[str, float] = {}
    for checkpoint in checkpoints:
        snap = outputs.get(checkpoint, {})
        tag = f"{int(round(checkpoint))}s"
        for key, value in snap.items():
            result[f"{key}_{tag}"] = float(value)
    return result


def parse_kinit(text: str) -> dict[str, float]:
    values = dict(DEFAULT_KINIT)
    for token in str(text).replace(",", " ").split():
        if "=" not in token:
            continue
        name, raw = token.split("=", 1)
        values[normalize_class_name(name)] = float(raw)
    return values


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--parameter-root", type=Path, default=Path("mpz_v9_11_parameters"))
    p.add_argument("--T-K", type=float, default=700.0)
    p.add_argument("--classes", default="ceramic weakT DBTT")
    p.add_argument(
        "--K-init",
        default="ceramic=11.820868 weakT=16.949365 DBTT=29.197177",
        help="own-initiation K values in MPa sqrt(m)",
    )
    p.add_argument("--r0-m", type=float, default=1.0e-6)
    p.add_argument("--mpz-length-um", type=float, default=100.0)
    p.add_argument("--mpz-n-bins", type=int, default=200)
    p.add_argument("--out", type=Path, default=Path("matched_stress_classes_v917"))
    args = p.parse_args()

    kinit = parse_kinit(args.K_init)
    classes = [normalize_class_name(x) for x in args.classes.replace(",", " ").split()]
    rows = []
    b_m = 2.74e-10
    for cls in classes:
        manifest = args.parameter_root / cls / "spatial_promotion_manifest.csv"
        row = load_selected_row(manifest, cls)
        K_MPa = float(kinit[cls])
        sigma = K_MPa * 1.0e6 / math.sqrt(2.0 * math.pi * max(args.r0_m, 1.0e-30))
        Gc = _exp_floor_barrier_eV(row, "cleave", sigma, args.T_K)
        Ge = _exp_floor_barrier_eV(row, "emit", sigma, args.T_K)
        lambda_c_raw = _rate(1.0e12, Gc, args.T_K)
        lambda_c = _multihit_rate(lambda_c_raw)
        lambda_e = _rate(1.0e11, Ge, args.T_K)

        cfg_args = SimpleNamespace(
            mpz_length_um=float(args.mpz_length_um),
            mpz_n_bins=int(args.mpz_n_bins),
            r_pz=1.0e-6,
        )
        cfg = build_mpz_config(cfg_args, row)
        state = MovingProcessZoneState(cfg)
        total_sites = float(np.sum(state.site_capacity))
        hold = _simulate_hold(state, args.T_K, sigma, lambda_e, b_m)
        t_open = 1.0 / max(lambda_c, 1.0e-300)
        t_emit_one = 1.0 / max(lambda_e * total_sites, 1.0e-300)
        recovery_rate = max(float(row["retained_recovery_rate_s"]), 0.0)
        retention_lifetime = 1.0 / recovery_rate if recovery_rate > 0.0 else math.inf
        rows.append({
            "class": cls,
            "T_K": float(args.T_K),
            "K_init_MPa_sqrt_m": K_MPa,
            "sigma_tip_Pa": sigma,
            "G_cleave_eV": Gc,
            "G_emit_eV": Ge,
            "lambda_c_raw_s-1": lambda_c_raw,
            "lambda_c_multihit_s-1": lambda_c,
            "self_consistent_opening_time_s": t_open,
            "lambda_e_per_site_s-1": lambda_e,
            "total_source_sites": total_sites,
            "time_to_one_emission_full_inventory_s": t_emit_one,
            "Pi_open_over_emit": t_open / t_emit_one,
            "retained_recovery_rate_s-1": recovery_rate,
            "retained_lifetime_s": retention_lifetime,
            "source_refresh_length_um": float(row["source_refresh_length_um"]),
            "source_refresh_fraction_per_5um": min(
                5.0 / max(float(row["source_refresh_length_um"]), 1.0e-30), 1.0
            ),
            "manifest_max_K_shield_MPa_sqrt_m": row.get("max_K_shield_MPa_sqrt_m"),
            **hold,
        })

    frame = pd.DataFrame(rows)
    args.out.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.out / "matched_stress_classes_v917.csv", index=False)
    (args.out / "matched_stress_classes_v917.json").write_text(
        json.dumps(rows, indent=2, default=str)
    )
    summary = {
        "schema": "matched_stress_class_screen_v917_v1",
        "T_K": float(args.T_K),
        "classes": classes,
        "opening_clock": "one absolute cleavage renewal hazard unit",
        "purpose": (
            "screen absolute cleavage/emission timescale separation and retention "
            "before FEM propagation"
        ),
        "rows": rows,
    }
    (args.out / "matched_stress_classes_v917_summary.json").write_text(
        json.dumps(summary, indent=2, default=str)
    )
    print(frame.to_string(index=False))
    print("wrote", args.out / "matched_stress_classes_v917.csv")


if __name__ == "__main__":
    main()
