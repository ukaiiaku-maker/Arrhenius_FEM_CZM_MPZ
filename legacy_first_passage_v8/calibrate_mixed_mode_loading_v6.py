#!/usr/bin/env python3
"""Exact adaptive-CZM production-backend calibration for mixed-mode v6.

Unlike v5 and earlier calibrators, this script does not assemble a separate
static damaged-notch model.  Every basis and verification point is a one-step
run through the same sharp-front/adaptive-CZM driver used by the physical
campaign.  The measured first production-step tractions therefore define the
loading response matrix.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import shlex
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from arrhenius_fracture.mixed_mode_first_passage_v6 import (
    MODEL_ID,
    angle_error_deg,
    loading_coefficients_from_response_basis,
    loading_z_from_coefficients,
    shear_sign_from_basis,
    traction_phase_deg,
)

CALIBRATION_ID = "mixed_mode_fem_czm_v6_exact_production_backend_calibration"


def vals(s):
    return [float(x) for x in str(s).replace(",", " ").split() if x]


def fs(x):
    return f"{float(x):.16g}"


def bool_col(s):
    return str(s).strip().lower() in {"1", "true", "yes"}


def probe_command(py, a, qo, qs, target, out, ref_c=1.0, ref_s=0.0, shear_sign=1.0):
    return [
        py, "-m", "arrhenius_fracture.mixed_mode_first_passage_v6",
        "--mixity-open-coeff", fs(qo),
        "--mixity-shear-coeff", fs(qs),
        "--target-traction-phase-deg", fs(target),
        "--traction-shear-sign", fs(shear_sign),
        "--traction-probe-radius-m", fs(a.traction_probe_radius_m),
        "--reference-cleavage-shape", fs(ref_c),
        "--reference-slip-shape", fs(ref_s),
        "--shear-emission-weight", fs(a.shear_emission_weight),
        "--directional-factor-max", fs(a.directional_factor_max),
        "--mode", "2d",
        "--nx", str(a.nx), "--ny", str(a.ny),
        "--tip-h-fine", fs(a.tip_h_fine), "--tip-ratio", fs(a.tip_ratio),
        "--dU", fs(a.U_cal_m), "--dt", fs(a.dt), "--steps", "1",
        "--n-stagger", "2", "--print-every", "1",
        "--stop-after-first-fire", "--max-fronts", "1",
        "--adaptive-events", "--adaptive-event-target", ".25",
        "--adaptive-min-frac", "1e-8", "--adaptive-grow", "4",
        "--da-phys", "5e-6",
        "--j-decomposition", "cluster", "--rJ-cluster", "20e-6",
        "--rJ-outer", "25e-6", "--temperatures", fs(a.T_K),
        "--crack-backend", "adaptive_czm", "--czm-max-angle-error-deg", "35",
        "--crystal-aniso", "--crystal-compete",
        "--crystal-theta-deg", fs(a.crystal_theta_deg),
        "--crystal-C11", fs(a.crystal_C11),
        "--crystal-C12", fs(a.crystal_C12),
        "--crystal-C44", fs(a.crystal_C44),
        "--cleave-gamma-aniso", fs(a.cleave_gamma_aniso),
        "--crystal-material", "w",
        # High barriers make the one-step calibration mechanically passive.
        "--emit-barrier-kind", "exp_floor",
        "--emit-G00-eV", "20", "--emit-gT-eV-per-K", "0",
        "--emit-sigc0-GPa", "5", "--emit-sT-GPa-per-K", "0",
        "--emit-exp-a", "1", "--emit-exp-n", "1",
        "--emit-floor-frac", "0.01", "--emit-Tref-K", "300",
        "--cleave-barrier-kind", "exp_floor", "--cleave-exp-T-mode", "linear",
        "--cleave-G00-eV", "20", "--cleave-gT-eV-per-K", "0",
        "--cleave-sigc0-GPa", "5", "--cleave-sT-GPa-per-K", "0",
        "--cleave-exp-a", "1", "--cleave-exp-n", "1",
        "--cleave-floor-frac", "0.01", "--cleave-S-hs-kB", "0",
        "--cleave-sigma-S-GPa", "6", "--cleave-S-hs-power", "2",
        "--cleave-S-hs-Tref-K", "300", "--cleave-Tref-K", "300",
        "--cleave-shield-chi", "0", "--n-sat", "inf",
        "--multihit-m", "3", "--multihit-tau", "1e-6",
        "--emb-sat-frac", "1", "--save-snapshots", "0", "--no-plots",
        "--out", str(out),
    ]


def read_probe(out: Path) -> dict:
    p = out / "anisotropic_calibrated_tip_calls.csv"
    if not p.exists():
        raise RuntimeError(f"missing production-backend probe output: {p}")
    df = pd.read_csv(p)
    if "traction_probe_reliable" in df:
        mask = df["traction_probe_reliable"].map(bool_col)
        if mask.any():
            df = df[mask]
    if df.empty:
        raise RuntimeError(f"no reliable production-backend probe rows in {p}")
    # Multiple stagger/J evaluations can occur at the same first load.  Use a
    # robust median for numeric quantities and the last categorical value.
    row = {}
    for col in df.columns:
        num = pd.to_numeric(df[col], errors="coerce")
        if num.notna().any():
            row[col] = float(num.median())
        else:
            row[col] = df[col].iloc[-1]
    row["n_probe_records"] = int(len(df))
    return row


def run_probe(py, a, qo, qs, target, out, ref_c=1.0, ref_s=0.0, shear_sign=1.0):
    out.mkdir(parents=True, exist_ok=True)
    cmd = probe_command(py, a, qo, qs, target, out, ref_c, ref_s, shear_sign)
    (out / "command.txt").write_text(shlex.join(cmd) + "\n")
    with (out / "run.log").open("w") as log:
        rc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT).returncode
    if rc:
        tail = "\n".join((out / "run.log").read_text(errors="replace").splitlines()[-30:])
        raise RuntimeError(f"production-backend probe failed rc={rc}:\n{tail}")
    return read_probe(out)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--out", default="runs/mixed_mode_fem_czm_v6_backend_calibration")
    p.add_argument("--target-psi-deg", default="-60 -45 -30 -15 0 15 30 45 60")
    p.add_argument("--U-cal-m", type=float, default=2e-7)
    p.add_argument("--dt", type=float, default=8.4)
    p.add_argument("--T-K", type=float, default=500.0)
    p.add_argument("--nx", type=int, default=24)
    p.add_argument("--ny", type=int, default=48)
    p.add_argument("--tip-h-fine", type=float, default=3e-6)
    p.add_argument("--tip-ratio", type=float, default=1.25)
    p.add_argument("--crystal-theta-deg", type=float, default=45.0)
    p.add_argument("--crystal-C11", type=float, default=523e9)
    p.add_argument("--crystal-C12", type=float, default=203e9)
    p.add_argument("--crystal-C44", type=float, default=160e9)
    p.add_argument("--cleave-gamma-aniso", type=float, default=0.3)
    p.add_argument("--traction-probe-radius-m", type=float, default=10e-6)
    p.add_argument("--psi-tol-deg", type=float, default=0.75)
    p.add_argument("--basis-condition-max", type=float, default=1e8)
    p.add_argument("--shear-emission-weight", type=float, default=1.0)
    p.add_argument("--directional-factor-max", type=float, default=5.0)
    p.add_argument("--force", action="store_true")
    a = p.parse_args()

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    probes = out / "production_backend_probes"
    py = sys.executable

    opening = run_probe(py, a, 1.0, 0.0, 0.0, probes / "basis_open")
    sliding = run_probe(py, a, 0.0, 1.0, 0.0, probes / "basis_slide")
    Mraw = np.array([
        [opening["reference_sigma_nn_Pa"], sliding["reference_sigma_nn_Pa"]],
        [opening["reference_tau_tn_Pa"], sliding["reference_tau_tn_Pa"]],
    ], float)
    sign = shear_sign_from_basis(Mraw)
    M = np.diag([1.0, sign]) @ Mraw
    cond = float(np.linalg.cond(M))
    if not np.isfinite(cond) or cond > a.basis_condition_max:
        raise SystemExit(f"production-backend traction basis invalid cond={cond}")

    q0 = loading_coefficients_from_response_basis(M, 0.0)
    ref = run_probe(py, a, *q0, 0.0, probes / "reference_mode_I", shear_sign=sign)
    ref_c = float(ref["cleavage_shape"])
    ref_s = max(float(ref["slip_shape"]), 0.0)
    if not np.isfinite(ref_c) or ref_c <= 1e-12:
        raise SystemExit(f"invalid production-backend Mode-I cleavage shape {ref_c}")

    rows = []
    for target in vals(a.target_psi_deg):
        qo, qs = loading_coefficients_from_response_basis(M, target)
        label = ("m" if target < 0 else "p") + f"{abs(target):05.1f}".replace(".", "p")
        r = run_probe(py, a, qo, qs, target, probes / f"verify_{label}", ref_c, ref_s, sign)
        phase = traction_phase_deg(r["reference_sigma_nn_Pa"], r["reference_tau_tn_Pa"], sign)
        error = angle_error_deg(phase, target)
        reliable = bool(bool_col(r.get("traction_probe_reliable", True)))
        row = {
            **r,
            "calibration_id": CALIBRATION_ID,
            "model": MODEL_ID,
            "target_psi_deg": target,
            "loading_open_coeff": qo,
            "loading_shear_coeff": qs,
            "loading_z": loading_z_from_coefficients(qo, qs),
            "loading_angle_deg": math.degrees(math.atan2(qs, qo)),
            "traction_shear_sign": sign,
            "achieved_traction_phase_deg": phase,
            "traction_phase_error_deg": error,
            "phase_converged": bool(reliable and abs(error) <= a.psi_tol_deg),
            "first_production_step_verified": bool(reliable and abs(error) <= a.psi_tol_deg),
            "reference_cleavage_shape": ref_c,
            "reference_slip_shape": ref_s,
            "basis_condition": cond,
            "response_11_Pa": M[0, 0], "response_12_Pa": M[0, 1],
            "response_21_Pa": M[1, 0], "response_22_Pa": M[1, 1],
            "crack_backend": "adaptive_czm",
            "crystal_theta_deg": a.crystal_theta_deg,
            "traction_probe_radius_m": a.traction_probe_radius_m,
        }
        rows.append(row)
        print({k: row[k] for k in (
            "target_psi_deg", "loading_open_coeff", "loading_shear_coeff",
            "loading_z", "achieved_traction_phase_deg", "traction_phase_error_deg",
            "phase_converged")})

    csv_path = out / "mixed_mode_loading_calibration_v6.csv"
    with csv_path.open("w", newline="") as fp:
        cols = sorted({k for r in rows for k in r})
        w = csv.DictWriter(fp, fieldnames=cols)
        w.writeheader(); w.writerows(rows)
    (out / "production_backend_basis_v6.json").write_text(json.dumps({
        "calibration_id": CALIBRATION_ID,
        "model": MODEL_ID,
        "raw_response_matrix_Pa": Mraw.tolist(),
        "normalized_response_matrix_Pa": M.tolist(),
        "traction_shear_sign": sign,
        "basis_condition": cond,
        "reference_loading_coefficients": list(q0),
        "reference_cleavage_shape": ref_c,
        "reference_slip_shape": ref_s,
        "backend": "adaptive_czm",
        "mesh": {"nx": a.nx, "ny": a.ny, "tip_h_fine": a.tip_h_fine, "tip_ratio": a.tip_ratio},
        "crystal": {"theta_deg": a.crystal_theta_deg, "C11": a.crystal_C11,
                    "C12": a.crystal_C12, "C44": a.crystal_C44},
        "probe_radius_m": a.traction_probe_radius_m,
    }, indent=2))
    print("raw production-backend traction basis [Pa]:\n", Mraw)
    print("normalized basis [Pa]:\n", M)
    print("basis condition:", cond)
    print("Mode-I reference shapes:", {"cleavage": ref_c, "slip": ref_s})
    bad = [r for r in rows if not r["phase_converged"]]
    if bad:
        raise SystemExit("v6 production-backend calibration failed: " +
                         ", ".join(f"{r['target_psi_deg']}:err={r['traction_phase_error_deg']:.4g}" for r in bad))
    print("wrote", csv_path)


if __name__ == "__main__":
    main()
