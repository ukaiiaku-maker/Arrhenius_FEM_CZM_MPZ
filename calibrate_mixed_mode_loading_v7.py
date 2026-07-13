#!/usr/bin/env python3
"""Exact adaptive-CZM, full-circle mixed-mode calibration, v7.

The calibration uses the exact production backend for every probe.  The linear
opening/sliding basis supplies only an initial guess.  Each requested phase is
then solved by direct root finding on the measured first-production-step
traction phase.  Boundary loading is represented by the full unit circle
(q_open, q_shear)=(cos alpha, sin alpha); negative q_open is retained and
reported rather than silently flipped.
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

from arrhenius_fracture.mixed_mode_first_passage_v7 import (
    MODEL_ID,
    angle_error_deg,
    loading_alpha_deg_from_coefficients,
    loading_coefficients_from_alpha_deg,
    loading_coefficients_from_response_basis,
    shear_sign_from_basis,
    traction_phase_deg,
    wrap_loading_angle_deg,
)

CALIBRATION_ID = "mixed_mode_fem_czm_v7_exact_backend_full_circle_calibration"


def vals(s):
    return [float(x) for x in str(s).replace(",", " ").split() if x]


def fs(x):
    return f"{float(x):.16g}"


def bool_col(s):
    return str(s).strip().lower() in {"1", "true", "yes"}


def probe_command(py, a, qo, qs, target, out, ref_c=1.0, ref_s=0.0, shear_sign=1.0):
    return [
        py, "-m", "arrhenius_fracture.mixed_mode_first_passage_v7",
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
        # Mechanically passive one-step calibration.
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
    phase_col = "traction_phase_probe_reliable"
    if phase_col in df:
        phase_mask = df[phase_col].map(bool_col)
        if phase_mask.any():
            df = df[phase_mask]
    elif "traction_probe_reliable" in df:
        phase_mask = df["traction_probe_reliable"].map(bool_col)
        if phase_mask.any():
            df = df[phase_mask]
    if df.empty:
        raise RuntimeError(f"no finite phase-probe rows in {p}")
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


def phase_sample_reliable(row: dict) -> bool:
    flag = row.get("traction_phase_probe_reliable", row.get("traction_probe_reliable", False))
    if not bool_col(flag):
        return False
    sn = float(row.get("reference_sigma_nn_Pa", np.nan))
    tt = float(row.get("reference_tau_tn_Pa", np.nan))
    return bool(np.isfinite(sn) and np.isfinite(tt) and math.hypot(sn, tt) > 1e-10)


def closest_unwrapped(alpha_deg: float, reference_deg: float) -> float:
    a = wrap_loading_angle_deg(alpha_deg)
    return float(a + 360.0 * round((float(reference_deg) - a) / 360.0))


def solve_exact_target(py, a, target, alpha0, probes, ref_c, ref_s, sign):
    """Directly solve the first-step phase target using exact backend probes."""
    cache: dict[float, dict] = {}
    serial = 0

    def evaluate(alpha_unwrapped):
        nonlocal serial
        alpha_unwrapped = float(alpha_unwrapped)
        canonical = wrap_loading_angle_deg(alpha_unwrapped)
        key = round(canonical, 10)
        if key in cache:
            old = dict(cache[key])
            old["loading_alpha_unwrapped_deg"] = closest_unwrapped(canonical, alpha_unwrapped)
            return old
        qo, qs = loading_coefficients_from_alpha_deg(canonical)
        label = f"target_{target:+06.1f}_probe_{serial:02d}_a_{canonical:+09.4f}".replace("+", "p").replace("-", "m").replace(".", "p")
        serial += 1
        raw = run_probe(py, a, qo, qs, target, probes / label, ref_c, ref_s, sign)
        phase = traction_phase_deg(raw["reference_sigma_nn_Pa"], raw["reference_tau_tn_Pa"], sign)
        err = angle_error_deg(phase, target)
        rec = {
            **raw,
            "loading_alpha_deg": canonical,
            "loading_alpha_unwrapped_deg": alpha_unwrapped,
            "loading_open_coeff": qo,
            "loading_shear_coeff": qs,
            "loading_open_is_tensile": bool(qo >= 0.0),
            "achieved_traction_phase_deg": phase,
            "traction_phase_error_deg": err,
            "phase_sample_reliable": phase_sample_reliable(raw),
        }
        cache[key] = dict(rec)
        return rec

    samples = [evaluate(alpha0)]
    if samples[-1]["phase_sample_reliable"] and abs(samples[-1]["traction_phase_error_deg"]) <= a.psi_tol_deg:
        return samples[-1], samples

    # Local probes first; cancellation-sensitive targets often converge within
    # a fraction of a degree in boundary-angle space.
    for delta in (0.25, -0.25):
        samples.append(evaluate(alpha0 + delta))
        good = [r for r in samples if r["phase_sample_reliable"] and abs(r["traction_phase_error_deg"]) <= a.psi_tol_deg]
        if good:
            return min(good, key=lambda r: abs(r["traction_phase_error_deg"])), samples

    expansion_offsets = [0.5, -0.5, 1, -1, 2, -2, 4, -4, 8, -8,
                         16, -16, 32, -32, 64, -64, 96, -96, 128, -128, 180]
    next_expand = 0

    for _ in range(max(1, a.max_root_iters)):
        valid = [r for r in samples if r["phase_sample_reliable"] and
                 np.isfinite(r["traction_phase_error_deg"])]
        good = [r for r in valid if abs(r["traction_phase_error_deg"]) <= a.psi_tol_deg]
        if good:
            return min(good, key=lambda r: abs(r["traction_phase_error_deg"])), samples

        # Narrowest sign-changing error bracket.
        bracket = None
        sv = sorted(valid, key=lambda r: r["loading_alpha_unwrapped_deg"])
        for r0, r1 in zip(sv[:-1], sv[1:]):
            e0 = r0["traction_phase_error_deg"]
            e1 = r1["traction_phase_error_deg"]
            if e0 * e1 < 0:
                width = abs(r1["loading_alpha_unwrapped_deg"] - r0["loading_alpha_unwrapped_deg"])
                if bracket is None or width < bracket[0]:
                    bracket = (width, r0, r1)
        if bracket is not None:
            _, r0, r1 = bracket
            a0, a1 = r0["loading_alpha_unwrapped_deg"], r1["loading_alpha_unwrapped_deg"]
            e0, e1 = r0["traction_phase_error_deg"], r1["traction_phase_error_deg"]
            anew = a1 - e1 * (a1-a0) / (e1-e0)
            if not (min(a0, a1) < anew < max(a0, a1)):
                anew = 0.5 * (a0+a1)
            samples.append(evaluate(anew))
            continue

        # Local empirical secant using the two best distinct samples.
        best = sorted(valid, key=lambda r: abs(r["traction_phase_error_deg"]))
        proposed = None
        if len(best) >= 2:
            r1 = best[0]
            r0 = next((r for r in best[1:] if abs(r["loading_alpha_unwrapped_deg"]-r1["loading_alpha_unwrapped_deg"]) > 1e-8), None)
            if r0 is not None:
                da = r1["loading_alpha_unwrapped_deg"] - r0["loading_alpha_unwrapped_deg"]
                de = r1["traction_phase_error_deg"] - r0["traction_phase_error_deg"]
                if abs(de) > 1e-6:
                    step = -r1["traction_phase_error_deg"] * da / de
                    step = float(np.clip(step, -a.max_alpha_step_deg, a.max_alpha_step_deg))
                    proposed = r1["loading_alpha_unwrapped_deg"] + step
        if proposed is not None:
            canonical = round(wrap_loading_angle_deg(proposed), 10)
            if canonical not in cache:
                samples.append(evaluate(proposed))
                continue

        if next_expand < len(expansion_offsets):
            samples.append(evaluate(alpha0 + expansion_offsets[next_expand]))
            next_expand += 1
            continue
        break

    valid = [r for r in samples if r["phase_sample_reliable"] and np.isfinite(r["traction_phase_error_deg"])]
    selected = min(valid, key=lambda r: abs(r["traction_phase_error_deg"])) if valid else samples[0]
    return selected, samples


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--out", default="runs/mixed_mode_fem_czm_v7_backend_calibration")
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
    p.add_argument("--max-root-iters", type=int, default=20)
    p.add_argument("--max-alpha-step-deg", type=float, default=20.0)
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
    all_hist = []
    for target in vals(a.target_psi_deg):
        qguess = loading_coefficients_from_response_basis(M, target)
        alpha0 = loading_alpha_deg_from_coefficients(*qguess)
        selected, history = solve_exact_target(
            py, a, target, alpha0, probes / f"solve_{target:+06.1f}".replace("+", "p").replace("-", "m").replace(".", "p"),
            ref_c, ref_s, sign)
        error = float(selected["traction_phase_error_deg"])
        phase_ok = bool(selected["phase_sample_reliable"] and abs(error) <= a.psi_tol_deg)
        row = {
            **selected,
            "calibration_id": CALIBRATION_ID,
            "model": MODEL_ID,
            "target_psi_deg": target,
            "phase_converged": phase_ok,
            "first_production_step_verified": phase_ok,
            "reference_cleavage_shape": ref_c,
            "reference_slip_shape": ref_s,
            "basis_condition": cond,
            "response_11_Pa": M[0, 0], "response_12_Pa": M[0, 1],
            "response_21_Pa": M[1, 0], "response_22_Pa": M[1, 1],
            "traction_shear_sign": sign,
            "crack_backend": "adaptive_czm",
            "crystal_theta_deg": a.crystal_theta_deg,
            "traction_probe_radius_m": a.traction_probe_radius_m,
            "calibration_probe_count": len(history),
            "basis_initial_alpha_deg": alpha0,
        }
        rows.append(row)
        for i, h in enumerate(history):
            all_hist.append({"target_psi_deg": target, "probe_index": i, **h})
        print({k: row.get(k) for k in (
            "target_psi_deg", "loading_alpha_deg", "loading_open_coeff",
            "loading_shear_coeff", "loading_open_is_tensile",
            "achieved_traction_phase_deg", "traction_phase_error_deg",
            "phase_converged", "calibration_probe_count")})

    csv_path = out / "mixed_mode_loading_calibration_v7.csv"
    with csv_path.open("w", newline="") as fp:
        cols = sorted({k for r in rows for k in r})
        w = csv.DictWriter(fp, fieldnames=cols)
        w.writeheader(); w.writerows(rows)
    pd.DataFrame(all_hist).to_csv(out / "mixed_mode_loading_calibration_history_v7.csv", index=False)
    (out / "production_backend_basis_v7.json").write_text(json.dumps({
        "calibration_id": CALIBRATION_ID,
        "model": MODEL_ID,
        "raw_response_matrix_Pa": Mraw.tolist(),
        "normalized_response_matrix_Pa": M.tolist(),
        "traction_shear_sign": sign,
        "basis_condition": cond,
        "reference_cleavage_shape": ref_c,
        "reference_slip_shape": ref_s,
        "coordinate": "full_circle_alpha_deg",
    }, indent=2))

    print("raw production-backend traction basis [Pa]:\n", Mraw)
    print("normalized basis [Pa]:\n", M)
    print("basis condition:", cond)
    print("Mode-I reference shapes:", {"cleavage": ref_c, "slip": ref_s})
    bad = [r for r in rows if not r["phase_converged"]]
    if bad and not a.force:
        raise SystemExit("v7 exact-backend full-circle calibration failed: " +
                         ", ".join(f"{r['target_psi_deg']}:err={r['traction_phase_error_deg']:.4g}" for r in bad))
    print("wrote", csv_path)


if __name__ == "__main__":
    main()
