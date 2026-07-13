#!/usr/bin/env python3
"""Deterministic anisotropic mixed-mode campaign with exact backend full-circle control, v8."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shlex
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

from arrhenius_fracture.mixed_mode_first_passage_v8 import (
    loading_coefficients_from_alpha_deg,
    safeguarded_event_alpha_update,
    wrap_loading_angle_deg,
)

RUNNER_ID = "mixed_mode_fem_czm_v8_exact_backend_full_circle_event_control"
REQ = {
    "target_class", "exp_G00_eV", "exp_gT_eV_per_K", "exp_sigc0_GPa",
    "exp_sT_MPa_per_K", "exp_a", "exp_n", "exp_floor_frac",
    "cleave_G00_eV", "cleave_gT_eV_per_K", "cleave_sigc0_GPa",
    "cleave_sT_MPa_per_K", "cleave_exp_a", "cleave_exp_n",
    "cleave_floor_frac", "cleave_S_hs_kB", "chi_shield", "N_sat",
}


def items(s, cast=str):
    return [cast(x) for x in str(s).replace(",", " ").split() if x]


def fs(x):
    return "inf" if math.isinf(float(x)) else f"{float(x):.16g}"


def tag(x):
    return ("m" if x < 0 else "p") + f"{abs(x):05.1f}".replace(".", "p")


def angle_error(value, target):
    return (float(value) - float(target) + 180.0) % 360.0 - 180.0


def pyenv(env):
    cp = subprocess.run(
        ["conda", "run", "-n", env, "python", "-c", "import sys;print(sys.executable)"],
        capture_output=True, text=True)
    if cp.returncode:
        raise SystemExit(cp.stderr)
    return [x for x in cp.stdout.splitlines() if x.strip()][-1]


def truthy(v):
    return str(v).strip().lower() in {"1", "true", "yes"}


def barrier_audit(row):
    keys = [
        "exp_G00_eV", "exp_gT_eV_per_K", "exp_sigc0_GPa",
        "exp_sT_MPa_per_K", "exp_a", "exp_n", "exp_floor_frac",
        "cleave_G00_eV", "cleave_gT_eV_per_K", "cleave_sigc0_GPa",
        "cleave_sT_MPa_per_K", "cleave_exp_a", "cleave_exp_n",
        "cleave_floor_frac", "cleave_S_hs_kB", "chi_shield", "N_sat",
    ]
    out = {k: float(row[k]) for k in keys}
    out["barrier_fingerprint_sha256"] = hashlib.sha256(
        json.dumps(out, sort_keys=True).encode()).hexdigest()
    return out


def response_matrix(cal_row):
    return np.array([
        [float(cal_row["response_11_Pa"]), float(cal_row["response_12_Pa"])],
        [float(cal_row["response_21_Pa"]), float(cal_row["response_22_Pa"])],
    ])


def build_case_command(py, args, row, cal_row, target, qo, qs, case_dir):
    emitG = 0.75 * float(row.exp_G00_eV)
    emitg = 0.75 * float(row.exp_gT_eV_per_K)
    return [
        py, "-m", "arrhenius_fracture.mixed_mode_first_passage_v8",
        "--mixity-open-coeff", fs(qo),
        "--mixity-shear-coeff", fs(qs),
        "--target-traction-phase-deg", fs(target),
        "--traction-shear-sign", fs(cal_row["traction_shear_sign"]),
        "--traction-probe-radius-m", fs(args.traction_probe_radius_m),
        "--reference-cleavage-shape", fs(cal_row["reference_cleavage_shape"]),
        "--reference-slip-shape", fs(cal_row["reference_slip_shape"]),
        "--shear-emission-weight", fs(args.shear_emission_weight),
        "--directional-factor-max", fs(args.directional_factor_max),
        "--mode", "2d",
        "--nx", str(args.nx), "--ny", str(args.ny),
        "--tip-h-fine", fs(args.tip_h_fine), "--tip-ratio", fs(args.tip_ratio),
        "--dU", fs(args.dU), "--dt", fs(args.dt), "--steps", str(args.steps),
        "--n-stagger", "2", "--print-every", str(args.print_every),
        "--stop-after-first-fire", "--max-fronts", "1",
        "--adaptive-events", "--adaptive-event-target", ".25",
        "--adaptive-min-frac", "1e-8", "--adaptive-grow", "4",
        "--da-phys", "5e-6",
        "--j-decomposition", "cluster", "--rJ-cluster", "20e-6",
        "--rJ-outer", "25e-6", "--temperatures", fs(args.T_K),
        "--crack-backend", "adaptive_czm", "--czm-max-angle-error-deg", "35",
        "--crystal-aniso", "--crystal-compete",
        "--crystal-theta-deg", fs(args.crystal_theta_deg),
        "--crystal-C11", fs(args.crystal_C11),
        "--crystal-C12", fs(args.crystal_C12),
        "--crystal-C44", fs(args.crystal_C44),
        "--cleave-gamma-aniso", fs(args.cleave_gamma_aniso),
        "--crystal-material", "w",
        "--emit-barrier-kind", "exp_floor",
        "--emit-G00-eV", fs(emitG), "--emit-gT-eV-per-K", fs(emitg),
        "--emit-sigc0-GPa", fs(row.exp_sigc0_GPa),
        "--emit-sT-GPa-per-K", fs(float(row.exp_sT_MPa_per_K) / 1000.0),
        "--emit-exp-a", fs(row.exp_a), "--emit-exp-n", fs(row.exp_n),
        "--emit-floor-frac", fs(row.exp_floor_frac), "--emit-Tref-K", "300",
        "--cleave-barrier-kind", "exp_floor", "--cleave-exp-T-mode", "linear",
        "--cleave-G00-eV", fs(row.cleave_G00_eV),
        "--cleave-gT-eV-per-K", fs(row.cleave_gT_eV_per_K),
        "--cleave-sigc0-GPa", fs(row.cleave_sigc0_GPa),
        "--cleave-sT-GPa-per-K", fs(float(row.cleave_sT_MPa_per_K) / 1000.0),
        "--cleave-exp-a", fs(row.cleave_exp_a),
        "--cleave-exp-n", fs(row.cleave_exp_n),
        "--cleave-floor-frac", fs(row.cleave_floor_frac),
        "--cleave-S-hs-kB", fs(row.cleave_S_hs_kB),
        "--cleave-sigma-S-GPa", "6", "--cleave-S-hs-power", "2",
        "--cleave-S-hs-Tref-K", "300", "--cleave-Tref-K", "300",
        "--cleave-shield-chi", fs(row.chi_shield), "--n-sat", fs(row.N_sat),
        "--multihit-m", "3", "--multihit-tau", "1e-6",
        "--emb-sat-frac", "1", "--save-snapshots", "0", "--no-plots",
        "--out", str(case_dir),
    ]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--parameter-table", default="four_class_exp_floor_exact_model_inputs.csv")
    p.add_argument("--calibration-csv", required=True)
    p.add_argument("--classes", default="ceramic DBTT")
    p.add_argument("--target-psi-deg", default="-60 -45 -30 -15 0 15 30 45 60")
    p.add_argument("--T-K", type=float, default=500)
    p.add_argument("--outroot", default="runs/mixed_mode_fem_czm_v8_production_backend_500K")
    p.add_argument("--conda-env", default="arrhenius-fem-czm")
    p.add_argument("--max-jobs", type=int, default=1)
    p.add_argument("--force", action="store_true")
    p.add_argument("--nx", type=int, default=24)
    p.add_argument("--ny", type=int, default=48)
    p.add_argument("--tip-h-fine", type=float, default=3e-6)
    p.add_argument("--tip-ratio", type=float, default=1.25)
    p.add_argument("--dU", type=float, default=2e-7)
    p.add_argument("--dt", type=float, default=8.4)
    p.add_argument("--steps", type=int, default=3000)
    p.add_argument("--print-every", type=int, default=50)
    p.add_argument("--crystal-theta-deg", type=float, default=45.0)
    p.add_argument("--crystal-C11", type=float, default=523e9)
    p.add_argument("--crystal-C12", type=float, default=203e9)
    p.add_argument("--crystal-C44", type=float, default=160e9)
    p.add_argument("--cleave-gamma-aniso", type=float, default=0.3)
    p.add_argument("--traction-probe-radius-m", type=float, default=10e-6)
    p.add_argument("--shear-emission-weight", type=float, default=1.0)
    p.add_argument("--directional-factor-max", type=float, default=5.0)
    p.add_argument("--event-psi-tol-deg", type=float, default=2.0)
    p.add_argument("--max-control-iters", type=int, default=6)
    p.add_argument("--max-alpha-step-deg", type=float, default=20.0)
    a = p.parse_args()

    py = pyenv(a.conda_env)
    probe = subprocess.run(
        [py, "-c", "from arrhenius_fracture.mixed_mode_first_passage_v8 import MODEL_ID;print(MODEL_ID)"],
        capture_output=True, text=True)
    if probe.returncode:
        raise SystemExit(probe.stderr)
    print("mixed-mode mechanics:", probe.stdout.strip())

    df = pd.read_csv(a.parameter_table)
    missing = REQ - set(df.columns)
    if missing:
        raise SystemExit(f"parameter table missing {sorted(missing)}")
    df = df.set_index("target_class", drop=False)
    cal_rows = {round(float(r["target_psi_deg"]), 6): r
                for r in csv.DictReader(open(a.calibration_csv))}
    out = Path(a.outroot)
    out.mkdir(parents=True, exist_ok=True)

    def run_one(class_name, target):
        cal = cal_rows[round(target, 6)]
        if not truthy(cal.get("phase_converged")) or not truthy(cal.get("first_production_step_verified")):
            raise RuntimeError(f"exact production-backend calibration not verified for {target}")
        row = df.loc[class_name]
        audit = barrier_audit(row)
        root = out / class_name / f"psi_{tag(target)}"
        root.mkdir(parents=True, exist_ok=True)
        final_path = root / "production_backend_control_final_summary.json"
        if final_path.exists() and not a.force:
            return json.loads(final_path.read_text())

        M = response_matrix(cal)
        alpha = float(cal.get("loading_alpha_unwrapped_deg", cal["loading_alpha_deg"]))
        history = []
        candidates = []
        for iteration in range(max(1, a.max_control_iters)):
            qo, qs = loading_coefficients_from_alpha_deg(alpha)
            case = root / f"iter_{iteration:02d}"
            case.mkdir(parents=True, exist_ok=True)
            summary_path = case / "anisotropic_calibrated_tip_first_passage_summary.json"
            cmd = build_case_command(py, a, row, cal, target, qo, qs, case)
            (case / "command.txt").write_text(shlex.join(cmd) + "\n")
            (case / "barrier_audit.json").write_text(json.dumps({"class": class_name, **audit}, indent=2))
            with (case / "run.log").open("w") as log:
                rc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT).returncode
            if rc:
                history.append({"iteration": iteration,
                                "loading_alpha_unwrapped_deg": alpha,
                                "loading_alpha_deg": wrap_loading_angle_deg(alpha),
                                "loading_open_coeff": qo, "loading_shear_coeff": qs,
                                "status": "failed", "return_code": rc})
                break
            zsum = json.loads(summary_path.read_text())
            achieved = float(zsum.get("traction_phase_first_deg", np.nan))
            phase_err = angle_error(achieved, target) if np.isfinite(achieved) else float("nan")
            event = zsum.get("control_state") == "first_passage"
            phase_reliable = bool(zsum.get("traction_phase_probe_reliable",
                                           zsum.get("traction_probe_reliable", False)))
            directional_reliable = bool(zsum.get("directional_metrics_reliable",
                                                 zsum.get("traction_probe_reliable", False)))
            rowhist = {
                "iteration": iteration,
                "loading_alpha_unwrapped_deg": alpha,
                "loading_alpha_deg": wrap_loading_angle_deg(alpha),
                "loading_open_coeff": qo,
                "loading_shear_coeff": qs,
                "loading_angle_deg": math.degrees(math.atan2(qs, qo)),
                "target_psi_deg": target,
                "achieved_psi_deg": achieved,
                "psi_error_deg": phase_err,
                "event_observed": event,
                "traction_phase_probe_reliable": phase_reliable,
                "directional_metrics_reliable": directional_reliable,
                "traction_probe_reliable": bool(phase_reliable and directional_reliable),
                "KJ_reference_first_MPa_sqrt_m": zsum.get("KJ_reference_first_MPa_sqrt_m"),
                "Kcleave_calibrated_first_MPa_sqrt_m": zsum.get("Kcleave_calibrated_first_MPa_sqrt_m"),
                "cleavage_factor_first": zsum.get("cleavage_factor_first"),
                "emission_factor_first": zsum.get("emission_factor_first"),
                "mode_classification": zsum.get("mode_classification"),
                "status": "event" if event else "right_censored",
            }
            history.append(rowhist)
            candidates.append((zsum, rowhist, case))
            if (phase_reliable and directional_reliable and np.isfinite(phase_err) and
                    abs(phase_err) <= a.event_psi_tol_deg):
                break
            if not np.isfinite(achieved):
                break
            anew = safeguarded_event_alpha_update(
                history, target, M, max_step_deg=a.max_alpha_step_deg)
            if abs(anew - alpha) < 1e-8:
                break
            alpha = anew

        pd.DataFrame(history).to_csv(root / "mixed_mode_control_history_v8.csv", index=False)
        if not candidates:
            raise RuntimeError(f"no completed iterations for {class_name} target {target}")
        # Prefer an event with controlled phase. Then any controlled endpoint,
        # then observed events, then minimum phase error.
        def rank(item):
            _, h, _ = item
            e = abs(float(h["psi_error_deg"])) if np.isfinite(float(h["psi_error_deg"])) else float("inf")
            ok = bool(h["traction_phase_probe_reliable"] and
                      h["directional_metrics_reliable"] and
                      e <= a.event_psi_tol_deg)
            ev = bool(h["event_observed"])
            return (0 if (ev and ok) else 1 if ok else 2 if ev else 3, e)
        selected = min(candidates, key=rank)
        zsum, h, case = selected
        phase_ok = bool(h["traction_phase_probe_reliable"] and
                        h["directional_metrics_reliable"] and
                        np.isfinite(h["psi_error_deg"]) and
                        abs(float(h["psi_error_deg"])) <= a.event_psi_tol_deg)
        event = bool(h["event_observed"])
        status = ("event" if event and phase_ok else
                  "event_phase_mismatch" if event else
                  "right_censored" if phase_ok else
                  "right_censored_phase_mismatch")
        final = {
            **zsum, **audit,
            "runner": RUNNER_ID,
            "class": class_name,
            "target_psi_deg": target,
            "status": status,
            "selected_iteration": int(h["iteration"]),
            "control_iterations_run": len(history),
            "event_phase_control_converged": phase_ok,
            "selected_case_dir": str(case),
            "selected_loading_alpha_unwrapped_deg": float(h["loading_alpha_unwrapped_deg"]),
            "selected_loading_alpha_deg": float(h["loading_alpha_deg"]),
            "selected_loading_open_coeff": float(h["loading_open_coeff"]),
            "selected_loading_shear_coeff": float(h["loading_shear_coeff"]),
            "calibration_id": cal.get("calibration_id"),
            "calibration_first_production_step_verified": truthy(cal.get("first_production_step_verified")),
            "calibration_loading_alpha_deg": float(cal["loading_alpha_deg"]),
            "calibration_loading_open_is_tensile": truthy(cal.get("loading_open_is_tensile", True)),
            "calibration_achieved_phase_deg": float(cal["achieved_traction_phase_deg"]),
            "calibration_first_step_phase_error_deg": float(cal["traction_phase_error_deg"]),
            "calibration_backend": cal.get("crack_backend"),
        }
        final_path.write_text(json.dumps(final, indent=2, default=str))
        with (root / "production_backend_control_final_summary.csv").open("w", newline="") as fp:
            w = csv.DictWriter(fp, fieldnames=list(final)); w.writeheader(); w.writerow(final)
        return final

    jobs = [(k, t) for k in items(a.classes) for t in items(a.target_psi_deg, float)]
    results = []
    with ThreadPoolExecutor(max_workers=max(1, a.max_jobs)) as ex:
        futures = {ex.submit(run_one, k, t): (k, t) for k, t in jobs}
        for future in as_completed(futures):
            k, t = futures[future]
            try:
                z = future.result()
            except Exception as exc:
                z = {"class": k, "target_psi_deg": t, "status": "failed", "error": repr(exc)}
            results.append(z)
            print({q: z.get(q) for q in (
                "class", "target_psi_deg", "status", "selected_loading_alpha_deg",
                "KJ_reference_first_MPa_sqrt_m", "Kcleave_calibrated_first_MPa_sqrt_m",
                "traction_phase_first_deg", "cleavage_factor_first", "B_final",
                "mode_classification")})

    results.sort(key=lambda z: (str(z.get("class")), float(z.get("target_psi_deg", 0))))
    pd.DataFrame(results).to_csv(out / "campaign_status_v8.csv", index=False)
    good = [x for x in results if x.get("status") != "failed"]
    if good:
        pd.DataFrame(good).to_csv(out / "mixed_mode_v8_anisotropic_all_cases.csv", index=False)
    if any(x.get("status") == "failed" for x in results):
        raise SystemExit("one or more v8 cases failed")


if __name__ == "__main__":
    main()
