#!/usr/bin/env python3
"""Deterministic mixed-mode 2-D campaign for the v9.11 MPZ integration branch."""
from __future__ import annotations

import argparse
import csv
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
from arrhenius_fracture.mpz_parameterization_v911 import (
    compact_audit,
    load_selected_row,
    normalize_class_name,
)

RUNNER_ID = "mixed_mode_fem_czm_v9_11_MPZ_independent_PT_event_control"


def items(text, cast=str):
    return [cast(x) for x in str(text).replace(",", " ").split() if x]


def fs(value):
    value = float(value)
    return "inf" if math.isinf(value) else f"{value:.16g}"


def tag(value):
    return ("m" if value < 0 else "p") + f"{abs(value):05.1f}".replace(".", "p")


def angle_error(value, target):
    return (float(value) - float(target) + 180.0) % 360.0 - 180.0


def truthy(value):
    return str(value).strip().lower() in {"1", "true", "yes"}


def python_in_env(env):
    cp = subprocess.run(
        ["conda", "run", "-n", env, "python", "-c", "import sys;print(sys.executable)"],
        capture_output=True,
        text=True,
    )
    if cp.returncode:
        raise SystemExit(cp.stderr)
    return [x for x in cp.stdout.splitlines() if x.strip()][-1]


def response_matrix(cal_row):
    return np.array([
        [float(cal_row["response_11_Pa"]), float(cal_row["response_12_Pa"])],
        [float(cal_row["response_21_Pa"]), float(cal_row["response_22_Pa"])],
    ])


def manifest_for(parameter_root: Path, class_name: str) -> Path:
    cls = normalize_class_name(class_name)
    path = parameter_root / cls / "spatial_promotion_manifest.csv"
    if not path.exists():
        raise FileNotFoundError(path)
    return path


def build_case_command(py, args, manifest, class_name, cal_row, target, qo, qs, case_dir):
    return [
        py, "-m", "arrhenius_fracture.mixed_mode_first_passage_v9_11",
        "--mpz-material-manifest", str(manifest),
        "--mpz-material-class", class_name,
        "--mpz-length-um", fs(args.mpz_length_um),
        "--mpz-n-bins", str(args.mpz_n_bins),
        "--mpz-profile-sector-half-angle-deg", fs(args.mpz_profile_sector_half_angle_deg),
        "--mpz-profile-damage-cutoff", fs(args.mpz_profile_damage_cutoff),
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
        "--da-phys", fs(args.da_phys_m),
        "--j-decomposition", "cluster",
        "--rJ-cluster", fs(args.rJ_cluster_m),
        "--rJ-outer", fs(args.rJ_outer_m),
        "--temperatures", fs(args.T_K),
        "--crack-backend", "adaptive_czm", "--czm-max-angle-error-deg", "35",
        "--crystal-aniso", "--crystal-compete",
        "--crystal-theta-deg", fs(args.crystal_theta_deg),
        "--crystal-C11", fs(args.crystal_C11),
        "--crystal-C12", fs(args.crystal_C12),
        "--crystal-C44", fs(args.crystal_C44),
        "--cleave-gamma-aniso", fs(args.cleave_gamma_aniso),
        "--crystal-material", "w",
        "--multihit-m", "3", "--multihit-tau", "1e-6",
        "--sigma-cap-GPa", "0",
        "--save-snapshots", str(args.save_snapshots),
        *( [] if args.make_plots else ["--no-plots"] ),
        "--out", str(case_dir),
    ]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--parameter-root", default="mpz_v9_11_parameters")
    p.add_argument("--calibration-csv", required=True)
    p.add_argument("--classes", default="ceramic weakT DBTT")
    p.add_argument("--target-psi-deg", default="-60 -45 -30 -15 0 15 30 45 60")
    p.add_argument("--T-K", type=float, default=500.0)
    p.add_argument("--outroot", default="runs/mixed_mode_fem_czm_v9_11_MPZ_500K")
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
    p.add_argument("--mpz-length-um", type=float, default=100.0)
    p.add_argument("--mpz-n-bins", type=int, default=200)
    p.add_argument("--mpz-profile-sector-half-angle-deg", type=float, default=45.0)
    p.add_argument("--mpz-profile-damage-cutoff", type=float, default=0.85)
    p.add_argument("--da-phys-m", type=float, default=5e-6)
    p.add_argument("--rJ-cluster-m", type=float, default=20e-6)
    p.add_argument("--rJ-outer-m", type=float, default=25e-6)
    p.add_argument("--save-snapshots", type=int, default=0)
    p.add_argument("--make-plots", action="store_true")
    a = p.parse_args()

    py = python_in_env(a.conda_env)
    probe = subprocess.run(
        [py, "-c", "from arrhenius_fracture.mixed_mode_first_passage_v9_11 import MODEL_ID;print(MODEL_ID)"],
        capture_output=True,
        text=True,
    )
    if probe.returncode:
        raise SystemExit(probe.stderr)
    print("mixed-mode MPZ mechanics:", probe.stdout.strip())

    parameter_root = Path(a.parameter_root).resolve()
    materials = {}
    manifests = {}
    for name in items(a.classes):
        cls = normalize_class_name(name)
        manifest = manifest_for(parameter_root, cls)
        materials[cls] = load_selected_row(manifest, cls)
        manifests[cls] = manifest

    cal_rows = {
        round(float(r["target_psi_deg"]), 6): r
        for r in csv.DictReader(open(a.calibration_csv))
    }
    out = Path(a.outroot)
    out.mkdir(parents=True, exist_ok=True)

    def run_one(class_name, target):
        cls = normalize_class_name(class_name)
        cal = cal_rows[round(target, 6)]
        if not truthy(cal.get("phase_converged")) or not truthy(cal.get("first_production_step_verified")):
            raise RuntimeError(f"exact production-backend calibration not verified for {target}")
        row = materials[cls]
        manifest = manifests[cls]
        audit = compact_audit(row)
        root = out / cls / f"psi_{tag(target)}"
        root.mkdir(parents=True, exist_ok=True)
        final_path = root / "production_backend_control_final_summary.json"
        if final_path.exists() and not a.force:
            return json.loads(final_path.read_text())

        M = response_matrix(cal)
        alpha = float(cal.get("loading_alpha_unwrapped_deg", cal["loading_alpha_deg"]))
        history, candidates = [], []
        for iteration in range(max(1, a.max_control_iters)):
            qo, qs = loading_coefficients_from_alpha_deg(alpha)
            case = root / f"iter_{iteration:02d}"
            case.mkdir(parents=True, exist_ok=True)
            summary_path = case / "anisotropic_calibrated_tip_first_passage_summary.json"
            cmd = build_case_command(py, a, manifest, cls, cal, target, qo, qs, case)
            (case / "command.txt").write_text(shlex.join(cmd) + "\n")
            (case / "parameter_audit.json").write_text(json.dumps(audit, indent=2, default=str))
            with (case / "run.log").open("w") as log:
                rc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT).returncode
            if rc:
                history.append({
                    "iteration": iteration,
                    "loading_alpha_unwrapped_deg": alpha,
                    "loading_alpha_deg": wrap_loading_angle_deg(alpha),
                    "loading_open_coeff": qo,
                    "loading_shear_coeff": qs,
                    "status": "failed",
                    "return_code": rc,
                })
                break
            zsum = json.loads(summary_path.read_text())
            achieved = float(zsum.get("traction_phase_first_deg", np.nan))
            phase_err = angle_error(achieved, target) if np.isfinite(achieved) else float("nan")
            event = zsum.get("control_state") == "first_passage"
            phase_reliable = bool(zsum.get("traction_phase_probe_reliable", zsum.get("traction_probe_reliable", False)))
            directional_reliable = bool(zsum.get("directional_metrics_reliable", zsum.get("traction_probe_reliable", False)))
            h = {
                "iteration": iteration,
                "loading_alpha_unwrapped_deg": alpha,
                "loading_alpha_deg": wrap_loading_angle_deg(alpha),
                "loading_open_coeff": qo,
                "loading_shear_coeff": qs,
                "target_psi_deg": target,
                "achieved_psi_deg": achieved,
                "psi_error_deg": phase_err,
                "event_observed": event,
                "traction_phase_probe_reliable": phase_reliable,
                "directional_metrics_reliable": directional_reliable,
                "KJ_reference_first_MPa_sqrt_m": zsum.get("KJ_reference_first_MPa_sqrt_m"),
                "Kcleave_calibrated_first_MPa_sqrt_m": zsum.get("Kcleave_calibrated_first_MPa_sqrt_m"),
                "status": "event" if event else "right_censored",
            }
            history.append(h)
            candidates.append((zsum, h, case))
            if phase_reliable and directional_reliable and np.isfinite(phase_err) and abs(phase_err) <= a.event_psi_tol_deg:
                break
            if not np.isfinite(achieved):
                break
            anew = safeguarded_event_alpha_update(history, target, M, max_step_deg=a.max_alpha_step_deg)
            if abs(anew - alpha) < 1e-8:
                break
            alpha = anew

        pd.DataFrame(history).to_csv(root / "mixed_mode_control_history_v9_11.csv", index=False)
        if not candidates:
            raise RuntimeError(f"no completed iterations for {cls} target {target}")

        def rank(item):
            _, h, _ = item
            err = abs(float(h["psi_error_deg"])) if np.isfinite(float(h["psi_error_deg"])) else float("inf")
            controlled = bool(h["traction_phase_probe_reliable"] and h["directional_metrics_reliable"] and err <= a.event_psi_tol_deg)
            event = bool(h["event_observed"])
            return (0 if event and controlled else 1 if controlled else 2 if event else 3, err)

        zsum, h, case = min(candidates, key=rank)
        phase_ok = bool(
            h["traction_phase_probe_reliable"] and h["directional_metrics_reliable"] and
            np.isfinite(h["psi_error_deg"]) and abs(float(h["psi_error_deg"])) <= a.event_psi_tol_deg
        )
        event = bool(h["event_observed"])
        status = (
            "event" if event and phase_ok else
            "event_phase_mismatch" if event else
            "right_censored" if phase_ok else
            "right_censored_phase_mismatch"
        )
        final = {
            **zsum,
            **audit,
            "runner": RUNNER_ID,
            "class": cls,
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
            "calibration_achieved_phase_deg": float(cal["achieved_traction_phase_deg"]),
            "calibration_first_step_phase_error_deg": float(cal["traction_phase_error_deg"]),
            "calibration_backend": cal.get("crack_backend"),
        }
        final_path.write_text(json.dumps(final, indent=2, default=str))
        with (root / "production_backend_control_final_summary.csv").open("w", newline="") as fp:
            writer = csv.DictWriter(fp, fieldnames=list(final))
            writer.writeheader(); writer.writerow(final)
        return final

    jobs = [(name, target) for name in materials for target in items(a.target_psi_deg, float)]
    results = []
    with ThreadPoolExecutor(max_workers=max(1, a.max_jobs)) as executor:
        futures = {executor.submit(run_one, name, target): (name, target) for name, target in jobs}
        for future in as_completed(futures):
            name, target = futures[future]
            try:
                result = future.result()
            except Exception as exc:
                result = {"class": name, "target_psi_deg": target, "status": "failed", "error": repr(exc)}
            results.append(result)
            print({k: result.get(k) for k in (
                "class", "target_psi_deg", "status",
                "KJ_reference_first_MPa_sqrt_m",
                "Kcleave_calibrated_first_MPa_sqrt_m",
                "traction_phase_first_deg",
            )})

    results.sort(key=lambda z: (str(z.get("class")), float(z.get("target_psi_deg", 0))))
    pd.DataFrame(results).to_csv(out / "campaign_status_v9_11.csv", index=False)
    good = [r for r in results if r.get("status") != "failed"]
    if good:
        pd.DataFrame(good).to_csv(out / "mixed_mode_v9_11_MPZ_all_cases.csv", index=False)
    if any(r.get("status") == "failed" for r in results):
        raise SystemExit("one or more v9.11 cases failed")


if __name__ == "__main__":
    main()
