#!/usr/bin/env python3
"""Compare the v10.0.5.18.4.0 theta=0 FEM control with PF v10.4.1.

The script aligns trajectories by projected crack extension rather than solver
step.  It reports R-curve shape, force-normalized J response, cleavage-clock
roughness, emitted/retained-state roughness, and exact run-contract differences.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SCHEMA = "theta0_pf_fem_parity_analysis_v10_0_5_18_4_0"


def _steps_file(case: Path) -> Path:
    direct = case / "steps_1000K.csv"
    if direct.is_file():
        return direct
    matches = sorted(case.rglob("steps_1000K.csv"))
    if len(matches) != 1:
        raise FileNotFoundError(
            f"expected one steps_1000K.csv below {case}; found {len(matches)}"
        )
    return matches[0]


def _json_if_present(path: Path) -> dict[str, Any]:
    if not path.is_file():
        return {}
    try:
        value = json.loads(path.read_text())
    except Exception:
        return {}
    return value if isinstance(value, dict) else {}


def _series(frame: pd.DataFrame, name: str, default: float = 0.0) -> np.ndarray:
    if name not in frame:
        return np.full(len(frame), float(default), dtype=float)
    return pd.to_numeric(frame[name], errors="coerce").fillna(default).to_numpy(float)


def _trajectory(case: Path, label: str) -> dict[str, Any]:
    steps = _steps_file(case)
    frame = pd.read_csv(steps)
    if frame.empty:
        raise ValueError(f"empty step table: {steps}")
    required = {"crack_extension_m", "KJ_Pa_sqrtm", "Ftop_N", "B", "N_em"}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"{steps} lacks required columns: {missing}")

    ext_um = 1.0e6 * _series(frame, "crack_extension_m")
    kj = 1.0e-6 * _series(frame, "KJ_Pa_sqrtm")
    force = _series(frame, "Ftop_N")
    B = _series(frame, "B")
    n_em = np.maximum(_series(frame, "N_em"), 0.0)
    n_fire = _series(frame, "n_fire")
    accepted = n_fire > 0.5
    if not np.any(accepted):
        accepted = np.r_[False, np.diff(ext_um) > 1.0e-9]
    positive_force = np.abs(force) > 1.0e-30
    kj_over_force = np.divide(
        kj,
        np.abs(force),
        out=np.full_like(kj, np.nan),
        where=positive_force,
    )

    order = np.argsort(ext_um, kind="stable")
    ext_sorted = ext_um[order]
    kj_sorted = kj[order]
    keep = np.r_[True, np.diff(ext_sorted) > 1.0e-12]
    ext_unique = ext_sorted[keep]
    kj_unique = kj_sorted[keep]

    accepted_index = np.where(accepted)[0]
    event_ext = ext_um[accepted_index]
    event_kj = kj[accepted_index]
    event_force = force[accepted_index]
    event_ratio = kj_over_force[accepted_index]

    if len(event_ext) >= 2 and np.ptp(event_ext) > 0.0:
        slope = float(np.polyfit(event_ext, event_kj, 1)[0])
    else:
        slope = float("nan")

    dB = np.diff(B)
    dlogN = np.diff(np.log10(1.0 + n_em))
    metrics = {
        "rows": int(len(frame)),
        "accepted_event_rows": int(np.count_nonzero(accepted)),
        "maximum_projected_extension_um": float(np.nanmax(ext_um)),
        "first_event_KJ_MPa_sqrt_m": (
            float(event_kj[0]) if len(event_kj) else None
        ),
        "final_event_KJ_MPa_sqrt_m": (
            float(event_kj[-1]) if len(event_kj) else None
        ),
        "mean_event_KJ_MPa_sqrt_m": (
            float(np.nanmean(event_kj)) if len(event_kj) else None
        ),
        "event_R_curve_linear_slope_MPa_sqrt_m_per_um": slope,
        "event_KJ_over_abs_force_initial": (
            float(event_ratio[0]) if len(event_ratio) else None
        ),
        "event_KJ_over_abs_force_final": (
            float(event_ratio[-1]) if len(event_ratio) else None
        ),
        "B_min": float(np.nanmin(B)),
        "B_max": float(np.nanmax(B)),
        "B_median_absolute_step_change": (
            float(np.nanmedian(np.abs(dB))) if len(dB) else 0.0
        ),
        "B_large_reset_count": int(np.count_nonzero(dB < -0.5)),
        "N_em_min": float(np.nanmin(n_em)),
        "N_em_max": float(np.nanmax(n_em)),
        "N_em_log10_dynamic_range": float(
            np.log10(1.0 + np.nanmax(n_em)) - np.log10(1.0 + np.nanmin(n_em))
        ),
        "N_em_median_absolute_log10_step_change": (
            float(np.nanmedian(np.abs(dlogN))) if len(dlogN) else 0.0
        ),
    }

    return {
        "label": label,
        "case": str(case.resolve()),
        "steps_file": str(steps.resolve()),
        "frame": frame,
        "ext_um": ext_um,
        "kj": kj,
        "force": force,
        "B": B,
        "n_em": n_em,
        "kj_over_force": kj_over_force,
        "accepted": accepted,
        "event_ext": event_ext,
        "event_kj": event_kj,
        "event_force": event_force,
        "event_ratio": event_ratio,
        "ext_unique": ext_unique,
        "kj_unique": kj_unique,
        "metrics": metrics,
    }


def _aligned(pf: dict[str, Any], fem: dict[str, Any], spacing_um: float) -> pd.DataFrame:
    maximum = min(
        float(np.nanmax(pf["ext_unique"])),
        float(np.nanmax(fem["ext_unique"])),
    )
    if maximum <= 0.0:
        raise ValueError("both trajectories must contain positive crack extension")
    spacing = max(float(spacing_um), 0.1)
    grid = np.arange(0.0, maximum + 0.5 * spacing, spacing)
    pf_k = np.interp(grid, pf["ext_unique"], pf["kj_unique"])
    fem_k = np.interp(grid, fem["ext_unique"], fem["kj_unique"])
    delta = fem_k - pf_k
    relative = np.divide(
        delta,
        np.maximum(np.abs(pf_k), 1.0e-30),
    )
    return pd.DataFrame(
        {
            "projected_extension_um": grid,
            "PF_KJ_MPa_sqrt_m": pf_k,
            "FEM_KJ_MPa_sqrt_m": fem_k,
            "FEM_minus_PF_KJ_MPa_sqrt_m": delta,
            "relative_error": relative,
        }
    )


def _contract_comparison(pf_case: Path, fem_case: Path) -> dict[str, Any]:
    pf_args = _json_if_present(pf_case / "run_args.json")
    fem_contract = _json_if_present(
        fem_case / "theta0_pf_parity_contract_v10_0_5_18_4_0.json"
    )
    matched = dict(fem_contract.get("matched_controls", {}) or {})
    pairs = {
        "theta_deg": (pf_args.get("crystal_theta_deg"), matched.get("theta_deg")),
        "dU_m": (pf_args.get("dU"), matched.get("dU_m")),
        "dt_s": (pf_args.get("dt"), matched.get("dt_s")),
        "n_stagger": (pf_args.get("n_stagger"), matched.get("n_stagger")),
        "tip_h_fine_m": (
            pf_args.get("tip_h_fine"),
            matched.get("tip_h_fine_m"),
        ),
        "tip_ratio": (pf_args.get("tip_ratio"), matched.get("tip_ratio")),
        "da_phys_m": (pf_args.get("da_phys"), matched.get("da_phys_m")),
        "adaptive_event_target": (
            pf_args.get("adaptive_event_target"),
            matched.get("adaptive_event_target"),
        ),
        "nx": (pf_args.get("nx"), matched.get("nx")),
        "ny": (pf_args.get("ny"), matched.get("ny")),
        "max_fronts": (pf_args.get("max_fronts"), matched.get("maximum_fronts")),
    }
    rows = {}
    for name, (left, right) in pairs.items():
        equal = False
        try:
            equal = math.isclose(float(left), float(right), rel_tol=1.0e-12, abs_tol=1.0e-15)
        except (TypeError, ValueError):
            equal = left == right
        rows[name] = {"PF": left, "FEM": right, "matched": bool(equal)}
    return {
        "fields": rows,
        "all_listed_numerical_controls_match": all(
            row["matched"] for row in rows.values()
        ),
        "intentional_model_form_differences": fem_contract.get(
            "intentional_model_form_differences", {}
        ),
    }


def _plot(out: Path, pf: dict[str, Any], fem: dict[str, Any]) -> None:
    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    ax.plot(pf["ext_um"], pf["kj"], label="PF v10.4.1")
    ax.plot(fem["ext_um"], fem["kj"], label="FEM/CZM v10.0.5.18.4.0")
    ax.set_xlabel("Projected crack extension (µm)")
    ax.set_ylabel(r"$K_J$ (MPa$\sqrt{m}$)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / "KJ_R_curve_comparison.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    ax.plot(pf["ext_um"], pf["B"], label="PF v10.4.1")
    ax.plot(fem["ext_um"], fem["B"], label="FEM/CZM v10.0.5.18.4.0")
    ax.set_xlabel("Projected crack extension (µm)")
    ax.set_ylabel("Cleavage progress B")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / "B_comparison.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    ax.plot(pf["ext_um"], 1.0 + pf["n_em"], label="PF v10.4.1")
    ax.plot(fem["ext_um"], 1.0 + fem["n_em"], label="FEM/CZM v10.0.5.18.4.0")
    ax.set_yscale("log")
    ax.set_xlabel("Projected crack extension (µm)")
    ax.set_ylabel(r"$1+N_{em}$")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / "N_em_comparison.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    ax.plot(pf["ext_um"], pf["kj_over_force"], label="PF v10.4.1")
    ax.plot(fem["ext_um"], fem["kj_over_force"], label="FEM/CZM v10.0.5.18.4.0")
    ax.set_xlabel("Projected crack extension (µm)")
    ax.set_ylabel(r"$K_J/|F|$")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / "KJ_over_force_comparison.png", dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pf-case", required=True, type=Path)
    parser.add_argument("--fem-case", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--spacing-um", type=float, default=5.0)
    args = parser.parse_args()

    pf_case = args.pf_case.expanduser().resolve()
    fem_case = args.fem_case.expanduser().resolve()
    out = args.out.expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)

    pf = _trajectory(pf_case, "PF")
    fem = _trajectory(fem_case, "FEM")
    aligned = _aligned(pf, fem, args.spacing_um)
    aligned.to_csv(out / "aligned_R_curve.csv", index=False)

    error = aligned["FEM_minus_PF_KJ_MPa_sqrt_m"].to_numpy(float)
    rel = aligned["relative_error"].to_numpy(float)
    pf_slope = pf["metrics"]["event_R_curve_linear_slope_MPa_sqrt_m_per_um"]
    fem_slope = fem["metrics"]["event_R_curve_linear_slope_MPa_sqrt_m_per_um"]
    slope_sign_same = bool(
        np.isfinite(pf_slope)
        and np.isfinite(fem_slope)
        and np.sign(pf_slope) == np.sign(fem_slope)
    )
    summary = {
        "schema": SCHEMA,
        "PF": {k: v for k, v in pf.items() if k in {"case", "steps_file", "metrics"}},
        "FEM": {k: v for k, v in fem.items() if k in {"case", "steps_file", "metrics"}},
        "aligned_comparison": {
            "maximum_common_extension_um": float(
                aligned["projected_extension_um"].max()
            ),
            "grid_spacing_um": float(args.spacing_um),
            "KJ_RMSE_MPa_sqrt_m": float(np.sqrt(np.mean(error * error))),
            "KJ_mean_absolute_error_MPa_sqrt_m": float(np.mean(np.abs(error))),
            "KJ_median_absolute_relative_error": float(np.median(np.abs(rel))),
            "R_curve_slope_sign_same": slope_sign_same,
            "PF_R_curve_slope_MPa_sqrt_m_per_um": pf_slope,
            "FEM_R_curve_slope_MPa_sqrt_m_per_um": fem_slope,
        },
        "run_contract": _contract_comparison(pf_case, fem_case),
        "interpretation_guards": {
            "theta0_removes_path_rotation": True,
            "FEM_adaptive_CZM_remeshing_still_active": True,
            "PF_and_FEM_bulk_constitutive_models_identical": False,
            "R_curve_difference_cannot_be_assigned_to_remeshing_alone_until_bulk_parity": True,
        },
    }
    (out / "parity_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True, default=str) + "\n"
    )
    _plot(out, pf, fem)

    print(json.dumps(summary["aligned_comparison"], indent=2, sort_keys=True))
    print(f"Wrote parity analysis to {out}")


if __name__ == "__main__":
    main()
