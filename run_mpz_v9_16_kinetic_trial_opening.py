#!/usr/bin/env python3
"""Run the v9.16 kinetic trial-cohesive/MPZ material-transfer gate."""
from __future__ import annotations

import json
import math
import os
from pathlib import Path
import shutil

import numpy as np
import pandas as pd

import run_mpz_v9_13_deterministic_material_transfer as _base
from arrhenius_fracture.coupled_event_audit_v916 import audit_case_v916
from arrhenius_fracture.rcurve_postprocess_v911 import cascade_metrics

_original_build_command = _base.build_command
_original_run_case = _base.run_case


def _replace_option(cmd: list[str], name: str, value: float | str) -> None:
    for i, token in enumerate(cmd):
        if token == name and i + 1 < len(cmd):
            cmd[i + 1] = str(value)
            return
        if token.startswith(name + "="):
            cmd[i] = f"{name}={value}"
            return
    raise RuntimeError(f"required option {name} not found in command: {cmd}")


def _solver_guard_target_um(args) -> tuple[float, float]:
    target = float(args.target_extension_um)
    guard = max(float(args.da_phys_um), 1.0e-9)
    if not math.isfinite(target) or target <= 0.0:
        return target, 0.0
    return target + guard, guard


def _bulk_mode() -> str:
    mode = os.environ.get("ARRHENIUS_BULK_PLASTICITY_MODE", "tip_only").strip()
    if mode not in {"tip_only", "bulk_same_pt_km"}:
        raise SystemExit(
            "ARRHENIUS_BULK_PLASTICITY_MODE must be tip_only or bulk_same_pt_km"
        )
    return mode


def _build_command_v916(args, class_name, run_root, force_rerun):
    cmd = _original_build_command(args, class_name, run_root, force_rerun)
    old = "run_mpz_v9_13_mode_i_rcurve.py"
    new = "run_mpz_v9_16_mode_i_rcurve.py"
    try:
        cmd[cmd.index(old)] = new
    except ValueError as exc:
        raise RuntimeError(f"v9.13 driver token not found in command: {cmd}") from exc
    solver_target, _ = _solver_guard_target_um(args)
    if math.isfinite(solver_target) and solver_target > 0.0:
        _replace_option(cmd, "--target-extension-um", solver_target)
    _replace_option(cmd, "--bulk-plasticity-mode", _bulk_mode())
    return cmd


def _copy_guard_full(path: Path) -> Path | None:
    if not path.exists():
        return None
    full = path.with_name(path.stem + "_solver_guard_full" + path.suffix)
    shutil.copy2(path, full)
    return full


def _trim_guard_rcurve_outputs(case_dir: Path, target_um: float) -> dict:
    if not math.isfinite(target_um) or target_um <= 0.0:
        return {"guard_trim_applied": False}
    tol = max(1.0e-8, 1.0e-8 * abs(target_um))
    raw_path = case_dir / "R_curve_topology_events_raw.csv"
    cluster_path = case_dir / "R_curve_load_events_clustered.csv"
    compat_path = case_dir / "R_curve_event_sampled.csv"
    metrics_path = case_dir / "R_curve_cascade_metrics.csv"
    plot_path = case_dir / "R_curve_cascade_aware.png"
    for path in (raw_path, cluster_path, compat_path, metrics_path, plot_path):
        _copy_guard_full(path)

    def read(path: Path) -> pd.DataFrame:
        if not path.exists() or path.stat().st_size == 0:
            return pd.DataFrame()
        try:
            return pd.read_csv(path)
        except (pd.errors.EmptyDataError, OSError, ValueError):
            return pd.DataFrame()

    raw = read(raw_path)
    clustered = read(cluster_path)
    if not raw.empty and "crack_extension_after_um" in raw:
        q = pd.to_numeric(raw["crack_extension_after_um"], errors="coerce")
        raw = raw[q <= target_um + tol].copy()
    if not clustered.empty:
        end_name = (
            "crack_extension_end_um"
            if "crack_extension_end_um" in clustered
            else "crack_extension_um"
        )
        if end_name in clustered:
            q = pd.to_numeric(clustered[end_name], errors="coerce")
            clustered = clustered[q <= target_um + tol].copy()
    raw.to_csv(raw_path, index=False)
    clustered.to_csv(cluster_path, index=False)
    clustered.to_csv(compat_path, index=False)
    metrics = cascade_metrics(raw, clustered)
    pd.DataFrame([metrics]).to_csv(metrics_path, index=False)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        if not clustered.empty:
            fig, ax = plt.subplots(figsize=(7.4, 4.9))
            ax.plot(
                clustered["crack_extension_um"],
                clustered["KJ_onset_MPa_sqrt_m"],
                marker="o", linewidth=1.0, markersize=3,
            )
            ax.set_xlabel(r"Projected crack extension $\Delta a_x$ ($\mu$m)")
            ax.set_ylabel(r"$K_J$ at load-event onset (MPa$\sqrt{m}$)")
            ax.set_title(fr"v9.16 trial-topology R-curve through {target_um:g} $\mu$m")
            ax.grid(alpha=0.25)
            fig.tight_layout()
            fig.savefig(plot_path, dpi=220)
            plt.close(fig)
    except Exception:
        pass

    return {
        "guard_trim_applied": True,
        "analysis_target_extension_um": float(target_um),
        "n_analysis_raw_topology_events": int(len(raw)),
        "n_analysis_independent_load_events": int(len(clustered)),
        **metrics,
    }


def _run_case_v916(args, base_seed, class_name, root):
    # The inherited completion marker is topology-based.  Never reuse it when a
    # prior v9.16 audit says the kinetic event was not physically committed.
    expected = (
        Path(root) / f"seed_{base_seed}" / "tip_only" / str(class_name)
        / f"T{int(round(args.T_K))}_th{args.crystal_theta_deg:g}"
    )
    marker = expected / ".long_growth_complete"
    prior_audit_path = expected / "kinetic_trial_event_audit_v916.json"
    if marker.exists():
        try:
            prior = json.loads(prior_audit_path.read_text())
        except Exception:
            prior = {}
        if not bool(prior.get("target_committed", False)):
            marker.unlink()

    row = _original_run_case(args, base_seed, class_name, root)
    case_dir = Path(row["case_dir"])
    target = float(args.target_extension_um)
    solver_target, guard = _solver_guard_target_um(args)
    trim = _trim_guard_rcurve_outputs(case_dir, target)
    audit = audit_case_v916(case_dir, args.T_K, target)

    try:
        solver_final = float(row.get("final_extension_um"))
    except (TypeError, ValueError):
        solver_final = np.nan
    committed_final = float(audit.get("committed_extension_um", 0.0) or 0.0)
    committed_complete = bool(audit.get("target_committed", False))
    if int(row.get("subprocess_returncode", row.get("returncode", 0)) or 0) == 0:
        row["status"] = "complete" if committed_complete else "right_censored_trial_uncommitted"

    row.update({
        "requested_bulk_plasticity_mode": _bulk_mode(),
        "solver_guard_target_extension_um": float(solver_target),
        "solver_guard_increment_um": float(guard),
        "solver_topology_extension_um": solver_final,
        "analysis_target_extension_um": target,
        "analysis_committed_extension_um": committed_final,
        "target_completed": bool(audit.get("target_committed", False)),
        **trim,
    })
    row["final_extension_um"] = committed_final
    row.update({f"v916_{k}": v for k, v in audit.items() if k not in {"case_dir", "T_K"}})

    protocol = {
        "schema": "v9_16_kinetic_trial_guard_v1",
        "analysis_target_extension_um": target,
        "solver_guard_target_extension_um": solver_target,
        "solver_guard_increment_um": guard,
        "solver_topology_extension_um": solver_final,
        "analysis_committed_extension_um": committed_final,
        "bulk_plasticity_mode": _bulk_mode(),
        "data_policy": (
            "solver topology outputs are retained; v9.16 physical completion is "
            "reported from committed kinetic trial events"
        ),
    }
    (case_dir / "analysis_window_v916.json").write_text(json.dumps(protocol, indent=2))
    (case_dir / "v9_16_case_summary.json").write_text(json.dumps(row, indent=2, default=str))
    pd.DataFrame([row]).to_csv(case_dir / "v9_16_case_summary.csv", index=False)
    return row


def main():
    original_build = _base.build_command
    original_run = _base.run_case
    _base.build_command = _build_command_v916
    _base.run_case = _run_case_v916
    try:
        return _base.main()
    finally:
        _base.build_command = original_build
        _base.run_case = original_run


if __name__ == "__main__":
    main()
