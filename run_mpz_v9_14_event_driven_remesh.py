#!/usr/bin/env python3
"""Run v9.14 through the conservative adaptive-CZM event path.

The inner solver advances one additional physical increment beyond the requested
analysis endpoint.  That guard increment allows the event exactly at the
requested endpoint to receive its zero-time/zero-load post-event equilibrium
solve before the inherited target-extension break is reached.  Guard data are
preserved separately and excluded from the reported R-curve files.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
import shutil

import numpy as np
import pandas as pd

import run_mpz_v9_13_deterministic_material_transfer as _base
from arrhenius_fracture.rcurve_postprocess_v911 import cascade_metrics
from arrhenius_fracture.remesh_audit_v914 import audit_case


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


def _build_command_v914(args, class_name, run_root, force_rerun):
    cmd = _original_build_command(args, class_name, run_root, force_rerun)
    old = "run_mpz_v9_13_mode_i_rcurve.py"
    new = "run_mpz_v9_14_mode_i_rcurve.py"
    try:
        cmd[cmd.index(old)] = new
    except ValueError as exc:
        raise RuntimeError(f"v9.13 driver token not found in command: {cmd}") from exc

    solver_target, _ = _solver_guard_target_um(args)
    if math.isfinite(solver_target) and solver_target > 0.0:
        _replace_option(cmd, "--target-extension-um", solver_target)
    return cmd


def _copy_guard_full(path: Path) -> Path | None:
    if not path.exists():
        return None
    full = path.with_name(path.stem + "_solver_guard_full" + path.suffix)
    shutil.copy2(path, full)
    return full


def _trim_guard_rcurve_outputs(case_dir: Path, target_um: float) -> dict:
    """Keep solver-guard data separately and expose only the analysis interval."""
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
        end_name = "crack_extension_end_um" if "crack_extension_end_um" in clustered else "crack_extension_um"
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
            ax.set_title(fr"Event-driven R-curve through {target_um:g} $\mu$m")
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


def _run_case_v914(args, base_seed, class_name, root):
    row = _original_run_case(args, base_seed, class_name, root)
    case_dir = Path(row["case_dir"])
    target = float(args.target_extension_um)
    solver_target, guard = _solver_guard_target_um(args)

    trim = _trim_guard_rcurve_outputs(case_dir, target)
    audit = audit_case(
        case_dir,
        args.T_K,
        analysis_target_extension_um=target,
    )

    solver_final = row.get("final_extension_um")
    try:
        solver_final_f = float(solver_final)
    except (TypeError, ValueError):
        solver_final_f = np.nan
    analysis_final = min(solver_final_f, target) if np.isfinite(solver_final_f) else np.nan
    row.update({
        "solver_guard_target_extension_um": float(solver_target),
        "solver_guard_increment_um": float(guard),
        "solver_final_extension_um": solver_final_f,
        "analysis_target_extension_um": target,
        "analysis_final_extension_um": analysis_final,
        "target_completed": bool(np.isfinite(solver_final_f) and solver_final_f + 1e-8 >= target),
        **trim,
    })
    # The public final_extension field describes the analyzed interval; the full
    # solver value remains available as solver_final_extension_um.
    if np.isfinite(analysis_final):
        row["final_extension_um"] = analysis_final
    row.update({f"v914_{k}": v for k, v in audit.items() if k not in {"case_dir", "T_K"}})

    protocol = {
        "schema": "v9_14_terminal_equilibrium_guard_v1",
        "analysis_target_extension_um": target,
        "solver_guard_target_extension_um": solver_target,
        "solver_guard_increment_um": guard,
        "solver_final_extension_um": solver_final_f,
        "analysis_final_extension_um": analysis_final,
        "guard_data_policy": (
            "full solver data preserved with _solver_guard_full suffix; standard "
            "R-curve CSV/plot files are trimmed to the requested analysis target"
        ),
    }
    (case_dir / "analysis_window_v914.json").write_text(json.dumps(protocol, indent=2))
    (case_dir / "v9_14_case_summary.json").write_text(json.dumps(row, indent=2, default=str))
    pd.DataFrame([row]).to_csv(case_dir / "v9_14_case_summary.csv", index=False)
    return row


def main():
    original_build = _base.build_command
    original_run = _base.run_case
    _base.build_command = _build_command_v914
    _base.run_case = _run_case_v914
    try:
        return _base.main()
    finally:
        _base.build_command = original_build
        _base.run_case = original_run


if __name__ == "__main__":
    main()
