#!/usr/bin/env python3
"""Analyze v10.0.5.17 audited paper-parameter FEM/CZM campaigns.

The analysis uses accepted physical cleavage events only. K_J is read from the
step history; J is read directly from the signed-positive directional domain
integral in the root-front history rather than reconstructed from K. Crack-path
CSVs provide committed segment directions and path-deflection diagnostics.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager, nullcontext
import hashlib
import json
from pathlib import Path
from typing import Any, Iterator

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

MANIFEST_NAME = "persistent_site_production_manifest_v10_0_5_17.json"
TRANSFER_NAME = "audited_PF_paper_parameter_transfer_v10_0_5_17.json"
EVENTS_NAME = "stochastic_geometry_events_v10_0_5_16.json"
SCHEMA = "v10.0.5.17.2_four_class_K_J_path_campaign_analysis"


def _load_json(path: Path, default: Any) -> Any:
    return json.loads(path.read_text()) if path.is_file() else default


def _first(case: Path, pattern: str) -> Path:
    matches = sorted(case.glob(pattern))
    if not matches:
        raise FileNotFoundError(f"missing {pattern} in {case}")
    return matches[0]


def _temperature(case: Path) -> int:
    name = case.name
    if len(name) != 6 or not name.startswith("T") or not name.endswith("K"):
        raise ValueError(f"invalid temperature directory {case}")
    return int(name[1:5])


def _accepted(events: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        dict(row)
        for row in events
        if bool(row.get("inserted", False)) and float(row.get("moved_m", 0.0)) > 0.0
    ]


def _fingerprint(events: list[dict[str, Any]]) -> str:
    payload = [
        {
            "event_index": int(row.get("hazard_event_index", index)),
            "threshold_action": float(row.get("threshold_action", 0.0)),
            "moved_m": float(row.get("moved_m", 0.0)),
        }
        for index, row in enumerate(events)
    ]
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _seed(manifest: dict[str, Any], events: list[dict[str, Any]]) -> int | None:
    physics = dict(manifest.get("physics_contract", {}) or {})
    value = physics.get("cleavage_hazard_seed")
    if value is None and events:
        value = events[0].get("hazard_seed")
    return None if value is None else int(value)


def _save(fig: plt.Figure, stem: Path) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def _short_label(option: str) -> str:
    labels = {
        "v913_paper_peak01_0242980_persistent_sites": "Peak 0242980",
        "v913_paper_peak02_0127508_persistent_sites": "Peak 0127508",
        "v913_paper_peak03_0115460_persistent_sites": "Peak 0115460",
        "v913_paper_dbtt01_0202500_persistent_sites": "DBTT 0202500",
        "v913_paper_dbtt02_0088403_persistent_sites": "DBTT 0088403",
        "v913_paper_weakT01_0257068_persistent_sites": "Weak-T 0257068",
        "v913_paper_ceramic01_0189364_persistent_sites": "Ceramic 0189364",
        "v913_paper_control01_0086420_persistent_sites": "Control 0086420",
    }
    return labels.get(option, option)


def _angle_wrap_deg(value: np.ndarray | float) -> np.ndarray | float:
    return (np.asarray(value) + 180.0) % 360.0 - 180.0


def _load_root_front_events(case: Path, event_steps: np.ndarray) -> pd.DataFrame:
    fronts = pd.read_csv(_first(case, "fronts_*K.csv"))
    required = {
        "step", "front_id", "n_fire", "J_signed_trial", "J_effective_trial",
        "x_m", "y_m",
    }
    missing = required.difference(fronts.columns)
    if missing:
        raise ValueError(f"{case}: missing root-front J/path columns {sorted(missing)}")
    root = fronts.loc[(fronts["front_id"] == 0) & (fronts["n_fire"] > 0)].copy()
    root = root.sort_values("step").drop_duplicates("step", keep="last").set_index("step")
    missing_steps = [int(step) for step in event_steps if int(step) not in root.index]
    if missing_steps:
        raise ValueError(f"{case}: no root-front direct J record for event steps {missing_steps[:8]}")
    return root.loc[[int(step) for step in event_steps]].reset_index()


def _load_path(case: Path, temperature: int, n_events: int) -> pd.DataFrame:
    exact = case / f"crack_path_{temperature}K.csv"
    path_file = exact if exact.is_file() else _first(case, "crack_path_[0-9]*K.csv")
    path = pd.read_csv(path_file)
    required = {"x_m", "y_m"}
    missing = required.difference(path.columns)
    if missing:
        raise ValueError(f"{case}: missing crack path columns {sorted(missing)}")
    if len(path) != n_events + 1:
        raise ValueError(
            f"{case}: crack path has {len(path)-1} segments but {n_events} accepted events"
        )
    return path


def _requested_angles(case: Path, event_steps: np.ndarray) -> np.ndarray:
    matches = sorted(case.glob("branch_diagnostics_*K.csv"))
    if not matches:
        return np.full(len(event_steps), np.nan)
    diag = pd.read_csv(matches[0])
    if not {"step", "angle1_deg"}.issubset(diag.columns):
        return np.full(len(event_steps), np.nan)
    diag = diag.sort_values("step").drop_duplicates("step", keep="last").set_index("step")
    return np.asarray(
        [float(diag.loc[int(step), "angle1_deg"]) if int(step) in diag.index else np.nan
         for step in event_steps],
        dtype=float,
    )


def load_case(case: Path, target_um: float) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    temperature = _temperature(case)
    steps = pd.read_csv(_first(case, "steps_*K.csv"))
    required = {"step", "n_fire", "KJ_Pa_sqrtm", "crack_extension_m"}
    missing = required.difference(steps.columns)
    if missing:
        raise ValueError(f"{case}: missing step columns {sorted(missing)}")
    rows = steps.loc[steps["n_fire"] > 0].copy().reset_index(drop=True)
    if rows.empty:
        raise ValueError(f"{case}: no accepted cleavage events")

    raw_events = _load_json(case / EVENTS_NAME, [])
    accepted = _accepted(raw_events)
    if len(accepted) != len(rows):
        raise ValueError(
            f"{case}: {len(rows)} n_fire rows but {len(accepted)} accepted geometry events"
        )

    manifest = _load_json(case / MANIFEST_NAME, {})
    transfer = _load_json(case / TRANSFER_NAME, {})
    option = str(
        transfer.get("parameter_option")
        or manifest.get("parameter_option")
        or case.parent.name
    )
    candidate = str(
        transfer.get("candidate_id") or manifest.get("candidate_id") or "unknown"
    )
    material_class = str(
        transfer.get("material_class") or manifest.get("material_class") or "unknown"
    )
    role = str(
        transfer.get("paper_role") or manifest.get("paper_role") or option
    )
    seed = _seed(manifest, accepted)

    moved_um = np.asarray([float(row["moved_m"]) * 1.0e6 for row in accepted])
    thresholds = np.asarray(
        [float(row.get("threshold_action", np.nan)) for row in accepted]
    )
    event_steps = rows["step"].astype(int).to_numpy()
    root = _load_root_front_events(case, event_steps)
    path = _load_path(case, temperature, len(rows))
    requested_angle = _requested_angles(case, event_steps)

    xy = path[["x_m", "y_m"]].to_numpy(float)
    dxy = np.diff(xy, axis=0)
    seg_len = np.linalg.norm(dxy, axis=1)
    committed_angle = np.degrees(np.arctan2(dxy[:, 1], dxy[:, 0]))
    turn_angle = np.zeros(len(committed_angle))
    if len(committed_angle) > 1:
        turn_angle[1:] = _angle_wrap_deg(np.diff(committed_angle))
    direction_error = _angle_wrap_deg(committed_angle - requested_angle)

    rows["parameter_option"] = option
    rows["candidate_id"] = candidate
    rows["material_class"] = material_class
    rows["paper_role"] = role
    rows["temperature_K"] = temperature
    rows["hazard_seed"] = seed
    rows["event_index"] = np.arange(len(rows), dtype=int)
    rows["K_MPa_sqrt_m"] = rows["KJ_Pa_sqrtm"] / 1.0e6
    rows["J_effective_J_m2"] = root["J_effective_trial"].to_numpy(float)
    rows["J_signed_J_m2"] = root["J_signed_trial"].to_numpy(float)
    rows["J_effective_kJ_m2"] = rows["J_effective_J_m2"] / 1.0e3
    rows["J_signed_kJ_m2"] = rows["J_signed_J_m2"] / 1.0e3
    rows["crack_extension_um"] = rows["crack_extension_m"] * 1.0e6
    rows["event_advance_um"] = moved_um
    rows["threshold_action"] = thresholds
    rows["requested_direction_deg"] = requested_angle
    rows["committed_direction_deg"] = committed_angle
    rows["direction_error_deg"] = direction_error
    rows["turn_angle_deg"] = turn_angle
    rows["path_x_um"] = (xy[1:, 0] - xy[0, 0]) * 1.0e6
    rows["path_y_um"] = (xy[1:, 1] - xy[0, 1]) * 1.0e6
    rows["committed_segment_length_um"] = seg_len * 1.0e6

    path_table = path.copy()
    path_table["parameter_option"] = option
    path_table["candidate_id"] = candidate
    path_table["temperature_K"] = temperature
    path_table["path_point_index"] = np.arange(len(path_table), dtype=int)
    path_table["x_extension_um"] = (path_table["x_m"] - path_table["x_m"].iloc[0]) * 1.0e6
    path_table["y_offset_um"] = (path_table["y_m"] - path_table["y_m"].iloc[0]) * 1.0e6

    snapshot_matches = sorted(case.glob("mpz_state_snapshots_*K.json"))
    snapshot_payload = _load_json(snapshot_matches[0], {}) if snapshot_matches else {}
    snapshot_count = len(snapshot_payload.get("snapshots", []))
    visual_snapshots = [
        path_item for path_item in case.rglob("*.png")
        if "snapshot" in path_item.name.lower() or "state" in path_item.name.lower()
    ]
    final_extension = float(rows["crack_extension_um"].iloc[-1])
    complete = bool(manifest.get("run_completed_without_exception", False))
    mechanics_complete = final_extension + 1.0e-9 >= target_um
    net = xy[-1] - xy[0]
    total_path = float(np.sum(seg_len))
    net_length = float(np.linalg.norm(net))
    finite_err = np.abs(direction_error[np.isfinite(direction_error)])

    summary = {
        "parameter_option": option,
        "candidate_id": candidate,
        "material_class": material_class,
        "paper_role": role,
        "parameter_entry": transfer.get("parameter_entry") or manifest.get("parameter_entry"),
        "selected_row_sha256": transfer.get("selected_row_sha256"),
        "temperature_K": temperature,
        "hazard_seed": seed,
        "event_sequence_sha256": _fingerprint(accepted),
        "n_geometry_events": len(accepted),
        "final_extension_um": final_extension,
        "mean_event_advance_um": float(np.mean(moved_um)),
        "minimum_event_advance_um": float(np.min(moved_um)),
        "maximum_event_advance_um": float(np.max(moved_um)),
        "K_initial_MPa_sqrt_m": float(rows["K_MPa_sqrt_m"].iloc[0]),
        "K_mean_MPa_sqrt_m": float(rows["K_MPa_sqrt_m"].mean()),
        "K_max_MPa_sqrt_m": float(rows["K_MPa_sqrt_m"].max()),
        "K_end_MPa_sqrt_m": float(rows["K_MPa_sqrt_m"].iloc[-1]),
        "J_initial_kJ_m2": float(rows["J_effective_kJ_m2"].iloc[0]),
        "J_mean_kJ_m2": float(rows["J_effective_kJ_m2"].mean()),
        "J_max_kJ_m2": float(rows["J_effective_kJ_m2"].max()),
        "J_end_kJ_m2": float(rows["J_effective_kJ_m2"].iloc[-1]),
        "J_signed_initial_kJ_m2": float(rows["J_signed_kJ_m2"].iloc[0]),
        "J_signed_mean_kJ_m2": float(rows["J_signed_kJ_m2"].mean()),
        "J_signed_end_kJ_m2": float(rows["J_signed_kJ_m2"].iloc[-1]),
        "path_initial_angle_deg": float(committed_angle[0]),
        "path_final_angle_deg": float(committed_angle[-1]),
        "path_net_deflection_deg": float(np.degrees(np.arctan2(net[1], net[0]))),
        "path_max_abs_turn_deg": float(np.max(np.abs(turn_angle))) if len(turn_angle) else 0.0,
        "path_n_direction_changes_gt_1deg": int(np.count_nonzero(np.abs(turn_angle) > 1.0)),
        "path_y_span_um": float((np.max(xy[:, 1]) - np.min(xy[:, 1])) * 1.0e6),
        "path_total_length_um": total_path * 1.0e6,
        "path_net_length_um": net_length * 1.0e6,
        "path_tortuosity": total_path / max(net_length, 1.0e-300),
        "direction_error_mean_abs_deg": float(np.mean(finite_err)) if len(finite_err) else np.nan,
        "direction_error_max_abs_deg": float(np.max(finite_err)) if len(finite_err) else np.nan,
        "state_snapshot_records": int(snapshot_count),
        "visual_snapshot_files": int(len(visual_snapshots)),
        "manifest_status": manifest.get("status", "missing"),
        "run_completed_without_exception": complete,
        "mechanics_reached_target": mechanics_complete,
        "runtime_error": manifest.get("runtime_error"),
        "case_directory": str(case.resolve()),
    }
    return rows, path_table, summary


def _plot_rcurves(
    events: pd.DataFrame,
    out: Path,
    option: str,
    y_column: str,
    ylabel: str,
    stem: str,
) -> None:
    subset = events.loc[events["parameter_option"] == option]
    temperatures = sorted(subset["temperature_K"].unique())
    cmap = plt.get_cmap("jet")
    denom = max(len(temperatures) - 1, 1)
    fig, ax = plt.subplots(figsize=(8.2, 5.8))
    for index, temperature in enumerate(temperatures):
        rows = subset.loc[subset["temperature_K"] == temperature]
        ax.plot(
            rows["crack_extension_um"], rows[y_column],
            marker="o", markersize=4.0, linewidth=1.2,
            color=cmap(index / denom), label=f"{temperature:g} K",
        )
    ax.set_xlabel(r"Crack extension, $\Delta a$ ($\mu$m)")
    ax.set_ylabel(ylabel)
    ax.legend(ncol=2, fontsize=9, frameon=False)
    ax.tick_params(direction="in", top=True, right=True)
    _save(fig, out / option / stem)


def _plot_temperature(
    summary: pd.DataFrame,
    out: Path,
    option: str,
    prefix: str,
    columns: tuple[str, str, str],
    ylabel: str,
) -> None:
    rows = summary.loc[summary["parameter_option"] == option].sort_values("temperature_K")
    fig, ax = plt.subplots(figsize=(7.6, 5.4))
    for column, label, marker in zip(
        columns,
        ("Initial", "Mean", "Final"),
        ("o", "s", "D"),
    ):
        ax.plot(
            rows["temperature_K"], rows[column],
            marker=marker, linewidth=1.4, label=label,
        )
    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel(ylabel)
    ax.legend(frameon=False)
    ax.tick_params(direction="in", top=True, right=True)
    _save(fig, out / option / f"{prefix}_vs_temperature_initial_mean_final")


def _plot_metric_comparison(
    summary: pd.DataFrame,
    out: Path,
    column: str,
    ylabel: str,
    stem: str,
) -> None:
    fig, ax = plt.subplots(figsize=(8.2, 5.8))
    markers = ("o", "s", "^", "D", "v", "P", "X", "<")
    for index, option in enumerate(sorted(summary["parameter_option"].unique())):
        rows = summary.loc[summary["parameter_option"] == option].sort_values("temperature_K")
        ax.plot(
            rows["temperature_K"], rows[column],
            marker=markers[index % len(markers)], linewidth=1.35,
            label=_short_label(option),
        )
    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel(ylabel)
    ax.legend(fontsize=8.5, frameon=False, ncol=2)
    ax.tick_params(direction="in", top=True, right=True)
    _save(fig, out / stem)


def _plot_paths(paths: pd.DataFrame, out: Path, option: str) -> None:
    subset = paths.loc[paths["parameter_option"] == option]
    temperatures = sorted(subset["temperature_K"].unique())
    cmap = plt.get_cmap("jet")
    denom = max(len(temperatures) - 1, 1)
    fig, ax = plt.subplots(figsize=(7.2, 5.8))
    for index, temperature in enumerate(temperatures):
        rows = subset.loc[subset["temperature_K"] == temperature]
        ax.plot(
            rows["x_extension_um"], rows["y_offset_um"],
            marker="o", markersize=2.8, linewidth=1.1,
            color=cmap(index / denom), label=f"{temperature:g} K",
        )
    ax.set_xlabel(r"Projected crack extension, $\Delta x$ ($\mu$m)")
    ax.set_ylabel(r"Crack-path offset, $\Delta y$ ($\mu$m)")
    ax.legend(ncol=2, fontsize=9, frameon=False)
    ax.tick_params(direction="in", top=True, right=True)
    ax.set_aspect("equal", adjustable="datalim")
    _save(fig, out / option / "crack_paths_all_temperatures")


def _plot_direction_history(events: pd.DataFrame, out: Path, option: str) -> None:
    subset = events.loc[events["parameter_option"] == option]
    temperatures = sorted(subset["temperature_K"].unique())
    cmap = plt.get_cmap("jet")
    denom = max(len(temperatures) - 1, 1)
    fig, ax = plt.subplots(figsize=(8.2, 5.8))
    for index, temperature in enumerate(temperatures):
        rows = subset.loc[subset["temperature_K"] == temperature]
        ax.plot(
            rows["crack_extension_um"], rows["committed_direction_deg"],
            marker="o", markersize=3.0, linewidth=1.0,
            color=cmap(index / denom), label=f"{temperature:g} K",
        )
    ax.set_xlabel(r"Crack extension, $\Delta a$ ($\mu$m)")
    ax.set_ylabel("Committed segment direction (deg)")
    ax.legend(ncol=2, fontsize=9, frameon=False)
    ax.tick_params(direction="in", top=True, right=True)
    _save(fig, out / option / "committed_crack_direction_vs_extension")


def _expected_cases(campaign_root: Path) -> pd.DataFrame:
    path = campaign_root / "case_matrix.tsv"
    if not path.is_file():
        return pd.DataFrame(columns=["parameter_entry", "parameter_option", "temperature_K", "hazard_seed"])
    return pd.read_csv(
        path, sep="\t", header=None,
        names=["parameter_entry", "parameter_option", "temperature_K", "hazard_seed"],
    )


@contextmanager
def _file_lock(path: Path) -> Iterator[None]:
    import fcntl
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a+") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def analyze(
    campaign_root: Path,
    out: Path,
    target_um: float,
    *,
    allow_partial: bool = False,
) -> dict[str, Any]:
    option_dirs = sorted(
        path for path in campaign_root.iterdir()
        if path.is_dir() and path.name.startswith("v913_paper_")
    )
    expected = _expected_cases(campaign_root)
    event_tables: list[pd.DataFrame] = []
    path_tables: list[pd.DataFrame] = []
    summaries: list[dict[str, Any]] = []
    skipped: list[dict[str, Any]] = []

    for option_dir in option_dirs:
        for case in sorted(path for path in option_dir.glob("T????K") if path.is_dir()):
            try:
                events, path_table, summary = load_case(case, target_um)
            except Exception as exc:
                skipped.append({
                    "parameter_option": option_dir.name,
                    "temperature_directory": case.name,
                    "error_type": type(exc).__name__,
                    "error": str(exc),
                })
                if not allow_partial:
                    raise
                continue
            event_tables.append(events)
            path_tables.append(path_table)
            summaries.append(summary)

    out.mkdir(parents=True, exist_ok=True)
    skipped_table = pd.DataFrame(skipped)
    skipped_table.to_csv(out / "analysis_skipped_cases.csv", index=False)

    if not summaries:
        payload = {
            "schema": SCHEMA,
            "campaign_root": str(campaign_root.resolve()),
            "analysis_root": str(out.resolve()),
            "allow_partial": allow_partial,
            "n_expected_cases": int(len(expected)),
            "n_analyzed_cases": 0,
            "n_skipped_cases": len(skipped),
            "analysis_complete": False,
        }
        (out / "analysis_manifest.json").write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n"
        )
        if allow_partial:
            return payload
        raise RuntimeError("no complete event histories were available for strict analysis")

    events = pd.concat(event_tables, ignore_index=True)
    paths = pd.concat(path_tables, ignore_index=True)
    summary = pd.DataFrame(summaries).sort_values(
        ["parameter_option", "temperature_K"]
    ).reset_index(drop=True)

    if summary["hazard_seed"].isna().any():
        raise RuntimeError("one or more analyzed cases do not report an explicit hazard seed")
    if summary["hazard_seed"].duplicated().any():
        raise RuntimeError("hazard seeds are not globally unique across analyzed cases")
    if summary["event_sequence_sha256"].duplicated().any():
        raise RuntimeError("two analyzed cases have identical stochastic event sequences")

    expected_count = int(len(expected)) if len(expected) else len(summary)
    all_expected_loaded = len(summary) == expected_count and not skipped
    all_target = bool(summary["mechanics_reached_target"].all())
    if not allow_partial:
        if not all_expected_loaded:
            raise RuntimeError(
                f"strict analysis expected {expected_count} cases but loaded {len(summary)}; "
                f"skipped={len(skipped)}"
            )
        if not all_target:
            failed = summary.loc[
                ~summary["mechanics_reached_target"],
                ["parameter_option", "temperature_K", "final_extension_um"],
            ]
            raise RuntimeError(f"one or more cases did not reach the target:\n{failed}")

    events.to_csv(out / "R_curve_event_points_all_parameterizations.csv", index=False)
    paths.to_csv(out / "crack_path_points_all_parameterizations.csv", index=False)
    summary.to_csv(out / "K_J_temperature_path_seed_completion_summary.csv", index=False)
    summary.to_csv(out / "K_temperature_seed_and_completion_summary.csv", index=False)
    summary[[
        "parameter_option", "candidate_id", "material_class", "temperature_K",
        "hazard_seed", "event_sequence_sha256", "state_snapshot_records",
        "visual_snapshot_files", "selected_row_sha256",
        "path_n_direction_changes_gt_1deg", "path_max_abs_turn_deg",
    ]].to_csv(out / "parameter_temperature_seed_snapshot_path_audit.csv", index=False)

    options = sorted(summary["parameter_option"].unique())
    for option in options:
        _plot_rcurves(
            events, out, option, "K_MPa_sqrt_m",
            r"$K_J$ (MPa$\sqrt{\mathrm{m}}$)", "K_R_curves_all_temperatures",
        )
        _plot_rcurves(
            events, out, option, "J_effective_kJ_m2",
            r"Directional $J$ (kJ m$^{-2}$)", "J_R_curves_all_temperatures",
        )
        _plot_temperature(
            summary, out, option, "K",
            ("K_initial_MPa_sqrt_m", "K_mean_MPa_sqrt_m", "K_end_MPa_sqrt_m"),
            r"$K_J$ (MPa$\sqrt{\mathrm{m}}$)",
        )
        _plot_temperature(
            summary, out, option, "J",
            ("J_initial_kJ_m2", "J_mean_kJ_m2", "J_end_kJ_m2"),
            r"Directional $J$ (kJ m$^{-2}$)",
        )
        _plot_paths(paths, out, option)
        _plot_direction_history(events, out, option)

    metric_groups = (
        (
            "K",
            r"$K_J$ (MPa$\sqrt{\mathrm{m}}$)",
            (
                ("K_initial_MPa_sqrt_m", "initial"),
                ("K_mean_MPa_sqrt_m", "mean"),
                ("K_end_MPa_sqrt_m", "final"),
            ),
        ),
        (
            "J",
            r"Directional $J$ (kJ m$^{-2}$)",
            (
                ("J_initial_kJ_m2", "initial"),
                ("J_mean_kJ_m2", "mean"),
                ("J_end_kJ_m2", "final"),
            ),
        ),
    )
    for prefix, ylabel, metrics in metric_groups:
        for column, label in metrics:
            _plot_metric_comparison(
                summary, out, column,
                ylabel, f"all_parameterizations_{prefix}_{label}_vs_temperature",
            )

    payload = {
        "schema": SCHEMA,
        "campaign_root": str(campaign_root.resolve()),
        "analysis_root": str(out.resolve()),
        "target_extension_um": target_um,
        "allow_partial": allow_partial,
        "n_expected_cases": expected_count,
        "n_analyzed_cases": len(summary),
        "n_cases": len(summary),
        "n_skipped_cases": len(skipped),
        "n_parameterizations": len(options),
        "parameter_options": options,
        "all_seeds_unique": True,
        "all_event_sequences_unique": True,
        "all_expected_cases_loaded": all_expected_loaded,
        "all_analyzed_cases_reached_target": all_target,
        "analysis_complete": all_expected_loaded and all_target,
        "K_definition": "event-level K_J reported by the FEM domain-integral output",
        "J_definition": (
            "direct signed-positive directional domain-integral J_effective_trial "
            "for root-front accepted events; not reconstructed from K"
        ),
        "initial_definition": "first accepted physical cleavage event",
        "mean_definition": "arithmetic mean over accepted physical cleavage events",
        "final_definition": "last accepted physical cleavage event",
        "path_definition": "committed adaptive-CZM physical endpoints",
    }
    (out / "analysis_manifest.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n"
    )
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("campaign_root", type=Path)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--target-extension-um", type=float, default=200.0)
    parser.add_argument("--allow-partial", action="store_true")
    parser.add_argument("--lock-file", type=Path)
    args = parser.parse_args()
    campaign_root = args.campaign_root.expanduser().resolve()
    out = (
        args.out.expanduser().resolve()
        if args.out is not None
        else campaign_root / "final_analysis"
    )
    lock = _file_lock(args.lock_file.expanduser().resolve()) if args.lock_file else nullcontext()
    with lock:
        payload = analyze(
            campaign_root, out, args.target_extension_um,
            allow_partial=bool(args.allow_partial),
        )
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
