#!/usr/bin/env python3
"""Analyze the v10.0.5.17 audited paper-parameter stochastic campaign."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

MANIFEST_NAME = "persistent_site_production_manifest_v10_0_5_17.json"
TRANSFER_NAME = "audited_PF_paper_parameter_transfer_v10_0_5_17.json"
EVENTS_NAME = "stochastic_geometry_events_v10_0_5_16.json"
SCHEMA = "v10.0.5.17_audited_paper_parameter_campaign_analysis"


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
    text = option.replace("v913_paper_", "").replace("_persistent_sites", "")
    return text.replace("weakT", "weak-T").replace("ceramic", "ceramic ")


def load_case(case: Path, target_um: float) -> tuple[pd.DataFrame, dict[str, Any]]:
    temperature = _temperature(case)
    steps = pd.read_csv(_first(case, "steps_*K.csv"))
    required = {"n_fire", "KJ_Pa_sqrtm", "crack_extension_m"}
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
    option = str(transfer.get("parameter_option") or manifest.get("parameter_option") or case.parent.name)
    candidate = str(transfer.get("candidate_id") or manifest.get("candidate_id") or "unknown")
    material_class = str(transfer.get("material_class") or manifest.get("material_class") or "unknown")
    role = str(transfer.get("paper_role") or manifest.get("paper_role") or option)
    seed = _seed(manifest, accepted)

    moved_um = np.asarray([float(row["moved_m"]) * 1.0e6 for row in accepted])
    thresholds = np.asarray([float(row.get("threshold_action", np.nan)) for row in accepted])
    rows["parameter_option"] = option
    rows["candidate_id"] = candidate
    rows["material_class"] = material_class
    rows["paper_role"] = role
    rows["temperature_K"] = temperature
    rows["hazard_seed"] = seed
    rows["event_index"] = np.arange(len(rows), dtype=int)
    rows["K_MPa_sqrt_m"] = rows["KJ_Pa_sqrtm"] / 1.0e6
    rows["crack_extension_um"] = rows["crack_extension_m"] * 1.0e6
    rows["event_advance_um"] = moved_um
    rows["threshold_action"] = thresholds

    snapshot_path = _first(case, "mpz_state_snapshots_*K.json")
    snapshot_payload = _load_json(snapshot_path, {})
    snapshot_count = len(snapshot_payload.get("snapshots", []))
    visual_snapshots = [
        path for path in case.rglob("*.png")
        if "snapshot" in path.name.lower() or "state" in path.name.lower()
    ]
    final_extension = float(rows["crack_extension_um"].iloc[-1])
    complete = bool(manifest.get("run_completed_without_exception", False))
    mechanics_complete = final_extension + 1.0e-9 >= target_um
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
        "K_rise_end_minus_initial_MPa_sqrt_m": float(
            rows["K_MPa_sqrt_m"].iloc[-1] - rows["K_MPa_sqrt_m"].iloc[0]
        ),
        "state_snapshot_records": int(snapshot_count),
        "visual_snapshot_files": int(len(visual_snapshots)),
        "manifest_status": manifest.get("status", "missing"),
        "run_completed_without_exception": complete,
        "mechanics_reached_target": mechanics_complete,
        "runtime_error": manifest.get("runtime_error"),
        "case_directory": str(case.resolve()),
    }
    return rows, summary


def plot_parameter_rcurves(events: pd.DataFrame, out: Path, option: str) -> None:
    subset = events.loc[events["parameter_option"] == option]
    temperatures = sorted(subset["temperature_K"].unique())
    cmap = plt.get_cmap("jet")
    denom = max(len(temperatures) - 1, 1)
    fig, ax = plt.subplots(figsize=(8.2, 5.8))
    for index, temperature in enumerate(temperatures):
        rows = subset.loc[subset["temperature_K"] == temperature]
        ax.plot(
            rows["crack_extension_um"],
            rows["K_MPa_sqrt_m"],
            marker="o",
            markersize=4.0,
            linewidth=1.2,
            color=cmap(index / denom),
            label=f"{temperature:g} K",
        )
    ax.set_xlabel(r"Crack extension, $\Delta a$ ($\mu$m)")
    ax.set_ylabel(r"$K_J$ (MPa$\sqrt{\mathrm{m}}$)")
    ax.legend(ncol=2, fontsize=9, frameon=False)
    ax.tick_params(direction="in", top=True, right=True)
    _save(fig, out / option / "R_curves_all_temperatures")


def plot_parameter_temperature(summary: pd.DataFrame, out: Path, option: str) -> None:
    rows = summary.loc[summary["parameter_option"] == option].sort_values("temperature_K")
    fig, ax = plt.subplots(figsize=(7.6, 5.4))
    for column, label, marker in (
        ("K_initial_MPa_sqrt_m", "Initial", "o"),
        ("K_mean_MPa_sqrt_m", "Mean", "s"),
        ("K_max_MPa_sqrt_m", "Maximum", "^"),
        ("K_end_MPa_sqrt_m", "Endpoint", "D"),
    ):
        ax.plot(rows["temperature_K"], rows[column], marker=marker, linewidth=1.4, label=label)
    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel(r"$K_J$ (MPa$\sqrt{\mathrm{m}}$)")
    ax.legend(frameon=False)
    ax.tick_params(direction="in", top=True, right=True)
    _save(fig, out / option / "K_vs_temperature_initial_mean_max_end")


def plot_metric_comparison(summary: pd.DataFrame, out: Path, column: str, stem: str) -> None:
    fig, ax = plt.subplots(figsize=(8.2, 5.8))
    markers = ("o", "s", "^", "D", "v", "P", "X", "<")
    for index, option in enumerate(sorted(summary["parameter_option"].unique())):
        rows = summary.loc[summary["parameter_option"] == option].sort_values("temperature_K")
        ax.plot(
            rows["temperature_K"],
            rows[column],
            marker=markers[index % len(markers)],
            linewidth=1.35,
            label=_short_label(option),
        )
    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel(r"$K_J$ (MPa$\sqrt{\mathrm{m}}$)")
    ax.legend(fontsize=8.5, frameon=False, ncol=2)
    ax.tick_params(direction="in", top=True, right=True)
    _save(fig, out / stem)


def analyze(campaign_root: Path, out: Path, target_um: float) -> dict[str, Any]:
    option_dirs = sorted(path for path in campaign_root.iterdir() if path.is_dir() and path.name.startswith("v913_paper_"))
    if not option_dirs:
        raise FileNotFoundError(f"no parameter-option directories found in {campaign_root}")
    event_tables: list[pd.DataFrame] = []
    summaries: list[dict[str, Any]] = []
    for option_dir in option_dirs:
        cases = sorted(path for path in option_dir.glob("T????K") if path.is_dir())
        if not cases:
            raise FileNotFoundError(f"no temperature cases in {option_dir}")
        for case in cases:
            events, summary = load_case(case, target_um)
            event_tables.append(events)
            summaries.append(summary)
    events = pd.concat(event_tables, ignore_index=True)
    summary = pd.DataFrame(summaries).sort_values(["parameter_option", "temperature_K"]).reset_index(drop=True)

    if summary["hazard_seed"].isna().any():
        raise RuntimeError("one or more cases do not report an explicit hazard seed")
    if summary["hazard_seed"].duplicated().any():
        raise RuntimeError("hazard seeds are not globally unique across option-temperature cases")
    if summary["event_sequence_sha256"].duplicated().any():
        raise RuntimeError("two option-temperature cases have identical stochastic event sequences")
    if not summary["mechanics_reached_target"].all():
        failed = summary.loc[~summary["mechanics_reached_target"], ["parameter_option", "temperature_K", "final_extension_um"]]
        raise RuntimeError(f"one or more cases did not reach the target:\n{failed}")

    out.mkdir(parents=True, exist_ok=True)
    events.to_csv(out / "R_curve_event_points_all_parameterizations.csv", index=False)
    summary.to_csv(out / "K_temperature_seed_and_completion_summary.csv", index=False)
    summary[[
        "parameter_option", "candidate_id", "material_class", "temperature_K",
        "hazard_seed", "event_sequence_sha256", "state_snapshot_records",
        "visual_snapshot_files", "selected_row_sha256",
    ]].to_csv(out / "parameter_temperature_seed_snapshot_audit.csv", index=False)

    options = sorted(summary["parameter_option"].unique())
    for option in options:
        plot_parameter_rcurves(events, out, option)
        plot_parameter_temperature(summary, out, option)
    for column, stem in (
        ("K_initial_MPa_sqrt_m", "all_parameterizations_K_initial_vs_temperature"),
        ("K_mean_MPa_sqrt_m", "all_parameterizations_K_mean_vs_temperature"),
        ("K_max_MPa_sqrt_m", "all_parameterizations_K_max_vs_temperature"),
        ("K_end_MPa_sqrt_m", "all_parameterizations_K_end_vs_temperature"),
    ):
        plot_metric_comparison(summary, out, column, stem)

    payload = {
        "schema": SCHEMA,
        "campaign_root": str(campaign_root.resolve()),
        "analysis_root": str(out.resolve()),
        "target_extension_um": target_um,
        "n_parameterizations": len(options),
        "n_cases": len(summary),
        "parameter_options": options,
        "all_seeds_unique": True,
        "all_event_sequences_unique": True,
        "all_cases_reached_target": True,
        "K_mean_definition": "arithmetic mean over accepted cleavage-event K_J values",
        "K_initial_definition": "K_J at first accepted event",
        "K_max_definition": "maximum accepted-event K_J",
        "K_end_definition": "K_J at final accepted event",
    }
    (out / "analysis_manifest.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("campaign_root", type=Path)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--target-extension-um", type=float, default=200.0)
    args = parser.parse_args()
    root = args.campaign_root.expanduser().resolve()
    out = (args.out or root / "final_analysis").expanduser().resolve()
    payload = analyze(root, out, args.target_extension_um)
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
