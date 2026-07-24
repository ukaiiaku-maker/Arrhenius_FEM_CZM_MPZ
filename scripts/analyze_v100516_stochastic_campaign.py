#!/usr/bin/env python3
"""Analyze and plot a v10.0.5.16 stochastic FEM/CZM temperature campaign.

The script consumes one campaign directory containing T####K case folders.  It
extracts only accepted cleavage-event rows (``n_fire > 0``) for the R-curves,
checks that every temperature used a unique stochastic seed and event sequence,
and writes publication-oriented figures plus auditable CSV/JSON summaries.

No simulation state is changed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

MANIFEST_NAME = "persistent_site_production_manifest_v10_0_5_16.json"
EVENTS_NAME = "stochastic_geometry_events_v10_0_5_16.json"
ANALYSIS_SCHEMA = "v10.0.5.16.1_stochastic_campaign_rcurve_analysis"


def _temperature_from_dir(path: Path) -> int:
    name = path.name
    if len(name) != 6 or not name.startswith("T") or not name.endswith("K"):
        raise ValueError(f"not a temperature case directory: {path}")
    return int(name[1:5])


def _load_json(path: Path, default: Any = None) -> Any:
    if not path.is_file():
        return default
    return json.loads(path.read_text())


def _first_existing(case: Path, pattern: str) -> Path | None:
    matches = sorted(case.glob(pattern))
    return matches[0] if matches else None


def _manifest_seed(manifest: dict[str, Any], events: list[dict[str, Any]]) -> int | None:
    physics = dict(manifest.get("physics_contract", {}) or {})
    seed = physics.get("cleavage_hazard_seed")
    if seed is None:
        front = dict(manifest.get("front_engine", {}) or {})
        stochastic = dict(front.get("stochastic_hazard", {}) or {})
        seed = stochastic.get("seed")
    if seed is None and events:
        seed = events[0].get("hazard_seed")
    return None if seed is None else int(seed)


def _accepted_events(events: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        dict(row)
        for row in events
        if bool(row.get("inserted", False)) and float(row.get("moved_m", 0.0)) > 0.0
    ]


def _sequence_fingerprint(events: list[dict[str, Any]]) -> str:
    payload = [
        {
            "event_index": int(row.get("hazard_event_index", index)),
            "threshold_action": float(row.get("threshold_action", 0.0)),
            "moved_m": float(row.get("moved_m", 0.0)),
        }
        for index, row in enumerate(events)
    ]
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def _case_status(manifest: dict[str, Any], extension_um: float, target_um: float | None) -> str:
    if bool(manifest.get("run_completed_without_exception", False)):
        return "complete"
    if target_um is not None and extension_um + 1.0e-9 >= target_um:
        return "mechanics_complete_postrun_audit_failed"
    return str(manifest.get("status", "unknown"))


def load_case(case: Path, target_um: float | None = None) -> tuple[pd.DataFrame, dict[str, Any]]:
    temperature = _temperature_from_dir(case)
    steps_path = _first_existing(case, "steps_*K.csv")
    if steps_path is None:
        raise FileNotFoundError(f"missing steps CSV in {case}")
    steps = pd.read_csv(steps_path)
    required = {"n_fire", "KJ_Pa_sqrtm", "crack_extension_m"}
    missing = required.difference(steps.columns)
    if missing:
        raise ValueError(f"{steps_path} is missing columns: {sorted(missing)}")

    event_rows = steps.loc[steps["n_fire"] > 0].copy()
    if event_rows.empty:
        raise ValueError(f"no accepted cleavage events in {steps_path}")
    event_rows = event_rows.reset_index(drop=True)
    event_rows["temperature_K"] = temperature
    event_rows["event_index"] = np.arange(len(event_rows), dtype=int)
    event_rows["K_MPa_sqrt_m"] = event_rows["KJ_Pa_sqrtm"] / 1.0e6
    event_rows["crack_extension_um"] = event_rows["crack_extension_m"] * 1.0e6

    manifest = _load_json(case / MANIFEST_NAME, {}) or {}
    raw_events = _load_json(case / EVENTS_NAME, []) or []
    accepted = _accepted_events(raw_events)
    if len(accepted) != len(event_rows):
        raise ValueError(
            f"{case.name}: accepted geometry events ({len(accepted)}) do not match "
            f"n_fire rows ({len(event_rows)})"
        )

    moved_um = np.asarray([float(row["moved_m"]) * 1.0e6 for row in accepted])
    thresholds = np.asarray(
        [float(row.get("threshold_action", np.nan)) for row in accepted]
    )
    event_rows["event_advance_um"] = moved_um
    event_rows["threshold_action"] = thresholds

    seed = _manifest_seed(manifest, accepted)
    event_rows["hazard_seed"] = seed
    extension_um = float(event_rows["crack_extension_um"].iloc[-1])

    snapshot_path = _first_existing(case, "mpz_state_snapshots_*K.json")
    snapshot_payload = _load_json(snapshot_path, {}) if snapshot_path else {}
    snapshot_records = len((snapshot_payload or {}).get("snapshots", []))
    visual_snapshot_files = sorted(
        path
        for path in case.rglob("*.png")
        if "snapshot" in path.name.lower() or "state" in path.name.lower()
    )

    summary = {
        "temperature_K": temperature,
        "hazard_seed": seed,
        "event_sequence_sha256": _sequence_fingerprint(accepted),
        "n_geometry_events": len(accepted),
        "final_extension_um": extension_um,
        "mean_event_advance_um": float(np.mean(moved_um)),
        "minimum_event_advance_um": float(np.min(moved_um)),
        "maximum_event_advance_um": float(np.max(moved_um)),
        "fraction_at_minimum_event_advance": float(
            np.mean(np.isclose(moved_um, np.min(moved_um), rtol=1.0e-10, atol=1.0e-12))
        ),
        "K_initial_MPa_sqrt_m": float(event_rows["K_MPa_sqrt_m"].iloc[0]),
        "K_mean_MPa_sqrt_m": float(event_rows["K_MPa_sqrt_m"].mean()),
        "K_max_MPa_sqrt_m": float(event_rows["K_MPa_sqrt_m"].max()),
        "K_end_MPa_sqrt_m": float(event_rows["K_MPa_sqrt_m"].iloc[-1]),
        "K_rise_end_minus_initial_MPa_sqrt_m": float(
            event_rows["K_MPa_sqrt_m"].iloc[-1]
            - event_rows["K_MPa_sqrt_m"].iloc[0]
        ),
        "state_snapshot_records": int(snapshot_records),
        "visual_snapshot_files": int(len(visual_snapshot_files)),
        "manifest_status": str(manifest.get("status", "missing")),
        "run_completed_without_exception": bool(
            manifest.get("run_completed_without_exception", False)
        ),
        "runtime_error_type": manifest.get("runtime_error_type"),
        "runtime_error": manifest.get("runtime_error"),
        "analysis_status": _case_status(manifest, extension_um, target_um),
        "case_directory": str(case.resolve()),
    }
    return event_rows, summary


def _save_figure(fig: plt.Figure, stem: Path) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_overlay(events: pd.DataFrame, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.2, 5.8))
    temperatures = sorted(events["temperature_K"].unique())
    cmap = plt.get_cmap("jet")
    denom = max(len(temperatures) - 1, 1)
    for index, temperature in enumerate(temperatures):
        rows = events.loc[events["temperature_K"] == temperature]
        ax.plot(
            rows["crack_extension_um"],
            rows["K_MPa_sqrt_m"],
            marker="o",
            markersize=4.2,
            linewidth=1.25,
            color=cmap(index / denom),
            label=f"{temperature:g} K",
        )
    ax.set_xlabel(r"Crack extension, $\Delta a$ ($\mu$m)")
    ax.set_ylabel(r"$K_J$ (MPa$\sqrt{\mathrm{m}}$)")
    ax.legend(ncol=2, fontsize=9, frameon=False)
    ax.tick_params(direction="in", top=True, right=True)
    _save_figure(fig, out / "R_curves_all_temperatures")


def plot_individual(events: pd.DataFrame, out: Path) -> None:
    target = out / "individual_R_curves"
    for temperature in sorted(events["temperature_K"].unique()):
        rows = events.loc[events["temperature_K"] == temperature]
        fig, ax = plt.subplots(figsize=(6.5, 4.8))
        ax.plot(
            rows["crack_extension_um"],
            rows["K_MPa_sqrt_m"],
            marker="o",
            markersize=4.8,
            linewidth=1.4,
        )
        ax.set_xlabel(r"Crack extension, $\Delta a$ ($\mu$m)")
        ax.set_ylabel(r"$K_J$ (MPa$\sqrt{\mathrm{m}}$)")
        ax.text(0.04, 0.94, f"{temperature:g} K", transform=ax.transAxes, va="top")
        ax.tick_params(direction="in", top=True, right=True)
        _save_figure(fig, target / f"R_curve_{int(temperature):04d}K")


def plot_temperature_summary(summary: pd.DataFrame, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.6, 5.4))
    styles = [
        ("K_initial_MPa_sqrt_m", "Initial", "o"),
        ("K_mean_MPa_sqrt_m", "Mean", "s"),
        ("K_max_MPa_sqrt_m", "Maximum", "^"),
        ("K_end_MPa_sqrt_m", "Endpoint", "D"),
    ]
    for column, label, marker in styles:
        ax.plot(
            summary["temperature_K"],
            summary[column],
            marker=marker,
            markersize=6,
            linewidth=1.4,
            label=label,
        )
    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel(r"$K_J$ (MPa$\sqrt{\mathrm{m}}$)")
    ax.legend(frameon=False)
    ax.tick_params(direction="in", top=True, right=True)
    _save_figure(fig, out / "K_vs_temperature_initial_mean_max_end")


def plot_event_lengths(events: pd.DataFrame, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.2, 5.8))
    temperatures = sorted(events["temperature_K"].unique())
    cmap = plt.get_cmap("jet")
    denom = max(len(temperatures) - 1, 1)
    for index, temperature in enumerate(temperatures):
        rows = events.loc[events["temperature_K"] == temperature]
        ax.plot(
            rows["crack_extension_um"],
            rows["event_advance_um"],
            marker="o",
            markersize=3.8,
            linewidth=1.0,
            color=cmap(index / denom),
            label=f"{temperature:g} K",
        )
    ax.set_xlabel(r"Crack extension after event ($\mu$m)")
    ax.set_ylabel(r"Accepted event advance ($\mu$m)")
    ax.legend(ncol=2, fontsize=9, frameon=False)
    ax.tick_params(direction="in", top=True, right=True)
    _save_figure(fig, out / "event_advance_vs_crack_extension")


def analyze_campaign(campaign_root: Path, out_dir: Path, target_um: float | None) -> dict[str, Any]:
    cases = sorted(
        path for path in campaign_root.glob("T????K") if path.is_dir()
    )
    if not cases:
        raise FileNotFoundError(f"no T####K directories found in {campaign_root}")

    all_events: list[pd.DataFrame] = []
    summaries: list[dict[str, Any]] = []
    for case in cases:
        events, summary = load_case(case, target_um=target_um)
        all_events.append(events)
        summaries.append(summary)

    events_table = pd.concat(all_events, ignore_index=True)
    summary_table = pd.DataFrame(summaries).sort_values("temperature_K").reset_index(drop=True)

    seeds = summary_table["hazard_seed"]
    if seeds.isna().any():
        raise RuntimeError("one or more cases do not report an explicit stochastic seed")
    if seeds.duplicated().any():
        duplicate_rows = summary_table.loc[seeds.duplicated(keep=False), ["temperature_K", "hazard_seed"]]
        raise RuntimeError(f"duplicate temperature seeds detected:\n{duplicate_rows}")
    fingerprints = summary_table["event_sequence_sha256"]
    if fingerprints.duplicated().any():
        duplicate_rows = summary_table.loc[
            fingerprints.duplicated(keep=False),
            ["temperature_K", "hazard_seed", "event_sequence_sha256"],
        ]
        raise RuntimeError(f"duplicate stochastic event sequences detected:\n{duplicate_rows}")

    out_dir.mkdir(parents=True, exist_ok=True)
    events_table.to_csv(out_dir / "R_curve_event_points_all_temperatures.csv", index=False)
    summary_table.to_csv(out_dir / "K_temperature_and_seed_summary.csv", index=False)
    summary_table[["temperature_K", "hazard_seed", "event_sequence_sha256"]].to_csv(
        out_dir / "temperature_seed_map.csv", index=False
    )

    plot_overlay(events_table, out_dir)
    plot_individual(events_table, out_dir)
    plot_temperature_summary(summary_table, out_dir)
    plot_event_lengths(events_table, out_dir)

    audit = {
        "schema": ANALYSIS_SCHEMA,
        "campaign_root": str(campaign_root.resolve()),
        "output_directory": str(out_dir.resolve()),
        "temperatures_K": summary_table["temperature_K"].astype(int).tolist(),
        "unique_temperature_seeds": True,
        "unique_event_sequence_fingerprints": True,
        "event_rows_match_accepted_geometry_events": True,
        "metrics": {
            "initial": "first accepted cleavage event K_J",
            "mean": "arithmetic mean K_J over accepted cleavage events",
            "maximum": "maximum K_J over accepted cleavage events",
            "endpoint": "K_J at final accepted cleavage event",
        },
        "n_cases": int(len(summary_table)),
        "n_total_events": int(len(events_table)),
        "all_cases_reached_target": (
            None
            if target_um is None
            else bool((summary_table["final_extension_um"] + 1.0e-9 >= target_um).all())
        ),
        "state_snapshot_records_by_temperature": {
            str(int(row.temperature_K)): int(row.state_snapshot_records)
            for row in summary_table.itertuples()
        },
        "visual_snapshot_files_by_temperature": {
            str(int(row.temperature_K)): int(row.visual_snapshot_files)
            for row in summary_table.itertuples()
        },
    }
    (out_dir / "analysis_audit.json").write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n"
    )
    return audit


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign-root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, default=None)
    parser.add_argument("--target-extension-um", type=float, default=200.0)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    root = args.campaign_root.expanduser().resolve()
    out = (
        args.out_dir.expanduser().resolve()
        if args.out_dir is not None
        else root / "final_analysis"
    )
    audit = analyze_campaign(root, out, args.target_extension_um)
    print(json.dumps(audit, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
