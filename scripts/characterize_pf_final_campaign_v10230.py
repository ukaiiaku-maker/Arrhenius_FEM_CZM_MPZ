#!/usr/bin/env python3
"""Characterize the PF final v10.2.30 four-class campaign from its manifest.

Reads PF_FINAL_CAMPAIGN_REFERENCE_MANIFEST.json (produced by
inventory_pf_final_campaign_v10230.py) and emits comparison-ready CSV
tables: first-passage vs. temperature, R-curve level summary, event
statistics, material-class ordering, and rate sensitivity. Does not touch
the PF repository -- reads only the already-extracted manifest.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path
from statistics import mean, pstdev

WORKSPACE_ROOT = Path(__file__).resolve().parents[1]
MANIFEST_PATH = WORKSPACE_ROOT / "PF_FINAL_CAMPAIGN_REFERENCE_MANIFEST.json"
OUT_DIR = WORKSPACE_ROOT / "runs" / "pf_final_campaign_v10230"

CLASS_ORDER = ["Peak", "DBTT", "weak-T", "ceramic"]
RATE_ORDER = ["rate0p01x", "rate1x", "rate100x"]


def _load_cases() -> list[dict]:
    with open(MANIFEST_PATH) as f:
        manifest = json.load(f)
    return manifest["cases"]


def _write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    print(f"wrote {path} ({len(rows)} rows)")


def first_passage_vs_temperature(cases: list[dict]) -> None:
    rows = []
    for c in sorted(cases, key=lambda c: (c["rate"], c.get("material_class", ""), c.get("temperature_K", 0))):
        rows.append({
            "rate": c["rate"],
            "material_class": c.get("material_class"),
            "option_dir": c.get("option_dir"),
            "temperature_K": c.get("temperature_K"),
            "seed": c.get("seed"),
            "first_passage_KJ_MPa_sqrt_m": c.get("first_passage_KJ_MPa_sqrt_m"),
            "first_passage_J_J_per_m2": c.get("first_passage_J_J_per_m2"),
            "first_passage_time_s": c.get("first_passage_time_s"),
            "first_passage_step": c.get("first_passage_step"),
        })
    _write_csv(
        OUT_DIR / "pf_first_passage_vs_temperature.csv",
        rows,
        ["rate", "material_class", "option_dir", "temperature_K", "seed",
         "first_passage_KJ_MPa_sqrt_m", "first_passage_J_J_per_m2",
         "first_passage_time_s", "first_passage_step"],
    )


def rcurve_summary(cases: list[dict]) -> None:
    rows = []
    for c in sorted(cases, key=lambda c: (c["rate"], c.get("material_class", ""), c.get("temperature_K", 0))):
        first_kj = c.get("first_passage_KJ_MPa_sqrt_m")
        final_kj = c.get("final_KJ_MPa_sqrt_m")
        rise = (final_kj - first_kj) if (first_kj is not None and final_kj is not None) else None
        ratio = (final_kj / first_kj) if (first_kj and final_kj is not None) else None
        rows.append({
            "rate": c["rate"],
            "material_class": c.get("material_class"),
            "temperature_K": c.get("temperature_K"),
            "seed": c.get("seed"),
            "first_passage_KJ_MPa_sqrt_m": first_kj,
            "final_KJ_MPa_sqrt_m": final_kj,
            "KJ_rise_MPa_sqrt_m": rise,
            "toughening_ratio_final_over_first": ratio,
            "final_crack_extension_m": c.get("final_crack_extension_m"),
            "projected_extension_um": c.get("projected_extension_um"),
            "mode": c.get("mode"),
            "branched": c.get("branched"),
        })
    _write_csv(
        OUT_DIR / "pf_rcurve_summary.csv",
        rows,
        ["rate", "material_class", "temperature_K", "seed",
         "first_passage_KJ_MPa_sqrt_m", "final_KJ_MPa_sqrt_m",
         "KJ_rise_MPa_sqrt_m", "toughening_ratio_final_over_first",
         "final_crack_extension_m", "projected_extension_um", "mode", "branched"],
    )


def event_statistics(cases: list[dict]) -> None:
    rows = []
    for c in sorted(cases, key=lambda c: (c["rate"], c.get("material_class", ""), c.get("temperature_K", 0))):
        rows.append({
            "rate": c["rate"],
            "material_class": c.get("material_class"),
            "temperature_K": c.get("temperature_K"),
            "seed": c.get("seed"),
            "n_geometry_events": c.get("n_geometry_events"),
            "n_advances": c.get("n_advances"),
            "n_advances_primary": c.get("n_advances_primary"),
            "n_advances_branch": c.get("n_advances_branch"),
            "event_length_count": c.get("event_length_count"),
            "event_length_mean_m": c.get("event_length_mean_m"),
            "event_length_min_m": c.get("event_length_min_m"),
            "event_length_max_m": c.get("event_length_max_m"),
            "N_em_init": c.get("N_em_init"),
            "N_em_final": c.get("N_em_final"),
            "W_emit_J_per_m": c.get("W_emit_J_per_m"),
        })
    _write_csv(
        OUT_DIR / "pf_event_statistics.csv",
        rows,
        ["rate", "material_class", "temperature_K", "seed",
         "n_geometry_events", "n_advances", "n_advances_primary", "n_advances_branch",
         "event_length_count", "event_length_mean_m", "event_length_min_m", "event_length_max_m",
         "N_em_init", "N_em_final", "W_emit_J_per_m"],
    )


def material_class_ordering(cases: list[dict]) -> None:
    """At each (rate, temperature) present for ALL four classes, rank
    classes by first-passage KJ to check ordering consistency."""
    by_rate_temp: dict[tuple[str, float], dict[str, float]] = {}
    for c in cases:
        key = (c["rate"], c.get("temperature_K"))
        cls = c.get("material_class")
        kj = c.get("first_passage_KJ_MPa_sqrt_m")
        if cls is None or kj is None:
            continue
        by_rate_temp.setdefault(key, {})[cls] = kj

    rows = []
    for (rate, temp), class_kj in sorted(by_rate_temp.items(), key=lambda kv: (kv[0][0], kv[0][1])):
        if set(class_kj.keys()) != set(CLASS_ORDER):
            continue
        ranked = sorted(class_kj.items(), key=lambda kv: kv[1])
        rows.append({
            "rate": rate,
            "temperature_K": temp,
            **{f"KJ_{cls}_MPa_sqrt_m": class_kj[cls] for cls in CLASS_ORDER},
            "ranking_low_to_high": ">".join(cls for cls, _ in ranked),
        })
    _write_csv(
        OUT_DIR / "pf_material_class_ordering.csv",
        rows,
        ["rate", "temperature_K"] + [f"KJ_{cls}_MPa_sqrt_m" for cls in CLASS_ORDER] + ["ranking_low_to_high"],
    )


def rate_sensitivity_summary(cases: list[dict]) -> None:
    by_class_temp: dict[tuple[str, float], dict[str, float]] = {}
    for c in cases:
        key = (c.get("material_class"), c.get("temperature_K"))
        kj = c.get("first_passage_KJ_MPa_sqrt_m")
        if kj is None:
            continue
        by_class_temp.setdefault(key, {})[c["rate"]] = kj

    rows = []
    for (cls, temp), rate_kj in sorted(by_class_temp.items(), key=lambda kv: (str(kv[0][0]), kv[0][1] or 0)):
        row = {"material_class": cls, "temperature_K": temp}
        for r in RATE_ORDER:
            row[f"KJ_{r}_MPa_sqrt_m"] = rate_kj.get(r)
        vals = [rate_kj[r] for r in RATE_ORDER if r in rate_kj]
        if len(vals) == len(RATE_ORDER):
            row["max_over_min_ratio"] = max(vals) / min(vals) if min(vals) else None
        rows.append(row)
    _write_csv(
        OUT_DIR / "pf_rate_sensitivity_summary.csv",
        rows,
        ["material_class", "temperature_K"] + [f"KJ_{r}_MPa_sqrt_m" for r in RATE_ORDER] + ["max_over_min_ratio"],
    )


def main() -> None:
    cases = _load_cases()
    first_passage_vs_temperature(cases)
    rcurve_summary(cases)
    event_statistics(cases)
    material_class_ordering(cases)
    rate_sensitivity_summary(cases)


if __name__ == "__main__":
    main()
