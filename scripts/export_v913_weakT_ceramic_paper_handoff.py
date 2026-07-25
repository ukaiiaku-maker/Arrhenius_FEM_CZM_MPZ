#!/usr/bin/env python3
"""Export the selected weak-T/FCC-like and ceramic-like rows for 2-D transfer.

The exporter reads the completed local selection output, preserves the exact 29
active v9.13 candidate parameters, and writes a compact CSV plus a hash-checked
JSON manifest.  No parameter is fitted, transformed, rounded, or defaulted during
this handoff.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import pandas as pd

from arrhenius_fracture.emergent_gnd_contract_v913 import (
    ACTIVE_CANDIDATE_PARAMETER_FIELDS,
    candidate_parameter_fingerprint,
    effective_candidate_parameters,
)


EXPECTED = {
    "weakT_FCC_like": {
        "candidate_id": "v913_zeroD_sobol_0257068",
        "option_key": "v913_paper_weakT01_0257068_persistent_sites",
        "selection_role": "paper primary weak-temperature/FCC-like candidate",
    },
    "ceramic_like": {
        "candidate_id": "v913_zeroD_sobol_0189364",
        "option_key": "v913_paper_ceramic01_0189364_persistent_sites",
        "selection_role": "paper primary ceramic-like candidate",
    },
}

METRIC_FIELDS = (
    "complete_temperature_grid",
    "K50_mean_MPa_sqrt_m",
    "K50_median_MPa_sqrt_m",
    "K50_temperature_span_MPa_sqrt_m",
    "K50_low_mid_span_700_1100_MPa_sqrt_m",
    "K50_low_mid_median_MPa_sqrt_m",
    "K50_high_median_1200_1400_MPa_sqrt_m",
    "high_temperature_toughness_loss_MPa_sqrt_m",
    "pre_high_temperature_change_K1100_minus_K700_MPa_sqrt_m",
    "terminal_change_K1400_minus_K1100_MPa_sqrt_m",
    "high_temperature_rebound_MPa_sqrt_m",
    "median_R_rise_first_to_50_MPa_sqrt_m",
    "median_abs_R_rise_first_to_50_MPa_sqrt_m",
    "median_R_rise_25_to_50_MPa_sqrt_m",
    "median_abs_R_rise_25_to_50_MPa_sqrt_m",
    "positive_R_rise_temperature_fraction",
    "K50_peak_temperature_K",
    "K50_peak_prominence_MPa_sqrt_m",
    "weakT_gate",
    "weakT_score",
    "ceramic_gate",
    "ceramic_score",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selection-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    return parser.parse_args()


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, bool):
        return value
    if hasattr(value, "item"):
        return json_safe(value.item())
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    return value


def _selected_metrics(summary: dict[str, Any], material_class: str) -> dict[str, Any]:
    selected = summary.get("selected", {})
    row = selected.get(material_class)
    if not isinstance(row, dict):
        raise RuntimeError(f"selection summary lacks selected.{material_class}")
    return row


def build_handoff(
    registry: pd.DataFrame,
    summary: dict[str, Any],
) -> tuple[pd.DataFrame, list[dict[str, Any]]]:
    if registry["candidate_id"].astype(str).duplicated().any():
        raise RuntimeError("selected registry contains duplicate candidate IDs")
    source = registry.set_index(registry["candidate_id"].astype(str), drop=False)

    rows: list[dict[str, Any]] = []
    selected_metadata: list[dict[str, Any]] = []
    for material_class, expected in EXPECTED.items():
        candidate_id = expected["candidate_id"]
        if candidate_id not in source.index:
            raise RuntimeError(
                f"selected registry does not contain required {material_class} row: {candidate_id}"
            )
        source_row = source.loc[candidate_id].to_dict()
        active = effective_candidate_parameters(source_row)
        metrics = _selected_metrics(summary, material_class)
        if str(metrics.get("candidate_id")) != candidate_id:
            raise RuntimeError(
                f"selection summary {material_class} candidate mismatch: "
                f"expected={candidate_id}, observed={metrics.get('candidate_id')}"
            )

        strict_field = "weakT_gate" if material_class == "weakT_FCC_like" else "ceramic_gate"
        score_field = "weakT_score" if material_class == "weakT_FCC_like" else "ceramic_score"
        strict = bool(metrics[strict_field])
        score = float(metrics[score_field])
        record = {
            "option_key": expected["option_key"],
            "candidate_id": candidate_id,
            "paper_material_class": material_class,
            "selection_role": expected["selection_role"],
            "oneD_strict_gate_passed": strict,
            "oneD_selection_score": score,
            **active,
        }
        rows.append(record)

        selected_metadata.append(
            {
                "option_key": expected["option_key"],
                "candidate_id": candidate_id,
                "paper_material_class": material_class,
                "selection_role": expected["selection_role"],
                "oneD_strict_gate_passed": strict,
                "oneD_selection_score": score,
                "oneD_metrics": {
                    key: metrics.get(key)
                    for key in METRIC_FIELDS
                    if key in metrics
                },
            }
        )

    handoff = pd.DataFrame(rows)
    expected_columns = [
        "option_key",
        "candidate_id",
        "paper_material_class",
        "selection_role",
        "oneD_strict_gate_passed",
        "oneD_selection_score",
        *ACTIVE_CANDIDATE_PARAMETER_FIELDS,
    ]
    handoff = handoff[expected_columns]
    return handoff, selected_metadata


def main() -> int:
    args = parse_args()
    selection_dir = args.selection_dir.expanduser().resolve()
    source_registry = selection_dir / "selected_weakT_ceramic_registry.csv"
    source_summary = selection_dir / "selection_summary.json"
    if not source_registry.is_file():
        raise FileNotFoundError(source_registry)
    if not source_summary.is_file():
        raise FileNotFoundError(source_summary)

    registry = pd.read_csv(source_registry)
    summary = json.loads(source_summary.read_text())
    handoff, selected_metadata = build_handoff(registry, summary)

    out = args.out_dir.expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)
    csv_path = out / "v9_13_weakT_ceramic_paper_handoff.csv"
    json_path = out / "v9_13_weakT_ceramic_paper_handoff.json"
    handoff.to_csv(csv_path, index=False)

    records = handoff.to_dict(orient="records")
    fingerprint = candidate_parameter_fingerprint(records)
    manifest = {
        "schema": "v9.13_weakT_ceramic_paper_handoff_v1",
        "source_selection_schema": summary.get("schema"),
        "source_selection_dir": str(selection_dir),
        "source_registry": str(source_registry),
        "source_registry_sha256": sha256_path(source_registry),
        "source_summary": str(source_summary),
        "source_summary_sha256": sha256_path(source_summary),
        "handoff_csv": str(csv_path),
        "handoff_csv_sha256": sha256_path(csv_path),
        "active_parameter_fields": list(ACTIVE_CANDIDATE_PARAMETER_FIELDS),
        "active_parameter_fingerprint_sha256": fingerprint,
        "selected": selected_metadata,
        "fixed_closure": {
            "persistent_sites": True,
            "finite_source_inventory": False,
            "source_depletion_on_emission": False,
            "source_refresh_on_crack_advance": False,
            "explicit_recovery": False,
            "dynamic_tip_radius": True,
            "dynamic_front_width": True,
        },
        "transfer_policy": (
            "Exact active-row transfer only; no fitting, transformation, rounding, "
            "or substitution of inactive legacy source/recovery coordinates."
        ),
    }
    json_path.write_text(json.dumps(json_safe(manifest), indent=2, sort_keys=True) + "\n")

    print(
        "V913_WEAKT_CERAMIC_HANDOFF_WRITTEN "
        f"rows={len(handoff)} fingerprint={fingerprint} csv={csv_path} manifest={json_path}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
