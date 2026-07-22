#!/usr/bin/env python3
"""Aggregate reduced work/energy diagnostics from a v9.12 1-D campaign."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


ENERGY_FIELDS = (
    "external_plastic_work_J_per_m",
    "nonlocal_shielding_work_J_per_m",
    "internal_stress_work_J_per_m",
    "effective_plastic_work_J_per_m",
    "effective_plastic_dissipation_J_per_m",
    "external_plastic_work_per_crack_area_J_m2",
    "effective_plastic_dissipation_per_crack_area_J_m2",
    "mobile_line_energy_J_per_m",
    "retained_line_energy_J_per_m",
    "total_line_energy_J_per_m",
)
EXPECTED_BOOKKEEPING = "reduced_1d_orowan_power_and_log_line_energy_v1"
EXPECTED_INTEGRATOR = (
    "coupled_mobile_retained_backward_euler_v2_post_emit_refresh"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True)
    parser.add_argument("--window-um", nargs=2, type=float, default=(10.0, 30.0))
    parser.add_argument(
        "--out-prefix",
        default=None,
        help="Default: <root>/energy_analysis",
    )
    return parser.parse_args()


def values_in_window(
    extension_um: np.ndarray,
    values: np.ndarray,
    window_um: tuple[float, float],
) -> np.ndarray:
    mask = (
        np.isfinite(extension_um)
        & np.isfinite(values)
        & (extension_um >= window_um[0])
        & (extension_um <= window_um[1])
    )
    if not np.any(mask):
        raise RuntimeError(f"no values inside developed window {window_um}")
    return values[mask]


def main() -> int:
    args = parse_args()
    root = Path(args.root)
    if not root.is_dir():
        raise RuntimeError(f"campaign root does not exist: {root}")
    window_um = tuple(float(value) for value in args.window_um)
    out_prefix = (
        Path(args.out_prefix)
        if args.out_prefix
        else root / "energy_analysis"
    )
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    temperature_rows: list[dict[str, object]] = []
    candidate_rows: list[dict[str, object]] = []
    candidate_dirs = sorted(
        path
        for path in root.iterdir()
        if path.is_dir() and not path.name.startswith("_")
    )
    for candidate_root in candidate_dirs:
        temperature_files = sorted(candidate_root.glob("T*K.json"))
        if not temperature_files:
            continue
        per_temperature: list[dict[str, object]] = []
        for path in temperature_files:
            payload = json.loads(path.read_text())
            metadata = payload.get("numerical_integration", {})
            if metadata.get("spatial_integrator") != EXPECTED_INTEGRATOR:
                raise RuntimeError(
                    f"wrong spatial integrator in {path}: "
                    f"{metadata.get('spatial_integrator')!r}"
                )
            if metadata.get("energy_bookkeeping") != EXPECTED_BOOKKEEPING:
                raise RuntimeError(
                    f"missing energy bookkeeping in {path}: "
                    f"{metadata.get('energy_bookkeeping')!r}"
                )
            if metadata.get("energy_bookkeeping_feedback_active") is not False:
                raise RuntimeError(f"energy diagnostics feed back in {path}")

            extension = np.asarray(payload["extensions_um"], dtype=float)
            delta_k = np.asarray(
                payload["delta_K_micro_MPa_sqrt_m"],
                dtype=float,
            )
            if extension.size == 0 or delta_k.shape != extension.shape:
                raise RuntimeError(f"invalid checkpoint arrays in {path}")
            row: dict[str, object] = {
                "candidate_id": str(payload["candidate_id"]),
                "temperature_K": float(payload["temperature_K"]),
                "developed_delta_K_micro_MPa_sqrt_m": float(
                    np.median(values_in_window(extension, delta_k, window_um))
                ),
                "final_extension_um": float(extension[-1]),
            }
            for field_name in ENERGY_FIELDS:
                if field_name not in payload:
                    raise RuntimeError(f"missing {field_name} in {path}")
                values = np.asarray(payload[field_name], dtype=float)
                if (
                    values.shape != extension.shape
                    or not np.all(np.isfinite(values))
                ):
                    raise RuntimeError(f"invalid {field_name} in {path}")
                developed_values = values_in_window(
                    extension,
                    values,
                    window_um,
                )
                row[f"developed_{field_name}"] = float(
                    np.median(developed_values)
                )
                row[f"final_{field_name}"] = float(values[-1])
            per_temperature.append(row)
            temperature_rows.append(row)

        table = pd.DataFrame(per_temperature).sort_values("temperature_K")
        temperatures = table["temperature_K"].to_numpy(dtype=float)
        low_mask = temperatures <= 700.0
        high_mask = temperatures >= 1000.0
        if not np.any(low_mask) or not np.any(high_mask):
            raise RuntimeError(
                f"candidate {candidate_root.name} lacks low/high temperature bands"
            )

        dissipation_field = (
            "developed_effective_plastic_dissipation_per_crack_area_J_m2"
        )
        external_field = "developed_external_plastic_work_per_crack_area_J_m2"
        retained_field = "developed_retained_line_energy_J_per_m"
        delta_field = "developed_delta_K_micro_MPa_sqrt_m"
        dissipation = table[dissipation_field].to_numpy(dtype=float)
        external_work = table[external_field].to_numpy(dtype=float)
        retained_energy = table[retained_field].to_numpy(dtype=float)
        delta_k_values = table[delta_field].to_numpy(dtype=float)

        candidate_rows.append(
            {
                "candidate_id": candidate_root.name,
                "temperature_count": int(len(table)),
                "peak_delta_K_micro_MPa_sqrt_m": float(np.max(delta_k_values)),
                "peak_delta_K_temperature_K": float(
                    temperatures[int(np.argmax(delta_k_values))]
                ),
                "maximum_developed_effective_dissipation_J_m2": float(
                    np.max(dissipation)
                ),
                "maximum_dissipation_temperature_K": float(
                    temperatures[int(np.argmax(dissipation))]
                ),
                "low_temperature_median_effective_dissipation_J_m2": float(
                    np.median(dissipation[low_mask])
                ),
                "high_temperature_median_effective_dissipation_J_m2": float(
                    np.median(dissipation[high_mask])
                ),
                "high_minus_low_effective_dissipation_J_m2": float(
                    np.median(dissipation[high_mask])
                    - np.median(dissipation[low_mask])
                ),
                "high_temperature_median_external_work_J_m2": float(
                    np.median(external_work[high_mask])
                ),
                "high_temperature_median_retained_line_energy_J_per_m": float(
                    np.median(retained_energy[high_mask])
                ),
                "dissipation_at_delta_K_peak_J_m2": float(
                    dissipation[int(np.argmax(delta_k_values))]
                ),
                "delta_K_at_dissipation_max_MPa_sqrt_m": float(
                    delta_k_values[int(np.argmax(dissipation))]
                ),
                "temperature_correlation_delta_K_vs_dissipation": float(
                    pd.Series(delta_k_values).corr(
                        pd.Series(dissipation),
                        method="spearman",
                    )
                ),
            }
        )

    if not temperature_rows:
        raise RuntimeError(f"no per-temperature records found under {root}")

    temperature_table = pd.DataFrame(temperature_rows).sort_values(
        ["candidate_id", "temperature_K"]
    )
    candidate_table = pd.DataFrame(candidate_rows).sort_values(
        [
            "high_temperature_median_effective_dissipation_J_m2",
            "peak_delta_K_micro_MPa_sqrt_m",
        ],
        ascending=False,
    )
    temperature_path = Path(f"{out_prefix}_temperature.csv")
    candidate_path = Path(f"{out_prefix}_candidate.csv")
    summary_path = Path(f"{out_prefix}_summary.json")
    temperature_table.to_csv(temperature_path, index=False)
    candidate_table.to_csv(candidate_path, index=False)

    summary = {
        "campaign_root": str(root.resolve()),
        "candidate_count": int(candidate_table["candidate_id"].nunique()),
        "temperature_record_count": int(len(temperature_table)),
        "developed_window_um": list(window_um),
        "spatial_integrator": EXPECTED_INTEGRATOR,
        "energy_bookkeeping": EXPECTED_BOOKKEEPING,
        "energy_bookkeeping_feedback_active": False,
        "interpretation": {
            "effective_plastic_dissipation": (
                "reduced 1-D nonnegative Orowan resolved-work proxy per crack area"
            ),
            "external_plastic_work": (
                "signed work from K_applied before shielding"
            ),
            "nonlocal_shielding_work": (
                "signed work associated with the K_shield projection"
            ),
            "internal_stress_work": (
                "signed work associated with the local tau_GND field"
            ),
            "line_energy": (
                "logarithmic dislocation line-energy proxy per unit front thickness"
            ),
            "not_equivalent_to": [
                "full FEM J-integral",
                "validated J-R curve",
                "Charpy impact energy",
            ],
        },
        "temperature_csv": str(temperature_path.resolve()),
        "candidate_csv": str(candidate_path.resolve()),
    }
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    print(
        "ENERGY_ANALYSIS_COMPLETE "
        f"candidates={summary['candidate_count']} "
        f"temperature_records={summary['temperature_record_count']} "
        f"candidate_csv={candidate_path}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
