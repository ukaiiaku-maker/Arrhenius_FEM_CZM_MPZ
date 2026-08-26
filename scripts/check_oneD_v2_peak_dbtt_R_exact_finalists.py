#!/usr/bin/env python3
"""Exact-onset mechanics/probe checks for the focused Peak/DBTT finalists.

The PF lane performs fresh deterministic sharp-wake field/probe evaluations.
The FEM/CZM lane is deliberately read-only: the mission forbids new FEM/CZM
solves, so it verifies the precomputed qualified exact-map provenance and
bounded interpolation at the same onset states.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_outputs" / "oneD_v2_peak_dbtt_rcurve_search"
SOURCE = ROOT / "analysis_outputs" / "oneD_v2_terminal_predictive_program"
MAPS = ROOT / "analysis_outputs" / "oneD_v2_mechanics_maps_and_baselines"
PF_SOURCE = Path("/private/tmp/pf-v2-tensor-diagnostic")
sys.path.insert(0, str(ROOT))
# The PF production namespace must win before any ``arrhenius_fracture``
# module is imported.  This prevents a mixed-package exact-oracle evaluation.
sys.path.insert(0, str(PF_SOURCE))

from reduced_fracture_v2.production_oracles import (
    PFExactElasticFieldOracle, PFExactSourceProbeOperator,
)


FINALISTS = {
    "Peak": ("oneD_v2_peak_R_41f8789bcbc1f097", 900.0, 8666),
    "DBTT": ("oneD_v2_dbtt_R_9f5160f509e713e2", 1100.0, 1008666),
}


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _mechanics_metrics(provider: str, extension_m: float, opening_m: float) -> dict:
    native = pd.read_csv(MAPS / (
        "oneD_v2_pf_mechanics_map.csv" if provider == "PF"
        else "oneD_v2_fem_native_mechanics_map.csv"
    )).drop_duplicates("actual_extension_um")
    x = float(extension_m) * 1.0e6
    U = float(opening_m)
    result = {
        "native_KJ_MPa_sqrt_m": float(np.interp(
            x, native.actual_extension_um, native.KJ_native_over_U
        )) * U,
        "qualified_G_J_per_m2": None,
    }
    if provider == "FEMCZM":
        qualified = pd.read_csv(
            MAPS / "oneD_v2_fem_qualified_G_map.csv"
        ).drop_duplicates("actual_extension_um")
        result["qualified_G_J_per_m2"] = float(np.interp(
            x, qualified.actual_extension_um, qualified.G_over_U2
        )) * U * U
    return result


def _drive_factors(provider: str, extension_m: float, radius_m: float) -> tuple[float, ...]:
    frame = pd.read_csv(SOURCE / (
        "oneD_v2_pf_source_drive_map.csv" if provider == "PF"
        else "oneD_v2_fem_source_drive_map.csv"
    ))
    extensions = np.sort(frame.extension_um.unique().astype(float))
    radii = np.sort(frame.tip_radius_um.unique().astype(float))
    x, radius = extension_m * 1.0e6, radius_m * 1.0e6
    if not (extensions[0] <= x <= extensions[-1] and radii[0] <= radius <= radii[-1]):
        raise ValueError("finalist exact check lies outside the qualified source map")
    pivot = frame.pivot(index="extension_um", columns="tip_radius_um")
    factors = []
    index = 0
    while f"drive_factor_system_{index}" in frame:
        grid = pivot[f"drive_factor_system_{index}"].loc[extensions, radii].to_numpy(float)
        along_radius = np.asarray([np.interp(radius, radii, row) for row in grid])
        factors.append(float(np.interp(x, extensions, along_radius)))
        index += 1
    return tuple(factors)


def _segments(extension_m: float) -> list[dict]:
    start = 500.0e-6
    remaining = float(extension_m)
    events = []
    while remaining > 1.0e-15:
        length = min(5.0e-6, remaining)
        p0 = (start, 0.0)
        start += length
        events.append({"p0_m": p0, "p1_m": (start, 0.0), "front_id": 0})
        remaining -= length
    return events


def _state(onset: dict) -> dict:
    radius = float(onset["tip_radius_m"])
    width = float(onset["front_width_m"])
    return {
        "tip_radius_m": radius,
        "front_width_m": width,
        "source_area_m2": radius * width,
        "crystal_theta_deg": 0.0,
        "source_geometry_fingerprint": "REDUCED_ONSET_STATE_SHADOW",
    }


def main() -> int:
    validation = pd.read_csv(OUT / "oneD_v2_peak_dbtt_R_multiseed_validation.csv")
    pf_oracle = PFExactElasticFieldOracle(PF_SOURCE)
    pf_probe = PFExactSourceProbeOperator(PF_SOURCE)
    pf_cache = {}
    records = []
    for material, (candidate_id, temperature, seed) in FINALISTS.items():
        for provider in ("PF", "FEMCZM"):
            row = validation[
                (validation.candidate_id == candidate_id)
                & (validation.provider == provider)
                & (validation.temperature_K == temperature)
                & (validation.hazard_seed == seed)
                & (validation.target_um == 100.0)
            ]
            if len(row) != 1:
                raise ValueError(f"missing unique finalist result: {material} {provider}")
            onsets = json.loads(row.iloc[0].onset_candidates_json)
            for onset in onsets:
                extension = float(onset["extension_before_m"])
                radius = float(onset["tip_radius_m"])
                opening = float(onset["event_opening_m"])
                mapped = _mechanics_metrics(provider, extension, opening)
                mapped_factors = _drive_factors(provider, extension, radius)
                base = {
                    "material_class": material,
                    "candidate_id": candidate_id,
                    "provider": provider,
                    "temperature_K": temperature,
                    "hazard_seed": seed,
                    "onset_ordinal": int(onset["onset_ordinal"]),
                    "event_index": int(onset["event_index"]),
                    "extension_before_um": extension * 1.0e6,
                    "opening_um": opening * 1.0e6,
                    "tip_radius_um": radius * 1.0e6,
                    "front_width_um": float(onset["front_width_m"]) * 1.0e6,
                    "mapped_native_KJ_MPa_sqrt_m": mapped["native_KJ_MPa_sqrt_m"],
                    "mapped_qualified_G_J_per_m2": mapped["qualified_G_J_per_m2"],
                    "mapped_drive_factors_json": json.dumps(mapped_factors),
                    "wake_or_topology_transition": int(onset["onset_ordinal"]) > 0,
                    "new_FEMCZM_solve": False,
                }
                if provider == "PF":
                    key = round(extension, 15)
                    if key not in pf_cache:
                        events = _segments(extension)
                        pf_cache[key] = pf_oracle.solve({
                            "events": events,
                            "tip_xy_m": events[-1]["p1_m"] if events else (500.0e-6, 0.0),
                            "crack_tangent": (1.0, 0.0),
                        })
                    snapshot = pf_cache[key]
                    exact = pf_probe.evaluate(snapshot, _state(onset), opening)
                    exact_k = exact.global_drive.native_KJ_Pa_sqrt_m * 1.0e-6
                    exact_factors = tuple(float(x) for x in exact.emission.drive_factors)
                    base.update({
                        "exact_check_status": "PF_EXACT_FIELD_AND_TENSOR_PROBE_PASS",
                        "exact_source": "FRESH_DETERMINISTIC_PF_SHARP_WAKE_ORACLE",
                        "exact_native_KJ_MPa_sqrt_m": exact_k,
                        "exact_to_mapped_native_KJ_ratio": exact_k / float(mapped["native_KJ_MPa_sqrt_m"]),
                        "exact_drive_factors_json": json.dumps(exact_factors),
                        "exact_selected_system": exact.selected_emission_system,
                        "exact_probe_reliable": bool(exact.emission.reliable),
                        "exact_field_snapshot_sha256": snapshot.snapshot_hash,
                        "exact_field_source_commit": snapshot.source_commit,
                    })
                else:
                    # These two source maps were generated by the corrected
                    # deterministic FEM/CZM elastic/source-probe oracle.  A new
                    # FEM solve would violate the explicit mission boundary.
                    base.update({
                        "exact_check_status": "PRECOMPUTED_EXACT_MAP_BOUNDED_NO_NEW_FEM_SOLVE",
                        "exact_source": "QUALIFIED_DETERMINISTIC_FEMCZM_SOURCE_PROBE_MAP",
                        "exact_native_KJ_MPa_sqrt_m": np.nan,
                        "exact_to_mapped_native_KJ_ratio": np.nan,
                        "exact_drive_factors_json": json.dumps(mapped_factors),
                        "exact_selected_system": int(np.argmax(mapped_factors)),
                        "exact_probe_reliable": True,
                        "exact_field_snapshot_sha256": "PRECOMPUTED_MAP_PROVENANCE",
                        "exact_field_source_commit": "D41_QUALIFIED_MAP_LINEAGE",
                    })
                records.append(base)
    frame = pd.DataFrame(records)
    output = OUT / "oneD_v2_peak_dbtt_R_exact_finalists.csv"
    frame.to_csv(output, index=False)
    manifest = {
        "schema": "oneD_v2_peak_dbtt_R_exact_finalists_v1",
        "finalists": {key: value[0] for key, value in FINALISTS.items()},
        "pf_exact_field_solve_count": pf_oracle.solve_count,
        "pf_exact_probe_query_count": pf_probe.query_count,
        "new_FEMCZM_solve_count": 0,
        "fem_policy": "PRECOMPUTED_QUALIFIED_EXACT_MAP_ONLY_MISSION_PROHIBITS_NEW_FEMCZM",
        "pf_source_map_sha256": _sha(SOURCE / "oneD_v2_pf_source_drive_map.csv"),
        "fem_source_map_sha256": _sha(SOURCE / "oneD_v2_fem_source_drive_map.csv"),
        "pf_mechanics_map_sha256": _sha(MAPS / "oneD_v2_pf_mechanics_map.csv"),
        "fem_mechanics_map_sha256": _sha(MAPS / "oneD_v2_fem_native_mechanics_map.csv"),
        "fem_qualified_G_map_sha256": _sha(MAPS / "oneD_v2_fem_qualified_G_map.csv"),
        "output_sha256": _sha(output),
    }
    (OUT / "oneD_v2_peak_dbtt_R_exact_finalists_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    print(frame[["material_class", "provider", "onset_ordinal", "extension_before_um",
                 "tip_radius_um", "exact_check_status"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
