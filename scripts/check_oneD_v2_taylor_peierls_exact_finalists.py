#!/usr/bin/env python3
"""Exact mechanics/probe audit for Taylor/Peierls finalists.

PF uses fresh deterministic sharp-wake elastic fields and tensor probes.
FEM/CZM remains read-only and uses only the qualified bounded exact maps.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_outputs" / "oneD_v2_taylor_peierls_rcurve_search"
SOURCE = ROOT / "analysis_outputs" / "oneD_v2_terminal_predictive_program"
MAPS = ROOT / "analysis_outputs" / "oneD_v2_mechanics_maps_and_baselines"
PF_SOURCE = Path("/private/tmp/pf-v2-tensor-diagnostic")
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(PF_SOURCE))

from reduced_fracture_v2.production_oracles import (
    PFExactElasticFieldOracle,
    PFExactSourceProbeOperator,
)


FINALISTS = {
    "Peak": {
        "temperature_K": 1000.0, "seed": 8666,
        "candidate_ids": (
            "v913_zeroD_sobol_0242980",
            "oneD_v2_peak_TP_6962e84b2eb78fbb",
        ),
    },
    "DBTT": {
        "temperature_K": 1100.0, "seed": 1008666,
        "candidate_ids": (
            "v913_zeroD_sobol_0202500",
            "oneD_v2_dbtt_TP_6ca03f05fbae34e9",
        ),
    },
}


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _mechanics_metrics(provider: str, extension_m: float, opening_m: float) -> dict:
    native = pd.read_csv(MAPS / (
        "oneD_v2_pf_mechanics_map.csv" if provider == "PF"
        else "oneD_v2_fem_native_mechanics_map.csv"
    )).drop_duplicates("actual_extension_um")
    extension_um = float(extension_m) * 1e6
    result = {
        "native_KJ_MPa_sqrt_m": float(np.interp(
            extension_um, native.actual_extension_um, native.KJ_native_over_U
        )) * float(opening_m),
        "qualified_G_J_per_m2": None,
    }
    if provider == "FEMCZM":
        qualified = pd.read_csv(MAPS / "oneD_v2_fem_qualified_G_map.csv").drop_duplicates(
            "actual_extension_um"
        )
        result["qualified_G_J_per_m2"] = float(np.interp(
            extension_um, qualified.actual_extension_um, qualified.G_over_U2
        )) * float(opening_m) ** 2
    return result


def _drive_factors(provider: str, extension_m: float, radius_m: float) -> tuple[float, ...]:
    frame = pd.read_csv(SOURCE / (
        "oneD_v2_pf_source_drive_map.csv" if provider == "PF"
        else "oneD_v2_fem_source_drive_map.csv"
    ))
    extensions = np.sort(frame.extension_um.unique().astype(float))
    radii = np.sort(frame.tip_radius_um.unique().astype(float))
    x, radius = extension_m * 1e6, radius_m * 1e6
    if not (extensions[0] <= x <= extensions[-1] and radii[0] <= radius <= radii[-1]):
        raise ValueError("finalist state lies outside the qualified source-drive map")
    pivot = frame.pivot(index="extension_um", columns="tip_radius_um")
    factors = []
    index = 0
    while f"drive_factor_system_{index}" in frame:
        grid = pivot[f"drive_factor_system_{index}"].loc[extensions, radii].to_numpy(float)
        at_radius = np.asarray([np.interp(radius, radii, row) for row in grid])
        factors.append(float(np.interp(x, extensions, at_radius)))
        index += 1
    return tuple(factors)


def _segments(extension_m: float) -> list[dict]:
    start, remaining, events = 500e-6, float(extension_m), []
    while remaining > 1e-15:
        length = min(5e-6, remaining)
        p0 = (start, 0.0)
        start += length
        events.append({"p0_m": p0, "p1_m": (start, 0.0), "front_id": 0})
        remaining -= length
    return events


def _state(onset: dict) -> dict:
    radius, width = float(onset["tip_radius_m"]), float(onset["front_width_m"])
    return {
        "tip_radius_m": radius, "front_width_m": width,
        "source_area_m2": radius * width, "crystal_theta_deg": 0.0,
        "source_geometry_fingerprint": "REDUCED_ONSET_STATE_SHADOW",
    }


def main() -> int:
    validation_path = OUT / "oneD_v2_taylor_peierls_multiseed_results.parquet"
    validation = pd.read_parquet(validation_path)
    pf_oracle = PFExactElasticFieldOracle(PF_SOURCE)
    pf_probe = PFExactSourceProbeOperator(PF_SOURCE)
    pf_cache: dict[float, object] = {}
    records = []
    for material, specification in FINALISTS.items():
        for candidate_id in specification["candidate_ids"]:
            for provider in ("PF", "FEMCZM"):
                row = validation[
                    validation.candidate_id.eq(candidate_id)
                    & validation.provider.eq(provider)
                    & validation.temperature_K.eq(specification["temperature_K"])
                    & validation.hazard_seed.eq(specification["seed"])
                    & validation.target_um.eq(300.0)
                ]
                if len(row) != 1:
                    raise RuntimeError(f"missing unique exact-finalist row {candidate_id} {provider}")
                for onset in json.loads(row.iloc[0].onset_candidates_json):
                    extension = float(onset["extension_before_m"])
                    opening = float(onset["event_opening_m"])
                    radius = float(onset["tip_radius_m"])
                    mapped = _mechanics_metrics(provider, extension, opening)
                    mapped_factors = _drive_factors(provider, extension, radius)
                    record = {
                        "material_class": material, "candidate_id": candidate_id,
                        "provider": provider, "temperature_K": specification["temperature_K"],
                        "hazard_seed": specification["seed"], "onset_ordinal": onset["onset_ordinal"],
                        "extension_before_um": extension * 1e6, "opening_um": opening * 1e6,
                        "tip_radius_um": radius * 1e6,
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
                                "tip_xy_m": events[-1]["p1_m"] if events else (500e-6, 0.0),
                                "crack_tangent": (1.0, 0.0),
                            })
                        snapshot = pf_cache[key]
                        exact = pf_probe.evaluate(snapshot, _state(onset), opening)
                        exact_k = exact.global_drive.native_KJ_Pa_sqrt_m * 1e-6
                        record.update({
                            "exact_check_status": "PF_EXACT_FIELD_AND_TENSOR_PROBE_PASS",
                            "exact_source": "FRESH_DETERMINISTIC_PF_SHARP_WAKE_ORACLE",
                            "exact_native_KJ_MPa_sqrt_m": exact_k,
                            "exact_to_mapped_native_KJ_ratio": exact_k / mapped["native_KJ_MPa_sqrt_m"],
                            "exact_drive_factors_json": json.dumps(tuple(float(x) for x in exact.emission.drive_factors)),
                            "exact_selected_system": exact.selected_emission_system,
                            "exact_probe_reliable": bool(exact.emission.reliable),
                            "exact_field_snapshot_sha256": snapshot.snapshot_hash,
                            "exact_field_source_commit": snapshot.source_commit,
                        })
                    else:
                        record.update({
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
                    records.append(record)
    frame = pd.DataFrame(records)
    output = OUT / "oneD_v2_taylor_peierls_exact_finalists.csv"
    frame.to_csv(output, index=False)
    manifest = {
        "schema": "oneD_v2_taylor_peierls_exact_finalists_v1",
        "finalists": FINALISTS, "pf_exact_field_solve_count": pf_oracle.solve_count,
        "pf_exact_probe_query_count": pf_probe.query_count, "new_FEMCZM_solve_count": 0,
        "validation_sha256": _sha(validation_path), "output_sha256": _sha(output),
        "pf_source_map_sha256": _sha(SOURCE / "oneD_v2_pf_source_drive_map.csv"),
        "fem_source_map_sha256": _sha(SOURCE / "oneD_v2_fem_source_drive_map.csv"),
    }
    (OUT / "oneD_v2_taylor_peierls_exact_finalists_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    print(frame[["material_class", "candidate_id", "provider", "onset_ordinal",
                 "extension_before_um", "exact_check_status"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
