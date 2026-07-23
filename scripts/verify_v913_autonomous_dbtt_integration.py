#!/usr/bin/env python3
"""Fail-closed preflight for the integrated v9.13 DBTT search.

The checks tie the candidate pool, active parameter contract, common physics,
stochastic loading map, accepted 50-case calibration, and optional executable
sentinel cases together before a long acquisition wave is allowed to start.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from arrhenius_fracture.emergent_gnd_campaign_v913 import (
    candidate_from_registry_row,
)
from arrhenius_fracture.emergent_gnd_contract_v913 import (
    ACTIVE_CANDIDATE_PARAMETER_FIELDS,
    PERSISTENT_INACTIVE_REGISTRY_FIELDS,
    candidate_feature_record,
    candidate_parameter_fingerprint,
    effective_candidate_parameters,
)
from arrhenius_fracture.emergent_gnd_rcurve_v913 import (
    RCurveLoadingMap,
    run_autonomous_rcurve,
)
from scripts.run_mpz_v9_13_persistent_top5 import load_physics


EXPECTED_POOL_FILE_SHA256 = (
    "1633851df78f4848a897d87e5ee9679f8e708095c1895369609bac2ab7c78efe"
)
EXPECTED_POOL_PARAMETER_SHA256 = (
    "7903e0c65bb0da45ddccbb5a74814adc5b9713b1de4092b8c24e583fbe550083"
)
EXPECTED_TOP5_PARAMETER_SHA256 = (
    "a8befb167c06289fc0b19bc1452f6744501fa09647ac2c242efd9fbb846e137f"
)

SENTINEL_CASES = (
    ("v912_targeted_local_peak_013476_0368", 900.0),
    ("v912_targeted_local_peak_005518_0118", 1200.0),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--candidate-registry",
        type=Path,
        default=Path("candidates/v9_12_targeted_local_4096_registry.csv"),
    )
    parser.add_argument(
        "--top5-registry",
        type=Path,
        default=Path("candidates/v9_13_persistent_sites_top5_registry.csv"),
    )
    parser.add_argument(
        "--base-physics-json",
        type=Path,
        default=Path("mpz_v9_13_v10222_transfer_common_physics.json"),
    )
    parser.add_argument(
        "--loading-map",
        type=Path,
        default=Path(
            "runs/v9_13_v10222_rcurve_targets_v1/"
            "v10_2_22_rcurve_loading_map.json"
        ),
    )
    parser.add_argument(
        "--target-manifest",
        type=Path,
        default=Path(
            "runs/v9_13_v10222_rcurve_targets_v1/"
            "v9_13_v10_2_22_rcurve_target_manifest.json"
        ),
    )
    parser.add_argument(
        "--accepted-summary",
        type=Path,
        default=Path(
            "runs/v9_13_v10222_rcurve_alpha0p95_all50_v1/"
            "R_curve_calibration_summary.json"
        ),
    )
    parser.add_argument(
        "--sentinel-reference",
        type=Path,
        default=Path(
            "calibration_data/"
            "v9_13_v10222_rcurve_sentinel_reference.json"
        ),
    )
    parser.add_argument(
        "--run-sentinels",
        action=argparse.BooleanOptionalAction,
        default=False,
    )
    return parser.parse_args()


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise RuntimeError(f"empty CSV: {path}")
    return rows


def _sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _assert_close(
    actual: float,
    expected: float,
    *,
    name: str,
    atol: float = 1.0e-10,
    rtol: float = 1.0e-10,
) -> None:
    if not np.isclose(actual, expected, atol=atol, rtol=rtol):
        raise RuntimeError(f"{name}: actual={actual!r}, expected={expected!r}")


def _validate_pool(
    pool_rows: list[dict[str, str]],
    top5_rows: list[dict[str, str]],
    *,
    pool_path: Path,
) -> None:
    file_sha = _sha256_path(pool_path)
    if file_sha != EXPECTED_POOL_FILE_SHA256:
        raise RuntimeError(
            "4096-row registry file hash mismatch: "
            f"{file_sha} != {EXPECTED_POOL_FILE_SHA256}"
        )
    if len(pool_rows) != 4096:
        raise RuntimeError(f"candidate pool has {len(pool_rows)} rows, expected 4096")
    ids = [row["candidate_id"] for row in pool_rows]
    if len(set(ids)) != len(ids):
        raise RuntimeError("candidate pool contains duplicate candidate_id values")
    for row in pool_rows:
        effective_candidate_parameters(row)
    pool_fingerprint = candidate_parameter_fingerprint(pool_rows)
    if pool_fingerprint != EXPECTED_POOL_PARAMETER_SHA256:
        raise RuntimeError("candidate pool active-parameter fingerprint mismatch")
    top5_fingerprint = candidate_parameter_fingerprint(top5_rows)
    if top5_fingerprint != EXPECTED_TOP5_PARAMETER_SHA256:
        raise RuntimeError("top-five active-parameter fingerprint mismatch")

    by_id = {row["candidate_id"]: row for row in pool_rows}
    for top5 in top5_rows:
        candidate_id = top5["candidate_id"]
        if candidate_id not in by_id:
            raise RuntimeError(f"top-five candidate absent from pool: {candidate_id}")
        expected = effective_candidate_parameters(top5)
        actual = effective_candidate_parameters(by_id[candidate_id])
        for field in ACTIVE_CANDIDATE_PARAMETER_FIELDS:
            _assert_close(
                actual[field],
                expected[field],
                name=f"{candidate_id}.{field}",
                atol=0.0,
                rtol=5.0e-15,
            )

    feature_names = set(candidate_feature_record(pool_rows[0]))
    expected_features = {
        f"x_raw__{field}" for field in ACTIVE_CANDIDATE_PARAMETER_FIELDS
    }
    if feature_names != expected_features:
        raise RuntimeError("surrogate active-feature contract is incomplete")
    forbidden = {
        f"x_raw__{field}" for field in PERSISTENT_INACTIVE_REGISTRY_FIELDS
    }
    if feature_names & forbidden:
        raise RuntimeError("inactive legacy fields entered surrogate features")


def _validate_physics(
    physics_path: Path,
    loading_map_path: Path,
) -> tuple[Any, RCurveLoadingMap]:
    physics, metadata = load_physics(physics_path)
    common = metadata.get("common_physics", metadata)
    fixed = metadata.get("fixed_contract", {})
    required_false = (
        "finite_source_inventory",
        "source_depletion_on_emission",
        "source_refresh_on_crack_advance",
        "explicit_recovery",
    )
    for name in required_false:
        if fixed.get(name) is not False:
            raise RuntimeError(f"physics contract does not disable {name}")
    if fixed.get("persistent_multiplicity") is not True:
        raise RuntimeError("persistent multiplicity is not active")
    if fixed.get("front_width_independent_of_ahead_tip_mesh") is not True:
        raise RuntimeError("front width is not declared mesh-independent")
    if not bool(common["coupled_moving_tip_enabled"]):
        raise RuntimeError("accepted autonomous transfer requires coupled moving tip")
    _assert_close(
        float(common["minimum_front_width_m"]),
        float(common["b_m"]),
        name="physical minimum front width",
        atol=1.0e-30,
        rtol=1.0e-14,
    )
    dx = float(common["mpz_length_m"]) / int(common["n_bins"])
    if not float(common["minimum_front_width_m"]) < 1.0e-3 * dx:
        raise RuntimeError("front-width minimum remains coupled to MPZ spacing")
    _assert_close(
        float(common["mobile_transport_velocity_scale"]),
        0.0,
        name="validated-scalar transport scale",
        atol=0.0,
        rtol=0.0,
    )
    _assert_close(
        float(common["taylor_phi_max"]),
        20.0,
        name="Taylor amplification cap",
    )
    if any(float(value) != 0.0 for value in common["shielding_orientation_factors"]):
        raise RuntimeError("analytical 1-D signed shielding must remain disabled")

    loading_map = RCurveLoadingMap.from_dict(
        json.loads(loading_map_path.read_text())
    )
    if loading_map.seed != 3621 or loading_map.n_events != 16:
        raise RuntimeError("loading map is not the accepted CRN3621 16-event map")
    if loading_map.reference_candidate_id != (
        "v912_targeted_local_peak_013476_0368"
    ):
        raise RuntimeError("loading-map reference candidate changed")
    _assert_close(
        loading_map.reference_temperature_K,
        300.0,
        name="loading-map reference temperature",
        atol=0.0,
        rtol=0.0,
    )
    return physics, loading_map


def _validate_calibration(
    *,
    target_manifest_path: Path,
    accepted_summary_path: Path,
) -> None:
    target_manifest = json.loads(target_manifest_path.read_text())
    if target_manifest.get("candidate_registry_sha256") != (
        EXPECTED_TOP5_PARAMETER_SHA256
    ):
        raise RuntimeError("2-D target manifest has the wrong candidate fingerprint")
    if target_manifest.get("candidate_parameter_fields") != list(
        ACTIVE_CANDIDATE_PARAMETER_FIELDS
    ):
        raise RuntimeError("2-D target manifest omits active candidate fields")
    if (
        target_manifest.get("n_cases") != 50
        or target_manifest.get("n_events_per_case") != 16
    ):
        raise RuntimeError("2-D target manifest is not the accepted 50-case set")

    summary = json.loads(accepted_summary_path.read_text())
    if summary.get("candidate_registry_sha256") != (
        EXPECTED_TOP5_PARAMETER_SHA256
    ):
        raise RuntimeError("accepted calibration has the wrong candidate fingerprint")
    shared = summary["shared_calibration"]
    for name, value in shared.items():
        expected = 0.95 if name == "translation_action_exponent" else 1.0
        _assert_close(float(value), expected, name=f"shared calibration {name}")
    metrics = summary["metrics"]["all_K_checkpoints"]
    _assert_close(
        float(metrics["MAE"]),
        0.1330560783932988,
        name="50-case checkpoint MAE",
        atol=1.0e-12,
        rtol=1.0e-12,
    )
    _assert_close(
        float(metrics["RMSE"]),
        0.5986756471546335,
        name="50-case checkpoint RMSE",
        atol=1.0e-12,
        rtol=1.0e-12,
    )
    event_metrics = summary["event_metrics"]
    if int(event_metrics["n_events"]) != 800:
        raise RuntimeError("accepted calibration does not contain 800 events")
    _assert_close(
        float(event_metrics["MAE"]),
        0.13082554834835303,
        name="800-event MAE",
        atol=1.0e-12,
        rtol=1.0e-12,
    )


def _run_sentinels(
    *,
    top5_rows: list[dict[str, str]],
    physics: Any,
    loading_map: RCurveLoadingMap,
    reference: Mapping[str, Any],
) -> None:
    if reference.get("schema") != (
        "v9.13_autonomous_R_curve_sentinel_reference_v1"
    ):
        raise RuntimeError("sentinel reference has the wrong schema")
    _assert_close(
        float(reference["translation_action_exponent"]),
        0.95,
        name="sentinel translation exponent",
    )
    _assert_close(
        float(reference["max_hazard_increment"]),
        0.05,
        name="sentinel hazard increment",
    )
    if int(reference["seed"]) != loading_map.seed:
        raise RuntimeError("sentinel reference seed differs from loading map")
    reference_by_case = {
        (str(row["candidate_id"]), float(row["temperature_K"])): row
        for row in reference["cases"]
    }
    if set(reference_by_case) != set(SENTINEL_CASES):
        raise RuntimeError("sentinel reference does not define the exact cases")

    by_id = {row["candidate_id"]: row for row in top5_rows}
    for index, (candidate_id, temperature_K) in enumerate(
        SENTINEL_CASES,
        start=1,
    ):
        print(
            "V913_DBTT_SENTINEL_START "
            f"index={index}/{len(SENTINEL_CASES)} "
            f"candidate={candidate_id} T={temperature_K:g}",
            flush=True,
        )
        result = run_autonomous_rcurve(
            candidate_from_registry_row(by_id[candidate_id]),
            physics,
            loading_map,
            temperature_K,
            target_projected_extension_m=25.0e-6,
            max_hazard_increment=0.05,
            translation_mode="hazard_coupled",
            translation_action_exponent=0.95,
        )
        if result.status != "complete":
            raise RuntimeError(
                f"sentinel {candidate_id}:{temperature_K:g} did not complete"
            )
        reference_row = reference_by_case[(candidate_id, temperature_K)]
        for label, actual in (
            ("K_first", result.checkpoint_K(0.0)),
            ("K_10um", result.checkpoint_K(10.0e-6)),
            ("K_25um", result.checkpoint_K(25.0e-6)),
        ):
            _assert_close(
                actual,
                float(reference_row["predicted"][label]),
                name=f"sentinel {candidate_id}:{temperature_K:g} {label}",
                atol=2.0e-9,
                rtol=2.0e-11,
            )

        expected_prefix = reference_row["events_through_25um"]
        if len(expected_prefix) != len(result.events):
            raise RuntimeError("accepted event table is missing a sentinel prefix")
        for actual, expected in zip(result.events, expected_prefix):
            if actual.event_index != int(expected["event_index"]):
                raise RuntimeError("sentinel event index changed")
            _assert_close(
                actual.threshold_action,
                float(expected["threshold_action"]),
                name=(
                    f"sentinel {candidate_id}:{temperature_K:g} "
                    f"event {actual.event_index} threshold"
                ),
                atol=1.0e-14,
                rtol=1.0e-14,
            )
            _assert_close(
                actual.cumulative_projected_extension_m,
                float(expected["cumulative_projected_extension_m"]),
                name=(
                    f"sentinel {candidate_id}:{temperature_K:g} "
                    f"event {actual.event_index} extension"
                ),
                atol=1.0e-18,
                rtol=1.0e-12,
            )
            _assert_close(
                actual.K_MPa_sqrt_m,
                float(expected["K_MPa_sqrt_m"]),
                name=(
                    f"sentinel {candidate_id}:{temperature_K:g} "
                    f"event {actual.event_index} K"
                ),
                atol=2.0e-9,
                rtol=2.0e-11,
            )
        print(
            "V913_DBTT_SENTINEL_COMPLETE "
            f"index={index}/{len(SENTINEL_CASES)} "
            f"candidate={candidate_id} T={temperature_K:g} "
            f"events={len(result.events)} "
            f"K25={result.checkpoint_K(25.0e-6):.8g}",
            flush=True,
        )


def main() -> int:
    args = parse_args()
    print("V913_DBTT_PREFLIGHT_START", flush=True)
    pool_rows = _read_csv(args.candidate_registry)
    top5_rows = _read_csv(args.top5_registry)
    _validate_pool(
        pool_rows,
        top5_rows,
        pool_path=args.candidate_registry,
    )
    print(
        "V913_DBTT_PREFLIGHT_CHECK "
        f"name=candidate_contract status=ok pool_rows={len(pool_rows)}",
        flush=True,
    )
    physics, loading_map = _validate_physics(
        args.base_physics_json,
        args.loading_map,
    )
    print(
        "V913_DBTT_PREFLIGHT_CHECK name=physics_loading_map status=ok",
        flush=True,
    )
    _validate_calibration(
        target_manifest_path=args.target_manifest,
        accepted_summary_path=args.accepted_summary,
    )
    print(
        "V913_DBTT_PREFLIGHT_CHECK name=accepted_calibration status=ok",
        flush=True,
    )
    if args.run_sentinels:
        _run_sentinels(
            top5_rows=top5_rows,
            physics=physics,
            loading_map=loading_map,
            reference=json.loads(args.sentinel_reference.read_text()),
        )
    print(
        "V913_DBTT_INTEGRATION_OK "
        f"pool_rows={len(pool_rows)} top5_rows={len(top5_rows)} "
        f"active_features={len(ACTIVE_CANDIDATE_PARAMETER_FIELDS)} "
        "checkpoint_MAE=0.133056078 event_MAE=0.130825548 "
        f"sentinels={len(SENTINEL_CASES) if args.run_sentinels else 0}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
