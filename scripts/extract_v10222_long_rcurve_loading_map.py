#!/usr/bin/env python3
"""Extract an autonomous v9.13 loading map from one completed v10.2.22 case.

The extractor fails closed unless the source case uses the calibrated stochastic
contract (exponential cleavage thresholds and threshold-scaled event lengths).
It can also require the new long map to reproduce an existing calibrated map as
an exact prefix before accepting additional crack-growth events.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any, Mapping, Sequence

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import pandas as pd

from arrhenius_fracture.dbtt_long_alignment_v913 import loading_map_coverage_um
from arrhenius_fracture.emergent_gnd_rcurve_v913 import RCurveLoadingMap


_MAP_ARRAY_KEYS = (
    "K_per_U_MPa_sqrt_m_per_m",
    "threshold_actions",
    "path_advances_m",
    "projected_advances_m",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case-dir", type=Path, required=True)
    parser.add_argument("--steps-csv", type=Path)
    parser.add_argument("--events-json", type=Path)
    parser.add_argument("--run-args-json", type=Path)
    parser.add_argument("--stack-json", type=Path)
    parser.add_argument("--expected-prefix-loading-map", type=Path)
    parser.add_argument("--reference-candidate-id", required=True)
    parser.add_argument("--reference-temperature-K", type=float)
    parser.add_argument("--minimum-coverage-um", type=float, default=100.0)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--audit-csv", type=Path)
    return parser.parse_args()


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def unique_match(case_dir: Path, pattern: str, label: str) -> Path:
    matches = sorted(case_dir.glob(pattern))
    if len(matches) != 1:
        raise RuntimeError(
            f"expected one {label} matching {pattern!r} under {case_dir}, "
            f"found {len(matches)}"
        )
    return matches[0]


def _require_equal(actual: Any, expected: Any, label: str) -> None:
    if actual != expected:
        raise RuntimeError(
            f"uncalibrated stochastic reference: {label}={actual!r}, "
            f"expected {expected!r}"
        )


def validate_calibrated_stochastic_contract(stack: Mapping[str, Any]) -> dict[str, Any]:
    """Require the exact stochastic controls used by the calibrated 52 um map."""
    hazard = stack.get("stochastic_hazard")
    avalanche = stack.get("stochastic_avalanche")
    if not isinstance(hazard, Mapping) or not isinstance(avalanche, Mapping):
        raise RuntimeError(
            "v10_2_17_final_signed_stochastic_stack.json lacks stochastic contract"
        )

    _require_equal(hazard.get("mode"), "exponential", "stochastic_hazard.mode")
    _require_equal(
        hazard.get("distribution"),
        "exponential_unit_mean",
        "stochastic_hazard.distribution",
    )
    _require_equal(
        avalanche.get("mode"),
        "threshold_scaled",
        "stochastic_avalanche.mode",
    )
    _require_equal(
        bool(stack.get("event_length_uses_same_integrated_hazard_threshold")),
        True,
        "event_length_uses_same_integrated_hazard_threshold",
    )

    numeric_contract = {
        "stochastic_avalanche.minimum_factor": (avalanche.get("minimum_factor"), 0.5),
        "stochastic_avalanche.maximum_factor": (avalanche.get("maximum_factor"), 4.0),
        "stochastic_avalanche.geometry_subsegment_fraction": (
            avalanche.get("geometry_subsegment_fraction"),
            0.1,
        ),
    }
    for label, (actual, expected) in numeric_contract.items():
        try:
            value = float(actual)
        except (TypeError, ValueError):
            raise RuntimeError(
                f"uncalibrated stochastic reference: {label}={actual!r}"
            ) from None
        if not math.isclose(value, expected, rel_tol=0.0, abs_tol=1.0e-14):
            raise RuntimeError(
                f"uncalibrated stochastic reference: {label}={value:.17g}, "
                f"expected {expected:.17g}"
            )

    return {
        "hazard_mode": str(hazard["mode"]),
        "hazard_distribution": str(hazard["distribution"]),
        "event_length_mode": str(avalanche["mode"]),
        "event_length_minimum_factor": float(avalanche["minimum_factor"]),
        "event_length_maximum_factor": float(avalanche["maximum_factor"]),
        "event_length_subsegment_fraction": float(
            avalanche["geometry_subsegment_fraction"]
        ),
    }


def validate_expected_prefix(
    current: Mapping[str, Sequence[float]],
    expected: Mapping[str, Any],
) -> dict[str, Any]:
    """Require every calibrated map array to be an exact prefix of the long map."""
    lengths: set[int] = set()
    maximum_errors: dict[str, float] = {}
    for key in _MAP_ARRAY_KEYS:
        if key not in expected:
            raise KeyError(f"expected prefix loading map is missing {key}")
        reference = [float(value) for value in expected[key]]
        observed = [float(value) for value in current[key]]
        lengths.add(len(reference))
        if len(reference) > len(observed):
            raise RuntimeError(
                f"expected prefix {key} has {len(reference)} entries but long map "
                f"has only {len(observed)}"
            )
        errors = [abs(a - b) for a, b in zip(observed, reference, strict=False)]
        maximum = max(errors, default=0.0)
        maximum_errors[key] = maximum
        for index, (actual, target) in enumerate(
            zip(observed, reference, strict=False)
        ):
            absolute_tolerance = 1.0e-9 if key.startswith("K_per_U") else 1.0e-14
            if not math.isclose(
                actual,
                target,
                rel_tol=1.0e-12,
                abs_tol=absolute_tolerance,
            ):
                raise RuntimeError(
                    "long loading map does not reproduce calibrated prefix: "
                    f"{key}[{index}]={actual:.17g}, expected {target:.17g}"
                )
    if len(lengths) != 1:
        raise RuntimeError(
            f"expected prefix loading-map arrays have inconsistent lengths: {lengths}"
        )

    prefix_events = next(iter(lengths))
    for key in ("seed", "nominal_dU_m", "nominal_dt_s"):
        if key not in expected:
            raise KeyError(f"expected prefix loading map is missing {key}")
        actual = current[key]
        target = expected[key]
        if key == "seed":
            if int(actual) != int(target):
                raise RuntimeError(
                    f"long loading-map seed {actual} does not match prefix seed {target}"
                )
        elif not math.isclose(
            float(actual), float(target), rel_tol=0.0, abs_tol=1.0e-15
        ):
            raise RuntimeError(
                f"long loading-map {key}={actual} does not match prefix {target}"
            )

    return {
        "prefix_events": prefix_events,
        "prefix_coverage_um": (
            sum(float(value) for value in expected["projected_advances_m"]) * 1.0e6
        ),
        "maximum_absolute_errors": maximum_errors,
    }


def main() -> int:
    args = parse_args()
    case_dir = args.case_dir.resolve()
    steps_path = args.steps_csv or unique_match(case_dir, "steps_*K.csv", "steps CSV")
    events_path = args.events_json or case_dir / "stochastic_avalanche_geometry_events.json"
    run_args_path = args.run_args_json or case_dir / "run_args.json"
    stack_path = args.stack_json or case_dir / "v10_2_17_final_signed_stochastic_stack.json"
    for path in (steps_path, events_path, run_args_path, stack_path):
        if not path.is_file():
            raise FileNotFoundError(path)
    if args.minimum_coverage_um <= 0.0:
        raise ValueError("minimum coverage must be positive")

    stack = json.loads(stack_path.read_text())
    stochastic_contract = validate_calibrated_stochastic_contract(stack)

    steps = pd.read_csv(steps_path)
    required_columns = {
        "step",
        "Uapp_m",
        "KJ_Pa_sqrtm",
        "da_block_m",
        "crack_extension_m",
    }
    missing = sorted(required_columns - set(steps.columns))
    if missing:
        raise KeyError(f"steps CSV is missing columns: {missing}")
    accepted = steps[steps["da_block_m"].astype(float) > 0.0].copy()
    accepted = accepted.sort_values("step", kind="stable").reset_index(drop=True)

    events_payload = json.loads(events_path.read_text())
    if not isinstance(events_payload, list) or not events_payload:
        raise RuntimeError("geometry event JSON must contain a nonempty list")
    events = sorted(events_payload, key=lambda row: int(row["event_index"]))
    expected_indices = list(range(len(events)))
    actual_indices = [int(row["event_index"]) for row in events]
    if actual_indices != expected_indices:
        raise RuntimeError(
            f"geometry event indices are not contiguous: {actual_indices[:10]}"
        )
    if len(accepted) != len(events):
        raise RuntimeError(
            f"accepted step count {len(accepted)} does not match geometry event "
            f"count {len(events)}"
        )

    run_args = json.loads(run_args_path.read_text())
    nominal_dU = float(run_args["dU"])
    nominal_dt = float(run_args["dt"])
    seed_values = {int(row["hazard_seed"]) for row in events}
    if len(seed_values) != 1:
        raise RuntimeError(f"geometry events contain multiple hazard seeds: {seed_values}")
    seed = next(iter(seed_values))
    if args.reference_temperature_K is not None:
        reference_temperature = float(args.reference_temperature_K)
    else:
        temperatures = run_args.get("temperatures", ())
        if not isinstance(temperatures, list) or len(temperatures) != 1:
            raise RuntimeError(
                "reference temperature is not explicit and run_args temperatures "
                "does not contain exactly one value"
            )
        reference_temperature = float(temperatures[0])

    geometry_factors: list[float] = []
    threshold_actions: list[float] = []
    path_advances: list[float] = []
    projected_advances: list[float] = []
    reference_K: list[float] = []
    reference_U: list[float] = []
    audit_rows: list[dict[str, Any]] = []

    for index, (event, (_, step)) in enumerate(
        zip(events, accepted.iterrows(), strict=True)
    ):
        U_m = float(step["Uapp_m"])
        K_MPa = float(step["KJ_Pa_sqrtm"]) * 1.0e-6
        if not math.isfinite(U_m) or U_m <= 0.0:
            raise ValueError(f"event {index} has invalid applied displacement {U_m}")
        if not math.isfinite(K_MPa) or K_MPa <= 0.0:
            raise ValueError(f"event {index} has invalid KJ {K_MPa}")
        factor = K_MPa / U_m
        threshold = float(event["threshold_action"])
        path_advance = float(event["event_advance_m"])
        projected_advance = float(event["x1"]) - float(event["x0"])
        if projected_advance <= 0.0:
            raise ValueError(
                f"event {index} has nonpositive projected x advance "
                f"{projected_advance}"
            )
        reported_projected = float(step["da_block_m"])
        projected_error = projected_advance - reported_projected
        if abs(projected_error) > max(1.0e-12, 1.0e-6 * projected_advance):
            raise RuntimeError(
                f"event {index} geometry/steps projected advance mismatch: "
                f"{projected_advance} versus {reported_projected}"
            )
        geometry_factors.append(factor)
        threshold_actions.append(threshold)
        path_advances.append(path_advance)
        projected_advances.append(projected_advance)
        reference_K.append(K_MPa)
        reference_U.append(U_m)
        audit_rows.append(
            {
                "event_index": index,
                "step": int(round(float(step["step"]))),
                "threshold_action": threshold,
                "U_m": U_m,
                "K_MPa_sqrt_m": K_MPa,
                "K_per_U_MPa_sqrt_m_per_m": factor,
                "path_advance_m": path_advance,
                "steps_da_block_m": reported_projected,
                "projected_advance_difference_m": projected_error,
                "projected_advance_m": projected_advance,
                "cumulative_projected_extension_m": sum(projected_advances),
            }
        )

    raw_map: dict[str, Any] = {
        "K_per_U_MPa_sqrt_m_per_m": geometry_factors,
        "threshold_actions": threshold_actions,
        "path_advances_m": path_advances,
        "projected_advances_m": projected_advances,
        "nominal_dU_m": nominal_dU,
        "nominal_dt_s": nominal_dt,
        "seed": seed,
    }
    prefix_audit: dict[str, Any] | None = None
    if args.expected_prefix_loading_map is not None:
        prefix_path = args.expected_prefix_loading_map.resolve()
        if not prefix_path.is_file():
            raise FileNotFoundError(prefix_path)
        expected_prefix = json.loads(prefix_path.read_text())
        prefix_audit = validate_expected_prefix(raw_map, expected_prefix)
    else:
        prefix_path = None

    loading_map = RCurveLoadingMap(
        K_per_U_MPa_sqrt_m_per_m=tuple(geometry_factors),
        threshold_actions=tuple(threshold_actions),
        path_advances_m=tuple(path_advances),
        projected_advances_m=tuple(projected_advances),
        nominal_dU_m=nominal_dU,
        nominal_dt_s=nominal_dt,
        seed=seed,
        reference_candidate_id=str(args.reference_candidate_id),
        reference_temperature_K=reference_temperature,
        reference_event_K_MPa_sqrt_m=tuple(reference_K),
        reference_event_U_m=tuple(reference_U),
        provenance={
            "schema": "v9.13_loading_map_from_v10.2.22_case_v2",
            "case_dir": str(case_dir),
            "steps_csv": str(steps_path.resolve()),
            "steps_csv_sha256": sha256_path(steps_path),
            "geometry_events_json": str(events_path.resolve()),
            "geometry_events_json_sha256": sha256_path(events_path),
            "run_args_json": str(run_args_path.resolve()),
            "run_args_json_sha256": sha256_path(run_args_path),
            "stochastic_stack_json": str(stack_path.resolve()),
            "stochastic_stack_json_sha256": sha256_path(stack_path),
            "calibrated_stochastic_contract": stochastic_contract,
            "expected_prefix_loading_map": (
                str(prefix_path) if prefix_path is not None else None
            ),
            "expected_prefix_loading_map_sha256": (
                sha256_path(prefix_path) if prefix_path is not None else None
            ),
            "expected_prefix_audit": prefix_audit,
            "geometry_map": "accepted_event_KJ_over_applied_displacement",
            "state_translation": "event_advance_m",
            "R_curve_abscissa": "projected_x_advance",
        },
    )
    loading_map.validate()
    payload = loading_map.as_dict()
    coverage_um = loading_map_coverage_um(payload)
    if coverage_um + 1.0e-9 < float(args.minimum_coverage_um):
        raise RuntimeError(
            f"extracted map covers {coverage_um:.9g} um, below required "
            f"{args.minimum_coverage_um:.9g} um"
        )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")
    audit_path = args.audit_csv or args.out.with_suffix(".audit.csv")
    audit_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(audit_rows).to_csv(audit_path, index=False)
    prefix_text = (
        f" prefix_events={prefix_audit['prefix_events']}"
        if prefix_audit is not None
        else ""
    )
    print(
        "V913_LONG_LOADING_MAP_EXTRACTED "
        f"events={len(events)} coverage_um={coverage_um:.9g} seed={seed}"
        f"{prefix_text} out={args.out}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
