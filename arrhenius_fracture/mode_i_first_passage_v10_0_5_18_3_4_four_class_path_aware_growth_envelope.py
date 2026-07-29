"""v10.0.5.18.3.4 production entry with path-aware growth support."""
from __future__ import annotations

import json
import math
import os
from pathlib import Path
import sys
from typing import Any

import numpy as np
from scipy.spatial import ConvexHull

from . import mode_i_first_passage_v9_18_5_3 as _v91853
from . import mode_i_first_passage_v10_0_5_13_2_barrier_only as _v1005132
from . import mode_i_first_passage_v10_0_5_18_3_2_four_class_joint_K_single_trial_stochastic_emission as _base
from . import path_aware_growth_envelope_v10051834 as _envelope
from . import physical_refinement_mesh_v100510 as _physical
from .path_aware_growth_envelope_v10051834 import (
    CORRIDOR_SCHEMA,
    MODEL_ID as ENVELOPE_MODEL,
    PATH_AWARE_PHYSICAL_POLICY,
    path_aware_long_growth_envelope_mesh,
)

POINT_RELEASE = "10.0.5.18.3.4"
MODEL_ID = (
    "FEM_CZM_four_class_joint_K_single_trial_path_aware_"
    "growth_envelope_v10_0_5_18_3_4"
)
PRODUCTION_MANIFEST = "persistent_site_production_manifest_v10_0_5_18_3_4.json"
SELECTION_MANIFEST = "persistent_site_parameter_selection_v10_0_5_18_3_4.json"
TRANSFER_MANIFEST = "four_class_parameter_transfer_v10_0_5_18_3_4.json"
SEED_MANIFEST = "stochastic_seed_manifest_v10_0_5_18_3_4.json"

_OLD_FILES = {
    "persistent_site_production_manifest_v10_0_5_18_3_2.json": PRODUCTION_MANIFEST,
    "persistent_site_parameter_selection_v10_0_5_18_3_2.json": SELECTION_MANIFEST,
    "four_class_parameter_transfer_v10_0_5_18_3_2.json": TRANSFER_MANIFEST,
    "stochastic_seed_manifest_v10_0_5_18_3_2.json": SEED_MANIFEST,
}


def _option_value(argv: list[str], name: str) -> str | None:
    for index, token in enumerate(argv):
        if token == name and index + 1 < len(argv):
            return argv[index + 1]
        if token.startswith(name + "="):
            return token.split("=", 1)[1]
    return None


def _required_positive_option(argv: list[str], name: str) -> float:
    raw = _option_value(argv, name)
    if raw is None:
        raise SystemExit(f"v10.0.5.18.3.4 requires {name}")
    value = float(raw)
    if value <= 0.0:
        raise SystemExit(f"v10.0.5.18.3.4 requires positive {name}")
    return value


def _out_path(argv: list[str]) -> Path | None:
    raw = _option_value(argv, "--out")
    return None if raw is None else Path(raw).expanduser().resolve()


def _active_physical_refinement_provider(_function) -> bool:
    """Detect the validated provider through nested point-release wrappers."""
    return _physical._ACTIVE_SPEC is not None


def _event_endpoint(row: dict[str, Any]) -> np.ndarray | None:
    for name in (
        "actual_p1_m",
        "p1_m",
        "stochastic_p1_requested_m",
        "requested_p1_m",
    ):
        value = row.get(name)
        if value is None:
            continue
        point = np.asarray(value, dtype=float).reshape(-1)
        if point.size >= 2 and np.all(np.isfinite(point[:2])):
            return point[:2]
    if all(name in row for name in ("x1", "y1")):
        point = np.asarray([row["x1"], row["y1"]], dtype=float)
        if np.all(np.isfinite(point)):
            return point
    return None


def _event_startpoint(row: dict[str, Any]) -> np.ndarray | None:
    for name in ("p0_m", "actual_p0_m", "requested_p0_m"):
        value = row.get(name)
        if value is None:
            continue
        point = np.asarray(value, dtype=float).reshape(-1)
        if point.size >= 2 and np.all(np.isfinite(point[:2])):
            return point[:2]
    if all(name in row for name in ("x0", "y0")):
        point = np.asarray([row["x0"], row["y0"]], dtype=float)
        if np.all(np.isfinite(point)):
            return point
    return None


def _committed_endpoint_audit(
    out: Path,
    corridor: dict[str, Any],
) -> dict[str, Any]:
    events_path = out / "stochastic_geometry_events_v10_0_5_16.json"
    vertices_raw = corridor.get("reachable_envelope_vertices_m")
    if not events_path.is_file() or not vertices_raw:
        return {
            "committed_endpoint_audit_available": False,
            "committed_endpoint_count": 0,
            "all_committed_endpoints_inside_certified_envelope": None,
            "committed_endpoint_support_certified": None,
        }

    try:
        events = list(json.loads(events_path.read_text()) or [])
        vertices = np.asarray(vertices_raw, dtype=float)
        hull = ConvexHull(vertices)
    except Exception as exc:
        return {
            "committed_endpoint_audit_available": False,
            "committed_endpoint_count": 0,
            "committed_endpoint_audit_error": f"{type(exc).__name__}: {exc}",
            "all_committed_endpoints_inside_certified_envelope": None,
            "committed_endpoint_support_certified": None,
        }

    accepted = [
        dict(row)
        for row in events
        if bool(row.get("inserted", False))
        and float(row.get("moved_m", 0.0) or 0.0) > 0.0
    ]
    endpoint_rows: list[np.ndarray] = []
    segment_forward_cosines: list[float] = []
    for row in accepted:
        endpoint = _event_endpoint(row)
        if endpoint is not None:
            endpoint_rows.append(endpoint)
        start = _event_startpoint(row)
        if start is not None and endpoint is not None:
            segment = endpoint - start
            length = float(np.linalg.norm(segment))
            if length > 1.0e-30:
                segment_forward_cosines.append(float(segment[0] / length))

    if not endpoint_rows:
        return {
            "committed_endpoint_audit_available": True,
            "committed_endpoint_count": 0,
            "accepted_geometry_event_count": len(accepted),
            "all_committed_endpoints_inside_certified_envelope": True,
            "minimum_committed_segment_global_forward_cosine": None,
            "all_committed_segments_pass_global_forward_gate": True,
            "committed_endpoint_support_certified": True,
        }

    endpoints = np.asarray(endpoint_rows, dtype=float)
    normals = np.asarray(hull.equations[:, :2], dtype=float)
    offsets = np.asarray(hull.equations[:, 2], dtype=float)
    norms = np.linalg.norm(normals, axis=1)
    excess = np.max(
        (endpoints @ normals.T + offsets[None, :])
        / np.maximum(norms[None, :], 1.0e-300),
        axis=1,
    )
    tolerance_m = max(
        float(corridor.get("physical_da_um", 5.0) or 5.0) * 1.0e-7,
        1.0e-9,
    )
    all_inside = bool(np.all(excess <= tolerance_m))
    minimum_forward = (
        None
        if not segment_forward_cosines
        else float(min(segment_forward_cosines))
    )
    required_forward = float(
        corridor.get("minimum_global_forward_cosine", 0.05) or 0.0
    )
    forward_ok = bool(
        minimum_forward is None
        or required_forward <= -1.0
        or minimum_forward >= required_forward - 1.0e-9
    )
    startup_support_ok = bool(
        float(
            corridor.get(
                "maximum_sampled_hbar_tip_over_L_pz",
                math.inf,
            )
        )
        <= float(
            corridor.get(
                "maximum_hbar_tip_over_L_pz_required",
                0.25,
            )
        )
        + 1.0e-12
    )
    return {
        "committed_endpoint_audit_available": True,
        "committed_endpoint_count": int(len(endpoints)),
        "accepted_geometry_event_count": int(len(accepted)),
        "committed_endpoint_tolerance_m": float(tolerance_m),
        "maximum_committed_endpoint_envelope_excess_m": float(np.max(excess)),
        "all_committed_endpoints_inside_certified_envelope": all_inside,
        "minimum_committed_segment_global_forward_cosine": minimum_forward,
        "all_committed_segments_pass_global_forward_gate": forward_ok,
        "startup_reachable_envelope_support_passed": startup_support_ok,
        "committed_endpoint_support_certified": bool(
            all_inside and forward_ok and startup_support_ok
        ),
    }


def _rewrite_outputs(
    out: Path | None,
    target_um: float,
    da_um: float,
    lpz_um: float,
    theta_deg: float,
    min_global_forward: float,
) -> None:
    if out is None:
        return
    corridor_path = out / "compact_corridor_mesh_v91852.json"
    corridor: dict[str, Any] = {}
    if corridor_path.is_file():
        try:
            corridor = dict(json.loads(corridor_path.read_text()) or {})
        except Exception:
            corridor = {}

    endpoint_audit = _committed_endpoint_audit(out, corridor)
    corridor.update(endpoint_audit)
    if corridor:
        corridor_path.write_text(
            json.dumps(corridor, indent=2, sort_keys=True, default=str) + "\n"
        )

    fields = {
        "path_aware_long_growth_envelope": True,
        "path_aware_envelope_schema": CORRIDOR_SCHEMA,
        "path_aware_envelope_model": ENVELOPE_MODEL,
        "path_aware_physical_refinement_policy": PATH_AWARE_PHYSICAL_POLICY,
        "directional_support_basis": (
            "complete_global_forward_gate_conservative_over_"
            "continuous_crack_direction_competition"
        ),
        "committed_target_extension_um_propagated_before_mesh": target_um,
        "physical_event_length_um_propagated_before_mesh": da_um,
        "process_zone_length_um_propagated_before_mesh": lpz_um,
        "crystal_theta_deg_propagated_before_mesh": theta_deg,
        "minimum_global_forward_cosine_propagated_before_mesh": (
            min_global_forward
        ),
        "corridor_acceptance_metric": (
            "triangle_quality_and_full_reachable_envelope_"
            "local_mean_edge_over_L_pz"
        ),
        "tip_h_over_da_role": "audit_warning_only",
        "physical_provider_detection": "active_validated_refinement_spec",
        "child_area_ratio_floor_relaxed": False,
        "constitutive_physics_changed": False,
        **endpoint_audit,
    }

    for old_name, new_name in _OLD_FILES.items():
        old = out / old_name
        if not old.is_file():
            continue
        payload = json.loads(old.read_text())
        payload["point_release"] = POINT_RELEASE
        payload["model"] = MODEL_ID
        if old_name.startswith("persistent_site_production_manifest"):
            payload["schema"] = (
                "persistent_site_production_manifest_v10_0_5_18_3_4"
            )
            physics = dict(payload.get("physics_contract", {}) or {})
            physics.update(fields)
            payload["physics_contract"] = physics
            payload["long_growth_corridor"] = corridor
        elif old_name.startswith("persistent_site_parameter_selection"):
            payload["schema"] = MODEL_ID
            policy = dict(payload.get("policy", {}) or {})
            policy.update(fields)
            payload["policy"] = policy
        else:
            payload.update(fields)
        new = out / new_name
        new.write_text(
            json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n"
        )
        old.unlink()


def main(argv: list[str] | None = None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    target_um = _required_positive_option(
        user_args,
        "--target-crack-extension-um",
    )
    da_m = _required_positive_option(user_args, "--da-phys")
    da_um = da_m * 1.0e6
    lpz_raw = _option_value(user_args, "--mpz-length-um")
    lpz_um = 50.0 if lpz_raw is None else float(lpz_raw)
    if lpz_um <= 0.0:
        raise SystemExit(
            "v10.0.5.18.3.4 requires positive --mpz-length-um"
        )
    theta_raw = _option_value(user_args, "--crystal-theta-deg")
    theta_deg = 45.0 if theta_raw is None else float(theta_raw)
    forward_raw = _option_value(user_args, "--min-global-forward")
    min_global_forward = 0.05 if forward_raw is None else float(forward_raw)
    if not math.isfinite(min_global_forward) or min_global_forward > 1.0:
        raise SystemExit(
            "v10.0.5.18.3.4 requires finite --min-global-forward <= 1"
        )
    out = _out_path(user_args)

    requested_env = {
        "ARRHENIUS_COMMITTED_TARGET_EXTENSION_UM": str(target_um),
        "ARRHENIUS_PHYSICAL_DA_UM": str(da_um),
        "ARRHENIUS_CORRIDOR_PROCESS_ZONE_UM": str(lpz_um),
        "ARRHENIUS_CRYSTAL_THETA_DEG": str(theta_deg),
        "ARRHENIUS_MIN_GLOBAL_FORWARD": str(min_global_forward),
        "ARRHENIUS_MAX_CORRIDOR_H_OVER_LPZ": os.environ.get(
            "ARRHENIUS_MAX_CORRIDOR_H_OVER_LPZ",
            "0.25",
        ),
        "ARRHENIUS_ENVELOPE_SUPPORT_H_OVER_LPZ": os.environ.get(
            "ARRHENIUS_ENVELOPE_SUPPORT_H_OVER_LPZ",
            "0.245",
        ),
        "ARRHENIUS_ENVELOPE_BUFFER_LPZ": os.environ.get(
            "ARRHENIUS_ENVELOPE_BUFFER_LPZ",
            "1.0",
        ),
    }
    saved_env = {key: os.environ.get(key) for key in requested_env}
    saved_v91853_corridor = _v91853._quality_selected_corridor_mesh
    saved_v1005132_corridor = (
        _v1005132._quality_selected_corridor_mesh_v1005132
    )
    saved_provider_detector = (
        _envelope._is_production_physical_constructor
    )

    _v91853._quality_selected_corridor_mesh = (
        path_aware_long_growth_envelope_mesh
    )
    _v1005132._quality_selected_corridor_mesh_v1005132 = (
        path_aware_long_growth_envelope_mesh
    )
    _envelope._is_production_physical_constructor = (
        _active_physical_refinement_provider
    )
    for key, value in requested_env.items():
        os.environ[key] = value

    try:
        return _base.main(user_args)
    finally:
        _envelope._is_production_physical_constructor = (
            saved_provider_detector
        )
        _v1005132._quality_selected_corridor_mesh_v1005132 = (
            saved_v1005132_corridor
        )
        _v91853._quality_selected_corridor_mesh = saved_v91853_corridor
        for key, value in saved_env.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value
        _rewrite_outputs(
            out,
            target_um,
            da_um,
            lpz_um,
            theta_deg,
            min_global_forward,
        )


if __name__ == "__main__":
    main()


__all__ = [
    "MODEL_ID",
    "POINT_RELEASE",
    "PRODUCTION_MANIFEST",
    "SEED_MANIFEST",
    "SELECTION_MANIFEST",
    "TRANSFER_MANIFEST",
    "main",
]
