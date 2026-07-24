"""v10.0.5.16: stochastic PF hazard and variable event-length parity.

This point release retains the exact v10.0.5.15 PF local state update and
continuous moving-tip coupling.  It restores the production PF renewal law:
one unit-mean exponential integrated-hazard threshold per cleavage event and a
bounded, mean-preserving crack-length reward derived from that same threshold.
"""
from __future__ import annotations

import json
import math
import os
from pathlib import Path
import statistics
import sys
from typing import Any

import numpy as np

from . import mode_i_first_passage_v10_0_5_14_1_persistent_site_family as _base
from .adaptive_czm_tip_support_v1005144 import (
    MODEL_ID as TIP_SUPPORT_MODEL,
    installed_tip_support_repair_v1005144,
    tip_support_audit_v1005144,
)
from .persistent_site_diagnostics_v1005144 import (
    installed_persistent_diagnostics_v1005144,
)
from .persistent_site_diagnostics_v100516 import (
    DIAGNOSTIC_MODEL,
    installed_persistent_diagnostics_v100516,
)
from .persistent_site_pf_update_v100515 import PF_REFERENCE_COMMIT, PF_UPDATE_MAP
from .persistent_site_stochastic_tip_v100516 import (
    EVENT_LENGTH_SCHEMA,
    HAZARD_SCHEMA,
    MODEL_ID as ENGINE_MODEL,
    PersistentSitePFStochasticMovingTipFrontEngineV100516,
)
from .refinement_audit_v100516 import (
    AUDIT_MODEL as REFINEMENT_AUDIT_MODEL,
    installed_refinement_audit_v100516,
    refinement_audit_payload_v100516,
)
from .stochastic_geometry_v100516 import (
    GEOMETRY_MODEL,
    installed_stochastic_geometry_v100516,
    realized_stochastic_geometry_events,
    write_stochastic_geometry_events,
)

POINT_RELEASE = "10.0.5.16"
MODEL_ID = (
    "FEM_CZM_full_2D_PF_v10_2_22_stochastic_hazard_"
    "threshold_correlated_event_length_v10_0_5_16"
)
PRODUCTION_MANIFEST = "persistent_site_production_manifest_v10_0_5_16.json"
SELECTION_MANIFEST = "persistent_site_parameter_selection_v10_0_5_16.json"


def _out_path(argv: list[str]) -> Path | None:
    if "--out" not in argv:
        return None
    index = argv.index("--out")
    if index + 1 >= len(argv):
        return None
    return Path(argv[index + 1]).expanduser().resolve()


def _require_stochastic_config() -> dict[str, Any]:
    hazard_mode = os.environ.get("CLEAVAGE_HAZARD_MODE", "exponential").strip().lower()
    length_mode = os.environ.get(
        "CLEAVAGE_EVENT_LENGTH_MODE", "threshold_scaled"
    ).strip().lower()
    if hazard_mode != "exponential":
        raise SystemExit(
            "v10.0.5.16 requires CLEAVAGE_HAZARD_MODE=exponential"
        )
    if length_mode != "threshold_scaled":
        raise SystemExit(
            "v10.0.5.16 requires CLEAVAGE_EVENT_LENGTH_MODE=threshold_scaled"
        )
    raw_seed = os.environ.get("CLEAVAGE_HAZARD_SEED", "").strip()
    if not raw_seed:
        raise SystemExit("v10.0.5.16 requires explicit CLEAVAGE_HAZARD_SEED")
    try:
        seed = int(raw_seed)
    except ValueError as exc:
        raise SystemExit("CLEAVAGE_HAZARD_SEED must be an integer") from exc
    if seed < 0:
        raise SystemExit("CLEAVAGE_HAZARD_SEED must be nonnegative")
    minimum_threshold = float(
        os.environ.get("CLEAVAGE_HAZARD_MIN_THRESHOLD", "1e-12")
    )
    minimum_factor = float(os.environ.get("CLEAVAGE_EVENT_MIN_FACTOR", "0.5"))
    maximum_factor = float(os.environ.get("CLEAVAGE_EVENT_MAX_FACTOR", "4.0"))
    if not all(
        math.isfinite(value) and value > 0.0
        for value in (minimum_threshold, minimum_factor, maximum_factor)
    ):
        raise SystemExit("stochastic threshold and event-length bounds must be positive")
    if maximum_factor < minimum_factor:
        raise SystemExit("CLEAVAGE_EVENT_MAX_FACTOR must be >= minimum factor")
    return {
        "hazard_mode": hazard_mode,
        "hazard_seed": seed,
        "hazard_minimum_threshold": minimum_threshold,
        "event_length_mode": length_mode,
        "event_minimum_factor": minimum_factor,
        "event_maximum_factor": maximum_factor,
    }


def _sum_nested(value: Any) -> float:
    try:
        return float(np.sum(np.asarray(value, dtype=float)))
    except Exception:
        return 0.0


def _final_line_content_accounting(out: Path) -> dict[str, Any] | None:
    records: list[dict[str, Any]] = []
    for path in sorted(out.glob("mpz_state_snapshots_*K.json")):
        try:
            payload = json.loads(path.read_text())
        except Exception:
            continue
        for front in payload.get("final_fronts", []):
            state = dict(front.get("state", {}) or {})
            active_mobile = _sum_nested(state.get("mobile_positive", [])) + _sum_nested(
                state.get("mobile_negative", [])
            )
            active_retained = _sum_nested(
                state.get("retained_positive", [])
            ) + _sum_nested(state.get("retained_negative", []))
            wake_mobile = _sum_nested(
                state.get("wake_mobile_positive", [])
            ) + _sum_nested(state.get("wake_mobile_negative", []))
            wake_retained = _sum_nested(
                state.get("wake_retained_positive", [])
            ) + _sum_nested(state.get("wake_retained_negative", []))
            emitted = float(state.get("emitted_total", 0.0) or 0.0)
            escaped = float(state.get("escaped_total", 0.0) or 0.0)
            recovered = float(state.get("recovered_total", 0.0) or 0.0)
            legacy_discarded = float(state.get("wake_discarded_total", 0.0) or 0.0)
            present_or_escaped = (
                active_mobile
                + active_retained
                + wake_mobile
                + wake_retained
                + escaped
                + recovered
            )
            discarded_line = max(emitted - present_or_escaped, 0.0)
            discarded_slip = max(legacy_discarded - discarded_line, 0.0)
            accounted = present_or_escaped + discarded_line
            balance = emitted - accounted
            scale = max(abs(emitted), abs(accounted), 1.0e-300)
            records.append(
                {
                    "temperature_K": payload.get("temperature_K"),
                    "front_id": front.get("front_id"),
                    "active_mobile_total": active_mobile,
                    "active_retained_total": active_retained,
                    "wake_mobile_total": wake_mobile,
                    "wake_retained_total": wake_retained,
                    "emitted_total": emitted,
                    "escaped_total": escaped,
                    "recovered_total": recovered,
                    "legacy_wake_discarded_total": legacy_discarded,
                    "wake_discarded_line_total": discarded_line,
                    "wake_discarded_slip_total": discarded_slip,
                    "accounted_line_content_total": accounted,
                    "line_content_balance_signed": balance,
                    "line_content_balance_error": abs(balance),
                    "line_content_balance_relative_error": abs(balance) / scale,
                    "advance_total_m": float(state.get("advance_total_m", 0.0) or 0.0),
                }
            )
    if not records:
        return None
    return {
        "schema": "final_persistent_line_content_accounting_v10_0_5_16",
        "legacy_N_em_semantics": "instantaneous_active_retained_line_content",
        "legacy_wake_discarded_total_includes_slip_history": True,
        "authoritative_discarded_line_field": "wake_discarded_line_total",
        "fronts": records,
        "maximum_relative_balance_error": max(
            row["line_content_balance_relative_error"] for row in records
        ),
    }


def _rewrite_summary_event_semantics(out: Path) -> None:
    summary_path = out / "summary.json"
    events_path = out / "stochastic_geometry_events_v10_0_5_16.json"
    if not summary_path.is_file() or not events_path.is_file():
        return
    try:
        summary = json.loads(summary_path.read_text())
        events = json.loads(events_path.read_text())
    except Exception:
        return
    accepted = [
        row
        for row in events
        if bool(row.get("inserted", False)) and float(row.get("moved_m", 0.0)) > 0.0
    ]
    if not isinstance(summary, list) or not summary or not accepted:
        return
    lengths = [float(row["moved_m"]) for row in accepted]
    nominal_values = [
        float(row.get("nominal_checkpoint_m", row.get("requested_fixed_length_m", 0.0)))
        for row in accepted
    ]
    nominal_values = [value for value in nominal_values if value > 0.0]
    if not nominal_values:
        return
    nominal = float(statistics.median(nominal_values))
    path_length = float(sum(lengths))
    first_p0 = np.asarray(accepted[0].get("p0_m", [0.0, 0.0]), dtype=float)
    last_p1 = np.asarray(
        accepted[-1].get(
            "actual_p1_m", accepted[-1].get("stochastic_p1_requested_m", first_p0)
        ),
        dtype=float,
    )
    row = summary[0]
    row.update(
        {
            "n_geometry_events": int(len(accepted)),
            "n_equivalent_checkpoints_exact": path_length / nominal,
            "n_equivalent_checkpoints_rounded": int(round(path_length / nominal)),
            "nominal_checkpoint_length_m": nominal,
            "geometry_path_length_m": path_length,
            "geometry_projected_extension_m": float(last_p1[0] - first_p0[0]),
            "n_advances_semantics": "rounded_path_length_over_nominal_checkpoint",
            "n_geometry_events_semantics": (
                "accepted_stochastic_cleavage_renewals_and_adaptive_CZM_commits"
            ),
        }
    )
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")


def _rewrite_release_metadata(
    out: Path | None, stochastic: dict[str, Any]
) -> None:
    if out is None or not out.exists():
        return
    _rewrite_summary_event_semantics(out)
    manifest = out / PRODUCTION_MANIFEST
    if manifest.is_file():
        payload = json.loads(manifest.read_text())
        payload["schema"] = "persistent_site_production_manifest_v10_0_5_16"
        payload["model"] = MODEL_ID
        payload["point_release"] = POINT_RELEASE
        physics = dict(payload.get("physics_contract", {}))
        physics.update(
            {
                "PF_reference_commit": PF_REFERENCE_COMMIT,
                "PF_local_update_map": PF_UPDATE_MAP,
                "cleavage_hazard_schema": HAZARD_SCHEMA,
                "cleavage_hazard_mode": stochastic["hazard_mode"],
                "cleavage_hazard_seed": stochastic["hazard_seed"],
                "cleavage_threshold_distribution": "exponential_unit_mean",
                "cleavage_event_length_schema": EVENT_LENGTH_SCHEMA,
                "cleavage_event_length_mode": stochastic["event_length_mode"],
                "cleavage_event_minimum_factor": stochastic[
                    "event_minimum_factor"
                ],
                "cleavage_event_maximum_factor": stochastic[
                    "event_maximum_factor"
                ],
                "event_length_uses_same_integrated_hazard_threshold": True,
                "mean_event_length_preserved": True,
                "noise_added_to_K": False,
                "noise_added_to_J": False,
                "noise_added_to_barriers": False,
                "geometry_realization": GEOMETRY_MODEL,
                "adaptive_czm_quality_veto_retained": True,
                "persistent_state_diagnostic_model": DIAGNOSTIC_MODEL,
                "refinement_audit_model": REFINEMENT_AUDIT_MODEL,
                "constitutive_parameters_changed_from_v10_0_5_15": False,
            }
        )
        payload["physics_contract"] = physics
        payload["front_engine"] = (
            PersistentSitePFStochasticMovingTipFrontEngineV100516.audit_payload()
        )
        payload["adaptive_czm_tip_support_audit"] = tip_support_audit_v1005144()
        payload["refinement_metadata_lifecycle_audit"] = (
            refinement_audit_payload_v100516()
        )
        geometry_path = out / "stochastic_geometry_events_v10_0_5_16.json"
        if geometry_path.is_file():
            events = json.loads(geometry_path.read_text())
            accepted = [row for row in events if bool(row.get("inserted", False))]
            payload["stochastic_geometry"] = {
                "schema": GEOMETRY_MODEL,
                "events_file": str(geometry_path),
                "n_attempted_events": len(events),
                "n_accepted_events": len(accepted),
                "event_lengths_m": [float(row.get("moved_m", 0.0)) for row in accepted],
                "thresholds": [
                    float(row.get("threshold_action", 0.0)) for row in accepted
                ],
            }
        accounting = _final_line_content_accounting(out)
        if accounting is not None:
            payload["final_line_content_accounting"] = accounting
        manifest.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

    old_selection = out / "persistent_site_parameter_selection_v10_0_5_14_1.json"
    if old_selection.is_file():
        selection = json.loads(old_selection.read_text())
        selection["schema"] = MODEL_ID
        selection["point_release"] = POINT_RELEASE
        policy = dict(selection.get("policy", {}))
        policy.update(
            {
                "PF_reference_commit": PF_REFERENCE_COMMIT,
                "PF_local_update_map": PF_UPDATE_MAP,
                "stochastic_hazard_schema": HAZARD_SCHEMA,
                "stochastic_hazard_mode": stochastic["hazard_mode"],
                "stochastic_hazard_seed": stochastic["hazard_seed"],
                "stochastic_event_length_schema": EVENT_LENGTH_SCHEMA,
                "stochastic_event_length_mode": stochastic["event_length_mode"],
                "event_length_uses_same_hazard_threshold": True,
                "mean_event_length_preserved": True,
                "persistent_state_diagnostic_model": DIAGNOSTIC_MODEL,
                "refinement_audit_model": REFINEMENT_AUDIT_MODEL,
                "candidate_source_or_shielding_closure_applied": True,
                "barrier_only_transfer": False,
            }
        )
        selection["policy"] = policy
        (out / SELECTION_MANIFEST).write_text(
            json.dumps(selection, indent=2, sort_keys=True) + "\n"
        )
        old_selection.unlink()


def main(argv: list[str] | None = None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    out = _out_path(user_args)
    stochastic = _require_stochastic_config()
    PersistentSitePFStochasticMovingTipFrontEngineV100516.configure_stochastic(
        **stochastic
    )
    saved = {
        "point_release": _base.POINT_RELEASE,
        "model_id": _base.MODEL_ID,
        "manifest": _base.PRODUCTION_MANIFEST,
        "engine_class": _base.PersistentSiteMovingProcessZoneFrontEngineV100514,
    }
    _base.POINT_RELEASE = POINT_RELEASE
    _base.MODEL_ID = MODEL_ID
    _base.PRODUCTION_MANIFEST = PRODUCTION_MANIFEST
    _base.PersistentSiteMovingProcessZoneFrontEngineV100514 = (
        PersistentSitePFStochasticMovingTipFrontEngineV100516
    )
    try:
        with (
            installed_refinement_audit_v100516(),
            installed_tip_support_repair_v1005144(),
            installed_persistent_diagnostics_v1005144(),
            installed_persistent_diagnostics_v100516(),
            installed_stochastic_geometry_v100516(),
        ):
            return _base.main(user_args)
    finally:
        if out is not None:
            write_stochastic_geometry_events(out)
        _rewrite_release_metadata(out, stochastic)
        _base.POINT_RELEASE = saved["point_release"]
        _base.MODEL_ID = saved["model_id"]
        _base.PRODUCTION_MANIFEST = saved["manifest"]
        _base.PersistentSiteMovingProcessZoneFrontEngineV100514 = saved[
            "engine_class"
        ]


if __name__ == "__main__":
    main()


__all__ = ["POINT_RELEASE", "MODEL_ID", "PRODUCTION_MANIFEST", "main"]
