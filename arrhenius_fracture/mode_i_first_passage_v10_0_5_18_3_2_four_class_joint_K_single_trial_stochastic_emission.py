"""v10.0.5.18.3.2 final joint-K single-trial stochastic FEM/CZM entry."""
from __future__ import annotations

import json
from pathlib import Path
import sys

from . import mode_i_first_passage_v10_0_5_18_3_2_four_class_joint_K_ramp_stochastic_emission as _base
from .persistent_site_joint_K_single_trial_v10051832 import (
    SINGLE_TRIAL_SCHEMA,
    PersistentSiteJointKSingleTrialFrontEngineV10051832,
)


class ProductionJointKSingleTrialFrontEngineV10051832(
    PersistentSiteJointKSingleTrialFrontEngineV10051832
):
    """Production adapter retaining the instantaneous aggregate hazard audit."""

    def step_drives(self, K_cleave, K_emit, T, dt, metadata=None):
        out = super().step_drives(
            K_cleave,
            K_emit,
            T,
            dt,
            metadata=metadata,
        )
        last = dict(self.mpz_state.last_emission or {})
        aggregate = last.get(
            "aggregate_hazard_final_by_system_s",
            last.get("aggregate_hazard_initial_by_system_s", [0.0, 0.0]),
        )
        out["lambda_e"] = float(sum(float(value) for value in aggregate))
        out["lambda_e_semantics"] = (
            "instantaneous_sum_exact_stochastic_signed_channel_aggregate_hazards"
        )
        return out


def _update_json(path: Path, update) -> None:
    if not path.is_file():
        return
    payload = json.loads(path.read_text())
    update(payload)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n")


def _augment_single_trial_metadata(out: Path | None) -> None:
    if out is None:
        return

    fields = {
        "single_trial_event_horizon_active": True,
        "single_trial_event_horizon_schema": SINGLE_TRIAL_SCHEMA,
        "state_dependent_log_rate_rejection_active": False,
        "log_rate_change_is_diagnostic_only": True,
        "PF_transport_operator_preserved": True,
        "PF_transport_recursive_composition_error_control": False,
        "one_endpoint_transport_trial_per_event_horizon": True,
        "one_committed_transport_update_to_localized_event": True,
        "emission_threshold_crossings_localized_individually": True,
        "emission_events_batched": False,
    }

    def production(payload: dict) -> None:
        physics = dict(payload.get("physics_contract", {}) or {})
        physics.update(fields)
        payload["physics_contract"] = physics
        payload["front_engine"] = (
            ProductionJointKSingleTrialFrontEngineV10051832.audit_payload()
        )

    def selection(payload: dict) -> None:
        policy = dict(payload.get("policy", {}) or {})
        policy.update(fields)
        payload["policy"] = policy

    def transfer(payload: dict) -> None:
        payload.update(fields)

    _update_json(out / _base.PRODUCTION_MANIFEST, production)
    _update_json(out / _base.SELECTION_MANIFEST, selection)
    _update_json(out / _base.TRANSFER_MANIFEST, transfer)


def main(argv: list[str] | None = None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    saved_engine = _base.ProductionJointKRampStochasticEmissionFrontEngineV10051832
    saved_rewrite = _base._rewrite_point_release_metadata

    def rewrite(original_rewrite, out, audit, solver_args):
        saved_rewrite(original_rewrite, out, audit, solver_args)
        _augment_single_trial_metadata(out)

    _base.ProductionJointKRampStochasticEmissionFrontEngineV10051832 = (
        ProductionJointKSingleTrialFrontEngineV10051832
    )
    _base._rewrite_point_release_metadata = rewrite
    try:
        return _base.main(user_args)
    finally:
        _base.ProductionJointKRampStochasticEmissionFrontEngineV10051832 = (
            saved_engine
        )
        _base._rewrite_point_release_metadata = saved_rewrite


if __name__ == "__main__":
    main()


__all__ = [
    "ProductionJointKSingleTrialFrontEngineV10051832",
    "main",
]
