"""v10.0.5.18.3.2 production entry with action-error event horizons."""
from __future__ import annotations

import sys

from . import mode_i_first_passage_v10_0_5_18_3_2_four_class_joint_K_ramp_stochastic_emission as _base
from .persistent_site_joint_K_action_error_v10051832 import (
    PersistentSiteJointKActionErrorFrontEngineV10051832,
)


class ProductionJointKActionErrorFrontEngineV10051832(
    PersistentSiteJointKActionErrorFrontEngineV10051832
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


def main(argv: list[str] | None = None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    saved = _base.ProductionJointKRampStochasticEmissionFrontEngineV10051832
    _base.ProductionJointKRampStochasticEmissionFrontEngineV10051832 = (
        ProductionJointKActionErrorFrontEngineV10051832
    )
    try:
        return _base.main(user_args)
    finally:
        _base.ProductionJointKRampStochasticEmissionFrontEngineV10051832 = saved


if __name__ == "__main__":
    main()


__all__ = [
    "ProductionJointKActionErrorFrontEngineV10051832",
    "main",
]
