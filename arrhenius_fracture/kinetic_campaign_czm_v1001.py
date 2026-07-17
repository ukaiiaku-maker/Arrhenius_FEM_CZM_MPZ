"""v10.0.1 reset-safe PF-equivalent kinetic CZM front state.

The v10.0 campaign engine replaces the v9.11 state after construction, but it
inherits the v9.11 ``reset`` implementation.  Any later explicit reset would
therefore rebuild the wrong state class.  This point release preserves the
v9.11 clock/control reset and then reinstantiates the campaign-calibrated state
from the immutable PF manifest.
"""
from __future__ import annotations

from typing import Any

from .kinetic_campaign_czm import (
    CampaignKineticMPZState,
    DevelopedStateDiagnosticCZMFrontEngine,
)


class ResetSafeDevelopedStateDiagnosticCZMFrontEngine(
    DevelopedStateDiagnosticCZMFrontEngine
):
    """Campaign engine whose virgin-state reset cannot fall back to v9.11."""

    state_model_detail = (
        "pf_v10_1_7_1_campaign_calibrated_continuous_tip_reset_safe_v1001"
    )
    campaign_state_reset_safe = True

    def reset(self) -> None:
        # Preserve threshold-stream and propagation-control initialization from
        # the established v9.11 front engine, then replace only its MPZ object.
        super().reset()
        self.mpz_state = CampaignKineticMPZState(
            self.mpz_config,
            self.manifest,
            b=self.b,
            G_Pa=self.G,
            kinetic_cfg=self.kinetic_config,
        )
        self.N_em = 0.0
        self.B = 0.0
        self.a_adv = 0.0
        self.n_adv = 0
        self.W_emit = 0.0
        self.t = 0.0
        self.K_prev = None
        self._lambda_c_prev = None
        self._K_cleave_prev = None
        self.micro_advance_total_m = 0.0
        self.checkpoint_advance_total_m = 0.0
        self._last_pre_checkpoint_snapshot = None
        self._last_channels: dict[str, Any] = {}
        self._sync_compat()

    def audit_payload(self) -> dict[str, Any]:
        payload = super().audit_payload()
        payload.update({
            "campaign_state_reset_safe": True,
            "reset_state_class": type(self.mpz_state).__name__,
            "reset_state_model": getattr(self.mpz_state, "state_model", None),
        })
        return payload


__all__ = ["ResetSafeDevelopedStateDiagnosticCZMFrontEngine"]
