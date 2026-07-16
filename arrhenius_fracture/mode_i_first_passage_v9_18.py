"""v9.18 Mode-I absolute-hazard opening with a persistent plastic wake.

This wrapper preserves the v9.17.1 one-fire routing and v9.17 absolute hazard
clock, replaces only the MPZ state instantiated by the v9.11 front engine, and
separates the nominal loading interval from adaptive event-localization dt.
"""
from __future__ import annotations

import json
import math
import os
from pathlib import Path
import sys
from typing import Any

from . import mode_i_first_passage_v9_17 as _v917
from . import mode_i_first_passage_v9_17_1 as _v9171
from . import mpz_front_engine_v911 as _front_module
from .moving_process_zone_v918 import MovingProcessZoneState as PersistentWakeState


class PersistentWakeHazardController(_v917.HazardClockTrialEventController):
    """Absolute-hazard controller with nominal-dt and wake-aware event records."""

    def __init__(self) -> None:
        super().__init__()
        self.nominal_loading_dt_s = max(
            float(os.environ.get("ARRHENIUS_NOMINAL_LOADING_DT_S", "8.4")),
            self.min_event_dt_s,
        )
        target_um = float(
            os.environ.get("ARRHENIUS_COMMITTED_TARGET_EXTENSION_UM", "inf")
        )
        self.committed_target_m = (
            target_um * 1.0e-6 if math.isfinite(target_um) and target_um > 0.0
            else math.inf
        )
        self.committed_target_reached = False
        self.post_target_renewals_suppressed = 0

    def _hold_cap_s(self) -> float:
        caps = [self.max_fixed_hold_s, self.nominal_loading_dt_s]
        finite = [x for x in caps if math.isfinite(x) and x > 0.0]
        return min(finite) if finite else math.inf

    def _state_snapshot(self) -> dict[str, float | None]:
        eng = self.active_engine
        state = getattr(eng, "mpz_state", None)
        if eng is None or state is None:
            return {}
        try:
            active_K = float(state.active_K_shielding(eng.G, eng.nu, eng.b))
            wake_K = float(state.wake_K_shielding(eng.G, eng.nu, eng.b))
        except Exception:
            active_K = wake_K = None
        return {
            "active_mobile_count": float(getattr(state, "mobile_count", 0.0)),
            "active_retained_count": float(getattr(state, "retained_count", 0.0)),
            "active_slip_count": float(state.local_slip_count()),
            "wake_mobile_count": float(getattr(state, "wake_mobile_count", 0.0)),
            "wake_retained_count": float(getattr(state, "wake_retained_count", 0.0)),
            "wake_slip_count": float(getattr(state, "wake_slip_count", 0.0)),
            "active_K_shield_Pa_sqrt_m": active_K,
            "wake_K_shield_Pa_sqrt_m": wake_K,
            "total_K_shield_Pa_sqrt_m": (
                None if active_K is None or wake_K is None else active_K + wake_K
            ),
        }

    def schedule_event(self, out: dict[str, Any]) -> int:
        event_id = super().schedule_event(out)
        if self._active_record is not None:
            self._active_record.update({
                f"{key}_at_nucleation": value
                for key, value in self._state_snapshot().items()
            })
            self._active_record["nominal_loading_dt_s"] = float(
                self.nominal_loading_dt_s
            )
            self._active_record["adaptive_dt_not_used_as_hold_cap"] = True
        return event_id

    def _commit_deferred_renewal(self) -> dict[str, float]:
        eng = self.active_engine
        before = self._state_snapshot()
        wake = super()._commit_deferred_renewal()
        after = self._state_snapshot()
        for key, value in before.items():
            if value is not None:
                wake[f"{key}_precommit"] = float(value)
        for key, value in after.items():
            if value is not None:
                wake[f"{key}_postcommit"] = float(value)

        if self.total_committed_distance_m + 1.0e-15 >= self.committed_target_m:
            self.committed_target_reached = True
            if eng is not None:
                eng.f.max_advances_per_step = 0.0
                eng.v918_committed_target_reached = True
                eng.v918_committed_target_m = float(self.committed_target_m)
            wake["committed_target_reached"] = 1.0
        else:
            wake["committed_target_reached"] = 0.0
        return wake

    def finish_substep(
        self,
        out: dict[str, Any],
        *,
        lambda_c_current: float,
        lambda_c_raw_current: float,
        K_cleave: float,
        T: float,
    ) -> None:
        rec = self._active_record
        precommit = self._state_snapshot()
        super().finish_substep(
            out,
            lambda_c_current=lambda_c_current,
            lambda_c_raw_current=lambda_c_raw_current,
            K_cleave=K_cleave,
            T=T,
        )
        if rec is not None and rec.get("substeps"):
            row = rec["substeps"][-1]
            row.update({
                f"{key}_before_commit": value
                for key, value in precommit.items()
            })
        if rec is not None and rec.get("status") == "complete_committed":
            wake = dict(rec.get("wake_on_commit", {}))
            rec["active_retained_precommit"] = wake.get(
                "active_retained_count_precommit"
            )
            rec["active_retained_postcommit"] = wake.get(
                "active_retained_count_postcommit"
            )
            rec["wake_retained_postcommit"] = wake.get(
                "wake_retained_count_postcommit"
            )
            rec["wake_K_shield_postcommit_Pa_sqrt_m"] = wake.get(
                "wake_K_shield_Pa_sqrt_m_postcommit"
            )
            rec["total_K_shield_postcommit_Pa_sqrt_m"] = wake.get(
                "total_K_shield_Pa_sqrt_m_postcommit"
            )
            rec["persistent_wake_state_committed"] = True

    def defer_engine_renewal(
        self,
        engine: Any,
        nfire: int,
        distance_m: float,
        wake_preview: dict[str, float] | None = None,
    ) -> None:
        if self.committed_target_reached:
            self.post_target_renewals_suppressed += max(int(nfire), 0)
            engine.f.max_advances_per_step = 0.0
            return
        super().defer_engine_renewal(engine, nfire, distance_m, wake_preview)

    def payload(self) -> dict[str, Any]:
        data = super().payload()
        data.update({
            "schema": "persistent_plastic_wake_hazard_event_v918_v1",
            "persistent_wake_enabled": True,
            "wake_toughening_mechanism": (
                "signed_dislocation_residual_stress_only; no bridging or transformation law"
            ),
            "nominal_loading_dt_s": float(self.nominal_loading_dt_s),
            "adaptive_dt_used_as_hold_cap": False,
            "committed_target_m": (
                None if not math.isfinite(self.committed_target_m)
                else float(self.committed_target_m)
            ),
            "committed_target_reached": bool(self.committed_target_reached),
            "post_target_renewals_suppressed": int(
                self.post_target_renewals_suppressed
            ),
        })
        return data


def main(argv=None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    original_state = _front_module.MovingProcessZoneState
    original_controller = _v917.HazardClockTrialEventController
    _front_module.MovingProcessZoneState = PersistentWakeState
    _v917.HazardClockTrialEventController = PersistentWakeHazardController
    try:
        results = _v9171.main(user_args)
    finally:
        _front_module.MovingProcessZoneState = original_state
        _v917.HazardClockTrialEventController = original_controller

    out_value = _v917._v916._option_value(user_args, "--out")
    if out_value is not None:
        out = Path(out_value)
        source = out / "absolute_hazard_event_relaxation_v917.json"
        if source.exists():
            payload = json.loads(source.read_text())
            payload["compatibility_source_filename"] = source.name
            (out / "persistent_wake_event_relaxation_v918.json").write_text(
                json.dumps(payload, indent=2, default=str)
            )
    return results


if __name__ == "__main__":
    main()
