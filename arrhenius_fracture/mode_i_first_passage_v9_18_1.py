"""v9.18.1 guard against accidental renewal during an active trial event.

The inherited emission-only continuation path temporarily sets the current
cleavage prefactor and action to zero, but the first trapezoidal hazard update can
still inherit ``_lambda_c_prev`` from the preceding loading step.  If that
interpolated half-step crosses a renewal threshold while a cohesive event is
already active, v9.18 raises ``cannot defer a second renewal``.

This wrapper treats such a crossing as a numerical continuation artifact:
restore the exact pre-renewal MPZ/action/threshold snapshot and continue the
already active cohesive event.  It does not discard a physically admissible
renewal when no event is active.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys
from typing import Any

from . import mode_i_first_passage_v9_18 as _v918


class RenewalRollbackPersistentWakeController(
    _v918.PersistentWakeHazardController
):
    """Persistent-wake controller with transactional active-event renewal veto."""

    def __init__(self) -> None:
        super().__init__()
        self.active_event_renewals_rolled_back = 0
        self.active_event_thresholds_rolled_back = 0

    @staticmethod
    def _restore_accidental_renewal(engine: Any) -> bool:
        snapshot = getattr(engine, "_last_pre_renewal_event_snapshot", None)
        state = getattr(engine, "_last_pre_renewal_state", None)
        if not isinstance(snapshot, dict):
            return False

        if state is not None:
            engine.mpz_state = state.copy()
        engine.B = float(snapshot.get("B", engine.B))
        threshold = snapshot.get("threshold")
        stream = getattr(engine, "_threshold_stream", None)
        if threshold is not None and stream is not None:
            stream.restore(threshold)
        engine.a_adv = float(snapshot.get("a_adv", engine.a_adv))
        engine.n_adv = int(snapshot.get("n_adv", engine.n_adv))
        if hasattr(engine, "_reload_until_U_m"):
            engine._reload_until_U_m = snapshot.get("reload_until_U_m")
        if hasattr(engine, "_reload_until_K"):
            engine._reload_until_K = snapshot.get("reload_until_K")
        engine._sync_compat()
        engine._last_pre_renewal_state = None
        engine._last_pre_renewal_event_snapshot = None
        return True

    def defer_engine_renewal(
        self,
        engine: Any,
        nfire: int,
        distance_m: float,
        wake_preview: dict[str, float] | None = None,
    ) -> None:
        if self.active:
            restored = self._restore_accidental_renewal(engine)
            if not restored:
                raise RuntimeError(
                    "renewal crossed while a cohesive event was active, but the "
                    "transactional pre-renewal snapshot was unavailable"
                )
            self.active_event_renewals_rolled_back += max(int(nfire), 0)
            self.active_event_thresholds_rolled_back += 1
            if self._active_record is not None:
                self._active_record["active_event_renewals_rolled_back"] = int(
                    self._active_record.get(
                        "active_event_renewals_rolled_back", 0
                    )
                    + max(int(nfire), 0)
                )
                self._active_record[
                    "active_event_renewal_rollback_reason"
                ] = "inherited_previous_lambda_during_emission_only_continuation"
            return
        super().defer_engine_renewal(
            engine, nfire, distance_m, wake_preview
        )

    def payload(self) -> dict[str, Any]:
        data = super().payload()
        data.update({
            "schema": "persistent_plastic_wake_hazard_event_v9181_v1",
            "active_event_renewal_transactional_rollback_enabled": True,
            "active_event_renewals_rolled_back": int(
                self.active_event_renewals_rolled_back
            ),
            "active_event_thresholds_rolled_back": int(
                self.active_event_thresholds_rolled_back
            ),
        })
        return data


def main(argv=None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    original_controller = _v918.PersistentWakeHazardController
    _v918.PersistentWakeHazardController = RenewalRollbackPersistentWakeController
    try:
        results = _v918.main(user_args)
    finally:
        _v918.PersistentWakeHazardController = original_controller

    out_value = _v918._v917._v916._option_value(user_args, "--out")
    if out_value is not None:
        out = Path(out_value)
        source = out / "persistent_wake_event_relaxation_v918.json"
        if source.exists():
            payload = json.loads(source.read_text())
            payload["compatibility_source_filename"] = source.name
            (out / "persistent_wake_event_relaxation_v9181.json").write_text(
                json.dumps(payload, indent=2, default=str)
            )
    return results


if __name__ == "__main__":
    main()
