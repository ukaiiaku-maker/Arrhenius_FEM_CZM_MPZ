"""Reset-safe Stage-C entry point for PF-equivalent kinetic CZM v10.0.1.

This point release remains limited to the one-segment progressive smoke.  It
fixes campaign-state reset wiring but deliberately does not authorize longer
runs until rejected-damage retry and unused-time carry are integrated into the
outer geometry lifecycle.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

from . import mode_i_first_passage_v10_0 as _foundation
from . import mode_i_first_passage_v10_0_progressive as _base
from .kinetic_campaign_czm_v1001 import (
    ResetSafeDevelopedStateDiagnosticCZMFrontEngine,
)

MODEL_ID = (
    "FEM_CZM_Mode_I_kinetic_campaign_czm_v10_0_1_"
    "progressive_clock_linear_reset_safe"
)


def _annotate(argv: list[str]) -> None:
    out_value = _foundation._option_value(argv, "--out")
    if out_value is None:
        return
    root = Path(out_value)
    path = root / "kinetic_campaign_czm_v10_0_audit.json"
    if path.exists():
        payload = json.loads(path.read_text())
        payload.update({
            "model": MODEL_ID,
            "point_release": "10.0.1",
            "campaign_state_reset_safe": True,
            "scope": "single_segment_progressive_smoke_only",
            "rejected_damage_retry_integrated": False,
            "unused_event_time_carry_integrated": False,
            "long_progressive_runs_authorized": False,
        })
        path.write_text(json.dumps(payload, indent=2, default=str))


def main(argv: list[str] | None = None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    original = _foundation.DevelopedStateDiagnosticCZMFrontEngine
    _foundation.DevelopedStateDiagnosticCZMFrontEngine = (
        ResetSafeDevelopedStateDiagnosticCZMFrontEngine
    )
    try:
        results = _base.main(user_args)
    finally:
        _foundation.DevelopedStateDiagnosticCZMFrontEngine = original
    _annotate(user_args)
    return results


if __name__ == "__main__":
    main()
