"""v10.0.2 progressive Mode-I entry point with full event-time lifecycle."""
from __future__ import annotations

import json
from pathlib import Path
import sys

from . import mode_i_first_passage_v10_0 as _foundation
from . import mode_i_first_passage_v10_0_progressive as _base
from .kinetic_campaign_czm_v1001 import (
    ResetSafeDevelopedStateDiagnosticCZMFrontEngine,
)
from .kinetic_progressive_2d_v1002 import (
    build_progressive_run_2d_v1002,
    progressive_runtime_payload_v1002,
    reset_progressive_runtime_v1002,
    write_progressive_runtime_audit_v1002,
)

MODEL_ID = (
    "FEM_CZM_Mode_I_kinetic_campaign_czm_v10_0_2_"
    "progressive_event_lifecycle"
)


def _annotate(argv: list[str]) -> None:
    out_value = _foundation._option_value(argv, "--out")
    if out_value is None:
        return
    root = Path(out_value)
    path = root / "kinetic_campaign_czm_v10_0_audit.json"
    if not path.exists():
        return
    payload = json.loads(path.read_text())
    runtime = progressive_runtime_payload_v1002()
    payload.update({
        "model": MODEL_ID,
        "point_release": "10.0.2",
        "campaign_state_reset_safe": True,
        "rejected_damage_retry_integrated": True,
        "unused_event_time_carry_integrated": True,
        "same_load_re_equilibration_after_commit": True,
        "dot_ep_transactionally_restored": True,
        "event_lifecycle_runtime": runtime,
        "scope": "one_segment_progressive_validation_before_penalty_convergence",
        "long_progressive_runs_authorized": False,
    })
    path.write_text(json.dumps(payload, indent=2, default=str))


def main(argv: list[str] | None = None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    original_engine = _foundation.DevelopedStateDiagnosticCZMFrontEngine
    original_build = _base.build_progressive_run_2d
    original_payload = _base.progressive_runtime_payload
    original_reset = _base.reset_progressive_runtime
    original_write = _base.write_progressive_runtime_audit

    _foundation.DevelopedStateDiagnosticCZMFrontEngine = (
        ResetSafeDevelopedStateDiagnosticCZMFrontEngine
    )
    _base.build_progressive_run_2d = build_progressive_run_2d_v1002
    _base.progressive_runtime_payload = progressive_runtime_payload_v1002
    _base.reset_progressive_runtime = reset_progressive_runtime_v1002
    _base.write_progressive_runtime_audit = write_progressive_runtime_audit_v1002
    try:
        results = _base.main(user_args)
    finally:
        _foundation.DevelopedStateDiagnosticCZMFrontEngine = original_engine
        _base.build_progressive_run_2d = original_build
        _base.progressive_runtime_payload = original_payload
        _base.reset_progressive_runtime = original_reset
        _base.write_progressive_runtime_audit = original_write
    _annotate(user_args)
    return results


if __name__ == "__main__":
    main()
