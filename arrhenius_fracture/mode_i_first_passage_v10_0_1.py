"""Reset-safe Stage-B entry point for PF-equivalent kinetic CZM v10.0.1."""
from __future__ import annotations

import json
from pathlib import Path
import sys

from . import mode_i_first_passage_v10_0 as _base
from .kinetic_campaign_czm_v1001 import (
    ResetSafeDevelopedStateDiagnosticCZMFrontEngine,
)

MODEL_ID = "FEM_CZM_Mode_I_kinetic_campaign_czm_v10_0_1_reset_safe"


def _annotate(argv: list[str]) -> None:
    out_value = _base._option_value(argv, "--out")
    if out_value is None:
        return
    path = Path(out_value) / "kinetic_campaign_czm_v10_0_audit.json"
    if not path.exists():
        return
    payload = json.loads(path.read_text())
    payload.update({
        "model": MODEL_ID,
        "point_release": "10.0.1",
        "campaign_state_reset_safe": True,
        "long_progressive_runs_authorized": False,
    })
    path.write_text(json.dumps(payload, indent=2, default=str))


def main(argv: list[str] | None = None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    original = _base.DevelopedStateDiagnosticCZMFrontEngine
    _base.DevelopedStateDiagnosticCZMFrontEngine = (
        ResetSafeDevelopedStateDiagnosticCZMFrontEngine
    )
    try:
        results = _base.main(user_args)
    finally:
        _base.DevelopedStateDiagnosticCZMFrontEngine = original
    _annotate(user_args)
    return results


if __name__ == "__main__":
    main()
