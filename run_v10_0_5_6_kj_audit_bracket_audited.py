#!/usr/bin/env python3
"""Audited launcher for the v10.0.5.6 KJ audit and bracket workflow."""
from __future__ import annotations

from pathlib import Path

import run_v10_0_5_6_kj_audit_bracket as _base


AUDITED_CAMPAIGN = (
    Path(__file__).resolve().parent
    / "run_v10_0_5_6_stochastic_delta_sigma_audited.py"
)


def main(argv=None) -> int:
    saved = _base.CAMPAIGN
    _base.CAMPAIGN = AUDITED_CAMPAIGN
    try:
        return int(_base.main(argv) or 0)
    finally:
        _base.CAMPAIGN = saved


if __name__ == "__main__":
    raise SystemExit(main())
