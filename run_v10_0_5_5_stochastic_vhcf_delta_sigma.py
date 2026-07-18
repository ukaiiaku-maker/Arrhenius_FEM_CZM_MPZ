#!/usr/bin/env python3
"""v10.0.5.5 stochastic VHCF remote-stress campaign.

The existing v10.0.5.4 campaign owns calibration, physical-horizon accounting,
case summaries and plots.  This wrapper selects the stochastic v10.0.5.5 entry
point and completion manifest while preserving those audited campaign tools.
"""
from __future__ import annotations

import os

import run_v10_0_5_4_vhcf_delta_sigma as _base

POINT_RELEASE = "10.0.5.5"
ENTRY_MODULE = (
    "arrhenius_fracture.mode_i_first_passage_v10_0_5_5_stochastic_vhcf"
)
COMPLETION_MANIFEST = "run_completion_v10_0_5_5_stochastic_vhcf.json"


def main(argv=None) -> int:
    saved_entry = _base.ENTRY_MODULE
    saved_release = _base.POINT_RELEASE
    saved_completion = _base.COMPLETION_MANIFEST
    _base.ENTRY_MODULE = ENTRY_MODULE
    _base.POINT_RELEASE = POINT_RELEASE
    _base.COMPLETION_MANIFEST = COMPLETION_MANIFEST

    # Stochastic is the production default for this point release.  Deterministic
    # mode remains an explicit regression ablation via EVENT_STATISTICS.
    os.environ.setdefault("ARRHENIUS_EVENT_STATISTICS", "stochastic")
    os.environ.setdefault("ARRHENIUS_STOCHASTIC_EMISSION", "1")
    os.environ.setdefault("ARRHENIUS_STOCHASTIC_SEED", "1")
    os.environ.setdefault("ARRHENIUS_STOCHASTIC_BLOCKS", "1")
    os.environ.setdefault("ARRHENIUS_RARE_EVENT_TARGET", "0.25")
    os.environ.setdefault("ARRHENIUS_TAU_LEAP_TARGET", "3.0")
    os.environ.setdefault("ARRHENIUS_TAU_SWITCH_EXPECTED_EVENTS", "10.0")
    os.environ.setdefault("ARRHENIUS_VHCF_FEM_CACHE", "1")
    try:
        return int(_base.main(argv) or 0)
    finally:
        _base.ENTRY_MODULE = saved_entry
        _base.POINT_RELEASE = saved_release
        _base.COMPLETION_MANIFEST = saved_completion


if __name__ == "__main__":
    raise SystemExit(main())
