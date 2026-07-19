#!/usr/bin/env python3
"""v10.0.5.7 contour audit and first-passage bracket using fixed grip."""
from __future__ import annotations

from pathlib import Path

import run_v10_0_5_6_kj_audit_bracket as _base
from arrhenius_fracture.kj_audit_v10057 import (
    POINT_RELEASE,
    select_contour_plateau,
)

CAMPAIGN = Path(__file__).resolve().parent / "run_v10_0_5_7_stochastic_delta_sigma.py"


def main(argv=None) -> int:
    saved_campaign = _base.CAMPAIGN
    saved_release = _base.POINT_RELEASE
    saved_selector = _base.select_contour_plateau
    try:
        _base.CAMPAIGN = CAMPAIGN
        _base.POINT_RELEASE = POINT_RELEASE
        _base.select_contour_plateau = select_contour_plateau
        return int(_base.main(argv) or 0)
    finally:
        _base.CAMPAIGN = saved_campaign
        _base.POINT_RELEASE = saved_release
        _base.select_contour_plateau = saved_selector


if __name__ == "__main__":
    raise SystemExit(main())
