#!/usr/bin/env python3
"""v10.0.5.7 stochastic delta-sigma runner with fixed-grip K audit."""
from __future__ import annotations

from typing import Iterable

import run_v10_0_5_6_stochastic_delta_sigma as _base
from arrhenius_fracture.kj_audit_v10057 import (
    POINT_RELEASE,
    build_kj_audit_row,
)

_ORIGINAL_BUILD_PARSER = _base.build_parser


def build_parser():
    parser = _ORIGINAL_BUILD_PARSER()
    parser.description = __doc__
    parser.set_defaults(KJ_LEFM_ratio_min=0.85, KJ_LEFM_ratio_max=1.15)
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    saved_parser = _base.build_parser
    saved_builder = _base.build_kj_audit_row
    saved_release = _base.POINT_RELEASE
    try:
        _base.build_parser = build_parser
        _base.build_kj_audit_row = build_kj_audit_row
        _base.POINT_RELEASE = POINT_RELEASE
        return int(_base.main(argv) or 0)
    finally:
        _base.build_parser = saved_parser
        _base.build_kj_audit_row = saved_builder
        _base.POINT_RELEASE = saved_release


if __name__ == "__main__":
    raise SystemExit(main())
