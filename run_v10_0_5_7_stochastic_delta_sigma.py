#!/usr/bin/env python3
"""v10.0.5.7 stochastic delta-sigma runner with fixed-grip K audit."""
from __future__ import annotations

import os
import sys
from typing import Iterable

import run_v10_0_5_6_stochastic_delta_sigma as _base
from arrhenius_fracture.kj_audit_v10057 import (
    POINT_RELEASE,
    build_kj_audit_row,
    load_fixed_grip_reference,
)

_ORIGINAL_BUILD_PARSER = _base.build_parser
_REFERENCE_ENV = "ARRHENIUS_FIXED_GRIP_REFERENCE_JSON"


def build_parser():
    parser = _ORIGINAL_BUILD_PARSER()
    parser.description = __doc__
    parser.set_defaults(KJ_LEFM_ratio_min=0.85, KJ_LEFM_ratio_max=1.15)
    parser.add_argument(
        "--fixed-grip-reference-json",
        default=os.environ.get(_REFERENCE_ENV),
        help=(
            "convergence-approved fixed-grip reference artifact; may also be set "
            f"through {_REFERENCE_ENV}"
        ),
    )
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    args_list = list(sys.argv[1:] if argv is None else argv)
    parsed = build_parser().parse_args(args_list)
    if not parsed.fixed_grip_reference_json:
        raise SystemExit(
            "v10.0.5.7 requires --fixed-grip-reference-json PATH or "
            f"{_REFERENCE_ENV}; no built-in geometry factor is permitted"
        )
    reference = load_fixed_grip_reference(parsed.fixed_grip_reference_json)

    def configured_audit_row(**kwargs):
        return build_kj_audit_row(
            **kwargs,
            fixed_grip_reference=reference,
        )

    saved_parser = _base.build_parser
    saved_builder = _base.build_kj_audit_row
    saved_release = _base.POINT_RELEASE
    try:
        _base.build_parser = build_parser
        _base.build_kj_audit_row = configured_audit_row
        _base.POINT_RELEASE = POINT_RELEASE
        return int(_base.main(args_list) or 0)
    finally:
        _base.build_parser = saved_parser
        _base.build_kj_audit_row = saved_builder
        _base.POINT_RELEASE = saved_release


if __name__ == "__main__":
    raise SystemExit(main())
