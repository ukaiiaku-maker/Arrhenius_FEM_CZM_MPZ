#!/usr/bin/env python3
"""v10.0.5.7 contour audit and first-passage bracket using fixed grip."""
from __future__ import annotations

import os
from pathlib import Path
import sys

import run_v10_0_5_6_kj_audit_bracket as _base
from arrhenius_fracture.kj_audit_v10057 import (
    POINT_RELEASE,
    load_fixed_grip_reference,
    select_contour_plateau,
)

CAMPAIGN = Path(__file__).resolve().parent / "run_v10_0_5_7_stochastic_delta_sigma.py"
_REFERENCE_ENV = "ARRHENIUS_FIXED_GRIP_REFERENCE_JSON"
_ORIGINAL_BUILD_PARSER = _base.build_parser


def build_parser():
    parser = _ORIGINAL_BUILD_PARSER()
    parser.add_argument(
        "--fixed-grip-reference-json",
        default=os.environ.get(_REFERENCE_ENV),
        help=(
            "convergence-approved fixed-grip reference artifact; may also be set "
            f"through {_REFERENCE_ENV}"
        ),
    )
    return parser


def main(argv=None) -> int:
    args_list = list(sys.argv[1:] if argv is None else argv)
    parsed = build_parser().parse_args(args_list)
    if not parsed.fixed_grip_reference_json:
        raise SystemExit(
            "v10.0.5.7 requires --fixed-grip-reference-json PATH or "
            f"{_REFERENCE_ENV}; no built-in geometry factor is permitted"
        )
    reference_path = str(Path(parsed.fixed_grip_reference_json).expanduser().resolve())
    load_fixed_grip_reference(reference_path)

    saved_campaign = _base.CAMPAIGN
    saved_release = _base.POINT_RELEASE
    saved_selector = _base.select_contour_plateau
    saved_parser = _base.build_parser
    saved_env = os.environ.get(_REFERENCE_ENV)
    try:
        os.environ[_REFERENCE_ENV] = reference_path
        _base.CAMPAIGN = CAMPAIGN
        _base.POINT_RELEASE = POINT_RELEASE
        _base.select_contour_plateau = select_contour_plateau
        _base.build_parser = build_parser
        return int(_base.main(args_list) or 0)
    finally:
        _base.CAMPAIGN = saved_campaign
        _base.POINT_RELEASE = saved_release
        _base.select_contour_plateau = saved_selector
        _base.build_parser = saved_parser
        if saved_env is None:
            os.environ.pop(_REFERENCE_ENV, None)
        else:
            os.environ[_REFERENCE_ENV] = saved_env


if __name__ == "__main__":
    raise SystemExit(main())
