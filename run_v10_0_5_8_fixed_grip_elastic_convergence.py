#!/usr/bin/env python3
"""Run the v10.0.5.8 fixed-grip elastic FEM convergence audit."""
from __future__ import annotations

import sys

from arrhenius_fracture.fixed_grip_elastic_audit_v10058 import main


def _production_args(argv: list[str]) -> list[str]:
    """Supply increments that are resolvable on the default 10 um coarse mesh."""
    out = list(argv)
    if "--crack-increment-um" not in out:
        out.extend(["--crack-increment-um", "40 20 10"])
    return out


if __name__ == "__main__":
    raise SystemExit(main(_production_args(sys.argv[1:])))
