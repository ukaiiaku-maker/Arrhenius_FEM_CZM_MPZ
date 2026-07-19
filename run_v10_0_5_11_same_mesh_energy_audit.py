#!/usr/bin/env python3
"""Run v10.0.5.11 same-production-mesh J versus fixed-grip energy audit."""
from __future__ import annotations

from typing import Iterable

import run_v10_0_5_10_refinement_support_audit as _base
from arrhenius_fracture.production_j_same_mesh_energy_v100511 import (
    CONTOUR_CSV,
    LAUNCH_FAILURE_JSON,
    PROBE_CSV,
    PROBE_JSON,
    RADIAL_CSV,
    SUMMARY_JSON,
    analyze_same_mesh_energy_v100511,
)

ENTRY_MODULE = "arrhenius_fracture.mode_i_first_passage_v10_0_5_11_same_mesh_probe"


def main(argv: Iterable[str] | None = None) -> int:
    saved = {
        "ENTRY_MODULE": _base.ENTRY_MODULE,
        "PROBE_JSON": _base.PROBE_JSON,
        "SUMMARY_JSON": _base.SUMMARY_JSON,
        "PROBE_CSV": _base.PROBE_CSV,
        "CONTOUR_CSV": _base.CONTOUR_CSV,
        "RADIAL_CSV": _base.RADIAL_CSV,
        "LAUNCH_FAILURE_JSON": _base.LAUNCH_FAILURE_JSON,
        "analyze_refinement_support_v100510": _base.analyze_refinement_support_v100510,
    }
    _base.ENTRY_MODULE = ENTRY_MODULE
    _base.PROBE_JSON = PROBE_JSON
    _base.SUMMARY_JSON = SUMMARY_JSON
    _base.PROBE_CSV = PROBE_CSV
    _base.CONTOUR_CSV = CONTOUR_CSV
    _base.RADIAL_CSV = RADIAL_CSV
    _base.LAUNCH_FAILURE_JSON = LAUNCH_FAILURE_JSON
    _base.analyze_refinement_support_v100510 = analyze_same_mesh_energy_v100511
    try:
        return _base.main(argv)
    finally:
        for name, value in saved.items():
            setattr(_base, name, value)


if __name__ == "__main__":
    raise SystemExit(main())
