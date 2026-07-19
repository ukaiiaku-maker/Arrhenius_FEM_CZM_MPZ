#!/usr/bin/env python3
"""Full-field v10.0.5.7 material R-curve campaign entry point.

Composes the repaired v10.0.5.7 summary/publication gate with the validated
v9.12 command substitution that selects the mapped MPZ field renderer.
"""
from __future__ import annotations

import run_mpz_v10_0_5_7_tip_only_material_rcurve as _repaired
import run_mpz_v9_12_tip_only_material_rcurve_fullfield as _fullfield

POINT_RELEASE = "10.0.5.7"


def main() -> None:
    campaign = _repaired._base
    saved_build_command = campaign.build_command
    campaign.build_command = _fullfield._build_command_fullfield
    try:
        _repaired.main()
    finally:
        campaign.build_command = saved_build_command


if __name__ == "__main__":
    main()


__all__ = ["POINT_RELEASE", "main"]
