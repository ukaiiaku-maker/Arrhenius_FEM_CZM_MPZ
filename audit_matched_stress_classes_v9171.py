#!/usr/bin/env python3
"""Run the v9.17 matched-stress audit with the actual v9.11 MPZ state.

The original v9.17 audit imported the package-level legacy
``moving_process_zone.MovingProcessZoneState``.  The FEM integration instead
uses ``moving_process_zone_v911.MovingProcessZoneState``, which inherits the
v9.10.2 unified Peierls encounter/Taylor release kinetics and independent
barrier shapes.  This wrapper makes the constitutive preflight use the same
state class as the FEM solver.
"""
from __future__ import annotations

import audit_matched_stress_classes_v917 as _base
from arrhenius_fracture.moving_process_zone_v911 import MovingProcessZoneState


def main():
    original = _base.MovingProcessZoneState
    _base.MovingProcessZoneState = MovingProcessZoneState
    try:
        return _base.main()
    finally:
        _base.MovingProcessZoneState = original


if __name__ == "__main__":
    main()
