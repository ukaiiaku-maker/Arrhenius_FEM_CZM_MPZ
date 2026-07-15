#!/usr/bin/env python3
"""Write raw and cascade-aware v9.11 propagation-event tables."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from arrhenius_fracture.rcurve_postprocess_v911 import write_cascade_aware_outputs


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("case_dir", type=Path)
    p.add_argument("--T-K", type=float, required=True)
    p.add_argument("--relative-load-tolerance", type=float, default=1.0e-4)
    p.add_argument("--absolute-load-tolerance-m", type=float, default=1.0e-12)
    args = p.parse_args()
    metrics = write_cascade_aware_outputs(
        args.case_dir,
        args.T_K,
        relative_load_tolerance=args.relative_load_tolerance,
        absolute_load_tolerance_m=args.absolute_load_tolerance_m,
    )
    print(json.dumps(metrics, indent=2, default=str))


if __name__ == "__main__":
    main()
