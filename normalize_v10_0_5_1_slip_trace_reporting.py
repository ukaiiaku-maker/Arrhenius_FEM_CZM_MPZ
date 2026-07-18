#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
import sys

from arrhenius_fracture.slip_trace_reporting_v10051 import normalize_output


def main(argv: list[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    if len(args) != 1:
        raise SystemExit(
            "usage: normalize_v10_0_5_1_slip_trace_reporting.py <v10.0.5-output-dir>"
        )
    payload = normalize_output(Path(args[0]))
    print("V10.0.5.1 REDUCED 2-D SLIP-TRACE REPORTING NORMALIZATION PASSED")
    print(json.dumps({
        "slip_trace_channel_count": payload["slip_trace_channel_count"],
        "channel_rows_written": payload["channel_rows_written"],
        "emission_observed_in_this_run": payload["emission_observed_in_this_run"],
        "emission_observation_required": payload[
            "emission_observation_required_for_implementation_certification"
        ],
        "physics_recomputed": payload["physics_recomputed"],
        "source_outputs_modified": payload["source_outputs_modified"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
