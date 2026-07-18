#!/usr/bin/env python3
"""Audited v10.0.5.6 runner with deterministic-calibration labels."""
from __future__ import annotations

from typing import Any, Mapping, Sequence

import run_v10_0_5_6_stochastic_delta_sigma as _base
from arrhenius_fracture.kj_audit_v10056 import LEGACY_LIMITER_LABELS


def _enrich_with_calibration_fallback(
    rows: Sequence[Mapping[str, Any]],
    scheduler_records: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    if scheduler_records:
        return _ORIGINAL_ENRICH(rows, scheduler_records)
    output = []
    for raw in rows:
        row = dict(raw)
        code = int(float(row.get("cycle_limiter_code", -1)))
        label = LEGACY_LIMITER_LABELS.get(code, f"unknown_code_{code}")
        row.update(
            {
                "stochastic_scheduler_mode": "deterministic_calibration",
                "stochastic_event_rate_per_cycle": 0.0,
                "stochastic_expected_state_events": 0.0,
                "base_cycle_limiter": label,
                "final_cycle_limiter": label,
                "cycle_limiter_label": label,
            }
        )
        output.append(row)
    return output


_ORIGINAL_ENRICH = _base.enrich_stochastic_block_rows


def main(argv=None) -> int:
    saved = _base.enrich_stochastic_block_rows
    _base.enrich_stochastic_block_rows = _enrich_with_calibration_fallback
    try:
        return int(_base.main(argv) or 0)
    finally:
        _base.enrich_stochastic_block_rows = saved


if __name__ == "__main__":
    raise SystemExit(main())
