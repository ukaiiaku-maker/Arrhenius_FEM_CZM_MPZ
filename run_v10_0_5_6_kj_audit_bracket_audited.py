#!/usr/bin/env python3
"""Audited launcher for the v10.0.5.6 KJ audit and bracket workflow."""
from __future__ import annotations

from pathlib import Path

import run_v10_0_5_6_kj_audit_bracket as _base


AUDITED_CAMPAIGN = (
    Path(__file__).resolve().parent
    / "run_v10_0_5_6_stochastic_delta_sigma_audited.py"
)


def _read_csv_with_booleans(path):
    rows = _ORIGINAL_READ_CSV(path)
    for row in rows:
        for key in (
            "first_passage_observed",
            "right_censored",
            "reached_cycle_horizon",
            "reached_target_extension",
        ):
            value = row.get(key)
            if isinstance(value, str):
                token = value.strip().lower()
                if token in {"true", "1", "yes"}:
                    row[key] = True
                elif token in {"false", "0", "no", ""}:
                    row[key] = False
    return rows


_ORIGINAL_READ_CSV = _base._read_csv


def main(argv=None) -> int:
    saved_campaign = _base.CAMPAIGN
    saved_reader = _base._read_csv
    _base.CAMPAIGN = AUDITED_CAMPAIGN
    _base._read_csv = _read_csv_with_booleans
    try:
        return int(_base.main(argv) or 0)
    finally:
        _base.CAMPAIGN = saved_campaign
        _base._read_csv = saved_reader


if __name__ == "__main__":
    raise SystemExit(main())
