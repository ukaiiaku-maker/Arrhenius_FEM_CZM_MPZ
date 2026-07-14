#!/usr/bin/env python3
"""Spatial promotion for v9.10.2 independent four-barrier shapes."""
from __future__ import annotations

import json
from pathlib import Path
import sys

import promote_mpz_v9_10_spatial as base
from arrhenius_fracture.moving_process_zone_v9102 import (
    MovingProcessZoneState as IndependentShapeMPZState,
)


_BASE_RUN_SPATIAL = base.run_spatial


def run_spatial(candidate, T_K, opt, *, mpz_length_m, mpz_n_bins):
    """Inject candidate-specific Peierls/Taylor alpha and n into v9.10 MPZ."""
    original_state = base.MovingProcessZoneState

    def state_factory(cfg):
        cfg.pt_peierls_exp_a = base.finite(
            candidate, "peierls_exp_a", base.finite(candidate, "emit_exp_a", 0.2)
        )
        cfg.pt_peierls_exp_n = base.finite(
            candidate, "peierls_exp_n", base.finite(candidate, "emit_exp_n", 1.0)
        )
        cfg.pt_taylor_exp_a = base.finite(
            candidate, "taylor_exp_a", base.finite(candidate, "emit_exp_a", 0.2)
        )
        cfg.pt_taylor_exp_n = base.finite(
            candidate, "taylor_exp_n", base.finite(candidate, "emit_exp_n", 1.0)
        )
        return IndependentShapeMPZState(cfg)

    base.MovingProcessZoneState = state_factory
    try:
        result, engine, extra = _BASE_RUN_SPATIAL(
            candidate,
            T_K,
            opt,
            mpz_length_m=mpz_length_m,
            mpz_n_bins=mpz_n_bins,
        )
    finally:
        base.MovingProcessZoneState = original_state
    extra.update(
        {
            "independent_shape_all_four_active": 1.0,
            "peierls_exp_a": base.finite(candidate, "peierls_exp_a", 0.2),
            "peierls_exp_n": base.finite(candidate, "peierls_exp_n", 1.0),
            "taylor_exp_a": base.finite(candidate, "taylor_exp_a", 0.2),
            "taylor_exp_n": base.finite(candidate, "taylor_exp_n", 1.0),
        }
    )
    return result, engine, extra


def _argument_value(flag: str, default: str) -> str:
    try:
        index = sys.argv.index(flag)
    except ValueError:
        return default
    if index + 1 >= len(sys.argv):
        return default
    return sys.argv[index + 1]


def _mark_outputs() -> None:
    outroot = Path(
        _argument_value(
            "--out", "runs/mpz_v9_10_2_independent_shape_spatial_promotion_v1"
        )
    ).resolve()
    report = outroot / "unified_spatial_report.json"
    if report.exists():
        data = json.loads(report.read_text())
        data.update(
            {
                "status": "V9_10_2_INDEPENDENT_SHAPE_SPATIAL_COMPLETE_NOT_2D_VALIDATED",
                "shape_mode": "INDEPENDENT_ALPHA_N_FOR_C_E_P_T",
            }
        )
        report.write_text(json.dumps(data, indent=2))


def main() -> None:
    base.run_spatial = run_spatial
    base.main()
    _mark_outputs()


if __name__ == "__main__":
    main()
