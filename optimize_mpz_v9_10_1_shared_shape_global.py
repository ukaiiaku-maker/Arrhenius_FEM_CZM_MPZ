#!/usr/bin/env python3
"""Full-space v9.10 search with one EXP-floor shape for all four barriers.

The v9.10 optimizer allowed independent alpha/n values for cleavage and
emission, while Peierls and Taylor inherited the emission shape.  This wrapper
implements the lower-dimensional physical hypothesis requested for v9.10.1:

    alpha_c = alpha_e = alpha_P = alpha_T = alpha_shared
    n_c     = n_e     = n_P     = n_T     = n_shared

The common shape is optimized globally.  Cleavage, emission, Peierls, and
Taylor retain independent barrier heights and thermal coordinates.  Peierls
and Taylor continue to satisfy H_P < H_T through
H_T = H_P + Delta H_PT, Delta H_PT > 0.

Every restart still begins from a fresh full Sobol population.  No prior
shortlist or class-specific parameter down-selection is used.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np

import optimize_mpz_v9_10_unified_global as base


PARAMETER_NAMES = (
    "cleave_G00_eV",
    "cleave_gT_eV_per_K",
    "cleave_sigc0_GPa",
    "cleave_sT_GPa_per_K",
    "cleave_floor_frac",
    "emit_G00_eV",
    "emit_gT_eV_per_K",
    "emit_sigc0_GPa",
    "emit_sT_GPa_per_K",
    "emit_floor_frac",
    "shared_exp_a",
    "shared_exp_n",
    "peierls_H0_eV",
    "delta_H_PT_eV",
    "peierls_activation_entropy_kB",
    "taylor_activation_entropy_kB",
    "log10_taylor_corr_rho_c_m2",
    "log10_taylor_corr_scale",
    "log10_source_sites_per_system",
    "log10_encounter_efficiency",
    "log10_retained_recovery_rate_s",
    "log10_source_refresh_length_um",
    "c_blunt",
)

BOUNDS = {
    name: base.BOUNDS[name]
    for name in PARAMETER_NAMES
    if name in base.BOUNDS
}
BOUNDS.update(
    {
        "shared_exp_a": (0.02, 2.0),
        "shared_exp_n": (0.45, 2.5),
    }
)


def bounds_list() -> list[tuple[float, float]]:
    return [BOUNDS[name] for name in PARAMETER_NAMES]


def decode(x: np.ndarray) -> dict[str, float]:
    p = {name: float(value) for name, value in zip(PARAMETER_NAMES, x)}
    alpha = p["shared_exp_a"]
    exponent = p["shared_exp_n"]
    p.update(
        {
            "cleave_exp_a": alpha,
            "cleave_exp_n": exponent,
            "emit_exp_a": alpha,
            "emit_exp_n": exponent,
            "peierls_exp_a": alpha,
            "peierls_exp_n": exponent,
            "taylor_exp_a": alpha,
            "taylor_exp_n": exponent,
            "taylor_H0_eV": p["peierls_H0_eV"] + p["delta_H_PT_eV"],
            "taylor_corr_rho_c_m2": 10.0
            ** p["log10_taylor_corr_rho_c_m2"],
            "taylor_corr_scale": 10.0 ** p["log10_taylor_corr_scale"],
            "source_sites_per_system": 10.0
            ** p["log10_source_sites_per_system"],
            "encounter_efficiency": 10.0 ** p["log10_encounter_efficiency"],
            "retained_recovery_rate_s": 10.0
            ** p["log10_retained_recovery_rate_s"],
            "source_refresh_length_um": 10.0
            ** p["log10_source_refresh_length_um"],
            "peierls_nu0_s": 1.0e12,
            "taylor_nu0_s": 1.0e11,
            "shared_shape_all_four_active": 1.0,
        }
    )
    return p


def _argument_value(flag: str, default: str) -> str:
    try:
        index = sys.argv.index(flag)
    except ValueError:
        return default
    if index + 1 >= len(sys.argv):
        return default
    return sys.argv[index + 1]


def _mark_v9101_outputs() -> None:
    target_class = _argument_value("--target-class", "ceramic")
    outroot = Path(
        _argument_value(
            "--out", "runs/mpz_v9_10_1_shared_shape_global_search_v1"
        )
    )
    class_dir = outroot.resolve() / target_class
    summary_path = class_dir / "unified_global_summary.json"
    config_path = class_dir / "unified_global_config.json"
    if summary_path.exists():
        summary = json.loads(summary_path.read_text())
        summary.update(
            {
                "status": "V9_10_1_SHARED_SHAPE_GLOBAL_SEARCH_COMPLETE",
                "shape_mode": "ONE_COMMON_ALPHA_N_FOR_ALL_FOUR_BARRIERS",
                "parameter_count": len(PARAMETER_NAMES),
            }
        )
        summary_path.write_text(json.dumps(summary, indent=2))
    if config_path.exists():
        config = json.loads(config_path.read_text())
        config.update(
            {
                "shape_mode": "ONE_COMMON_ALPHA_N_FOR_ALL_FOUR_BARRIERS",
                "parameter_names": list(PARAMETER_NAMES),
                "bounds": BOUNDS,
                "derived_shape_equalities": [
                    "cleave_exp_a = emit_exp_a = peierls_exp_a = taylor_exp_a = shared_exp_a",
                    "cleave_exp_n = emit_exp_n = peierls_exp_n = taylor_exp_n = shared_exp_n",
                ],
            }
        )
        config_path.write_text(json.dumps(config, indent=2))


def main() -> None:
    # UnifiedObjective, simulate_zero_d_rcurve, and base.main resolve these names
    # from the base module at evaluation time.  Patching them here preserves the
    # tested v9.10 unified transport/retention implementation while changing only
    # the global parameterization.
    base.PARAMETER_NAMES = PARAMETER_NAMES
    base.BOUNDS = BOUNDS
    base.bounds_list = bounds_list
    base.decode = decode
    base.main()
    _mark_v9101_outputs()


if __name__ == "__main__":
    main()
