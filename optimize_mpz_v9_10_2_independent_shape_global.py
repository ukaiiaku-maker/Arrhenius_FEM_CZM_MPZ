#!/usr/bin/env python3
"""Full Sobol global search with independent alpha/n for all four barriers.

The v9.10 formulation represented the originally intended hierarchy:
cleavage had its own EXP-floor shape while emission, Peierls, and Taylor shared
one shape.  The v9.10.1 audit tied all four shapes together.  Version 9.10.2
releases the complete shape space:

    (alpha_c, n_c), (alpha_e, n_e), (alpha_P, n_P), (alpha_T, n_T).

All other v9.10 unified mobile/retained physics is retained.  Each restart uses
a fresh full Sobol population and never starts from a previous shortlist.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np

import optimize_mpz_v9_10_unified_global as base
from arrhenius_fracture.emission_derived_plasticity import (
    CorrelatedTaylorConfig,
    EmissionDerivedPeierlsTaylorConfig,
    ExpFloorSurface,
)
from arrhenius_fracture.emission_derived_plasticity_v9102 import (
    EmissionDerivedPeierlsTaylorModel,
    IndependentShapeEntropyMechanismScale,
)


PARAMETER_NAMES = (
    "cleave_G00_eV",
    "cleave_gT_eV_per_K",
    "cleave_sigc0_GPa",
    "cleave_sT_GPa_per_K",
    "cleave_exp_a",
    "cleave_exp_n",
    "cleave_floor_frac",
    "emit_G00_eV",
    "emit_gT_eV_per_K",
    "emit_sigc0_GPa",
    "emit_sT_GPa_per_K",
    "emit_exp_a",
    "emit_exp_n",
    "emit_floor_frac",
    "peierls_H0_eV",
    "peierls_exp_a",
    "peierls_exp_n",
    "delta_H_PT_eV",
    "taylor_exp_a",
    "taylor_exp_n",
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

BOUNDS = dict(base.BOUNDS)
BOUNDS.update(
    {
        "peierls_exp_a": (0.02, 2.0),
        "peierls_exp_n": (0.45, 2.5),
        "taylor_exp_a": (0.02, 2.0),
        "taylor_exp_n": (0.45, 2.5),
    }
)


def bounds_list() -> list[tuple[float, float]]:
    return [BOUNDS[name] for name in PARAMETER_NAMES]


def decode(x: np.ndarray) -> dict[str, float]:
    p = {name: float(value) for name, value in zip(PARAMETER_NAMES, x)}
    p.update(
        {
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
            "independent_shape_all_four_active": 1.0,
        }
    )
    return p


def build_model(
    p: dict[str, float], Tref_K: float
) -> EmissionDerivedPeierlsTaylorModel:
    emit0 = max(p["emit_G00_eV"], 1.0e-12)
    parent = ExpFloorSurface(
        G00_eV=p["emit_G00_eV"],
        gT_eV_per_K=p["emit_gT_eV_per_K"],
        sigc0_Pa=p["emit_sigc0_GPa"] * 1.0e9,
        sT_Pa_per_K=p["emit_sT_GPa_per_K"] * 1.0e9,
        Tref_K=Tref_K,
        a=p["emit_exp_a"],
        n=p["emit_exp_n"],
        floor_fraction=p["emit_floor_frac"],
        floor_min_eV=1.0e-4,
        floor_max_fraction=0.95,
    )
    return EmissionDerivedPeierlsTaylorModel(
        EmissionDerivedPeierlsTaylorConfig(
            parent=parent,
            peierls=IndependentShapeEntropyMechanismScale(
                energy_ratio=p["peierls_H0_eV"] / emit0,
                activation_entropy_kB=p["peierls_activation_entropy_kB"],
                exp_a=p["peierls_exp_a"],
                exp_n=p["peierls_exp_n"],
                stress_ratio=1.0,
                rate_prefactor_s=p["peierls_nu0_s"],
            ),
            taylor=IndependentShapeEntropyMechanismScale(
                energy_ratio=p["taylor_H0_eV"] / emit0,
                activation_entropy_kB=p["taylor_activation_entropy_kB"],
                exp_a=p["taylor_exp_a"],
                exp_n=p["taylor_exp_n"],
                stress_ratio=1.0,
                rate_prefactor_s=p["taylor_nu0_s"],
            ),
            correlated_taylor=CorrelatedTaylorConfig(
                rho_c_m2=p["taylor_corr_rho_c_m2"],
                renewal_time_s=1.0,
                m_exponent=1.0,
                m_scale=p["taylor_corr_scale"],
                m_cap=float("inf"),
            ),
            mobile_saturation_density_m2=float("inf"),
            mobile_density_floor_m2=0.0,
            jump_length_min_m=0.0,
            taylor_phi_max=float("inf"),
            rate_cap_s=float("inf"),
        )
    )


def _argument_value(flag: str, default: str) -> str:
    try:
        index = sys.argv.index(flag)
    except ValueError:
        return default
    if index + 1 >= len(sys.argv):
        return default
    return sys.argv[index + 1]


def _mark_outputs() -> None:
    target_class = _argument_value("--target-class", "ceramic")
    outroot = Path(
        _argument_value(
            "--out", "runs/mpz_v9_10_2_independent_shape_global_search_v1"
        )
    )
    class_dir = outroot.resolve() / target_class
    summary_path = class_dir / "unified_global_summary.json"
    config_path = class_dir / "unified_global_config.json"
    if summary_path.exists():
        summary = json.loads(summary_path.read_text())
        summary.update(
            {
                "status": "V9_10_2_INDEPENDENT_SHAPE_GLOBAL_SEARCH_COMPLETE",
                "shape_mode": "INDEPENDENT_ALPHA_N_FOR_C_E_P_T",
                "parameter_count": len(PARAMETER_NAMES),
            }
        )
        summary_path.write_text(json.dumps(summary, indent=2))
    if config_path.exists():
        config = json.loads(config_path.read_text())
        config.update(
            {
                "shape_mode": "INDEPENDENT_ALPHA_N_FOR_C_E_P_T",
                "parameter_names": list(PARAMETER_NAMES),
                "bounds": BOUNDS,
                "derived_shape_equalities": [],
                "shape_coordinates": [
                    ["cleave_exp_a", "cleave_exp_n"],
                    ["emit_exp_a", "emit_exp_n"],
                    ["peierls_exp_a", "peierls_exp_n"],
                    ["taylor_exp_a", "taylor_exp_n"],
                ],
            }
        )
        config_path.write_text(json.dumps(config, indent=2))


def main() -> None:
    base.PARAMETER_NAMES = PARAMETER_NAMES
    base.BOUNDS = BOUNDS
    base.bounds_list = bounds_list
    base.decode = decode
    base.build_model = build_model
    base.main()
    _mark_outputs()


if __name__ == "__main__":
    main()
