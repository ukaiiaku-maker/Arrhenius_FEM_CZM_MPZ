#!/usr/bin/env python3
"""Fast constitutive verification for the full v9.11 integration."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from arrhenius_fracture.bulk_plasticity_v9102 import independent_config_from_dislocation_config
from arrhenius_fracture.emission_derived_plasticity_v9102 import EmissionDerivedPeierlsTaylorModel
from arrhenius_fracture.moving_process_zone_v911 import MovingProcessZoneState
from arrhenius_fracture.mpz_parameterization_v911 import (
    apply_pt_dislocation_config,
    build_mpz_config,
    load_selected_row,
)
from arrhenius_fracture.process_zone_2d_v911 import ProcessZone2DProfile


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--parameter-root", type=Path, default=Path("mpz_v9_11_parameters"))
    args = ap.parse_args()
    rows = []
    for cls in ("ceramic", "weakT", "DBTT"):
        row = load_selected_row(args.parameter_root / cls / "spatial_promotion_manifest.csv", cls)
        disl = SimpleNamespace()
        apply_pt_dislocation_config(disl, row)
        cfg = independent_config_from_dislocation_config(disl)
        model = EmissionDerivedPeierlsTaylorModel(cfg)
        margin = float("inf")
        finite = True
        for T in (300.0, 700.0, 900.0, 1200.0):
            for sigma in np.linspace(0.0, 30.0e9, 31):
                gp = float(np.asarray(model.barrier_eV("peierls", sigma, T)))
                gt = float(np.asarray(model.barrier_eV("taylor", sigma, T)))
                finite = finite and np.isfinite(gp) and np.isfinite(gt) and gp > 0 and gt > 0
                margin = min(margin, gt - gp)

        ns = SimpleNamespace(mpz_length_um=100.0, mpz_n_bins=200, r_pz=1.0e-6)
        mpz_cfg = build_mpz_config(ns, row)
        state = MovingProcessZoneState(mpz_cfg)
        state.retained[0, 0] = 1.0
        Ksh = state.shielding_K(160.0e9, 0.28, 2.74e-10)
        profile = ProcessZone2DProfile(
            state.x.copy(),
            np.full(state.n_bins, 1.0e14),
            np.linspace(1.0, 0.2, state.n_bins),
            True, 1.0, state.n_bins, "verification",
        )
        state.set_2d_profile(profile)
        local_rho = state.local_forest_density_m2()
        local_stress = state.local_stress_profile_Pa(5.0e9)
        passed = bool(
            finite and margin >= -1.0e-8 and Ksh > 0.0 and
            np.all(local_rho >= 1.0e14) and
            np.isclose(local_stress[0], 5.0e9) and
            local_stress[-1] < local_stress[0]
        )
        rows.append({
            "class": cls,
            "candidate_id": row["candidate_id"],
            "PT_surface_order_margin_eV": margin,
            "finite_positive_barriers": finite,
            "K_shield_for_one_retained_line_Pa_sqrt_m": Ksh,
            "profile_rho_min_m2": float(np.min(local_rho)),
            "profile_stress_first_Pa": float(local_stress[0]),
            "profile_stress_last_Pa": float(local_stress[-1]),
            "passed": passed,
        })
    result = {"model": "mpz_v9_11", "rows": rows, "passed": all(r["passed"] for r in rows)}
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if result["passed"] else 1)


if __name__ == "__main__":
    main()
