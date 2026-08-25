#!/usr/bin/env python3
"""Reduced loading-rate trends with authoritative PF rate-extreme checks."""
from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_outputs" / "oneD_v2_terminal_predictive_program"
PF_RUN = Path(
    "/Volumes/Data/Data/Nanopillar_calculation/"
    "PF-fracture-fatigue_v10_2_21_persistent_sites_top1/"
    "runs/1_final_v10_2_30_theta0_rate_extremes_four_class_1000um_base3621_v1"
)
sys.path.insert(0, str(ROOT))

from scripts.run_oneD_v2_predictive_campaign import inputs, run_case, summary

RATE_DIRECTORIES = {0.01: "rate0p01x", 1.0: "rate1x", 100.0: "rate100x"}
TEMPERATURES = (600.0, 1000.0, 1200.0)


def main() -> int:
    registry = pd.read_csv(OUT / "oneD_v2_pf_transfer_registry.csv")
    rows = {"Peak": registry.iloc[0], "DBTT": registry.iloc[1]}
    physics, _, providers = inputs()
    direct = {}
    for factor, directory in RATE_DIRECTORIES.items():
        path = PF_RUN / directory / "temperature_response/v10_2_28_four_class_KJ_temperature_response.csv"
        frame = pd.read_csv(path)
        for material in rows:
            q = frame[(frame.plot_label == material) & (frame.temperature_K == 1000.0)].iloc[0]
            direct[(material, factor)] = float(q.initial_K_MPa_sqrt_m)
    records = []
    for material, row in rows.items():
        for backend, (mechanics, drive) in providers.items():
            for temperature in TEMPERATURES:
                for factor in RATE_DIRECTORIES:
                    result = summary(run_case(
                        row, material, temperature, backend, mechanics, drive,
                        physics, 100.0, nominal_dt_s=8.4 / factor,
                    ))
                    result.update({
                        "loading_rate_factor": factor,
                        "nominal_opening_rate_m_s": 0.2e-6 / (8.4 / factor),
                        "authoritative_PF_2D_initial_KJ_MPa_sqrt_m": (
                            direct[(material, factor)]
                            if backend == "PF" and temperature == 1000.0 else None
                        ),
                    })
                    if result["authoritative_PF_2D_initial_KJ_MPa_sqrt_m"] is not None:
                        result["PF_2D_initial_relative_error"] = (
                            result["first_event_native_KJ_MPa_sqrt_m"]
                            - result["authoritative_PF_2D_initial_KJ_MPa_sqrt_m"]
                        ) / result["authoritative_PF_2D_initial_KJ_MPa_sqrt_m"]
                    records.append(result)
    frame = pd.DataFrame(records)
    frame.to_csv(OUT / "oneD_v2_loading_rate_sensitivity.csv", index=False)
    print(
        "RATE_SENSITIVITY_COMPLETE "
        f"cases={len(frame)} complete={int(frame.status.eq('TARGET_RIGHT_CENSORED').sum())}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
