#!/usr/bin/env python3
"""Compare final reduced rows with existing and bounded direct PF evidence."""
from __future__ import annotations

import hashlib
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_outputs" / "oneD_v2_terminal_predictive_program"
PF_ROOT = Path(
    "/Volumes/Data/Data/Nanopillar_calculation/"
    "PF-fracture-fatigue_v10_2_21_persistent_sites_top1"
)
PF_RESPONSE = PF_ROOT / (
    "runs/1_final_v10_2_30_theta0_rate_extremes_four_class_1000um_base3621_v1/"
    "rate1x/temperature_response/v10_2_28_four_class_KJ_temperature_response.csv"
)
NEW_ROOT = Path("/private/tmp/oneD-v2-terminal-pf-transfer-runs")
TEMPERATURES = (300.0, 1000.0, 1200.0)
sys.path.insert(0, str(ROOT))

from scripts.run_oneD_v2_predictive_campaign import inputs, run_case, summary


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def step_metrics(path: Path, target_um: float = 100.0) -> dict:
    steps = pd.read_csv(path)
    within = steps.crack_extension_m.to_numpy(float) <= target_um * 1.0e-6 + 20.0e-6
    steps = steps.loc[within].copy()
    positions = np.flatnonzero(steps.n_fire.to_numpy(float) > 0.0)
    if not len(positions):
        raise RuntimeError(f"PF transfer trajectory contains no event: {path}")
    avalanche = 0
    counts = {0: 1}
    for index in range(1, len(positions)):
        between = steps.iloc[positions[index - 1] + 1 : positions[index] + 1]
        if (between.adaptive_frac.to_numpy(float) >= 1.0 - 1.0e-12).any():
            avalanche += 1
        counts[avalanche] = counts.get(avalanche, 0) + 1
    event_rows = steps.iloc[positions]
    initial = float(event_rows.iloc[0].KJ_Pa_sqrtm) * 1.0e-6
    terminal = float(steps.crack_extension_m.max()) * 1.0e6
    return {
        "pf_2D_initial_native_KJ_MPa_sqrt_m": initial,
        "pf_2D_event_count": int(len(event_rows)),
        "pf_2D_physical_avalanche_count": int(len(counts)),
        "pf_2D_precursor_reinitiation_count": int(len(counts) - 1),
        "pf_2D_largest_avalanche_fraction": max(counts.values()) / len(event_rows),
        "pf_2D_mean_event_size_um": float(event_rows.da_block_m.mean()) * 1.0e6,
        "pf_2D_max_backstress_GPa": float(steps.sigma_back_Pa.max()) * 1.0e-9,
        "pf_2D_max_mobile_count": float(steps.mpz_mobile_count.max()),
        "pf_2D_max_retained_count": float(steps.mpz_retained_count.max()),
        "pf_2D_terminal_extension_um": terminal,
        "pf_2D_target_right_censored": terminal >= target_um - 1.0e-6,
        "source_steps_file": str(path),
        "source_steps_sha256": sha(path),
    }


def main() -> int:
    registry = pd.read_csv(OUT / "oneD_v2_pf_transfer_registry.csv")
    physics, _, providers = inputs()
    mechanics, drive = providers["PF"]
    response = pd.read_csv(PF_RESPONSE)
    labels = {"Peak": "Peak", "DBTT": "DBTT"}
    records = []
    for material, row in zip(("Peak", "DBTT", "weak-T", "ceramic-like"), registry.itertuples(index=False)):
        row = pd.Series(row._asdict())
        for temperature in TEMPERATURES:
            reduced = summary(run_case(
                row, material, temperature, "PF", mechanics, drive, physics,
                100.0,
            ))
            if material in labels:
                source = response[
                    (response.plot_label == labels[material])
                    & (response.temperature_K == temperature)
                ].iloc[0]
                path = Path(source.steps_file)
                role = "EXISTING_AUTHORITATIVE_2D_PF"
                new_run = False
                seed = int(source.seed)
            else:
                seed = 2008666 if material == "weak-T" else 3008666
                tag = int(temperature)
                stem = f"steps_{tag:04d}K.csv" if tag < 1000 else f"steps_{tag}K.csv"
                path = NEW_ROOT / material / f"T{tag}K_seed{seed}" / stem
                role = "NEW_BOUNDED_2D_PF_FINALIST_TRANSFER"
                new_run = True
            direct = step_metrics(path)
            initial_error = (
                float(reduced["first_event_native_KJ_MPa_sqrt_m"])
                - direct["pf_2D_initial_native_KJ_MPa_sqrt_m"]
            ) / direct["pf_2D_initial_native_KJ_MPa_sqrt_m"]
            records.append({
                "material_class": material,
                "candidate_id": row.candidate_id,
                "temperature_K": temperature,
                "transfer_data_role": role,
                "pf_2D_seed": seed,
                **direct,
                "oneD_PF_status": reduced["status"],
                "oneD_PF_initial_native_KJ_MPa_sqrt_m": reduced["first_event_native_KJ_MPa_sqrt_m"],
                "oneD_PF_onset_envelope_max_MPa_sqrt_m": reduced["onset_envelope_max_MPa_sqrt_m"],
                "oneD_PF_event_count": reduced["event_count"],
                "oneD_PF_physical_avalanche_count": reduced["physical_avalanche_count"],
                "oneD_PF_precursor_reinitiation_count": reduced["precursor_reinitiation_count"],
                "oneD_PF_largest_avalanche_fraction": reduced["largest_avalanche_fraction"],
                "oneD_PF_mean_event_size_um": reduced["mean_event_size_um"],
                "oneD_PF_max_tip_radius_um": reduced["max_tip_radius_um"],
                "oneD_PF_max_backstress_GPa": reduced["max_backstress_GPa"],
                "initial_onset_relative_error": initial_error,
                "topology_match": int(reduced["physical_avalanche_count"]) == int(direct["pf_2D_physical_avalanche_count"]),
                "new_PF_run_launched": new_run,
                "new_FEMCZM_run_launched": False,
            })
    result = pd.DataFrame(records)
    result.to_csv(OUT / "oneD_v2_pf_transfer_results.csv", index=False)
    print(
        "PF_TRANSFER_COMPLETE "
        f"cases={len(result)} new={int(result.new_PF_run_launched.sum())} "
        f"max_abs_onset_error={result.initial_onset_relative_error.abs().max():.6f} "
        f"topology_matches={int(result.topology_match.sum())}/{len(result)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
