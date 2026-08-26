#!/usr/bin/env python3
"""V2 transaction/physical-avalanche analysis of bounded Peak/DBTT PF transfer."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_outputs" / "oneD_v2_peak_dbtt_rcurve_search"
PF_ROOT = Path(
    "/Volumes/Data/Data/Nanopillar_calculation/"
    "PF-fracture-fatigue_v10_2_21_persistent_sites_top1"
)
PF_RESPONSE = PF_ROOT / (
    "runs/1_final_v10_2_30_theta0_rate_extremes_four_class_1000um_base3621_v1/"
    "rate1x/temperature_response/v10_2_28_four_class_KJ_temperature_response.csv"
)
NEW_ROOT = Path("/private/tmp/oneD-v2-peak-dbtt-R-pf-runs")
TARGET_UM = 100.0

FINALISTS = {
    "Peak": ("oneD_v2_peak_R_41f8789bcbc1f097", 8666, (600.0, 900.0, 1200.0)),
    "DBTT": ("oneD_v2_dbtt_R_9f5160f509e713e2", 1008666, (600.0, 1100.0, 1200.0)),
}
CONTROLS = {
    "Peak": "v913_zeroD_sobol_0242980",
    "DBTT": "v913_zeroD_sobol_0202500",
}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def cases() -> list[dict]:
    response = pd.read_csv(PF_RESPONSE)
    records = []
    for material, (finalist, seed, temperatures) in FINALISTS.items():
        for temperature in temperatures:
            tag = int(temperature)
            stem = f"steps_{tag:04d}K.csv" if tag < 1000 else f"steps_{tag}K.csv"
            records.append({
                "material_class": material,
                "candidate_role": "R_FINALIST",
                "candidate_id": finalist,
                "temperature_K": temperature,
                "hazard_seed": seed,
                "seed_comparison_status": "CANONICAL_FINALIST_SEED",
                "new_PF_run": True,
                "steps_file": NEW_ROOT / f"{material}_R" / f"T{tag}K_seed{seed}" / stem,
            })
            old = response[
                (response.plot_label == material)
                & (response.temperature_K == temperature)
                & (response.candidate_id == CONTROLS[material])
            ]
            if len(old) != 1:
                raise ValueError(f"missing unique archived {material} control at {temperature:g} K")
            old = old.iloc[0]
            records.append({
                "material_class": material,
                "candidate_role": "CONTROL",
                "candidate_id": CONTROLS[material],
                "temperature_K": temperature,
                "hazard_seed": int(old.seed),
                "seed_comparison_status": "AUTHORITATIVE_TEMPERATURE_SWEEP_SEED_NOT_PAIRED",
                "new_PF_run": False,
                "steps_file": Path(old.steps_file),
            })
    return records


def _target_prefix(steps: pd.DataFrame) -> pd.DataFrame:
    reached = np.flatnonzero(steps.crack_extension_m.to_numpy(float) * 1.0e6 >= TARGET_UM - 1.0e-9)
    if not len(reached):
        raise RuntimeError("PF trajectory did not reach the 100 um target")
    return steps.iloc[: int(reached[0]) + 1].copy()


def classify(case: dict) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    path = Path(case["steps_file"])
    if not path.is_file():
        raise FileNotFoundError(path)
    all_steps = pd.read_csv(path)
    audit_path = path.parent / "anisotropic_emission_audit_v10174.json"
    audit_payload = json.loads(audit_path.read_text())
    audit = pd.DataFrame(audit_payload["records"])
    if len(audit) != len(all_steps):
        raise RuntimeError(f"PF step/audit record count mismatch: {path}")
    steps = _target_prefix(all_steps)
    audit = audit.iloc[:len(steps)].copy()
    steps["physical_time_s"] = steps.dt_cur_s.to_numpy(float).cumsum()
    event_positions = np.flatnonzero(steps.n_fire.to_numpy(float) > 0.0)
    if not len(event_positions):
        raise RuntimeError(f"PF trajectory has no cleavage transactions: {path}")
    rows = []
    avalanche = 0
    def audit_number(record: pd.Series, name: str) -> float:
        value = record.get(name, np.nan)
        return float(value) if value is not None else np.nan

    for transaction, pos in enumerate(event_positions):
        event = steps.iloc[pos]
        event_audit = audit.iloc[pos]
        before = steps.iloc[max(pos - 1, 0)]
        if transaction:
            prior = event_positions[transaction - 1]
            reload_slice = steps.iloc[prior + 1 : pos + 1]
            if (reload_slice.adaptive_frac.to_numpy(float) >= 1.0 - 1.0e-12).any():
                avalanche += 1
        next_pos = event_positions[transaction + 1] if transaction + 1 < len(event_positions) else None
        reload_slice = steps.iloc[pos + 1 : next_pos + 1] if next_pos is not None else steps.iloc[0:0]
        certified_reload = bool(
            next_pos is not None
            and (reload_slice.adaptive_frac.to_numpy(float) >= 1.0 - 1.0e-12).any()
        )
        rows.append({
            **{key: case[key] for key in (
                "material_class", "candidate_role", "candidate_id", "temperature_K",
                "hazard_seed", "seed_comparison_status", "new_PF_run",
            )},
            "event_transaction_index": transaction,
            "physical_avalanche_index": avalanche,
            "pre_event_step": int(event.step),
            "pre_event_time_s": float(event.physical_time_s),
            "pre_event_opening_m": float(event.Uapp_m),
            "pre_event_reaction_N_per_m": float(event.Ftop_N),
            "pre_event_native_J_J_per_m2": float(event.J_effective_direct_J_per_m2),
            "pre_event_native_KJ_MPa_sqrt_m": float(event.KJ_Pa_sqrtm) * 1.0e-6,
            "pre_event_projected_extension_um": float(before.crack_extension_m) * 1.0e6,
            "post_event_projected_extension_um": float(event.crack_extension_m) * 1.0e6,
            "event_extension_um": float(event.crack_extension_m - before.crack_extension_m) * 1.0e6,
            "event_driver_da_block_um": float(event.da_block_m) * 1.0e6,
            "pre_event_tip_radius_um": audit_number(event_audit, "persistent_tip_radius_m") * 1.0e6,
            "pre_event_front_width_um": audit_number(event_audit, "persistent_site_front_width_m") * 1.0e6,
            "pre_event_backstress_GPa": float(event.sigma_back_Pa) * 1.0e-9,
            "pre_event_mobile_count": float(event.mpz_mobile_count),
            "pre_event_retained_count": float(event.mpz_retained_count),
            "pre_event_source_multiplicity": audit_number(event_audit, "persistent_site_multiplicity_per_system"),
            "pre_event_signed_shielding_MPa_sqrt_m": float(event.mpz_K_shield_Pa_sqrt_m) * 1.0e-6,
            "reload_time_to_next_event_s": np.nan if next_pos is None else float(
                steps.iloc[next_pos].physical_time_s - event.physical_time_s
            ),
            "reload_opening_to_next_event_m": np.nan if next_pos is None else float(
                steps.iloc[next_pos].Uapp_m - event.Uapp_m
            ),
            "certified_reload_before_next_event": certified_reload,
            "right_censored_at_target": transaction == len(event_positions) - 1,
            "in_avalanche_drive_interpretation": "PF_MODEL_NATIVE_DRIVING_TRAJECTORY_NOT_RESISTANCE",
            "source_steps_file": str(path.resolve()),
            "source_steps_sha256": sha(path),
            "source_tensor_audit_file": str(audit_path.resolve()),
            "source_tensor_audit_sha256": sha(audit_path),
            "source_process_zone_state_availability": (
                "FULL_PERSISTENT_SITE_OBSERVER"
                if "persistent_tip_radius_m" in audit.columns
                else "ARCHIVED_SCALAR_STATE_ONLY"
            ),
        })
    tx = pd.DataFrame(rows)
    avalanches = []
    for avalanche_id, group in tx.groupby("physical_avalanche_index", sort=True):
        avalanches.append({
            **{key: case[key] for key in (
                "material_class", "candidate_role", "candidate_id", "temperature_K",
                "hazard_seed", "seed_comparison_status", "new_PF_run",
            )},
            "physical_avalanche_index": int(avalanche_id),
            "first_event_transaction_index": int(group.event_transaction_index.min()),
            "last_event_transaction_index": int(group.event_transaction_index.max()),
            "event_transaction_count": int(len(group)),
            "start_extension_um": float(group.pre_event_projected_extension_um.iloc[0]),
            "end_extension_um": float(group.post_event_projected_extension_um.iloc[-1]),
            "avalanche_extension_um": float(group.event_extension_um.sum()),
            "onset_opening_m": float(group.pre_event_opening_m.iloc[0]),
            "onset_native_KJ_MPa_sqrt_m": float(group.pre_event_native_KJ_MPa_sqrt_m.iloc[0]),
            "target_right_censored": bool(group.right_censored_at_target.iloc[-1]),
            "grouping_rule": "CONTIGUOUS_EVENTS_WITHOUT_INTERVENING_FULL_ACCEPTED_LOADING_INTERVAL",
            "source_steps_file": str(path.resolve()),
            "source_steps_sha256": sha(path),
        })
    ava = pd.DataFrame(avalanches)
    onset = tx.groupby("physical_avalanche_index", sort=True).head(1).copy()
    onset["onset_role"] = np.where(
        onset.physical_avalanche_index == 0,
        "INITIAL_ONSET_PRE_EVENT", "REINITIATION_ONSET_PRE_EVENT",
    )
    onset["resistance_candidate_policy"] = "RELOAD_SEPARATED_PRE_EVENT_ONLY"
    return tx, ava, onset


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    grouped = [classify(case) for case in cases()]
    tx = pd.concat([item[0] for item in grouped], ignore_index=True)
    ava = pd.concat([item[1] for item in grouped], ignore_index=True)
    onset = pd.concat([item[2] for item in grouped], ignore_index=True)
    tx_path = OUT / "pf_2d_peak_dbtt_R_event_transactions.csv"
    ava_path = OUT / "pf_2d_peak_dbtt_R_physical_avalanches.csv"
    onset_path = OUT / "pf_2d_peak_dbtt_R_onset_candidates.csv"
    tx.to_csv(tx_path, index=False)
    ava.to_csv(ava_path, index=False)
    onset.to_csv(onset_path, index=False)

    summaries = []
    for keys, local in tx.groupby(
        ["material_class", "candidate_role", "candidate_id", "temperature_K",
         "hazard_seed", "seed_comparison_status", "new_PF_run"], sort=False,
    ):
        material, role, candidate, temperature, seed, seed_status, new_run = keys
        starts = local.groupby("physical_avalanche_index", sort=True).head(1)
        avalanche_extension = local.groupby("physical_avalanche_index").event_extension_um.sum()
        first = float(starts.pre_event_native_KJ_MPa_sqrt_m.iloc[0])
        maximum = float(starts.pre_event_native_KJ_MPa_sqrt_m.max())
        reinit_values = starts.pre_event_native_KJ_MPa_sqrt_m.iloc[1:].to_numpy(float)
        increments = np.diff(starts.pre_event_native_KJ_MPa_sqrt_m.to_numpy(float))
        if len(increments) == 0:
            monotonicity = "NOT_APPLICABLE"
        elif np.all(increments >= -1.0e-12):
            monotonicity = "NONDECREASING"
        elif np.all(increments <= 1.0e-12):
            monotonicity = "NONINCREASING"
        else:
            monotonicity = "MIXED"
        summaries.append({
            "material_class": material,
            "candidate_role": role,
            "candidate_id": candidate,
            "temperature_K": temperature,
            "hazard_seed": seed,
            "seed_comparison_status": seed_status,
            "new_PF_run": new_run,
            "target_um": TARGET_UM,
            "terminal_semantics": "TARGET_RIGHT_CENSORED",
            "event_transaction_count": int(len(local)),
            "physical_avalanche_count": int(local.physical_avalanche_index.nunique()),
            "N_reinit": int(local.physical_avalanche_index.nunique() - 1),
            "initial_onset_native_KJ_MPa_sqrt_m": first,
            "maximum_onset_native_KJ_MPa_sqrt_m": maximum,
            "deltaK_reinit_MPa_sqrt_m": maximum - first,
            "relative_deltaK_reinit": (maximum - first) / first,
            "first_reinitiation_onset_native_KJ_MPa_sqrt_m": (
                np.nan if not len(reinit_values) else float(reinit_values[0])
            ),
            "maximum_reinitiation_onset_native_KJ_MPa_sqrt_m": (
                np.nan if not len(reinit_values) else float(np.max(reinit_values))
            ),
            "signed_max_reinitiation_minus_initial_K_MPa_sqrt_m": (
                np.nan if not len(reinit_values) else float(np.max(reinit_values) - first)
            ),
            "onset_sequence_monotonicity": monotonicity,
            "largest_avalanche_fraction": float(avalanche_extension.max() / avalanche_extension.sum()),
            "maximum_tip_radius_um": float(local.pre_event_tip_radius_um.max()),
            "maximum_backstress_GPa": float(local.pre_event_backstress_GPa.max()),
            "maximum_mobile_count": float(local.pre_event_mobile_count.max()),
            "maximum_retained_count": float(local.pre_event_retained_count.max()),
            "target_right_censored": True,
            "source_steps_file": local.source_steps_file.iloc[0],
            "source_steps_sha256": local.source_steps_sha256.iloc[0],
        })
    summary = pd.DataFrame(summaries)
    summary_path = OUT / "pf_2d_peak_dbtt_R_transfer_summary.csv"
    summary.to_csv(summary_path, index=False)
    manifest = {
        "schema": "pf_2d_peak_dbtt_R_transfer_analysis_v1",
        "case_count": len(summary),
        "new_PF_run_count": int(summary.new_PF_run.sum()),
        "new_FEMCZM_run_count": 0,
        "maximum_concurrent_heavy_PF_workers": 2,
        "target_um": TARGET_UM,
        "grouping_rule": "CONTIGUOUS_EVENTS_WITHOUT_INTERVENING_FULL_ACCEPTED_LOADING_INTERVAL",
        "resistance_candidate_policy": "RELOAD_SEPARATED_PRE_EVENT_ONLY",
        "control_pairing_qualification": "TEMPERATURE_MATCHED_AUTHORITATIVE_UNPAIRED_SEEDS",
        "artifacts": {
            path.name: sha(path)
            for path in (tx_path, ava_path, onset_path, summary_path)
        },
    }
    (OUT / "pf_2d_peak_dbtt_R_transfer_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    print(summary[["material_class", "candidate_role", "temperature_K", "hazard_seed",
                   "event_transaction_count", "physical_avalanche_count", "N_reinit",
                   "initial_onset_native_KJ_MPa_sqrt_m",
                   "maximum_onset_native_KJ_MPa_sqrt_m"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
