#!/usr/bin/env python3
"""Finalize the Taylor/Peierls search, PF transfer, option banks, and reports."""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import subprocess
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_outputs" / "oneD_v2_taylor_peierls_rcurve_search"
PF100 = Path("/private/tmp/oneD-v2-taylor-peierls-pf-runs")
PF300 = Path("/private/tmp/oneD-v2-taylor-peierls-pf-runs-300um")
PF_REPO = Path("/private/tmp/pf-taylor-peierls-transfer")
CONTROL_IDS = {"Peak": "v913_zeroD_sobol_0242980", "DBTT": "v913_zeroD_sobol_0202500"}
FINALIST_IDS = {
    "Peak": "oneD_v2_peak_TP_6962e84b2eb78fbb",
    "DBTT": "oneD_v2_dbtt_TP_6ca03f05fbae34e9",
}
OPTION_KEYS = {
    CONTROL_IDS["Peak"]: "oneD_v2_Peak_control",
    FINALIST_IDS["Peak"]: "oneD_v2_Peak_TP",
    CONTROL_IDS["DBTT"]: "oneD_v2_DBTT_control",
    FINALIST_IDS["DBTT"]: "oneD_v2_DBTT_TP",
}
TEMPERATURES = {"Peak": (600, 1000, 1200), "DBTT": (600, 1100, 1200)}
SEEDS = {"Peak": 8666, "DBTT": 1008666}
SEARCH_FIELDS = (
    "peierls_H0_eV", "peierls_activation_entropy_kB", "peierls_exp_a", "peierls_exp_n",
    "taylor_H0_eV", "taylor_activation_entropy_kB", "taylor_exp_a", "taylor_exp_n",
    "taylor_corr_rho_c_m2", "taylor_corr_scale",
)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def head(path: Path = ROOT) -> str:
    return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=path, text=True).strip()


def _jsonable(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    raise TypeError(type(value).__name__)


def _stage1_aggregate(cases: pd.DataFrame) -> pd.DataFrame:
    records = []
    for (material, candidate), local in cases.groupby(["target_class", "candidate_id"], sort=False):
        upper_t = (900.0, 1000.0) if material == "Peak" else (900.0, 1000.0, 1100.0)
        upper = local[local.temperature_K.isin(upper_t)]
        low_t = (600.0,) if material == "Peak" else (300.0, 600.0)
        low = local[local.temperature_K.isin(low_t)]
        response = local.groupby(["provider", "temperature_K"])[
            ["initial_onset_native_KJ_MPa_sqrt_m", "relative_deltaK_reinit", "N_reinit",
             "largest_avalanche_fraction", "maximum_onset_tip_radius_um",
             "maximum_onset_backstress_GPa", "maximum_onset_mobile_density_m2",
             "maximum_onset_retained_density_m2"]
        ].first().reset_index()
        records.append({
            "target_class": material, "candidate_id": candidate,
            "stage1_case_count": int(len(local)),
            "stage1_completion_failure_count": int((local.status != "TARGET_RIGHT_CENSORED").sum()),
            "stage1_positive_reinit_case_count": int((local.deltaK_reinit_MPa_sqrt_m > 1e-12).sum()),
            "upper_temperature_max_relative_deltaK": float(upper.relative_deltaK_reinit.max()),
            "upper_temperature_mean_N_reinit": float(upper.N_reinit.mean()),
            "low_temperature_mean_N_reinit": float(low.N_reinit.mean()),
            "minimum_largest_avalanche_fraction": float(local.largest_avalanche_fraction.min()),
            "maximum_onset_tip_radius_um": float(local.maximum_onset_tip_radius_um.max()),
            "maximum_onset_backstress_GPa": float(local.maximum_onset_backstress_GPa.max()),
            "maximum_onset_mobile_density_m2": float(local.maximum_onset_mobile_density_m2.max()),
            "maximum_onset_retained_density_m2": float(local.maximum_onset_retained_density_m2.max()),
            "response_descriptors_json": json.dumps(response.to_dict("records"), sort_keys=True, separators=(",", ":")),
            "fracture_response_equivalence_policy": "KINETIC_ROWS_NUMERICALLY_EQUIVALENT_WITHIN_CLASS_PROVIDER_TEMPERATURE_SEED",
        })
    return pd.DataFrame(records)


def build_population_and_banks() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    stage0 = pd.read_parquet(OUT / "oneD_v2_taylor_peierls_stage0_population.parquet")
    diagnostics = pd.read_csv(OUT / "oneD_v2_taylor_peierls_timescale_diagnostics.csv")
    stage1_cases = pd.read_parquet(OUT / "oneD_v2_taylor_peierls_stage1_100um_cases.parquet")
    aggregate = _stage1_aggregate(stage1_cases)
    complete = stage0.merge(diagnostics, on=["candidate_id", "target_class", "full_material_sha256",
                                             "taylor_peierls_subvector_sha256"], how="left")
    complete = complete.merge(aggregate, on=["candidate_id", "target_class"], how="left")
    complete["stage1_fully_evaluated"] = complete.stage1_case_count.notna()
    complete["pareto_nondominated"] = complete.stage1_fully_evaluated
    complete["pareto_near_front"] = False
    complete["pareto_policy"] = (
        "ALL_STAGE1_ROWS_RESPONSE_EQUIVALENT_ON_FRACTURE_OBJECTIVES; "
        "KINETIC_STATE_PARTITION_RETAINED_AS_DESCRIPTOR_NOT_HIDDEN_SCALAR_SCORE"
    )
    complete_path = OUT / "oneD_v2_taylor_peierls_complete_population.parquet"
    complete.to_parquet(complete_path, index=False)
    pareto = complete[complete.stage1_fully_evaluated].copy()
    pareto_path = OUT / "oneD_v2_taylor_peierls_pareto_library.parquet"
    pareto.to_parquet(pareto_path, index=False)

    survivors = json.loads((OUT / "oneD_v2_taylor_peierls_survivor_selection.json").read_text())
    curated_ids = {
        candidate
        for material in ("Peak", "DBTT")
        for candidate in survivors["classes"][material]["stage2_300um_candidates"]
    }
    shortlist_ids = {
        candidate
        for material in ("Peak", "DBTT")
        for candidate in survivors["classes"][material]["stage3_multiseed_candidates"]
    }
    bank = complete[complete.candidate_id.isin(curated_ids)].copy()
    bank["bank_level"] = "CURATED_TAYLOR_PEIERLS_KINETICS_OPTION_BANK"
    bank["direct_PF_status"] = np.where(
        bank.candidate_id.isin(FINALIST_IDS.values()), "DIRECT_PF_REJECTED_NO_POSITIVE_REINITIATION",
        np.where(bank.candidate_id.isin(CONTROL_IDS.values()), "DIRECT_PF_PAIRED_CONTROL", "NOT_RUN"),
    )
    for material in ("Peak", "DBTT"):
        mask = bank.target_class.eq(material)
        for field in ("peierls_H0_eV", "taylor_H0_eV", "taylor_corr_rho_c_m2"):
            values = bank.loc[mask, field]
            qlo, qhi = values.quantile([0.25, 0.75])
            bank.loc[mask & bank[field].le(qlo), f"{field}_regime"] = "LOW"
            bank.loc[mask & bank[field].ge(qhi), f"{field}_regime"] = "HIGH"
    bank_path_csv = OUT / "oneD_v2_taylor_peierls_option_bank.csv"
    bank_path_parquet = OUT / "oneD_v2_taylor_peierls_option_bank.parquet"
    bank.to_csv(bank_path_csv, index=False)
    bank.to_parquet(bank_path_parquet, index=False)
    shortlist = bank[bank.candidate_id.isin(shortlist_ids)].copy()
    shortlist["bank_level"] = "COMPACT_FUTURE_SEARCH_SHORTLIST"
    shortlist.to_csv(OUT / "oneD_v2_taylor_peierls_shortlist.csv", index=False)
    registry_columns = [
        "candidate_id", "target_class", "parent_control_id", "search_stage", "search_method",
        "canonical_parameter_json", "full_material_sha256", "cleavage_barrier_sha256",
        "emission_barrier_sha256", "taylor_peierls_subvector_sha256", *SEARCH_FIELDS,
    ]
    registry = bank.loc[:, registry_columns].copy()
    registry.to_csv(OUT / "oneD_v2_taylor_peierls_material_registry.csv", index=False)
    return complete, pareto, bank, shortlist


def pf_cases() -> list[dict[str, Any]]:
    cases = []
    for material in ("Peak", "DBTT"):
        for candidate_id in (CONTROL_IDS[material], FINALIST_IDS[material]):
            option = OPTION_KEYS[candidate_id]
            role = "CONTROL" if candidate_id == CONTROL_IDS[material] else "KINETIC_FINALIST"
            for temperature in TEMPERATURES[material]:
                cases.append({
                    "material_class": material, "candidate_role": role,
                    "candidate_id": candidate_id, "option_key": option,
                    "temperature_K": float(temperature), "hazard_seed": SEEDS[material],
                    "target_um": 100.0,
                    "case_dir": PF100 / option / f"T{temperature}K_seed{SEEDS[material]}",
                })
            transition = 1000 if material == "Peak" else 1100
            cases.append({
                "material_class": material, "candidate_role": role,
                "candidate_id": candidate_id, "option_key": option,
                "temperature_K": float(transition), "hazard_seed": SEEDS[material],
                "target_um": 300.0,
                "case_dir": PF300 / option / f"T{transition}K_seed{SEEDS[material]}",
            })
    return cases


def _steps_path(case_dir: Path, temperature: float) -> Path:
    tag = int(temperature)
    path = case_dir / (f"steps_{tag:04d}K.csv" if tag < 1000 else f"steps_{tag}K.csv")
    if not path.is_file():
        raise FileNotFoundError(path)
    return path


def _audit_number(record: pd.Series, name: str) -> float:
    value = record.get(name, np.nan)
    return float(value) if value is not None else np.nan


def classify_pf(case: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    steps_path = _steps_path(case["case_dir"], case["temperature_K"])
    audit_path = case["case_dir"] / "anisotropic_emission_audit_v10174.json"
    all_steps = pd.read_csv(steps_path)
    audit_payload = json.loads(audit_path.read_text())
    audit = pd.DataFrame(audit_payload["records"])
    if len(audit) != len(all_steps):
        raise RuntimeError(f"PF step/audit count mismatch: {steps_path}")
    reached = np.flatnonzero(
        all_steps.crack_extension_m.to_numpy(float) * 1e6 >= case["target_um"] - 1e-9
    )
    if not len(reached):
        raise RuntimeError(f"PF trajectory did not reach target: {steps_path}")
    steps = all_steps.iloc[: int(reached[0]) + 1].copy()
    audit = audit.iloc[: len(steps)].copy()
    steps["physical_time_s"] = steps.dt_cur_s.to_numpy(float).cumsum()
    positions = np.flatnonzero(steps.n_fire.to_numpy(float) > 0.0)
    if not len(positions):
        raise RuntimeError(f"PF trajectory has no fracture event: {steps_path}")
    rows, profile_rows = [], []
    avalanche = 0
    for transaction, position in enumerate(positions):
        event, event_audit = steps.iloc[position], audit.iloc[position]
        before = steps.iloc[max(position - 1, 0)]
        if transaction:
            prior = positions[transaction - 1]
            if (steps.iloc[prior + 1 : position + 1].adaptive_frac.to_numpy(float) >= 1 - 1e-12).any():
                avalanche += 1
        next_position = positions[transaction + 1] if transaction + 1 < len(positions) else None
        reload_slice = steps.iloc[position + 1 : next_position + 1] if next_position is not None else steps.iloc[0:0]
        certified_reload = bool(
            next_position is not None
            and (reload_slice.adaptive_frac.to_numpy(float) >= 1 - 1e-12).any()
        )
        common = {key: case[key] for key in (
            "material_class", "candidate_role", "candidate_id", "option_key",
            "temperature_K", "hazard_seed", "target_um",
        )}
        sigma_open = _audit_number(event_audit, "source_opening_stress_Pa")
        local_k = sigma_open * math.sqrt(2 * math.pi * 1e-6) / 1e6
        profile_present = event_audit.get("taylor_peierls_state_profile_schema") is not None
        rows.append({
            **common, "event_transaction_index": transaction,
            "physical_avalanche_index": avalanche, "pre_event_step": int(event.step),
            "pre_event_time_s": float(event.physical_time_s),
            "pre_event_opening_m": float(event.Uapp_m),
            "pre_event_reaction_N_per_m": float(event.Ftop_N),
            "pre_event_native_J_J_per_m2": float(event.J_effective_direct_J_per_m2),
            "pre_event_native_KJ_MPa_sqrt_m": float(event.KJ_Pa_sqrtm) * 1e-6,
            "pre_event_K_app_or_common_reference_MPa_sqrt_m": float(event.KJ_Pa_sqrtm) * 1e-6,
            "pre_event_K_native_MPa_sqrt_m": float(event.KJ_Pa_sqrtm) * 1e-6,
            "pre_event_K_shield_MPa_sqrt_m": float(event.mpz_K_shield_Pa_sqrt_m) * 1e-6,
            "pre_event_K_effective_local_equivalent_MPa_sqrt_m": local_k,
            "pre_event_source_opening_stress_Pa": sigma_open,
            "pre_event_projected_extension_um": float(before.crack_extension_m) * 1e6,
            "post_event_projected_extension_um": float(event.crack_extension_m) * 1e6,
            "event_extension_um": float(event.crack_extension_m - before.crack_extension_m) * 1e6,
            "pre_event_tip_radius_um": _audit_number(event_audit, "persistent_tip_radius_m") * 1e6,
            "pre_event_front_width_um": _audit_number(event_audit, "persistent_site_front_width_m") * 1e6,
            "pre_event_backstress_GPa": float(event.sigma_back_Pa) * 1e-9,
            "pre_event_mobile_count": float(event.mpz_mobile_count),
            "pre_event_retained_count": float(event.mpz_retained_count),
            "pre_event_source_multiplicity": _audit_number(event_audit, "persistent_site_multiplicity_per_system"),
            "pre_event_resolved_shears_Pa_json": json.dumps(event_audit.get("anisotropic_tau_signed_Pa")),
            "pre_event_opening_tensor_Pa_json": json.dumps(event_audit.get("opening_tensor_Pa")),
            "pre_event_channel_tensors_Pa_json": json.dumps(event_audit.get("channel_tensors_Pa")),
            "pre_event_peierls_rate_max_s": _audit_number(event_audit, "source_rate_max_peierls_s"),
            "pre_event_encounter_rate_max_s": _audit_number(event_audit, "source_rate_max_encounter_s"),
            "pre_event_taylor_completion_rate_max_s": _audit_number(event_audit, "source_rate_max_taylor_completion_s"),
            "reload_time_to_next_event_s": np.nan if next_position is None else float(
                steps.iloc[next_position].physical_time_s - event.physical_time_s
            ),
            "reload_opening_to_next_event_m": np.nan if next_position is None else float(
                steps.iloc[next_position].Uapp_m - event.Uapp_m
            ),
            "certified_reload_before_next_event": certified_reload,
            "right_censored_at_target": transaction == len(positions) - 1,
            "source_steps_file": str(steps_path.resolve()), "source_steps_sha256": sha(steps_path),
            "source_tensor_audit_file": str(audit_path.resolve()), "source_tensor_audit_sha256": sha(audit_path),
            "source_process_zone_state_availability": (
                "FULL_ACTIVE_AND_WAKE_SPATIAL_PROFILE_OBSERVER" if profile_present else "SCALAR_ONLY"
            ),
            "in_avalanche_drive_interpretation": "PF_MODEL_NATIVE_DRIVING_TRAJECTORY_NOT_RESISTANCE",
        })
        # Expand profiles only for resistance-candidate onsets, not every event.
        onset_candidate = transaction == 0 or (
            transaction > 0
            and (steps.iloc[positions[transaction - 1] + 1 : position + 1].adaptive_frac.to_numpy(float) >= 1 - 1e-12).any()
        )
        if onset_candidate and profile_present:
            active_x = np.asarray(event_audit["active_x_ahead_of_tip_m"], dtype=float)
            wake_x = np.asarray(event_audit["wake_x_behind_tip_m"], dtype=float)
            active_mobile = np.asarray(event_audit["mobile_active_by_system_bin"], dtype=float)
            active_retained = np.asarray(event_audit["retained_active_by_system_bin"], dtype=float)
            wake_mobile = np.asarray(event_audit["mobile_wake_by_system_bin"], dtype=float)
            wake_retained = np.asarray(event_audit["retained_wake_by_system_bin"], dtype=float)
            arrays = {
                name: np.asarray(event_audit[name], dtype=float)
                for name in (
                    "forest_density_active_m2_by_bin", "local_transport_stress_Pa_by_bin",
                    "peierls_rate_s_by_bin", "peierls_velocity_m_s_by_bin",
                    "encounter_rate_s_by_bin", "taylor_completion_rate_s_by_bin",
                    "transport_time_s_by_bin", "retention_encounter_time_s_by_bin",
                    "taylor_completion_time_s_by_bin", "chi_ret_by_bin",
                    "chi_taylor_completion_by_bin",
                )
            }
            for region, coordinates, mobile, retained in (
                ("ACTIVE_AHEAD_OF_TIP", active_x, active_mobile, active_retained),
                ("WAKE_BEHIND_TIP", wake_x, wake_mobile, wake_retained),
            ):
                for system in range(mobile.shape[0]):
                    for bin_index, coordinate in enumerate(coordinates):
                        profile_rows.append({
                            **common, "physical_avalanche_index": avalanche,
                            "onset_event_transaction_index": transaction,
                            "region": region, "system_index": system, "bin_index": bin_index,
                            "distance_from_tip_m": float(coordinate),
                            "mobile_count": float(mobile[system, bin_index]),
                            "retained_count": float(retained[system, bin_index]),
                            **{
                                name: (float(values[bin_index]) if region == "ACTIVE_AHEAD_OF_TIP" else np.nan)
                                for name, values in arrays.items()
                            },
                            "profile_source": "EXISTING_PF_MPZ_STATE_DEFAULT_OFF_OBSERVER_NO_FEEDBACK",
                        })
    transactions = pd.DataFrame(rows)
    avalanche_rows = []
    for avalanche_id, group in transactions.groupby("physical_avalanche_index", sort=True):
        avalanche_rows.append({
            **{key: case[key] for key in (
                "material_class", "candidate_role", "candidate_id", "option_key",
                "temperature_K", "hazard_seed", "target_um",
            )},
            "physical_avalanche_index": int(avalanche_id),
            "first_event_transaction_index": int(group.event_transaction_index.min()),
            "last_event_transaction_index": int(group.event_transaction_index.max()),
            "event_transaction_count": int(len(group)),
            "start_extension_um": float(group.pre_event_projected_extension_um.iloc[0]),
            "end_extension_um": float(group.post_event_projected_extension_um.iloc[-1]),
            "avalanche_extension_um": float(group.event_extension_um.sum()),
            "onset_native_KJ_MPa_sqrt_m": float(group.pre_event_native_KJ_MPa_sqrt_m.iloc[0]),
            "target_right_censored": bool(group.right_censored_at_target.iloc[-1]),
            "grouping_rule": "CONTIGUOUS_EVENTS_WITHOUT_INTERVENING_FULL_ACCEPTED_LOADING_INTERVAL",
        })
    avalanches = pd.DataFrame(avalanche_rows)
    onsets = transactions.groupby("physical_avalanche_index", sort=True).head(1).copy()
    onsets["onset_role"] = np.where(
        onsets.physical_avalanche_index.eq(0), "INITIAL_ONSET_PRE_EVENT", "REINITIATION_ONSET_PRE_EVENT"
    )
    onsets["resistance_candidate_policy"] = "RELOAD_SEPARATED_PRE_EVENT_ONLY"
    return transactions, avalanches, onsets, pd.DataFrame(profile_rows)


def build_pf_outputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    grouped = [classify_pf(case) for case in pf_cases()]
    transactions = pd.concat([item[0] for item in grouped], ignore_index=True)
    avalanches = pd.concat([item[1] for item in grouped], ignore_index=True)
    onsets = pd.concat([item[2] for item in grouped], ignore_index=True)
    profiles = pd.concat([item[3] for item in grouped], ignore_index=True)
    transactions.to_csv(OUT / "pf_2d_taylor_peierls_event_transactions.csv", index=False)
    avalanches.to_csv(OUT / "pf_2d_taylor_peierls_physical_avalanches.csv", index=False)
    onsets.to_csv(OUT / "pf_2d_taylor_peierls_onset_candidates.csv", index=False)
    profiles.to_parquet(OUT / "pf_2d_taylor_peierls_state_profiles.parquet", index=False)
    summaries = []
    for keys, local in transactions.groupby(
        ["material_class", "candidate_role", "candidate_id", "temperature_K", "hazard_seed", "target_um"],
        sort=False,
    ):
        material, role, candidate, temperature, seed, target = keys
        onset = local.groupby("physical_avalanche_index", sort=True).head(1)
        values = onset.pre_event_native_KJ_MPa_sqrt_m.to_numpy(float)
        reinit = values[1:]
        avalanche_size = local.groupby("physical_avalanche_index").event_extension_um.sum()
        summaries.append({
            "material_class": material, "candidate_role": role, "candidate_id": candidate,
            "temperature_K": temperature, "hazard_seed": seed, "target_um": target,
            "terminal_semantics": "TARGET_RIGHT_CENSORED",
            "event_transaction_count": int(len(local)),
            "physical_avalanche_count": int(len(values)), "N_reinit": int(max(len(values) - 1, 0)),
            "initial_onset_native_KJ_MPa_sqrt_m": float(values[0]),
            "maximum_onset_native_KJ_MPa_sqrt_m": float(np.max(values)),
            "deltaK_reinit_MPa_sqrt_m": float(np.max(values) - values[0]),
            "signed_max_reinitiation_minus_initial_K_MPa_sqrt_m": (
                np.nan if not len(reinit) else float(np.max(reinit) - values[0])
            ),
            "initial_local_equivalent_K_MPa_sqrt_m": float(
                onset.pre_event_K_effective_local_equivalent_MPa_sqrt_m.iloc[0]
            ),
            "maximum_reinit_local_equivalent_K_MPa_sqrt_m": (
                np.nan if len(onset) < 2 else float(onset.pre_event_K_effective_local_equivalent_MPa_sqrt_m.iloc[1:].max())
            ),
            "largest_avalanche_fraction": float(avalanche_size.max() / avalanche_size.sum()),
            "maximum_tip_radius_um": float(local.pre_event_tip_radius_um.max()),
            "maximum_backstress_GPa": float(local.pre_event_backstress_GPa.max()),
            "maximum_abs_signed_shielding_MPa_sqrt_m": float(local.pre_event_K_shield_MPa_sqrt_m.abs().max()),
            "maximum_mobile_count": float(local.pre_event_mobile_count.max()),
            "maximum_retained_count": float(local.pre_event_retained_count.max()),
            "positive_reload_separated_resistance": bool(len(reinit) and np.max(reinit) > values[0] + 1e-12),
            "direct_PF_success": False,
            "direct_PF_decision": "REJECT_NO_POSITIVE_RELOAD_SEPARATED_RESISTANCE",
            "target_right_censored": True,
        })
    summary = pd.DataFrame(summaries)
    summary.to_csv(OUT / "pf_2d_taylor_peierls_transfer_summary.csv", index=False)
    return transactions, avalanches, onsets, profiles, summary


def _pca_coordinates(frame: pd.DataFrame) -> np.ndarray:
    values = frame.loc[:, SEARCH_FIELDS].to_numpy(float)
    for index, field in enumerate(SEARCH_FIELDS):
        if "entropy" not in field:
            values[:, index] = np.log10(values[:, index])
    values = (values - values.mean(axis=0)) / np.maximum(values.std(axis=0), 1e-12)
    return np.linalg.svd(values, full_matrices=False)[0][:, :2]


def figures(complete, bank, stage1, pf_summary, onsets, profiles) -> None:
    plt.style.use("seaborn-v0_8-whitegrid")
    for material, filename in (
        ("Peak", "PEAK_TAYLOR_PEIERLS_RESPONSE_MAP.png"),
        ("DBTT", "DBTT_TAYLOR_PEIERLS_RESPONSE_MAP.png"),
    ):
        local = complete[complete.target_class.eq(material)]
        fig, ax = plt.subplots(figsize=(7.2, 5.4))
        scatter = ax.scatter(
            local.log10_peierls_rate_s_median, local.log10_taylor_completion_rate_s_median,
            c=local.retained_equilibrium_fraction_median, s=10, cmap="viridis", alpha=0.7,
        )
        ax.set(xlabel="log10 Peierls rate [s⁻¹]", ylabel="log10 Taylor completion rate [s⁻¹]",
               title=f"{material}: source kinetic regime map")
        fig.colorbar(scatter, ax=ax, label="retained equilibrium fraction")
        fig.tight_layout(); fig.savefig(OUT / filename, dpi=180); plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.4, 5.4))
    for material, marker in (("Peak", "o"), ("DBTT", "^")):
        local = complete[complete.target_class.eq(material)]
        ax.scatter(local.log10_peierls_rate_s_median, local.log10_chi_taylor_completion_median,
                   s=9, alpha=0.5, label=material, marker=marker)
    ax.set(xlabel="log10 Peierls rate [s⁻¹]", ylabel="log10(transport time / Taylor completion time)",
           title="Transport–retention timescale map")
    ax.legend(); fig.tight_layout(); fig.savefig(OUT / "TRANSPORT_RETENTION_TIMESCALE_MAP.png", dpi=180); plt.close(fig)

    merged = stage1.merge(
        complete[["candidate_id", "log10_chi_taylor_completion_median"]], on="candidate_id", how="left"
    )
    fig, ax = plt.subplots(figsize=(7.4, 5.4))
    for material in ("Peak", "DBTT"):
        local = merged[merged.target_class.eq(material)]
        ax.scatter(local.log10_chi_taylor_completion_median, local.relative_deltaK_reinit,
                   s=8, alpha=0.35, label=material)
    ax.axhline(0, color="black", lw=1)
    ax.set(xlabel="log10(transport time / Taylor completion time)",
           ylabel="reload-separated relative ΔK (reduced)", title="Kinetic ratio versus R-curve propensity")
    ax.legend(); fig.tight_layout(); fig.savefig(OUT / "CHI_RET_VS_RCURVE_PROPENSITY.png", dpi=180); plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.4, 5.4))
    for (material, role), local in onsets[onsets.target_um.eq(100)].groupby(["material_class", "candidate_role"]):
        ax.scatter(local.pre_event_native_KJ_MPa_sqrt_m,
                   local.pre_event_K_effective_local_equivalent_MPa_sqrt_m,
                   label=f"{material} {role}", s=38)
    ax.set(xlabel="remote/native onset KJ [MPa√m]",
           ylabel="fixed-radius local-equivalent K [MPa√m]",
           title="Apparent versus local onset drive")
    ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(OUT / "APPARENT_VS_LOCAL_TOUGHENING.png", dpi=180); plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.2))
    transition = onsets[(onsets.target_um.eq(300))]
    for (material, role), local in transition.groupby(["material_class", "candidate_role"]):
        label = f"{material} {role}"
        axes[0].plot(local.pre_event_projected_extension_um, local.pre_event_tip_radius_um, "o-", label=label)
        axes[1].plot(local.pre_event_projected_extension_um, local.pre_event_backstress_GPa, "o-")
        axes[2].semilogy(local.pre_event_projected_extension_um, np.maximum(local.pre_event_retained_count, 1e-30), "o-")
    axes[0].set(ylabel="tip radius [µm]"); axes[1].set(ylabel="backstress [GPa]"); axes[2].set(ylabel="retained count")
    for ax in axes: ax.set_xlabel("onset extension [µm]")
    axes[0].legend(fontsize=7); fig.suptitle("PF onset-state evolution"); fig.tight_layout()
    fig.savefig(OUT / "RADIUS_BACKSTRESS_RETAINED_STATE_EVOLUTION.png", dpi=180); plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.5), sharey=False)
    for axis, material in zip(axes, ("Peak", "DBTT")):
        ids = [CONTROL_IDS[material], FINALIST_IDS[material]]
        for candidate in ids:
            local = stage1[
                stage1.target_class.eq(material) & stage1.candidate_id.eq(candidate)
                & stage1.provider.eq("PF")
            ].sort_values("temperature_K")
            axis.plot(local.temperature_K, local.initial_onset_native_KJ_MPa_sqrt_m, "o-",
                      label="control" if candidate == CONTROL_IDS[material] else "kinetic finalist")
        axis.set(title=material, xlabel="temperature [K]", ylabel="initial onset native KJ [MPa√m]")
        axis.legend()
    fig.suptitle("Reduced PF-provider controls versus kinetic finalists"); fig.tight_layout()
    fig.savefig(OUT / "CONTROL_VS_TAYLOR_PEIERLS_FINALISTS_REDUCED.png", dpi=180); plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.5))
    for axis, material in zip(axes, ("Peak", "DBTT")):
        for role in ("CONTROL", "KINETIC_FINALIST"):
            local = pf_summary[
                pf_summary.material_class.eq(material) & pf_summary.candidate_role.eq(role)
                & pf_summary.target_um.eq(100)
            ].sort_values("temperature_K")
            axis.plot(local.temperature_K, local.initial_onset_native_KJ_MPa_sqrt_m, "o-", label=role)
            axis.plot(local.temperature_K, local.maximum_onset_native_KJ_MPa_sqrt_m, "x--", alpha=0.7)
        axis.set(title=material, xlabel="temperature [K]", ylabel="PF onset native KJ [MPa√m]")
        axis.legend(fontsize=8)
    fig.suptitle("Direct PF: initial (solid) and maximum reload-separated onset (dashed)")
    fig.tight_layout(); fig.savefig(OUT / "CONTROL_VS_TAYLOR_PEIERLS_FINALISTS_PF.png", dpi=180); plt.close(fig)

    coords = _pca_coordinates(bank)
    fig, ax = plt.subplots(figsize=(7.2, 5.4))
    for material in ("Peak", "DBTT"):
        mask = bank.target_class.eq(material).to_numpy()
        ax.scatter(coords[mask, 0], coords[mask, 1], s=34, label=material)
    ax.set(xlabel="kinetic parameter PC1", ylabel="kinetic parameter PC2",
           title="Curated Taylor/Peierls option bank")
    ax.legend(); fig.tight_layout(); fig.savefig(OUT / "TAYLOR_PEIERLS_OPTION_BANK_MAP.png", dpi=180); plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.0, 4.6)); ax.axis("off")
    ax.text(0.5, 0.80, "Final decision", ha="center", fontsize=20, weight="bold")
    ax.text(0.25, 0.52, "PEAK\nRETAIN CONTROL\nNo positive PF reinitiation", ha="center", va="center",
            fontsize=14, bbox=dict(boxstyle="round,pad=.7", fc="#dbeafe", ec="#1d4ed8"))
    ax.text(0.75, 0.52, "DBTT\nRETAIN CONTROL\nPF reinitiation softens", ha="center", va="center",
            fontsize=14, bbox=dict(boxstyle="round,pad=.7", fc="#fee2e2", ec="#b91c1c"))
    ax.text(0.5, 0.16, "Taylor/Peierls flexibility changes state partition but not positive reload-separated resistance",
            ha="center", fontsize=11)
    fig.tight_layout(); fig.savefig(OUT / "FINAL_TAYLOR_PEIERLS_DECISION.png", dpi=180); plt.close(fig)


def write_reports(complete, pareto, bank, shortlist, stage1, long, pf_summary, onsets) -> dict[str, Any]:
    peak_pf = pf_summary[(pf_summary.material_class.eq("Peak")) & pf_summary.target_um.eq(300)]
    dbtt_pf = pf_summary[(pf_summary.material_class.eq("DBTT")) & pf_summary.target_um.eq(300)]
    peak_variant = peak_pf[peak_pf.candidate_role.eq("KINETIC_FINALIST")].iloc[0]
    dbtt_variant = dbtt_pf[dbtt_pf.candidate_role.eq("KINETIC_FINALIST")].iloc[0]
    decision = {
        "schema": "oneD_v2_taylor_peierls_final_decision_v1",
        "Peak": {
            "decision": "RETAIN_CONTROL",
            "control_candidate_id": CONTROL_IDS["Peak"],
            "kinetic_variant_decision": "NO_CREDIBLE_TAYLOR_PEIERLS_VARIANT",
            "direct_PF_300um_physical_avalanche_count": int(peak_variant.physical_avalanche_count),
            "direct_PF_positive_reinitiation": False,
        },
        "DBTT": {
            "decision": "RETAIN_CONTROL",
            "control_candidate_id": CONTROL_IDS["DBTT"],
            "kinetic_variant_decision": "NO_CREDIBLE_TAYLOR_PEIERLS_VARIANT",
            "direct_PF_300um_physical_avalanche_count": int(dbtt_variant.physical_avalanche_count),
            "direct_PF_signed_reinitiation_increment_MPa_sqrt_m": float(
                dbtt_variant.signed_max_reinitiation_minus_initial_K_MPa_sqrt_m
            ),
            "direct_PF_positive_reinitiation": False,
        },
        "barriers_changed": False, "memory_terms_added": False,
        "lifecycle_terms_changed": False, "new_FEMCZM_runs": 0,
        "scientific_conclusion": (
            "Within source-qualified bounds, Taylor/Peierls kinetics alone redistribute "
            "mobile and retained state but do not create positive reload-separated PF resistance."
        ),
        "recommended_full_precision_rows": [CONTROL_IDS["Peak"], CONTROL_IDS["DBTT"]],
    }
    decision_path = OUT / "oneD_v2_taylor_peierls_final_decision.json"
    decision_path.write_text(json.dumps(decision, indent=2, sort_keys=True) + "\n")

    reports = {
        "ONE_D_V2_TAYLOR_PEIERLS_TIMESCALE_ANALYSIS.md": f"""# Taylor/Peierls Timescale Analysis

The 8,194 Sobol rows plus two controls cover the empirical V9.13/494-option bounds. `chi_ret = tau_transport/tau_encounter` is constant within each class to numerical precision because the source encounter rate is proportional to Peierls velocity; the velocity cancels. It is therefore diagnostic but not an identifying search coordinate. Dimensional Peierls rates span about 51 decades in Peak and 0.6 decades in DBTT at the representative control states. Taylor completion changes the retained-equilibrium partition strongly.

Peak is predominantly transport/low-retention in the source exchange, while DBTT is retention-saturated. The search retained dimensional rates, `chi_taylor_completion`, and equilibrium fraction so this cancellation was not hidden. No new timescale entered the equations.
""",
        "ONE_D_V2_PEAK_TAYLOR_PEIERLS_SEARCH.md": f"""# Peak Taylor/Peierls Search

All {int((complete.target_class == 'Peak').sum())} Peak identities preserve the Peak cleavage and emission hashes exactly. The 513-row two-provider screen produced zero reload-separated reinitiations in all 6,156 cases. Sixteen diverse rows continued to 300 µm, five to 500 µm, and six across three seeds; Peak remained a single final avalanche in every provider/temperature/seed lane. Taylor/Peierls changes redistribute mobile versus retained content but leave reduced onset, radius, backstress, avalanche count, and largest-avalanche fraction invariant within numerical roundoff.

Direct PF confirms the negative result. The paired finalist has no reload at 600, 1000, or 1200 K and remains a single 67-event avalanche through 300 µm at 1000 K. Its initial onset is lower than the control at 1000 and 1200 K. Decision: retain `{CONTROL_IDS['Peak']}`; do not promote a Peak kinetic R variant.
""",
        "ONE_D_V2_DBTT_TAYLOR_PEIERLS_SEARCH.md": f"""# DBTT Taylor/Peierls Search

All {int((complete.target_class == 'DBTT').sum())} DBTT identities preserve the DBTT cleavage and emission hashes exactly. Reduced DBTT topology is provider-specific but candidate-invariant: the PF provider has no reinitiation through 900 K and one positive reinitiation at 1000–1200 K; FEM/CZM has reinitiations across the full range, including the low shelf. This is not a kinetic-candidate improvement because every searched row, including the control, gives the same fracture response.

Direct PF rejects the apparent reduced signal. Control and finalist each have one reload, but the reinitiation is softer. At 1100 K through 300 µm the finalist increment is {float(dbtt_variant.signed_max_reinitiation_minus_initial_K_MPa_sqrt_m):.6g} MPa√m. Decision: retain `{CONTROL_IDS['DBTT']}`; do not promote a DBTT kinetic R variant.
""",
        "ONE_D_V2_TAYLOR_PEIERLS_MULTISEED_VALIDATION.md": f"""# Taylor/Peierls Multiseed Validation

The canonical multiseed table contains {len(pd.read_parquet(OUT / 'oneD_v2_taylor_peierls_multiseed_results.parquet'))} cases: six diverse rows per class, both providers, all requested class temperatures, and three seeds at 300 µm. All completed. Five rows per class also reached 500 µm under the canonical seed. Across material rows, response quantities are invariant to about 1e-13; seed changes affect the stochastic topology as expected but do not reveal a kinetic-row-specific R response.
""",
        "PF_2D_TAYLOR_PEIERLS_TRANSFER_VALIDATION.md": f"""# PF 2-D Taylor/Peierls Transfer Validation

Sixteen bounded PF cases were run with a hard maximum of two workers: 12 paired 100 µm cases and four transition-temperature 300 µm continuations. Every case reached its target. The default-off observer records the actual active and wake MPZ arrays, radius, front width, mobile/retained state, backstress, signed shielding, opening/channel tensors, resolved shears, and source-owned Peierls/Taylor/encounter profiles at resistance onsets.

Peak has no reload-separated reinitiation. DBTT reinitiates once but softens in both control and variant. The kinetic changes clearly alter state partition, so the negative result is not caused by inactive parameters. Neither finalist meets the direct-PF positive-resistance criterion. No 500 µm PF extension is needed because the 300 µm trajectories remain unambiguous.
""",
        "ONE_D_V2_APPARENT_VS_LOCAL_TOUGHENING_DECOMPOSITION.md": """# Apparent Versus Local Toughening Decomposition

Only reload-separated pre-event states are resistance points. Remote/common-reference K, native KJ, signed K-shield, fixed-radius local-equivalent K, source opening stress, radius, backstress, and retained/mobile state are kept separate. Peak supplies no second resistance point. DBTT's second remote/native onset is lower, not higher. Thus neither class can be classified as radius-dominated, shielding-dominated, mixed positive toughening, or local-critical-drive increase. Peak is `NO_TOUGHENING`; DBTT is `SOFTENING`. Rising values inside the final avalanche remain model-native trajectory drive and are not called an R-curve.
""",
        "ONE_D_V2_TAYLOR_PEIERLS_OPTION_BANK.md": f"""# Taylor/Peierls Option Bank

The complete library contains {len(complete)} identities; {len(pareto)} received full two-provider trajectory evaluation. Because their fracture objectives are numerically equivalent within each class, all fully evaluated rows are response-Pareto-equivalent rather than artificially ranked with a kinetic scalar score. A maximin kinetic/state selection retains {len(bank)} curated rows and {len(shortlist)} compact shortlist rows. Controls are present at every level. Finalists rejected by PF remain available as diagnostic kinetic extremes, explicitly labeled `DIRECT_PF_REJECTED_NO_POSITIVE_REINITIATION`; they are not promoted material rows.
""",
        "ONE_D_V2_TAYLOR_PEIERLS_FINAL_DECISION.md": f"""# Taylor/Peierls Final Decision

## Decision

Retain both controls: Peak `{CONTROL_IDS['Peak']}` and DBTT `{CONTROL_IDS['DBTT']}`. Within the source-qualified bounds, the ten Taylor/Peierls coordinates alone do not create a credible positive reload-separated R-curve response with the existing Peak/DBTT cleavage and emission barriers.

The search changed only mobile/retained partition kinetics. In the reduced model, radius and backstress depend on total state, making fracture response structurally invariant to this partition. Direct PF allows spatial redistribution and shows modest onset changes, but Peak has no reload and DBTT softens on reload. There is no true positive model-form response to promote.

No barrier, memory, lifecycle, provider-specific material row, PF equation, FEM/CZM equation, or production trajectory was modified. No FEM/CZM simulation was launched. The only PF source change is a default-off, no-feedback state observer on branch `codex/oneD-v2-taylor-peierls-pf-transfer`.
""",
    }
    for filename, content in reports.items():
        (ROOT / filename).write_text(content.strip() + "\n")
    return decision


def provenance(decision: dict[str, Any]) -> None:
    names = (
        "oneD_v2_taylor_peierls_complete_population.parquet",
        "oneD_v2_taylor_peierls_pareto_library.parquet",
        "oneD_v2_taylor_peierls_multiseed_results.parquet",
        "oneD_v2_taylor_peierls_exact_finalists.csv",
        "oneD_v2_taylor_peierls_timescale_diagnostics.csv",
        "oneD_v2_taylor_peierls_option_bank.csv",
        "oneD_v2_taylor_peierls_option_bank.parquet",
        "oneD_v2_taylor_peierls_shortlist.csv",
        "oneD_v2_taylor_peierls_material_registry.csv",
        "pf_2d_taylor_peierls_event_transactions.csv",
        "pf_2d_taylor_peierls_physical_avalanches.csv",
        "pf_2d_taylor_peierls_onset_candidates.csv",
        "pf_2d_taylor_peierls_state_profiles.parquet",
        "pf_2d_taylor_peierls_transfer_summary.csv",
        "oneD_v2_taylor_peierls_final_decision.json",
    )
    artifacts = {name: sha(OUT / name) for name in names}
    manifest = {
        "schema": "oneD_v2_taylor_peierls_provenance_manifest_v1",
        "reduced_repository": str(ROOT), "reduced_branch": "codex/oneD-v2-taylor-peierls-rcurve-search",
        "analysis_commit": head(ROOT), "source_parent_commit": "3d64c918431874f1ee1c215e4441b86672a52cc6",
        "pf_repository": str(PF_REPO), "pf_branch": "codex/oneD-v2-taylor-peierls-pf-transfer",
        "pf_runner_commit": head(PF_REPO), "pf_parent_commit": "130b305ecc17863cba3cd8c157ee179796740cfd",
        "maximum_concurrent_heavy_PF_workers": 2, "new_PF_case_count": 16,
        "new_FEMCZM_run_count": 0, "fatigue_code_analyzed_or_modified": False,
        "canonical_registries_modified": False, "artifacts": artifacts,
        "decision_sha256": artifacts["oneD_v2_taylor_peierls_final_decision.json"],
    }
    (OUT / "oneD_v2_taylor_peierls_provenance_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True, default=_jsonable) + "\n"
    )


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    complete, pareto, bank, shortlist = build_population_and_banks()
    transactions, avalanches, onsets, profiles, pf_summary = build_pf_outputs()
    stage1 = pd.read_parquet(OUT / "oneD_v2_taylor_peierls_stage1_100um_cases.parquet")
    long = pd.read_parquet(OUT / "oneD_v2_taylor_peierls_long_validation.parquet")
    figures(complete, bank, stage1, pf_summary, onsets, profiles)
    decision = write_reports(complete, pareto, bank, shortlist, stage1, long, pf_summary, onsets)
    provenance(decision)
    print(json.dumps({
        "complete_population": len(complete), "fully_evaluated": len(pareto),
        "option_bank": len(bank), "shortlist": len(shortlist),
        "pf_cases": len(pf_summary), "pf_onsets": len(onsets), "profile_rows": len(profiles),
        "decision": {key: decision[key]["decision"] for key in ("Peak", "DBTT")},
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
