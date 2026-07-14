#!/usr/bin/env python3
"""Promote v9.10 broad-search candidates to the unified spatial MPZ."""
from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
import pandas as pd

from arrhenius_fracture import sharp_front as sf
from arrhenius_fracture.moving_process_zone_v910 import MovingProcessZoneState
from fit_mpz_three_classes import simulate
import promote_mpz_v9_9_spatial as v99


def parse_floats(text: str) -> list[float]:
    return [float(x) for x in str(text).replace(",", " ").split() if x]


def finite(row: pd.Series, name: str, default: float) -> float:
    try:
        value = float(row.get(name, default))
    except (TypeError, ValueError):
        return float(default)
    return value if np.isfinite(value) else float(default)


def material_row(candidate: pd.Series, mpz_length_m: float, mpz_n_bins: int) -> pd.Series:
    proxy = candidate.copy()
    proxy["continuation_candidate_id"] = str(candidate.get("candidate_id", "candidate"))
    proxy["barrier_scale"] = 1.0
    proxy["recovery_rate_s"] = finite(candidate, "retained_recovery_rate_s", 1.0e-5)
    proxy["mobile_fraction"] = 0.01
    return v99.material_row(proxy, mpz_length_m, mpz_n_bins)


def run_spatial(
    candidate: pd.Series,
    T_K: float,
    opt: SimpleNamespace,
    *,
    mpz_length_m: float,
    mpz_n_bins: int,
) -> tuple[dict[str, Any], Any, dict[str, float]]:
    row = material_row(candidate, mpz_length_m, mpz_n_bins)
    original_build = sf.build_engine
    holder: dict[str, Any] = {}

    def patched_build(args, material):
        eng = original_build(args, material)
        cfg = copy.deepcopy(eng.mpz_config)
        cfg.length_m = float(mpz_length_m)
        cfg.n_bins = int(mpz_n_bins)
        cfg.n_systems = 2
        cfg.source_sites_per_system = finite(candidate, "source_sites_per_system", 200.0)
        cfg.source_recovery_rate_s = 0.0
        cfg.source_refresh_length_m = finite(candidate, "source_refresh_length_um", 0.25) * 1.0e-6
        cfg.pt_emit_G00_eV = finite(candidate, "emit_G00_eV", 1.5)
        cfg.pt_emit_gT_eV_per_K = finite(candidate, "emit_gT_eV_per_K", 0.0)
        cfg.pt_emit_sigc0_Pa = finite(candidate, "emit_sigc0_GPa", 2.5) * 1.0e9
        cfg.pt_emit_sT_Pa_per_K = finite(candidate, "emit_sT_GPa_per_K", 0.0) * 1.0e9
        cfg.pt_emit_exp_a = finite(candidate, "emit_exp_a", 0.2)
        cfg.pt_emit_exp_n = finite(candidate, "emit_exp_n", 1.0)
        cfg.pt_emit_floor_frac = finite(candidate, "emit_floor_frac", 0.02)
        emit0 = max(cfg.pt_emit_G00_eV, 1.0e-12)
        hp = finite(candidate, "peierls_H0_eV", 1.0)
        ht = max(finite(candidate, "taylor_H0_eV", hp), hp)
        cfg.pt_peierls_energy_ratio = hp / emit0
        cfg.pt_peierls_entropy_ratio = finite(candidate, "peierls_activation_entropy_kB", 0.0)
        cfg.pt_peierls_nu0_s = 1.0e12
        cfg.pt_taylor_energy_ratio = ht / emit0
        cfg.pt_taylor_entropy_ratio = finite(candidate, "taylor_activation_entropy_kB", 0.0)
        cfg.pt_taylor_nu0_s = 1.0e11
        cfg.pt_taylor_corr_rho_c = finite(candidate, "taylor_corr_rho_c_m2", 1.0e14)
        cfg.pt_taylor_renewal_time_s = 1.0
        cfg.pt_taylor_m_exponent = 1.0
        cfg.pt_taylor_m_scale = finite(candidate, "taylor_corr_scale", 1.0)
        cfg.pt_taylor_m_cap = float("inf")
        cfg.pt_mobile_fraction = 0.01
        cfg.pt_mobile_saturation_density_m2 = float("inf")
        cfg.pt_mobile_density_floor_m2 = 0.0
        cfg.pt_jump_length_min_m = 0.0
        cfg.pt_taylor_phi_max = float("inf")
        cfg.pt_encounter_efficiency = finite(candidate, "encounter_efficiency", 1.0)
        cfg.retained_recovery_nu0_s = finite(candidate, "retained_recovery_rate_s", 0.0)
        cfg.retained_recovery_barrier_eV = 0.0
        cfg.retained_recovery_activation_volume_b3 = 0.0
        cfg.mobile_recovery_rate_s = 0.0
        cfg.pair_annihilation_rate_per_count_s = 0.0
        cfg.blunting_length_m = max(0.5e-6, 0.5 * finite(row, "r_pz_m", 1.0e-6))
        cfg.mobile_shield_fraction = 0.0
        eng.mpz_config = cfg
        eng.mpz_state = MovingProcessZoneState(cfg)
        eng._sync_compat()
        holder["engine"] = eng
        return eng

    sf.build_engine = patched_build
    try:
        result = simulate(row, float(T_K), opt)
    finally:
        sf.build_engine = original_build
    eng = holder.get("engine")
    if eng is None:
        raise RuntimeError("v9.10 unified spatial engine was not captured")
    Kdiag = result.get("K_plateau", result.get("K_init", 0.0))
    if not np.isfinite(Kdiag):
        Kdiag = result.get("K_init", 0.0)
    sigma = eng.sigma_tip(max(float(Kdiag), 0.0) * 1.0e6)
    diag = eng.mpz_state.evolve(0.0, float(T_K), sigma, eng.b, 0.0)
    load_time = max(float(result.get("K_init", 0.0)) / max(opt.Kdot, 1.0e-30), 0.0)
    traverse = float(diag.get("glide_velocity_m_s", 0.0)) * load_time / max(float(eng.f.r0), 1.0e-30)
    extra = {
        "final_peierls_rate_s": float(diag.get("peierls_rate_s", 0.0)),
        "final_taylor_completion_rate_s": float(diag.get("taylor_completion_rate_s", 0.0)),
        "final_encounter_rate_s": float(diag.get("encounter_rate_s", 0.0)),
        "final_glide_velocity_m_s": float(diag.get("glide_velocity_m_s", 0.0)),
        "final_peierls_traverse_number": traverse,
        "final_retained_count": float(eng.mpz_state.retained_count),
        "final_mobile_count": float(eng.mpz_state.mobile_count),
        "final_K_shield_MPa_sqrt_m": float(eng.mpz_state.shielding_K(eng.G, eng.nu, eng.b) / 1.0e6),
        "final_available_site_fraction": float(eng.mpz_state.available_site_fraction),
        "unified_transport_retention_active": float(diag.get("unified_transport_retention_active", 0.0)),
        "legacy_trap_barrier_active": float(diag.get("legacy_trap_barrier_active", 1.0)),
    }
    return result, eng, extra


def candidate_summary(group: pd.DataFrame) -> dict[str, Any]:
    ordered = group.sort_values("T_K")
    target_class = str(ordered.target_class.iloc[0])
    low = ordered[ordered.T_K <= 700.0]
    high = ordered[ordered.T_K >= 900.0]
    plateau_rise = float(high.K_plateau.median() - low.K_plateau.median())
    init_rise = float(high.K_init.median() - low.K_init.median())
    high_dkr = float(high.delta_KR.median())
    low_dkr = float(low.delta_KR.median())
    if target_class == "ceramic":
        accepted = bool(ordered.completed.all() and ordered.delta_KR.max() <= 3.0)
        reason = "ceramic_small_Rcurve_persists" if accepted else "ceramic_spatial_failure"
    elif target_class == "weakT":
        flat = float(ordered.K_plateau.max() - ordered.K_plateau.min()) <= 6.0
        fast = float(ordered.final_peierls_traverse_number.min()) >= 1.0
        moderate = 2.0 <= float(ordered.delta_KR.median()) <= 9.0
        accepted = bool(ordered.completed.all() and flat and fast and moderate)
        reason = "weakT_FCC_unified_persists" if accepted else "weakT_spatial_failure"
    else:
        accepted = bool(ordered.completed.all() and plateau_rise >= 15.0 and high_dkr >= 3.0)
        reason = "DBTT_unified_persists_with_crack_growth" if accepted else "DBTT_spatial_failure"
    return {
        "candidate_id": str(ordered.candidate_id.iloc[0]),
        "target_class": target_class,
        "all_temperatures_completed": bool(ordered.completed.all()),
        "K_init_high_minus_low": init_rise,
        "K_plateau_high_minus_low": plateau_rise,
        "high_temperature_delta_KR": high_dkr,
        "low_temperature_delta_KR": low_dkr,
        "minimum_peierls_traverse_number": float(ordered.final_peierls_traverse_number.min()),
        "maximum_final_K_shield_MPa_sqrt_m": float(ordered.final_K_shield_MPa_sqrt_m.max()),
        "accepted_after_spatial_growth": accepted,
        "acceptance_reason": reason,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest-root", type=Path, default=Path("runs/mpz_v9_10_unified_global_search_v1"))
    ap.add_argument("--classes", default="ceramic weakT DBTT")
    ap.add_argument("--temperatures", default="300 700 900 1200")
    ap.add_argument("--max-per-class", type=int, default=3)
    ap.add_argument("--target-extension-um", type=float, default=500.0)
    ap.add_argument("--da-um", type=float, default=5.0)
    ap.add_argument("--dK", type=float, default=0.25)
    ap.add_argument("--Kdot", type=float, default=0.005)
    ap.add_argument("--Kmax", type=float, default=80.0)
    ap.add_argument("--mpz-length-um", type=float, default=100.0)
    ap.add_argument("--mpz-n-bins", type=int, default=200)
    ap.add_argument("--out", type=Path, default=Path("runs/mpz_v9_10_unified_spatial_promotion_v1"))
    a = ap.parse_args()

    classes = str(a.classes).replace(",", " ").split()
    manifests = []
    for class_name in classes:
        path = a.manifest_root / class_name / "spatial_promotion_manifest.csv"
        if not path.exists():
            raise SystemExit(f"v9.10 promotion manifest not found: {path}")
        manifests.append(pd.read_csv(path))
    manifest = pd.concat(manifests, ignore_index=True)
    selected = (
        manifest.sort_values(["target_class", "objective"])
        .groupby("target_class", as_index=False)
        .head(a.max_per_class)
        .reset_index(drop=True)
    )
    if selected.empty:
        raise SystemExit("no v9.10 candidates selected for spatial promotion")
    temperatures = parse_floats(a.temperatures)
    opt = SimpleNamespace(
        dK=float(a.dK), Kdot=float(a.Kdot),
        n_advances=int(round(a.target_extension_um / a.da_um)) + 1,
        Kmax=float(a.Kmax), da_um=float(a.da_um),
        early_window_um=(20.0, min(100.0, 0.5 * a.target_extension_um)),
        plateau_window_um=(0.70 * a.target_extension_um, a.target_extension_um),
        target_dB_substep=0.25, target_emission_hazard_substep=1.0,
        source_active_fraction_min=1.0e-4, min_substep_fraction=1.0e-8,
        max_substeps=2_000_000, objective_mode="rcurve",
    )
    out = a.out.resolve(); out.mkdir(parents=True, exist_ok=True)
    metric_rows: list[dict[str, Any]] = []
    event_rows: list[dict[str, Any]] = []
    total = len(selected) * len(temperatures); count = 0
    for _, candidate in selected.iterrows():
        for T in temperatures:
            count += 1
            result, _, extra = run_spatial(
                candidate, T, opt,
                mpz_length_m=a.mpz_length_um * 1.0e-6,
                mpz_n_bins=a.mpz_n_bins,
            )
            rec = {
                "candidate_id": str(candidate.candidate_id),
                "target_class": str(candidate.target_class),
                "analytical_objective": float(candidate.objective),
                **{k: v for k, v in result.items() if k != "events"},
                **extra,
                "status": "V9_10_UNIFIED_SPATIAL_NOT_2D_VALIDATED",
            }
            metric_rows.append(rec)
            for event in result.get("events", []):
                event_rows.append({
                    "candidate_id": str(candidate.candidate_id),
                    "target_class": str(candidate.target_class),
                    "T_K": float(T), **event,
                })
            print(
                f"evaluated {count}/{total} candidate={candidate.candidate_id} T={T:g} "
                f"completed={result['completed']} Kinit={result['K_init']:.4g} "
                f"Kplateau={result['K_plateau']:.4g} dKR={result['delta_KR']:.4g}",
                flush=True,
            )
    metrics = pd.DataFrame(metric_rows)
    events = pd.DataFrame(event_rows)
    summaries = pd.DataFrame([
        candidate_summary(group)
        for _, group in metrics.groupby("candidate_id", sort=False)
    ])
    metrics.to_csv(out / "unified_spatial_temperature_metrics.csv", index=False)
    events.to_csv(out / "unified_spatial_event_history.csv", index=False)
    summaries.to_csv(out / "unified_spatial_summary.csv", index=False)
    selected.to_csv(out / "selected_candidates.csv", index=False)
    report = {
        "n_candidates": int(len(summaries)),
        "n_temperature_runs": int(len(metrics)),
        "n_accepted": int(summaries.accepted_after_spatial_growth.sum()),
        "accepted_by_class": summaries[summaries.accepted_after_spatial_growth].groupby("target_class").size().to_dict(),
        "target_extension_um": float(a.target_extension_um),
        "status": "V9_10_UNIFIED_SPATIAL_COMPLETE_NOT_2D_VALIDATED",
        "output": str(out),
    }
    (out / "unified_spatial_report.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
