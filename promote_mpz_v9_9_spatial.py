#!/usr/bin/env python3
"""Promote v9.9 continuation candidates to finite crack-growth MPZ runs.

The analytical continuation is only a screening stage.  This script replaces
the engine's active state with the v9.9 independent-entropy spatial MPZ adapter
and runs the existing reduced sharp-front crack-growth calculation.  It records
full event histories and evaluates whether the response class persists as the
crack advances.
"""
from __future__ import annotations

import argparse
import copy
import json
import math
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
import pandas as pd

from arrhenius_fracture import sharp_front as sf
from arrhenius_fracture.moving_process_zone_v99 import MovingProcessZoneState
from fit_mpz_three_classes import simulate


def parse_floats(text: str) -> list[float]:
    return [float(x) for x in str(text).replace(",", " ").split() if x]


def finite(row: pd.Series, name: str, default: float) -> float:
    try:
        value = float(row.get(name, default))
    except (TypeError, ValueError):
        return float(default)
    return value if np.isfinite(value) else float(default)


def shape_value(row: pd.Series, name: str, default: float) -> float:
    return finite(row, f"shape_{name}", finite(row, name, default))


def material_row(row: pd.Series, mpz_length_m: float, mpz_n_bins: int) -> pd.Series:
    """Construct the complete legacy row consumed by ``simulate``.

    The parser-created MPZ state is replaced immediately after engine build;
    the legacy PT scale fields here only keep parser construction valid.
    """
    emit0 = max(finite(row, "emit_G00_eV", 1.5), 1.0e-12)
    hp = finite(row, "peierls_H0_eV", 1.0)
    ht = max(finite(row, "taylor_H0_eV", hp), hp)
    source_refresh_m = finite(row, "source_refresh_length_um", 0.25) * 1.0e-6
    values: dict[str, Any] = {
        "candidate_id": str(row.get("continuation_candidate_id", "candidate")),
        "r_pz_m": 1.0e-6,
        "c_blunt": finite(row, "c_blunt", 1.0),
        "mpz_length_m": float(mpz_length_m),
        "mpz_n_bins": int(mpz_n_bins),
        "mpz_n_systems": 2,
        "mpz_source_sites_per_system": finite(row, "source_sites_per_system", 200.0),
        "mpz_source_recovery_rate_s": 0.0,
        "mpz_source_refresh_length_m": source_refresh_m,
        "mpz_shielding_factors": "1 1",
        "mpz_glide_barrier_eV": 0.8,
        "mpz_glide_activation_volume_b3": 8.0,
        "mpz_trap_barrier_eV": 0.65,
        "mpz_detrap_barrier_eV": 1.2,
        "mpz_retained_recovery_barrier_eV": 1.5,
        "mpz_pair_annihilation_rate_per_count_s": 0.0,
        "pt_peierls_energy_ratio": hp / emit0,
        "pt_peierls_entropy_ratio": finite(row, "peierls_activation_entropy_kB", -20.0),
        "pt_taylor_energy_ratio": ht / emit0,
        "pt_taylor_entropy_ratio": finite(row, "taylor_activation_entropy_kB", -20.0),
        "pt_taylor_corr_rho_c": finite(row, "taylor_corr_rho_c_m2", 1.0e14),
        "pt_taylor_renewal_time_s": 1.0,
        "pt_taylor_m_exponent": 1.0,
        "pt_taylor_m_scale": finite(row, "taylor_corr_scale", 1.0),
        "pt_taylor_m_cap": float("inf"),
        "pt_mobile_fraction": finite(row, "mobile_fraction", 0.01),
        "pt_mobile_saturation_density_m2": float("inf"),
        "cleave_G00_eV": finite(row, "cleave_G00_eV", 2.0),
        "cleave_gT_eV_per_K": finite(row, "cleave_gT_eV_per_K", 0.0),
        "cleave_sigc0_GPa": finite(row, "cleave_sigc0_GPa", 4.0),
        "cleave_sT_GPa_per_K": shape_value(row, "cleave_sT_GPa_per_K", 0.0),
        "cleave_exp_a": shape_value(row, "cleave_exp_a", 0.2),
        "cleave_exp_n": shape_value(row, "cleave_exp_n", 1.0),
        "cleave_floor_frac": shape_value(row, "cleave_floor_frac", 0.02),
        "cleave_S_hs_kB": 0.0,
        "emit_G00_eV": emit0,
        "emit_gT_eV_per_K": finite(row, "emit_gT_eV_per_K", 0.0),
        "emit_sigc0_GPa": finite(row, "emit_sigc0_GPa", 2.5),
        "emit_sT_GPa_per_K": shape_value(row, "emit_sT_GPa_per_K", 0.0),
        "emit_exp_a": shape_value(row, "emit_exp_a", 0.2),
        "emit_exp_n": shape_value(row, "emit_exp_n", 1.0),
        "emit_floor_frac": shape_value(row, "emit_floor_frac", 0.02),
    }
    return pd.Series(values)


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
        cfg.source_sites_per_system = finite(
            candidate, "source_sites_per_system", 200.0
        )
        cfg.source_recovery_rate_s = 0.0
        cfg.source_refresh_length_m = finite(
            candidate, "source_refresh_length_um", 0.25
        ) * 1.0e-6
        cfg.pt_emit_G00_eV = finite(candidate, "emit_G00_eV", 1.5)
        cfg.pt_emit_gT_eV_per_K = finite(candidate, "emit_gT_eV_per_K", 0.0)
        cfg.pt_emit_sigc0_Pa = finite(candidate, "emit_sigc0_GPa", 2.5) * 1.0e9
        cfg.pt_emit_sT_Pa_per_K = shape_value(
            candidate, "emit_sT_GPa_per_K", 0.0
        ) * 1.0e9
        cfg.pt_emit_exp_a = shape_value(candidate, "emit_exp_a", 0.2)
        cfg.pt_emit_exp_n = shape_value(candidate, "emit_exp_n", 1.0)
        cfg.pt_emit_floor_frac = shape_value(candidate, "emit_floor_frac", 0.02)
        emit0 = max(cfg.pt_emit_G00_eV, 1.0e-12)
        hp = finite(candidate, "peierls_H0_eV", 1.0)
        ht = max(finite(candidate, "taylor_H0_eV", hp), hp)
        cfg.pt_peierls_energy_ratio = hp / emit0
        # v9.9 adapter interprets these two slots as S*/k_B.
        cfg.pt_peierls_entropy_ratio = finite(
            candidate, "peierls_activation_entropy_kB", -20.0
        )
        cfg.pt_peierls_nu0_s = finite(candidate, "peierls_nu0_s", 1.0e12)
        cfg.pt_taylor_energy_ratio = ht / emit0
        cfg.pt_taylor_entropy_ratio = finite(
            candidate, "taylor_activation_entropy_kB", -20.0
        )
        cfg.pt_taylor_nu0_s = finite(candidate, "taylor_nu0_s", 1.0e11)
        cfg.pt_taylor_corr_rho_c = finite(
            candidate, "taylor_corr_rho_c_m2", 1.0e14
        )
        cfg.pt_taylor_renewal_time_s = 1.0
        cfg.pt_taylor_m_exponent = 1.0
        cfg.pt_taylor_m_scale = finite(candidate, "taylor_corr_scale", 1.0)
        cfg.pt_taylor_m_cap = float("inf")
        cfg.pt_mobile_fraction = finite(candidate, "mobile_fraction", 0.01)
        cfg.pt_mobile_saturation_density_m2 = float("inf")
        cfg.pt_mobile_density_floor_m2 = 0.0
        cfg.pt_jump_length_min_m = 0.0
        cfg.pt_taylor_phi_max = float("inf")
        cfg.retained_recovery_nu0_s = finite(candidate, "recovery_rate_s", 1.0e-5)
        cfg.retained_recovery_barrier_eV = 0.0
        cfg.retained_recovery_activation_volume_b3 = 0.0
        cfg.mobile_recovery_rate_s = 0.0
        cfg.pair_annihilation_rate_per_count_s = 0.0
        cfg.blunting_length_m = max(0.5 * finite(row, "r_pz_m", 1.0e-6), 1.0e-9)
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
        raise RuntimeError("v9.9 spatial MPZ engine was not captured")

    Kdiag = result.get("K_plateau", result.get("K_init", 0.0))
    if not np.isfinite(Kdiag):
        Kdiag = result.get("K_init", 0.0)
    sigma = eng.sigma_tip(max(float(Kdiag), 0.0) * 1.0e6)
    diag = eng.mpz_state.evolve(0.0, float(T_K), sigma, eng.b, 0.0)
    load_time = max(float(result.get("K_init", 0.0)) / max(opt.Kdot, 1.0e-30), 0.0)
    traverse = (
        float(diag.get("glide_velocity_m_s", 0.0))
        * load_time
        / max(float(eng.f.r0), 1.0e-30)
    )
    extra = {
        "final_peierls_rate_s": float(diag.get("peierls_rate_s", 0.0)),
        "final_taylor_completion_rate_s": float(
            diag.get("taylor_completion_rate_s", 0.0)
        ),
        "final_glide_velocity_m_s": float(diag.get("glide_velocity_m_s", 0.0)),
        "final_peierls_traverse_number": float(traverse),
        "final_retained_count": float(eng.mpz_state.retained_count),
        "final_mobile_count": float(eng.mpz_state.mobile_count),
        "final_K_shield_MPa_sqrt_m": float(
            eng.mpz_state.shielding_K(eng.G, eng.nu, eng.b) / 1.0e6
        ),
        "final_available_site_fraction": float(
            eng.mpz_state.available_site_fraction
        ),
    }
    return result, eng, extra


def candidate_summary(group: pd.DataFrame) -> dict[str, Any]:
    ordered = group.sort_values("T_K")
    target_class = str(ordered.target_class.iloc[0])
    low = ordered[ordered.T_K <= 700.0]
    high = ordered[ordered.T_K >= 900.0]
    plateau_rise = (
        float(high.K_plateau.median() - low.K_plateau.median())
        if not low.empty and not high.empty
        else float("nan")
    )
    init_rise = (
        float(high.K_init.median() - low.K_init.median())
        if not low.empty and not high.empty
        else float("nan")
    )
    high_dKR = float(high.delta_KR.median()) if not high.empty else float("nan")
    low_dKR = float(low.delta_KR.median()) if not low.empty else float("nan")
    if target_class == "ceramic":
        accepted = bool(ordered.completed.all() and ordered.delta_KR.max() <= 3.0)
        reason = "ceramic_small_Rcurve_persists" if accepted else "ceramic_spatial_failure"
    elif target_class == "weakT":
        flat = float(ordered.K_plateau.max() - ordered.K_plateau.min()) <= 6.0
        fast = float(ordered.final_peierls_traverse_number.min()) >= 1.0
        moderate = 2.0 <= float(ordered.delta_KR.median()) <= 9.0
        accepted = bool(ordered.completed.all() and flat and fast and moderate)
        reason = "weakT_fast_Peierls_persists" if accepted else "weakT_spatial_failure"
    else:
        accepted = bool(
            ordered.completed.all()
            and np.isfinite(plateau_rise)
            and plateau_rise >= 15.0
            and high_dKR >= 3.0
        )
        reason = "DBTT_persists_with_crack_growth" if accepted else "DBTT_spatial_failure"
    return {
        "continuation_candidate_id": str(ordered.continuation_candidate_id.iloc[0]),
        "target_class": target_class,
        "barrier_scale": float(ordered.barrier_scale.iloc[0]),
        "all_temperatures_completed": bool(ordered.completed.all()),
        "K_init_high_minus_low": init_rise,
        "K_plateau_high_minus_low": plateau_rise,
        "high_temperature_delta_KR": high_dKR,
        "low_temperature_delta_KR": low_dKR,
        "minimum_peierls_traverse_number": float(
            ordered.final_peierls_traverse_number.min()
        ),
        "maximum_final_K_shield_MPa_sqrt_m": float(
            ordered.final_K_shield_MPa_sqrt_m.max()
        ),
        "accepted_after_spatial_growth": accepted,
        "acceptance_reason": reason,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--manifest",
        type=Path,
        default=Path(
            "runs/mpz_v9_9_barrier_continuation_v1/spatial_promotion_manifest.csv"
        ),
    )
    ap.add_argument("--classes", default="ceramic weakT DBTT")
    ap.add_argument("--temperatures", default="300 700 900 1200")
    ap.add_argument("--max-per-class", type=int, default=2)
    ap.add_argument("--target-extension-um", type=float, default=500.0)
    ap.add_argument("--da-um", type=float, default=5.0)
    ap.add_argument("--dK", type=float, default=0.25)
    ap.add_argument("--Kdot", type=float, default=0.005)
    ap.add_argument("--Kmax", type=float, default=80.0)
    ap.add_argument("--mpz-length-um", type=float, default=100.0)
    ap.add_argument("--mpz-n-bins", type=int, default=200)
    ap.add_argument(
        "--out", type=Path, default=Path("runs/mpz_v9_9_spatial_promotion_v1")
    )
    a = ap.parse_args()

    if not a.manifest.exists():
        raise SystemExit(f"promotion manifest not found: {a.manifest}")
    manifest = pd.read_csv(a.manifest)
    classes = str(a.classes).replace(",", " ").split()
    manifest = manifest[manifest.target_class.astype(str).isin(classes)].copy()
    selected = (
        manifest.sort_values(["target_class", "objective"])
        .groupby("target_class", as_index=False)
        .head(a.max_per_class)
        .reset_index(drop=True)
    )
    if selected.empty:
        raise SystemExit("no v9.9 candidates selected for spatial promotion")
    temperatures = parse_floats(a.temperatures)
    opt = SimpleNamespace(
        dK=float(a.dK),
        Kdot=float(a.Kdot),
        n_advances=int(round(a.target_extension_um / a.da_um)) + 1,
        Kmax=float(a.Kmax),
        da_um=float(a.da_um),
        early_window_um=(20.0, min(100.0, 0.5 * a.target_extension_um)),
        plateau_window_um=(0.70 * a.target_extension_um, a.target_extension_um),
        target_dB_substep=0.25,
        target_emission_hazard_substep=1.0,
        source_active_fraction_min=1.0e-4,
        min_substep_fraction=1.0e-8,
        max_substeps=2_000_000,
        objective_mode="rcurve",
    )
    out = a.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    metric_rows: list[dict[str, Any]] = []
    event_rows: list[dict[str, Any]] = []
    total = len(selected) * len(temperatures)
    count = 0
    for _, candidate in selected.iterrows():
        for T in temperatures:
            count += 1
            result, _, extra = run_spatial(
                candidate,
                T,
                opt,
                mpz_length_m=a.mpz_length_um * 1.0e-6,
                mpz_n_bins=a.mpz_n_bins,
            )
            rec = {
                "continuation_candidate_id": str(candidate.continuation_candidate_id),
                "target_class": str(candidate.target_class),
                "barrier_scale": float(candidate.barrier_scale),
                "analytical_objective": float(candidate.objective),
                **{k: v for k, v in result.items() if k != "events"},
                **extra,
                "status": "V9_9_SPATIAL_MPZ_NOT_2D_VALIDATED",
            }
            metric_rows.append(rec)
            for event in result.get("events", []):
                event_rows.append(
                    {
                        "continuation_candidate_id": str(
                            candidate.continuation_candidate_id
                        ),
                        "target_class": str(candidate.target_class),
                        "barrier_scale": float(candidate.barrier_scale),
                        "T_K": float(T),
                        **event,
                    }
                )
            print(
                f"evaluated {count}/{total} candidate={candidate.continuation_candidate_id} "
                f"T={T:g} completed={result['completed']} "
                f"Kinit={result['K_init']:.4g} Kplateau={result['K_plateau']:.4g} "
                f"dKR={result['delta_KR']:.4g}",
                flush=True,
            )

    metrics = pd.DataFrame(metric_rows)
    events = pd.DataFrame(event_rows)
    summary = pd.DataFrame(
        [candidate_summary(group) for _, group in metrics.groupby("continuation_candidate_id")]
    ).sort_values(
        ["target_class", "accepted_after_spatial_growth", "K_plateau_high_minus_low"],
        ascending=[True, False, False],
    )
    metrics.to_csv(out / "spatial_promotion_metrics.csv", index=False)
    events.to_csv(out / "spatial_promotion_events.csv", index=False)
    summary.to_csv(out / "spatial_promotion_summary.csv", index=False)
    accepted = summary[summary.accepted_after_spatial_growth.astype(bool)]
    accepted.to_csv(out / "spatial_promotion_accepted.csv", index=False)
    report = {
        "n_candidates": int(summary.shape[0]),
        "n_temperature_runs": int(metrics.shape[0]),
        "n_accepted": int(accepted.shape[0]),
        "accepted_by_class": accepted.target_class.value_counts().to_dict(),
        "target_extension_um": float(a.target_extension_um),
        "status": "V9_9_SPATIAL_PROMOTION_COMPLETE_NOT_2D_VALIDATED",
        "output": str(out),
    }
    (out / "spatial_promotion_report.json").write_text(json.dumps(report, indent=2))
    config = vars(a).copy()
    for key, value in list(config.items()):
        if isinstance(value, Path):
            config[key] = str(value)
    config.update(
        {
            "active_state_for_this_workflow": "arrhenius_fracture.moving_process_zone_v99",
            "global_package_state_unchanged": True,
            "independent_entropy_slots": {
                "pt_peierls_entropy_ratio": "activation_entropy_kB",
                "pt_taylor_entropy_ratio": "activation_entropy_kB",
            },
        }
    )
    (out / "spatial_promotion_config.json").write_text(json.dumps(config, indent=2))
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
