#!/usr/bin/env python3
"""Physical-interpretation audit of the eight DBTT Taylor/Peierls PF cases.

This program is analysis-only.  It reads immutable accepted trajectories and
default-off observer records, evaluates the committed production source laws
on frozen states, and writes a self-contained paper/audit/option-bank bundle.
It never calls the stochastic trajectory driver.
"""
from __future__ import annotations

import hashlib
from itertools import combinations, product
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys
from types import SimpleNamespace
from typing import Any, Iterable

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from scipy.special import gammainc


ROOT = Path(__file__).resolve().parents[1]
PF_REPO = Path("/private/tmp/pf-taylor-peierls-transfer")
PF_RUN_ROOT = Path("/private/tmp/oneD-v2-taylor-peierls-spatial-pf-runs")
PF_MATRIX_MANIFEST = PF_RUN_ROOT / "spatial_pf_matrix_manifest.json"
BASE = ROOT / "analysis_outputs" / "oneD_v2_taylor_peierls_spatial_transfer"
OUT = ROOT / "analysis_outputs" / "taylor_peierls_spatial_coupling_paper_audit"
FIG = OUT / "figures"
FIGDATA = OUT / "figure_source_data"
KERNEL_PATH = Path(
    "/Volumes/Data/Data/Nanopillar_calculation/"
    "PF-fracture-fatigue_v10_2_21_persistent_sites_top1/runs/"
    "v10_2_28_kernel_cache/1447653d199f0b43cb475951092d69444c9b785f6fdf518c723792abb3b1f5e5/family.json"
)
TEMPERATURE_K = 1100.0
HAZARD_SEED = 1008666
TARGET_UM = 300.0
CHECKPOINTS_UM = (0.0, 25.0, 50.0, 100.0, 200.0, 300.0)
TP_FIELDS = (
    "peierls_H0_eV", "peierls_activation_entropy_kB", "peierls_exp_a", "peierls_exp_n",
    "taylor_H0_eV", "taylor_activation_entropy_kB", "taylor_exp_a", "taylor_exp_n",
    "taylor_corr_rho_c_m2", "taylor_corr_scale",
)
META_FIELDS = ("option_key", "candidate_id", "material_class", "role", "mechanism_summary",
               "validation_status", "target_class", "L_pz_um_recommended", "n_bins_recommended")
G_PA = 410.0e9 / (2.0 * (1.0 + 0.28))
B_M = 2.74e-10
NU = 0.28
R0_M = 1.0e-6
BLUNTING_LENGTH_M = 0.5e-6
DX_M = 50.0e-6 / 80.0
MULTIHIT_M = 3.0
MULTIHIT_TAU_S = 1.0e-6
KB_EV_PER_K = 8.617333262145e-5
COUNTERFACTUAL = "COUNTERFACTUAL_DIAGNOSTIC_NOT_PRODUCTION_PHYSICS"

EXPECTED = {
    "v913_zeroD_sobol_0202500": "CONTROL",
    "oneD_v2_dbtt_TP_f07fe99faea93479": "TRANSPORT_EXTREME",
    "oneD_v2_dbtt_TP_4895f9e5b44deea5": "RETENTION_EXTREME",
    "oneD_v2_dbtt_TP_7e668ee637fc3ac5": "NEAR_TIP_RETENTION_EXTREME",
    "oneD_v2_dbtt_TP_b60d78111740b058": "PERSISTENT_WAKE_EXTREME",
    "oneD_v2_dbtt_TP_bd5be1610f6e1bce": "MOBILE_TRANSPORT_EXTREME",
    "oneD_v2_dbtt_TP_9555ff54d637c974": "BACKSTRESS_SHIELDING_EXTREME",
    "oneD_v2_dbtt_TP_f2817e7998cb7be6": "BALANCED_MAXIMIN_MEDOID",
}
SHORT_LABEL = {
    "CONTROL": "Control", "TRANSPORT_EXTREME": "Transport",
    "RETENTION_EXTREME": "Retention", "NEAR_TIP_RETENTION_EXTREME": "Near-tip",
    "PERSISTENT_WAKE_EXTREME": "Wake", "MOBILE_TRANSPORT_EXTREME": "Mobile",
    "BACKSTRESS_SHIELDING_EXTREME": "Backstress", "BALANCED_MAXIMIN_MEDOID": "Balanced",
}
PALETTE = {
    role: color for role, color in zip(EXPECTED.values(), [
        "#0072B2", "#E69F00", "#009E73", "#D55E00",
        "#CC79A7", "#56B4E9", "#F0E442", "#000000",
    ])
}

sys.path.insert(0, str(PF_REPO))
from arrhenius_fracture.campaign_calibrated_tip import _campaign_backstress  # noqa: E402
from arrhenius_fracture.fractional_moving_frame import _translate_toward_tip  # noqa: E402
from arrhenius_fracture.material_manifest import MaterialManifest  # noqa: E402
from arrhenius_fracture.persistent_site_source_v10221 import (  # noqa: E402
    effective_front_width_m, persistent_site_multiplicity,
)
from arrhenius_fracture.signed_kernel_family_v10214 import (  # noqa: E402
    ActiveOnlySigned2DShieldingKernelFamily,
)
from arrhenius_fracture.unified_mpz import MPZConfig, UnifiedMPZState  # noqa: E402


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def tree_hash(paths: Iterable[Path]) -> str:
    digest = hashlib.sha256()
    for path in sorted(paths, key=lambda p: p.name):
        digest.update(path.name.encode())
        digest.update(path.read_bytes())
    return digest.hexdigest()


def git_head(path: Path) -> str:
    return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=path, text=True).strip()


def canonical_hash(mapping: dict[str, Any]) -> str:
    raw = json.dumps(mapping, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    return hashlib.sha256(raw).hexdigest()


def builtin_scalar(value: Any) -> Any:
    """Return a JSON-safe scalar without losing numeric type or precision."""
    return value.item() if isinstance(value, np.generic) else value


def load_inputs() -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    selection = pd.read_csv(BASE / "oneD_v2_spatial_transfer_selected_candidates.csv")
    registry = pd.read_csv(BASE / "oneD_v2_spatial_transfer_pf_registry.csv")
    matrix = json.loads(PF_MATRIX_MANIFEST.read_text())
    observed = dict(zip(selection.candidate_id, selection.selection_role))
    if observed != EXPECTED:
        raise RuntimeError("authoritative candidate IDs/roles do not match the mission contract")
    if set(registry.candidate_id) != set(EXPECTED):
        raise RuntimeError("PF registry does not contain exactly the eight authoritative cases")
    if len(matrix["cases"]) != 8 or matrix["maximum_concurrent_heavy_workers"] != 2:
        raise RuntimeError("PF execution manifest violates the bounded eight-case contract")
    return selection, registry, matrix


def steps_path(case_dir: Path) -> Path:
    path = case_dir / "steps_1100K.csv"
    if not path.is_file():
        raise FileNotFoundError(path)
    return path


def load_case(case: dict[str, Any]) -> dict[str, Any]:
    root = Path(case["case_path"])
    sp = steps_path(root)
    ap = root / "anisotropic_emission_audit_v10174.json"
    gp = root / "stochastic_avalanche_geometry_events.json"
    cp = root / "crack_path_1100K.csv"
    mp = root / "selected_material_manifest_v10_2_22.csv"
    args = root / "run_args.json"
    steps = pd.read_csv(sp)
    audit_payload = json.loads(ap.read_text())
    audit = audit_payload["records"]
    if len(steps) != len(audit):
        raise RuntimeError(f"step/observer length mismatch for {case['candidate_id']}")
    reached = np.flatnonzero(steps.crack_extension_m.to_numpy(float) * 1e6 >= TARGET_UM - 1e-9)
    if not len(reached):
        raise RuntimeError(f"case did not reach {TARGET_UM} um: {case['candidate_id']}")
    stop = int(reached[0])
    return {
        "meta": case, "root": root, "steps_path": sp, "audit_path": ap,
        "geometry_path": gp, "crack_path": cp, "material_path": mp, "args_path": args,
        "steps": steps.iloc[: stop + 1].copy(), "audit": audit[: stop + 1],
        "audit_payload": audit_payload, "material": MaterialManifest.from_csv(mp),
        "material_row": pd.read_csv(mp).iloc[0], "args": json.loads(args.read_text()),
    }


def event_positions(steps: pd.DataFrame) -> np.ndarray:
    return np.flatnonzero(steps.n_fire.to_numpy(float) > 0.0)


def physical_onset_positions(steps: pd.DataFrame) -> list[int]:
    events = event_positions(steps)
    positions = [int(events[0])]
    for i in range(1, len(events)):
        between = steps.iloc[int(events[i - 1]) + 1 : int(events[i]) + 1]
        if (between.adaptive_frac.to_numpy(float) >= 1.0 - 1e-12).any():
            positions.append(int(events[i]))
    return positions


def array(record: dict[str, Any], name: str) -> np.ndarray:
    return np.asarray(record[name], dtype=float)


def state_arrays(record: dict[str, Any]) -> dict[str, np.ndarray]:
    return {
        "mobile": array(record, "mobile_active_by_system_bin"),
        "retained": array(record, "retained_active_by_system_bin"),
        "mobile_positive": array(record, "mobile_positive_by_system_bin"),
        "mobile_negative": array(record, "mobile_negative_by_system_bin"),
        "retained_positive": array(record, "retained_positive_by_system_bin"),
        "retained_negative": array(record, "retained_negative_by_system_bin"),
        "wake_mobile": array(record, "mobile_wake_by_system_bin"),
        "wake_retained": array(record, "retained_wake_by_system_bin"),
        "wake_mobile_positive": array(record, "wake_mobile_positive_by_system_bin"),
        "wake_mobile_negative": array(record, "wake_mobile_negative_by_system_bin"),
        "wake_retained_positive": array(record, "wake_retained_positive_by_system_bin"),
        "wake_retained_negative": array(record, "wake_retained_negative_by_system_bin"),
        "x": array(record, "active_x_ahead_of_tip_m"),
        "wake_x": array(record, "wake_x_behind_tip_m"),
    }


def moments(x: np.ndarray, weight: np.ndarray) -> tuple[float, float]:
    w = np.maximum(np.asarray(weight, dtype=float), 0.0)
    total = float(w.sum())
    if total <= 0.0:
        return 0.0, 0.0
    center = float(np.sum(x * w) / total)
    width = float(np.sqrt(max(np.sum((x - center) ** 2 * w) / total, 0.0)))
    return center, width


def make_backstress_state(record: dict[str, Any]) -> SimpleNamespace:
    fields = state_arrays(record)
    return SimpleNamespace(
        mobile=fields["mobile"], retained=fields["retained"], x=fields["x"],
        dx=float(fields["x"][1] - fields["x"][0]),
        cfg=SimpleNamespace(blunting_length_m=BLUNTING_LENGTH_M,
                            taylor_stress_fraction=1.0 / math.sqrt(3.0)),
        _campaign_b=B_M, _campaign_G_Pa=G_PA, _campaign_backstress_scale=1.0,
    )


def exact_backstress(record: dict[str, Any]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    return _campaign_backstress(make_backstress_state(record))


def source_geometry(record: dict[str, Any], material: MaterialManifest) -> dict[str, float]:
    # These four diagnostics were evaluated by the production engine itself.
    # Preserve them rather than silently rebuilding with defaults: in the
    # neutral archive they are the source-exact geometry at the accepted row.
    # (The accumulated-slip profile needed to reconstruct radius is not saved.)
    width = float(record["persistent_site_front_width_m"])
    rho_width = float(record.get(
        "persistent_site_width_density_m2",
        5.0e12 * (10.0e-6 / max(width, 1.0e-30)) ** 2,
    ))
    radius = float(record["persistent_tip_radius_m"])
    multiplicity = float(record["persistent_site_multiplicity_per_system"])
    return {"radius_m": radius, "front_width_m": width,
            "rho_width_m2": rho_width, "multiplicity": multiplicity}


class KernelProbe:
    def __init__(self) -> None:
        family = ActiveOnlySigned2DShieldingKernelFamily.from_json(KERNEL_PATH)
        dummy = SimpleNamespace(
            n_systems=2,
            x=(np.arange(80, dtype=float) + 0.5) * (50.0e-6 / 80.0),
            wake_x=(np.arange(160, dtype=float) + 0.5) * (100.0e-6 / 160.0),
        )
        self.family = family.bind_to_state_grid(dummy)

    def active_kernel(self, extension_m: float) -> np.ndarray:
        kernel = self.family.resolve(
            r_eff_over_r0=1.0, opening_strength_fraction=0.0,
            crack_extension_m=float(extension_m),
        )[0]
        return np.asarray(kernel, dtype=float)

    def shielding(self, record: dict[str, Any], extension_m: float,
                  support: str = "full_active") -> float:
        fields = state_arrays(record)
        signed = fields["retained_positive"] - fields["retained_negative"]
        if support == "near_tip":
            signed = signed.copy(); signed[:, fields["x"] > 2.0e-6] = 0.0
        elif support not in {"full_active", "full_production_with_wake"}:
            raise ValueError(support)
        # The production family has an identically zero wake operator.
        return float(np.sum(self.active_kernel(extension_m) * signed))


def cleavage_rate(material: MaterialManifest, stress_Pa: float, T_K: float) -> tuple[float, float, float]:
    barrier = float(np.asarray(material.cleavage.values_eV(max(stress_Pa, 0.0), T_K)))
    raw = float(np.asarray(material.cleavage.rate(max(stress_Pa, 0.0), T_K)))
    effective = float(gammainc(MULTIHIT_M, min(max(raw, 0.0) * MULTIHIT_TAU_S, 1.0e12)) / MULTIHIT_TAU_S)
    return barrier, raw, effective


def transport_probe(record: dict[str, Any], material: MaterialManifest,
                    emission_stress_by_system: np.ndarray, T_K: float) -> dict[str, np.ndarray]:
    cfg = MPZConfig(length_m=50.0e-6, n_bins=80, n_systems=2, source_bin_count=2,
                    blunting_length_m=BLUNTING_LENGTH_M, forest_density_floor_m2=5.0e12,
                    peierls_stress_fraction=1.0 / math.sqrt(3.0),
                    taylor_stress_fraction=1.0 / math.sqrt(3.0),
                    wake_length_m=100.0e-6, wake_n_bins=160, wake_shielding=False)
    state = UnifiedMPZState(material, cfg)
    fields = state_arrays(record)
    state.mobile = fields["mobile"].copy(); state.retained = fields["retained"].copy()
    rho_shared = state.local_forest_density_m2(False)
    by_system = []
    for system in range(2):
        stress_profile = state.local_stress_profile_Pa(float(emission_stress_by_system[system]))
        by_system.append(state._transport_rates(stress_profile, rho_shared, T_K, B_M))
    return {name: np.asarray([item[name] for item in by_system]) for name in by_system[0]}


def frozen_probe(record: dict[str, Any], material: MaterialManifest, kernel: KernelProbe,
                 *, extension_m: float, native_K_Pa_sqrt_m: float,
                 radius_m: float | None = None, shielding_mode: str = "actual",
                 backstress_mode: str = "actual", multiplicity: float | None = None,
                 state_support: str = "full_active", T_K: float = TEMPERATURE_K) -> dict[str, Any]:
    geom = source_geometry(record, material)
    radius = geom["radius_m"] if radius_m is None else float(radius_m)
    shield = kernel.shielding(record, extension_m, state_support)
    if shielding_mode == "zero": shield = 0.0
    rho, tau_back, recomputed_sigma_back = exact_backstress(record)
    sigma_back = np.asarray(record["anisotropic_sigma_back_by_system_Pa"], dtype=float)
    if backstress_mode == "zero": sigma_back = np.zeros_like(sigma_back)
    mult = geom["multiplicity"] if multiplicity is None else float(multiplicity)
    factors = array(record, "anisotropic_drive_factors")
    opening = max(float(native_K_Pa_sqrt_m), 0.0) / math.sqrt(2.0 * math.pi * radius)
    cleavage_stress = max(float(native_K_Pa_sqrt_m) - shield, 0.0) / math.sqrt(2.0 * math.pi * radius)
    emission_stress = np.maximum(factors * opening - sigma_back, 0.0)
    emission_barrier = np.asarray(material.emission.values_eV(emission_stress, T_K), dtype=float)
    emission_rate = np.asarray(material.emission.rate(emission_stress, T_K), dtype=float)
    aggregate_emission = float(mult * np.sum(emission_rate))
    cleave_barrier, cleave_raw, cleave_effective = cleavage_rate(material, cleavage_stress, T_K)
    transport = transport_probe(record, material, emission_stress, T_K)
    return {
        "native_K_Pa_sqrt_m": float(native_K_Pa_sqrt_m), "extension_m": float(extension_m),
        "radius_m": radius, "front_width_m": geom["front_width_m"],
        "shielding_Pa_sqrt_m": shield, "backstress_mean_Pa": float(np.mean(sigma_back)),
        "accepted_snapshot_recomputed_backstress_mean_Pa": float(np.mean(recomputed_sigma_back)),
        "multiplicity_per_system": mult, "opening_stress_Pa": opening,
        "cleavage_stress_Pa": cleavage_stress,
        "cleavage_barrier_eV": cleave_barrier, "cleavage_raw_rate_s": cleave_raw,
        "cleavage_effective_rate_s": cleave_effective,
        "emission_stress_by_system_Pa_json": json.dumps(emission_stress.tolist()),
        "emission_barrier_by_system_eV_json": json.dumps(emission_barrier.tolist()),
        "emission_rate_by_system_s_json": json.dumps(emission_rate.tolist()),
        "aggregate_emission_rate_s": aggregate_emission,
        "log10_cleavage_to_emission": math.log10(max(cleave_effective, 1e-300) / max(aggregate_emission, 1e-300)),
        "peierls_rate_max_s": float(np.max(transport["peierls"])),
        "encounter_rate_max_s": float(np.max(transport["encounter"])),
        "taylor_completion_rate_max_s": float(np.max(transport["taylor"])),
        "state_support": state_support,
    }


def normalized_table_hash(path: Path) -> str:
    frame = pd.read_csv(path)
    numeric = frame.select_dtypes(include=[np.number]).copy()
    return canonical_hash({c: numeric[c].astype(float).round(15).tolist() for c in numeric})


def build_provenance(selection: pd.DataFrame, registry: pd.DataFrame,
                     matrix: dict[str, Any], cases: list[dict[str, Any]]) -> pd.DataFrame:
    rows = []
    prior_onsets = pd.read_csv(BASE / "pf_2d_spatial_transfer_onsets.csv")
    prior_avalanches = pd.read_csv(BASE / "pf_2d_spatial_transfer_physical_avalanches.csv")
    registry_index = registry.set_index("candidate_id")
    selection_index = selection.set_index("candidate_id")
    physical_columns = [c for c in registry.columns if c not in META_FIELDS]
    non_tp = [c for c in physical_columns if c not in TP_FIELDS]
    fixed_non_tp_hashes = {
        canonical_hash({c: registry_index.loc[candidate, c].item() if isinstance(registry_index.loc[candidate, c], np.generic)
                        else registry_index.loc[candidate, c] for c in non_tp})
        for candidate in registry_index.index
    }
    if len(fixed_non_tp_hashes) != 1:
        raise RuntimeError("a non-Taylor/non-Peierls material coordinate differs")
    for data in cases:
        case = data["meta"]; cid = case["candidate_id"]
        records = data["audit"]
        args = data["args"]
        observer_neutral = all(
            r.get("taylor_peierls_state_profile_feedback") is False
            for r in records if r.get("taylor_peierls_state_profile_schema") is not None
        )
        no_clip = all(not bool(r.get("anisotropic_factors_clipped", False)) for r in records)
        stack_path = data["root"] / "v10_2_17_final_signed_stochastic_stack.json"
        stack = json.loads(stack_path.read_text())
        family = stack["state_resolved_kernel_family"]
        no_extrap = bool(family["interpolation"].get("extrapolation_allowed") is False)
        mesh_contract = {
            key: args[key] for key in (
                "nx", "ny", "tip_h_fine", "tip_ratio", "mesh_levels", "crack_backend",
                "wake_length_um", "wake_n_bins", "wake_shielding", "mpz_n_bins",
                "mpz_length_um", "mpz_blunting_length_um", "mpz_source_bins",
            )
        }
        loading_contract = {key: args[key] for key in (
            "Kdot", "dU", "dt", "adaptive_events", "adaptive_event_target",
            "adaptive_min_frac", "adaptive_safety", "target_crack_extension_um",
            "temperatures", "crystal_theta_deg", "bulk_kinetics_model",
        )}
        event_count = int((data["steps"].n_fire > 0).sum())
        local_onsets = prior_onsets[prior_onsets.candidate_id.eq(cid)]
        local_avalanches = prior_avalanches[prior_avalanches.candidate_id.eq(cid)]
        geometry_files = [data["geometry_path"], data["crack_path"],
                          data["root"] / "sharp_wake_advance_log.csv",
                          data["root"] / "fronts_1100K.csv"]
        row = selection_index.loc[cid]
        rows.append({
            "candidate_id": cid, "original_1d_selection_role": case["selection_role"],
            "full_material_sha256": row.full_material_sha256,
            "taylor_peierls_subvector_sha256": row.taylor_peierls_subvector_sha256,
            "cleavage_barrier_sha256": row.cleavage_barrier_sha256,
            "emission_barrier_sha256": row.emission_barrier_sha256,
            "fixed_non_taylor_peierls_coordinates_sha256": next(iter(fixed_non_tp_hashes)),
            "source_repository": "ukaiiaku-maker/PF-fracture-fatigue",
            "source_commit": matrix["pf_runner_commit"], "run_directory": str(data["root"]),
            "steps_file": str(data["steps_path"]), "steps_sha256": sha(data["steps_path"]),
            "observer_file": str(data["audit_path"]), "observer_sha256": sha(data["audit_path"]),
            "geometry_events_sha256": sha(data["geometry_path"]),
            "crack_path_numeric_sha256": normalized_table_hash(data["crack_path"]),
            "geometry_path_bundle_sha256": tree_hash(geometry_files),
            "mesh_wake_fingerprint": canonical_hash(mesh_contract),
            "loading_contract_sha256": canonical_hash(loading_contract),
            "kernel_family_sha256": sha(KERNEL_PATH),
            "event_transaction_count": event_count,
            "physical_avalanche_count": int(len(local_avalanches)),
            "physical_onset_step_indices_json": json.dumps(local_onsets.pre_event_step.astype(int).tolist()),
            "extension_checkpoints_um_json": json.dumps(CHECKPOINTS_UM),
            "target_right_censored": bool(local_avalanches.target_right_censored.iloc[-1]),
            "observer_neutral_no_feedback": observer_neutral,
            "anisotropic_factor_clipping_absent": no_clip,
            "kernel_extrapolation_absent": no_extrap,
            "geometry_files_complete": all(path.is_file() and path.stat().st_size > 0 for path in geometry_files),
            "production_trajectory_rerun_for_audit": False,
            "fatigue_evaluated": False,
        })
    result = pd.DataFrame(rows)
    for field in ("cleavage_barrier_sha256", "emission_barrier_sha256",
                  "fixed_non_taylor_peierls_coordinates_sha256", "mesh_wake_fingerprint",
                  "loading_contract_sha256", "crack_path_numeric_sha256"):
        if result[field].nunique() != 1:
            raise RuntimeError(f"shared contract field differs: {field}")
    result.to_csv(OUT / "PF_TP_CASE_PROVENANCE_TABLE.csv", index=False)
    audit_manifest = {
        "schema": "PF_TP_SPATIAL_COUPLING_AUDIT_MANIFEST_v1",
        "immutable_dataset": True, "material_class": "DBTT", "temperature_K": TEMPERATURE_K,
        "hazard_seed": HAZARD_SEED, "theta_deg": 0.0, "bulk_mode": "tip_only",
        "target_projected_extension_um": TARGET_UM, "candidate_count": 8,
        "result_source_commit": "ff9e12f7c4c00e1171c1dead75f29920b089f935",
        "pf_source_commit": matrix["pf_runner_commit"],
        "new_stochastic_trajectories": 0, "new_femczm_runs": 0,
        "memory_term_added": False, "fatigue_evaluated": False,
        "observer_contract": "DEFAULT_OFF_NEUTRAL_NO_FEEDBACK",
        "cases": result.to_dict("records"),
    }
    (OUT / "PF_TP_SPATIAL_COUPLING_AUDIT_MANIFEST.json").write_text(
        json.dumps(audit_manifest, indent=2, sort_keys=True) + "\n"
    )
    return result


def build_dependency_audit() -> pd.DataFrame:
    rows = [
        ("mobile population", "campaign_calibrated_tip.py", "_campaign_local_density_m2", "active mobile[system,bin]", "tip-relative x>0", "all active bins; exp(-x/L), L=max(0.5 um, dx)=0.625 um", "line count", "unsigned", "translated toward tip; crossing enters wake", False),
        ("retained/tangled population", "campaign_calibrated_tip.py", "_campaign_local_density_m2", "active retained[system,bin]", "tip-relative x>0", "all active bins; exp(-x/L), L=max(0.5 um, dx)=0.625 um", "line count", "unsigned for backstress", "translated; crossing enters wake", False),
        ("near-tip retained content", "campaign_calibrated_tip.py", "_campaign_backstress", "mobile+retained active", "tip-relative", "exponential; first 0.625-2 um dominate", "m^-2 after volume conversion", "nonnegative", "recentered continuously", False),
        ("source multiplicity", "persistent_site_physical_width_v10222.py", "physical_source_geometry/persistent_site_multiplicity", "radius, density-limited along-front width, rho_site0", "tip-local scalar", "near-tip density through mesh-independent width; no wake", "sites/system", "positive", "recomputed continuously", False),
        ("active accumulated slip", "unified_mpz.py", "accumulated_slip", "accepted inner-step slip increments", "tip-relative", "all active bins", "line count", "nonnegative", "translated with active populations; profile not serialized by neutral observer", False),
        ("tip radius", "unified_mpz.py", "local_slip_count/blunted_radius", "active accumulated slip", "tip-relative", "all active bins; exp(-x/L), L=max(0.5 um, dx)=0.625 um", "m", "positive blunting", "slip convects into wake, natural resharpening", False),
        ("front width", "persistent_site_physical_width_v10222.py", "physical_source_geometry/effective_front_width_m", "active unsigned near-tip density", "orthogonal along-front scalar", "rho^-1/2; bounded by explicit physical minimum/b and 50 um; ahead-tip dx is deliberately not a floor", "m", "positive", "recomputed after translation", False),
        ("backstress", "campaign_calibrated_tip.py", "_campaign_backstress", "active unsigned mobile+retained", "tip-relative", "all active bins with exp[-x/max(0.5um,dx)]; L=0.625um", "Pa", "positive, subtracts from emission drive", "recomputed after every state update", False),
        ("signed shielding", "signed_burgers_shared_v1025.py", "_active_K", "active signed retained; mobile fraction=0", "tip-relative physical-x atlas", "all 80 active bins, extension-dependent measured kernel", "Pa sqrt(m)", "signed; K_cleave=K_native-K_shield", "active population translates; kernel follows crack extension", False),
        ("retained wake shielding", "state_resolved_signed_engine_v10214.py", "K_shield/_wake_K", "wake signed retained", "tip-relative wake", "identically zero production wake operator", "Pa sqrt(m)", "zero", "wake persists/advects but no shielding feedback", False),
        ("opening tensor", "anisotropic_emission_v10174.py", "probe_tensor_ahead", "2-D FEM stress tensor", "laboratory mesh, sector ahead of tip", "0.25-1.75 times 10um probe, adaptive expansion", "Pa", "standard tension/shear", "reprobed after geometry solve", False),
        ("channel tensors/resolved shears", "anisotropic_emission_v10174.py", "resolved_channel_drives", "opening tensor and two reduced slip traces", "laboratory tensor projected into tip frame", "same finite-radius sector", "Pa", "signed shear chooses Burgers species", "reprobed after geometry solve", False),
        ("Peierls transport", "unified_mpz.py", "_transport_rates", "local emission stress profile and forest density", "tip-relative active bins", "all active bins", "s^-1 and m/s", "forward magnitude", "advects mobile away from tip", False),
        ("encounter retention", "unified_mpz.py", "_transport_rates/_exchange", "Peierls velocity and local retained forest", "tip-relative active bins", "all active bins", "s^-1", "mobile to retained", "local exchange before transport", False),
        ("Taylor completion/release", "unified_mpz.py", "_transport_rates/_exchange", "Taylor barrier, stress, density multiplicity", "tip-relative active bins", "all active bins", "s^-1", "retained to mobile", "local exchange before transport", False),
        ("cleavage barrier/rate", "unified_front.py", "lambda_cleave", "K_native, active shielding, radius", "tip-local scalar", "no direct dislocation bins beyond radius/shielding", "eV, s^-1", "shielding subtracts K", "integrated continuously; first passage B", False),
        ("emission barrier/rate", "persistent_site_source_v10221.py", "_persistent_emit", "unshielded opening, tensor factor, backstress, multiplicity", "tip-local scalar per channel", "near-tip active state through backstress/width/radius", "eV, s^-1", "backstress subtracts drive", "persistent sites not consumed/refreshed", False),
        ("cleavage first-passage action", "kinetic_tip_cell.py", "_integrate_coupled", "multihit cleavage rate", "scalar front state", "time integral; max action substep 0.02", "dimensionless", "increasing", "B subtracts one at checkpoint; MPZ not reset", False),
        ("PF structural solve", "sharp_front.py", "sharp_wake 2-D solve", "geometry, opening, elasticity, damage wake", "laboratory frame", "whole specimen", "N/m, J/m2, Pa sqrt(m)", "root-signed J", "geometry remesh/advance", False),
    ]
    columns = ["quantity", "source_file", "class_or_function", "input_state", "coordinate_frame",
               "spatial_support_or_weighting", "units", "sign_convention", "event_translation_behavior",
               "wake_behind_current_tip_contributes"]
    table = pd.DataFrame(rows, columns=columns)
    table.to_csv(OUT / "pf_tp_state_support_table.csv", index=False)
    graph = {
        "schema": "pf_tp_state_to_tip_dependency_graph_v1",
        "nodes": sorted(set(table.quantity) | {"native PF KJ", "fracture event", "wake ledger"}),
        "edges": [
            {"from": "mobile population", "to": "backstress", "coupling": "unsigned exponential near-tip density"},
            {"from": "retained/tangled population", "to": "backstress", "coupling": "unsigned exponential near-tip density"},
            {"from": "retained/tangled population", "to": "signed shielding", "coupling": "signed active measured kernel"},
            {"from": "retained wake shielding", "to": "signed shielding", "coupling": "disabled: zero kernel"},
            {"from": "active accumulated slip", "to": "tip radius", "coupling": "exponential near-tip weight"},
            {"from": "backstress", "to": "emission barrier/rate", "coupling": "subtracted from tensor-weighted opening"},
            {"from": "tip radius", "to": "emission barrier/rate", "coupling": "K/sqrt(2*pi*r)"},
            {"from": "tip radius", "to": "cleavage barrier/rate", "coupling": "(K-Kshield)/sqrt(2*pi*r)"},
            {"from": "signed shielding", "to": "cleavage barrier/rate", "coupling": "K_native-Kshield"},
            {"from": "cleavage barrier/rate", "to": "cleavage first-passage action", "coupling": "multihit rate integration"},
            {"from": "cleavage first-passage action", "to": "fracture event", "coupling": "B reaches one"},
            {"from": "fracture event", "to": "wake ledger", "coupling": "continuous moving-frame translation"},
        ],
        "key_findings": {
            "wake_feedback": "retained wake is persisted but wake shielding is identically zero; it does not enter backstress/radius/multiplicity",
            "radius_structural_coupling": "tip radius is a local hazard/source scalar and is not passed into the sharp-wake structural PF solve or kernel interpolation",
            "commit_reset": "outer fracture commit resets/subtracts first-passage action only; persistent sites and MPZ populations are not renewed",
        },
    }
    (OUT / "pf_tp_state_to_tip_dependency_graph.json").write_text(json.dumps(graph, indent=2) + "\n")
    return table


def state_totals(record: dict[str, Any]) -> dict[str, float]:
    f = state_arrays(record)
    return {
        "active_mobile": float(f["mobile"].sum()),
        "active_retained": float(f["retained"].sum()),
        "wake_mobile": float(f["wake_mobile"].sum()),
        "wake_retained": float(f["wake_retained"].sum()),
        "mobile_total": float(f["mobile"].sum() + f["wake_mobile"].sum()),
        "retained_total": float(f["retained"].sum() + f["wake_retained"].sum()),
        "all_total": float(f["mobile"].sum() + f["wake_mobile"].sum() +
                           f["retained"].sum() + f["wake_retained"].sum()),
    }


def regional_state(record: dict[str, Any]) -> dict[str, float]:
    f = state_arrays(record)
    rows: dict[str, float] = {}
    active_mobile = f["mobile"].sum(axis=0); active_retained = f["retained"].sum(axis=0)
    wake_mobile = f["wake_mobile"].sum(axis=0); wake_retained = f["wake_retained"].sum(axis=0)
    wake_distance = np.abs(f["wake_x"])
    for name, lo, hi in (("wake_0_10um", 0.0, 10e-6), ("wake_10_50um", 10e-6, 50e-6),
                         ("wake_50_100um", 50e-6, 100e-6), ("wake_beyond_100um", 100e-6, math.inf)):
        mask = (wake_distance >= lo) & (wake_distance < hi)
        rows[f"{name}_mobile"] = float(wake_mobile[mask].sum())
        rows[f"{name}_retained"] = float(wake_retained[mask].sum())
    for name, limit in (("near_tip_0_2um", 2e-6), ("ahead_0_10um", 10e-6),
                        ("ahead_full_active", math.inf)):
        mask = f["x"] < limit
        rows[f"{name}_mobile"] = float(active_mobile[mask].sum())
        rows[f"{name}_retained"] = float(active_retained[mask].sum())
    return rows


def expanded_state_rows(common: dict[str, Any], record: dict[str, Any], tip_m: float,
                        boundary: str) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    f = state_arrays(record); tip_rows, lab_rows = [], []
    for region, x, mobile, retained in (
        ("ACTIVE_AHEAD", f["x"], f["mobile"], f["retained"]),
        ("WAKE_BEHIND", f["wake_x"], f["wake_mobile"], f["wake_retained"]),
    ):
        for system in range(mobile.shape[0]):
            for i, distance in enumerate(x):
                base = {**common, "event_boundary": boundary, "region": region,
                        "system_index": system, "bin_index": i,
                        "mobile_count": float(mobile[system, i]),
                        "retained_count": float(retained[system, i])}
                # The neutral observer stores active x as positive ahead of the
                # tip and wake x as negative behind it.
                relative = float(distance)
                tip_rows.append({**base, "tip_relative_position_m": relative})
                lab_rows.append({**base, "laboratory_position_m": tip_m + relative,
                                 "current_tip_laboratory_m": tip_m})
    return tip_rows, lab_rows


def translated_prior_ledger(before: dict[str, Any], distance_m: float) -> dict[str, float]:
    f = state_arrays(before); result: dict[str, float] = {}
    for family in ("mobile", "retained"):
        active = f[family]
        wake = f[f"wake_{family}"]
        shifted_wake, old_lost = UnifiedMPZState._advect_forward(
            wake, distance_m, abs(float(f["wake_x"][1] - f["wake_x"][0]))
        )
        shifted_active, crossed, active_lost = _translate_toward_tip(
            active, distance_m, float(f["x"][1] - f["x"][0]),
            len(f["wake_x"]), abs(float(f["wake_x"][1] - f["wake_x"][0])),
        )
        result[f"operator_active_{family}_after"] = float(shifted_active.sum())
        result[f"operator_wake_{family}_after"] = float(shifted_wake.sum() + crossed.sum())
        result[f"operator_{family}_entering_wake"] = float(crossed.sum())
        result[f"operator_{family}_discarded"] = float(old_lost + active_lost)
        result[f"operator_{family}_conservation_residual"] = float(
            active.sum() + wake.sum() - shifted_active.sum() - shifted_wake.sum() -
            crossed.sum() - old_lost - active_lost
        )
    return result


def observed_role(record: dict[str, Any], original_role: str) -> str:
    totals = state_totals(record)
    fraction = totals["retained_total"] / max(totals["all_total"], 1e-300)
    if original_role == "CONTROL":
        return "MOBILE_RICH_REFERENCE_AT_REINITIATION"
    if fraction >= 0.9:
        return "RETENTION_DOMINATED_LONG_WAKE_AT_REINITIATION"
    return "MOBILE_RICH_LONG_WAKE_AT_REINITIATION"


def build_translation_audit(cases: list[dict[str, Any]]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    event_rows: list[dict[str, Any]] = []
    tip_rows: list[dict[str, Any]] = []
    lab_rows: list[dict[str, Any]] = []
    role_rows: list[dict[str, Any]] = []
    for data in cases:
        cid = data["meta"]["candidate_id"]; role = data["meta"]["selection_role"]
        steps, audit = data["steps"], data["audit"]
        events = event_positions(steps); onsets = physical_onset_positions(steps)
        reinit_record = audit[onsets[1]]
        role_rows.append({"candidate_id": cid, "original_1d_selection_role": role,
                          "pf_observed_microstructure_role": observed_role(reinit_record, role)})
        prior_cumulative = {name: 0.0 for name in ("emitted", "escaped", "recovered")}
        for transaction, pos in enumerate(events):
            pos = int(pos); before_pos = max(pos - 1, 0)
            before, after = audit[before_pos], audit[pos]
            before_totals, after_totals = state_totals(before), state_totals(after)
            extension_before = float(steps.iloc[before_pos].crack_extension_m)
            extension_after = float(steps.iloc[pos].crack_extension_m)
            distance = max(extension_after - extension_before, 0.0)
            cumulative = {
                "emitted": float(after.get("developed_state_cumulative_emitted", 0.0)),
                "escaped": float(after.get("developed_state_cumulative_escaped", 0.0)),
                "recovered": float(after.get("developed_state_cumulative_recovered", 0.0)),
            }
            # The previous accepted record, rather than prior fracture event, is the
            # exact bracket available in the neutral observer.
            before_cumulative = {
                key: float(before.get(f"developed_state_cumulative_{key}", 0.0))
                for key in cumulative
            }
            delta = {key: cumulative[key] - before_cumulative[key] for key in cumulative}
            unaccounted_removed = (
                before_totals["all_total"] + delta["emitted"] - delta["escaped"] -
                delta["recovered"] - after_totals["all_total"]
            )
            bf = state_arrays(before); af = state_arrays(after)
            before_x = np.r_[bf["x"], bf["wake_x"]]
            after_x = np.r_[af["x"], af["wake_x"]]
            before_mobile = np.r_[bf["mobile"].sum(axis=0), bf["wake_mobile"].sum(axis=0)]
            after_mobile = np.r_[af["mobile"].sum(axis=0), af["wake_mobile"].sum(axis=0)]
            before_ret = np.r_[bf["retained"].sum(axis=0), bf["wake_retained"].sum(axis=0)]
            after_ret = np.r_[af["retained"].sum(axis=0), af["wake_retained"].sum(axis=0)]
            mc0, mw0 = moments(before_x, before_mobile); mc1, mw1 = moments(after_x, after_mobile)
            rc0, rw0 = moments(before_x, before_ret); rc1, rw1 = moments(after_x, after_ret)
            avalanche_index = 0
            if pos >= onsets[1]: avalanche_index = 1
            common = {"candidate_id": cid, "original_1d_selection_role": role,
                      "pf_observed_microstructure_role": role_rows[-1]["pf_observed_microstructure_role"],
                      "event_transaction_index": transaction, "physical_avalanche_index": avalanche_index,
                      "accepted_pre_step_index": int(steps.iloc[before_pos].step),
                      "accepted_event_step_index": int(steps.iloc[pos].step),
                      "physical_onset_boundary": pos in onsets,
                      "projected_extension_before_um": extension_before * 1e6,
                      "projected_extension_after_um": extension_after * 1e6,
                      "moving_tip_translation_m": distance}
            event_rows.append({
                **common,
                **{f"before_{k}": v for k, v in before_totals.items()},
                **{f"after_{k}": v for k, v in after_totals.items()},
                "interval_emitted": delta["emitted"], "interval_escaped": delta["escaped"],
                "interval_recovered": delta["recovered"],
                "inferred_wake_discard_or_unobserved_removal": unaccounted_removed,
                "state_conservation_residual": float(
                    after_totals["all_total"] - (before_totals["all_total"] + delta["emitted"] -
                    delta["escaped"] - delta["recovered"] - unaccounted_removed)
                ),
                "mobile_centroid_shift_tip_relative_m": mc1 - mc0,
                "mobile_width_change_m": mw1 - mw0,
                "retained_centroid_shift_tip_relative_m": rc1 - rc0,
                "retained_width_change_m": rw1 - rw0,
                **translated_prior_ledger(before, distance),
                **({f"onset_after_{k}": v for k, v in regional_state(after).items()} if pos in onsets else {}),
                "observer_bracket_semantics": "PREVIOUS_ACCEPTED_ROW_TO_EVENT_CROSSING_ROW; INNER_MICROSTEP_PRECOMMIT_NOT_SERIALIZED",
                "outer_commit_state_reset": "NONE_FOR_MPZ; FIRST_PASSAGE_ACTION_SUBTRACTS_ONE",
            })
            for boundary, record, row_pos in (("PREVIOUS_ACCEPTED", before, before_pos),
                                               ("EVENT_CROSSING_ACCEPTED", after, pos)):
                detail = {**common, "projected_extension_um": float(steps.iloc[row_pos].crack_extension_m) * 1e6}
                tip, lab = expanded_state_rows(detail, record, float(steps.iloc[row_pos].a_tip_m), boundary)
                tip_rows.extend(tip); lab_rows.extend(lab)
    events_frame = pd.DataFrame(event_rows)
    tip_frame = pd.DataFrame(tip_rows); lab_frame = pd.DataFrame(lab_rows)
    roles = pd.DataFrame(role_rows)
    events_frame.to_csv(OUT / "pf_tp_event_state_translation.csv", index=False)
    flow_summary = events_frame.groupby(
        ["candidate_id", "original_1d_selection_role", "pf_observed_microstructure_role"],
        as_index=False,
    ).agg(
        event_transactions=("event_transaction_index", "count"),
        total_projected_extension_um=("projected_extension_after_um", "max"),
        operator_mobile_left_active_sum=("operator_active_mobile_after", "sum"),
        operator_retained_left_active_sum=("operator_active_retained_after", "sum"),
        operator_mobile_entering_wake_sum=("operator_mobile_entering_wake", "sum"),
        operator_retained_entering_wake_sum=("operator_retained_entering_wake", "sum"),
        operator_mobile_discarded_sum=("operator_mobile_discarded", "sum"),
        operator_retained_discarded_sum=("operator_retained_discarded", "sum"),
        accepted_interval_emitted_sum=("interval_emitted", "sum"),
        accepted_interval_escaped_sum=("interval_escaped", "sum"),
        accepted_interval_recovered_sum=("interval_recovered", "sum"),
        inferred_discard_or_unobserved_removal_sum=("inferred_wake_discard_or_unobserved_removal", "sum"),
        max_operator_mobile_conservation_residual=("operator_mobile_conservation_residual", lambda x: float(np.max(np.abs(x)))),
        max_operator_retained_conservation_residual=("operator_retained_conservation_residual", lambda x: float(np.max(np.abs(x)))),
        max_accepted_balance_residual=("state_conservation_residual", lambda x: float(np.max(np.abs(x)))),
    )
    flow_summary["persistent_source_site_renewal"] = 0.0
    flow_summary.to_csv(OUT / "pf_tp_translation_flow_summary.csv", index=False)
    tip_frame.to_parquet(OUT / "pf_tp_tip_relative_state_ledger.parquet", index=False)
    lab_frame.to_parquet(OUT / "pf_tp_lab_frame_state_ledger.parquet", index=False)
    roles.to_csv(OUT / "pf_tp_observed_microstructure_roles.csv", index=False)
    return events_frame, tip_frame, lab_frame, roles


def state_descriptor(record: dict[str, Any]) -> dict[str, float]:
    f = state_arrays(record); totals = state_totals(record)
    wake_ret = f["wake_retained"].sum(axis=0)
    wake_center, wake_width = moments(np.abs(f["wake_x"]), wake_ret)
    near = f["x"] <= 2.0e-6 + 1e-15
    return {
        **totals,
        "retained_fraction": totals["retained_total"] / max(totals["all_total"], 1e-300),
        "near_tip_mobile": float(f["mobile"][:, near].sum()),
        "near_tip_retained": float(f["retained"][:, near].sum()),
        "wake_retained_centroid_behind_tip_m": wake_center,
        "wake_retained_width_m": wake_width,
    }


def onset_state_table(cases: list[dict[str, Any]], kernel: KernelProbe,
                      roles: pd.DataFrame) -> tuple[pd.DataFrame, dict[tuple[str, str], tuple[dict[str, Any], pd.Series, dict[str, Any]]]]:
    rows = []; states: dict[tuple[str, str], tuple[dict[str, Any], pd.Series, dict[str, Any]]] = {}
    role_map = roles.set_index("candidate_id").pf_observed_microstructure_role.to_dict()
    for data in cases:
        cid = data["meta"]["candidate_id"]; original = data["meta"]["selection_role"]
        for onset_number, pos in enumerate(physical_onset_positions(data["steps"])):
            label = "INITIAL_ONSET_PRE_EVENT" if onset_number == 0 else "REINITIATION_ONSET_PRE_EVENT"
            step = data["steps"].iloc[pos]; record = data["audit"][pos]
            # The observer row carries the accepted event-crossing state and the
            # production shielding scalar is resolved on that row's committed
            # crack extension.  There is no MPZ reset at the outer commit.
            extension_m = float(step.crack_extension_m)
            native_K = float(step.KJ_Pa_sqrtm)
            probe = frozen_probe(record, data["material"], kernel,
                                 extension_m=extension_m, native_K_Pa_sqrt_m=native_K)
            rho, tau, recomputed_sigma_back = exact_backstress(record)
            sigma_back = np.asarray(record["anisotropic_sigma_back_by_system_Pa"], dtype=float)
            geom = source_geometry(record, data["material"])
            archived_shield = float(step.mpz_K_shield_Pa_sqrt_m)
            archived_back = float(step.sigma_back_Pa)
            archived_mult = float(record["persistent_site_multiplicity_per_system"])
            if not np.isclose(probe["shielding_Pa_sqrt_m"], archived_shield, rtol=5e-9, atol=5e-3):
                raise RuntimeError(f"source-exact shielding replay mismatch for {cid} {label}")
            if not np.isclose(float(np.mean(sigma_back)), archived_back, rtol=5e-9, atol=5.0):
                raise RuntimeError(f"source-exact backstress replay mismatch for {cid} {label}")
            if not np.isclose(geom["multiplicity"], archived_mult, rtol=5e-9, atol=5e-6):
                raise RuntimeError(f"source-exact multiplicity replay mismatch for {cid} {label}")
            descriptors = state_descriptor(record)
            emission_barrier = json.loads(probe["emission_barrier_by_system_eV_json"])
            emission_rates = json.loads(probe["emission_rate_by_system_s_json"])
            rows.append({
                "candidate_id": cid, "original_1d_selection_role": original,
                "pf_observed_microstructure_role": role_map[cid], "onset_role": label,
                "event_step_index": int(step.step), "physical_time_s": float(data["steps"].dt_cur_s.iloc[:pos + 1].sum()),
                "external_opening_m": float(step.Uapp_m), "reaction_N_per_m": float(step.Ftop_N),
                "native_PF_J_J_per_m2": float(step.J_effective_direct_J_per_m2),
                "native_PF_KJ_MPa_sqrt_m": native_K * 1e-6,
                "common_reference_or_apparent_K_MPa_sqrt_m": native_K * 1e-6,
                "K_shield_MPa_sqrt_m": probe["shielding_Pa_sqrt_m"] * 1e-6,
                "K_cleavage_effective_MPa_sqrt_m": (native_K - probe["shielding_Pa_sqrt_m"]) * 1e-6,
                "source_opening_equivalent_K_at_1um_MPa_sqrt_m": probe["opening_stress_Pa"] * math.sqrt(2 * math.pi * 1e-6) * 1e-6,
                "tip_radius_um": geom["radius_m"] * 1e6, "front_width_um": geom["front_width_m"] * 1e6,
                "backstress_GPa": float(np.mean(sigma_back)) * 1e-9,
                "accepted_snapshot_recomputed_backstress_GPa": float(np.mean(recomputed_sigma_back)) * 1e-9,
                "backstress_archive_phase": "LAST_SOURCE_EVALUATION; ACCEPTED_POPULATION_SNAPSHOT_RECOMPUTATION_REPORTED_SEPARATELY",
                "source_multiplicity_per_system": geom["multiplicity"],
                "opening_tensor_Pa_json": json.dumps(record.get("opening_tensor_Pa")),
                "channel_tensors_Pa_json": json.dumps(record.get("channel_tensors_Pa")),
                "resolved_signed_shears_Pa_json": json.dumps(record.get("anisotropic_tau_signed_Pa")),
                "cleavage_barrier_eV": probe["cleavage_barrier_eV"],
                "emission_barriers_eV_json": json.dumps(emission_barrier),
                "cleavage_raw_rate_s": probe["cleavage_raw_rate_s"],
                "cleavage_multihit_rate_s": probe["cleavage_effective_rate_s"],
                "emission_rates_per_site_s_json": json.dumps(emission_rates),
                "aggregate_emission_rate_s": probe["aggregate_emission_rate_s"],
                "cleavage_action_current": float(record.get("hazard_action_current", np.nan)),
                "cleavage_threshold_current": float(record.get("hazard_threshold_current_action", np.nan)),
                "peierls_transport_rate_max_s": probe["peierls_rate_max_s"],
                "encounter_retention_rate_max_s": probe["encounter_rate_max_s"],
                "taylor_completion_rate_max_s": probe["taylor_completion_rate_max_s"],
                **descriptors,
            })
            states[(cid, label)] = (record, step, data)
    table = pd.DataFrame(rows)
    table.to_csv(OUT / "pf_tp_onset_state_table.csv", index=False)
    return table, states


def swap_outcome(record: dict[str, Any], material: MaterialManifest, kernel: KernelProbe,
                 extension_m: float, native_K: float, components: dict[str, Any]) -> dict[str, Any]:
    return frozen_probe(
        record, material, kernel, extension_m=extension_m, native_K_Pa_sqrt_m=native_K,
        radius_m=components.get("radius_m"), shielding_mode=components.get("shielding_mode", "actual"),
        backstress_mode=components.get("backstress_mode", "actual"),
        multiplicity=components.get("multiplicity"), state_support=components.get("state_support", "full_active"),
    )


def build_frozen_swaps(cases: list[dict[str, Any]], kernel: KernelProbe,
                       onset_table: pd.DataFrame,
                       states: dict[tuple[str, str], tuple[dict[str, Any], pd.Series, dict[str, Any]]],
                       roles: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    control = "v913_zeroD_sobol_0202500"
    role_map = roles.set_index("candidate_id").pf_observed_microstructure_role.to_dict()
    state_swap_rows = []
    for onset_role in ("INITIAL_ONSET_PRE_EVENT", "REINITIATION_ONSET_PRE_EVENT"):
        control_row = onset_table[(onset_table.candidate_id.eq(control)) & onset_table.onset_role.eq(onset_role)].iloc[0]
        common_U = float(control_row.external_opening_m)
        common_K = float(control_row.native_PF_KJ_MPa_sqrt_m) * 1e6
        # Source state tuple retains the accepted event-crossing row.
        _control_record, cstep, control_data = states[(control, onset_role)]
        cpos = physical_onset_positions(control_data["steps"])[0 if onset_role.startswith("INITIAL") else 1]
        extension = float(cstep.crack_extension_m)
        target_cleave = float(control_row.cleavage_multihit_rate_s)
        for cid, original in EXPECTED.items():
            record, step, data = states[(cid, onset_role)]
            probe = frozen_probe(record, data["material"], kernel, extension_m=extension,
                                 native_K_Pa_sqrt_m=common_K)
            target_stress = float(states[(control, onset_role)][2]["audit"][cpos].get("sigma_cleave_eff_Pa",
                                   control_row.K_cleavage_effective_MPa_sqrt_m * 1e6 /
                                   math.sqrt(2 * math.pi * control_row.tip_radius_um * 1e-6)))
            required_K = target_stress * math.sqrt(2 * math.pi * probe["radius_m"]) + probe["shielding_Pa_sqrt_m"]
            state_swap_rows.append({
                "experiment": "STATE_SWAP_AT_FIXED_GEOMETRY_AND_OPENING",
                "geometry_owner": "CONTROL", "state_owner_candidate_id": cid,
                "original_1d_selection_role": original,
                "pf_observed_microstructure_role": role_map[cid], "onset_role": onset_role,
                "common_external_opening_m": common_U,
                "native_K_from_common_geometry_opening_MPa_sqrt_m": common_K * 1e-6,
                "opening_required_equivalent_K_for_control_cleavage_rate_MPa_sqrt_m": required_K * 1e-6,
                "control_target_cleavage_rate_s": target_cleave,
                **probe, "diagnostic_semantics": COUNTERFACTUAL,
            })
            for support in ("near_tip", "full_active", "full_production_with_wake"):
                supported = frozen_probe(record, data["material"], kernel, extension_m=extension,
                                         native_K_Pa_sqrt_m=common_K, state_support=support)
                state_swap_rows.append({
                    "experiment": "NEAR_TIP_VERSUS_FULL_WAKE_STATE",
                    "geometry_owner": "CONTROL", "state_owner_candidate_id": cid,
                    "original_1d_selection_role": original,
                    "pf_observed_microstructure_role": role_map[cid], "onset_role": onset_role,
                    "common_external_opening_m": common_U,
                    "native_K_from_common_geometry_opening_MPa_sqrt_m": common_K * 1e-6,
                    **supported, "diagnostic_semantics": COUNTERFACTUAL,
                })
    state_swap = pd.DataFrame(state_swap_rows)
    state_swap.to_csv(OUT / "pf_tp_state_swap_matrix.csv", index=False)

    # Geometry swap: exact archived linear K/U coefficient; all cases have the
    # same physical path, but retaining every owner makes that equivalence auditable.
    geometry_rows = []
    fixed_record, _, fixed_data = states[(control, "INITIAL_ONSET_PRE_EVENT")]
    fixed_U = float(onset_table[(onset_table.candidate_id.eq(control)) &
                                onset_table.onset_role.eq("INITIAL_ONSET_PRE_EVENT")].external_opening_m.iloc[0])
    for geometry_owner in EXPECTED:
        for onset_role in ("INITIAL_ONSET_PRE_EVENT", "REINITIATION_ONSET_PRE_EVENT"):
            geometry_row = onset_table[(onset_table.candidate_id.eq(geometry_owner)) &
                                       onset_table.onset_role.eq(onset_role)].iloc[0]
            coeff = float(geometry_row.native_PF_KJ_MPa_sqrt_m) * 1e6 / float(geometry_row.external_opening_m)
            native_K = coeff * fixed_U
            data = states[(geometry_owner, onset_role)][2]
            pos = physical_onset_positions(data["steps"])[0 if onset_role.startswith("INITIAL") else 1]
            extension = float(data["steps"].iloc[pos].crack_extension_m)
            probe = frozen_probe(fixed_record, fixed_data["material"], kernel,
                                 extension_m=extension, native_K_Pa_sqrt_m=native_K)
            geometry_rows.append({"experiment": "GEOMETRY_SWAP_AT_FIXED_CONTROL_INITIAL_STATE",
                                  "geometry_owner_candidate_id": geometry_owner,
                                  "geometry_onset_role": onset_role, "fixed_state_owner": control,
                                  "fixed_external_opening_m": fixed_U,
                                  "native_K_per_opening_Pa_sqrt_m_per_m": coeff,
                                  **probe, "diagnostic_semantics": COUNTERFACTUAL})
    geometry_swap = pd.DataFrame(geometry_rows)
    geometry_swap.to_csv(OUT / "pf_tp_geometry_swap_matrix.csv", index=False)

    ablation_rows = []
    for cid, original in EXPECTED.items():
        for onset_role in ("INITIAL_ONSET_PRE_EVENT", "REINITIATION_ONSET_PRE_EVENT"):
            record, step, data = states[(cid, onset_role)]
            pos = physical_onset_positions(data["steps"])[0 if onset_role.startswith("INITIAL") else 1]
            extension = float(data["steps"].iloc[pos].crack_extension_m)
            native_K = float(step.KJ_Pa_sqrtm)
            control_record = states[(control, onset_role)][0]
            ref_geom = source_geometry(control_record, data["material"])
            actual_geom = source_geometry(record, data["material"])
            configs = {
                "ACTUAL_RADIUS_ACTUAL_SHIELD_ACTUAL_BACKSTRESS_ACTUAL_MULTIPLICITY": {},
                "REFERENCE_RADIUS_ACTUAL_SHIELD_ACTUAL_BACKSTRESS": {"radius_m": ref_geom["radius_m"]},
                "ACTUAL_RADIUS_ZERO_SHIELD_ACTUAL_BACKSTRESS": {"shielding_mode": "zero"},
                "ACTUAL_RADIUS_ACTUAL_SHIELD_ZERO_BACKSTRESS": {"backstress_mode": "zero"},
                "REFERENCE_RADIUS_ZERO_SHIELD_ZERO_BACKSTRESS": {"radius_m": ref_geom["radius_m"], "shielding_mode": "zero", "backstress_mode": "zero"},
                "ACTUAL_RADIUS_ACTUAL_SHIELD_ACTUAL_BACKSTRESS_REFERENCE_MULTIPLICITY": {"multiplicity": ref_geom["multiplicity"]},
            }
            for name, config in configs.items():
                probe = swap_outcome(record, data["material"], kernel, extension, native_K, config)
                ablation_rows.append({"candidate_id": cid, "original_1d_selection_role": original,
                                      "onset_role": onset_role, "ablation": name,
                                      "factorial_member": False,
                                      **probe, "diagnostic_semantics": COUNTERFACTUAL})
            # Complete 2^4 diagnostic factorial.  This is needed to distinguish
            # one-at-a-time effects, pair interactions, and higher-order residuals
            # without pretending that nonlinear hazards form an additive energy.
            for radius_actual, shield_actual, backstress_actual, mult_actual in product((0, 1), repeat=4):
                bits = (radius_actual, shield_actual, backstress_actual, mult_actual)
                config = {
                    "radius_m": actual_geom["radius_m"] if radius_actual else ref_geom["radius_m"],
                    "shielding_mode": "actual" if shield_actual else "zero",
                    "backstress_mode": "actual" if backstress_actual else "zero",
                    "multiplicity": actual_geom["multiplicity"] if mult_actual else ref_geom["multiplicity"],
                }
                probe = swap_outcome(record, data["material"], kernel, extension, native_K, config)
                ablation_rows.append({
                    "candidate_id": cid, "original_1d_selection_role": original,
                    "onset_role": onset_role,
                    "ablation": "FACTORIAL_R%d_S%d_B%d_M%d" % bits,
                    "factorial_member": True,
                    "radius_actual": bool(radius_actual), "shielding_actual": bool(shield_actual),
                    "backstress_actual": bool(backstress_actual), "multiplicity_actual": bool(mult_actual),
                    **probe, "diagnostic_semantics": COUNTERFACTUAL,
                })
    ablations = pd.DataFrame(ablation_rows)
    ablations.to_csv(OUT / "pf_tp_component_ablation_matrix.csv", index=False)

    # Onset changes and a nonlinear hazard-competition factorial decomposition.
    decomposition_rows = []
    for cid, original in EXPECTED.items():
        first = onset_table[(onset_table.candidate_id.eq(cid)) & onset_table.onset_role.str.startswith("INITIAL")].iloc[0]
        second = onset_table[(onset_table.candidate_id.eq(cid)) & onset_table.onset_role.str.startswith("REINITIATION")].iloc[0]
        native_delta = float(second.native_PF_KJ_MPa_sqrt_m - first.native_PF_KJ_MPa_sqrt_m)
        local_delta = float(second.source_opening_equivalent_K_at_1um_MPa_sqrt_m -
                            first.source_opening_equivalent_K_at_1um_MPa_sqrt_m)
        radius_req = float(first.K_cleavage_effective_MPa_sqrt_m) * (
            math.sqrt(second.tip_radius_um / first.tip_radius_um) - 1.0
        )
        shield_delta = float(second.K_shield_MPa_sqrt_m - first.K_shield_MPa_sqrt_m)
        decomposition_rows.extend([
            {"candidate_id": cid, "original_1d_selection_role": original, "component": "NATIVE_APPARENT_DELTA_K_REINIT",
             "value": native_delta, "units": "MPa_sqrt_m", "additivity": "OBSERVED_NOT_DECOMPOSED"},
            {"candidate_id": cid, "original_1d_selection_role": original, "component": "LOCAL_SOURCE_OPENING_EQUIVALENT_DELTA_AT_REFERENCE_1UM",
             "value": local_delta, "units": "MPa_sqrt_m", "additivity": "DIAGNOSTIC_NOT_NATIVE_K"},
            {"candidate_id": cid, "original_1d_selection_role": original, "component": "RADIUS_CONSTANT_CLEAVAGE_STRESS_K_REQUIREMENT",
             "value": radius_req, "units": "MPa_sqrt_m", "additivity": "ONE_AT_A_TIME_DIAGNOSTIC"},
            {"candidate_id": cid, "original_1d_selection_role": original, "component": "SHIELDING_CHANGE",
             "value": shield_delta, "units": "MPa_sqrt_m", "additivity": "ONE_AT_A_TIME_DIAGNOSTIC"},
            {"candidate_id": cid, "original_1d_selection_role": original, "component": "GLOBAL_STRUCTURAL_MINUS_LOCAL_REFERENCE_RESIDUAL",
             "value": native_delta - local_delta, "units": "MPa_sqrt_m", "additivity": "CONTRAST_NOT_CAUSAL_SUM"},
        ])
        factorial = ablations[(ablations.candidate_id.eq(cid)) &
                              ablations.onset_role.str.startswith("REINITIATION") &
                              ablations.factorial_member.eq(True)].set_index("ablation")
        def response(bits: tuple[int, int, int, int]) -> float:
            return float(factorial.loc["FACTORIAL_R%d_S%d_B%d_M%d" % bits,
                                       "log10_cleavage_to_emission"])
        zero_bits = (0, 0, 0, 0); full_bits = (1, 1, 1, 1)
        zero = response(zero_bits); full = response(full_bits)
        names = ("radius", "shielding", "backstress", "multiplicity")
        one_effects = {}
        for i, name in enumerate(names):
            bits = [0, 0, 0, 0]; bits[i] = 1
            one_effects[name] = response(tuple(bits)) - zero
        for name, value in one_effects.items():
            decomposition_rows.append({"candidate_id": cid, "original_1d_selection_role": original,
                                       "component": f"{name.upper()}_HAZARD_COMPETITION_MAIN_EFFECT",
                                       "value": value, "units": "log10_cleavage_to_emission",
                                       "additivity": "ONE_AT_A_TIME_DIAGNOSTIC"})
        pair_effects = {}
        for i, j in combinations(range(4), 2):
            both = [0, 0, 0, 0]; both[i] = both[j] = 1
            only_i = [0, 0, 0, 0]; only_i[i] = 1
            only_j = [0, 0, 0, 0]; only_j[j] = 1
            value = response(tuple(both)) - response(tuple(only_i)) - response(tuple(only_j)) + zero
            pair_name = f"{names[i].upper()}_X_{names[j].upper()}_PAIR_INTERACTION"
            pair_effects[pair_name] = value
            decomposition_rows.append({"candidate_id": cid, "original_1d_selection_role": original,
                                       "component": pair_name, "value": value,
                                       "units": "log10_cleavage_to_emission",
                                       "additivity": "SECOND_ORDER_FACTORIAL_INTERACTION"})
        decomposition_rows.append({"candidate_id": cid, "original_1d_selection_role": original,
                                   "component": "NONADDITIVE_HIGHER_ORDER_INTERACTION_RESIDUAL",
                                   "value": full - zero - sum(one_effects.values()) - sum(pair_effects.values()),
                                   "units": "log10_cleavage_to_emission",
                                   "additivity": "THIRD_AND_FOURTH_ORDER_FACTORIAL_RESIDUAL"})
    decomposition = pd.DataFrame(decomposition_rows)
    decomposition["diagnostic_semantics"] = COUNTERFACTUAL
    decomposition.to_csv(OUT / "pf_tp_resistance_component_decomposition.csv", index=False)
    return state_swap, geometry_swap, ablations, decomposition


def selected_state_positions(data: dict[str, Any]) -> list[tuple[str, int]]:
    steps = data["steps"]; ext = steps.crack_extension_m.to_numpy(float) * 1e6
    positions: list[tuple[str, int]] = []
    for checkpoint in CHECKPOINTS_UM:
        hit = np.flatnonzero(ext >= checkpoint - 1e-9)
        positions.append((f"CHECKPOINT_{int(checkpoint):03d}UM", int(hit[0]) if len(hit) else len(steps) - 1))
    for i, pos in enumerate(physical_onset_positions(steps)):
        positions.append(("INITIAL_ONSET" if i == 0 else "REINITIATION_ONSET", int(pos)))
    return positions


def barrier_floor(surface, T_K: float) -> float:
    G0 = max(surface.G00_eV + surface.gT_eV_per_K * (T_K - surface.Tref_K), 1e-12)
    return min(surface.floor_max_fraction * G0,
               max(surface.floor_min_eV, surface.floor_fraction * G0))


def rate_flags(barrier: np.ndarray, rate: np.ndarray, surface, T_K: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    floor = barrier_floor(surface, T_K)
    on_floor = barrier <= floor + max(1e-8, 1e-6 * max(floor, 1e-12))
    saturated = rate >= 0.1 * surface.attempt_frequency_s
    inactive = rate <= 1.0e-12
    return on_floor, saturated, inactive


def state_channel_arrays(record: dict[str, Any], material: MaterialManifest, kernel: KernelProbe,
                         extension_m: float, native_K: float, T_K: float) -> dict[str, dict[str, np.ndarray]]:
    probe = frozen_probe(record, material, kernel, extension_m=extension_m,
                         native_K_Pa_sqrt_m=native_K, T_K=T_K)
    emission_stress = np.asarray(json.loads(probe["emission_stress_by_system_Pa_json"]), dtype=float)
    fields = state_arrays(record)
    cfg = MPZConfig(length_m=50e-6, n_bins=80, n_systems=2, source_bin_count=2,
                    blunting_length_m=BLUNTING_LENGTH_M, forest_density_floor_m2=5e12,
                    peierls_stress_fraction=1/math.sqrt(3), taylor_stress_fraction=1/math.sqrt(3),
                    wake_length_m=100e-6, wake_n_bins=160, wake_shielding=False)
    state = UnifiedMPZState(material, cfg); state.mobile = fields["mobile"]; state.retained = fields["retained"]
    rho = state.local_forest_density_m2(False)
    p_surface = material.peierls.as_surface(material.emission)
    t_surface = material.taylor.as_surface(material.emission)
    results: dict[str, list[np.ndarray]] = {name: [] for name in ("p_stress", "p_barrier", "p_rate", "t_stress", "t_barrier", "t_single", "t_complete", "encounter", "m", "rho")}
    for system in range(2):
        stress_profile = state.local_stress_profile_Pa(float(emission_stress[system]))
        tau_p = np.maximum(stress_profile, 0.0) / math.sqrt(3)
        spacing = 1.0 / (2.0 * np.sqrt(np.maximum(rho, 1.0)))
        phi = spacing / B_M
        tau_t = np.maximum(stress_profile, 0.0) / math.sqrt(3) * phi
        pbar = np.asarray(p_surface.values_eV(tau_p, T_K)); prate = np.asarray(p_surface.rate(tau_p, T_K))
        tbar = np.asarray(t_surface.values_eV(tau_t, T_K)); tsingle = np.asarray(t_surface.rate(tau_t, T_K))
        m = 1.0 + max(material.taylor_corr_scale, 0.0) * np.maximum(
            np.sqrt(np.maximum(rho, 0.0) / max(material.taylor_corr_rho_c_m2, 1.0)) - 1.0, 0.0)
        tcomplete = gammainc(np.maximum(m, 1.0), np.minimum(np.maximum(tsingle, 0.0), 1e12))
        velocity = spacing * prate
        encounter = max(material.encounter_efficiency, 0.0) * velocity * np.sqrt(np.maximum(rho, 0.0))
        for name, value in (("p_stress", tau_p), ("p_barrier", pbar), ("p_rate", prate),
                            ("t_stress", tau_t), ("t_barrier", tbar), ("t_single", tsingle),
                            ("t_complete", tcomplete), ("encounter", encounter), ("m", m), ("rho", rho)):
            results[name].append(np.asarray(value, dtype=float))
    out = {name: np.asarray(values) for name, values in results.items()}
    out["emission_stress"] = emission_stress
    out["emission_barrier"] = np.asarray(material.emission.values_eV(emission_stress, T_K))
    out["emission_rate"] = np.asarray(material.emission.rate(emission_stress, T_K))
    out["cleavage_stress"] = np.asarray([probe["cleavage_stress_Pa"]])
    out["cleavage_barrier"] = np.asarray([probe["cleavage_barrier_eV"]])
    out["cleavage_rate"] = np.asarray([probe["cleavage_effective_rate_s"]])
    return out


def build_barrier_audit(cases: list[dict[str, Any]], kernel: KernelProbe,
                        roles: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    actual_rows = []; grid_rows = []; landscape_rows = []; ratio_rows = []
    role_map = roles.set_index("candidate_id").pf_observed_microstructure_role.to_dict()
    for data in cases:
        cid = data["meta"]["candidate_id"]; material = data["material"]
        p_surface = material.peierls.as_surface(material.emission)
        t_surface = material.taylor.as_surface(material.emission)
        state_cache = []
        for label, pos in selected_state_positions(data):
            step = data["steps"].iloc[pos]; record = data["audit"][pos]
            extension = float(step.crack_extension_m); native_K = float(step.KJ_Pa_sqrtm)
            channels = state_channel_arrays(record, material, kernel, extension, native_K, TEMPERATURE_K)
            state_cache.append((label, record, extension, native_K, channels))
            for channel, stress, barrier, rate, surface in (
                ("PEIERLS", channels["p_stress"], channels["p_barrier"], channels["p_rate"], p_surface),
                ("TAYLOR_COMPLETION_SINGLE_HIT", channels["t_stress"], channels["t_barrier"], channels["t_single"], t_surface),
                ("EMISSION", channels["emission_stress"], channels["emission_barrier"], channels["emission_rate"], material.emission),
                ("CLEAVAGE_MULTIHIT", channels["cleavage_stress"], channels["cleavage_barrier"], channels["cleavage_rate"], material.cleavage),
            ):
                floor, sat, inactive = rate_flags(np.asarray(barrier), np.asarray(rate), surface, TEMPERATURE_K)
                for index in np.ndindex(np.asarray(stress).shape):
                    actual_rows.append({"candidate_id": cid, "original_1d_selection_role": data["meta"]["selection_role"],
                                        "pf_observed_microstructure_role": role_map[cid], "state_label": label,
                                        "projected_extension_um": extension * 1e6, "temperature_K": TEMPERATURE_K,
                                        "channel": channel, "system_index": index[0] if len(index) > 1 else (index[0] if channel == "EMISSION" else -1),
                                        "bin_index": index[1] if len(index) > 1 else -1,
                                        "effective_stress_Pa": float(np.asarray(stress)[index]),
                                        "effective_barrier_eV": float(np.asarray(barrier)[index]),
                                        "rate_s": float(np.asarray(rate)[index]),
                                        "on_barrier_floor": bool(np.asarray(floor)[index]),
                                        "rate_saturated": bool(np.asarray(sat)[index]),
                                        "rate_inactive": bool(np.asarray(inactive)[index]),
                                        "source_formula": f"{surface.__class__.__module__}.{surface.__class__.__name__}"})
            for system in range(2):
                for bin_index in range(80):
                    actual_rows.append({"candidate_id": cid, "original_1d_selection_role": data["meta"]["selection_role"],
                                        "pf_observed_microstructure_role": role_map[cid], "state_label": label,
                                        "projected_extension_um": extension * 1e6, "temperature_K": TEMPERATURE_K,
                                        "channel": "ENCOUNTER_RETENTION", "system_index": system, "bin_index": bin_index,
                                        "effective_stress_Pa": float(channels["p_stress"][system, bin_index]),
                                        "effective_barrier_eV": np.nan, "rate_s": float(channels["encounter"][system, bin_index]),
                                        "on_barrier_floor": False, "rate_saturated": False,
                                        "rate_inactive": bool(channels["encounter"][system, bin_index] <= 1e-12),
                                        "source_formula": "UnifiedMPZState._transport_rates: encounter_efficiency*v*sqrt(rho)"})
                    actual_rows.append({"candidate_id": cid, "original_1d_selection_role": data["meta"]["selection_role"],
                                        "pf_observed_microstructure_role": role_map[cid], "state_label": label,
                                        "projected_extension_um": extension * 1e6, "temperature_K": TEMPERATURE_K,
                                        "channel": "TAYLOR_COMPLETION_MULTIPLICITY", "system_index": system, "bin_index": bin_index,
                                        "effective_stress_Pa": float(channels["t_stress"][system, bin_index]),
                                        "effective_barrier_eV": float(channels["t_barrier"][system, bin_index]),
                                        "rate_s": float(channels["t_complete"][system, bin_index]),
                                        "on_barrier_floor": bool(channels["t_barrier"][system, bin_index] <= barrier_floor(t_surface, TEMPERATURE_K) + 1e-8),
                                        "rate_saturated": bool(channels["t_complete"][system, bin_index] >= .1),
                                        "rate_inactive": bool(channels["t_complete"][system, bin_index] <= 1e-12),
                                        "source_formula": "gammainc(m,taylor_single_hit_rate), m=rho-dependent"})

        # Source-exact temperature audit at characteristic observed stresses.
        for label, record, extension, native_K, channel1100 in state_cache:
            stress_stats = {
                "PEIERLS": channel1100["p_stress"], "TAYLOR_COMPLETION_SINGLE_HIT": channel1100["t_stress"],
                "EMISSION": channel1100["emission_stress"], "CLEAVAGE_MULTIHIT": channel1100["cleavage_stress"],
            }
            for T in np.arange(300.0, 1200.0 + 1e-9, 50.0):
                channels = state_channel_arrays(record, material, kernel, extension, native_K, T)
                for system in range(2):
                    emission_barrier = float(channels["emission_barrier"][system])
                    emission_rate = float(channels["emission_rate"][system])
                    for bin_index in range(80):
                        pbar = float(channels["p_barrier"][system, bin_index])
                        tbar = float(channels["t_barrier"][system, bin_index])
                        prate = float(channels["p_rate"][system, bin_index])
                        encounter = float(channels["encounter"][system, bin_index])
                        tcomplete = float(channels["t_complete"][system, bin_index])
                        ratio_rows.append({
                            "candidate_id": cid, "original_1d_selection_role": data["meta"]["selection_role"],
                            "pf_observed_microstructure_role": role_map[cid], "state_label": label,
                            "projected_extension_um": extension * 1e6, "temperature_K": T,
                            "system_index": system, "bin_index": bin_index,
                            "peierls_effective_barrier_eV": pbar,
                            "taylor_effective_barrier_eV": tbar,
                            "emission_effective_barrier_eV": emission_barrier,
                            "deltaG_peierls_over_emission": pbar / max(emission_barrier, 1e-300),
                            "deltaG_taylor_over_emission": tbar / max(emission_barrier, 1e-300),
                            "log10_rate_peierls_over_emission": math.log10(max(prate, 1e-300) / max(emission_rate, 1e-300)),
                            "log10_rate_encounter_over_emission": math.log10(max(encounter, 1e-300) / max(emission_rate, 1e-300)),
                            "log10_rate_taylor_completion_over_encounter": math.log10(max(tcomplete, 1e-300) / max(encounter, 1e-300)),
                            "same_physical_state": True,
                        })
                channel_map = {
                    "PEIERLS": (channels["p_stress"], channels["p_barrier"], channels["p_rate"], p_surface),
                    "TAYLOR_COMPLETION_SINGLE_HIT": (channels["t_stress"], channels["t_barrier"], channels["t_single"], t_surface),
                    "EMISSION": (channels["emission_stress"], channels["emission_barrier"], channels["emission_rate"], material.emission),
                    "CLEAVAGE_MULTIHIT": (channels["cleavage_stress"], channels["cleavage_barrier"], channels["cleavage_rate"], material.cleavage),
                    "ENCOUNTER_RETENTION": (channels["p_stress"], np.full_like(channels["p_stress"], np.nan), channels["encounter"], None),
                    "TAYLOR_COMPLETION_MULTIPLICITY": (channels["t_stress"], channels["t_barrier"], channels["t_complete"], t_surface),
                }
                emission_median = float(np.median(channels["emission_rate"]))
                encounter_median = float(np.median(channels["encounter"]))
                for channel, (stress, barrier, rate, surface) in channel_map.items():
                    for stat, func in (("MIN", np.min), ("MEDIAN", np.median), ("MAX", np.max)):
                        grid_rows.append({"candidate_id": cid, "state_label": label, "temperature_K": T,
                                          "channel": channel, "stress_statistic": stat,
                                          "effective_stress_Pa": float(func(stress)),
                                          "effective_barrier_eV": (float(func(barrier)) if np.isfinite(barrier).any() else np.nan),
                                          "rate_s": float(func(rate)),
                                          "log10_rate_relative_to_emission_median": math.log10(max(float(func(rate)), 1e-300) / max(emission_median, 1e-300)),
                                          "log10_taylor_completion_to_encounter_median": math.log10(max(float(np.median(channels["t_complete"])), 1e-300) / max(encounter_median, 1e-300))})

        # Smooth 1100-K barrier curves over the actually observed channel stress domain.
        all_actual = pd.DataFrame([r for r in actual_rows if r["candidate_id"] == cid])
        for channel, surface in (("PEIERLS", p_surface), ("TAYLOR_COMPLETION_SINGLE_HIT", t_surface),
                                 ("EMISSION", material.emission), ("CLEAVAGE_MULTIHIT", material.cleavage)):
            values = all_actual[all_actual.channel.eq(channel)].effective_stress_Pa.to_numpy(float)
            high = max(float(np.max(values)), 1.0)
            for stress in np.linspace(0.0, high, 161):
                barrier = float(np.asarray(surface.values_eV(stress, TEMPERATURE_K)))
                rate = float(np.asarray(surface.rate(stress, TEMPERATURE_K)))
                if channel == "CLEAVAGE_MULTIHIT":
                    rate = float(gammainc(MULTIHIT_M, min(rate * MULTIHIT_TAU_S, 1e12)) / MULTIHIT_TAU_S)
                landscape_rows.append({"candidate_id": cid, "original_1d_selection_role": data["meta"]["selection_role"],
                                       "channel": channel, "temperature_K": TEMPERATURE_K,
                                       "effective_stress_Pa": stress, "effective_barrier_eV": barrier,
                                       "rate_s": rate, "observed_stress_min_Pa": float(np.min(values)),
                                       "observed_stress_max_Pa": float(np.max(values)),
                                       "stress_free_barrier_eV": float(np.asarray(surface.values_eV(0.0, TEMPERATURE_K)))})
    actual = pd.DataFrame(actual_rows); grid = pd.DataFrame(grid_rows); landscape = pd.DataFrame(landscape_rows)
    ratios = pd.DataFrame(ratio_rows)
    actual.to_parquet(OUT / "pf_tp_actual_state_barriers_and_rates.parquet", index=False)
    grid.to_parquet(OUT / "pf_tp_barrier_landscape_temperature_grid.parquet", index=False)
    landscape.to_csv(OUT / "pf_tp_barrier_landscape_1100K.csv", index=False)
    ratios.to_parquet(OUT / "pf_tp_source_relevant_barrier_rate_ratios.parquet", index=False)

    classifications = []
    for cid, local in actual.groupby("candidate_id", sort=False):
        flags = []
        metrics: dict[str, float] = {}
        for channel in ("PEIERLS", "TAYLOR_COMPLETION_SINGLE_HIT", "TAYLOR_COMPLETION_MULTIPLICITY", "ENCOUNTER_RETENTION"):
            sub = local[local.channel.eq(channel)]
            metrics[f"{channel.lower()}_floor_fraction"] = float(sub.on_barrier_floor.mean())
            metrics[f"{channel.lower()}_saturated_fraction"] = float(sub.rate_saturated.mean())
            metrics[f"{channel.lower()}_inactive_fraction"] = float(sub.rate_inactive.mean())
            positive = sub.rate_s.to_numpy(float); positive = positive[positive > 0]
            metrics[f"{channel.lower()}_log10_rate_span"] = (
                float(np.log10(positive.max()) - np.log10(positive.min())) if len(positive) else 0.0
            )
            metrics[f"{channel.lower()}_observed_stress_min_Pa"] = float(sub.effective_stress_Pa.min())
            metrics[f"{channel.lower()}_observed_stress_max_Pa"] = float(sub.effective_stress_Pa.max())
            loaded = sub[~sub.state_label.eq("CHECKPOINT_000UM")]
            metrics[f"{channel.lower()}_loaded_floor_fraction"] = float(loaded.on_barrier_floor.mean())
            metrics[f"{channel.lower()}_loaded_saturated_fraction"] = float(loaded.rate_saturated.mean())
            metrics[f"{channel.lower()}_loaded_inactive_fraction"] = float(loaded.rate_inactive.mean())
            metrics[f"{channel.lower()}_observed_regime"] = (
                "LOADED_FLOOR_OR_RATE_ASYMPTOTE"
                if max(metrics[f"{channel.lower()}_loaded_floor_fraction"],
                       metrics[f"{channel.lower()}_loaded_saturated_fraction"]) > 0.5
                else "DYNAMIC_TRANSITION_RESOLVED"
            )
        if max(metrics["peierls_floor_fraction"], metrics["taylor_completion_single_hit_floor_fraction"]) > 0.5:
            flags.append("FLOOR_DOMINATED")
        if max(metrics["peierls_saturated_fraction"], metrics["taylor_completion_multiplicity_saturated_fraction"]) > 0.5:
            flags.append("RATE_SATURATED")
        if metrics["peierls_inactive_fraction"] > 0.5 and metrics["encounter_retention_inactive_fraction"] > 0.5:
            flags.append("RATE_INERT")
        if flags:
            primary = "EFFECTIVE_MODEL_EXTREME"
        elif cid == "v913_zeroD_sobol_0202500":
            primary = "PHYSICALLY_INTERPRETABLE_PRIMARY"
        else:
            primary = "PHYSICALLY_INTERPRETABLE_EXTREME"
        classifications.append({"candidate_id": cid,
                                "original_1d_selection_role": EXPECTED[cid],
                                "pf_observed_microstructure_role": role_map[cid],
                                "primary_plausibility_class": primary,
                                "secondary_flags_json": json.dumps(flags),
                                "observed_stress_domain_only": True,
                                "classification_uses_source_evaluated_barriers_not_raw_H0": True,
                                **metrics})
    classification = pd.DataFrame(classifications)
    classification.to_csv(OUT / "pf_tp_candidate_plausibility_classification.csv", index=False)
    return landscape, grid, actual, classification


def build_option_bank(selection: pd.DataFrame, registry: pd.DataFrame, provenance: pd.DataFrame,
                      roles: pd.DataFrame, onset: pd.DataFrame, classification: pd.DataFrame,
                      barrier_landscape: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    summary = pd.read_csv(BASE / "pf_2d_spatial_transfer_summary.csv")
    states = pd.read_parquet(BASE / "pf_2d_spatial_transfer_state_features.parquet")
    roles_i = roles.set_index("candidate_id"); class_i = classification.set_index("candidate_id")
    prov_i = provenance.set_index("candidate_id"); sel_i = selection.set_index("candidate_id")
    reg_i = registry.set_index("candidate_id")
    candidates = []
    for cid, original in EXPECTED.items():
        first = onset[(onset.candidate_id.eq(cid)) & onset.onset_role.str.startswith("INITIAL")].iloc[0]
        second = onset[(onset.candidate_id.eq(cid)) & onset.onset_role.str.startswith("REINITIATION")].iloc[0]
        final = states[(states.candidate_id.eq(cid)) & states.state_label.eq("CHECKPOINT_300UM")].iloc[0]
        material_coordinates = {
            field: builtin_scalar(reg_i.loc[cid, field])
            for field in reg_i.columns if field not in META_FIELDS
        }
        classification_row = class_i.loc[cid]
        flags = json.loads(classification_row.secondary_flags_json)
        candidates.append({
            "option_bank_schema": "Taylor_Peierls_Microstructure_Option_Bank_v1",
            "candidate_id": cid, "original_1d_selection_role": original,
            "pf_observed_microstructure_role": roles_i.loc[cid, "pf_observed_microstructure_role"],
            "full_material_vector_json": json.dumps(material_coordinates, sort_keys=True, separators=(",", ":")),
            "full_material_sha256": sel_i.loc[cid, "full_material_sha256"],
            "taylor_peierls_subvector_sha256": sel_i.loc[cid, "taylor_peierls_subvector_sha256"],
            "cleavage_barrier_sha256": sel_i.loc[cid, "cleavage_barrier_sha256"],
            "emission_barrier_sha256": sel_i.loc[cid, "emission_barrier_sha256"],
            **{field: reg_i.loc[cid, field] for field in TP_FIELDS},
            "source_evaluated_barrier_curve_artifact": "pf_tp_barrier_landscape_1100K.csv",
            "source_evaluated_barrier_curve_sha256": sha(OUT / "pf_tp_barrier_landscape_1100K.csv"),
            "source_evaluated_temperature_grid_artifact": "pf_tp_barrier_landscape_temperature_grid.parquet",
            "source_evaluated_temperature_grid_sha256": sha(OUT / "pf_tp_barrier_landscape_temperature_grid.parquet"),
            "physical_plausibility_class": classification_row.primary_plausibility_class,
            "physical_plausibility_secondary_flags_json": json.dumps(flags),
            "initial_onset_state_descriptors_json": json.dumps(first.to_dict(), default=str, sort_keys=True),
            "reinitiation_onset_state_descriptors_json": json.dumps(second.to_dict(), default=str, sort_keys=True),
            "state_300um_descriptors_json": json.dumps(final.to_dict(), default=str, sort_keys=True),
            "event_transaction_count": int(summary.set_index("candidate_id").loc[cid, "event_transaction_count"]),
            "physical_avalanche_count": int(summary.set_index("candidate_id").loc[cid, "physical_avalanche_count"]),
            "N_reinit": int(summary.set_index("candidate_id").loc[cid, "N_reinit"]),
            "largest_avalanche_fraction": float(summary.set_index("candidate_id").loc[cid, "largest_avalanche_fraction"]),
            "deltaK_reinit_MPa_sqrt_m": float(summary.set_index("candidate_id").loc[cid, "deltaK_reinit_MPa_sqrt_m"]),
            "local_equivalent_deltaK_at_reference_1um_MPa_sqrt_m": float(
                second.source_opening_equivalent_K_at_1um_MPa_sqrt_m - first.source_opening_equivalent_K_at_1um_MPa_sqrt_m),
            "PF_steps_sha256": prov_i.loc[cid, "steps_sha256"],
            "PF_observer_sha256": prov_i.loc[cid, "observer_sha256"],
            "PF_geometry_sha256": prov_i.loc[cid, "geometry_events_sha256"],
            "direct_PF_completed": True, "fatigue_evaluated": False,
            "fatigue_validation_status": "NOT_EVALUATED",
            "future_fatigue_hypothesis_status": "HYPOTHESIS_NOT_EVALUATED",
        })
    bank = pd.DataFrame(candidates)
    interpretable = bank.physical_plausibility_class.str.startswith("PHYSICALLY_INTERPRETABLE")
    shortlist_ids = bank.loc[interpretable, "candidate_id"].tolist()
    # Keep at most six while preserving actual PF-state groups and the balanced row.
    preferred = [
        "v913_zeroD_sobol_0202500", "oneD_v2_dbtt_TP_f07fe99faea93479",
        "oneD_v2_dbtt_TP_4895f9e5b44deea5", "oneD_v2_dbtt_TP_bd5be1610f6e1bce",
        "oneD_v2_dbtt_TP_9555ff54d637c974", "oneD_v2_dbtt_TP_f2817e7998cb7be6",
    ]
    shortlist_ids = [cid for cid in preferred if cid in shortlist_ids][:6]
    minimal_preferred = [
        "v913_zeroD_sobol_0202500", "oneD_v2_dbtt_TP_bd5be1610f6e1bce",
        "oneD_v2_dbtt_TP_4895f9e5b44deea5", "oneD_v2_dbtt_TP_f2817e7998cb7be6",
    ]
    minimal_ids = [cid for cid in minimal_preferred if cid in shortlist_ids]
    bank["complete_mechanism_set"] = True
    bank["physically_interpretable_shortlist"] = bank.candidate_id.isin(shortlist_ids)
    bank["minimal_future_test_set"] = bank.candidate_id.isin(minimal_ids)
    bank["tier_membership_json"] = bank.apply(lambda r: json.dumps([
        "COMPLETE_MECHANISM_SET",
        *(["PHYSICALLY_INTERPRETABLE_SHORTLIST"] if r.physically_interpretable_shortlist else []),
        *(["MINIMAL_FUTURE_TEST_SET"] if r.minimal_future_test_set else []),
    ]), axis=1)
    bank.to_csv(OUT / "taylor_peierls_microstructure_option_bank.csv", index=False)
    bank.to_parquet(OUT / "taylor_peierls_microstructure_option_bank.parquet", index=False)
    long_rows = []
    for row in bank.itertuples(index=False):
        material = json.loads(row.full_material_vector_json)
        for field, value in material.items():
            long_rows.append({"candidate_id": row.candidate_id, "parameter": field,
                              "value": value, "is_taylor_peierls_coordinate": field in TP_FIELDS,
                              "coordinate_type": "MATERIAL_COORDINATE"})
    pd.DataFrame(long_rows).to_csv(OUT / "taylor_peierls_microstructure_option_parameters_long.csv", index=False)
    feature_columns = [
        "candidate_id", "original_1d_selection_role", "pf_observed_microstructure_role",
        "physical_plausibility_class", "physical_plausibility_secondary_flags_json",
        "event_transaction_count", "physical_avalanche_count", "N_reinit",
        "largest_avalanche_fraction", "deltaK_reinit_MPa_sqrt_m",
        "local_equivalent_deltaK_at_reference_1um_MPa_sqrt_m",
        "complete_mechanism_set", "physically_interpretable_shortlist", "minimal_future_test_set",
    ]
    bank[feature_columns].to_csv(OUT / "taylor_peierls_microstructure_option_features.csv", index=False)
    bank[bank.physically_interpretable_shortlist].to_csv(
        OUT / "taylor_peierls_future_cyclic_shortlist.csv", index=False
    )
    manifest = {
        "schema": "Taylor_Peierls_Microstructure_Option_Bank_v1",
        "producer_code_commit": git_head(ROOT), "source_result_commit": "ff9e12f7c4c00e1171c1dead75f29920b089f935",
        "source_pf_commit": "f8f76435a3509553197e8d28a0e8b3cd2b9ca7ce",
        "complete_mechanism_set_candidate_ids": bank.candidate_id.tolist(),
        "physically_interpretable_shortlist_candidate_ids": shortlist_ids,
        "minimal_future_test_set_candidate_ids": minimal_ids,
        "selection_basis": "SOURCE_EVALUATED_BARRIER_RATE_PLAUSIBILITY_AND_ACTUAL_2D_PF_MICROSTRUCTURE_DIVERSITY",
        "fatigue_evaluated": False, "fatigue_validation_status": "NOT_EVALUATED",
        "future_hypothesis": (
            "Under monotonic crack growth, rapid crack advance and the current future-tip coupling make the fracture "
            "response insensitive to large variations in the mobile/retained dislocation partition. Under cyclic loading, "
            "repeated subcritical loading at a slowly advancing or temporarily stationary tip may allow those same transport "
            "and retention differences to accumulate into different cyclic blunting, shielding, return, or hazard histories."
        ),
        "future_hypothesis_status": "HYPOTHESIS_NOT_EVALUATED",
        "forbidden_non_material_coordinates_absent": [
            "PF lifecycle settings", "PF mesh/wake settings", "reduced-model event rules",
            "observer settings", "numerical tolerances", "cache settings", "avalanche grouping rules",
        ],
    }
    (OUT / "taylor_peierls_microstructure_option_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return bank, bank[bank.physically_interpretable_shortlist].copy()


def figure_style() -> None:
    mpl.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 8.5, "axes.titlesize": 9,
        "axes.labelsize": 8.5, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
        "legend.fontsize": 7, "pdf.fonttype": 42, "ps.fonttype": 42,
        "svg.fonttype": "none", "axes.spines.top": False, "axes.spines.right": False,
    })


def save_figure(fig: plt.Figure, name: str, source: pd.DataFrame) -> None:
    source.to_parquet(FIGDATA / f"{name}_source.parquet", index=False)
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(FIG / f"{name}.svg", bbox_inches="tight")
    fig.savefig(FIG / f"{name}.png", dpi=600, bbox_inches="tight")
    plt.close(fig)


def generate_figures(onset: pd.DataFrame, decomposition: pd.DataFrame,
                     state_swap: pd.DataFrame, geometry_swap: pd.DataFrame,
                     translation: pd.DataFrame, landscape: pd.DataFrame,
                     actual_rates: pd.DataFrame, classification: pd.DataFrame,
                     bank: pd.DataFrame) -> None:
    figure_style()
    summary = pd.read_csv(BASE / "pf_2d_spatial_transfer_summary.csv")
    reinit = onset[onset.onset_role.str.startswith("REINITIATION")].copy()
    source1 = reinit.merge(summary[["candidate_id", "deltaK_reinit_MPa_sqrt_m", "largest_avalanche_fraction"]], on="candidate_id")
    source1["label"] = source1.original_1d_selection_role.map(SHORT_LABEL)
    colors = [PALETTE[r] for r in source1.original_1d_selection_role]
    x = np.arange(len(source1))
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.4), constrained_layout=True)
    ax = axes[0, 0]
    ax.bar(x - .18, source1.mobile_total, .36, label="mobile", color="#56B4E9")
    ax.bar(x + .18, source1.retained_total, .36, label="retained", color="#D55E00")
    ax.set_yscale("symlog", linthresh=1); ax.set_ylabel("line content"); ax.set_title("a  Re-initiation populations")
    ax.legend(frameon=False, ncol=2); ax.set_xticks(x, [])
    ax = axes[0, 1]
    ax.bar(x, source1.wake_retained, color=colors); ax.set_yscale("log")
    ax.set_ylabel("retained line content"); ax.set_title("b  Retained wake"); ax.set_xticks(x, [])
    ax = axes[1, 0]
    radius = source1.tip_radius_um.to_numpy(); back = source1.backstress_GPa.to_numpy(); shield = source1.K_shield_MPa_sqrt_m.to_numpy()
    ax.plot(x, radius / radius.mean(), "o-", label="radius / mean", color="#0072B2")
    ax.plot(x, back / back.mean(), "s-", label="backstress / mean", color="#009E73")
    ax.plot(x, np.abs(shield) / np.abs(shield).mean(), "^-", label="|shielding| / mean", color="#CC79A7")
    ax.set_ylabel("normalized local quantity"); ax.set_title("c  Coupled near-tip scalars"); ax.legend(frameon=False)
    ax.set_xticks(x, source1.label, rotation=35, ha="right")
    ax = axes[1, 1]
    ax.bar(x, source1.deltaK_reinit_MPa_sqrt_m, color=colors)
    ax.axhline(0, color="black", lw=.7); ax.set_ylabel(r"native $\Delta K_{reinit}$ [MPa$\sqrt{m}$]")
    ax.set_title("d  Invariant two-avalanche softening")
    ax.set_xticks(x, source1.label, rotation=35, ha="right")
    save_figure(fig, "PF_TP_STATE_DIVERSITY_FRACTURE_INVARIANCE", source1)

    profiles = pd.read_parquet(BASE / "pf_2d_spatial_transfer_state_profiles.parquet")
    wanted = ["INITIAL_ONSET", "REINITIATION_ONSET_01", "CHECKPOINT_100UM", "CHECKPOINT_300UM"]
    source2 = profiles[profiles.state_label.isin(wanted)].copy()
    # Correct the legacy visualization-only profile coordinates: its wake
    # distances were stored as positive magnitudes and added to the tip.  The
    # immutable counts are unchanged; the paper source data records both the
    # recovered tip and the corrected signed laboratory coordinate.
    source2["current_tip_laboratory_m"] = source2.laboratory_position_m - source2.tip_relative_position_m
    source2["corrected_tip_relative_position_m"] = np.where(
        source2.region.eq("WAKE_BEHIND_TIP"),
        -np.abs(source2.tip_relative_position_m), np.abs(source2.tip_relative_position_m),
    )
    source2["corrected_laboratory_position_m"] = (
        source2.current_tip_laboratory_m + source2.corrected_tip_relative_position_m
    )
    source2["total"] = source2.mobile_count + source2.retained_count
    fig, axes = plt.subplots(2, 4, figsize=(7.4, 4.2), constrained_layout=True, sharey=True)
    for col, state_label in enumerate(wanted):
        local = source2[source2.state_label.eq(state_label)]
        tip_um = float(local.current_tip_laboratory_m.median() * 1e6)
        for row, field in enumerate(("mobile_count", "retained_count")):
            pivot = local.groupby(["selection_role", "corrected_tip_relative_position_m"])[field].sum().unstack(fill_value=0.0)
            pivot = pivot.reindex(EXPECTED.values())
            image = axes[row, col].imshow(np.log10(pivot.to_numpy() + 1e-3), aspect="auto", cmap="cividis",
                                          extent=[tip_um + pivot.columns.min()*1e6,
                                                  tip_um + pivot.columns.max()*1e6, 7.5, -.5])
            axes[row, col].axvline(tip_um, color="white", ls="--", lw=.9)
            axes[row, col].set_title(state_label.replace("CHECKPOINT_", "").replace("_ONSET_01", " onset").replace("INITIAL_ONSET", "initial onset"))
            axes[row, col].set_xlabel("laboratory x [um]")
            if col == 0:
                axes[row, col].set_yticks(np.arange(8), [SHORT_LABEL[r] for r in EXPECTED.values()])
                axes[row, col].set_ylabel("mobile" if row == 0 else "retained")
    cbar = fig.colorbar(image, ax=axes, shrink=.75, pad=.01); cbar.set_label(r"log$_{10}$(line content + $10^{-3}$)")
    save_figure(fig, "PF_TP_LAB_FRAME_MICROSTRUCTURE_PROFILES", source2)

    source3 = decomposition.copy()
    k_components = ["NATIVE_APPARENT_DELTA_K_REINIT", "LOCAL_SOURCE_OPENING_EQUIVALENT_DELTA_AT_REFERENCE_1UM",
                    "RADIUS_CONSTANT_CLEAVAGE_STRESS_K_REQUIREMENT", "SHIELDING_CHANGE"]
    log_components = ["RADIUS_HAZARD_COMPETITION_MAIN_EFFECT", "SHIELDING_HAZARD_COMPETITION_MAIN_EFFECT",
                      "BACKSTRESS_HAZARD_COMPETITION_MAIN_EFFECT", "MULTIPLICITY_HAZARD_COMPETITION_MAIN_EFFECT",
                      "NONADDITIVE_HIGHER_ORDER_INTERACTION_RESIDUAL"]
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.2), constrained_layout=True)
    left = source3[source3.component.isin(k_components)].pivot(index="original_1d_selection_role", columns="component", values="value").reindex(EXPECTED.values())
    left = left.reindex(columns=k_components)
    left.plot.bar(ax=axes[0], width=.8, color=["#0072B2", "#D55E00", "#009E73", "#CC79A7"])
    axes[0].axhline(0, color="black", lw=.7); axes[0].set_ylabel(r"diagnostic change [MPa$\sqrt{m}$]")
    axes[0].set_title("a  Apparent, local, radius, shielding")
    axes[0].set_xticklabels([SHORT_LABEL[r] for r in left.index], rotation=35, ha="right")
    axes[0].legend(
        handles=[Patch(color=color, label=label) for color, label in zip(
            ["#0072B2", "#D55E00", "#009E73", "#CC79A7"],
            ["native apparent", "local opening-equivalent", "radius K requirement", "shielding"],
        )], frameon=False, fontsize=6,
    )
    right = source3[source3.component.isin(log_components)].pivot(index="original_1d_selection_role", columns="component", values="value").reindex(EXPECTED.values())
    right = right.reindex(columns=log_components)
    right.plot.bar(ax=axes[1], width=.8, color=["#0072B2", "#CC79A7", "#D55E00", "#009E73", "#777777"])
    axes[1].axhline(0, color="black", lw=.7); axes[1].set_ylabel(r"change in log$_{10}(\lambda_c/\Lambda_e)$")
    axes[1].set_title("b  Nonlinear diagnostics (not additive)")
    axes[1].set_xticklabels([SHORT_LABEL[r] for r in right.index], rotation=35, ha="right")
    axes[1].legend(
        handles=[Patch(color=color, label=label) for color, label in zip(
            ["#0072B2", "#CC79A7", "#D55E00", "#009E73", "#777777"],
            ["radius", "shielding", "backstress", "multiplicity", "higher-order residual"],
        )], frameon=False, fontsize=5.8,
    )
    save_figure(fig, "PF_TP_APPARENT_LOCAL_RESISTANCE_DECOMPOSITION", source3)

    source4 = landscape.copy()
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.2), constrained_layout=True)
    for ax, channel in zip(axes.flat, ("PEIERLS", "TAYLOR_COMPLETION_SINGLE_HIT", "EMISSION", "CLEAVAGE_MULTIHIT")):
        local = source4[source4.channel.eq(channel)]
        for role in EXPECTED.values():
            row = local[local.original_1d_selection_role.eq(role)]
            if channel in {"EMISSION", "CLEAVAGE_MULTIHIT"} and role != "CONTROL": continue
            ax.plot(row.effective_stress_Pa * 1e-9, row.effective_barrier_eV,
                    color=PALETTE[role], label=SHORT_LABEL[role], lw=1.2)
            ax.axvspan(row.observed_stress_min_Pa.iloc[0] * 1e-9, row.observed_stress_max_Pa.iloc[0] * 1e-9,
                       color=PALETTE[role], alpha=.035)
        for state_label, marker in (("INITIAL_ONSET", "o"), ("REINITIATION_ONSET", "D")):
            points = actual_rates[(actual_rates.channel.eq(channel)) &
                                  actual_rates.state_label.eq(state_label)]
            points = points.groupby("original_1d_selection_role", as_index=False).agg(
                effective_stress_Pa=("effective_stress_Pa", "median"),
                effective_barrier_eV=("effective_barrier_eV", "median"),
            )
            for point in points.itertuples(index=False):
                ax.scatter(point.effective_stress_Pa * 1e-9, point.effective_barrier_eV,
                           marker=marker, s=18, facecolor=PALETTE[point.original_1d_selection_role],
                           edgecolor="white", linewidth=.35, zorder=4)
        ax.set(xlabel="source-effective stress [GPa]", ylabel="free-energy barrier [eV]", title=channel.replace("_", " ").title())
    axes[1, 1].scatter([], [], marker="o", color="0.4", label="initial onset")
    axes[1, 1].scatter([], [], marker="D", color="0.4", label="re-initiation onset")
    axes[1, 1].legend(frameon=False, fontsize=6)
    axes[0, 0].legend(frameon=False, ncol=2, fontsize=5.8)
    save_figure(fig, "PF_TP_BARRIER_LANDSCAPES_1100K", source4)

    source5 = state_swap[state_swap.experiment.eq("STATE_SWAP_AT_FIXED_GEOMETRY_AND_OPENING")].copy()
    pivot = source5.pivot(index="original_1d_selection_role", columns="onset_role", values="log10_cleavage_to_emission").reindex(EXPECTED.values())
    fig, ax = plt.subplots(figsize=(4.7, 3.7), constrained_layout=True)
    im = ax.imshow(pivot.to_numpy(), aspect="auto", cmap="coolwarm")
    ax.set_yticks(np.arange(8), [SHORT_LABEL[r] for r in pivot.index]); ax.set_xticks([0, 1], ["initial", "re-initiation"])
    ax.set_title("Frozen state swaps on control geometry/opening")
    cbar = fig.colorbar(im, ax=ax); cbar.set_label(r"log$_{10}(\lambda_c/\Lambda_e)$")
    save_figure(fig, "PF_TP_FROZEN_STATE_SWAP_MATRIX", source5)

    # SI: tip-relative profiles, translation, distance transfer, rates, plausibility, bank.
    source6 = profiles[profiles.state_label.eq("REINITIATION_ONSET_01")].copy()
    source6["corrected_tip_relative_position_m"] = np.where(
        source6.region.eq("WAKE_BEHIND_TIP"),
        -np.abs(source6.tip_relative_position_m), np.abs(source6.tip_relative_position_m),
    )
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), constrained_layout=True)
    for role in EXPECTED.values():
        local = source6[source6.selection_role.eq(role)]
        for ax, field, title in ((axes[0], "mobile_count", "mobile"), (axes[1], "retained_count", "retained")):
            curve = local.groupby("corrected_tip_relative_position_m")[field].sum()
            ax.plot(curve.index * 1e6, curve.values, color=PALETTE[role], label=SHORT_LABEL[role])
            ax.set_yscale("symlog", linthresh=1); ax.axvline(0, color="black", lw=.6)
            ax.set(xlabel="tip-relative x [um]", ylabel="line content", title=title)
    axes[1].legend(frameon=False, fontsize=5.8, ncol=2)
    save_figure(fig, "PF_TP_TIP_RELATIVE_PROFILES", source6)

    source7 = translation.copy()
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), constrained_layout=True)
    for role in EXPECTED.values():
        local = source7[source7.original_1d_selection_role.eq(role)]
        axes[0].plot(local.projected_extension_after_um, local.operator_retained_entering_wake,
                     color=PALETTE[role], label=SHORT_LABEL[role])
        axes[1].plot(local.projected_extension_after_um, local.inferred_wake_discard_or_unobserved_removal,
                     color=PALETTE[role])
    axes[0].set(xlabel="projected extension [um]", ylabel="operator-projected retained entering wake", title="moving-frame translation")
    axes[1].set(xlabel="projected extension [um]", ylabel="inferred removal over accepted bracket", title="discard/escape accounting")
    axes[0].legend(frameon=False, fontsize=5.5, ncol=2)
    save_figure(fig, "PF_TP_STATE_TRANSLATION_LEDGER", source7)

    source8 = pd.read_csv(BASE / "oneD_to_pf_spatial_pairwise_rank_transfer.csv")
    fig, ax = plt.subplots(figsize=(4.4, 3.3), constrained_layout=True)
    ax.scatter(source8.distance_1d, source8.distance_2d, color="#0072B2", s=18)
    ax.set(xlabel="1-D response-feature distance", ylabel="2-D PF state distance",
           title="Poor 1-D to 2-D distance transfer")
    save_figure(fig, "PF_TP_1D_TO_2D_DISTANCE_TRANSFER", source8)

    source9 = actual_rates.groupby(["candidate_id", "original_1d_selection_role", "state_label", "projected_extension_um", "channel"], as_index=False).rate_s.max()
    fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.8), constrained_layout=True)
    for ax, channel in zip(axes, ("PEIERLS", "ENCOUNTER_RETENTION", "TAYLOR_COMPLETION_MULTIPLICITY")):
        for role in EXPECTED.values():
            local = source9[(source9.original_1d_selection_role.eq(role)) & source9.channel.eq(channel) & source9.state_label.str.startswith("CHECKPOINT")].sort_values("projected_extension_um")
            ax.plot(local.projected_extension_um, np.log10(np.maximum(local.rate_s, 1e-300)), color=PALETTE[role])
        ax.set(xlabel="extension [um]", ylabel=r"log$_{10}$(max rate [s$^{-1}$])", title=channel.replace("_", " ").title())
    save_figure(fig, "PF_TP_RATE_TIMESCALE_EVOLUTION", source9)

    source10 = classification.copy()
    fig, ax = plt.subplots(figsize=(5.2, 3.3), constrained_layout=True)
    for row in source10.itertuples():
        ax.scatter(row.peierls_log10_rate_span, row.taylor_completion_multiplicity_log10_rate_span,
                   color=PALETTE[row.original_1d_selection_role], label=SHORT_LABEL[row.original_1d_selection_role], s=45)
    ax.set(xlabel="Peierls log-rate span", ylabel="Taylor-completion log-rate span", title="Source-evaluated plausibility map")
    ax.legend(frameon=False, fontsize=5.8, ncol=2)
    save_figure(fig, "PF_TP_CANDIDATE_PLAUSIBILITY_MAP", source10)

    source11 = bank.copy()
    reinit_map = reinit.set_index("candidate_id")
    fig, ax = plt.subplots(figsize=(5.2, 3.3), constrained_layout=True)
    for row in bank.itertuples():
        state = reinit_map.loc[row.candidate_id]
        marker = "o" if row.physically_interpretable_shortlist else "x"
        ax.scatter(state.mobile_total, state.wake_retained, color=PALETTE[row.original_1d_selection_role],
                   marker=marker, s=50, label=SHORT_LABEL[row.original_1d_selection_role])
    ax.set_xscale("symlog", linthresh=1); ax.set_yscale("log")
    ax.set(xlabel="re-initiation mobile content", ylabel="re-initiation retained wake", title="Option-bank actual PF microstructure map")
    ax.legend(frameon=False, fontsize=5.8, ncol=2)
    save_figure(fig, "PF_TP_OPTION_BANK_MICROSTRUCTURE_MAP", source11)


def write_reports(provenance: pd.DataFrame, support: pd.DataFrame, translation: pd.DataFrame,
                  onset: pd.DataFrame, state_swap: pd.DataFrame, geometry_swap: pd.DataFrame,
                  ablations: pd.DataFrame, decomposition: pd.DataFrame,
                  landscape: pd.DataFrame, actual: pd.DataFrame,
                  classification: pd.DataFrame, bank: pd.DataFrame) -> None:
    initial = onset[onset.onset_role.str.startswith("INITIAL")]
    reinit = onset[onset.onset_role.str.startswith("REINITIATION")]
    summary = pd.read_csv(BASE / "pf_2d_spatial_transfer_summary.csv")
    near_sum = reinit.near_tip_mobile + reinit.near_tip_retained
    near_rel_span = float((near_sum.max() - near_sum.min()) / near_sum.mean())
    back_rel_span = float((reinit.backstress_GPa.max() - reinit.backstress_GPa.min()) / reinit.backstress_GPa.mean())
    radius_rel_span = float((reinit.tip_radius_um.max() - reinit.tip_radius_um.min()) / reinit.tip_radius_um.mean())
    shield_rel_span = float((reinit.K_shield_MPa_sqrt_m.max() - reinit.K_shield_MPa_sqrt_m.min()) /
                            abs(reinit.K_shield_MPa_sqrt_m.mean()))
    dK_span = float(summary.deltaK_reinit_MPa_sqrt_m.max() - summary.deltaK_reinit_MPa_sqrt_m.min())
    local_delta = reinit.set_index("candidate_id").source_opening_equivalent_K_at_1um_MPa_sqrt_m - \
                  initial.set_index("candidate_id").source_opening_equivalent_K_at_1um_MPa_sqrt_m
    structural_coeff_initial = initial.native_PF_KJ_MPa_sqrt_m / initial.external_opening_m
    structural_coeff_reinit = reinit.native_PF_KJ_MPa_sqrt_m / reinit.external_opening_m
    coeff_drop = 1.0 - float((structural_coeff_reinit / structural_coeff_initial.to_numpy()).mean())

    dependency_report = f"""# PF Taylor/Peierls state-to-tip dependency audit

## Production source conclusion

The eight cases use the v10.2.21 persistent-site source with the v10.2.22 mesh-independent physical along-front-width overlay, over the v10.2.14 active-only signed shielding atlas. The production fracture law does **not** consume total or retained-wake content directly. It consumes a near-tip projection of the active state through radius, density-limited source width/multiplicity, unsigned Taylor backstress, and active signed shielding. The retained wake remains serialized and advected, but its measured shielding operator is identically zero and it is absent from backstress, radius, and multiplicity.

## Direct answers

1. Backstress uses all 80 active bins with `exp[-x/max(0.5 um, dx)]` weighting of unsigned mobile plus retained content. Here `dx=0.625 um`, so the actual weighting length is 0.625 um and the first bins dominate.
2. Shielding uses signed retained content in all 80 active bins and the extension-dependent measured physical-x kernel. Mobile shielding is zero.
3. Retained material behind the advanced tip remains in the wake ledger but is mechanically decoupled in this stack: wake shielding is disabled and wake state does not enter local source geometry.
4. Backstress and shielding are evaluated in tip-relative coordinates. The tensor probe begins in the laboratory FEM mesh and is projected into the current tip frame.
5. Active mobile, retained, signed species, and accumulated slip translate toward the moving tip.
6. Fractions crossing the tip are deposited in the laboratory-frame wake; old wake advects farther behind the new tip.
7. Persistent source sites are neither depleted nor refreshed. Population beyond the 100-um wake support is explicitly discarded; escape is separately accumulated. The outer event commit only subtracts one from first-passage action.
8. Radius uses exponentially weighted active accumulated slip, not total/wake state or net signed slip.
9. Multiplicity uses radius and a mesh-independent along-front width computed from active unsigned near-tip density; it does not use total or wake content. The v10.2.22 overlay explicitly prevents the ahead-tip MPZ `dx` from acting as the along-front width floor.
10. The sharp-wake PF structural solve does not see analytical tip radius. Radius acts only in local opening/cleavage/emission conversion and is deliberately disabled as a shielding-kernel interpolation coordinate.

At re-initiation the candidate span of the active near-tip unsigned population is only {near_rel_span:.3%}; backstress spans {back_rel_span:.3%}, radius {radius_rel_span:.3%}, and signed shielding {shield_rel_span:.3%}. Thus the orders-of-magnitude total-state diversity is largely orthogonal to the quantities consumed by the current fracture source.
"""
    (OUT / "PF_TP_STATE_TO_TIP_DEPENDENCY_AUDIT.md").write_text(dependency_report)

    trans_onset = translation[translation.physical_onset_boundary]
    max_operator_residual = max(translation.operator_mobile_conservation_residual.abs().max(),
                                translation.operator_retained_conservation_residual.abs().max())
    max_state_residual = translation.state_conservation_residual.abs().max()
    translated_mobile = translation.operator_mobile_entering_wake.groupby(translation.candidate_id).sum()
    translated_retained = translation.operator_retained_entering_wake.groupby(translation.candidate_id).sum()
    discarded = (translation.operator_mobile_discarded + translation.operator_retained_discarded).groupby(translation.candidate_id).sum()
    emitted = translation.interval_emitted.groupby(translation.candidate_id).sum()
    escaped = translation.interval_escaped.groupby(translation.candidate_id).sum()
    recovered = translation.interval_recovered.groupby(translation.candidate_id).sum()
    translation_report = f"""# PF Taylor/Peierls state translation and conservation audit

The moving process zone is translated continuously inside the Strang-coupled kinetic integrator; the outer `n_fire` commit does not apply a second population reset. The neutral observer serializes accepted-row endpoints, not every inner microstep. Therefore `pf_tp_event_state_translation.csv` distinguishes the exact accepted-row population balance from an exact application of the production moving-frame operator to the prior accepted profile.

For every one of the 472 physical transactions, the operator-projected conservation residual closes to at most {max_operator_residual:.3e} line-count units. The accepted-row source ledger closes to {max_state_residual:.3e} after explicitly accounting for emission, escape, recovery, and inferred finite-wake discard/unobserved removal.

Summed over the 59 transactions in each candidate, the exact moving-frame operator deposits {translated_mobile.min():.3g}–{translated_mobile.max():.3g} mobile and {translated_retained.min():.3g}–{translated_retained.max():.3g} retained line-count units into the wake. Finite-window discard totals {discarded.min():.3g}–{discarded.max():.3g}. The accepted-row brackets add {emitted.min():.3g}–{emitted.max():.3g} by emission, report {escaped.min():.3g}–{escaped.max():.3g} escaped and {recovered.min():.3g}–{recovered.max():.3g} recovered, and renew no persistent source sites. These sums are flow-through totals, not the final stored population.

At the two physical onset boundaries, regional columns report the active 0–2 um and 0–10 um state and wake bands 0–10, 10–50, 50–100, and beyond 100 um. The production wake support ends at 100 um, so the farther-wake row is identically zero; older state is recorded as discarded rather than silently retained. Candidate diversity at re-initiation resides chiefly in mobile-versus-retained partitioning over the active/wake ledger. The exponentially weighted near-tip unsigned sum is nearly invariant.

No stochastic trajectory was replayed. `PREVIOUS_ACCEPTED_ROW_TO_EVENT_CROSSING_ROW` is the strongest time bracket supported by the immutable observer; it is not mislabelled as an inner-microstep snapshot.
"""
    (OUT / "PF_TP_STATE_TRANSLATION_AND_CONSERVATION_AUDIT.md").write_text(translation_report)

    onset_report = f"""# PF Taylor/Peierls onset resistance decomposition

All values use V2 reload-separated pre-event onsets. Native PF KJ decreases by {summary.deltaK_reinit_MPa_sqrt_m.min():.5f} to {summary.deltaK_reinit_MPa_sqrt_m.max():.5f} MPa sqrt(m), with only {dK_span:.5f} MPa sqrt(m) candidate spread. In contrast, the explicitly labelled source-opening equivalent at a fixed 1-um reference radius rises by {local_delta.min():.5f} to {local_delta.max():.5f} MPa sqrt(m).

The sign opposition is real but the quantities are not additive. The local measure maps unshielded opening stress to a reporting reference radius. Native KJ is a sharp-wake structural diagnostic. Radius, signed shielding, backstress, and multiplicity enter nonlinear source equations. The component table therefore reports one-at-a-time frozen diagnostics, paired/higher interactions in log hazard competition, and a non-additive residual rather than claiming an energy decomposition.

The source audit reproduces archived active shielding and uses the source-evaluated backstress and persistent-site geometry saved at both onsets. The accepted population snapshot recomputes a slightly later backstress separately; the observer does not serialize the inner-microstep precommit state, so those phases are not conflated. The modest local hardening-like change is produced mainly by tip resharpening: radius falls from about {initial.tip_radius_um.mean():.4f} to {reinit.tip_radius_um.mean():.4f} um, raising K/sqrt(2*pi*r) even though native KJ falls. Signed shielding is negative (anti-shielding in the production convention) and becomes less negative; this partly opposes the resharpening increase. Backstress and multiplicity are nearly candidate invariant at the physical onsets.

The mean structural K-per-opening coefficient falls by {coeff_drop:.2%} between initial and reload-separated geometries. That sharp-wake/global-to-local geometry change dominates the remote/native response, producing apparent softening despite the local stress increase.
"""
    (OUT / "PF_TP_ONSET_RESISTANCE_DECOMPOSITION.md").write_text(onset_report)

    primary_swaps = state_swap[state_swap.experiment.eq("STATE_SWAP_AT_FIXED_GEOMETRY_AND_OPENING")]
    swap_spans = primary_swaps.groupby("onset_role").log10_cleavage_to_emission.agg(lambda x: x.max()-x.min())
    geometry_initial = geometry_swap[geometry_swap.geometry_onset_role.str.startswith("INITIAL")]
    geometry_reinit = geometry_swap[geometry_swap.geometry_onset_role.str.startswith("REINITIATION")]
    frozen_report = f"""# PF Taylor/Peierls frozen-state swap audit

All rows are `{COUNTERFACTUAL}` and were generated by deterministic evaluation of the exact committed production barrier, multihit, backstress, persistent-site, and signed-kernel functions. Production CSV/JSON trajectories were opened read-only and were never modified.

At fixed control geometry and opening, swapping the eight initial states changes log10(cleavage/aggregate-emission hazard) over a span of {swap_spans.iloc[0]:.4g}; the re-initiation-state span is {swap_spans.iloc[1]:.4g}. This confirms that state partition can change local emission competition, but not enough to alter the common two-avalanche topology.

The near-tip/full-active/full-production-with-wake rows show that adding the retained wake changes no production source quantity: the wake shielding operator is exactly zero. Differences arise only when active support beyond 2 um is restored, principally through the signed shielding kernel; backstress remains dominated by the near-tip exponential weight.

Geometry swaps use the exact archived linear K/U coefficient at each fixed sharp-wake geometry. Candidate geometries are numerically equivalent, while initial-to-post-first-avalanche geometry changes the coefficient strongly. The complete 2^4 component factorial reports one-at-a-time effects, all six pair interactions, and a third/fourth-order non-additive residual for radius, shielding, backstress, and multiplicity; these are diagnostics, not proposed production physics.
"""
    (OUT / "PF_TP_FROZEN_STATE_SWAP_AUDIT.md").write_text(frozen_report)

    class_lines = "\n".join(
        f"- {SHORT_LABEL[r.original_1d_selection_role]}: {r.primary_plausibility_class}; flags {r.secondary_flags_json}."
        for r in classification.itertuples()
    )
    barrier_report = f"""# PF Taylor/Peierls barrier plausibility audit

The audit evaluates `ExpFloorBarrier.values_eV/rate`, `TransportBarrier.as_surface`, the rho-dependent Taylor multiplicity, encounter law, fixed emission surface, fixed cleavage surface, and cleavage multihit transform. It does not classify raw H0 or activation entropy in isolation.

For every candidate and checkpoint/onset, the 1100-K table records source-effective stress, free-energy barrier, rate, barrier-floor flag, rate-saturation flag, and inactivity flag. The temperature artifact repeats characteristic observed stresses from 300 to 1200 K in 50-K increments. `pf_tp_source_relevant_barrier_rate_ratios.parquet` evaluates Peierls/emission and Taylor/emission barriers plus all requested rate ratios on the same candidate, temperature, state, system, and bin.

{class_lines}

`OUTSIDE_OBSERVED_STRESS_DOMAIN` is not assigned because classification uses only saved source-effective stress states. Diagnostic landscape curves extend only from zero to each candidate's observed maximum and do not extrapolate the shielding atlas. Large or negative activation entropy alone was not used as a rejection rule.

All eight loaded trajectories spend most sampled states on the Peierls/Taylor floor or rate asymptote. The physically interpretable shortlist is therefore empty and no four-row minimal future set is certified. This is a fail-closed handoff decision, not evidence that the eight full-set rows are unusable as diagnostic mechanism extremes.
"""
    (OUT / "PF_TP_BARRIER_PLAUSIBILITY_AUDIT.md").write_text(barrier_report)

    shortlist = bank[bank.physically_interpretable_shortlist]
    minimal = bank[bank.minimal_future_test_set]
    handoff = f"""# Fracture-to-future-cyclic Taylor/Peierls handoff

`Taylor_Peierls_Microstructure_Option_Bank_v1` contains eight fracture-derived material options. Their fatigue response is unknown: every row has `fatigue_evaluated=false` and `fatigue_validation_status=NOT_EVALUATED`.

The rows were retained because Taylor/Peierls kinetics produce strongly different 2-D mobile/retained and wake states while giving nearly identical monotonic fracture topology and reload-separated softening. This makes them controlled options for a future test of whether cyclic loading is more sensitive to transport/retention partitioning than monotonic fracture.

- Complete mechanism set: all eight candidates.
- Physically interpretable shortlist ({len(shortlist)}): {', '.join(shortlist.original_1d_selection_role.map(SHORT_LABEL)) or 'none passed the source-derived gate'}.
- Minimal future test set ({len(minimal)}): {', '.join(minimal.original_1d_selection_role.map(SHORT_LABEL)) or 'not formed because fewer than four rows passed'}.

Future hypothesis — `HYPOTHESIS_NOT_EVALUATED`: under monotonic growth, rapid advance and current future-tip coupling suppress sensitivity to mobile/retained partitioning. Under repeated subcritical loading at a slowly advancing or stationary tip, the same transport/retention differences may accumulate into different blunting, shielding, return, or hazard histories. The future result may show strong, weak, or no divergence.

No PF lifecycle setting, mesh/wake control, reduced event rule, observer option, numerical tolerance, cache key, or avalanche-grouping rule is stored as a material coordinate. No fatigue repository or run was touched.
"""
    (OUT / "FRACTURE_TO_FUTURE_CYCLIC_TAYLOR_PEIERLS_HANDOFF.md").write_text(handoff)

    claims = [
        ("Only Taylor/Peierls material coordinates differ", "SUPPORTED", "PF_TP_CASE_PROVENANCE_TABLE.csv", "all eight rows", "DBTT/1100K/seed1008666", "one material class and seed", "Only Taylor/Peierls coordinates differed in this controlled set.", "The candidates differ only in all conceivable physics."),
        ("2-D dislocation states are strongly nonunique", "SUPPORTED", "pf_tp_onset_state_table.csv", "re-initiation rows", "saved observer state", "line content is model state", "The candidates developed orders-of-magnitude differences in mobile/retained partition and wake content.", "All crack-tip physics differed by orders of magnitude."),
        ("Active near-tip fracture inputs are nearly invariant", "SUPPORTED", "pf_tp_onset_state_table.csv", "re-initiation near-tip/radius/backstress/shielding", "current source law", "signed shielding differs modestly", "Near-tip unsigned state, radius, backstress, and shielding were nearly invariant at re-initiation.", "The dislocation state has no physical effect."),
        ("Retained wake is mechanically decoupled", "SUPPORTED", "pf_tp_state_support_table.csv", "retained wake shielding row", "v10.2.14 active-only stack", "not universal to other kernels", "The serialized wake does not feed the present backstress/radius/multiplicity laws and has a zero shielding operator.", "Wake dislocations can never affect fracture."),
        ("Monotonic fracture response is nearly invariant", "SUPPORTED", "pf_2d_spatial_transfer_summary.csv", "all eight rows", "one seed to 300um", "not a universal kinetics claim", "All eight cases shared two avalanches and nearly identical reload-separated softening.", "Taylor/Peierls kinetics are universally irrelevant."),
        ("Local and native resistance changes have opposite signs", "SUPPORTED", "pf_tp_resistance_component_decomposition.csv", "native and local rows", "local quantity uses 1um reporting radius", "not additive", "Local opening-equivalent hardening and native apparent softening had opposite signs.", "The local quantity is applied remote K."),
        ("Tip resharpening creates the local hardening-like change", "SUPPORTED", "PF_TP_ONSET_RESISTANCE_DECOMPOSITION.md", "initial/reinit radius rows", "production local conversion", "interacts with anti-shielding", "Tip resharpening raises local stress despite lower native KJ.", "Dislocation hardening alone raises the remote resistance."),
        ("Sharp-wake geometry dominates native softening", "SUPPORTED", "pf_tp_geometry_swap_matrix.csv", "initial versus reinit geometry", "archived linear K/U maps", "single geometry family", "The sharp-wake structural K/U decrease dominates the native apparent response.", "Geometry is the sole possible cause in all models."),
        ("Fatigue may be more sensitive", "HYPOTHESIS_NOT_EVALUATED", "taylor_peierls_microstructure_option_manifest.json", "future_hypothesis", "future work", "fatigue not run", "These options enable a future test of cyclic sensitivity.", "Taylor/Peierls kinetics will improve fatigue resistance."),
    ]
    claim_columns = ["claim", "status", "supporting_artifact", "supporting_rows_checkpoints", "scope", "limitation", "allowed_manuscript_wording", "wording_to_avoid"]
    pd.DataFrame(claims, columns=claim_columns).to_csv(OUT / "paper_claims_and_evidence.csv", index=False)

    paper = f"""# Taylor–Peierls microstructure diversity and fracture invariance

## Methods

Eight DBTT material rows differing only in Taylor/Peierls transport and completion coordinates were evaluated at 1100 K, seed 1008666, theta=0, tip-only loading to 300 um projected extension. Cleavage and emission barriers, elasticity, source density, geometry, mesh, wake, loading, and stochastic seed were fixed. A default-off observer serialized signed mobile/retained active and wake profiles without feedback. V2 reload-separated pre-event onsets defined initial and re-initiation resistance. Frozen-state swaps evaluated committed source functions without advancing the stochastic trajectory.

## Results

The eight candidates developed two distinct re-initiation microstructure families. Total mobile content ranged from {reinit.mobile_total.min():.3g} to {reinit.mobile_total.max():.3g}, while retained wake content ranged from {reinit.wake_retained.min():.3g} to {reinit.wake_retained.max():.3g}. Nevertheless, the exponentially weighted active near-tip unsigned content varied by only {near_rel_span:.3%}; backstress, radius, and shielding varied by less than {max(back_rel_span, radius_rel_span, shield_rel_span):.3%}. Every trajectory contained 59 event transactions, two physical avalanches, one reload-separated re-initiation, and largest-avalanche fraction 0.833262. Native DeltaK_reinit ranged only from {summary.deltaK_reinit_MPa_sqrt_m.min():.4f} to {summary.deltaK_reinit_MPa_sqrt_m.max():.4f} MPa sqrt(m).

The local source-opening equivalent increased by {local_delta.min():.3f}–{local_delta.max():.3f} MPa sqrt(m), whereas native apparent KJ decreased by about 6.33 MPa sqrt(m). Source and frozen-state audits attribute the local increase chiefly to moving-frame resharpening. The structural K-per-opening coefficient fell by {coeff_drop:.2%} after the first avalanche, overwhelming the local increase. Wake-rich states did not alter the production source because the active-only measured kernel sets wake shielding to zero and the other local laws consume active near-tip state. At fixed control geometry/opening, candidate state swaps span only {swap_spans.iloc[1]:.4g} in log10(cleavage/emission hazard) at re-initiation; loaded Peierls/Taylor states are predominantly floor/rate-asymptotic.

## Physical interpretation

Taylor–Peierls kinetics generated strongly nonunique active/wake partitions, but rapid moving-frame translation and exchange drove a nearly invariant unsigned near-tip population at the reload-separated onset. The fracture source therefore sampled similar radius, multiplicity, backstress, and active signed shielding even when total and wake populations differed by orders of magnitude. Local microstructural resharpening increased opening stress, but sharp-wake structural evolution reduced model-native KJ more strongly. Predominantly asymptotic transport/completion rates further compress candidate sensitivity. The invariant fracture topology is consequently dominated by near-tip-state convergence, zero wake feedback, source-rate saturation, and structural geometry—not by proof that dislocations have no effect. Frozen swaps show finite local hazard effects and nonlinear interactions, but those effects are too small to change the observed topology.

## Limitations

This is one DBTT temperature, one seed, theta=0, tip-only mode, and 300-um monotonic growth. The v10.2.14 wake operator is zero by construction. The neutral archive stores accepted endpoints rather than every inner precommit state. All eight rows are source-rate/floor extremes at loaded states, so no physically interpretable cyclic shortlist is certified. Results do not establish universal Taylor/Peierls irrelevance and do not predict fatigue. Repeated subcritical loading could expose transport/retention sensitivity absent here; that remains `HYPOTHESIS_NOT_EVALUATED`.

## Abstract-ready statement

At fixed cleavage and emission physics, Taylor–Peierls kinetics produced orders-of-magnitude variation in 2-D mobile/retained and wake populations but nearly identical two-avalanche monotonic fracture responses. The discrepancy arises because the production fracture source samples a nearly invariant active near-tip projection while retained wake content is mechanically decoupled, and a local resharpening-driven stress increase is overwhelmed by sharp-wake structural softening.

## Figure captions

1. **Microstructure diversity with fracture-response invariance.** Re-initiation mobile/retained populations and wake content vary strongly, while near-tip scalars and native reload-separated softening remain nearly invariant.
2. **Laboratory-frame profiles.** Mobile and retained heatmaps at initial onset, re-initiation, 100 um, and 300 um show wake persistence and moving-tip redistribution.
3. **Apparent versus local decomposition.** Native KJ softening is separated from the reference-radius local opening equivalent and from counterfactual radius, shielding, backstress, multiplicity, and nonlinear residual diagnostics.
4. **Barrier landscapes.** Source-exact Peierls, Taylor, emission, and cleavage free-energy surfaces at 1100 K with observed stress domains and initial/re-initiation markers.
5. **Frozen state-swap matrix.** Local cleavage/emission hazard competition for each saved state on common control geometries/openings.

SI figures document tip-relative profiles, translation ledgers, 1-D-to-2-D distance transfer, kinetic rates, source-derived plausibility, and option-bank state coverage.
"""
    (OUT / "PAPER_TAYLOR_PEIERLS_MICROSTRUCTURE_DIVERSITY_AND_FRACTURE_INVARIANCE.md").write_text(paper)


def final_provenance() -> None:
    files = [path for path in OUT.rglob("*") if path.is_file() and path.name != "PF_TP_PAPER_AUDIT_PROVENANCE.json"]
    scientific_files = [
        path for path in files
        if path.suffix in {".csv", ".parquet", ".json"}
        and "figure_source_data" not in path.parts
    ]
    payload = {
        "schema": "PF_TP_PAPER_AUDIT_PROVENANCE_v1",
        "producer_code_commit": git_head(ROOT),
        "source_result_commit": "ff9e12f7c4c00e1171c1dead75f29920b089f935",
        "source_pf_commit": "f8f76435a3509553197e8d28a0e8b3cd2b9ca7ce",
        "kernel_family_sha256": sha(KERNEL_PATH),
        "audit_scope": {"material_class": "DBTT", "temperature_K": TEMPERATURE_K,
                        "hazard_seed": HAZARD_SEED, "theta_deg": 0.0,
                        "bulk_mode": "tip_only", "target_extension_um": TARGET_UM},
        "candidate_ids": list(EXPECTED),
        "new_memory_law": False, "cleavage_barrier_changed": False,
        "emission_barrier_changed": False, "production_trajectory_changed": False,
        "new_stochastic_pf_trajectory_count": 0, "new_femczm_run_count": 0,
        "fatigue_code_modified": False, "fatigue_calculations_run": False,
        "frozen_diagnostic_semantics": COUNTERFACTUAL,
        "scientific_fingerprint_sha256": tree_hash(scientific_files),
        "scientific_fingerprint_artifacts": [str(path.relative_to(OUT)) for path in sorted(scientific_files)],
        "artifact_sha256": {str(path.relative_to(OUT)): sha(path) for path in sorted(files)},
    }
    (OUT / "PF_TP_PAPER_AUDIT_PROVENANCE.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True); FIG.mkdir(parents=True, exist_ok=True); FIGDATA.mkdir(parents=True, exist_ok=True)
    selection, registry, matrix = load_inputs()
    cases = [load_case(case) for case in matrix["cases"]]
    provenance = build_provenance(selection, registry, matrix, cases)
    support = build_dependency_audit()
    translation, tip_ledger, lab_ledger, roles = build_translation_audit(cases)
    kernel = KernelProbe()
    onset, states = onset_state_table(cases, kernel, roles)
    state_swap, geometry_swap, ablations, decomposition = build_frozen_swaps(
        cases, kernel, onset, states, roles
    )
    landscape, temperature_grid, actual, classification = build_barrier_audit(cases, kernel, roles)
    bank, shortlist = build_option_bank(
        selection, registry, provenance, roles, onset, classification, landscape
    )
    generate_figures(
        onset, decomposition, state_swap, geometry_swap, translation,
        landscape, actual, classification, bank,
    )
    write_reports(
        provenance, support, translation, onset, state_swap, geometry_swap,
        ablations, decomposition, landscape, actual, classification, bank,
    )
    final_provenance()
    print(
        "PF_TP_PAPER_AUDIT_COMPLETE "
        f"candidates={len(bank)} shortlist={int(bank.physically_interpretable_shortlist.sum())} "
        f"minimal={int(bank.minimal_future_test_set.sum())} artifacts={sum(1 for p in OUT.rglob('*') if p.is_file())}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
