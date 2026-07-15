"""Exact v9.10.2/v9.10.3 material-row loading for the 2-D MPZ validation branch.

This module deliberately treats the selected manifest as the source of truth.  It
contains no fitted defaults for the four EXP-floor surfaces.  A run either loads
one unambiguous selected row or fails before mechanics starts.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping

TREF_K = 481.33
PRIMARY_IDS = {
    "ceramic": "ceramic_restart02_candidate00",
    "weakT": "weakT_restart00_candidate00",
    "DBTT": "DBTT_restart01_candidate05",
}

REQUIRED_COLUMNS = {
    "candidate_id", "target_class",
    "cleave_G00_eV", "cleave_gT_eV_per_K", "cleave_sigc0_GPa",
    "cleave_sT_GPa_per_K", "cleave_exp_a", "cleave_exp_n",
    "cleave_floor_frac",
    "emit_G00_eV", "emit_gT_eV_per_K", "emit_sigc0_GPa",
    "emit_sT_GPa_per_K", "emit_exp_a", "emit_exp_n", "emit_floor_frac",
    "peierls_H0_eV", "peierls_activation_entropy_kB",
    "peierls_exp_a", "peierls_exp_n",
    "taylor_H0_eV", "taylor_activation_entropy_kB",
    "taylor_exp_a", "taylor_exp_n",
    "taylor_corr_rho_c_m2", "taylor_corr_scale",
    "source_sites_per_system", "encounter_efficiency",
    "retained_recovery_rate_s", "source_refresh_length_um", "c_blunt",
}


def _coerce(value: str) -> Any:
    text = str(value).strip()
    if text == "":
        return text
    low = text.lower()
    if low in {"true", "false"}:
        return low == "true"
    try:
        return float(text)
    except ValueError:
        return text


def normalize_class_name(value: str) -> str:
    key = str(value).strip().lower().replace("-", "").replace("_", "")
    if key in {"ceramic", "ceramiclike"}:
        return "ceramic"
    if key in {"weakt", "weaktfcclike", "fcclike", "weaktemperature"}:
        return "weakT"
    if key == "dbtt":
        return "DBTT"
    raise ValueError(f"unknown MPZ class: {value!r}")


def load_selected_row(path: str | Path, class_name: str | None = None) -> dict[str, Any]:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    with path.open(newline="") as fp:
        reader = csv.DictReader(fp)
        fields = set(reader.fieldnames or [])
        missing = REQUIRED_COLUMNS - fields
        if missing:
            raise KeyError(f"{path} is missing v9.10.2/v9.10.3 columns: {sorted(missing)}")
        rows = [{k: _coerce(v) for k, v in row.items()} for row in reader]
    if not rows:
        raise RuntimeError(f"parameter manifest is empty: {path}")

    expected_class = normalize_class_name(class_name) if class_name else None
    if expected_class:
        rows = [r for r in rows if normalize_class_name(str(r["target_class"])) == expected_class]
    if not rows:
        raise RuntimeError(f"no row for class={class_name!r} in {path}")

    primary = [r for r in rows if str(r.get("selection_role", "")).lower() == "primary"]
    if len(primary) == 1:
        row = primary[0]
    else:
        cls = expected_class or normalize_class_name(str(rows[0]["target_class"]))
        cid = PRIMARY_IDS[cls]
        match = [r for r in rows if str(r["candidate_id"]) == cid]
        if len(match) != 1:
            raise RuntimeError(
                f"expected exactly one primary row ({cid}) in {path}; found {len(match)}"
            )
        row = match[0]
    row = dict(row)
    row["target_class"] = normalize_class_name(str(row["target_class"]))
    row["parameter_manifest"] = str(path.resolve())
    row["parameter_fingerprint_sha256"] = parameter_fingerprint(row)
    return row


def parameter_fingerprint(row: Mapping[str, Any]) -> str:
    payload = {k: row[k] for k in sorted(REQUIRED_COLUMNS) if k in row}
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str).encode()
    ).hexdigest()


def apply_exact_barrier_args(args: Any, row: Mapping[str, Any]) -> None:
    """Overwrite parsed sharp-front arguments with the selected exact barriers."""
    values = {
        "cleave_barrier_kind": "exp_floor",
        "cleave_exp_T_mode": "linear",
        "cleave_G00_eV": row["cleave_G00_eV"],
        "cleave_gT_eV_per_K": row["cleave_gT_eV_per_K"],
        "cleave_sigc0_GPa": row["cleave_sigc0_GPa"],
        "cleave_sT_GPa_per_K": row["cleave_sT_GPa_per_K"],
        "cleave_exp_a": row["cleave_exp_a"],
        "cleave_exp_n": row["cleave_exp_n"],
        "cleave_floor_frac": row["cleave_floor_frac"],
        "cleave_floor_min_eV": 1.0e-4,
        "cleave_floor_max_frac": 0.95,
        "cleave_Tref_K": TREF_K,
        "cleave_S_hs_kB": 0.0,
        "emit_barrier_kind": "exp_floor",
        "emit_G00_eV": row["emit_G00_eV"],
        "emit_gT_eV_per_K": row["emit_gT_eV_per_K"],
        "emit_sigc0_GPa": row["emit_sigc0_GPa"],
        "emit_sT_GPa_per_K": row["emit_sT_GPa_per_K"],
        "emit_exp_a": row["emit_exp_a"],
        "emit_exp_n": row["emit_exp_n"],
        "emit_floor_frac": row["emit_floor_frac"],
        "emit_floor_min_eV": 1.0e-4,
        "emit_floor_max_frac": 0.95,
        "emit_Tref_K": TREF_K,
        "sigma_cap_GPa": 0.0,
        "multihit_m": 3.0,
        "multihit_tau": 1.0e-6,
        "cleave_shield_chi": 0.0,
        "N_sat": math.inf,
        "n_sat": math.inf,
        "dN_cap": math.inf,
        "emb_sat_frac": 1.0,
    }
    for name, value in values.items():
        setattr(args, name, value)


def apply_pt_dislocation_config(disl_cfg: Any, row: Mapping[str, Any]) -> None:
    """Install independent v9.10.2 PT parameters on a bulk or MPZ config object."""
    emit0 = max(float(row["emit_G00_eV"]), 1.0e-30)
    values = {
        "bulk_kinetics_model": "emission_derived_peierls_taylor_multihit",
        "bulk_kinetics_model_detail": "v9102_independent_shapes",
        "pt_emit_G00_eV": float(row["emit_G00_eV"]),
        "pt_emit_gT_eV_per_K": float(row["emit_gT_eV_per_K"]),
        "pt_emit_sigc0_Pa": float(row["emit_sigc0_GPa"]) * 1.0e9,
        "pt_emit_sT_Pa_per_K": float(row["emit_sT_GPa_per_K"]) * 1.0e9,
        "pt_emit_Tref_K": TREF_K,
        "pt_emit_exp_a": float(row["emit_exp_a"]),
        "pt_emit_exp_n": float(row["emit_exp_n"]),
        "pt_emit_floor_frac": float(row["emit_floor_frac"]),
        "pt_emit_floor_min_eV": 1.0e-4,
        "pt_emit_floor_max_frac": 0.95,
        "pt_peierls_energy_ratio": float(row["peierls_H0_eV"]) / emit0,
        "pt_peierls_entropy_ratio": float(row["peierls_activation_entropy_kB"]),
        "pt_peierls_activation_entropy_kB": float(row["peierls_activation_entropy_kB"]),
        "pt_peierls_exp_a": float(row["peierls_exp_a"]),
        "pt_peierls_exp_n": float(row["peierls_exp_n"]),
        "pt_peierls_stress_ratio": 1.0,
        "pt_peierls_stress_fraction": 1.0 / math.sqrt(3.0),
        "pt_peierls_nu0_s": 1.0e12,
        "pt_taylor_energy_ratio": float(row["taylor_H0_eV"]) / emit0,
        "pt_taylor_entropy_ratio": float(row["taylor_activation_entropy_kB"]),
        "pt_taylor_activation_entropy_kB": float(row["taylor_activation_entropy_kB"]),
        "pt_taylor_exp_a": float(row["taylor_exp_a"]),
        "pt_taylor_exp_n": float(row["taylor_exp_n"]),
        "pt_taylor_stress_ratio": 1.0,
        "pt_taylor_stress_fraction": 1.0 / math.sqrt(3.0),
        "pt_taylor_nu0_s": 1.0e11,
        "pt_taylor_corr_rho_c": float(row["taylor_corr_rho_c_m2"]),
        "pt_taylor_renewal_time_s": 1.0,
        "pt_taylor_m_exponent": 1.0,
        "pt_taylor_m_scale": float(row["taylor_corr_scale"]),
        "pt_taylor_m_cap": math.inf,
        "pt_encounter_efficiency": float(row["encounter_efficiency"]),
        "pt_forest_density_floor_m2": 5.0e12,
        "pt_mobile_fraction": 0.01,
        "pt_mobile_saturation_density_m2": math.inf,
        "pt_mobile_density_floor_m2": 0.0,
        "pt_jump_fraction": 1.0,
        "pt_jump_length_min_m": 0.0,
        "pt_taylor_phi_max": math.inf,
        "rate_cap_s": math.inf,
        "dot_ep_max": math.inf,
        "rho_cap": math.inf,
    }
    for name, value in values.items():
        setattr(disl_cfg, name, value)


def build_mpz_config(args: Any, row: Mapping[str, Any]):
    from .moving_process_zone import build_mpz_config_from_namespace

    length_m = float(getattr(args, "mpz_length_um", 100.0)) * 1.0e-6
    cfg = build_mpz_config_from_namespace(args, default_length_m=length_m)
    cfg.length_m = length_m
    cfg.n_bins = int(getattr(args, "mpz_n_bins", 200))
    cfg.n_systems = 2
    cfg.source_sites_per_system = float(row["source_sites_per_system"])
    cfg.source_recovery_rate_s = 0.0
    cfg.source_refresh_length_m = float(row["source_refresh_length_um"]) * 1.0e-6
    cfg.source_bin_count = max(2, int(round(0.02 * cfg.n_bins)))
    cfg.shielding_orientation_factors = (1.0, 1.0)
    cfg.mobile_shield_fraction = 0.0
    cfg.retained_recovery_nu0_s = float(row["retained_recovery_rate_s"])
    cfg.retained_recovery_barrier_eV = 0.0
    cfg.retained_recovery_activation_volume_b3 = 0.0
    cfg.mobile_recovery_rate_s = 0.0
    cfg.pair_annihilation_rate_per_count_s = 0.0
    cfg.blunting_length_m = max(0.5e-6, 0.5 * float(getattr(args, "r_pz", 1.0e-6)))
    cfg.blunting_slip_fraction = 1.0
    apply_pt_dislocation_config(cfg, row)
    return cfg


def compact_audit(row: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "candidate_id": row["candidate_id"],
        "target_class": row["target_class"],
        "parameter_manifest": row.get("parameter_manifest"),
        "parameter_fingerprint_sha256": row.get("parameter_fingerprint_sha256") or parameter_fingerprint(row),
        "Tref_K": TREF_K,
        "independent_shape_all_four_active": True,
        "source_sites_per_system": row["source_sites_per_system"],
        "source_recovery_rate_s": 0.0,
        "source_refresh_length_um": row["source_refresh_length_um"],
        "mobile_shield_fraction": 0.0,
        "N_sat_active": False,
        "stored_energy_cleavage_active": False,
        "rate_caps_active": False,
    }
