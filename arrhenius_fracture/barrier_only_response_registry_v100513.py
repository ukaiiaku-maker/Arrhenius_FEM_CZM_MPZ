"""Barrier-only four-option transfer for the full 2-D FEM/CZM model.

The v9.11.1 response registry contains both Arrhenius barrier parameters and
reduced-model state/closure parameters.  This adapter transfers only the four
barrier surfaces and their attempt frequencies.  Source inventory, source
refresh, encounter/retention, shielding, blunting, MPZ discretization, and
initial state remain one common 2-D solver policy for every option.
"""
from __future__ import annotations

from dataclasses import dataclass
import csv
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from .mpz_response_registry_v100512 import (
    PRIMARY_OPTION_KEYS,
    default_registry_path,
    load_option,
    normalize_option_key,
)

POINT_RELEASE = "10.0.5.13"
PARAMETER_SOURCE = "mpz_v9_11_1_barriers_only"

BARRIER_FIELDS = (
    "Tref_K",
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
    "peierls_activation_entropy_kB",
    "peierls_exp_a",
    "peierls_exp_n",
    "peierls_nu0_s",
    "taylor_H0_eV",
    "taylor_activation_entropy_kB",
    "taylor_exp_a",
    "taylor_exp_n",
    "taylor_nu0_s",
)

IGNORED_CANDIDATE_STATE_FIELDS = (
    "source_sites_per_system",
    "encounter_efficiency",
    "retained_recovery_rate_s",
    "source_refresh_length_um",
    "c_blunt",
    "mobile_shield_fraction",
    "source_recovery_rate_s",
    "L_pz_um_recommended",
    "n_bins_recommended",
    "rho_forest_floor_m2",
    "peierls_stress_fraction",
    "taylor_stress_fraction",
    "taylor_corr_rho_c_m2",
    "taylor_corr_scale",
)

# Common solver-native state policy.  These values are independent of response
# option and reproduce the pre-parameter-overlay 2-D defaults rather than the
# candidate-specific 0-D/1-D state fit.  The state itself still evolves during
# the FEM/CZM solve.
TWO_D_STATE_POLICY: dict[str, Any] = {
    "policy_id": "full_2d_common_state_v100513",
    "bulk_plasticity_mode": "bulk_same_pt_km",
    "mpz_length_um": 100.0,
    "mpz_n_bins": 80,
    "source_sites_per_system": 200.0,
    "source_recovery_rate_s": 0.0,
    "source_refresh_length_um": 0.25,
    "source_bin_count": 2,
    "shielding_orientation_factors": [1.0, 1.0],
    "mobile_shield_fraction": 0.0,
    "shielding_core_m": 2.5e-10,
    "retained_recovery_nu0_s": 1.0e9,
    "retained_recovery_barrier_eV": 1.50,
    "retained_recovery_activation_volume_b3": 0.0,
    "mobile_recovery_rate_s": 0.0,
    "pair_annihilation_rate_per_count_s": 0.0,
    "blunting_length_um": 0.5,
    "blunting_slip_fraction": 1.0,
    "c_blunt": 1.0,
    "encounter_efficiency": 1.0,
    "forest_density_floor_m2": 5.0e12,
    "peierls_stress_fraction": 1.0 / (3.0 ** 0.5),
    "taylor_stress_fraction": 1.0 / (3.0 ** 0.5),
    "taylor_corr_rho_c_m2": 1.0e14,
    "taylor_renewal_time_s": 1.0e-9,
    "taylor_m_exponent": 1.0,
    "taylor_m_scale": 1.0,
    "taylor_m_cap": float("inf"),
}


def _fingerprint(payload: Mapping[str, Any]) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str).encode()
    ).hexdigest()


@dataclass(frozen=True)
class BarrierOnlyOptionV100513:
    option_key: str
    candidate_id: str
    canonical_class: str
    barrier_row: dict[str, Any]
    barrier_fingerprint_sha256: str
    source_registry_path: str
    ignored_candidate_state: dict[str, Any]

    def legacy_row(self, manifest_path: str | None = None) -> dict[str, Any]:
        """Return a compatibility row for the unchanged v9.11 parser.

        State fields are populated from the common 2-D policy solely because the
        legacy parser requires those column names.  They never come from the
        candidate row.
        """
        row = {
            **self.barrier_row,
            "option_key": self.option_key,
            "candidate_id": self.candidate_id,
            "target_class": self.canonical_class,
            "selection_role": "primary",
            "source_sites_per_system": TWO_D_STATE_POLICY["source_sites_per_system"],
            "encounter_efficiency": TWO_D_STATE_POLICY["encounter_efficiency"],
            # Compatibility-only scalar.  The v10.0.5.13 config builder retains
            # the solver's activated recovery law instead of consuming this.
            "retained_recovery_rate_s": 0.0,
            "source_refresh_length_um": TWO_D_STATE_POLICY["source_refresh_length_um"],
            "c_blunt": TWO_D_STATE_POLICY["c_blunt"],
            "taylor_corr_rho_c_m2": TWO_D_STATE_POLICY["taylor_corr_rho_c_m2"],
            "taylor_corr_scale": TWO_D_STATE_POLICY["taylor_m_scale"],
            "parameter_source": PARAMETER_SOURCE,
            "barrier_fingerprint_sha256": self.barrier_fingerprint_sha256,
            "two_d_state_policy_id": TWO_D_STATE_POLICY["policy_id"],
        }
        if manifest_path is not None:
            row["parameter_manifest"] = str(manifest_path)
        row["parameter_fingerprint_sha256"] = _fingerprint(
            {
                "barriers": self.barrier_row,
                "two_d_state_policy": TWO_D_STATE_POLICY,
                "candidate_id": self.candidate_id,
            }
        )
        return row

    def write_barrier_csv(self, path: str | Path) -> Path:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        row = {
            "option_key": self.option_key,
            "candidate_id": self.candidate_id,
            "target_class": self.canonical_class,
            "parameter_source": PARAMETER_SOURCE,
            "barrier_fingerprint_sha256": self.barrier_fingerprint_sha256,
            **self.barrier_row,
        }
        with target.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(row))
            writer.writeheader()
            writer.writerow(row)
        return target

    def audit_payload(self) -> dict[str, Any]:
        return {
            "point_release": POINT_RELEASE,
            "parameter_source": PARAMETER_SOURCE,
            "option_key": self.option_key,
            "candidate_id": self.candidate_id,
            "canonical_class": self.canonical_class,
            "barrier_fields_transferred": list(BARRIER_FIELDS),
            "barrier_fingerprint_sha256": self.barrier_fingerprint_sha256,
            "source_registry_path": self.source_registry_path,
            "candidate_state_fields_ignored": self.ignored_candidate_state,
            "two_d_state_policy": TWO_D_STATE_POLICY,
        }


def load_barrier_option(
    value: str,
    path: str | Path | None = None,
) -> BarrierOnlyOptionV100513:
    option_key = normalize_option_key(value)
    source_path = Path(path) if path is not None else default_registry_path()
    source = load_option(option_key, source_path)
    missing = [name for name in BARRIER_FIELDS if name not in source.row]
    if missing:
        raise RuntimeError(f"{option_key} lacks barrier-only fields: {missing}")
    barrier_row = {name: source.row[name] for name in BARRIER_FIELDS}
    ignored = {
        name: source.row.get(name)
        for name in IGNORED_CANDIDATE_STATE_FIELDS
        if name in source.row
    }
    fingerprint = _fingerprint(
        {
            "option_key": option_key,
            "candidate_id": source.candidate_id,
            "barriers": barrier_row,
        }
    )
    return BarrierOnlyOptionV100513(
        option_key=option_key,
        candidate_id=source.candidate_id,
        canonical_class=source.canonical_class,
        barrier_row=barrier_row,
        barrier_fingerprint_sha256=fingerprint,
        source_registry_path=str(source_path.resolve()),
        ignored_candidate_state=ignored,
    )


__all__ = [
    "POINT_RELEASE",
    "PARAMETER_SOURCE",
    "PRIMARY_OPTION_KEYS",
    "BARRIER_FIELDS",
    "IGNORED_CANDIDATE_STATE_FIELDS",
    "TWO_D_STATE_POLICY",
    "BarrierOnlyOptionV100513",
    "load_barrier_option",
]
