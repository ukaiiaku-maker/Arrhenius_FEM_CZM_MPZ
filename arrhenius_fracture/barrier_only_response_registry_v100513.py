"""Barrier-only four-option transfer for the full 2-D FEM/CZM model.

The v9.11.1 response registry contains both Arrhenius barrier parameters and
reduced-model state/closure parameters.  This adapter transfers only the four
barrier surfaces and their attempt frequencies.  Source inventory, source
refresh, encounter/retention, shielding, blunting, and initial state remain
owned by the existing 2-D solver configuration and are not reassigned here.
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

POINT_RELEASE = "10.0.5.13.3"
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

# Campaign controls common to all materials. The validated monotonic model is
# tip/source only: the surrounding FEM bulk remains elastic, while source
# activity, transport, retention, clearing, shielding, and renewal evolve in
# the moving crack-tip MPZ. Candidate rows do not configure those state laws.
TWO_D_STATE_POLICY: dict[str, Any] = {
    "policy_id": "preserve_existing_tip_only_moving_mpz_v1005133",
    "bulk_plasticity_mode": "tip_only",
    "mpz_length_um": 100.0,
    "mpz_n_bins": 80,
    "candidate_source_inventory_applied": False,
    "candidate_source_refresh_applied": False,
    "candidate_encounter_recovery_applied": False,
    "candidate_shielding_blunting_applied": False,
    "candidate_initial_state_applied": False,
    "state_configuration_source": "existing_tip_only_2d_solver_and_explicit_cli",
    "state_evolution_source": "existing_moving_crack_tip_MPZ",
    "continuum_bulk_role": "elastic_fem_only",
    "uniform_bulk_mobile_retained_state_active": False,
}

# The legacy v9.11 parser requires these columns even though v10.0.5.13
# intercepts the row before the old state-loading functions consume them. The
# values are compatibility placeholders only and are marked non-authoritative.
LEGACY_COMPATIBILITY_PLACEHOLDERS: dict[str, float] = {
    "source_sites_per_system": 200.0,
    "encounter_efficiency": 1.0,
    "retained_recovery_rate_s": 0.0,
    "source_refresh_length_um": 0.25,
    "c_blunt": 1.0,
    "taylor_corr_rho_c_m2": 1.0e14,
    "taylor_corr_scale": 1.0,
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
        """Return a compatibility row for the unchanged v9.11 call contract."""
        row = {
            **self.barrier_row,
            "option_key": self.option_key,
            "candidate_id": self.candidate_id,
            "target_class": self.canonical_class,
            "selection_role": "primary",
            **LEGACY_COMPATIBILITY_PLACEHOLDERS,
            "parameter_source": PARAMETER_SOURCE,
            "barrier_fingerprint_sha256": self.barrier_fingerprint_sha256,
            "two_d_state_policy_id": TWO_D_STATE_POLICY["policy_id"],
            "legacy_state_columns_are_non_authoritative_placeholders": True,
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
            "legacy_compatibility_placeholders": LEGACY_COMPATIBILITY_PLACEHOLDERS,
            "legacy_state_placeholders_consumed": False,
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
    "LEGACY_COMPATIBILITY_PLACEHOLDERS",
    "BarrierOnlyOptionV100513",
    "load_barrier_option",
]
