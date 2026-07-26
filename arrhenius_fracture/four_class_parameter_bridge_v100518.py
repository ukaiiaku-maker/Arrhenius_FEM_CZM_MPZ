"""Exact local four-class parameter transfer for FEM/CZM v10.0.5.18.

This loader reads the vendored v10.2.27 four-row registry as data only.  It does
not import or execute the PF solver.  Both the source registry SHA-256 and the
active-parameter fingerprints are verified before a row is converted to the
existing FEM/CZM persistent-site state contract.
"""
from __future__ import annotations

import csv
from dataclasses import asdict
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping

from .persistent_site_registry_v100514 import PersistentSiteRowV100514

BRIDGE_SCHEMA = "four_class_exact_local_parameter_bridge_v10_0_5_18"
REGISTRY_NAME = "v10_2_27_v913_four_class_paper_registry.csv"
MANIFEST_NAME = "four_class_transfer_manifest_v10_0_5_18.json"
SOURCE_COMMIT = "0a340f6feb8de47200b7733fe8f1d9663c75a43e"
SOURCE_REGISTRY_SHA256 = (
    "760ce58801c978b5d9f06be151b8d3e8be46debfa0ef3a8df8dee390cbcc4bf3"
)
SOURCE_WEAKT_CERAMIC_FINGERPRINT_SHA256 = (
    "eeb992cb8f956351833ff148d7611c8fc9e210f2a02eb44a00caa85558001fc5"
)
INSTALLED_FOUR_ROW_FINGERPRINT_SHA256 = (
    "923463f936ddbdb5a997de7d1859b82774df2257b1faf1dde9fef44283816010"
)

ACTIVE_FIELDS = (
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
    "rho_source0_m2",
    "taylor_corr_rho_c_m2",
    "taylor_corr_scale",
    "c_blunt",
)

EXPECTED_OPTIONS = {
    "v913_paper_peak01_0242980_persistent_sites": (
        "v913_zeroD_sobol_0242980",
        "peak",
    ),
    "v913_paper_dbtt01_0202500_persistent_sites": (
        "v913_zeroD_sobol_0202500",
        "DBTT",
    ),
    "v913_paper_weakT01_0129902_persistent_sites": (
        "v913_zeroD_sobol_0129902",
        "weakT",
    ),
    "v913_paper_ceramic01_0077080_persistent_sites": (
        "v913_zeroD_sobol_0077080",
        "ceramic",
    ),
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _float(row: Mapping[str, Any], name: str) -> float:
    try:
        value = float(row[name])
    except KeyError as exc:
        raise ValueError(f"registry row is missing required field {name!r}") from exc
    if not math.isfinite(value):
        raise ValueError(f"registry field {name!r} must be finite")
    return value


def _int(row: Mapping[str, Any], name: str) -> int:
    return int(round(_float(row, name)))


def _active_payload(row: Mapping[str, Any]) -> dict[str, Any]:
    candidate_id = str(row.get("candidate_id", ""))
    if not candidate_id:
        raise ValueError("registry row lacks candidate_id")
    payload: dict[str, Any] = {"candidate_id": candidate_id}
    for field in ACTIVE_FIELDS:
        payload[field] = _float(row, field)
    return payload


def active_fingerprint(rows: list[Mapping[str, Any]]) -> str:
    payload = [
        _active_payload(row)
        for row in sorted(rows, key=lambda item: str(item["candidate_id"]))
    ]
    encoded = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _validate_closure(row: Mapping[str, Any]) -> None:
    zero_fields = (
        "source_recovery_rate_s",
        "retained_recovery_rate_s",
        "source_refresh_length_um",
        "recovery_nu0_s",
        "recovery_H0_eV",
        "recovery_activation_entropy_kB",
        "legacy_source_sites_active",
        "legacy_source_refresh_active",
        "explicit_recovery_active",
    )
    nonzero = {
        field: _float(row, field)
        for field in zero_fields
        if abs(_float(row, field)) > 1.0e-30
    }
    if nonzero:
        raise ValueError(
            "four-class row violates the persistent-site closure: "
            + json.dumps(nonzero, sort_keys=True)
        )
    if _int(row, "n_slip_channels") != 2:
        raise ValueError("four-class FEM/CZM transfer requires two signed channels")
    for field in (
        "rho_source0_m2",
        "reference_source_area_um2",
        "reference_front_width_um",
        "source_zone_length_um",
    ):
        if _float(row, field) <= 0.0:
            raise ValueError(f"{field} must be positive")


def _candidate_from_row(row: Mapping[str, Any]) -> PersistentSiteRowV100514:
    return PersistentSiteRowV100514(
        option_key=str(row["option_key"]),
        candidate_id=str(row["candidate_id"]),
        role=str(row["role"]),
        emit_G00_eV=_float(row, "emit_G00_eV"),
        emit_gT_eV_per_K=_float(row, "emit_gT_eV_per_K"),
        peierls_H0_eV=_float(row, "peierls_H0_eV"),
        peierls_activation_entropy_kB=_float(
            row, "peierls_activation_entropy_kB"
        ),
        taylor_H0_eV=_float(row, "taylor_H0_eV"),
        taylor_activation_entropy_kB=_float(
            row, "taylor_activation_entropy_kB"
        ),
        taylor_corr_rho_c_m2=_float(row, "taylor_corr_rho_c_m2"),
        rho_source0_m2=_float(row, "rho_source0_m2"),
        source_refresh_length_um_provenance=_float(
            row, "source_refresh_length_um"
        ),
        Tref_K=_float(row, "Tref_K"),
        n_slip_channels=_int(row, "n_slip_channels"),
        rho_forest_floor_m2=_float(row, "rho_forest_floor_m2"),
        peierls_stress_fraction=_float(row, "peierls_stress_fraction"),
        taylor_stress_fraction=_float(row, "taylor_stress_fraction"),
        mobile_shield_fraction=_float(row, "mobile_shield_fraction"),
        source_recovery_rate_s=_float(row, "source_recovery_rate_s"),
        L_pz_um_recommended=_float(row, "L_pz_um_recommended"),
        n_bins_recommended=_int(row, "n_bins_recommended"),
        cleave_G00_eV=_float(row, "cleave_G00_eV"),
        cleave_gT_eV_per_K=_float(row, "cleave_gT_eV_per_K"),
        cleave_sigc0_GPa=_float(row, "cleave_sigc0_GPa"),
        cleave_sT_GPa_per_K=_float(row, "cleave_sT_GPa_per_K"),
        cleave_exp_a=_float(row, "cleave_exp_a"),
        cleave_exp_n=_float(row, "cleave_exp_n"),
        cleave_floor_frac=_float(row, "cleave_floor_frac"),
        emit_sigc0_GPa=_float(row, "emit_sigc0_GPa"),
        emit_sT_GPa_per_K=_float(row, "emit_sT_GPa_per_K"),
        emit_exp_a=_float(row, "emit_exp_a"),
        emit_exp_n=_float(row, "emit_exp_n"),
        emit_floor_frac=_float(row, "emit_floor_frac"),
        peierls_exp_a=_float(row, "peierls_exp_a"),
        peierls_exp_n=_float(row, "peierls_exp_n"),
        peierls_nu0_s=_float(row, "peierls_nu0_s"),
        taylor_exp_a=_float(row, "taylor_exp_a"),
        taylor_exp_n=_float(row, "taylor_exp_n"),
        taylor_nu0_s=_float(row, "taylor_nu0_s"),
        taylor_corr_scale=_float(row, "taylor_corr_scale"),
        source_sites_per_system_provenance=_float(
            row, "source_sites_per_system"
        ),
        encounter_efficiency=_float(row, "encounter_efficiency"),
        retained_recovery_rate_s=_float(row, "retained_recovery_rate_s"),
        c_blunt=_float(row, "c_blunt"),
        recovery_nu0_s=_float(row, "recovery_nu0_s"),
        reference_source_area_um2=_float(row, "reference_source_area_um2"),
        reference_front_width_um=_float(row, "reference_front_width_um"),
        source_zone_length_um=_float(row, "source_zone_length_um"),
        minimum_front_width_um=0.0,
    ).validate()


def load_four_class_parameter_option(
    source_root: str | Path,
    parameter_option: str,
) -> tuple[PersistentSiteRowV100514, dict[str, Any]]:
    root = Path(source_root).expanduser().resolve()
    registry_path = root / REGISTRY_NAME
    manifest_path = root / MANIFEST_NAME
    for path in (registry_path, manifest_path):
        if not path.is_file():
            raise FileNotFoundError(path)

    registry_sha = _sha256(registry_path)
    if registry_sha != SOURCE_REGISTRY_SHA256:
        raise ValueError(
            "vendored four-class registry SHA-256 mismatch: "
            f"{registry_sha} != {SOURCE_REGISTRY_SHA256}"
        )
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("source_commit") != SOURCE_COMMIT:
        raise ValueError("four-class source commit mismatch")
    if manifest.get("source_registry_sha256") != SOURCE_REGISTRY_SHA256:
        raise ValueError("four-class manifest registry SHA mismatch")

    with registry_path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 4:
        raise ValueError(f"four-class registry must contain four rows, found {len(rows)}")
    options = {str(row["option_key"]): row for row in rows}
    if set(options) != set(EXPECTED_OPTIONS):
        raise ValueError("four-class option set does not match the canonical selection")
    if len({str(row["candidate_id"]) for row in rows}) != 4:
        raise ValueError("four-class candidate IDs must be unique")
    if {str(row["material_class"]) for row in rows} != {
        "peak",
        "DBTT",
        "weakT",
        "ceramic",
    }:
        raise ValueError("four-class material labels are invalid")

    for option, (candidate_id, material_class) in EXPECTED_OPTIONS.items():
        row = options[option]
        if str(row["candidate_id"]) != candidate_id:
            raise ValueError(f"candidate mapping mismatch for {option}")
        if str(row["material_class"]) != material_class:
            raise ValueError(f"material class mismatch for {option}")
        _validate_closure(row)

    full_fingerprint = active_fingerprint(rows)
    if full_fingerprint != INSTALLED_FOUR_ROW_FINGERPRINT_SHA256:
        raise ValueError("installed four-row active-parameter fingerprint mismatch")
    transferred_rows = [
        row for row in rows if str(row["material_class"]) in {"weakT", "ceramic"}
    ]
    transfer_fingerprint = active_fingerprint(transferred_rows)
    if transfer_fingerprint != SOURCE_WEAKT_CERAMIC_FINGERPRINT_SHA256:
        raise ValueError("weakT/ceramic source active-parameter fingerprint mismatch")

    try:
        selected = options[str(parameter_option)]
    except KeyError as exc:
        raise KeyError(
            f"unknown four-class option {parameter_option!r}; "
            f"allowed={sorted(options)}"
        ) from exc
    candidate = _candidate_from_row(selected)
    selected_row_sha = hashlib.sha256(
        json.dumps(selected, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    audit = {
        "schema": BRIDGE_SCHEMA,
        "parameter_source_root": str(root),
        "registry_path": str(registry_path),
        "manifest_path": str(manifest_path),
        "source_repository": manifest.get("source_repository"),
        "source_commit": SOURCE_COMMIT,
        "source_registry_sha256": SOURCE_REGISTRY_SHA256,
        "source_weakT_ceramic_active_parameter_fingerprint_sha256": (
            SOURCE_WEAKT_CERAMIC_FINGERPRINT_SHA256
        ),
        "installed_four_row_active_parameter_fingerprint_sha256": (
            INSTALLED_FOUR_ROW_FINGERPRINT_SHA256
        ),
        "parameter_option": str(parameter_option),
        "candidate_id": candidate.candidate_id,
        "material_class": str(selected["material_class"]),
        "paper_role": str(selected["role"]),
        "selected_registry_row": dict(selected),
        "converted_candidate": asdict(candidate),
        "selected_row_sha256": selected_row_sha,
        "parameter_values_manually_reconstructed": False,
        "parameter_values_refit": False,
        "mechanics_changed_by_parameter_transfer": False,
        "source_closure_changed_by_parameter_transfer": False,
        "PF_runtime_dependency": False,
        "persistent_sites": True,
        "finite_source_inventory": False,
        "source_refresh": False,
        "explicit_recovery": False,
        "dynamic_tip_radius": True,
        "dynamic_front_width": True,
    }
    return candidate, audit


__all__ = [
    "ACTIVE_FIELDS",
    "BRIDGE_SCHEMA",
    "EXPECTED_OPTIONS",
    "INSTALLED_FOUR_ROW_FINGERPRINT_SHA256",
    "SOURCE_REGISTRY_SHA256",
    "SOURCE_WEAKT_CERAMIC_FINGERPRINT_SHA256",
    "active_fingerprint",
    "load_four_class_parameter_option",
]
