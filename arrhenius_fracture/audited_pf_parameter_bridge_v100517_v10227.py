"""Load the exact four-row PF v10.2.27 paper registry into FEM/CZM.

This module is a parameter-only bridge layered on the validated v10.0.5.17
mechanics.  It does not import the PF package into the FEM environment.  It
reads the audited PF files as data, validates their hashes, class identities,
active-parameter fingerprint, and persistent-site closure, then converts one
selected row into the existing FEM/CZM persistent-site state contract.
"""
from __future__ import annotations

import csv
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from .persistent_site_registry_v100514 import PersistentSiteRowV100514

BRIDGE_SCHEMA = "audited_PF_v10_2_27_four_class_bridge_v10_0_5_17"
PF_BRANCH_REQUIRED = "v10.2.22-physical-front-width-top5-dbtt-screen"
PF_COMMIT_PREFIX_REQUIRED = "0a340f6"
REGISTRY_NAME = "v10_2_27_v913_four_class_paper_registry.csv"
SELECTION_NAME = "v10_2_27_v913_four_class_paper_selection.json"
MAPPING_MODULE = "sharp_front_v10_2_27.py"
AUDITED_MODULE = "sharp_front_v10_2_27_audited.py"
EXPECTED_REGISTRY_SHA256 = "760ce58801c978b5d9f06be151b8d3e8be46debfa0ef3a8df8dee390cbcc4bf3"
EXPECTED_ACTIVE_FINGERPRINT_SHA256 = "eeb992cb8f956351833ff148d7611c8fc9e210f2a02eb44a00caa85558001fc5"

EXPECTED_OPTIONS = {
    "v913_paper_peak01_0242980_persistent_sites": ("v913_zeroD_sobol_0242980", "peak"),
    "v913_paper_dbtt01_0202500_persistent_sites": ("v913_zeroD_sobol_0202500", "DBTT"),
    "v913_paper_weakT01_0129902_persistent_sites": ("v913_zeroD_sobol_0129902", "weakT"),
    "v913_paper_ceramic01_0077080_persistent_sites": ("v913_zeroD_sobol_0077080", "ceramic"),
}

ACTIVE_FIELDS = (
    "Tref_K",
    "cleave_G00_eV", "cleave_gT_eV_per_K", "cleave_sigc0_GPa",
    "cleave_sT_GPa_per_K", "cleave_exp_a", "cleave_exp_n", "cleave_floor_frac",
    "emit_G00_eV", "emit_gT_eV_per_K", "emit_sigc0_GPa",
    "emit_sT_GPa_per_K", "emit_exp_a", "emit_exp_n", "emit_floor_frac",
    "peierls_H0_eV", "peierls_activation_entropy_kB", "peierls_exp_a",
    "peierls_exp_n", "peierls_nu0_s",
    "taylor_H0_eV", "taylor_activation_entropy_kB", "taylor_exp_a",
    "taylor_exp_n", "taylor_nu0_s", "rho_source0_m2",
    "taylor_corr_rho_c_m2", "taylor_corr_scale", "c_blunt",
)

ZERO_CLOSURE_FIELDS = (
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


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _float(row: Mapping[str, Any], name: str) -> float:
    try:
        return float(row[name])
    except KeyError as exc:
        raise ValueError(f"registry row is missing required field {name!r}") from exc


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
    payload = [_active_payload(row) for row in sorted(rows, key=lambda item: str(item["candidate_id"]))]
    text = json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _validate_audited_files(mapping_path: Path, audited_path: Path) -> None:
    mapping = mapping_path.read_text()
    for token in (REGISTRY_NAME, SELECTION_NAME, "load_valid_options", "VALID_OPTIONS"):
        if token not in mapping:
            raise ValueError(f"{mapping_path} does not satisfy the v10.2.27 mapping contract: missing {token}")
    audited = audited_path.read_text()
    required = (
        "from . import sharp_front_v10_2_27 as _entry",
        "AuditedPersistentSiteStateResolvedTipEngine",
        "install_backstress_complementarity_fix",
        "install_physical_front_width",
    )
    missing = [token for token in required if token not in audited]
    if missing:
        raise ValueError(f"{audited_path} does not satisfy the audited entry contract; missing={missing}")


def _validate_rows(rows: list[dict[str, str]], selection: dict[str, Any]) -> None:
    if len(rows) != 4:
        raise ValueError(f"v10.2.27 registry must contain exactly four rows, found {len(rows)}")
    by_option = {str(row.get("option_key", "")): row for row in rows}
    if set(by_option) != set(EXPECTED_OPTIONS):
        raise ValueError(
            "v10.2.27 option set mismatch: "
            f"expected={sorted(EXPECTED_OPTIONS)} observed={sorted(by_option)}"
        )
    candidates = set()
    classes = set()
    for option, (candidate_expected, class_expected) in EXPECTED_OPTIONS.items():
        row = by_option[option]
        candidate = str(row.get("candidate_id", ""))
        material_class = str(row.get("material_class", ""))
        if candidate != candidate_expected:
            raise ValueError(f"candidate mismatch for {option}: expected={candidate_expected} observed={candidate}")
        if material_class != class_expected:
            raise ValueError(f"class mismatch for {option}: expected={class_expected} observed={material_class}")
        candidates.add(candidate)
        classes.add(material_class)
        for field in ZERO_CLOSURE_FIELDS:
            if abs(_float(row, field)) > 1.0e-30:
                raise ValueError(f"persistent-site closure requires {field}=0 for {option}")
        if _int(row, "n_slip_channels") != 2:
            raise ValueError(f"{option} must use exactly two reduced slip channels")
        if abs(_float(row, "mobile_shield_fraction")) > 1.0e-30:
            raise ValueError(f"{option} requires mobile_shield_fraction=0")
    if len(candidates) != 4 or classes != {"peak", "DBTT", "weakT", "ceramic"}:
        raise ValueError("registry candidate IDs or class labels are not the required four-class set")
    if int(selection.get("candidate_count", -1)) != 4:
        raise ValueError("selection manifest candidate_count must equal four")
    selected = {
        str(item.get("option_key")): str(item.get("candidate_id"))
        for item in selection.get("primary_candidates", [])
        if isinstance(item, dict)
    }
    expected_map = {key: value[0] for key, value in EXPECTED_OPTIONS.items()}
    if selected != expected_map:
        raise ValueError(f"selection manifest primary map mismatch: expected={expected_map} observed={selected}")


def _selection_metadata(selection: dict[str, Any], candidate_id: str) -> dict[str, Any] | None:
    for item in selection.get("primary_candidates", []):
        if isinstance(item, dict) and item.get("candidate_id") == candidate_id:
            return dict(item)
    return None


def _candidate_from_row(row: Mapping[str, Any]) -> PersistentSiteRowV100514:
    return PersistentSiteRowV100514(
        option_key=str(row["option_key"]), candidate_id=str(row["candidate_id"]), role=str(row["role"]),
        emit_G00_eV=_float(row, "emit_G00_eV"), emit_gT_eV_per_K=_float(row, "emit_gT_eV_per_K"),
        peierls_H0_eV=_float(row, "peierls_H0_eV"),
        peierls_activation_entropy_kB=_float(row, "peierls_activation_entropy_kB"),
        taylor_H0_eV=_float(row, "taylor_H0_eV"),
        taylor_activation_entropy_kB=_float(row, "taylor_activation_entropy_kB"),
        taylor_corr_rho_c_m2=_float(row, "taylor_corr_rho_c_m2"),
        rho_source0_m2=_float(row, "rho_source0_m2"),
        source_refresh_length_um_provenance=_float(row, "source_refresh_length_um"),
        Tref_K=_float(row, "Tref_K"), n_slip_channels=_int(row, "n_slip_channels"),
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
        cleave_exp_a=_float(row, "cleave_exp_a"), cleave_exp_n=_float(row, "cleave_exp_n"),
        cleave_floor_frac=_float(row, "cleave_floor_frac"),
        emit_sigc0_GPa=_float(row, "emit_sigc0_GPa"),
        emit_sT_GPa_per_K=_float(row, "emit_sT_GPa_per_K"),
        emit_exp_a=_float(row, "emit_exp_a"), emit_exp_n=_float(row, "emit_exp_n"),
        emit_floor_frac=_float(row, "emit_floor_frac"),
        peierls_exp_a=_float(row, "peierls_exp_a"), peierls_exp_n=_float(row, "peierls_exp_n"),
        peierls_nu0_s=_float(row, "peierls_nu0_s"),
        taylor_exp_a=_float(row, "taylor_exp_a"), taylor_exp_n=_float(row, "taylor_exp_n"),
        taylor_nu0_s=_float(row, "taylor_nu0_s"),
        taylor_corr_scale=_float(row, "taylor_corr_scale"),
        source_sites_per_system_provenance=_float(row, "source_sites_per_system"),
        encounter_efficiency=_float(row, "encounter_efficiency"),
        retained_recovery_rate_s=_float(row, "retained_recovery_rate_s"),
        c_blunt=_float(row, "c_blunt"), recovery_nu0_s=_float(row, "recovery_nu0_s"),
        reference_source_area_um2=_float(row, "reference_source_area_um2"),
        reference_front_width_um=_float(row, "reference_front_width_um"),
        source_zone_length_um=_float(row, "source_zone_length_um"), minimum_front_width_um=0.0,
    ).validate()


def load_audited_parameter_option(
    pf_repo_root: str | Path,
    parameter_option: str,
) -> tuple[PersistentSiteRowV100514, dict[str, Any]]:
    root = Path(pf_repo_root).expanduser().resolve()
    package = root / "arrhenius_fracture"
    materials = package / "data" / "materials"
    mapping_path = package / MAPPING_MODULE
    audited_path = package / AUDITED_MODULE
    registry_path = materials / REGISTRY_NAME
    selection_path = materials / SELECTION_NAME
    required = (mapping_path, audited_path, registry_path, selection_path)
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing audited PF v10.2.27 files: " + ", ".join(missing))
    _validate_audited_files(mapping_path, audited_path)
    registry_sha = _sha256(registry_path)
    if registry_sha != EXPECTED_REGISTRY_SHA256:
        raise ValueError(
            "PF v10.2.27 registry SHA-256 mismatch: "
            f"expected={EXPECTED_REGISTRY_SHA256} observed={registry_sha}"
        )
    with registry_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    selection = json.loads(selection_path.read_text())
    _validate_rows(rows, selection)
    fingerprint = active_fingerprint(rows)
    if fingerprint != EXPECTED_ACTIVE_FINGERPRINT_SHA256:
        raise ValueError(
            "PF v10.2.27 active-parameter fingerprint mismatch: "
            f"expected={EXPECTED_ACTIVE_FINGERPRINT_SHA256} observed={fingerprint}"
        )
    manifest_fp = str(selection.get("source_active_parameter_fingerprint_sha256", ""))
    if manifest_fp != EXPECTED_ACTIVE_FINGERPRINT_SHA256:
        raise ValueError(
            "PF selection manifest active fingerprint mismatch: "
            f"expected={EXPECTED_ACTIVE_FINGERPRINT_SHA256} observed={manifest_fp}"
        )
    if parameter_option not in EXPECTED_OPTIONS:
        raise KeyError(f"unknown v10.2.27 paper option {parameter_option!r}; allowed={sorted(EXPECTED_OPTIONS)}")
    selected = next(row for row in rows if row["option_key"] == parameter_option)
    candidate = _candidate_from_row(selected)
    audit = {
        "schema": BRIDGE_SCHEMA,
        "PF_branch_required": PF_BRANCH_REQUIRED,
        "PF_commit_prefix_required": PF_COMMIT_PREFIX_REQUIRED,
        "PF_repo_root": str(root),
        "parameter_entry": "v10.2.27",
        "audited_model_entry": str(audited_path),
        "stable_option_mapping_module": str(mapping_path),
        "parameter_option": parameter_option,
        "candidate_id": candidate.candidate_id,
        "material_class": selected.get("material_class"),
        "paper_role": selected.get("role"),
        "mechanism_summary": selected.get("mechanism_summary"),
        "validation_status": selected.get("validation_status"),
        "selection_metadata": _selection_metadata(selection, candidate.candidate_id),
        "registry_path": str(registry_path),
        "selection_path": str(selection_path),
        "audited_entry_sha256": _sha256(audited_path),
        "mapping_module_sha256": _sha256(mapping_path),
        "registry_sha256": registry_sha,
        "selection_sha256": _sha256(selection_path),
        "active_parameter_fingerprint_sha256": fingerprint,
        "selected_registry_row": dict(selected),
        "converted_candidate": asdict(candidate),
        "parameter_values_manually_reconstructed": False,
        "selected_through_audited_v10_2_27_registry": True,
        "mechanics_changed": False,
        "source_closure_changed": False,
        "stochastic_cleavage_law_changed": False,
        "persistent_sites": True,
        "finite_source_inventory": False,
        "source_refresh": False,
        "explicit_recovery": False,
    }
    audit["selected_row_sha256"] = hashlib.sha256(
        json.dumps(selected, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    return candidate, audit


__all__ = [
    "ACTIVE_FIELDS", "BRIDGE_SCHEMA", "EXPECTED_ACTIVE_FINGERPRINT_SHA256",
    "EXPECTED_OPTIONS", "EXPECTED_REGISTRY_SHA256", "PF_BRANCH_REQUIRED",
    "PF_COMMIT_PREFIX_REQUIRED", "active_fingerprint", "load_audited_parameter_option",
]
