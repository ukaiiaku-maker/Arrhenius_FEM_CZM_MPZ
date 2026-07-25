"""Load exact paper-campaign rows from audited PF model entries.

The FEM/CZM repository does not duplicate the v10.2.25/v10.2.26 numerical
registries.  Instead this bridge validates the requested audited PF entry and its
stable option map, reads the selected installed CSV row from the supplied PF
source tree, verifies the persistent-site closure, and converts that exact row to
the existing FEM/CZM persistent-site state contract.
"""
from __future__ import annotations

import ast
import csv
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
from typing import Any

from .persistent_site_registry_v100514 import PersistentSiteRowV100514

BRIDGE_SCHEMA = "audited_PF_v10_2_25_v10_2_26_parameter_bridge_v10_0_5_17"
PF_BRANCH_REQUIRED = "v10.2.22-physical-front-width-top5-dbtt-screen"


@dataclass(frozen=True)
class AuditedEntrySpec:
    cli_key: str
    audited_module: str
    mapping_module: str
    registry_name: str
    selection_name: str


ENTRY_SPECS = {
    "v10.2.25": AuditedEntrySpec(
        cli_key="v10.2.25",
        audited_module="sharp_front_v10_2_25_audited.py",
        mapping_module="sharp_front_v10_2_25.py",
        registry_name="v10_2_25_v913_paper_campaign_registry.csv",
        selection_name="v10_2_25_v913_paper_campaign_selection.json",
    ),
    "v10.2.26": AuditedEntrySpec(
        cli_key="v10.2.26",
        audited_module="sharp_front_v10_2_26_audited.py",
        mapping_module="sharp_front_v10_2_26.py",
        registry_name="v10_2_26_v913_weakT_ceramic_registry.csv",
        selection_name="v10_2_26_v913_weakT_ceramic_selection.json",
    ),
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _literal_assignment(source: str, name: str) -> Any:
    tree = ast.parse(source)
    for node in tree.body:
        if isinstance(node, ast.Assign):
            if any(isinstance(target, ast.Name) and target.id == name for target in node.targets):
                return ast.literal_eval(node.value)
    raise ValueError(f"{name} is not a literal assignment")


def _float(row: dict[str, str], name: str) -> float:
    try:
        return float(row[name])
    except KeyError as exc:
        raise ValueError(f"registry row is missing required field {name!r}") from exc


def _int(row: dict[str, str], name: str) -> int:
    return int(round(_float(row, name)))


def _selection_metadata(payload: dict[str, Any], candidate_id: str) -> dict[str, Any] | None:
    pools: list[Any] = []
    for key in (
        "primary_candidates",
        "secondary_candidates",
        "selected_candidates",
        "candidates",
    ):
        value = payload.get(key, [])
        if isinstance(value, list):
            pools.extend(value)
    for item in pools:
        if isinstance(item, dict) and item.get("candidate_id") == candidate_id:
            return dict(item)
    return None


def _verify_audited_entry(audited_path: Path, mapping_stem: str) -> None:
    source = audited_path.read_text()
    required_tokens = (
        f"from . import {mapping_stem} as _entry",
        "AuditedPersistentSiteStateResolvedTipEngine",
        "install_backstress_complementarity_fix",
        "install_physical_front_width",
    )
    missing = [token for token in required_tokens if token not in source]
    if missing:
        raise ValueError(
            f"{audited_path} does not satisfy the audited entry contract; missing={missing}"
        )


def _validate_closure(row: dict[str, str]) -> None:
    zero_fields = (
        "source_recovery_rate_s",
        "retained_recovery_rate_s",
        "recovery_nu0_s",
        "source_refresh_length_um",
        "legacy_source_sites_active",
        "legacy_source_refresh_active",
        "explicit_recovery_active",
    )
    nonzero = {name: _float(row, name) for name in zero_fields if abs(_float(row, name)) > 1.0e-30}
    if nonzero:
        raise ValueError(
            "selected PF row violates the audited persistent-site closure: "
            + json.dumps(nonzero, sort_keys=True)
        )
    if _int(row, "n_slip_channels") != 2:
        raise ValueError("audited FEM/CZM parity requires exactly two reduced slip channels")
    for name in ("reference_source_area_um2", "reference_front_width_um", "source_zone_length_um"):
        if _float(row, name) <= 0.0:
            raise ValueError(f"{name} must be positive")


def _candidate_from_row(row: dict[str, str]) -> PersistentSiteRowV100514:
    return PersistentSiteRowV100514(
        option_key=row["option_key"],
        candidate_id=row["candidate_id"],
        role=row["role"],
        emit_G00_eV=_float(row, "emit_G00_eV"),
        emit_gT_eV_per_K=_float(row, "emit_gT_eV_per_K"),
        peierls_H0_eV=_float(row, "peierls_H0_eV"),
        peierls_activation_entropy_kB=_float(row, "peierls_activation_entropy_kB"),
        taylor_H0_eV=_float(row, "taylor_H0_eV"),
        taylor_activation_entropy_kB=_float(row, "taylor_activation_entropy_kB"),
        taylor_corr_rho_c_m2=_float(row, "taylor_corr_rho_c_m2"),
        rho_source0_m2=_float(row, "rho_source0_m2"),
        source_refresh_length_um_provenance=_float(row, "source_refresh_length_um"),
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
        source_sites_per_system_provenance=_float(row, "source_sites_per_system"),
        encounter_efficiency=_float(row, "encounter_efficiency"),
        retained_recovery_rate_s=_float(row, "retained_recovery_rate_s"),
        c_blunt=_float(row, "c_blunt"),
        recovery_nu0_s=_float(row, "recovery_nu0_s"),
        reference_source_area_um2=_float(row, "reference_source_area_um2"),
        reference_front_width_um=_float(row, "reference_front_width_um"),
        source_zone_length_um=_float(row, "source_zone_length_um"),
        minimum_front_width_um=0.0,
    ).validate()


def load_audited_parameter_option(
    pf_repo_root: str | Path,
    parameter_entry: str,
    parameter_option: str,
) -> tuple[PersistentSiteRowV100514, dict[str, Any]]:
    root = Path(pf_repo_root).expanduser().resolve()
    try:
        spec = ENTRY_SPECS[str(parameter_entry)]
    except KeyError as exc:
        raise KeyError(
            f"unknown audited parameter entry {parameter_entry!r}; allowed={sorted(ENTRY_SPECS)}"
        ) from exc

    package = root / "arrhenius_fracture"
    materials = package / "data" / "materials"
    audited_path = package / spec.audited_module
    mapping_path = package / spec.mapping_module
    registry_path = materials / spec.registry_name
    selection_path = materials / spec.selection_name
    required = (audited_path, mapping_path, registry_path, selection_path)
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing audited PF parameter files: " + ", ".join(missing))

    _verify_audited_entry(audited_path, Path(spec.mapping_module).stem)
    mapping_source = mapping_path.read_text()
    valid_options = _literal_assignment(mapping_source, "VALID_OPTIONS")
    if spec.registry_name not in mapping_source or spec.selection_name not in mapping_source:
        raise ValueError(f"{mapping_path} does not reference the expected registry/selection files")
    if parameter_option not in valid_options:
        raise KeyError(
            f"option {parameter_option!r} is not exposed by {spec.audited_module}; "
            f"allowed={sorted(valid_options)}"
        )

    with registry_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    selected = next((row for row in rows if row.get("option_key") == parameter_option), None)
    if selected is None:
        raise KeyError(f"option {parameter_option!r} is absent from {registry_path}")
    if selected.get("candidate_id") != valid_options[parameter_option]:
        raise ValueError(
            "stable option mapping disagrees with registry candidate: "
            f"mapping={valid_options[parameter_option]!r} registry={selected.get('candidate_id')!r}"
        )
    _validate_closure(selected)
    candidate = _candidate_from_row(selected)

    selection_payload = json.loads(selection_path.read_text())
    metadata = _selection_metadata(selection_payload, candidate.candidate_id)
    audit = {
        "schema": BRIDGE_SCHEMA,
        "PF_branch_required": PF_BRANCH_REQUIRED,
        "PF_repo_root": str(root),
        "parameter_entry": spec.cli_key,
        "audited_model_entry": str(audited_path),
        "stable_option_mapping_module": str(mapping_path),
        "parameter_option": parameter_option,
        "candidate_id": candidate.candidate_id,
        "material_class": selected.get("material_class"),
        "paper_role": selected.get("role"),
        "mechanism_summary": selected.get("mechanism_summary"),
        "validation_status": selected.get("validation_status"),
        "selection_metadata": metadata,
        "registry_path": str(registry_path),
        "selection_path": str(selection_path),
        "audited_entry_sha256": _sha256(audited_path),
        "mapping_module_sha256": _sha256(mapping_path),
        "registry_sha256": _sha256(registry_path),
        "selection_sha256": _sha256(selection_path),
        "selected_registry_row": dict(selected),
        "converted_candidate": asdict(candidate),
        "parameter_values_manually_reconstructed": False,
        "selected_through_audited_entry_option_map": True,
        "persistent_sites": True,
        "finite_source_inventory": False,
        "source_refresh": False,
        "explicit_recovery": False,
    }
    audit["selected_row_sha256"] = hashlib.sha256(
        json.dumps(selected, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    return candidate, audit


__all__ = [
    "BRIDGE_SCHEMA",
    "ENTRY_SPECS",
    "PF_BRANCH_REQUIRED",
    "load_audited_parameter_option",
]
