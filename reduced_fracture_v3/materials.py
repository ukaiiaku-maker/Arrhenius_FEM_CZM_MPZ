"""Versioned material ownership and fail-closed transfer audit for OneD V3."""
from __future__ import annotations

from dataclasses import asdict, dataclass
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Iterable, Mapping


MATERIAL_BUNDLE_SCHEMA = "oneD.v3.material-bundle/1"
MATERIAL_MAPPING_SCHEMA = "oneD.v3.material-field-mapping-audit/1"

FRACTURE_ROWS = {
    "Peak": "v913_zeroD_sobol_0242980",
    "DBTT": "v913_zeroD_sobol_0202500",
    "weak-T": "oneD_v2_focused_weak_T_0016",
    "ceramic-like": "oneD_v2_focused_ceramic_like_0018",
}

COMMON_VOID_KINETICS_ROW_ID = "v5.voiding-config.reference-one-void/1"
COMMON_ELASTIC_ROW_ID = "v5.tungsten-plane-strain-E210GPa-nu0p3/1"
COMMON_SITE_POPULATION_ROW_ID = "v5.single-centerline-site-seed3621/1"
COMMON_SPECIMEN_LOADING_ROW_ID = "v5.square-1mm-fixed-opening-G0-companion/1"


@dataclass(frozen=True)
class MaterialBundle:
    """Five independently owned rows needed by a paired transfer case."""

    fracture_material_row_id: str
    void_kinetics_row_id: str
    elastic_row_id: str
    site_population_row_id: str
    specimen_loading_row_id: str
    schema: str = MATERIAL_BUNDLE_SCHEMA

    def __post_init__(self) -> None:
        for name, value in asdict(self).items():
            if name != "schema" and (not isinstance(value, str) or not value):
                raise ValueError(f"{name} must be a nonempty versioned identity")
        if self.fracture_material_row_id == self.void_kinetics_row_id:
            raise ValueError("fracture and void parameter ownership must remain separate")


def pilot_material_bundles() -> dict[str, MaterialBundle]:
    return {
        family: MaterialBundle(
            fracture_material_row_id=row_id,
            void_kinetics_row_id=COMMON_VOID_KINETICS_ROW_ID,
            elastic_row_id=COMMON_ELASTIC_ROW_ID,
            site_population_row_id=COMMON_SITE_POPULATION_ROW_ID,
            specimen_loading_row_id=COMMON_SPECIMEN_LOADING_ROW_ID,
        )
        for family, row_id in FRACTURE_ROWS.items()
    }


_UNITS = {
    "Tref_K": "K",
    "rho_forest_floor_m2": "m^-2",
    "peierls_stress_fraction": "1",
    "taylor_stress_fraction": "1",
    "mobile_shield_fraction": "1",
    "source_recovery_rate_s": "s^-1",
    "L_pz_um_recommended": "um",
    "cleave_G00_eV": "eV",
    "cleave_gT_eV_per_K": "eV K^-1",
    "cleave_sigc0_GPa": "GPa",
    "cleave_sT_GPa_per_K": "GPa K^-1",
    "cleave_exp_a": "1",
    "cleave_exp_n": "1",
    "cleave_floor_frac": "1",
    "emit_G00_eV": "eV",
    "emit_gT_eV_per_K": "eV K^-1",
    "emit_sigc0_GPa": "GPa",
    "emit_sT_GPa_per_K": "GPa K^-1",
    "emit_exp_a": "1",
    "emit_exp_n": "1",
    "emit_floor_frac": "1",
    "peierls_H0_eV": "eV",
    "peierls_activation_entropy_kB": "kB",
    "peierls_exp_a": "1",
    "peierls_exp_n": "1",
    "taylor_H0_eV": "eV",
    "taylor_activation_entropy_kB": "kB",
    "taylor_exp_a": "1",
    "taylor_exp_n": "1",
    "taylor_corr_rho_c_m2": "m^-2",
    "taylor_corr_scale": "1",
    "source_sites_per_system": "1",
    "encounter_efficiency": "1",
    "retained_recovery_rate_s": "s^-1",
    "source_refresh_length_um": "um",
    "c_blunt": "1",
    "peierls_nu0_s": "s^-1",
    "taylor_nu0_s": "s^-1",
    "rho_source0_m2": "m^-2",
    "recovery_nu0_s": "s^-1",
    "recovery_H0_eV": "eV",
    "recovery_activation_entropy_kB": "kB",
    "reference_source_area_um2": "um^2",
    "reference_front_width_um": "um",
    "source_zone_length_um": "um",
}

_METADATA_FIELDS = {
    "material_class",
    "candidate_id",
    "source_or_search_provenance",
    "parameter_decision",
    "reduced_selection_status",
    "same_material_row_both_providers",
}

_ZERO_DISABLED_FIELDS = {
    "mobile_shield_fraction",
    "source_recovery_rate_s",
    "retained_recovery_rate_s",
    "source_refresh_length_um",
    "recovery_nu0_s",
    "recovery_H0_eV",
    "recovery_activation_entropy_kB",
}

_ONE_D_TARGETS = {
    "Tref_K": "CanonicalParameters.values['Tref_K']",
    **{
        name: f"CanonicalParameters.values['{name}']"
        for name in _UNITS
    },
}


def _two_d_structural_target(field: str) -> str | None:
    if field == "Tref_K":
        return "MaterialManifest.cleavage/emission.Tref_K"
    if field.startswith("cleave_"):
        return "FrontEngine.cb." + field.removeprefix("cleave_")
    if field.startswith("emit_"):
        return "FrontEngine.eb." + field.removeprefix("emit_")
    if field.startswith("peierls_"):
        return "MaterialManifest.peierls." + field.removeprefix("peierls_")
    if field.startswith("taylor_"):
        return "MaterialManifest.taylor_or_correlation." + field.removeprefix("taylor_")
    if field == "c_blunt":
        return "FrontEngine.f.c_blunt"
    if field == "L_pz_um_recommended":
        return "FrontEngine.f.L_pz"
    return None


def _conversion(field: str) -> str:
    if field.endswith("_GPa") or "_GPa_per_K" in field:
        return "multiply by 1e9 to Pa"
    if field.endswith("_um_recommended") or field.endswith("_length_um") or field.endswith("_width_um"):
        return "multiply by 1e-6 to m"
    if field.endswith("_area_um2"):
        return "multiply by 1e-12 to m^2"
    return "identity"


def _full_precision(value: str) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    if not math.isfinite(number):
        raise ValueError("material rows must contain finite numeric values")
    return str(value)


def load_exact_fracture_rows(path: str | Path) -> dict[str, dict[str, str]]:
    source = Path(path)
    with source.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    by_id = {str(row["candidate_id"]): row for row in rows}
    expected = set(FRACTURE_ROWS.values())
    if not expected.issubset(by_id):
        raise ValueError(f"registry is missing exact fracture rows: {sorted(expected - set(by_id))}")
    selected = {family: by_id[row_id] for family, row_id in FRACTURE_ROWS.items()}
    for family, row in selected.items():
        if row.get("material_class") != family:
            raise ValueError(f"material family identity mismatch for {family}")
    return selected


def fracture_rows_sha256(rows: Mapping[str, Mapping[str, str]]) -> str:
    payload = {
        family: {key: str(value) for key, value in row.items()}
        for family, row in sorted(rows.items())
    }
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def material_field_mapping_audit(
    rows: Mapping[str, Mapping[str, str]],
    *,
    source_registry: str,
    v5_driver_identity: str,
) -> dict[str, Any]:
    """Audit every source column and reject active unbound V5 fields.

    The pinned V5 driver creates default ``FrontEngine`` objects internally and
    has no fracture-row argument. Structural destinations elsewhere in the 2-D
    repository therefore do not constitute a runtime binding.
    """

    records: list[dict[str, Any]] = []
    unsupported_active: list[dict[str, str]] = []
    for family in FRACTURE_ROWS:
        row = rows[family]
        for field, raw_value in row.items():
            one_d_target = None if field in _METADATA_FIELDS else _ONE_D_TARGETS.get(field)
            two_d_target = None if field in _METADATA_FIELDS else _two_d_structural_target(field)
            inactive = field in _METADATA_FIELDS or (
                field in _ZERO_DISABLED_FIELDS and float(raw_value) == 0.0
            )
            if inactive:
                classification = "inactive"
                reason = "identity/provenance field" if field in _METADATA_FIELDS else "zero disables this pathway in every selected row"
            else:
                # Registry storage and a structural class member are insufficient:
                # both trajectories need an exact runtime consumer.
                classification = "unsupported"
                reason = (
                    "pinned V5 one-void runtime exposes no fracture-material-row binding; "
                    "default FrontEngine/material state would be substituted"
                )
                unsupported_active.append({"material_class": family, "source_field": field})
            records.append({
                "material_class": family,
                "fracture_material_row_id": row["candidate_id"],
                "source_field": field,
                "units": _UNITS.get(field, "string"),
                "full_precision_value": _full_precision(raw_value),
                "one_d_target": one_d_target,
                "two_d_target": two_d_target,
                "conversion": "identity" if field in _METADATA_FIELDS else _conversion(field),
                "classification": classification,
                "source_activity": "inactive" if inactive else "active",
                "reason": reason,
            })
    return {
        "schema": MATERIAL_MAPPING_SCHEMA,
        "source_registry": source_registry,
        "source_rows_sha256": fracture_rows_sha256(rows),
        "v5_driver_identity": v5_driver_identity,
        "records": records,
        "summary": {
            "material_families": len(rows),
            "records": len(records),
            "unsupported_active_records": len(unsupported_active),
            "unsupported_active_fields": sorted({row["source_field"] for row in unsupported_active}),
            "all_active_fields_mapped": not unsupported_active,
            "gate": (
                "PASS_EXACT_RUNTIME_BINDING"
                if not unsupported_active
                else "BLOCKED_UNMAPPED_ACTIVE_FIELD"
            ),
        },
    }


def require_complete_material_mapping(audit: Mapping[str, Any]) -> None:
    if not bool(audit["summary"]["all_active_fields_mapped"]):
        fields = ", ".join(audit["summary"]["unsupported_active_fields"])
        raise RuntimeError("BLOCKED_UNMAPPED_ACTIVE_FIELD: " + fields)


def paired_case_ledger(
    anchors_by_family: Mapping[str, Iterable[float | None]],
    *,
    gate: str,
) -> list[dict[str, Any]]:
    rows = []
    for family in FRACTURE_ROWS:
        anchors = tuple(anchors_by_family[family])
        if len(anchors) != 3:
            raise ValueError("each material family requires exactly three anchor slots")
        for label, temperature in zip(("low", "feature", "high"), anchors):
            rows.append({
                "material_class": family,
                "fracture_material_row_id": FRACTURE_ROWS[family],
                "temperature_anchor": label,
                "temperature_K": temperature,
                "no_void_1d_status": "NOT_RUN_M2_GATE_BLOCKED",
                "no_void_2d_status": "NOT_RUN_M2_GATE_BLOCKED",
                "void_1d_status": "NOT_RUN_M2_GATE_BLOCKED",
                "void_2d_status": "NOT_RUN_M2_GATE_BLOCKED",
                "delta_observable_status": "UNAVAILABLE_NO_PAIRED_RUN",
                "material_mapping_gate": gate,
            })
    return rows


__all__ = [
    "COMMON_ELASTIC_ROW_ID",
    "COMMON_SITE_POPULATION_ROW_ID",
    "COMMON_SPECIMEN_LOADING_ROW_ID",
    "COMMON_VOID_KINETICS_ROW_ID",
    "FRACTURE_ROWS",
    "MATERIAL_BUNDLE_SCHEMA",
    "MATERIAL_MAPPING_SCHEMA",
    "MaterialBundle",
    "fracture_rows_sha256",
    "load_exact_fracture_rows",
    "material_field_mapping_audit",
    "paired_case_ledger",
    "pilot_material_bundles",
    "require_complete_material_mapping",
]
