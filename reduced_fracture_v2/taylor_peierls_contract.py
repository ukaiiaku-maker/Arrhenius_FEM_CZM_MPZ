"""Fail-closed material contract for the Taylor/Peierls-only R search."""
from __future__ import annotations

import hashlib
import json
import struct
from typing import Mapping

from arrhenius_fracture.emergent_gnd_contract_v913 import (
    ACTIVE_CANDIDATE_PARAMETER_FIELDS,
)


CLEAVAGE_FIELDS = (
    "cleave_G00_eV", "cleave_gT_eV_per_K", "cleave_sigc0_GPa",
    "cleave_sT_GPa_per_K", "cleave_exp_a", "cleave_exp_n",
    "cleave_floor_frac",
)
EMISSION_FIELDS = (
    "emit_G00_eV", "emit_gT_eV_per_K", "emit_sigc0_GPa",
    "emit_sT_GPa_per_K", "emit_exp_a", "emit_exp_n", "emit_floor_frac",
)
PEIERLS_FIELDS = (
    "peierls_H0_eV", "peierls_activation_entropy_kB",
    "peierls_exp_a", "peierls_exp_n",
)
TAYLOR_FIELDS = (
    "taylor_H0_eV", "taylor_activation_entropy_kB",
    "taylor_exp_a", "taylor_exp_n", "taylor_corr_rho_c_m2",
    "taylor_corr_scale",
)
SEARCH_WHITELIST = PEIERLS_FIELDS + TAYLOR_FIELDS
FIXED_ATTEMPT_FREQUENCIES = ("peierls_nu0_s", "taylor_nu0_s")
FROZEN_ACTIVE_FIELDS = tuple(
    field for field in ACTIVE_CANDIDATE_PARAMETER_FIELDS
    if field not in SEARCH_WHITELIST
)
SOURCE_FIXED_PROJECTION_FACTORS = {
    "peierls_stress_fraction": 1.0 / (3.0 ** 0.5),
    "taylor_stress_fraction": 1.0 / (3.0 ** 0.5),
}


def _hex_payload(row: Mapping[str, object], fields: tuple[str, ...]) -> str:
    return json.dumps(
        {field: float(row[field]).hex() for field in fields},
        sort_keys=True,
        separators=(",", ":"),
    )


def coordinate_hash(row: Mapping[str, object], fields: tuple[str, ...]) -> str:
    return hashlib.sha256(_hex_payload(row, fields).encode()).hexdigest()


def full_material_json(row: Mapping[str, object]) -> str:
    return json.dumps(
        {field: float(row[field]) for field in ACTIVE_CANDIDATE_PARAMETER_FIELDS},
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )


def full_material_hash(row: Mapping[str, object]) -> str:
    return hashlib.sha256(full_material_json(row).encode()).hexdigest()


def taylor_peierls_hash(row: Mapping[str, object]) -> str:
    return coordinate_hash(row, SEARCH_WHITELIST)


def barrier_hashes(row: Mapping[str, object]) -> dict[str, str]:
    return {
        "cleavage_barrier_sha256": coordinate_hash(row, CLEAVAGE_FIELDS),
        "emission_barrier_sha256": coordinate_hash(row, EMISSION_FIELDS),
    }


def bit_identical(left: object, right: object) -> bool:
    return struct.pack(">d", float(left)) == struct.pack(">d", float(right))


def validate_candidate(candidate: Mapping[str, object], control: Mapping[str, object]) -> None:
    changed = [
        field for field in FROZEN_ACTIVE_FIELDS
        if not bit_identical(candidate[field], control[field])
    ]
    if changed:
        raise ValueError(f"non-whitelisted material coordinates changed: {changed}")
    for fields, label in ((CLEAVAGE_FIELDS, "cleavage"), (EMISSION_FIELDS, "emission")):
        if coordinate_hash(candidate, fields) != coordinate_hash(control, fields):
            raise ValueError(f"{label} barrier hash differs from class control")


__all__ = [
    "CLEAVAGE_FIELDS", "EMISSION_FIELDS", "PEIERLS_FIELDS", "TAYLOR_FIELDS",
    "SEARCH_WHITELIST", "FIXED_ATTEMPT_FREQUENCIES", "FROZEN_ACTIVE_FIELDS",
    "SOURCE_FIXED_PROJECTION_FACTORS", "barrier_hashes", "bit_identical",
    "coordinate_hash", "full_material_hash", "full_material_json",
    "taylor_peierls_hash", "validate_candidate",
]
