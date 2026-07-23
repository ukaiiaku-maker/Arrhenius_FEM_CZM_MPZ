"""Immutable candidate-parameter contract for the v9.13 persistent model.

The v10.2.22 transfer disables finite source refresh and explicit recovery.
Those legacy registry coordinates may remain in historical candidate CSVs for
provenance, but they do not define the persistent-site constitutive state and
must not be used as surrogate features.  This module provides one canonical
list of the candidate values that do remain active.
"""

from __future__ import annotations

import hashlib
import json
import math
from typing import Any, Mapping, Sequence


ACTIVE_CANDIDATE_PARAMETER_FIELDS = (
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

ACTIVE_CANDIDATE_PARAMETER_DEFAULTS = {
    "Tref_K": 481.33,
    "peierls_nu0_s": 1.0e12,
    "taylor_nu0_s": 1.0e11,
}

PERSISTENT_INACTIVE_REGISTRY_FIELDS = (
    "source_refresh_length_um",
    "recovery_nu0_s",
    "recovery_H0_eV",
    "recovery_activation_entropy_kB",
)

# Backward-compatible name used by the target and calibration scripts.
CANDIDATE_PARAMETER_FIELDS = ACTIVE_CANDIDATE_PARAMETER_FIELDS


def effective_candidate_parameters(row: Mapping[str, Any]) -> dict[str, float]:
    """Return the finite values actually consumed by the v9.13 model."""
    values: dict[str, float] = {}
    missing: list[str] = []
    for field in ACTIVE_CANDIDATE_PARAMETER_FIELDS:
        raw = row.get(field)
        if raw in (None, ""):
            if field not in ACTIVE_CANDIDATE_PARAMETER_DEFAULTS:
                missing.append(field)
                continue
            raw = ACTIVE_CANDIDATE_PARAMETER_DEFAULTS[field]
        value = float(raw)
        if not math.isfinite(value):
            raise ValueError(
                f"candidate {row.get('candidate_id', '<unknown>')!r} has "
                f"nonfinite {field}={raw!r}"
            )
        values[field] = value
    if missing:
        raise KeyError(
            f"candidate {row.get('candidate_id', '<unknown>')!r} is missing "
            f"active fields: {missing}"
        )
    return values


def candidate_parameter_payload(row: Mapping[str, Any]) -> dict[str, Any]:
    candidate_id = row.get("candidate_id")
    if candidate_id in (None, ""):
        raise KeyError("candidate row is missing candidate_id")
    return {
        "candidate_id": str(candidate_id),
        **effective_candidate_parameters(row),
    }


def candidate_parameter_fingerprint(
    rows: Sequence[Mapping[str, Any]],
) -> str:
    """Hash candidate IDs and normalized active values in stable order."""
    payload = [
        candidate_parameter_payload(row)
        for row in sorted(rows, key=lambda item: str(item["candidate_id"]))
    ]
    text = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def candidate_feature_record(row: Mapping[str, Any]) -> dict[str, float]:
    """Return the complete active feature vector used by the surrogate."""
    return {
        f"x_raw__{field}": value
        for field, value in effective_candidate_parameters(row).items()
    }


__all__ = [
    "ACTIVE_CANDIDATE_PARAMETER_DEFAULTS",
    "ACTIVE_CANDIDATE_PARAMETER_FIELDS",
    "CANDIDATE_PARAMETER_FIELDS",
    "PERSISTENT_INACTIVE_REGISTRY_FIELDS",
    "candidate_feature_record",
    "candidate_parameter_fingerprint",
    "candidate_parameter_payload",
    "effective_candidate_parameters",
]
