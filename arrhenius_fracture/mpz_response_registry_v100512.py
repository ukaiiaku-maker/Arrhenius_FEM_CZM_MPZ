"""Authoritative v9.11.1 response-option registry for Phase-C FEM/CZM runs.

This module is a loader and fingerprinting layer, not a fitter.  It preserves the
exact selected candidate rows supplied by the completed 0-D/1-D/2-D
parameterization campaign and materializes one-row CSV manifests for the
existing v9.11 constitutive parser.
"""
from __future__ import annotations

from dataclasses import dataclass
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping

from .pf_equivalent_material_manifest import PFEquivalentMaterialManifest

POINT_RELEASE = "10.0.5.12"
PARAMETER_SOURCE = "mpz_v9_11_1"
REGISTRY_FILENAME = "response_registry.json"
PRIMARY_OPTION_KEYS = (
    "ceramic_primary",
    "weakT_primary",
    "dbtt_primary",
    "peak_primary",
)
EXPECTED_CANDIDATE_IDS = {
    "ceramic_primary": "ceramic_restart02_candidate00",
    "weakT_primary": "weakT_restart00_candidate00",
    "dbtt_primary": "DBTT_restart04_candidate03",
    "peak_primary": "DBTT_restart05_candidate61",
}
CANONICAL_CLASSES = {
    "ceramic_primary": "ceramic",
    "weakT_primary": "weakT",
    "dbtt_primary": "DBTT",
    "peak_primary": "DBTT",
}
ALIASES = {
    "ceramic": "ceramic_primary",
    "ceramicprimary": "ceramic_primary",
    "ceramiclike": "ceramic_primary",
    "weakt": "weakT_primary",
    "weaktprimary": "weakT_primary",
    "weaktemperature": "weakT_primary",
    "dbtt": "dbtt_primary",
    "dbttprimary": "dbtt_primary",
    "peak": "peak_primary",
    "peakprimary": "peak_primary",
    "toughnesspeak": "peak_primary",
}

_REQUIRED_NUMERIC_FIELDS = (
    "Tref_K",
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
    "peierls_nu0_s", "taylor_nu0_s",
    "L_pz_um_recommended", "n_bins_recommended",
)


def default_registry_path() -> Path:
    return Path(__file__).resolve().parent / "data" / "mpz_v9_11_1" / REGISTRY_FILENAME


def _default_registry_path() -> Path:
    """Backward-compatible private alias for the public path helper."""
    return default_registry_path()


def normalize_option_key(value: str) -> str:
    raw = str(value).strip()
    if raw in PRIMARY_OPTION_KEYS:
        return raw
    key = raw.lower().replace("-", "").replace("_", "").replace(" ", "")
    if key in ALIASES:
        return ALIASES[key]
    raise ValueError(
        f"unknown Phase-C response option {value!r}; expected one of {PRIMARY_OPTION_KEYS}"
    )


def _finite(row: Mapping[str, Any], key: str) -> float:
    try:
        value = float(row[key])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"registry field {key!r} is missing or nonnumeric") from exc
    if not math.isfinite(value):
        raise ValueError(f"registry field {key!r} is not finite")
    return value


def _fingerprint_payload(row: Mapping[str, Any]) -> dict[str, Any]:
    excluded = {
        "validation_status",
        "mechanism_summary",
        "role",
        "parameter_manifest",
        "parameter_fingerprint_sha256",
    }
    return {key: row[key] for key in sorted(row) if key not in excluded}


def parameter_fingerprint(row: Mapping[str, Any]) -> str:
    return hashlib.sha256(
        json.dumps(
            _fingerprint_payload(row),
            sort_keys=True,
            separators=(",", ":"),
            default=str,
        ).encode("utf-8")
    ).hexdigest()


@dataclass(frozen=True)
class ResponseOptionV100512:
    option_key: str
    candidate_id: str
    canonical_class: str
    mpz_length_um: float
    mpz_n_bins: int
    row: dict[str, Any]
    fingerprint_sha256: str
    registry_path: str

    def manifest_row(self) -> dict[str, Any]:
        row = dict(self.row)
        row.update(
            {
                "option_key": self.option_key,
                "candidate_id": self.candidate_id,
                "target_class": self.canonical_class,
                "selection_role": "primary",
                # v10.0.5 uses raw signed state-derived shielding; no fitted clip.
                "max_K_shield_MPa_sqrt_m": 0.0,
                "parameter_source": PARAMETER_SOURCE,
                "parameter_fingerprint_sha256": self.fingerprint_sha256,
            }
        )
        return row

    def material_manifest(self, source_path: str | None = None) -> PFEquivalentMaterialManifest:
        return PFEquivalentMaterialManifest.from_row(
            self.manifest_row(),
            parameter_source=PARAMETER_SOURCE,
            source_path=source_path,
        )

    def write_selected_csv(self, path: str | Path) -> Path:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        row = self.manifest_row()
        fields = sorted(row)
        with target.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
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
            "mpz_length_um": self.mpz_length_um,
            "mpz_n_bins": self.mpz_n_bins,
            "fingerprint_sha256": self.fingerprint_sha256,
            "registry_path": self.registry_path,
            "mobile_shield_fraction": float(self.row.get("mobile_shield_fraction", 0.0)),
            "source_recovery_rate_s": float(self.row.get("source_recovery_rate_s", 0.0)),
        }


def _validate_primary_row(row: Mapping[str, Any], registry_path: Path) -> ResponseOptionV100512:
    option_key = normalize_option_key(str(row.get("option_key", "")))
    candidate_id = str(row.get("candidate_id", ""))
    expected = EXPECTED_CANDIDATE_IDS[option_key]
    if candidate_id != expected:
        raise RuntimeError(
            f"{option_key} candidate mismatch: registry={candidate_id!r}, expected={expected!r}"
        )
    for key in _REQUIRED_NUMERIC_FIELDS:
        _finite(row, key)
    mobile_fraction = _finite(row, "mobile_shield_fraction")
    source_recovery = _finite(row, "source_recovery_rate_s")
    if abs(mobile_fraction) > 1.0e-15:
        raise RuntimeError(f"{option_key} requires mobile_shield_fraction=0, got {mobile_fraction}")
    if abs(source_recovery) > 1.0e-15:
        raise RuntimeError(f"{option_key} requires source_recovery_rate_s=0, got {source_recovery}")
    mpz_length = _finite(row, "L_pz_um_recommended")
    mpz_bins_float = _finite(row, "n_bins_recommended")
    mpz_bins = int(round(mpz_bins_float))
    if mpz_length <= 0.0 or mpz_bins < 2 or abs(mpz_bins_float - mpz_bins) > 1.0e-9:
        raise RuntimeError(f"invalid MPZ grid for {option_key}: {mpz_length} um, {mpz_bins_float} bins")
    normalized = dict(row)
    normalized["option_key"] = option_key
    normalized["candidate_id"] = candidate_id
    fingerprint = parameter_fingerprint(normalized)
    return ResponseOptionV100512(
        option_key=option_key,
        candidate_id=candidate_id,
        canonical_class=CANONICAL_CLASSES[option_key],
        mpz_length_um=mpz_length,
        mpz_n_bins=mpz_bins,
        row=normalized,
        fingerprint_sha256=fingerprint,
        registry_path=str(registry_path.resolve()),
    )


def load_registry(path: str | Path | None = None) -> dict[str, ResponseOptionV100512]:
    registry_path = Path(path) if path is not None else _default_registry_path()
    if not registry_path.is_file():
        raise FileNotFoundError(registry_path)
    payload = json.loads(registry_path.read_text())
    if not isinstance(payload, list):
        raise ValueError(f"registry must contain a JSON list: {registry_path}")
    rows: dict[str, ResponseOptionV100512] = {}
    for source in payload:
        if not isinstance(source, dict) or "option_key" not in source:
            continue
        raw_key = str(source["option_key"])
        try:
            key = normalize_option_key(raw_key)
        except ValueError:
            continue
        if key not in PRIMARY_OPTION_KEYS:
            continue
        if key in rows:
            raise RuntimeError(f"duplicate primary option {key} in {registry_path}")
        rows[key] = _validate_primary_row(source, registry_path)
    missing = [key for key in PRIMARY_OPTION_KEYS if key not in rows]
    if missing:
        raise RuntimeError(f"registry lacks primary options: {missing}")
    return rows


def load_option(value: str, path: str | Path | None = None) -> ResponseOptionV100512:
    return load_registry(path)[normalize_option_key(value)]


__all__ = [
    "POINT_RELEASE",
    "PARAMETER_SOURCE",
    "PRIMARY_OPTION_KEYS",
    "EXPECTED_CANDIDATE_IDS",
    "ResponseOptionV100512",
    "default_registry_path",
    "normalize_option_key",
    "parameter_fingerprint",
    "load_registry",
    "load_option",
]
