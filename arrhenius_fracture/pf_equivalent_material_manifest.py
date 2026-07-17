"""PF v10.1.7.1 material import for the kinetic CZM parity branch.

This module is intentionally a loader, not a fitter. The packaged PF manifests
are copied verbatim from PF-fracture-fatigue branch
``v10.1.7.1-final-production-temperature-sweep``.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import csv
import math
from pathlib import Path
from typing import Any, Mapping

import numpy as np

KB_EV_PER_K = 8.617333262145e-5
TREF_K = 481.33
PF_SOURCE = "pf_v10_1_7_1"
LEGACY_SOURCE = "czm_legacy_v9"
VALID_PARAMETER_SOURCES = (PF_SOURCE, LEGACY_SOURCE)


def _finite(row: Mapping[str, Any], key: str) -> float:
    try:
        value = float(row[key])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"material field {key!r} is missing or nonnumeric") from exc
    if not math.isfinite(value):
        raise ValueError(f"material field {key!r} is not finite")
    return value


def normalize_material_class(value: str) -> str:
    key = str(value).strip().lower().replace("-", "").replace("_", "")
    aliases = {
        "ceramic": "ceramic",
        "ceramiclike": "ceramic",
        "weakt": "weakT",
        "fcclike": "weakT",
        "weaktemperature": "weakT",
        "dbtt": "DBTT",
    }
    if key not in aliases:
        raise ValueError(f"unknown material class {value!r}")
    return aliases[key]


@dataclass(frozen=True)
class ExpFloorBarrier:
    G00_eV: float
    gT_eV_per_K: float
    sigc0_Pa: float
    sT_Pa_per_K: float
    alpha: float
    exponent: float
    floor_fraction: float
    attempt_frequency_s: float
    floor_min_eV: float = 1.0e-4
    floor_max_fraction: float = 0.95
    Tref_K: float = TREF_K

    def values_eV(self, stress_Pa: np.ndarray | float, T_K: float) -> np.ndarray:
        sigma = np.maximum(np.asarray(stress_Pa, dtype=float), 0.0)
        dT = float(T_K) - self.Tref_K
        G0 = max(self.G00_eV + self.gT_eV_per_K * dT, 1.0e-12)
        sigc = max(self.sigc0_Pa + self.sT_Pa_per_K * dT, 1.0)
        raw_floor = max(self.floor_min_eV, self.floor_fraction * G0)
        floor = min(self.floor_max_fraction * G0, raw_floor)
        return np.maximum(
            floor + (G0 - floor) * np.exp(
                -max(self.alpha, 0.0)
                * np.power(sigma / sigc, max(self.exponent, 1.0e-12))
            ),
            0.0,
        )

    def rate(self, stress_Pa: np.ndarray | float, T_K: float) -> np.ndarray:
        barrier = self.values_eV(stress_Pa, T_K)
        return self.attempt_frequency_s * np.exp(
            np.clip(-barrier / max(KB_EV_PER_K * float(T_K), 1.0e-30), -700.0, 0.0)
        )


@dataclass(frozen=True)
class TransportBarrier:
    H0_eV: float
    activation_entropy_kB: float
    alpha: float
    exponent: float
    attempt_frequency_s: float


@dataclass(frozen=True)
class PFEquivalentMaterialManifest:
    name: str
    candidate_id: str
    cleavage: ExpFloorBarrier
    emission: ExpFloorBarrier
    peierls: TransportBarrier
    taylor: TransportBarrier
    taylor_corr_rho_c_m2: float
    taylor_corr_scale: float
    source_sites_per_system: float
    encounter_efficiency: float
    retained_recovery_rate_s: float
    source_refresh_length_m: float
    c_blunt: float
    max_K_shield_MPa_sqrt_m: float
    orientation_factors: tuple[float, ...] = (1.0, 1.0)
    parameter_source: str = PF_SOURCE
    source_path: str | None = None

    @classmethod
    def from_row(
        cls,
        row: Mapping[str, Any],
        *,
        parameter_source: str,
        source_path: str | None = None,
    ) -> "PFEquivalentMaterialManifest":
        emission = ExpFloorBarrier(
            G00_eV=_finite(row, "emit_G00_eV"),
            gT_eV_per_K=_finite(row, "emit_gT_eV_per_K"),
            sigc0_Pa=_finite(row, "emit_sigc0_GPa") * 1.0e9,
            sT_Pa_per_K=_finite(row, "emit_sT_GPa_per_K") * 1.0e9,
            alpha=_finite(row, "emit_exp_a"),
            exponent=_finite(row, "emit_exp_n"),
            floor_fraction=_finite(row, "emit_floor_frac"),
            attempt_frequency_s=float(row.get("taylor_nu0_s", 1.0e11) or 1.0e11),
        )
        cleavage = ExpFloorBarrier(
            G00_eV=_finite(row, "cleave_G00_eV"),
            gT_eV_per_K=_finite(row, "cleave_gT_eV_per_K"),
            sigc0_Pa=_finite(row, "cleave_sigc0_GPa") * 1.0e9,
            sT_Pa_per_K=_finite(row, "cleave_sT_GPa_per_K") * 1.0e9,
            alpha=_finite(row, "cleave_exp_a"),
            exponent=_finite(row, "cleave_exp_n"),
            floor_fraction=_finite(row, "cleave_floor_frac"),
            attempt_frequency_s=float(row.get("peierls_nu0_s", 1.0e12) or 1.0e12),
        )
        return cls(
            name=normalize_material_class(str(row.get("target_class", ""))),
            candidate_id=str(row.get("candidate_id", "UNKNOWN")),
            cleavage=cleavage,
            emission=emission,
            peierls=TransportBarrier(
                H0_eV=_finite(row, "peierls_H0_eV"),
                activation_entropy_kB=_finite(row, "peierls_activation_entropy_kB"),
                alpha=_finite(row, "peierls_exp_a"),
                exponent=_finite(row, "peierls_exp_n"),
                attempt_frequency_s=float(row.get("peierls_nu0_s", 1.0e12) or 1.0e12),
            ),
            taylor=TransportBarrier(
                H0_eV=_finite(row, "taylor_H0_eV"),
                activation_entropy_kB=_finite(row, "taylor_activation_entropy_kB"),
                alpha=_finite(row, "taylor_exp_a"),
                exponent=_finite(row, "taylor_exp_n"),
                attempt_frequency_s=float(row.get("taylor_nu0_s", 1.0e11) or 1.0e11),
            ),
            taylor_corr_rho_c_m2=_finite(row, "taylor_corr_rho_c_m2"),
            taylor_corr_scale=_finite(row, "taylor_corr_scale"),
            source_sites_per_system=_finite(row, "source_sites_per_system"),
            encounter_efficiency=_finite(row, "encounter_efficiency"),
            retained_recovery_rate_s=_finite(row, "retained_recovery_rate_s"),
            source_refresh_length_m=_finite(row, "source_refresh_length_um") * 1.0e-6,
            c_blunt=_finite(row, "c_blunt"),
            max_K_shield_MPa_sqrt_m=_finite(row, "max_K_shield_MPa_sqrt_m"),
            orientation_factors=(1.0, 1.0),
            parameter_source=parameter_source,
            source_path=source_path,
        )

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)


def _one_csv_row(path: Path) -> dict[str, str]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 1:
        raise ValueError(f"expected exactly one material row in {path}; found {len(rows)}")
    return rows[0]


def pf_manifest_path(material_class: str) -> Path:
    cls = normalize_material_class(material_class)
    return Path(__file__).resolve().parent / "data" / "pf_v10_1_7_1" / f"{cls}.csv"


def load_material_manifest(
    material_class: str,
    *,
    parameter_source: str = PF_SOURCE,
    legacy_manifest_path: str | Path | None = None,
) -> PFEquivalentMaterialManifest:
    source = str(parameter_source).strip()
    if source not in VALID_PARAMETER_SOURCES:
        raise ValueError(
            f"unknown parameter source {source!r}; expected one of {VALID_PARAMETER_SOURCES}"
        )
    cls = normalize_material_class(material_class)
    if source == PF_SOURCE:
        path = pf_manifest_path(cls)
        if not path.exists():
            raise FileNotFoundError(path)
        return PFEquivalentMaterialManifest.from_row(
            _one_csv_row(path), parameter_source=source, source_path=str(path)
        )

    if legacy_manifest_path is None:
        raise ValueError("czm_legacy_v9 requires legacy_manifest_path")
    from .mpz_parameterization_v911 import load_selected_row

    path = Path(legacy_manifest_path)
    row = load_selected_row(path, cls)
    row = dict(row)
    row.setdefault("max_K_shield_MPa_sqrt_m", 0.0)
    row.setdefault("peierls_nu0_s", 1.0e12)
    row.setdefault("taylor_nu0_s", 1.0e11)
    return PFEquivalentMaterialManifest.from_row(
        row, parameter_source=source, source_path=str(path.resolve())
    )


__all__ = [
    "ExpFloorBarrier",
    "TransportBarrier",
    "PFEquivalentMaterialManifest",
    "PF_SOURCE",
    "LEGACY_SOURCE",
    "VALID_PARAMETER_SOURCES",
    "normalize_material_class",
    "pf_manifest_path",
    "load_material_manifest",
]
