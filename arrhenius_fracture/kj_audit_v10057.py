"""v10.0.5.7 fixed-grip K-reference utilities.

The FEM specimen is driven by symmetric prescribed displacement, not uniform
remote traction.  The legacy single-edge-tension polynomial is retained as a
secondary diagnostic, while the publication gate uses the fixed-grip geometry
factor for the exact 2 x 4 mm, a=0.5 mm benchmark geometry.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from .kj_audit_v10056 import (
    SpecimenGeometryV10056,
    build_kj_audit_row as _legacy_build_row,
    select_contour_plateau as _legacy_select_plateau,
)

POINT_RELEASE = "10.0.5.7"
REFERENCE_SCHEMA = "fixed_grip_edge_crack_reference_v10057"


@dataclass(frozen=True)
class FixedGripReferenceV10057:
    width_m: float = 2.0e-3
    height_m: float = 4.0e-3
    initial_crack_m: float = 0.5e-3
    geometry_factor_Y: float = 1.2003
    relative_geometry_tolerance: float = 1.0e-10
    provenance: str = (
        "fixed-grip elastic FEM convergence benchmark; regenerate with "
        "scripts/generate_fixed_grip_reference_v10057.py"
    )

    def validate_geometry(
        self, geometry: SpecimenGeometryV10056
    ) -> "FixedGripReferenceV10057":
        geometry.validate()
        expected = np.asarray(
            [self.width_m, self.height_m, self.initial_crack_m], dtype=float
        )
        supplied = np.asarray(
            [geometry.width_m, geometry.height_m, geometry.initial_crack_m],
            dtype=float,
        )
        if not np.allclose(
            supplied,
            expected,
            rtol=self.relative_geometry_tolerance,
            atol=1.0e-15,
        ):
            raise ValueError(
                "fixed-grip reference is geometry specific; regenerate it for "
                f"width/height/crack={supplied.tolist()}"
            )
        if self.geometry_factor_Y <= 0.0 or not math.isfinite(self.geometry_factor_Y):
            raise ValueError("fixed-grip geometry factor must be positive and finite")
        return self

    def to_dict(self) -> dict[str, Any]:
        return {"schema": REFERENCE_SCHEMA, **asdict(self)}


DEFAULT_FIXED_GRIP_REFERENCE = FixedGripReferenceV10057()


def load_fixed_grip_reference(
    path: str | Path | None = None,
) -> FixedGripReferenceV10057:
    if path is None:
        return DEFAULT_FIXED_GRIP_REFERENCE
    payload = json.loads(Path(path).expanduser().resolve().read_text())
    if payload.get("schema") != REFERENCE_SCHEMA:
        raise ValueError(f"fixed-grip reference schema must be {REFERENCE_SCHEMA}")
    if not bool(payload.get("convergence_passed", False)):
        raise ValueError("fixed-grip reference has not passed convergence")
    return FixedGripReferenceV10057(
        width_m=float(payload["width_m"]),
        height_m=float(payload["height_m"]),
        initial_crack_m=float(payload["initial_crack_m"]),
        geometry_factor_Y=float(payload["geometry_factor_Y"]),
        relative_geometry_tolerance=float(
            payload.get("relative_geometry_tolerance", 1.0e-10)
        ),
        provenance=str(payload.get("provenance", "generated fixed-grip reference")),
    )


def fixed_grip_reference_K_Pa_sqrt_m(
    sigma_gross_Pa: float,
    geometry: SpecimenGeometryV10056,
    reference: FixedGripReferenceV10057 = DEFAULT_FIXED_GRIP_REFERENCE,
) -> float:
    reference.validate_geometry(geometry)
    return (
        max(float(sigma_gross_Pa), 0.0)
        * math.sqrt(math.pi * geometry.initial_crack_m)
        * float(reference.geometry_factor_Y)
    )


def build_kj_audit_row(
    *,
    Ftop_N_per_thickness: float,
    KJ_Pa_sqrt_m: float,
    outer_radius_m: float,
    geometry: SpecimenGeometryV10056,
    n_active_elements: int | None = None,
    safety_fraction: float = 0.80,
    fixed_grip_reference: FixedGripReferenceV10057 = DEFAULT_FIXED_GRIP_REFERENCE,
) -> dict[str, Any]:
    row = _legacy_build_row(
        Ftop_N_per_thickness=Ftop_N_per_thickness,
        KJ_Pa_sqrt_m=KJ_Pa_sqrt_m,
        outer_radius_m=outer_radius_m,
        geometry=geometry,
        n_active_elements=n_active_elements,
        safety_fraction=safety_fraction,
    )
    legacy_K = float(row["K_LEFM_gross_MPa_sqrt_m"])
    legacy_ratio = float(row["KJ_over_K_LEFM_gross"])
    legacy_Y = float(row["edge_geometry_factor_Y"])
    sigma = float(row["sigma_gross_MPa"]) * 1.0e6
    fixed_K = fixed_grip_reference_K_Pa_sqrt_m(
        sigma, geometry, fixed_grip_reference
    )
    KJ = max(float(KJ_Pa_sqrt_m), 0.0)
    row.update(
        {
            "K_reference_boundary_condition": "symmetric_fixed_grip_displacement",
            "fixed_grip_reference_schema": REFERENCE_SCHEMA,
            "fixed_grip_geometry_factor_Y": float(
                fixed_grip_reference.geometry_factor_Y
            ),
            "fixed_grip_reference_provenance": fixed_grip_reference.provenance,
            "K_fixed_grip_reference_MPa_sqrt_m": fixed_K / 1.0e6,
            "KJ_over_K_fixed_grip_reference": KJ / fixed_K
            if fixed_K > 0.0
            else math.nan,
            "uniform_tension_edge_geometry_factor_Y": legacy_Y,
            "K_uniform_tension_edge_MPa_sqrt_m": legacy_K,
            "KJ_over_K_uniform_tension_edge": legacy_ratio,
            # Compatibility keys consumed by v10.0.5.6 runners now point to the
            # boundary-condition-consistent primary reference.
            "edge_geometry_factor_Y": float(fixed_grip_reference.geometry_factor_Y),
            "K_LEFM_gross_MPa_sqrt_m": fixed_K / 1.0e6,
            "KJ_over_K_LEFM_gross": KJ / fixed_K if fixed_K > 0.0 else math.nan,
        }
    )
    return row


def select_contour_plateau(
    rows: Sequence[Mapping[str, Any]],
    *,
    relative_tolerance: float = 0.10,
    minimum_points: int = 3,
    minimum_active_elements: int = 12,
    fixed_grip_ratio_min: float = 0.85,
    fixed_grip_ratio_max: float = 1.15,
) -> dict[str, Any]:
    selection = _legacy_select_plateau(
        rows,
        relative_tolerance=relative_tolerance,
        minimum_points=minimum_points,
        minimum_active_elements=minimum_active_elements,
    )
    if selection.get("status") != "plateau_selected":
        selection["fixed_grip_reference_gate_passed"] = False
        return selection
    radii = set(float(value) for value in selection["plateau_outer_radii_m"])
    selected_rows = [
        dict(row) for row in rows if float(row.get("outer_radius_m", -1.0)) in radii
    ]
    ratios = [
        float(row.get("KJ_over_K_fixed_grip_reference", math.nan))
        for row in selected_rows
    ]
    finite = [value for value in ratios if math.isfinite(value)]
    passed = bool(
        len(finite) == len(selected_rows)
        and finite
        and min(finite) >= float(fixed_grip_ratio_min)
        and max(finite) <= float(fixed_grip_ratio_max)
    )
    selection.update(
        {
            "reference_boundary_condition": "symmetric_fixed_grip_displacement",
            "fixed_grip_reference_gate_passed": passed,
            "fixed_grip_ratio_acceptance": [
                float(fixed_grip_ratio_min),
                float(fixed_grip_ratio_max),
            ],
            "plateau_KJ_over_K_fixed_grip_min": min(finite)
            if finite
            else math.nan,
            "plateau_KJ_over_K_fixed_grip_max": max(finite)
            if finite
            else math.nan,
            "status": "plateau_selected"
            if passed
            else "plateau_rejected_fixed_grip_mismatch",
        }
    )
    return selection


__all__ = [
    "POINT_RELEASE",
    "REFERENCE_SCHEMA",
    "FixedGripReferenceV10057",
    "DEFAULT_FIXED_GRIP_REFERENCE",
    "load_fixed_grip_reference",
    "fixed_grip_reference_K_Pa_sqrt_m",
    "build_kj_audit_row",
    "select_contour_plateau",
]
