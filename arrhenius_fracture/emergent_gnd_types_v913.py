"""Constitutive types for the v9.13 persistent-site 1-D transfer."""
from __future__ import annotations

from dataclasses import dataclass
import math
import numpy as np

from .emergent_gnd_types_v912 import *  # noqa: F401,F403
from .emergent_gnd_types_v912 import (
    CandidateParameters as _CandidateParametersV912,
    CommonPhysics as _CommonPhysicsV912,
)


@dataclass(frozen=True)
class CandidateParameters(_CandidateParametersV912):
    """v9.12 candidate kinetics plus the transferred crack-tip blunting factor."""

    c_blunt: float = 0.0


@dataclass(frozen=True)
class CommonPhysics(_CommonPhysicsV912):
    """Fixed persistent-site geometry and normalization shared with v10.2.21."""

    reference_source_area_m2: float = 25.0e-12
    reference_front_width_m: float = 10.0e-6
    reference_density_m2: float = 5.0e12
    blunting_length_m: float = 0.5e-6
    blunting_slip_fraction: float = 1.0
    persistent_backstress_scale: float = 1.0
    minimum_front_width_m: float = 0.0
    maximum_front_width_m: float = 0.0
    activation_to_line_content_per_system: tuple[float, ...] = (1.0, 1.0)
    implicit_tolerance: float = 1.0e-10
    implicit_max_iterations: int = 96

    def validate(self) -> None:
        super().validate()
        positive_values = (
            (self.reference_source_area_m2, "reference_source_area_m2"),
            (self.reference_front_width_m, "reference_front_width_m"),
            (self.reference_density_m2, "reference_density_m2"),
            (self.blunting_length_m, "blunting_length_m"),
            (self.implicit_tolerance, "implicit_tolerance"),
        )
        for value, name in positive_values:
            if not math.isfinite(float(value)) or float(value) <= 0.0:
                raise ValueError(f"{name} must be positive and finite")
        if self.blunting_slip_fraction < 0.0:
            raise ValueError("blunting_slip_fraction must be nonnegative")
        if self.persistent_backstress_scale <= 0.0:
            raise ValueError("persistent_backstress_scale must be positive")
        if self.minimum_front_width_m < 0.0:
            raise ValueError("minimum_front_width_m must be nonnegative")
        if self.maximum_front_width_m < 0.0:
            raise ValueError("maximum_front_width_m must be nonnegative")
        if int(self.implicit_max_iterations) < 8:
            raise ValueError("implicit_max_iterations must be at least 8")
        conversion = np.asarray(
            self.activation_to_line_content_per_system,
            dtype=float,
        )
        if conversion.shape != (self.n_systems,):
            raise ValueError(
                "activation_to_line_content_per_system length must equal n_systems"
            )
        if np.any(~np.isfinite(conversion)) or np.any(conversion <= 0.0):
            raise ValueError(
                "activation_to_line_content_per_system must be positive and finite"
            )


__all__ = [name for name in globals() if not name.startswith("_")]
