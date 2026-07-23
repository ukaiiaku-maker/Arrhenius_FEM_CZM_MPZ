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
    # Along-front/out-of-plane correlation length.  This is deliberately
    # independent of the ahead-of-tip MPZ spacing ``dx``.
    minimum_front_width_m: float = 10.0e-9
    maximum_front_width_m: float = 0.0
    # Keep the existing event-translation contract by default.  The coupled
    # moving-tip integrator is available as an explicit audit because it
    # materially changes the signed 1-D shielding response.
    coupled_moving_tip_enabled: bool = False
    # Maximum moving-frame translation per coupled kinetic substep, expressed
    # as a fraction of the ahead-of-tip MPZ cell width.
    moving_tip_cfl: float = 0.25
    activation_to_line_content_per_system: tuple[float, ...] = (1.0, 1.0)
    # Geometric mobile/obstacle encounter multiplier used by the 2-D MPZ
    # transport law: k_enc = eta * v_P * sqrt(rho_forest).
    encounter_efficiency: float = 1.0
    # Maximum forest-spacing/Burgers-vector amplification of the local Taylor
    # obstacle stress.  ``inf`` preserves the original v9.13 expression; the
    # v10.2.22 2-D runs use the shared physical value 20.
    taylor_phi_max: float = float("inf")
    # Scale only spatial Peierls advection/escape.  Encounter storage continues
    # to use the unscaled Peierls velocity.  This separation is required for
    # the v10.2.22 validated-scalar transport mode, whose channel-resolved
    # transport velocity and cumulative escape are effectively zero.
    mobile_transport_velocity_scale: float = 1.0
    # Optional reduced representation of a measured 2-D slip-projection path.
    # The factor row whose breakpoint is the greatest value not exceeding the
    # current cumulative crack extension replaces the constant
    # ``emission_schmid_factors``.  An empty table preserves the legacy law.
    emission_geometry_extension_m: tuple[float, ...] = ()
    emission_geometry_factors: tuple[tuple[float, ...], ...] = ()
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
        if (
            not math.isfinite(float(self.encounter_efficiency))
            or self.encounter_efficiency < 0.0
        ):
            raise ValueError("encounter_efficiency must be finite and nonnegative")
        if math.isnan(float(self.taylor_phi_max)) or self.taylor_phi_max < 1.0:
            raise ValueError("taylor_phi_max must be at least one")
        if (
            not math.isfinite(float(self.mobile_transport_velocity_scale))
            or self.mobile_transport_velocity_scale < 0.0
        ):
            raise ValueError(
                "mobile_transport_velocity_scale must be finite and nonnegative"
            )
        if (
            not math.isfinite(float(self.minimum_front_width_m))
            or self.minimum_front_width_m <= 0.0
        ):
            raise ValueError(
                "minimum_front_width_m must be a positive physical length"
            )
        if (
            not math.isfinite(float(self.maximum_front_width_m))
            or self.maximum_front_width_m < 0.0
        ):
            raise ValueError("maximum_front_width_m must be finite and nonnegative")
        if (
            self.maximum_front_width_m > 0.0
            and self.maximum_front_width_m < self.minimum_front_width_m
        ):
            raise ValueError(
                "maximum_front_width_m must not be smaller than "
                "minimum_front_width_m"
            )
        if (
            not math.isfinite(float(self.moving_tip_cfl))
            or not 0.0 < self.moving_tip_cfl <= 1.0
        ):
            raise ValueError("moving_tip_cfl must lie in (0, 1]")
        if not isinstance(self.coupled_moving_tip_enabled, (bool, np.bool_)):
            raise ValueError("coupled_moving_tip_enabled must be boolean")
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
        geometry_extension = np.asarray(
            self.emission_geometry_extension_m,
            dtype=float,
        )
        geometry_factors = np.asarray(
            self.emission_geometry_factors,
            dtype=float,
        )
        if geometry_extension.size == 0:
            if geometry_factors.size != 0:
                raise ValueError(
                    "emission_geometry_factors requires extension breakpoints"
                )
        else:
            if geometry_extension.ndim != 1:
                raise ValueError(
                    "emission_geometry_extension_m must be one-dimensional"
                )
            if (
                np.any(~np.isfinite(geometry_extension))
                or np.any(geometry_extension < 0.0)
                or np.any(np.diff(geometry_extension) <= 0.0)
            ):
                raise ValueError(
                    "emission geometry breakpoints must be finite, nonnegative, "
                    "and strictly increasing"
                )
            if geometry_factors.shape != (
                geometry_extension.size,
                self.n_systems,
            ):
                raise ValueError(
                    "emission_geometry_factors must have one n_systems row "
                    "per extension breakpoint"
                )
            if (
                np.any(~np.isfinite(geometry_factors))
                or np.any(geometry_factors < 0.0)
            ):
                raise ValueError(
                    "emission geometry factors must be finite and nonnegative"
                )


__all__ = [name for name in globals() if not name.startswith("_")]
