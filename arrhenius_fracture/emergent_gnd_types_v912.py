"""Constitutive types for the v9.12 emergent-GND campaign."""
from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Any

import numpy as np

KB_J_PER_K = 1.380649e-23
EV_TO_J = 1.602176634e-19
KB_EV_PER_K = KB_J_PER_K / EV_TO_J


def positive(value: float, floor: float = 0.0) -> float:
    return max(float(value), float(floor))


@dataclass(frozen=True)
class ExpFloorSurface:
    G00_eV: float
    gT_eV_per_K: float
    sigc0_Pa: float
    sT_Pa_per_K: float
    exp_a: float
    exp_n: float
    floor_fraction: float
    Tref_K: float = 481.33
    floor_min_eV: float = 1.0e-4
    floor_max_fraction: float = 0.95

    def zero_stress_eV(self, T_K: float) -> float:
        return positive(
            self.G00_eV + self.gT_eV_per_K * (float(T_K) - self.Tref_K),
            1.0e-12,
        )

    def characteristic_stress_Pa(self, T_K: float) -> float:
        return positive(
            self.sigc0_Pa + self.sT_Pa_per_K * (float(T_K) - self.Tref_K),
            1.0,
        )

    def barrier_eV(self, stress_Pa: np.ndarray | float, T_K: float) -> np.ndarray:
        sigma = np.maximum(np.asarray(stress_Pa, dtype=float), 0.0)
        G0 = self.zero_stress_eV(T_K)
        floor = min(
            self.floor_max_fraction * G0,
            max(self.floor_min_eV, self.floor_fraction * G0),
        )
        x = sigma / self.characteristic_stress_Pa(T_K)
        return np.maximum(
            floor
            + (G0 - floor)
            * np.exp(
                -positive(self.exp_a)
                * np.power(x, positive(self.exp_n, 1.0e-9))
            ),
            0.0,
        )


@dataclass(frozen=True)
class PTMechanism:
    H0_eV: float
    activation_entropy_kB: float
    exp_a: float
    exp_n: float
    nu0_s: float
    stress_fraction: float = 1.0 / math.sqrt(3.0)
    stress_scale: float = 1.0

    def surface(self, parent: ExpFloorSurface) -> ExpFloorSurface:
        ratio = positive(self.H0_eV) / positive(parent.G00_eV, 1.0e-30)
        return ExpFloorSurface(
            G00_eV=self.H0_eV,
            gT_eV_per_K=-self.activation_entropy_kB * KB_EV_PER_K,
            sigc0_Pa=self.stress_scale * parent.sigc0_Pa,
            sT_Pa_per_K=self.stress_scale * parent.sT_Pa_per_K,
            exp_a=self.exp_a,
            exp_n=self.exp_n,
            floor_fraction=parent.floor_fraction,
            Tref_K=parent.Tref_K,
            floor_min_eV=parent.floor_min_eV * max(ratio, 1.0e-12),
            floor_max_fraction=parent.floor_max_fraction,
        )


@dataclass(frozen=True)
class CandidateParameters:
    candidate_id: str
    cleavage: ExpFloorSurface
    emission: ExpFloorSurface
    peierls: PTMechanism
    taylor: PTMechanism
    rho_source0_m2: float
    source_refresh_length_m: float
    taylor_corr_rho_c_m2: float
    taylor_corr_scale: float
    recovery_nu0_s: float = 0.0
    recovery_H0_eV: float = 0.0
    recovery_activation_entropy_kB: float = 0.0


@dataclass(frozen=True)
class CommonPhysics:
    G_Pa: float = 160.15625e9
    nu: float = 0.28
    b_m: float = 2.74e-10
    r0_m: float = 1.0e-6
    mpz_length_m: float = 50.0e-6
    active_strip_width_m: float = 10.0e-6
    n_bins: int = 80
    n_systems: int = 2
    source_zone_length_m: float = 2.0e-6
    rho_forest_floor_m2: float = 5.0e12
    mean_free_path_coefficient: float = 1.0
    jump_fraction_of_forest_spacing: float = 1.0
    annihilation_capture_radius_b: float = 5.0
    core_regularization_b: float = 5.0
    cleavage_nu0_s: float = 1.0e12
    cleavage_hits: float = 3.0
    cleavage_correlation_time_s: float = 1.0e-6
    emission_nu0_s: float = 1.0e11
    emission_signs: tuple[int, ...] = (1, -1)
    emission_schmid_factors: tuple[float, ...] = (1.0, 1.0)
    shielding_orientation_factors: tuple[float, ...] = (1.0, -1.0)
    forest_interaction_matrix: tuple[tuple[float, ...], ...] = (
        (1.0, 1.0),
        (1.0, 1.0),
    )
    gnd_stress_projection_matrix: tuple[tuple[float, ...], ...] = (
        (1.0, 0.0),
        (0.0, 1.0),
    )
    max_fractional_state_change: float = 0.05
    min_substep_s: float = 1.0e-12

    def validate(self) -> None:
        if self.n_bins < 1 or self.n_systems < 1:
            raise ValueError("n_bins and n_systems must be positive")
        for values, name in (
            (self.emission_signs, "emission_signs"),
            (self.emission_schmid_factors, "emission_schmid_factors"),
            (self.shielding_orientation_factors, "shielding_orientation_factors"),
        ):
            if len(values) != self.n_systems:
                raise ValueError(f"{name} length must equal n_systems")
        forest = np.asarray(self.forest_interaction_matrix, dtype=float)
        if forest.shape != (self.n_systems, self.n_systems):
            raise ValueError("forest_interaction_matrix has wrong shape")
        if np.any(forest < 0.0):
            raise ValueError("forest_interaction_matrix must be nonnegative")
        gnd = np.asarray(self.gnd_stress_projection_matrix, dtype=float)
        if gnd.shape != (self.n_systems, self.n_systems):
            raise ValueError("gnd_stress_projection_matrix has wrong shape")
        if self.active_strip_width_m <= 0.0:
            raise ValueError("active_strip_width_m must be positive")


@dataclass(frozen=True)
class ProtocolSegment:
    extension_start_m: float
    extension_end_m: float
    K_start_MPa_sqrt_m: float
    K_end_MPa_sqrt_m: float
    duration_s: float

    @property
    def da_m(self) -> float:
        return self.extension_end_m - self.extension_start_m


@dataclass
class TemperatureResult:
    candidate_id: str
    temperature_K: float
    extensions_um: list[float] = field(default_factory=list)
    K_applied_MPa_sqrt_m: list[float] = field(default_factory=list)
    delta_K_micro_MPa_sqrt_m: list[float] = field(default_factory=list)
    K_shield_MPa_sqrt_m: list[float] = field(default_factory=list)
    tau_gnd_tip_MPa: list[float] = field(default_factory=list)
    retained_line_count_per_unit_thickness: list[float] = field(default_factory=list)
    gnd_abs_line_count_per_unit_thickness: list[float] = field(default_factory=list)
    source_available_fraction: list[float] = field(default_factory=list)
    pi_store_max: list[float] = field(default_factory=list)
    pi_release_max: list[float] = field(default_factory=list)

    def as_dict(self) -> dict[str, Any]:
        return dict(vars(self))
