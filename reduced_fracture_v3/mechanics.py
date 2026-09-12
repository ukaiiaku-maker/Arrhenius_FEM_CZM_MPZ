"""Fail-closed aligned response maps for the V3 mechanics reduction."""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
import math
from typing import Iterable, Mapping

import numpy as np


class MechanicsMapDomainError(ValueError):
    """A response-map query is outside its qualified rectangular domain."""


class TopologyRegime(str, Enum):
    PRECONNECTION = "PRECONNECTION"
    CONNECTED_CAVITY = "CONNECTED_CAVITY"
    DOWNSTREAM_CHILD = "DOWNSTREAM_CHILD"


@dataclass(frozen=True)
class ResponseMapNode:
    regime: TopologyRegime
    ligament_over_radius: float
    radius_over_scale: float
    delta_G_per_load2: float | None
    delta_sigma_nn_per_load: float | None
    delta_sigma_tt_per_load: float | None
    delta_tau_nt_per_load: float | None
    delta_sigma_m_per_load: float | None
    delta_compliance: float
    delta_potential_energy_per_load2: float
    source_identity: str


@dataclass(frozen=True)
class MechanicsResponse:
    regime: TopologyRegime
    G_kinetic: float | None
    cavity_tensor: tuple[float, float, float, float] | None
    compliance: float
    potential_energy: float
    map_source_identities: tuple[str, ...]
    qualification: str = "QUALIFIED_BOUNDED_INTERPOLATION"


class TensorProductResponseMap:
    """Two-axis piecewise-linear map with complete-grid qualification."""

    _OUTPUTS = (
        "delta_G_per_load2",
        "delta_sigma_nn_per_load",
        "delta_sigma_tt_per_load",
        "delta_tau_nt_per_load",
        "delta_sigma_m_per_load",
        "delta_compliance",
        "delta_potential_energy_per_load2",
    )

    def __init__(self, nodes: Iterable[ResponseMapNode], regime: TopologyRegime):
        selected = tuple(node for node in nodes if node.regime == regime)
        if not selected:
            raise ValueError(f"no nodes supplied for {regime.value}")
        self.regime = regime
        self.ligament_axis = np.array(sorted({node.ligament_over_radius for node in selected}))
        self.radius_axis = np.array(sorted({node.radius_over_scale for node in selected}))
        expected = len(self.ligament_axis) * len(self.radius_axis)
        keys = {(node.ligament_over_radius, node.radius_over_scale) for node in selected}
        if len(selected) != expected or len(keys) != expected:
            raise ValueError("response map requires one complete Cartesian geometry grid")
        by_key = {(node.ligament_over_radius, node.radius_over_scale): node for node in selected}
        self.grids: dict[str, np.ndarray | None] = {}
        for output in self._OUTPUTS:
            values = [getattr(node, output) for node in selected]
            if all(value is None for value in values):
                self.grids[output] = None
                continue
            if any(value is None or not math.isfinite(float(value)) for value in values):
                raise ValueError(f"{output} must be complete and finite or uniformly absent")
            self.grids[output] = np.array(
                [
                    [float(getattr(by_key[(ligament, radius)], output)) for radius in self.radius_axis]
                    for ligament in self.ligament_axis
                ]
            )
        if regime == TopologyRegime.CONNECTED_CAVITY and self.grids["delta_G_per_load2"] is not None:
            raise ValueError("the connected-cavity map cannot expose a sharp-tip G")
        if regime in (TopologyRegime.PRECONNECTION, TopologyRegime.DOWNSTREAM_CHILD):
            if self.grids["delta_G_per_load2"] is None:
                raise ValueError("sharp-front map requires candidate kinetic G correction")
        stress_presence = [
            self.grids[name] is not None
            for name in (
                "delta_sigma_nn_per_load",
                "delta_sigma_tt_per_load",
                "delta_tau_nt_per_load",
                "delta_sigma_m_per_load",
            )
        ]
        if any(stress_presence) and not all(stress_presence):
            raise ValueError("cavity tensor components must be complete or uniformly absent")
        self.source_identities = tuple(sorted({node.source_identity for node in selected}))

    def _at(self, grid: np.ndarray | None, ligament: float, radius: float) -> float | None:
        if grid is None:
            return None
        if (
            ligament < self.ligament_axis[0]
            or ligament > self.ligament_axis[-1]
            or radius < self.radius_axis[0]
            or radius > self.radius_axis[-1]
        ):
            raise MechanicsMapDomainError(
                "mechanics query outside qualified domain: "
                f"L/R={ligament:.12g}, R/l*={radius:.12g}"
            )
        along_radius = np.array(
            [np.interp(radius, self.radius_axis, row) for row in grid]
        )
        return float(np.interp(ligament, self.ligament_axis, along_radius))

    def evaluate(
        self,
        *,
        ligament_over_radius: float,
        radius_over_scale: float,
        load: float,
        no_void_G_kinetic: float | None,
        no_void_compliance: float,
        no_void_potential_energy: float,
        isolated_cavity_tensor_per_load: tuple[float, float, float, float] | None,
    ) -> MechanicsResponse:
        load_value = float(load)
        if not math.isfinite(load_value):
            raise ValueError("load must be finite")
        if self.regime == TopologyRegime.CONNECTED_CAVITY and no_void_G_kinetic is not None:
            raise ValueError("connected-cavity mechanics cannot receive or return sharp-tip G")
        if radius_over_scale == 0.0 or math.isinf(ligament_over_radius):
            tensor = None
            if isolated_cavity_tensor_per_load is not None:
                tensor = tuple(float(value) * load_value for value in isolated_cavity_tensor_per_load)
            return MechanicsResponse(
                regime=self.regime,
                G_kinetic=no_void_G_kinetic,
                cavity_tensor=tensor,
                compliance=float(no_void_compliance),
                potential_energy=float(no_void_potential_energy),
                map_source_identities=self.source_identities,
                qualification="EXACT_FAR_OR_SMALL_VOID_REDUCTION_LIMIT",
            )
        values = {
            name: self._at(grid, ligament_over_radius, radius_over_scale)
            for name, grid in self.grids.items()
        }
        delta_g = values["delta_G_per_load2"]
        G = None
        if delta_g is not None:
            if no_void_G_kinetic is None:
                raise ValueError("sharp-front response requires the existing no-void kinetic G")
            G = float(no_void_G_kinetic) + delta_g * load_value**2
        stress_values = tuple(values[name] for name in (
            "delta_sigma_nn_per_load",
            "delta_sigma_tt_per_load",
            "delta_tau_nt_per_load",
            "delta_sigma_m_per_load",
        ))
        tensor = None
        if not all(value is None for value in stress_values):
            if isolated_cavity_tensor_per_load is None:
                raise ValueError("cavity-tensor correction requires the isolated-cavity baseline")
            tensor = tuple(
                (float(base) + float(correction)) * load_value
                for base, correction in zip(isolated_cavity_tensor_per_load, stress_values)
            )
        return MechanicsResponse(
            regime=self.regime,
            G_kinetic=G,
            cavity_tensor=tensor,
            compliance=float(no_void_compliance) + float(values["delta_compliance"]),
            potential_energy=(
                float(no_void_potential_energy)
                + float(values["delta_potential_energy_per_load2"]) * load_value**2
            ),
            map_source_identities=self.source_identities,
        )


@dataclass(frozen=True)
class AlignedVoidMechanicsBackend:
    maps: Mapping[TopologyRegime, TensorProductResponseMap]

    def evaluate(self, regime: TopologyRegime, **query: object) -> MechanicsResponse:
        if regime not in self.maps:
            raise MechanicsMapDomainError(f"no qualified map for {regime.value}")
        return self.maps[regime].evaluate(**query)

    @staticmethod
    def no_void_response(
        regime: TopologyRegime,
        *,
        G_kinetic: float | None,
        compliance: float,
        potential_energy: float,
    ) -> MechanicsResponse:
        if regime == TopologyRegime.CONNECTED_CAVITY and G_kinetic is not None:
            raise ValueError("connected cavity has no sharp-tip G")
        return MechanicsResponse(
            regime=regime,
            G_kinetic=G_kinetic,
            cavity_tensor=None,
            compliance=float(compliance),
            potential_energy=float(potential_energy),
            map_source_identities=(),
            qualification="EXACT_NO_VOID_LIMIT",
        )


def relative_error(observed: float, expected: float, floor: float = 1.0e-30) -> float:
    return abs(float(observed) - float(expected)) / max(abs(float(expected)), float(floor))


def held_out_geometry_family(
    nodes: Iterable[ResponseMapNode],
    *,
    regime: TopologyRegime,
    radius_over_scale: float,
) -> dict[str, float | int | str]:
    """Hold out a complete radius family and score bounded interpolation."""

    nodes = tuple(node for node in nodes if node.regime == regime)
    held = tuple(node for node in nodes if node.radius_over_scale == radius_over_scale)
    training = tuple(node for node in nodes if node.radius_over_scale != radius_over_scale)
    if not held:
        raise ValueError("held-out radius family is absent")
    model = TensorProductResponseMap(training, regime)
    errors: dict[str, list[float]] = {name: [] for name in TensorProductResponseMap._OUTPUTS}
    for node in held:
        for name, grid in model.grids.items():
            expected = getattr(node, name)
            if expected is not None:
                observed = model._at(grid, node.ligament_over_radius, node.radius_over_scale)
                errors[name].append(relative_error(float(observed), float(expected)))
    return {
        "regime": regime.value,
        "held_out_radius_over_scale": float(radius_over_scale),
        "training_nodes": len(training),
        "held_out_nodes": len(held),
        **{
            f"maximum_{name}_relative_error": max(values) if values else 0.0
            for name, values in errors.items()
        },
    }


__all__ = [
    "AlignedVoidMechanicsBackend",
    "MechanicsMapDomainError",
    "MechanicsResponse",
    "ResponseMapNode",
    "TensorProductResponseMap",
    "TopologyRegime",
    "held_out_geometry_family",
]
