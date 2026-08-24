from __future__ import annotations

from dataclasses import dataclass, replace


@dataclass(frozen=True)
class ReducedState:
    time_s: float = 0.0
    external_control: float = 0.0
    crack_extension_m: float = 0.0
    cleavage_action: float = 0.0
    cleavage_threshold: float = 1.0
    emission_action: float = 0.0
    emission_threshold: float = 1.0
    mobile_density_m2: float = 0.0
    retained_density_m2: float = 0.0
    local_slip: float = 0.0
    tip_radius_m: float = 1e-9
    front_width_m: float = 1e-9
    source_area_m2: float = 1e-18
    source_multiplicity: float = 1.0
    backstress_Pa: float = 0.0
    signed_shielding_Pa_sqrt_m: float = 0.0
    event_count: int = 0
    physical_avalanche_index: int = 0

    def with_hazard(self, cleavage_increment: float, emission_increment: float, dt_s: float) -> "ReducedState":
        return replace(self, time_s=self.time_s + dt_s, cleavage_action=self.cleavage_action + cleavage_increment, emission_action=self.emission_action + emission_increment)
