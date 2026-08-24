from __future__ import annotations
from dataclasses import dataclass, asdict

@dataclass(frozen=True)
class OneDEventTransaction:
    implementation: str; material_class: str; candidate_id: str; option_key: str
    temperature_K: float; event_transaction_index: int; physical_avalanche_index: int | None
    pre_event_time_s: float; pre_event_load_kind: str; pre_event_load: float
    pre_event_displacement_m: float; pre_event_native_drive_kind: str
    pre_event_native_drive_MPa_sqrt_m: float; post_event_extension_um: float
    event_extension_um: float; post_event_native_drive_at_fixed_load_MPa_sqrt_m: float | None
    event_induced_native_drive_change_MPa_sqrt_m: float | None
    reload_to_next_event_MPa_sqrt_m: float | None; reload_to_next_event_displacement_m: float | None
    reload_to_next_event_time_s: float | None; physical_avalanche_role: str
    right_censored_at_target: bool; post_event_state: str

@dataclass(frozen=True)
class OneDPhysicalAvalanche:
    implementation: str; material_class: str; candidate_id: str; temperature_K: float
    physical_avalanche_index: int; first_event_transaction_index: int
    last_event_transaction_index: int; event_transaction_count: int
    start_extension_um: float; end_extension_um: float; extension_um: float
    onset_native_drive_MPa_sqrt_m: float; classification: str
    grouping_method: str; grouping_log_gap_decades: float | None
    right_censored_at_target: bool

@dataclass(frozen=True)
class OneDOnsetCandidate:
    implementation: str; material_class: str; candidate_id: str; temperature_K: float
    event_transaction_index: int; physical_avalanche_index: int; onset_role: str
    pre_event_extension_um: float; pre_event_displacement_m: float
    onset_native_drive_kind: str; onset_native_drive_MPa_sqrt_m: float
    is_structural_G: bool; right_censored_avalanche_after_onset: bool

@dataclass(frozen=True)
class OneDInAvalancheDriveState:
    implementation: str; material_class: str; candidate_id: str; temperature_K: float
    event_transaction_index: int; physical_avalanche_index: int
    physical_avalanche_role: str; post_event_extension_um: float
    native_drive_kind: str; post_event_native_drive_MPa_sqrt_m: float
    is_resistance_point: bool; right_censored_at_target: bool

def row(value): return asdict(value)
