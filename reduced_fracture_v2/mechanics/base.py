from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from enum import Enum


class ObservableKind(str, Enum):
    PF_NATIVE_DRIVE="PF_NATIVE_DRIVE"
    PF_NATIVE_KJ="PF_NATIVE_KJ"
    FEM_STRUCTURAL_G="FEM_STRUCTURAL_G"


class Qualification(str, Enum):
    ARCHIVED_REPLAY="ARCHIVED_REPLAY"
    QUALIFIED_ENERGY_COMPLIANCE_VCCT="QUALIFIED_ENERGY_COMPLIANCE_VCCT"
    MODEL_NATIVE_NOT_STRUCTURAL_G="MODEL_NATIVE_NOT_STRUCTURAL_G"
    UNQUALIFIED="UNQUALIFIED"


@dataclass(frozen=True)
class MechanicsState:
    external_control_kind: str
    external_control: float
    crack_extension_m: float
    observable_kind: ObservableKind
    native_drive: float | None
    structural_G_J_m2: float | None
    K_G_Pa_sqrt_m: float | None
    local_tip_stress_Pa: float | None
    reaction_equivalent_N: float | None
    qualification: Qualification
    source_event_index: int | None = None

    def __post_init__(self):
        if self.observable_kind != ObservableKind.FEM_STRUCTURAL_G and self.structural_G_J_m2 is not None:
            raise ValueError("model-native PF drive cannot be labelled structural G")
        if self.observable_kind == ObservableKind.FEM_STRUCTURAL_G and self.structural_G_J_m2 is not None and self.qualification != Qualification.QUALIFIED_ENERGY_COMPLIANCE_VCCT and self.qualification != Qualification.ARCHIVED_REPLAY:
            raise ValueError("FEM G requires qualified mechanics provenance")


class MechanicsProvider(ABC):
    @abstractmethod
    def evaluate_pre_event(self,index:int)->MechanicsState: ...
    @abstractmethod
    def evaluate_post_event(self,index:int)->MechanicsState: ...
    @abstractmethod
    def evaluate_reload(self,index:int)->MechanicsState: ...
