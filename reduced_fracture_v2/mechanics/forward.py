from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from .base import MechanicsProvider, MechanicsState, ObservableKind, Qualification


@dataclass
class TabulatedForwardProvider(MechanicsProvider):
    extensions_m: tuple[float,...]
    coefficients: tuple[float,...]
    observable_kind: ObservableKind
    qualification: Qualification
    map_fingerprint: str
    candidate_independent: bool=True

    def evaluate(self,extension_m:float,opening_m:float)->MechanicsState:
        if not self.candidate_independent: raise ValueError("mechanics map must be candidate independent")
        if len(self.extensions_m)<2 or extension_m<self.extensions_m[0] or extension_m>self.extensions_m[-1]: raise ValueError("forward mechanics map does not cover requested extension")
        c=float(np.interp(extension_m,self.extensions_m,self.coefficients))
        if self.observable_kind==ObservableKind.FEM_STRUCTURAL_G:
            g=c*opening_m**2; return MechanicsState("OPENING_M",opening_m,extension_m,self.observable_kind,None,g,None,None,None,self.qualification)
        return MechanicsState("OPENING_M",opening_m,extension_m,self.observable_kind,c*opening_m,None,None,None,None,self.qualification)
    def evaluate_pre_event(self,index): raise RuntimeError("forward provider requires extension and opening")
    def evaluate_post_event(self,index): raise RuntimeError("forward provider requires extension and opening")
    def evaluate_reload(self,index): raise RuntimeError("forward provider requires extension and opening")
