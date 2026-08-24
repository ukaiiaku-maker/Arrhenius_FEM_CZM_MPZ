from __future__ import annotations

from dataclasses import dataclass
import csv
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


class PFMechanicsProvider(TabulatedForwardProvider):
    def __post_init__(self):
        if self.observable_kind not in (ObservableKind.PF_NATIVE_DRIVE,ObservableKind.PF_NATIVE_KJ): raise ValueError("PF provider requires PF-native observable")
        if self.qualification==Qualification.QUALIFIED_ENERGY_COMPLIANCE_VCCT: raise ValueError("PF discrete drive is not qualified continuum G")

    @classmethod
    def from_csv(cls,path,map_fingerprint):
        with open(path,newline="") as f: rows=list(csv.DictReader(f))
        return _TypedPFProvider(
            tuple(float(r["actual_extension_um"])*1e-6 for r in rows),
            tuple(float(r["KJ_native_over_U"])*1e6 for r in rows),
            ObservableKind.PF_NATIVE_KJ,Qualification.MODEL_NATIVE_NOT_STRUCTURAL_G,
            map_fingerprint,True,
            tuple(float(r["F_over_U_N_per_m2"]) for r in rows),
            tuple(float(r["J_native_over_U2"]) for r in rows),
        )


class FEMCZMMechanicsProvider(TabulatedForwardProvider):
    def __post_init__(self):
        if self.observable_kind!=ObservableKind.FEM_STRUCTURAL_G: raise ValueError("FEM/CZM primary provider requires structural G")
        if self.qualification!=Qualification.QUALIFIED_ENERGY_COMPLIANCE_VCCT: raise ValueError("FEM/CZM structural G map must be qualified")

    @classmethod
    def from_csv(cls,native_path,g_path,map_fingerprint):
        with open(native_path,newline="") as f: n=list(csv.DictReader(f))
        with open(g_path,newline="") as f: g=list(csv.DictReader(f))
        if [r["actual_extension_um"] for r in n] != [r["actual_extension_um"] for r in g]:
            raise ValueError("native-J and qualified-G nodes differ")
        if any(r["node_qualification"]!="QUALIFIED" for r in g):
            raise ValueError("unqualified structural-G node")
        return _TypedFEMProvider(
            tuple(float(r["actual_extension_um"])*1e-6 for r in n),
            tuple(float(r["G_over_U2"]) for r in g),
            ObservableKind.FEM_STRUCTURAL_G,Qualification.QUALIFIED_ENERGY_COMPLIANCE_VCCT,
            map_fingerprint,True,
            tuple(float(r["F_over_U_N_per_m2"]) for r in n),
            tuple(float(r["J_native_over_U2"]) for r in n),
            tuple(float(r["KJ_native_over_U"])*1e6 for r in n),
            tuple(float(r["KG_over_U"])*1e6 for r in g),
            tuple(float(r["G_uncertainty_J_per_m2"])/(float(r["reference_opening_m"])**2) for r in g),
            tuple(r["J_status"] for r in n),
        )


@dataclass
class _TypedPFProvider(PFMechanicsProvider):
    reaction_over_U: tuple[float,...]=()
    J_over_U2: tuple[float,...]=()

    def evaluate(self,extension_m,opening_m):
        s=super().evaluate(extension_m,opening_m)
        j=float(np.interp(extension_m,self.extensions_m,self.J_over_U2))*opening_m**2
        return MechanicsState("OPENING_M",opening_m,extension_m,ObservableKind.PF_NATIVE_KJ,
            s.native_drive,None,None,None,float(np.interp(extension_m,self.extensions_m,self.reaction_over_U))*opening_m,
            Qualification.MODEL_NATIVE_NOT_STRUCTURAL_G,pf_native_J_J_m2=j,
            pf_native_KJ_Pa_sqrt_m=s.native_drive,map_qualification="PF_PRODUCTION_DISCRETE_MAP_QUALIFIED",
            veto_policy="FAIL_CLOSED_TERMINATE")


@dataclass
class _TypedFEMProvider(FEMCZMMechanicsProvider):
    reaction_over_U: tuple[float,...]=()
    J_over_U2: tuple[float,...]=()
    KJ_over_U: tuple[float,...]=()
    KG_over_U: tuple[float,...]=()
    G_uncertainty_over_U2: tuple[float,...]=()
    J_statuses: tuple[str,...]=()

    def evaluate(self,extension_m,opening_m):
        s=super().evaluate(extension_m,opening_m)
        interp=lambda xs:float(np.interp(extension_m,self.extensions_m,xs))
        # Status is conservative across an interpolation interval.
        i=max(0,min(len(self.extensions_m)-2,int(np.searchsorted(self.extensions_m,extension_m)-1)))
        status=self.J_statuses[i] if self.J_statuses[i]==self.J_statuses[i+1] else "UNQUALIFIED"
        return MechanicsState("OPENING_M",opening_m,extension_m,ObservableKind.FEM_STRUCTURAL_G,
            interp(self.KJ_over_U)*opening_m,s.structural_G_J_m2,interp(self.KG_over_U)*opening_m,None,
            interp(self.reaction_over_U)*opening_m,Qualification.QUALIFIED_ENERGY_COMPLIANCE_VCCT,
            fem_native_J_J_m2=interp(self.J_over_U2)*opening_m**2,
            fem_native_KJ_Pa_sqrt_m=interp(self.KJ_over_U)*opening_m,
            qualified_KG_Pa_sqrt_m=interp(self.KG_over_U)*opening_m,
            G_uncertainty_J_m2=interp(self.G_uncertainty_over_U2)*opening_m**2,
            J_status=status,map_qualification="FEM_NATIVE_AND_STRUCTURAL_G_MAPS_QUALIFIED")
