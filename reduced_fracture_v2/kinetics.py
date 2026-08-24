from __future__ import annotations

from dataclasses import dataclass
import math

from .parameters import CanonicalParameters
from .state import ReducedState

KB_EV_K = 8.617333262145e-5


@dataclass(frozen=True)
class KernelResult:
    cleavage_barrier_eV: float
    emission_barrier_eV: float
    cleavage_rate_s: float
    emission_rate_s: float
    peierls_rate_s: float
    taylor_rate_s: float
    cleavage_hazard_increment: float
    emission_hazard_increment: float


class SharedBarrierHazardCore:
    """Backend-free canonical EXP-floor barrier and hazard evaluation.

    Mechanics supplies local stresses. This class neither knows nor labels J/G/K.
    Process-zone transport remains an explicit qualification boundary until PF
    and FEM/CZM frozen-state source functions are proven equivalent.
    """

    def __init__(self, parameters: CanonicalParameters): self.p=parameters.values

    def barrier_eV(self, prefix: str, stress_Pa: float, temperature_K: float) -> float:
        p=self.p; tref=p["Tref_K"]
        g0=max(p[f"{prefix}_G00_eV"]+p[f"{prefix}_gT_eV_per_K"]*(temperature_K-tref),1e-8)
        sigc=max((p[f"{prefix}_sigc0_GPa"]+p[f"{prefix}_sT_GPa_per_K"]*(temperature_K-tref))*1e9,1e6)
        floor=max(0.0,min(g0,p[f"{prefix}_floor_frac"]*g0))
        x=(abs(stress_Pa)/sigc)**max(p[f"{prefix}_exp_n"],1e-8)
        return floor+(g0-floor)*math.exp(-max(p[f"{prefix}_exp_a"],0.0)*x)

    def _rate(self, barrier: float, temperature_K: float, prefactor: float) -> float:
        return prefactor*math.exp(max(-745.0,min(700.0,-barrier/(KB_EV_K*temperature_K))))

    def evaluate(self, state: ReducedState, temperature_K: float, cleavage_stress_Pa: float, emission_stress_Pa: float, dt_s: float) -> KernelResult:
        c=self.barrier_eV("cleave",cleavage_stress_Pa,temperature_K); e=self.barrier_eV("emit",emission_stress_Pa-state.backstress_Pa,temperature_K)
        cr=self._rate(c,temperature_K,1e13)*state.source_multiplicity; er=self._rate(e,temperature_K,1e13)*state.source_multiplicity
        pr=self._rate(self.p["peierls_H0_eV"],temperature_K,self.p["peierls_nu0_s"]); tr=self._rate(self.p["taylor_H0_eV"],temperature_K,self.p["taylor_nu0_s"])
        return KernelResult(c,e,cr,er,pr,tr,cr*dt_s,er*dt_s)


ReducedFractureKernel = SharedBarrierHazardCore
