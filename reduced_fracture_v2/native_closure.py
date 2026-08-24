from __future__ import annotations
from dataclasses import dataclass,replace
from typing import Any
import hashlib,json
import numpy as np

from .kinetics import KernelResult,SharedBarrierHazardCore
from .lifecycle import PFLifecyclePolicy,FEMCZMLifecyclePolicy


def _json_hash(value:Any)->str:
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(",",":"),default=str).encode()).hexdigest()


@dataclass(frozen=True)
class PFNativeForwardState:
    time_s:float;opening_m:float;extension_m:float
    native_J_J_m2:float;native_KJ_Pa_sqrt_m:float
    mobile_count:float;retained_count:float;tip_radius_m:float;front_width_m:float
    source_area_m2:float;source_multiplicity:float;backstress_Pa:float
    signed_shielding_Pa_sqrt_m:float;effective_K_Pa_sqrt_m:float;tip_stress_Pa:float
    cleavage_action:float;cleavage_threshold:float;threshold_identity:str
    process_zone_fingerprint:str;event_count:int;physical_avalanche_index:int
    last_committed_geometry_fingerprint:str;terminal_veto_state:str|None=None


@dataclass(frozen=True)
class FEMCZMNativeForwardState:
    time_s:float;opening_m:float;extension_m:float
    native_J_J_m2:float;native_KJ_Pa_sqrt_m:float
    qualified_G_J_m2:float;qualified_KG_Pa_sqrt_m:float;G_uncertainty_J_m2:float;J_status:str
    mobile_count:float;retained_count:float;tip_radius_m:float;front_width_m:float
    source_area_m2:float;source_multiplicity:float;backstress_Pa:float
    signed_shielding_Pa_sqrt_m:float;effective_K_Pa_sqrt_m:float;tip_stress_Pa:float
    cleavage_action:float;cleavage_threshold:float;threshold_RNG_identity:str
    process_zone_fingerprint:str;transaction_snapshot_identity:str
    event_count:int;physical_avalanche_index:int;last_committed_geometry_fingerprint:str


@dataclass(frozen=True)
class NativeIntervalResult:
    pre_state:PFNativeForwardState|FEMCZMNativeForwardState
    post_state:PFNativeForwardState|FEMCZMNativeForwardState
    shared_core:KernelResult
    fired:bool;event_length_m:float;geometry_commits:int;post_event_renewals:int


class _NativeClosure:
    def __init__(self,mechanics,engine,adapter,parameters):
        self.mechanics=mechanics;self.engine=engine;self.adapter=adapter
        self.core=SharedBarrierHazardCore(parameters);self.extension_m=0.;self.event_count=0
        self.geometry_fingerprint=_json_hash({"extension_m":0.,"event_count":0})
    def _pz(self):
        return getattr(self.engine,"mpz",getattr(self.engine,"mpz_state",None))
    def _counts(self):
        p=self._pz();return float(getattr(p,"mobile_count",np.sum(getattr(p,"mobile",0.)))),float(getattr(p,"retained_count",np.sum(getattr(p,"retained",0.))))
    def _threshold(self):
        target=float(getattr(self.engine,"B_target",1.));stream=getattr(self.engine,"_threshold_stream",None)
        identity=_json_hash(stream.snapshot() if stream is not None else {"deterministic_target":target})
        return float(getattr(self.engine,"B",0.)),target,identity
    def _scalars(self,K):
        p=self._pz();mobile,retained=self._counts();r=float(self.engine.r_eff())
        kshield=float(self.engine.K_shield()) if hasattr(self.engine,"K_shield") else max(K-float(self.engine.sigma_tip(K))*np.sqrt(2*np.pi*r),0.)
        sig=float(self.engine.sigma_tip(K));back=float(self.engine.sigma_back()) if hasattr(self.engine,"sigma_back") else kshield/np.sqrt(2*np.pi*r)
        diag={}
        if hasattr(p,"diagnostics"):
            try: diag=p.diagnostics(self.engine.G,self.engine.nu,self.engine.b,self.engine.f.r0,self.engine.f.c_blunt)
            except TypeError: diag=p.diagnostics(self.engine.G,self.engine.nu,self.engine.b,self.engine.f.r0)
        width=float(diag.get("persistent_front_width_m",diag.get("front_width_m",self.engine.f.L_pz)))
        area=float(diag.get("persistent_source_area_m2",max(r*width,0.)))
        mult=float(diag.get("persistent_site_multiplicity_total",getattr(p,"source_multiplicity",1.)))
        return mobile,retained,r,width,area,mult,back,kshield,max(K-kshield,0.),sig
    def _core(self,K,T,dt):
        *_,back,_,_,sig=self._scalars(K)
        from .state import ReducedState
        state=ReducedState(source_multiplicity=self._scalars(K)[5],backstress_Pa=back)
        return self.core.evaluate(state,T,sig,max(sig-back,0.),dt)


class PFNativeStateClosure(_NativeClosure):
    def __init__(self,mechanics,engine,adapter,parameters):super().__init__(mechanics,engine,adapter,parameters);self.lifecycle=PFLifecyclePolicy()
    def state(self,opening):
        m=self.mechanics.evaluate(self.extension_m,opening);vals=self._scalars(m.pf_native_KJ_Pa_sqrt_m);B,target,tid=self._threshold();snap=self.adapter.snapshot(self._pz())
        return PFNativeForwardState(float(self.engine.t),opening,self.extension_m,m.pf_native_J_J_m2,m.pf_native_KJ_Pa_sqrt_m,*vals[:6],vals[6],vals[7],vals[8],vals[9],B,target,tid,snap.fingerprint,self.event_count,0,self.geometry_fingerprint)
    def advance_interval(self,opening,T,dt):
        pre=self.state(opening);audit=self._core(pre.native_KJ_Pa_sqrt_m,T,dt);res=self.lifecycle.step(self.engine,pre.native_KJ_Pa_sqrt_m,T,dt)
        commits=renewals=0
        if res.fired:
            self.extension_m+=res.geometry_advance_m;self.event_count+=res.event_count;commits=renewals=1
            self.geometry_fingerprint=_json_hash({"extension_m":self.extension_m,"event_count":self.event_count})
        return NativeIntervalResult(pre,self.state(opening),audit,res.fired,res.geometry_advance_m,commits,renewals)
    def evaluate_reload(self,opening):return self.state(opening)


class FEMCZMNativeStateClosure(_NativeClosure):
    def __init__(self,mechanics,engine,adapter,parameters):super().__init__(mechanics,engine,adapter,parameters);self.lifecycle=FEMCZMLifecyclePolicy()
    def state(self,opening):
        m=self.mechanics.evaluate(self.extension_m,opening);vals=self._scalars(m.fem_native_KJ_Pa_sqrt_m);B,target,tid=self._threshold();snap=self.adapter.snapshot(self._pz())
        return FEMCZMNativeForwardState(float(self.engine.t),opening,self.extension_m,m.fem_native_J_J_m2,m.fem_native_KJ_Pa_sqrt_m,m.structural_G_J_m2,m.qualified_KG_Pa_sqrt_m,m.G_uncertainty_J_m2,m.J_status,*vals[:6],vals[6],vals[7],vals[8],vals[9],B,target,tid,snap.fingerprint,_json_hash({"B":B,"threshold":tid}),self.event_count,0,self.geometry_fingerprint)
    def advance_interval(self,opening,T,dt):
        pre=self.state(opening);audit=self._core(pre.native_KJ_Pa_sqrt_m,T,dt);res=self.lifecycle.step_engine(self.engine,pre.native_KJ_Pa_sqrt_m,T,dt)
        commits=renewals=0
        if res.fired:
            self.extension_m+=res.geometry_advance_m;self.event_count+=res.event_count;commits=renewals=1
            self.geometry_fingerprint=_json_hash({"extension_m":self.extension_m,"event_count":self.event_count})
        return NativeIntervalResult(pre,self.state(opening),audit,res.fired,res.geometry_advance_m,commits,renewals)
    def evaluate_reload(self,opening):return self.state(opening)
