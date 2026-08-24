from __future__ import annotations
from dataclasses import dataclass
from typing import Any
from enum import Enum

FINAL_POST_EVENT_RIGHT_CENSORED="FINAL_POST_EVENT_RIGHT_CENSORED"

class PFLateGeometryVetoPolicy(str,Enum):
 FAIL_CLOSED_TERMINATE="FAIL_CLOSED_TERMINATE"

@dataclass(frozen=True)
class PFTerminalVeto:
 terminal_class:str="NUMERICAL_BACKEND_FAIL_CLOSED_VETO"
 rollback_and_continue_supported:bool=False
 tentative_geometry_published:bool=False
 physical_arrest:bool=False
 successful_target_completion:bool=False

def pf_fail_closed_veto(last_authoritative_state, tentative_path_fingerprint, message):
 return {"policy":PFLateGeometryVetoPolicy.FAIL_CLOSED_TERMINATE.value,"last_authoritative_state":last_authoritative_state,"tentative_path_fingerprint":tentative_path_fingerprint,"terminal":PFTerminalVeto(),"exception_message":str(message)}

@dataclass(frozen=True)
class BackendEventResult:
 backend:str;fired:bool;event_count:int;threshold_before:float|None
 threshold_after:float|None;geometry_advance_m:float;snapshot_restored:bool

class PFLifecyclePolicy:
 backend="PF"
 def snapshot(self,engine):return engine.mpz.copy(),engine.B,engine.a_adv,engine.n_adv,engine.W_emit,engine.t,engine.K_prev,engine._lambda_c_prev
 def restore(self,engine,snapshot):engine.mpz,engine.B,engine.a_adv,engine.n_adv,engine.W_emit,engine.t,engine.K_prev,engine._lambda_c_prev=snapshot
 def step(self,engine,K,T,dt):
  b=float(engine.B);a=float(engine.a_adv);out=engine.step(K,T,dt)
  return BackendEventResult(self.backend,bool(out["fired"]),int(out["n_fire"]),b,float(engine.B),float(engine.a_adv-a),False)

class FEMCZMLifecyclePolicy:
 backend="FEMCZM"
 def snapshot(self,state):return state.copy()
 def restore(self,state,snapshot):return snapshot.copy()
 def no_event_update(self,state,dt,T,tip_stress,b):return state.evolve(dt,T,tip_stress,b)
 def advance_after_cleavage(self,state,distance_m):return state.advance(distance_m)
 def step_engine(self,engine,K,T,dt):
  b=float(engine.B);a=float(engine.a_adv);before=int(engine.n_adv);out=engine.step_drives(K,K,T,dt)
  return BackendEventResult(self.backend,bool(out["fired"]),int(engine.n_adv-before),b,float(engine.B),float(engine.a_adv-a),False)

# Backward-compatible names used by the lifecycle qualification artifacts.
PFEventLifecyclePolicy=PFLifecyclePolicy
FEMCZMEventLifecyclePolicy=FEMCZMLifecyclePolicy
