from __future__ import annotations
from dataclasses import dataclass
import copy,hashlib
import numpy as np

def _digest(fields):
 h=hashlib.sha256()
 for k,v in sorted(fields.items()):
  h.update(k.encode());a=np.asarray(v);h.update(str(a.dtype).encode());h.update(str(a.shape).encode());h.update(a.tobytes())
 return h.hexdigest()

@dataclass(frozen=True)
class ProcessZoneSnapshot:
 backend:str;fields:dict;fingerprint:str

class PFProcessZoneStateAdapter:
 owner="PF UnifiedMPZState";names=("mobile","retained","accumulated_slip","wake_mobile","wake_retained","wake_slip","available_sites","site_capacity","time_s","emitted_total","escaped_total","recovered_total","advance_total_m")
 def snapshot(self,state):
  f={n:copy.deepcopy(getattr(state,n)) for n in self.names};return ProcessZoneSnapshot("PF",f,_digest(f))
 def restore(self,state,snapshot):
  if snapshot.backend!="PF":raise ValueError("backend mismatch")
  for n,v in snapshot.fields.items():setattr(state,n,copy.deepcopy(v))

class FEMCZMProcessZoneStateAdapter:
 owner="FEM/CZM MovingProcessZoneState v9.11";names=("mobile","retained","accumulated_slip","time_s","emitted_total","escaped_total","recovered_total","advance_total_m")
 def snapshot(self,state):
  f={n:copy.deepcopy(getattr(state,n)) for n in self.names};return ProcessZoneSnapshot("FEMCZM",f,_digest(f))
 def restore(self,state,snapshot):
  if snapshot.backend!="FEMCZM":raise ValueError("backend mismatch")
  for n,v in snapshot.fields.items():setattr(state,n,copy.deepcopy(v))
