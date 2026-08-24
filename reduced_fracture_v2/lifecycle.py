from __future__ import annotations
from dataclasses import dataclass, replace
from .state import ReducedState


@dataclass(frozen=True)
class EventTransaction:
    event_transaction_index:int
    physical_avalanche_index:int|None
    pre_event_state:ReducedState
    post_event_geometry_extension_m:float
    post_event_process_zone_state:ReducedState|None
    reload_before_next_event:float|None
    role:str
    right_censored_at_target:bool


class CommonEventLifecyclePolicy:
    """Transactional first-passage policy; backend differences stay in mechanics."""
    def fire(self,state:ReducedState,event_length_m:float,target_m:float,avalanche_index:int|None)->EventTransaction:
        if event_length_m<=0: raise ValueError("event length must be positive")
        end=state.crack_extension_m+event_length_m
        post=replace(state,crack_extension_m=end,event_count=state.event_count+1,cleavage_action=0.0)
        return EventTransaction(state.event_count,avalanche_index,state,end,post,None,"FINAL_POST_EVENT_RIGHT_CENSORED" if end>=target_m else "INTERIOR_POST_EVENT",end>=target_m)
