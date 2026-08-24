from __future__ import annotations
from dataclasses import dataclass,replace
from .analysis import SharedEventAvalancheAnalysis
from .kinetics import SharedBarrierHazardCore
from .state import ReducedState


@dataclass(frozen=True)
class ControlledParityStep:
    opening_m:float
    dt_s:float
    cleavage_stress_Pa:float
    emission_stress_Pa:float
    event_extension_m:float=0.0
    reload_before_event:bool=False
    right_censored_at_target:bool=False


class ControlledParityLifecycle:
    """Analysis-only common lifecycle; never used as a production policy."""
    def apply(self,state,step,result,event_index):
        s=state.with_hazard(result.cleavage_hazard_increment,result.emission_hazard_increment,step.dt_s)
        s=replace(s,external_control=step.opening_m)
        row=None
        if step.event_extension_m>0:
            row={"event_transaction_index":event_index,"pre_event_time_s":s.time_s,
                 "pre_event_opening_m":step.opening_m,"pre_event_extension_m":s.crack_extension_m,
                 "event_extension_m":step.event_extension_m,"reload_before_event":step.reload_before_event,
                 "right_censored_at_target":step.right_censored_at_target}
            s=replace(s,crack_extension_m=s.crack_extension_m+step.event_extension_m,
                      cleavage_action=0.0,cleavage_threshold=1.0,event_count=s.event_count+1)
        return s,row


class _ParityDriverShell:
    backend="UNSET"
    def __init__(self,parameters):
        self.kernel=SharedBarrierHazardCore(parameters);self.lifecycle=ControlledParityLifecycle()
    def run(self,initial_state,temperature_K,steps):
        state=initial_state;rows=[]
        for step in steps:
            z=self.kernel.evaluate(state,temperature_K,step.cleavage_stress_Pa,step.emission_stress_Pa,step.dt_s)
            state,row=self.lifecycle.apply(state,step,z,len(rows))
            if row is not None:
                row.update({"backend_shell":self.backend,"post_event_extension_m":state.crack_extension_m,
                            "post_event_mobile_density_m2":state.mobile_density_m2,
                            "post_event_retained_density_m2":state.retained_density_m2})
                rows.append(row)
        return {"terminal_state":state,"event_transactions":rows,
                "physical_avalanches":SharedEventAvalancheAnalysis.group(rows)}


class PFParityDriverShell(_ParityDriverShell): backend="PF"
class FEMCZMParityDriverShell(_ParityDriverShell): backend="FEMCZM"
