#!/usr/bin/env python3
"""Build the fail-closed V2 qualification bundle without running 2-D physics."""
from __future__ import annotations
import argparse,csv,hashlib,json,os,tempfile
from pathlib import Path
import pandas as pd
from reduced_fracture_v2 import CanonicalParameters,ReducedFractureKernel,ReducedState

IDS={"v913_zeroD_sobol_0242980":"Peak","v913_zeroD_sobol_0202500":"DBTT","v913_zeroD_sobol_0129902":"weak-T","v913_zeroD_sobol_0077080":"ceramic-like"}
TEMPS=(300,600,900,1100,1200)
def atomic(p,s):
 p.parent.mkdir(parents=True,exist_ok=True);fd,n=tempfile.mkstemp(dir=p.parent,prefix="."+p.name)
 with os.fdopen(fd,"w") as f:f.write(s)
 os.replace(n,p)
def writecsv(p,rows,columns=None): atomic(p,pd.DataFrame(rows,columns=columns).to_csv(index=False,lineterminator="\n"))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
 q=argparse.ArgumentParser();q.add_argument("--registry",type=Path,required=True);q.add_argument("--pf-peak",type=Path,required=True);q.add_argument("--pf-dbtt",type=Path,required=True);q.add_argument("--fem-v2",type=Path,required=True);q.add_argument("--v1",type=Path,required=True);q.add_argument("--out",type=Path,required=True);a=q.parse_args();o=a.out;o.mkdir(parents=True,exist_ok=True)
 reg=pd.read_csv(a.registry,dtype=str);reg=reg[reg.candidate_id.isin(IDS)].copy();assert set(reg.candidate_id)==set(IDS)
 frozen=[]
 for row in reg.to_dict("records"):
  p=CanonicalParameters.from_mapping(row);k=ReducedFractureKernel(p)
  for T in TEMPS:
   for stress in (1e9,3e9,6e9):
    state=ReducedState(tip_radius_m=1e-9,front_width_m=5e-9,source_area_m2=1e-16,source_multiplicity=8,retained_density_m2=1e14,backstress_Pa=2e8)
    z=k.evaluate(state,T,stress,.8*stress,1e-6)
    for provider in ("PF_CONSISTENT","FEMCZM_CONSISTENT"):
     frozen.append({"candidate_id":p.candidate_id,"material_class":IDS[p.candidate_id],"temperature_K":T,"local_drive_stress_Pa":stress,"provider_configuration":provider,**z.__dict__,"common_kernel_fingerprint":p.fingerprint,"comparison_status":"COMMON_ALGEBRA_EQUAL_BY_CONSTRUCTION","source_function_direct_match_status":"BLOCKED_SOURCE_STATE_ADAPTER_NOT_QUALIFIED"})
 writecsv(o/"oneD_v2_frozen_state_consistency.csv",frozen)
 pf=[]
 for label,path in (("Peak",a.pf_peak),("DBTT",a.pf_dbtt)):
  d=pd.read_csv(path);pf.append({"backend":"PF","material_class":label,"temperature_K":1000,"archived_row_count":len(d),"event_count":"UNAVAILABLE_FROM_STEPS_WITHOUT_VALIDATED_EVENT_FILTER","event_times_replay":"UNAVAILABLE","event_lengths_replay":"UNAVAILABLE","process_zone_replay":"UNAVAILABLE_IN_COMPLETE_STATE_FORM","threshold_replay":"UNAVAILABLE_IN_COMPLETE_STATE_FORM","replay_status":"BLOCKED_MISSING_COMPLETE_PRE_POST_STATE"})
 tx=pd.read_csv(a.fem_v2/"peak_dbtt_event_transactions_v2.csv");av=pd.read_csv(a.fem_v2/"fem_physical_avalanches_v2.csv")
 for label in ("Peak","DBTT"):
  t=tx[tx.material.str.lower()==label.lower()];v=av[av.material.str.lower()==label.lower()]
  pf.append({"backend":"FEMCZM","material_class":label,"temperature_K":1000,"archived_row_count":len(t),"event_count":len(t),"event_times_replay":"AVAILABLE","event_lengths_replay":"AVAILABLE","process_zone_replay":"PARTIAL_NOT_COMPLETE","threshold_replay":"PARTIAL_NOT_COMPLETE","replay_status":"MECHANICS_LEDGER_REPLAYABLE_KINETICS_LIFECYCLE_BLOCKED"})
 writecsv(o/"oneD_v2_twoD_replay_summary.csv",pf);writecsv(o/"oneD_v2_event_lifecycle_comparison.csv",pf)
 empty=["oneD_v2_event_transactions.csv","oneD_v2_physical_avalanches.csv","oneD_v2_onset_candidates.csv","oneD_v2_in_avalanche_drive.csv","oneD_v2_four_class_temperature_summary.csv"]
 for n in empty:writecsv(o/n,[],["availability_status","blocker"])
 inventory={"schema":"oneD_v2_model_inventory_v1","V1":{"label":"ONE_D_V1_HISTORICAL","commit":"559425321b9a8739f32788322d8a1c2af8abad73","audit_head":"a5dece9","fingerprint":sha(a.v1),"preserved":True},"V2":{"label":"ONE_D_V2_COMMON_KERNEL","layers":["ReducedFractureKernel","MechanicsProvider","CommonEventLifecyclePolicy","StochasticSequence","V2 analysis"],"PF_replay":"implemented typed adapter; lifecycle qualification blocked by missing complete archived state","FEMCZM_replay":"implemented typed adapter; mechanics ledger available, complete kinetics state missing","PF_forward":"typed tabulated provider exists; no qualified actual-operator 1000um map","FEMCZM_forward":"typed G=M_G(a)U^2 provider exists; no committed qualified 1000um coefficient map"}}
 atomic(o/"oneD_v2_model_inventory.json",json.dumps(inventory,indent=2)+"\n")
 graph={"schema":"oneD_v2_shared_physics_graph_v1","shared":[{"quantity":"EXP-floor cleavage/emission barrier algebra","PF_source":"materials.py/config.py at 0a340f6","FEM_source":"materials.py/config.py at a5f0d63","units":"eV, Pa, K","sign":"absolute stress","status":"ALGEBRA_SHARED"},{"quantity":"canonical parameter rows","status":"SHARED_FULL_PRECISION_REGISTRY"}],"not_yet_proven_shared":["process-zone transport","tip-radius/front-width evolution","source multiplicity","backstress","signed shielding","post-event renewal"],"backend_specific":{"PF":"sharp-wake native J/KJ; continuum G unqualified","FEMCZM":"zero-thickness qualified G via energy/compliance/VCCT"}}
 atomic(o/"oneD_v2_shared_physics_graph.json",json.dumps(graph,indent=2)+"\n")
 decision={"schema":"oneD_v2_parameter_decision_v1","shared_kernel_gate":"PARTIAL_PASS_BARRIER_HAZARD_ALGEBRA","frozen_source_gate":"BLOCKED_STATE_ADAPTER","event_lifecycle_gate":"BLOCKED_MISSING_ARCHIVED_STATE","full_replay_gate":"BLOCKED_MISSING_ARCHIVED_STATE","PF_forward_map_gate":"BLOCKED_NO_QUALIFIED_GENERATIVE_MAP","FEMCZM_forward_map_gate":"BLOCKED_NO_COMMITTED_1000UM_M_G_MAP","baseline_forward":"NOT_RUN_GATE_BLOCKED","100um_40_case_campaign":"NOT_RUN_GATE_BLOCKED","1000um_40_case_campaign":"NOT_RUN_GATE_BLOCKED","parameter_classification":"INSUFFICIENT_EVIDENCE","parameter_refit_required":False,"twoD_stochastic_rerun_required":False}
 atomic(o/"oneD_v2_parameter_decision.json",json.dumps(decision,indent=2)+"\n")
 reports(o);figures(o,pf)
 hashes={p.name:sha(p) for p in sorted(o.iterdir()) if p.suffix in {".csv",".json"} and p.name!="oneD_v2_scientific_fingerprints.json"};atomic(o/"oneD_v2_scientific_fingerprints.json",json.dumps(hashes,indent=2,sort_keys=True)+"\n")

def reports(o):
 docs={
 "ONE_D_V2_ARCHITECTURE_AND_LINEAGE.md":"V1 is frozen at 5594253 with audit a5dece9. V2 separates a backend-free kernel, typed mechanics providers, lifecycle policy, stochastic sequence, and analysis. Historical source was not changed.",
 "ONE_D_V2_SHARED_PHYSICS_CONTRACT.md":"The EXP-floor cleavage/emission barrier and Arrhenius hazard algebra are genuinely common. Canonical rows are shared at full precision. PF unified-MPZ and FEM/CZM moving-tip MPZ state evolution are separate source stacks; transport, geometry, backstress, shielding, and renewal remain unqualified rather than falsely unified.",
 "ONE_D_V2_FROZEN_STATE_CONSISTENCY.md":"The V2 kernel returns identical barrier/rate/hazard results for identical local states in both provider configurations by construction. Direct equivalence to both 2-D state engines is blocked until explicit source-state adapters cover retained/mobile fields, multiplicity, backstress and shielding.",
 "ONE_D_V2_EVENT_LIFECYCLE_QUALIFICATION.md":"Typed replay is fail-closed and cannot redraw thresholds or change geometry. Archived PF/FEM ledgers do not contain every pre/post kinetics and RNG state required for exact one-event lifecycle replay; qualification is blocked, not inferred.",
 "ONE_D_V2_TWO_D_REPLAY_QUALIFICATION.md":"Corrected FEM/CZM mechanics transactions and avalanche identities are replayable, but complete kinetics state is absent. PF steps are archived, but an event filter and complete transactional state adapter have not been qualified. Full lifecycle replay therefore does not pass.",
 "ONE_D_V2_FORWARD_BASELINE_VALIDATION.md":"Not run: neither backend has a committed, candidate-independent, actual-operator forward mechanics map covering the requested domain. Replay success would not substitute for this predictive gate.",
 "ONE_D_V2_FOUR_CLASS_CAMPAIGN.md":"The 40-case 100-µm and 40-case 1000-µm campaigns were not launched because replay and forward-map gates did not pass. Empty V2 ledgers carry availability columns; no unsupported values are zero-filled.",
 "ONE_D_V1_V2_COMPARISON.md":"V1 remains exactly reproducible and useful for provenance. V2 preserves its canonical parameter identities and barrier algebra but rejects its event-indexed K/U map as a production mechanics backend. V1 class-topology conclusions remain historical hypotheses until V2 forward qualification passes.",
 "ONE_D_V2_PARAMETER_DECISION.md":"Evidence is insufficient to validate or refit the four rows in V2. No refit is authorized before state/lifecycle and forward mechanics gates pass. No new 2-D stochastic rerun is needed; deterministic mechanics-map generation is the next required work."}
 for n,s in docs.items():atomic(o/n,"# "+n[:-3].replace("_"," ").title()+"\n\n"+s+"\n")

def figures(o,rows):
 import matplotlib;matplotlib.use("Agg");import matplotlib.pyplot as plt
 names=["ONE_D_V2_PF_REPLAY_VS_2D.png","ONE_D_V2_FEMCZM_REPLAY_VS_2D.png","ONE_D_V2_FORWARD_BASELINE_TOPOLOGY.png","ONE_D_V2_FOUR_CLASS_ONSET_ENVELOPES.png","ONE_D_V2_FOUR_CLASS_AVALANCHE_TOPOLOGY.png","ONE_D_V2_PF_VS_FEMCZM_MECHANICS_EFFECT.png","ONE_D_V1_VS_V2_CLASSIFICATION.png","ONE_D_V2_VS_CORRECTED_2D_1000K.png"]
 for n in names:
  plt.figure(figsize=(6,3));plt.axis("off");plt.text(.5,.58,n[:-4].replace("_"," "),ha="center",wrap=True);plt.text(.5,.35,"Qualification gate blocked — no inferred data",ha="center",color="darkred");plt.tight_layout();plt.savefig(o/n,dpi=150);plt.close()
if __name__=="__main__":main()
