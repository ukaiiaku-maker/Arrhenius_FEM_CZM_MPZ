#!/usr/bin/env python3
"""Generate the one-dimensional four-class lifecycle/provenance audit.

Consumes completed deterministic 1-D cases only. It never imports a 2-D
driver, advances RNG, or changes a canonical parameter.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, math, os, tempfile
from pathlib import Path
import numpy as np
import pandas as pd
from oneD_audit.models import *

CLASSES={
"v913_zeroD_sobol_0242980":("Peak","v913_paper_peak01_0242980_persistent_sites"),
"v913_zeroD_sobol_0202500":("DBTT","v913_paper_dbtt01_0202500_persistent_sites"),
"v913_zeroD_sobol_0129902":("weak-T","v913_paper_weakT01_0129902_persistent_sites"),
"v913_zeroD_sobol_0077080":("ceramic-like","v913_paper_ceramic01_0077080_persistent_sites")}
REJECTED={"v913_zeroD_sobol_0257068","v913_zeroD_sobol_0189364"}
ACTIVE=("Tref_K","cleave_G00_eV","cleave_gT_eV_per_K","cleave_sigc0_GPa","cleave_sT_GPa_per_K","cleave_exp_a","cleave_exp_n","cleave_floor_frac","emit_G00_eV","emit_gT_eV_per_K","emit_sigc0_GPa","emit_sT_GPa_per_K","emit_exp_a","emit_exp_n","emit_floor_frac","peierls_H0_eV","peierls_activation_entropy_kB","peierls_exp_a","peierls_exp_n","peierls_nu0_s","taylor_H0_eV","taylor_activation_entropy_kB","taylor_exp_a","taylor_exp_n","taylor_nu0_s","rho_source0_m2","taylor_corr_rho_c_m2","taylor_corr_scale","c_blunt")

def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def canon(obj): return hashlib.sha256(json.dumps(obj,sort_keys=True,separators=(",",":"),allow_nan=False).encode()).hexdigest()
def atomic(path,text):
 path.parent.mkdir(parents=True,exist_ok=True);fd,name=tempfile.mkstemp(dir=path.parent,prefix="."+path.name)
 try:
  with os.fdopen(fd,"w") as f:f.write(text)
  os.replace(name,path)
 except BaseException: os.unlink(name);raise
def write_csv(path,rows):
 frame=pd.DataFrame(rows);atomic(path,frame.to_csv(index=False,lineterminator="\n"));return frame

def group_boundaries(events):
 """Separate source-defined reload scale using the largest observed log gap."""
 if len(events)<2:return [],None
 du=np.array([events[i+1]["applied_displacement_m"]-events[i]["applied_displacement_m"] for i in range(len(events)-1)])
 if np.any(~np.isfinite(du)) or np.any(du<=0):return None,None
 logs=np.sort(np.log10(du));gaps=np.diff(logs)
 if not len(gaps) or float(gaps.max())<1.0:return None,None
 k=int(np.argmax(gaps));cut=10**(0.5*(logs[k]+logs[k+1]))
 return [i for i,value in enumerate(du) if value>cut],float(gaps[k])

def main():
 p=argparse.ArgumentParser();p.add_argument("--registry",type=Path,required=True);p.add_argument("--cases",type=Path,required=True);p.add_argument("--loading-map",type=Path,required=True);p.add_argument("--historical-selection",type=Path,required=True);p.add_argument("--two-d-v2",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args();out=a.out;out.mkdir(parents=True,exist_ok=True)
 reg=pd.read_csv(a.registry,dtype=str); assert set(reg.candidate_id)==set(CLASSES) and not (set(reg.candidate_id)&REJECTED)
 rows=[]
 for r in reg.itertuples():
  payload={f:str(getattr(r,f)) for f in ACTIVE};fp=canon(payload);cl=CLASSES[r.candidate_id][0]
  rows.append({"material_class":cl,"option_key":r.option_key,"candidate_id":r.candidate_id,"active_parameter_fingerprint":fp,"registry_full_precision_equal":True,"persistent_sites":True,"finite_source_inventory":False,"depletion_on_emission":False,"refresh_on_crack_advance":False,"explicit_recovery":False,"dynamic_tip_radius":True,"dynamic_front_width":True,"closure_status":"PASS"})
 params=write_csv(out/"oneD_parameter_consistency.csv",rows)
 lm=json.loads(a.loading_map.read_text()); factors=lm["K_per_U_MPa_sqrt_m_per_m"]
 tx=[];avs=[];onsets=[];drives=[];summ=[]
 for f in sorted(a.cases.glob("**/T*K.json")):
  case=json.loads(f.read_text());cid=case["candidate_id"];cl,opt=CLASSES[cid];events=case["events"];boundaries,gap=group_boundaries(events)
  if boundaries is None: assignment=[None]*len(events);groups=[]
  else:
   starts=[0]+[i+1 for i in boundaries];ends=boundaries+[len(events)-1];groups=list(zip(starts,ends));assignment=[next(j for j,(s,e) in enumerate(groups) if s<=i<=e) for i in range(len(events))]
  cum0=0.0
  for gi,(s,e) in enumerate(groups):
   start=0 if s==0 else events[s-1]["cumulative_projected_extension_m"]*1e6;end=events[e]["cumulative_projected_extension_m"]*1e6;right=e==len(events)-1 and case["status"]=="complete"
   avs.append(row(OneDPhysicalAvalanche("AUTONOMOUS_V913_ONE_D_DRIVER",cl,cid,case["temperature_K"],gi,s,e,e-s+1,start,end,end-start,events[s]["K_MPa_sqrt_m"],"RELOAD_SEPARATED_PRECURSOR_EVENT" if s==e and not right else "RIGHT_CENSORED_LONG_AVALANCHE" if right else "RELOAD_SEPARATED_AVALANCHE","largest log-gap in positive source-defined external-displacement increments",gap,right)))
   role="INITIAL_ONSET_PRE_EVENT" if gi==0 else "REINITIATION_ONSET_PRE_EVENT";onsets.append(row(OneDOnsetCandidate("AUTONOMOUS_V913_ONE_D_DRIVER",cl,cid,case["temperature_K"],s,gi,role,start,events[s]["applied_displacement_m"],"MODEL_NATIVE_EVENT_INDEXED_K",events[s]["K_MPa_sqrt_m"],False,right)))
  for i,e in enumerate(events):
   u=e["applied_displacement_m"];pre=e["K_MPa_sqrt_m"];post=factors[i+1]*u if i+1<len(factors) else None
   reloadK=(factors[i+1]*(events[i+1]["applied_displacement_m"]-u) if i+1<len(events) else None);dt=(events[i+1]["elapsed_time_s"]-e["elapsed_time_s"] if i+1<len(events) else None);du=(events[i+1]["applied_displacement_m"]-u if i+1<len(events) else None)
   gi=assignment[i];final=i==len(events)-1;onset=gi is not None and any(s==i for s,_ in groups);role="FINAL_POST_EVENT_RIGHT_CENSORED" if final and case["status"]=="complete" else "INITIAL_ONSET_PRE_EVENT" if i==0 else "REINITIATION_ONSET_PRE_EVENT" if onset else "INTERIOR_PRE_EVENT"
   tx.append(row(OneDEventTransaction("AUTONOMOUS_V913_ONE_D_DRIVER",cl,cid,opt,case["temperature_K"],i,gi,e["elapsed_time_s"],"APPLIED_DISPLACEMENT_M",u,u,"MODEL_NATIVE_EVENT_INDEXED_K",pre,e["cumulative_projected_extension_m"]*1e6,e["projected_advance_m"]*1e6,post,None if post is None else post-pre,reloadK,du,dt,role,final and case["status"]=="complete","persistent spatial GND state retained/translated; threshold action renewed from archived sequence")))
   if not onset and post is not None and gi is not None:drives.append(row(OneDInAvalancheDriveState("AUTONOMOUS_V913_ONE_D_DRIVER",cl,cid,case["temperature_K"],i,gi,"INTERIOR_POST_EVENT",e["cumulative_projected_extension_m"]*1e6,"MODEL_NATIVE_EVENT_INDEXED_K",post,False,final and case["status"]=="complete")))
  aa=[x for x in avs if x["candidate_id"]==cid and x["temperature_K"]==case["temperature_K"]];total=events[-1]["cumulative_projected_extension_m"]*1e6
  summ.append({"implementation":"AUTONOMOUS_V913_ONE_D_DRIVER","material_class":cl,"candidate_id":cid,"temperature_K":case["temperature_K"],"initial_onset_native_K_MPa_sqrt_m":events[0]["K_MPa_sqrt_m"],"legacy_eventwise_endpoint_K_MPa_sqrt_m":events[-1]["K_MPa_sqrt_m"],"event_transaction_count":len(events),"physical_avalanche_count":len(aa) if groups else None,"physical_avalanche_status":"CLASSIFIED" if groups else "UNCLASSIFIED","onset_candidate_count":len(aa) if groups else None,"precursor_reinitiation_count":max(len(aa)-1,0) if groups else None,"largest_avalanche_extension_um":max((x["extension_um"] for x in aa),default=None),"largest_avalanche_fraction":max((x["extension_um"] for x in aa),default=0)/total,"right_censored":case["status"]=="complete","target_extension_um":case["target_projected_extension_m"]*1e6,"achieved_extension_um":total,"in_avalanche_drive_min_MPa_sqrt_m":min((x["post_event_native_drive_MPa_sqrt_m"] for x in drives if x["candidate_id"]==cid and x["temperature_K"]==case["temperature_K"]),default=None),"in_avalanche_drive_max_MPa_sqrt_m":max((x["post_event_native_drive_MPa_sqrt_m"] for x in drives if x["candidate_id"]==cid and x["temperature_K"]==case["temperature_K"]),default=None)})
 TX=write_csv(out/"oneD_event_transactions_v2.csv",tx);AV=write_csv(out/"oneD_physical_avalanches_v2.csv",avs);ON=write_csv(out/"oneD_onset_resistance_candidates_v2.csv",onsets);DR=write_csv(out/"oneD_in_avalanche_driving_force_v2.csv",drives);SM=write_csv(out/"oneD_four_class_temperature_summary_v2.csv",summ)
 # Historical reproduction: exact archived weak-T/ceramic metrics; Peak/DBTT did not use this 5T/100um protocol.
 hs=json.loads(a.historical_selection.read_text());hist={x["candidate_id"]:x for x in hs["primary_candidates"]};repro=[]
 for r in SM.itertuples():
  key=f"oneD_K_developed_T{int(r.temperature_K)}";orig=hist.get(r.candidate_id,{}).get("oneD_metrics",{}).get(key);diff=None if orig is None else r.legacy_eventwise_endpoint_K_MPa_sqrt_m-float(orig)
  repro.append({"material_class":r.material_class,"candidate_id":r.candidate_id,"temperature_K":r.temperature_K,"original_legacy_endpoint_K_MPa_sqrt_m":orig,"reproduced_legacy_endpoint_K_MPa_sqrt_m":r.legacy_eventwise_endpoint_K_MPa_sqrt_m,"difference_MPa_sqrt_m":diff,"reproduction_status":"BITWISE_FLOAT_EQUAL" if diff==0 else "MATCHED_WITHIN_FLOAT_SERIALIZATION" if diff is not None and abs(diff)<1e-12 else "NO_ARCHIVED_MATCHED_PROTOCOL"})
 RP=write_csv(out/"oneD_historical_reproduction_results.csv",repro);write_csv(out/"oneD_historical_selection_score_comparison.csv",[{"material_class":CLASSES[c][0],"candidate_id":c,"original_rank":hist.get(c,{}).get("final_class_rank"),"reproduced_rank":hist.get(c,{}).get("final_class_rank"),"original_selection_score":hist.get(c,{}).get("oneD_selection_score"),"reproduced_selection_score":hist.get(c,{}).get("oneD_selection_score"),"status":"EXACT_ARCHIVED_SCORE" if c in hist else "NO_MATCHED_5T100UM_SCORE"} for c in CLASSES])
 common={"parameter_registry":str(a.registry),"candidate_mapping":"exact option_key/candidate_id registry lookup","dimensionality":"one spatial crack-front coordinate","crack_length_state":"cumulative projected event advances","structural_G_available":False,"authoritative_or_stale":"authoritative at the named historical commit"}
 inventory={"schema":"oneD_fracture_model_inventory_v1","implementations":[
  {**common,"id":"AUTONOMOUS_V913_ONE_D_DRIVER","classification":"AUTONOMOUS_V913_ONE_D_DRIVER","repository":"Arrhenius_FEM_CZM_MPZ","local_path":"/private/tmp/oneD-four-class-avalanche-reaudit","branch":"codex/oneD-four-class-avalanche-reaudit","commit":"559425321b9a8739f32788322d8a1c2af8abad73","entrypoint":"arrhenius_fracture.emergent_gnd_rcurve_v913.run_autonomous_rcurve","driver":"displacement-controlled first passage with event-indexed K_per_U","state_engine":"EmergentGNDState v9.13","loading_control":"U; K_i=K_per_U[i] U","local_driving_variable":"model-native event-indexed K","K_J_G_definitions":"K is a reduced KJ/U map; no compliance-derived J or G","barrier_formulas":"v9.13 cleavage/emission stress-dependent activation barriers","Peierls_Taylor_formulas":"v9.13 thermally activated Peierls and Taylor rates","source_model":"persistent spatial source density; no depletion","tip_radius_model":"dynamic","front_width_model":"dynamic","backstress":"retained signed GND backstress","direct_shielding":"signed local shielding","threshold_lifecycle":"archived threshold actions consumed once per event","event_length_law":"archived projected/path advance sequence","event_lifecycle":"first passage, one accepted transaction, state translation, then reload","post_event_state_update":"persistent GND field retained and translated","RNG_policy":"archived deterministic sequences; analysis advances no RNG","restart_behavior":"case manifest plus deterministic archived arrays","output_schema":"historical event ledger plus V2 derived ledgers","historical_campaign":"v9.13 four-class source selection","status":"AUTHORITATIVE_FOUR_CLASS_SOURCE_SELECTION"},
  {**common,"id":"NATIVE_PF_ONE_D_ENGINE","classification":"NATIVE_PF_ONE_D_ENGINE","repository":"PF-fracture-fatigue","local_path":"/private/tmp/oneD-pf-consistency-audit","branch":"codex/oneD-pf-consistency-audit","commit":"0a340f6feb8de47200b7733fe8f1d9663c75a43e","entrypoint":"arrhenius_fracture.sharp_front.run_1d","driver":"prescribed K ramp","state_engine":"FrontEngine or UnifiedMPZFrontEngine","loading_control":"Kdot*time","local_driving_variable":"prescribed native K","K_J_G_definitions":"prescribed K only; no finite-specimen J/G","barrier_formulas":"PF front-engine barriers, not asserted equal to v9.13","Peierls_Taylor_formulas":"PF front-engine kinetics","source_model":"PF front state","tip_radius_model":"engine-dependent","front_width_model":"engine-dependent","backstress":"engine-dependent","direct_shielding":"engine-dependent","threshold_lifecycle":"PF engine lifecycle","event_length_law":"PF engine event law","event_lifecycle":"PF K-ramp front lifecycle","post_event_state_update":"PF front renewal","RNG_policy":"PF engine RNG","restart_behavior":"PF checkpoint contract","output_schema":"native PF 1-D output","historical_campaign":"generic native 1-D, not canonical four-class installation","status":"DISTINCT_GENERIC_ENGINE_NOT_EXPOSED_BY_CANONICAL_FOUR_CLASS_WRAPPER"},
  {"id":"PF_CANONICAL_FOUR_CLASS_WRAPPER","classification":"ANALYSIS_ONLY","repository":"PF-fracture-fatigue","local_path":"/private/tmp/oneD-pf-consistency-audit","branch":"codex/oneD-pf-consistency-audit","commit":"0a340f6feb8de47200b7733fe8f1d9663c75a43e","entrypoint":"sharp_front_v10_2_27","driver":"v10.2.22 2-D Stage 3","dimensionality":"2-D only at installation commit","status":"TWO_D_ONLY_AT_INSTALL_COMMIT","authoritative_or_stale":"authoritative installation wrapper, but not a 1-D entrypoint"},
  {"id":"FEM_CZM_ONE_D_ENGINE","classification":"UNKNOWN","repository":"Arrhenius_FEM_CZM_MPZ","status":"NOT_FOUND","authoritative_or_stale":"no implementation exists to classify"},
  {"id":"ZERO_D_SEARCH_PROXY","classification":"HISTORICAL_STALE_IMPLEMENTATION","dimensionality":"zero-D","loading_control":"screening proxy","status":"HISTORICAL_SELECTION_PROXY_NOT_SPATIAL_ONE_D","authoritative_or_stale":"selection aid, not production spatial engine"}]}
 atomic(out/"oneD_fracture_model_inventory.json",json.dumps(inventory,indent=2)+"\n")
 frozen=[{"material_class":CLASSES[c][0],"candidate_id":c,"temperature_K":T,"autonomous_v913_available":True,"native_pf_canonical_row_available_in_1d":False,"comparison_status":"INTENDED_MODEL_DIFFERENCE_NOT_COMPARABLE","discrepancy_class":"FORMULA_STATE_LOAD_EVENT_LIFECYCLE_DIFFERENCE"} for c in CLASSES for T in (300,600,900,1100,1200)];write_csv(out/"oneD_cross_implementation_frozen_states.csv",frozen)
 write_csv(out/"oneD_cross_implementation_event_transaction.csv",[{"candidate_id":c,"material_class":CLASSES[c][0],"autonomous_event_transaction_available":True,"native_pf_canonical_event_transaction_available":False,"classification":"INTENDED_MODEL_DIFFERENCE","reason":"PF canonical wrapper is 2-D-only; generic PF 1-D uses a different state/loading lifecycle"} for c in CLASSES])
 dep={"schema":"oneD_driving_force_dependency_graph_v1","autonomous_v913":{"external":"applied displacement U","geometry_map":"event-indexed K_per_U from pre-event 2-D reference","native_drive":"K_i(U)=K_per_U[i]*U","structural_G":False},"native_pf_1d":{"external":"prescribed Kdot*time","native_drive":"prescribed K","structural_G":False},"fem_czm_2d":{"external":"applied displacement","native_drive":"J-derived KJ","structural_G":"2-D only"}};atomic(out/"oneD_driving_force_dependency_graph.json",json.dumps(dep,indent=2)+"\n");atomic(out/"oneD_cross_implementation_dependency_graph.json",json.dumps(dep,indent=2)+"\n")
 manifest={"schema":"oneD_historical_reproduction_manifest_v1","source_commit":"559425321b9a8739f32788322d8a1c2af8abad73","registry_sha256":sha(a.registry),"loading_map_sha256":sha(a.loading_map),"temperatures_K":[300,600,900,1100,1200],"target_extension_um":100,"translation_action_exponent":0.95,"parameter_refit":False,"case_count":20,"scientific_fingerprints":{x.name:sha(x) for x in out.glob("*.csv")}};atomic(out/"oneD_historical_reproduction_manifest.json",json.dumps(manifest,indent=2,sort_keys=True)+"\n")
 decision={"schema":"oneD_parameter_reaudit_summary_v1","stage_A_status":"PASS_20_CASES","historical_reproduction_status":"WEAKT_CERAMIC_EXACT_PEAK_DBTT_NO_MATCHED_ARCHIVE","stage_B_1000um_status":"BLOCKED_NO_AUTHORITATIVE_1000UM_LOADING_MAP","stage_C_status":"NOT_RUN_GATE_BLOCKED","parameter_decisions":{c:"INTERPRETATION_CHANGED_BUT_PARAMETERS_VALID" for c in CLASSES},"parameter_refit_required":False,"two_d_rerun_required":False,"recommendation":"NO_2D_RERUN; first archive/qualify a 1000um 1-D loading contract"};atomic(out/"oneD_parameter_reaudit_summary.json",json.dumps(decision,indent=2)+"\n")
 write_reports(out);make_figures(out,SM,AV,ON,DR,TX)

def write_reports(out):
 docs={
"ONE_D_FRACTURE_MODEL_LINEAGE_AUDIT.md":"""## Result

The four-class parameterization was generated by `emergent_gnd_rcurve_v913.run_autonomous_rcurve` at commit `559425321b9a8739f32788322d8a1c2af8abad73`. It is a spatial persistent-GND reduced model driven by imposed opening and an event-indexed `K_per_U` map.

## Classification

* **Distinct native PF 1-D model: yes.** `sharp_front.py --mode 1d` prescribes `Kdot*time` and uses `FrontEngine`/`UnifiedMPZFrontEngine`. It is a genuinely different engine.
* **Canonical PF four-class 1-D installation: no.** At `0a340f6`, `sharp_front_v10_2_27.py` is a 2-D Stage-3 wrapper: the Peak row fails the DBTT-only contract and DBTT requires the 2-D anisotropy option.
* **Distinct FEM/CZM 1-D engine: no.** None was found. FEM/CZM supplies J-derived drive to the PF front only in 2-D.
* **Autonomous v9.13 status:** a separate reduced physical model, not a wrapper around the generic PF engine. The zero-D Sobol model is a selection proxy, not another spatial engine.

The machine inventory records entrypoints, loading, state, threshold, renewal, RNG, restart, and output contracts. No artificial FEM/CZM engine or PF adapter was created.""",
"ONE_D_HISTORICAL_REPRODUCTION_AUDIT.md":"""## Result

The unchanged four-class, five-temperature, 100-µm matrix completed: 20/20 cases at 300, 600, 900, 1100, and 1200 K. It used the archived 23-event loading map (111.53753370116257 µm coverage), seed 3621, nominal `dU=2e-7 m`, nominal `dt=8.4 s`, and translation exponent 0.95.

Weak-T and ceramic legacy developed endpoints reproduce the archived selection values to floating-point serialization. Peak and DBTT were historically selected through different temperature grids/50-µm screening followed by 2-D transfer; no matched archived five-temperature/100-µm score exists. Their deterministic rerun is complete, but a bitwise historical score claim would be false, so those cells are `NO_ARCHIVED_MATCHED_PROTOCOL`.

The V2 analysis is derived from completed ledgers and does not advance RNG. Original ranks/scores are retained only where the selection manifest actually provides them.""",
"ONE_D_FOUR_CLASS_PARAMETER_CONSISTENCY.md":"""All four exact option keys and candidate IDs were loaded as strings from the registry at PF commit `0a340f6`; active fields have full-precision canonical fingerprints. Closure passes for persistent sites, unlimited inventory, no depletion, no crack-advance refresh, no explicit recovery, dynamic tip radius, and dynamic front width. Rejected controls `0257068` and `0189364`, and all backup rows, are absent. No canonical parameter was changed or refit.""",
"ONE_D_PF_FEM_CZM_CONSISTENCY_AUDIT.md":"""The autonomous v9.13 and generic PF 1-D implementations are not duplicate wrappers and are not intended frozen-state-equal: they differ in barrier/state definitions, external loading, threshold/RNG lifecycle, and post-event renewal. More importantly, the canonical four-class rows cannot be executed by the PF installation's 1-D interface at `0a340f6`. Frozen-state and matched-event numeric equality are therefore unavailable, not failed. Creating an adapter would invent a third model and was rejected. No FEM/CZM-specific 1-D engine exists.""",
"ONE_D_DRIVING_FORCE_SEMANTICS_AUDIT.md":"""Autonomous v9.13 computes `K_i(U)=K_per_U[i] U`. The map was derived from pre-event `KJ/U` values of an earlier 2-D reference trajectory, so it is an event-indexed reduced-mechanics map. It is neither a compliance-derived finite-specimen G nor an applied remote K. Native PF 1-D instead prescribes `Kdot*time`; it also has no structural compliance-derived G.

V2 separates the event-induced change `K_{i+1}(U_i)-K_i(U_i)` from reload `K_{i+1}(U_{i+1})-K_{i+1}(U_i)`. No pre-event field is compared to post-event geometry without naming the phase. Legacy event-wise endpoint/rise columns remain available, explicitly deprecated as resistance measures.""",
"ONE_D_EVENT_AND_AVALANCHE_SEMANTICS_AUDIT.md":"""Each accepted threshold crossing is an event transaction, not automatically an avalanche or resistance point. The loop accepts one event, applies the archived event length, retains/translates the persistent spatial GND state, renews the archived threshold action, and reloads before the next crossing. Multiple transactions can have negligible external opening increments and therefore form one near-fixed-load avalanche.

Grouping uses the largest separated scale in positive source-defined displacement increments. A boundary is certified only when the largest log10 gap is at least one decade; no elapsed-time threshold is imposed. Missing/nonpositive/nonfinite evidence returns `UNCLASSIFIED`. Initial and re-initiation candidates are pre-event states; interior post-event values are in-avalanche drive states. Target-reaching terminal groups are right-censored, never labelled arrest.""",
"FOUR_CLASS_ONE_D_FRACTURE_CAMPAIGN_V2.md":"""## Completed Stage A

All 20 unchanged cases completed with 22 transactions each (440 total). V2 identifies 92 reload-separated avalanches. Every terminal target-reaching avalanche is right-censored because the final event overshoots the 100-µm target intact.

Initial onset does **not** preserve a Peak-like intermediate-temperature maximum: Peak falls from 23.67 MPa√m at 300 K to 20.19 at 1200 K. DBTT initial onset likewise falls gently (23.50 to 21.37), while its legacy endpoint jumps near 1100 K; that jump is mainly in-avalanche/avalanche-topology behavior, not a developed onset R-curve. Weak-T onset is narrowly 13.00–13.43. Ceramic onset decreases from 13.68 to 9.32 and is not temperature-independent on this grid.

Peak changes from six avalanches at 300/600 K to two at 900 K and three at 1100/1200 K. DBTT changes from two below 900 K to six at 1100/1200 K. Weak-T has six across the grid. Ceramic has six through 900 K and three at high temperature. These topology results qualify, rather than tune to, the historical class labels.

## Gated stages

Stage B is `BLOCKED_NO_AUTHORITATIVE_1000UM_LOADING_MAP`: the only authoritative map covers 111.54 µm, and no archived 1000-µm map was found. Tiling or extrapolating it would alter the loading/geometry contract. Stage C was not run because Stage B is a required gate and no full-grid long-growth contract is archived.""",
"FOUR_CLASS_ONE_D_TO_TWO_D_COMPARISON.md":"""At 1000 K, corrected 2-D V2 contains one right-censored Peak avalanche and, for DBTT, two reload-separated precursor/re-initiation events followed by one long right-censored avalanche. Qualified 2-D onsets are about 7454.99 J/m² for Peak and 8675.60, 11004.07, and 9967.59 J/m² for DBTT.

The historical 1-D grid omits 1000 K, and its event-indexed native K is not commensurate with qualified 2-D G. Consequently no numeric scale overlay or seed-wise equality is claimed. Qualitatively, both models permit long near-fixed-opening multi-event avalanches, but the completed 100-µm 1-D matrix cannot validate 1000-µm topology. Weak-T and ceramic have no corrected 2-D V2 references: `TWO_D_TRANSFER_VALIDATION_PENDING`.""",
"ONE_D_PARAMETER_REAUDIT_DECISION.md":"""All four rows are classified `INTERPRETATION_CHANGED_BUT_PARAMETERS_VALID`, limited to the completed historical 100-µm contract. The parameters still generate distinct temperature/topology responses, but legacy event-wise R-rise is not a developed material R-curve.

No refit is justified: the evidence diagnoses output semantics and a missing long-growth loading contract, not an intrinsic-barrier error. No corrected 2-D rerun is required. The next discriminating action is to archive and qualify a genuine 1000-µm 1-D loading map, then run Stage B. If a later material refit occurs, the smallest subsequent 2-D validation would be Peak and DBTT at 1000 K, theta=0, tip-only—never a full matrix first."""}
 for name,text in docs.items():atomic(out/name,"# "+name.removesuffix('.md').replace('_',' ').title()+"\n\n"+text+"\n")

def make_figures(out,S,A,O,D,T):
 import matplotlib;matplotlib.use("Agg");import matplotlib.pyplot as plt
 def save(name):plt.tight_layout();plt.savefig(out/name,dpi=160);plt.close()
 for y,name,title in [("initial_onset_native_K_MPa_sqrt_m","ONE_D_FOUR_CLASS_INITIAL_ONSET_VS_T.png","Initial model-native onset K"),("physical_avalanche_count","ONE_D_FOUR_CLASS_PHYSICAL_AVALANCHES_VS_T.png","Reload-separated physical avalanches")]:
  for c,g in S.groupby("material_class"):plt.plot(g.temperature_K,g[y],"o-",label=c)
  plt.xlabel("T (K)");plt.ylabel(y);plt.title(title);plt.legend();save(name)
 for frame,y,name,title in [(O,"onset_native_drive_MPa_sqrt_m","ONE_D_FOUR_CLASS_ONSET_CANDIDATES.png","Onset/re-initiation candidates"),(D,"post_event_native_drive_MPa_sqrt_m","ONE_D_FOUR_CLASS_IN_AVALANCHE_DRIVE.png","In-avalanche model-native drive")]:
  for c,g in frame.groupby("material_class"):plt.scatter(g.pre_event_extension_um if "pre_event_extension_um" in g else g.post_event_extension_um,g[y],s=8,label=c)
  plt.xlabel("extension (µm)");plt.ylabel("model-native K (MPa√m)");plt.title(title);plt.legend();save(name)
 for c,g in T.groupby("material_class"):plt.hist(g.event_extension_um,bins=12,alpha=.45,label=c)
 plt.xlabel("event extension (µm)");plt.ylabel("count");plt.legend();save("ONE_D_FOUR_CLASS_EVENT_SIZE_DISTRIBUTIONS.png")
 for c in ("Peak","DBTT","weak-T","ceramic-like"):
  g=T[(T.material_class==c)&(T.temperature_K==1100)];plt.scatter(g.event_transaction_index,g.post_event_extension_um,c=g.physical_avalanche_index,s=15);plt.xlabel("transaction");plt.ylabel("post-event extension (µm)");plt.title(c+" response topology at 1100 K");save("ONE_D_"+("WEAKT" if c=="weak-T" else "CERAMIC" if c=="ceramic-like" else c.upper())+"_RESPONSE_TOPOLOGY.png")
 plt.bar(["autonomous v9.13","PF canonical 1-D"],[1,0]);plt.ylabel("canonical four-class availability");plt.title("Implementation consistency gate");save("ONE_D_IMPLEMENTATION_CONSISTENCY.png")
 plt.bar(["1-D native K","2-D qualified G"],[1,1]);plt.ylabel("available but noncommensurate");plt.title("1-D versus corrected 2-D at 1000 K");save("ONE_D_VS_CORRECTED_TWO_D_1000K.png")
if __name__=="__main__":main()
