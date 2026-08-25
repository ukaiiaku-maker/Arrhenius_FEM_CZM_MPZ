#!/usr/bin/env python3
"""Materialize the fail-closed local-drive audit without launching mechanics."""
from __future__ import annotations
import csv, hashlib, json, subprocess
from pathlib import Path
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
OLD=ROOT/'analysis_outputs/oneD_v2_mechanics_maps_and_baselines'
OUT=ROOT/'analysis_outputs/oneD_v2_local_drive'
PF=Path('/Volumes/Data/Data/Nanopillar_calculation/PF-fracture-fatigue_v10_2_21_persistent_sites_top1')
FEM=Path('/Volumes/Data/Data/Nanopillar_calculation/Arrhenius_FEM_CZM_MPZ_theta0_pf_parity_claude')
RUNS={
 'Peak':PF/'runs/1_final_v10_2_30_theta0_rate_extremes_four_class_1000um_base3621_v1/rate1x/v913_paper_peak01_0242980_persistent_sites/T1000K_th0_seed8666',
 'DBTT':PF/'runs/1_final_v10_2_30_theta0_rate_extremes_four_class_1000um_base3621_v1/rate1x/v913_paper_dbtt01_0202500_persistent_sites/T1000K_th0_seed1008666'}

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(name,x): (OUT/name).write_text(json.dumps(x,indent=2,sort_keys=True)+'\n')
def write(name,text): (OUT/name).write_text(text.strip()+'\n')
def commit(repo): return subprocess.check_output(['git','-C',str(repo),'rev-parse','HEAD'],text=True).strip()

def main():
 OUT.mkdir(parents=True,exist_ok=True)
 old_manifest=json.loads((OLD/'oneD_v2_native_map_inputs_manifest.json').read_text())
 frozen={}
 for key in ('oneD_v2_pf_mechanics_map.csv','oneD_v2_fem_native_mechanics_map.csv','oneD_v2_fem_qualified_G_map.csv','oneD_v2_common_parity_results.csv','oneD_v2_source_shadow_results.csv'):
  p=OLD/key; frozen[key]={'content_sha256':sha(p),'bytes':p.stat().st_size}
 dump('oneD_v2_scalar_map_and_parity_freeze.json',{
  'schema':'oneD_v2_scalar_map_and_parity_freeze_v1','freeze_parent_commit':'8c7a4cb',
  'qualification':'PASS_HASHED_IMMUTABLE_INPUTS','valid_extension_interval_um':[0,1000],
  'level1':'PASS_EXACT','level2':'PASS_EXACT','level3_scalar':'UNQUALIFIED_INSUFFICIENT_DRIVE_DIMENSION',
  'map_manifest_sha256':sha(OLD/'oneD_v2_native_map_inputs_manifest.json'),'artifacts':frozen,
  'source_commits':{k:v['source_commit'] for k,v in old_manifest['maps'].items()},
  'generation_commits':{k:v['generation_commit'] for k,v in old_manifest['maps'].items()}})

 archive={}
 for material,run in RUNS.items():
  p=run/'anisotropic_emission_audit_v10174.json'; x=json.loads(p.read_text()); rows=x['records']
  archive[material]={
   'path':str(p),'sha256':sha(p),'record_count':len(rows),
   'reliable_count':sum(bool(r['anisotropic_drive_reliable']) for r in rows),
   'derived_fields_archived':['anisotropic_tau_signed_Pa','anisotropic_drive_factors','anisotropic_drive_reliable','anisotropic_drive_serial','anisotropic_mechanics_serial'],
   'raw_required_fields_missing':['opening_tensor_Pa','channel_tensors_Pa','probe_element_coordinates_m','probe_element_weights','front_direction','front_normal','tip_xy_m'],
   'first_tau_signed_Pa':rows[0]['anisotropic_tau_signed_Pa'],'last_tau_signed_Pa':rows[-1]['anisotropic_tau_signed_Pa'],
   'first_drive_factors':rows[0]['anisotropic_drive_factors'],'last_drive_factors':rows[-1]['anisotropic_drive_factors']}

 graph={
  'schema':'oneD_v2_pf_tensor_drive_dependency_graph_v1','source_commit_inspected':commit(PF),
  'nodes':[
   {'id':'accepted_2d_mechanics','provides':['mesh','sigma_gp=[sxx,syy,sxy]','damage'],'owner':'assemble_mechanics'},
   {'id':'front_probe','depends_on':['accepted_2d_mechanics','tip_xy','probe config'],'provides':['opening_tensor','two channel tensors','basis','quality'],'owner':'anisotropic_emission_v10174.build_front_drive'},
   {'id':'projection','depends_on':['front_probe','bcc_slip_traces(theta=0)'],'provides':['tau_signed[2]','drive_factors[2]'],'owner':'resolve_channel_drives'},
   {'id':'persistent_emission','depends_on':['projection','reliable','opening stress','dynamic source geometry','backstress','state'],'provides':['signed line content','emission rates'],'owner':'persistent_site_source_v10221._persistent_emit'}],
  'minimum_runtime_state':['reliable','drive_factors[2]','tau_signed_Pa[2]','opening_stress_Pa','persistent process-zone state'],
  'minimum_qualification_state':['opening_tensor_Pa[2,2]','channel_tensors_Pa[2,2,2]','crack-local basis','tip_xy_m','probe coordinates/weights or exact source selection metadata','resolved projections','quality metadata'],
  'archive_audit':archive}
 dump('oneD_v2_pf_tensor_drive_dependency_graph.json',graph)

 schema={
  '$schema':'https://json-schema.org/draft/2020-12/schema','title':'NativeDriveBundle','type':'object',
  'model_kind':'ONE_D_FRONT_MODEL_WITH_LOCAL_TENSOR_DRIVE',
  'required':['global_structural','crack_geometry','local_profile','emission_drive','qualification'],
  'properties':{
   'global_structural':{'required':['opening_m','reaction_N','native_J_J_m2','native_KJ_Pa_sqrt_m']},
   'crack_geometry':{'required':['tip_xy_m','tangent','normal','tip_radius_m','front_width_m']},
   'local_profile':{'required':['coordinates_crack_local_m','stress_tensor_Pa'],'symmetry_required':True},
   'emission_drive':{'required':['opening_tensor_Pa','channel_tensors_Pa','resolved_tau_signed_Pa','drive_factors','reliable']},
   'qualification':{'required':['source_commit','map_sha256','bounds','topology_fingerprint'],'extrapolation':'FORBIDDEN'}}}
 dump('oneD_v2_native_drive_schema.json',schema)

 pf_manifest={'schema':'oneD_v2_pf_tensor_drive_map_manifest_v1','status':'UNQUALIFIED_RAW_SOURCE_FIELDS_NOT_ARCHIVED',
  'candidate_independent_in_principle':True,'map_generated':False,'archive_audit':archive,
  'diagnostics_required':True,'diagnostics_completed':False,'new_2d_pf_diagnostics_launched':0,
  'reason':'Archived derived projections cannot establish raw tensor/profile parity or profile interpolation bounds.'}
 dump('oneD_v2_pf_tensor_drive_map_manifest.json',pf_manifest)
 pd.DataFrame(columns=['extension_m','sample_s_m','sample_n_m','sigma_xx_over_U','sigma_yy_over_U','sigma_xy_over_U','source_quality']).to_parquet(OUT/'oneD_v2_pf_tensor_drive_nodes.parquet',index=False)
 pd.DataFrame([{'quantity':'raw_tensor_profile','status':'NOT_RUN_MAP_UNAVAILABLE','max_abs_error':float('nan'),'sign_mismatch_count':0}]).to_csv(OUT/'oneD_v2_pf_tensor_interpolation_validation.csv',index=False)
 pd.DataFrame([{'material_class':m,'temperature_K':T,'status':'OUTSIDE_MAP_DOMAIN_REQUIRED_FIELD_UNAVAILABLE'} for m in ('Peak','DBTT','weak-T','ceramic-like') for T in (300,1000,1200)]).to_csv(OUT/'oneD_v2_pf_exact_wrapper_shadow.csv',index=False)

 fem_composition={
  'schema':'oneD_v2_fem_factory_composition_v1','source_commit':commit(FEM),
  'factory':'mode_i_first_passage_v10_0_5_14_1_persistent_site_family.main',
  'engine':'PersistentSiteMovingProcessZoneFrontEngineV100514','base_entry':'mode_i_first_passage_v10_0_5_13_5_barrier_only.main',
  'composition':['select persistent-site row','load/validate signed kernel family','apply row to barrier and MPZ options','replace MovingProcessZone2DFrontEngine','augment J wrapper with tensor probes','trial/commit state','renewal/rollback'],
  'instantiation_status':'SOURCE_COMPOSITION_TRACED_NOT_YET_FIELDWISE_FACTORY_SHADOWED','new_fem_2d_runs':0}
 dump('oneD_v2_fem_factory_composition.json',fem_composition)
 fem_drive={'schema':'oneD_v2_fem_drive_requirement_v1','classification':'FULL_LOCAL_TENSOR_REQUIRED',
  'cleavage':'SCALAR_K_TO_SIGMA_TIP','emission':'OPENING_STRESS_PLUS_TWO_SIGNED_TENSOR_CHANNELS',
  'source':['anisotropic_two_channel_drive_v100514.augmented_j_wrapper_factory','persistent_site_front_engine_v100514._two_channel_drive'],
  'fail_closed_message':'persistent signed emission requires a reliable two-channel FEM tensor drive','map_status':'NOT_BUILT_NO_NEW_FEM_2D_RUN_AUTHORIZED'}
 dump('oneD_v2_fem_drive_requirement.json',fem_drive)
 pd.DataFrame(columns=['backend','material_class','event_index','time_s','extension_m','event_length_m','status']).to_csv(OUT/'oneD_v2_level3_native_microtrajectories_v2.csv',index=False)
 decision={'schema':'oneD_v2_level3_decision_v2','Level1':'PASS_EXACT','Level2':'PASS_EXACT','Level3':'NATIVE_CLOSURE_UNQUALIFIED',
  'blockers':['PF raw tensor/profile map unavailable','PF exact tensor shadow not run','FEM exact factory fieldwise shadow incomplete','FEM tensor map unavailable'],
  'common_mode_tensor_parity':'NOT_RUN_PREREQUISITES_FAILED','native_microtrajectories_launched':0,'forward_baselines_authorized':False,
  'canonical_parameters_changed':False,'parameter_status':'INSUFFICIENT_EVIDENCE_FOR_V2_PARAMETER_DECISION','new_2d_PF_diagnostics_launched':0,'new_2d_FEMCZM_runs_launched':0}
 dump('oneD_v2_level3_decision_v2.json',decision)

 write('ONE_D_V2_LOCAL_DRIVE_ARCHITECTURE.md','''# One-dimensional V2 local-drive architecture

Status: **ONE_D_FRONT_MODEL_WITH_LOCAL_TENSOR_DRIVE**.

The scalar crack coordinate remains the only geometric evolution coordinate. Global reaction/J/KJ/G values are retained unchanged, but they are explicitly separated from local kinetic fields. `NativeDriveBundle` carries the crack-local basis, a bounded symmetric tensor profile, source-resolved signed emission channels, geometry, quality, and provenance. Scalar KJ alone fails closed for PF signed emission. No singular-field extrapolation is permitted.''')
 write('ONE_D_V2_PF_TENSOR_DRIVE_CONTRACT.md',f'''# PF tensor-drive contract

The exact source chain is `assemble_mechanics` → `build_front_drive` → `resolve_channel_drives` → `_persistent_emit`. The accepted plane-strain solve exposes Gauss-point `[sxx, syy, sxy]`; damage-filtered, area-weighted probes sample one opening tensor ahead of the inferred front and one tensor along each of the two theta=0 BCC traces. Signed shear is `t @ sigma @ n`. Magnitude factors are `abs(tau)/(0.5 sigma_amplitude)`; sign is retained separately. The persistent source also consumes opening stress, evolving radius/front width/source multiplicity, backstress, and signed MPZ state.

Minimum runtime closure is reliable + two factors + two signed shears. Minimum **qualification** additionally requires the raw tensors, basis, tip and source selection/weight metadata. Archive audit: Peak {archive['Peak']['record_count']} reliable derived rows; DBTT {archive['DBTT']['record_count']}; neither archive contains the raw qualification fields. Fallback is fail-closed in production.''')
 write('ONE_D_V2_PF_TENSOR_MAP_QUALIFICATION.md','''# PF tensor-map qualification

**UNQUALIFIED_RAW_SOURCE_FIELDS_NOT_ARCHIVED.** Existing archives preserve derived signed shears and factors but not the raw opening/channel tensors or probe coordinates/weights. An empty schema-valid parquet is emitted to make the absence machine-readable; it is not a claimed map. The authorized bounded Peak/DBTT observer diagnostics remain required. No interpolation, source shadow, or Level-3 result is claimed.''')
 write('ONE_D_V2_PF_EXACT_WRAPPER_SHADOW.md','''# PF exact-wrapper tensor shadow

**NOT RUN — REQUIRED TENSOR MAP UNAVAILABLE.** The previous scalar attempt remains correctly classified as `UNQUALIFIED_INSUFFICIENT_DRIVE_DIMENSION`. Twelve requested class/temperature cells are recorded as unavailable, not matches.''')
 write('ONE_D_V2_FEMCZM_FACTORY_COMPOSITION_AUDIT.md','''# FEM/CZM factory composition audit

The corrected entry composes the v10.0.5.14.1 persistent-site family wrapper, row adapter, signed kernel-family loader, barrier-only base entry, persistent signed MPZ engine replacement, and augmented tensor-probe J wrapper. This closes the source graph but is not yet a fieldwise isolated factory shadow; it remains fail-closed. No FEM/CZM trajectory or mechanics solve was launched.''')
 write('ONE_D_V2_FEMCZM_DRIVE_REQUIREMENT.md','''# FEM/CZM drive requirement

Classification: **FULL_LOCAL_TENSOR_REQUIRED** for exact native persistent signed emission. Cleavage uses scalar K-to-tip-stress conversion, while emission requires the opening stress plus reliable two-channel tensor projections (signed tau and magnitude factors). Therefore the earlier lower-level scalar shadow was incomplete. No analogous FEM tensor map was generated because new FEM/CZM 2-D mechanics are explicitly prohibited.''')
 write('ONE_D_V2_LEVEL3_NATIVE_CLOSURE_V2.md','''# Level-3 native closure V2

Level 1 and Level 2 remain `PASS_EXACT`. Level 3 remains **NATIVE_CLOSURE_UNQUALIFIED**. Common tensor parity and all four native microtrajectories are gated behind source-qualified PF and FEM local-drive providers and exact wrapper/factory shadows. The remaining discrepancy is currently classified as **MECHANICS_REDUCTION / MISSING_LOCAL_DRIVE**, not a parameter, interpolation, lifecycle, or model-form conclusion. Forward baselines remain prohibited and canonical parameters are unchanged.''')

if __name__=='__main__': main()
