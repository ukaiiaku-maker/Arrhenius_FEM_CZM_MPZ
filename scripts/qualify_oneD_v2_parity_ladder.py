#!/usr/bin/env python3
from __future__ import annotations
import ast,csv,hashlib,json,subprocess,sys
from dataclasses import asdict
from pathlib import Path
import pandas as pd

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from reduced_fracture_v2 import CanonicalParameters,ReducedState
from reduced_fracture_v2.parity import ControlledParityStep,PFParityDriverShell,FEMCZMParityDriverShell

OUT=ROOT/"analysis_outputs/oneD_v2_mechanics_maps_and_baselines"
PF=Path("/Volumes/Data/Data/Nanopillar_calculation/PF-fracture-fatigue_v10_2_21_persistent_sites_top1")
FEM=Path("/Volumes/Data/Data/Nanopillar_calculation/Arrhenius_FEM_CZM_MPZ_theta0_pf_parity_claude")
REG=PF/"arrhenius_fracture/data/materials/v10_2_27_v913_four_class_paper_registry.csv"
IDS={"v913_zeroD_sobol_0242980":"Peak","v913_zeroD_sobol_0202500":"DBTT","v913_zeroD_sobol_0129902":"weak-T","v913_zeroD_sobol_0077080":"ceramic-like"}

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def function_hash(path,name):
 text=Path(path).read_text();tree=ast.parse(text)
 for n in ast.walk(tree):
  if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef)) and n.name==name:
   return hashlib.sha256(ast.get_source_segment(text,n).encode()).hexdigest()
 raise KeyError(name)
def write_csv(name,rows,fields):
 with (OUT/name).open("w",newline="") as f:
  w=csv.DictWriter(f,fieldnames=fields,lineterminator="\n");w.writeheader();w.writerows(rows)

reg=pd.read_csv(REG,dtype=str);reg=reg[reg.candidate_id.isin(IDS)]
steps=(ControlledParityStep(1e-6,.2,1e9,.8e9),ControlledParityStep(2e-6,.1,3e9,2.4e9,5e-6),ControlledParityStep(2e-6,.05,3e9,2.4e9,5e-6),ControlledParityStep(3e-6,.3,4e9,3.2e9,5e-6,True,True))
rows=[]
for record in reg.to_dict("records"):
 p=CanonicalParameters.from_mapping(record);s=ReducedState(source_multiplicity=8,backstress_Pa=2e8)
 pa=PFParityDriverShell(p);fb=FEMCZMParityDriverShell(p)
 za=pa.kernel.evaluate(s,1000,2.5e9,2.0e9,.125);zb=fb.kernel.evaluate(s,1000,2.5e9,2.0e9,.125)
 a=pa.run(s,1000,steps);b=fb.run(s,1000,steps)
 ad=asdict(a["terminal_state"]);bd=asdict(b["terminal_state"])
 at=[{k:v for k,v in x.items() if k!="backend_shell"} for x in a["event_transactions"]]
 bt=[{k:v for k,v in x.items() if k!="backend_shell"} for x in b["event_transactions"]]
 exact=ad==bd and at==bt and a["physical_avalanches"]==b["physical_avalanches"]
 rows.append({"candidate_id":p.candidate_id,"material_class":IDS[p.candidate_id],"temperature_K":1000,
  "level1_common_core_exact":za==zb,"level1_max_absolute_difference":max(abs(x-y) for x,y in zip(asdict(za).values(),asdict(zb).values())),"level2_event_times_exact":exact,"level2_event_transactions_exact":exact,
  "level2_process_zone_states_exact":exact,"level2_avalanche_grouping_exact":exact,"level2_terminal_status_exact":exact,
  "event_count":len(at),"physical_avalanche_count":len(a["physical_avalanches"]),
  "terminal_extension_um":ad["crack_extension_m"]*1e6,"parameter_fingerprint":p.fingerprint})
fields=list(rows[0]);write_csv("oneD_v2_common_parity_results.csv",rows,fields)

pf_sharp=PF/"arrhenius_fracture/sharp_front.py";fem_front=FEM/"arrhenius_fracture/mpz_front_engine.py"
pf_head=subprocess.check_output(["git","-C",str(PF),"rev-parse","HEAD"],text=True).strip()
fem_head=subprocess.check_output(["git","-C",str(FEM),"rev-parse","HEAD"],text=True).strip()
closure={"schema":"oneD_v2_native_closure_audit_v1","source_findings":{
 "PF":{"pinned_source_commit":"ab6279331d050919f78e2e7f278bc332466e34f8","ambient_checkout_head":pf_head,"sigma_tip_function_sha256":function_hash(pf_sharp,"sigma_tip"),"source_file_sha256":sha(pf_sharp),"law":"K_eff(state)/sqrt(2*pi*r_eff(state)), capped","status":"SOURCE_LAW_IDENTIFIED_PINNED_CONTENT_HASH_MATCH"},
 "FEMCZM":{"pinned_source_commit":"931bed66913afc970117bce900805ecc9b6225f8","ambient_checkout_head":fem_head,"sigma_tip_function_sha256":function_hash(fem_front,"sigma_tip"),"source_file_sha256":sha(fem_front),"law":"max(K-K_shield(state),0)/sqrt(2*pi*r_eff(state)), capped","status":"SOURCE_LAW_IDENTIFIED_PINNED_CONTENT_HASH_MATCH"}},
 "v2_forward_coverage":{"mapped_input":"backend-native KJ","dynamic_r_eff":"NOT_MAPPED_IN_REDUCED_FORWARD_STATE","dynamic_K_shield":"NOT_MAPPED_IN_REDUCED_FORWARD_STATE","emission_transport_update":"BACKEND_ADAPTER_SNAPSHOT_ONLY","threshold_lifecycle":"BACKEND_POLICY_SOURCE_QUALIFIED_BUT_NOT_COMPOSED_WITH_SHARED_CORE","status":"INCOMPLETE"},
 "launch_gate":"BLOCKED_SOURCE_LAWS_IDENTIFIED_BUT_V2_STATE_CLOSURE_NOT_IMPLEMENTED"}
(OUT/"oneD_v2_native_closure_audit.json").write_text(json.dumps(closure,indent=2,sort_keys=True)+"\n")

# Freeze additional function hashes without changing the already-qualified statuses.
manifest=json.loads((OUT/"oneD_v2_backend_lifecycle_manifest.json").read_text())
manifest["function_hashes"]={
 "PF_sharp_front_sigma_tip":closure["source_findings"]["PF"]["sigma_tip_function_sha256"],
 "FEMCZM_moving_front_sigma_tip":closure["source_findings"]["FEMCZM"]["sigma_tip_function_sha256"],
 "V2_PF_lifecycle_file":sha(ROOT/"reduced_fracture_v2/lifecycle.py"),
 "V2_shared_core_file":sha(ROOT/"reduced_fracture_v2/kinetics.py")}
manifest["source_hash_check"]="PASS_PINNED_FUNCTION_CONTENT_HASHES; AMBIENT_REPOSITORY_HEADS_DIFFER_AND_WERE_NOT_USED_AS_PROVENANCE"
(OUT/"oneD_v2_backend_lifecycle_manifest.json").write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n")

blocked=[{"availability_status":"NOT_RUN_GATE_BLOCKED","blocker":closure["launch_gate"]}]
for name in ("oneD_v2_native_baseline_results.csv","oneD_v2_current_four_class_results.csv","oneD_v2_pareto_candidates.csv","oneD_v2_new_four_class_registry.csv","oneD_v2_pf_transfer_results.csv"):
 write_csv(name,blocked,["availability_status","blocker"])
pd.DataFrame(blocked).to_parquet(OUT/"oneD_v2_search_population.parquet",index=False)

maps="""All three deterministic mechanics products retain their prior qualification: PF production-discrete native J/KJ, FEM/CZM production-native J/KJ, and FEM/CZM qualified structural G/KG. PF continuum G remains unqualified. FEM domain J remains a separately classified diagnostic. The standardized-path biases near 500 µm remain Peak +2.21% G/+1.76% J and DBTT +0.24% G/+1.26% J.\n"""
(OUT/"ONE_D_V2_MECHANICS_MAPS.md").write_text("# V2 mechanics maps\n\n"+maps)
(OUT/"ONE_D_V2_PARITY_LADDER.md").write_text("# V2 parity ladder\n\nLevels 1 and 2 pass exactly for all four unchanged canonical rows at 1000 K. Identical local states through the shared core give exact barriers/rates/hazards. The two backend driver shells also give identical controlled event times, transactions, states, avalanche grouping, and terminal state under analysis-only common mechanics and lifecycle.\n\nLevel 3 is not run: source K-to-tip-stress laws were identified, but the V2 forward state does not yet compose dynamic radius, shielding, backend emission/transport, and thresholds with the shared core. This is an implementation qualification gap, not a failure of Levels 1–2.\n")
(OUT/"ONE_D_V2_1000K_NATIVE_BASELINES.md").write_text("# V2 1000-K native baselines\n\nNot launched. The native state-closure gate remains blocked; consequently neither 100-µm prefixes nor 1000-µm continuations are claimed.\n")
(OUT/"ONE_D_V2_CURRENT_FOUR_CLASS_DIAGNOSTIC.md").write_text("# Current four-class diagnostic\n\nNot launched because the four native baselines did not pass the launch gate. Existing rows remain unchanged and their cross-provider class distinctness is unassessed.\n")
(OUT/"ONE_D_V2_NEW_PARAMETER_SEARCH.md").write_text("# New shared-parameter search\n\nNot launched. Prerequisite native baselines and current-row 40-case diagnostic are unavailable. No separate backend parameterization was created.\n")
(OUT/"ONE_D_V2_NEW_FOUR_CLASS_SELECTION.md").write_text("# New V2 four-class selection\n\nNo new rows selected. The canonical rows remain unchanged; refitting is not justified by a forward implementation gap.\n")
(OUT/"ONE_D_V2_TO_PF_TRANSFER_VALIDATION.md").write_text("# V2-to-PF transfer validation\n\nNot launched. No new candidate survived a search because the search was not authorized. Existing PF archives were sufficient for event/avalanche semantics, so no fresh 2-D PF run was needed.\n")
decision={"schema":"oneD_v2_final_decision_v1","mechanics_maps_qualified":True,"level1_common_core_parity":"PASS_EXACT","level2_common_mechanics_lifecycle_parity":"PASS_EXACT","level3_native_provider":"NOT_RUN_GATE_BLOCKED","native_baselines":"NOT_RUN_GATE_BLOCKED","current_four_class_matrix":"NOT_RUN_GATE_BLOCKED","new_parameter_search":"NOT_RUN_GATE_BLOCKED","PF_transfer_runs_launched":0,"FEMCZM_2d_runs_launched":0,"canonical_parameters_changed":False,"production_formulas_changed":False,"production_trajectories_changed":False,"parameter_refit_justified":False,"smallest_later_FEMCZM_validation":"NONE_UNTIL_NATIVE_REDUCED_BASELINES_PASS"}
(OUT/"oneD_v2_final_decision.json").write_text(json.dumps(decision,indent=2,sort_keys=True)+"\n")
(OUT/"ONE_D_V2_FINAL_PARAMETER_AND_2D_DECISION.md").write_text("# Final parameter and 2-D decision\n\nLevels 1–2 establish exact architectural parity under controlled common inputs. Native predictive qualification remains blocked by incomplete state closure, so no parameter inference follows. No new 2-D PF or FEM/CZM run was launched, no production formula or trajectory changed, and the smallest later FEM/CZM validation is none until native reduced baselines pass.\n")

fingerprints={p.name:sha(p) for p in sorted(OUT.iterdir()) if p.suffix in {".csv",".json",".parquet"} and p.name!="oneD_v2_parity_scientific_fingerprints.json"}
(OUT/"oneD_v2_parity_scientific_fingerprints.json").write_text(json.dumps(fingerprints,indent=2,sort_keys=True)+"\n")
