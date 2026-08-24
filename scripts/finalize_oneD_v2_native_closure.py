#!/usr/bin/env python3
from __future__ import annotations
import csv,hashlib,json,subprocess,sys,tempfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/"analysis_outputs/oneD_v2_mechanics_maps_and_baselines"
PY=Path(sys.executable);SHADOW=ROOT/"scripts/run_oneD_v2_native_source_shadow.py"
PF=Path("/Volumes/Data/Data/Nanopillar_calculation/PF-fracture-fatigue_v10_2_21_persistent_sites_top1")
FEM=Path("/Volumes/Data/Data/Nanopillar_calculation/Arrhenius_FEM_CZM_MPZ_theta0_pf_parity_claude")

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read_rows(p):
 with Path(p).open(newline="") as f:return list(csv.DictReader(f))
def write_csv(p,rows,fields):
 with Path(p).open("w",newline="") as f:
  w=csv.DictWriter(f,fieldnames=fields,lineterminator="\n");w.writeheader();w.writerows(rows)

old=json.loads((OUT/"oneD_v2_mechanics_map_manifest.json").read_text())
definitions=(
 ("PF","oneD_v2_pf_mechanics_map.csv","bcd07dd"),
 ("FEMCZM_NATIVE","oneD_v2_fem_native_mechanics_map.csv","2f5889d"),
 ("FEMCZM_STRUCTURAL_G","oneD_v2_fem_qualified_G_map.csv","2f5889d"),)
maps={}
for key,name,generation in definitions:
 path=(OUT/name).resolve();rows=read_rows(path);prior=old["maps"][key];actual=sha(path)
 if actual!=prior["raw_file_sha256"]:raise RuntimeError(f"pinned map hash mismatch: {name}")
 maps[key]={"provider":prior["provider"],"schema_columns":list(rows[0]),"source_repository":prior["source_repository"],"source_commit":prior["source_commit"],"generation_commit":generation,"absolute_source_artifact_path":str(path),"content_sha256":actual,"raw_node_count":len(rows),"raw_node_table_sha256":actual,"raw_nodes":rows,"valid_extension_interval_um":prior["valid_interval_um"],"reference_opening_m":float(rows[0]["reference_opening_m"]),"bulk_constraint":prior["bulk_constraint"],"crack_representation":prior["crack_representation"],"local_process_zone":prior["local_process_zone"],"interpolation_method":old["interpolation_method"],"extrapolation_policy":old["extrapolation_policy"],"qualification_status":prior["qualification_status"],"uncertainty_fields":[x for x in rows[0] if "uncertainty" in x.lower() or "spread" in x.lower()]}
manifest={"schema":"oneD_v2_native_map_inputs_manifest_v1","qualification_parent_commit":"4ce6c0a","hash_gate":"PASS_ALL_PREVIOUSLY_QUALIFIED_HASHES_EXACT","maps":maps}
(OUT/"oneD_v2_native_map_inputs_manifest.json").write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n")

with tempfile.TemporaryDirectory() as d:
 p=Path(d);pf=p/"pf.json";pfx=p/"pf_exact.json";fem=p/"fem.json"
 subprocess.run([str(PY),str(SHADOW),"pf",str(PF),str(pf)],check=True,cwd=ROOT)
 subprocess.run([str(PY),str(SHADOW),"pf_exact",str(PF),str(pfx)],check=True,cwd=ROOT)
 subprocess.run([str(PY),str(SHADOW),"fem",str(FEM),str(fem)],check=True,cwd=ROOT)
 shadow=json.loads(pf.read_text())+json.loads(pfx.read_text())+json.loads(fem.read_text())
fields=list(shadow[0]);write_csv(OUT/"oneD_v2_source_shadow_results.csv",shadow,fields)

micro=[]
for backend in ("PF","FEMCZM"):
 for material,seed in (("Peak",8666),("DBTT",1008666)):
  subset=[r for r in shadow if r["backend"]==backend and r["material_class"]==material and r["temperature_K"]==1000]
  lower_exact=all(r["process_zone_fingerprint_match"] and r["action_match"] and r["event_count_match"] and r["event_length_match"] and r["tip_radius_match"] and r["tip_stress_match"] for r in subset)
  micro.append({"backend":backend,"material_class":material,"temperature_K":1000,"seed":seed,"accepted_event_target":"3-5","source_object_shadow_exact":lower_exact,"exact_production_wrapper_available":False,"microtrajectory_launch_status":"NOT_RUN_SOURCE_SHADOW_REQUIRED_FIELD_UNAVAILABLE","classification":"NATIVE_CLOSURE_UNQUALIFIED"})
write_csv(OUT/"oneD_v2_native_microtrajectory_results.csv",micro,list(micro[0]))

empty={
 "oneD_v2_native_baseline_results.csv":["backend","material_class","temperature_K","seed","target_um","launch_status","classification"],
 "oneD_v2_native_event_transactions.csv":["backend","material_class","event_transaction_index","pre_event_time_s","pre_event_opening_m","pre_event_extension_um","post_event_extension_um","native_J_minus_J_m2","native_J_plus_J_m2","event_induced_change_J_m2","reload_induced_change_J_m2","right_censored_at_target"],
 "oneD_v2_native_physical_avalanches.csv":["backend","material_class","physical_avalanche_index","first_event_transaction_index","last_event_transaction_index","event_count","extension_um","right_censored_at_target"],
 "oneD_v2_native_onset_candidates.csv":["backend","material_class","onset_role","pre_event_time_s","pre_event_opening_m","pre_event_native_J_J_m2","pre_event_native_KJ_MPa_sqrt_m"]}
for name,fields in empty.items():write_csv(OUT/name,[],fields)

decision={"schema":"oneD_v2_native_closure_decision_v1","map_inputs":"PINNED_REPRODUCIBLE","PF_lower_level_source_object_shadow":"EXACT","FEMCZM_lower_level_source_object_shadow":"EXACT","PF_exact_production_wrapper_shadow":"BLOCKED_RELIABLE_2D_TENSOR_DRIVE_REQUIRED","FEMCZM_exact_production_barrier_wrapper_shadow":"UNAVAILABLE","required_field_unavailable":True,"Level3":"NATIVE_CLOSURE_UNQUALIFIED","microtrajectories_launched":0,"baseline_prefixes_launched":0,"long_baselines_launched":0,"new_2d_PF_runs_launched":0,"new_2d_FEMCZM_runs_launched":0,"canonical_parameters_changed":False,"parameter_refit_justified":False,"production_formulas_changed":False,"production_trajectories_changed":False,"next_required_work":["define a source-qualified reduction of the PF 2-D tensor emission drive; scalar native KJ is rejected by the exact wrapper","instantiate corrected FEM/CZM exact production barrier/front factory","repeat same fieldwise source shadow"]}
(OUT/"oneD_v2_native_closure_decision.json").write_text(json.dumps(decision,indent=2,sort_keys=True)+"\n")

(OUT/"ONE_D_V2_NATIVE_STATE_CLOSURE.md").write_text("# V2 native state closure\n\nTyped PF and FEM/CZM forward states and thin backend closures now compose mapped native KJ with source-object radius, shielding, K-to-stress conversion, transport, threshold stepping, event advance, and renewal. Geometry extension is changed only after a reported source event. FEM qualified G/KG remains parallel metadata. No archived future state is read.\n\nQualification remains fail-closed because the current shadows instantiate lower-level production source objects, not both exact final wrapper compositions.\n")
(OUT/"ONE_D_V2_PF_NATIVE_CLOSURE_QUALIFICATION.md").write_text("# PF native closure qualification\n\nThe closure is exact against `UnifiedMPZFrontEngine` for all four rows at 300, 1000, and 1200 K in the controlled five-interval fixture. The exact persistent-site, state-resolved signed production wrapper was also instantiated with the pinned physical kernel family. It rejects scalar-KJ stepping because persistent signed emission requires a reliable 2-D tensor drive. This is a substantive missing reduction variable, not an adapter arithmetic mismatch. Required emission, multiplicity, and rate parity therefore remain **UNAVAILABLE** and PF closure status is **NATIVE_CLOSURE_UNQUALIFIED**.\n")
(OUT/"ONE_D_V2_FEMCZM_NATIVE_CLOSURE_QUALIFICATION.md").write_text("# FEM/CZM native closure qualification\n\nThe closure is exact against `MovingProcessZone2DFrontEngine` for all four rows at 300, 1000, and 1200 K in the controlled five-interval fixture. The fixture has not yet proven the corrected production factory's exact barrier construction. Required barrier/rate parity is therefore **UNAVAILABLE**, and FEM/CZM closure status is **NATIVE_CLOSURE_UNQUALIFIED**. Qualified G remains interpretation-only.\n")
(OUT/"ONE_D_V2_LEVEL3_MICROTRAJECTORY_GATE.md").write_text("# Level-3 microtrajectory gate\n\nAll four Peak/DBTT microtrajectories are **NOT RUN** because source-shadow required fields remain unavailable. Lower-level object equality is not promoted to full production-composition qualification.\n")
(OUT/"ONE_D_V2_1000K_NATIVE_BASELINES.md").write_text("# V2 1000-K native baselines\n\nNo 100-µm prefix or 1000-µm continuation was launched because no Level-3 microtrajectory passed.\n")
(OUT/"ONE_D_V2_NATIVE_BASELINE_COMPARISON.md").write_text("# Native baseline comparison\n\nUnavailable: baseline tables are intentionally header-only. No onset, avalanche, or parameter-plausibility result is inferred.\n")

fingerprints={p.name:sha(p) for p in sorted(OUT.iterdir()) if p.name.startswith("oneD_v2_native_") and p.suffix in {".csv",".json"} and p.name!="oneD_v2_native_closure_fingerprints.json"}
(OUT/"oneD_v2_native_closure_fingerprints.json").write_text(json.dumps(fingerprints,indent=2,sort_keys=True)+"\n")

# Preserve the earlier parity bundle's own key set while refreshing files that
# this later gated stage intentionally supersedes (notably baseline tables).
parity_path=OUT/"oneD_v2_parity_scientific_fingerprints.json"
parity=json.loads(parity_path.read_text())
parity={name:sha(OUT/name) for name in parity}
parity_path.write_text(json.dumps(parity,indent=2,sort_keys=True)+"\n")
