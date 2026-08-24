import hashlib,json
from pathlib import Path
import pandas as pd

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/"analysis_outputs/oneD_v2_mechanics_maps_and_baselines"

def test_actual_map_hashes_match_pinned_artifacts():
 d=json.loads((OUT/"oneD_v2_native_map_inputs_manifest.json").read_text());assert d["hash_gate"].startswith("PASS_")
 for m in d["maps"].values():assert hashlib.sha256(Path(m["absolute_source_artifact_path"]).read_bytes()).hexdigest()==m["content_sha256"] and len(m["raw_nodes"])==m["raw_node_count"]
def test_native_observable_separation_is_pinned():
 d=json.loads((OUT/"oneD_v2_native_map_inputs_manifest.json").read_text())
 assert "G" not in d["maps"]["PF"]["uncertainty_fields"]
 assert d["maps"]["FEMCZM_NATIVE"]["content_sha256"]!=d["maps"]["FEMCZM_STRUCTURAL_G"]["content_sha256"]
def test_lower_level_source_object_shadow_is_exact_but_not_overclaimed():
 x=pd.read_csv(OUT/"oneD_v2_source_shadow_results.csv");lower=x[x.backend.isin(["PF","FEM"])]
 exact=["process_zone_fingerprint_match","action_match","event_count_match","event_length_match","tip_radius_match","tip_stress_match","fired_match"]
 assert len(lower)==144 and lower[exact].all().all() and not lower.production_wrapper_exact.any()
 exact_pf=x[x.backend=="PF_EXACT"];assert len(exact_pf)==12 and exact_pf.production_wrapper_exact.all()
 assert exact_pf.closure_status.str.contains("reliable 2-D tensor drive",case=False).all()
 assert (x.required_field_status=="UNAVAILABLE").all()
def test_geometry_commit_and_reload_order_are_explicit():
 x=pd.read_csv(OUT/"oneD_v2_source_shadow_results.csv");x=x[(x.backend=="PF")&(x.material_class=="Peak")&(x.temperature_K==1000)]
 assert (x.loc[x.geometry_commits==0,"pre_extension_m"]==x.loc[x.geometry_commits==0,"post_extension_m"]).all()
 assert (x.geometry_commits==x.post_event_renewals).all() and x.geometry_commits.sum()==3
 reload=x[x.interval_index==3].iloc[0];assert reload.geometry_commits==0 and reload.pre_extension_m==reload.post_extension_m
def test_level3_and_baselines_fail_closed():
 d=json.loads((OUT/"oneD_v2_native_closure_decision.json").read_text())
 assert d["Level3"]=="NATIVE_CLOSURE_UNQUALIFIED" and d["microtrajectories_launched"]==d["baseline_prefixes_launched"]==0
 assert not d["canonical_parameters_changed"] and not d["parameter_refit_justified"]
def test_no_2d_runs_or_production_changes():
 d=json.loads((OUT/"oneD_v2_native_closure_decision.json").read_text())
 assert d["new_2d_PF_runs_launched"]==d["new_2d_FEMCZM_runs_launched"]==0
 assert not d["production_formulas_changed"] and not d["production_trajectories_changed"]
def test_native_closure_fingerprints_reproduce():
 d=json.loads((OUT/"oneD_v2_native_closure_fingerprints.json").read_text())
 assert all(hashlib.sha256((OUT/n).read_bytes()).hexdigest()==h for n,h in d.items())
