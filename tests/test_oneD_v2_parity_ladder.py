import hashlib,json
from pathlib import Path
import pandas as pd
from reduced_fracture_v2 import SharedEventAvalancheAnalysis
from reduced_fracture_v2.lifecycle import PFLifecyclePolicy,FEMCZMLifecyclePolicy

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/"analysis_outputs/oneD_v2_mechanics_maps_and_baselines"

def test_required_architecture_names_are_distinct():
 assert PFLifecyclePolicy is not FEMCZMLifecyclePolicy
 assert SharedEventAvalancheAnalysis.__name__=="SharedEventAvalancheAnalysis"
def test_common_parity_levels_are_exact():
 x=pd.read_csv(OUT/"oneD_v2_common_parity_results.csv")
 cols=[c for c in x if c.startswith("level") and c.endswith("exact")]
 assert len(x)==4 and x[cols].all().all() and (x.level1_max_absolute_difference==0).all()
def test_avalanche_grouping_uses_reload_not_time_threshold():
 rows=[{"event_transaction_index":0,"event_extension_m":1,"reload_before_event":False},{"event_transaction_index":1,"event_extension_m":1,"reload_before_event":False},{"event_transaction_index":2,"event_extension_m":1,"reload_before_event":True,"right_censored_at_target":True}]
 g=SharedEventAvalancheAnalysis.group(rows);assert [x["event_count"] for x in g]==[2,1] and g[-1]["right_censored_at_target"]
def test_native_and_downstream_gates_fail_closed():
 d=json.loads((OUT/"oneD_v2_final_decision.json").read_text())
 assert d["level3_native_provider"].startswith("NOT_RUN") and d["PF_transfer_runs_launched"]==0 and d["FEMCZM_2d_runs_launched"]==0
 assert not d["canonical_parameters_changed"] and not d["parameter_refit_justified"]
def test_scientific_fingerprints_reproduce():
 d=json.loads((OUT/"oneD_v2_parity_scientific_fingerprints.json").read_text())
 assert all(hashlib.sha256((OUT/n).read_bytes()).hexdigest()==h for n,h in d.items())
