#!/usr/bin/env python3
"""Close the V2 forward gate without inventing a native-drive/local-state law."""
from __future__ import annotations

import csv, hashlib, json, sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import matplotlib.pyplot as plt
from reduced_fracture_v2.mechanics import FEMCZMMechanicsProvider, PFMechanicsProvider

OUT=ROOT/"analysis_outputs/oneD_v2_mechanics_maps_and_baselines"

def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def empty_csv(name, fields):
    with (OUT/name).open("w",newline="") as f: csv.DictWriter(f,fieldnames=fields).writeheader()

pf_path=OUT/"oneD_v2_pf_mechanics_map.csv"
fn_path=OUT/"oneD_v2_fem_native_mechanics_map.csv"
fg_path=OUT/"oneD_v2_fem_qualified_G_map.csv"
pf=PFMechanicsProvider.from_csv(pf_path,sha(pf_path))
fem=FEMCZMMechanicsProvider.from_csv(fn_path,fg_path,sha(fn_path)+":"+sha(fg_path))

# Exercise the typed data contract at an interior state.  These are interface
# checks only, not stochastic trajectory calculations.
ps=pf.evaluate(100e-6,10e-6); fs=fem.evaluate(100e-6,10e-6)
assert ps.pf_native_J_J_m2 is not None and ps.pf_native_KJ_Pa_sqrt_m is not None
assert ps.structural_G_J_m2 is None and ps.local_tip_stress_Pa is None
assert fs.fem_native_J_J_m2 is not None and fs.fem_native_KJ_Pa_sqrt_m is not None
assert fs.structural_G_J_m2 is not None and fs.qualified_KG_Pa_sqrt_m is not None
assert fs.local_tip_stress_Pa is None

event_fields=["provider","material_class","seed","event_transaction_index","pre_event_time_s","pre_event_opening_m","pre_event_reaction_N_per_m","pre_event_native_J_J_per_m2","pre_event_native_KJ_MPa_sqrt_m","pre_event_projected_extension_um","post_event_projected_extension_um","event_extension_um","post_event_process_zone_fingerprint","reload_time_to_next_event_s","reload_opening_to_next_event_m","right_censored_at_target"]
av_fields=["provider","material_class","seed","physical_avalanche_index","first_event_transaction_index","last_event_transaction_index","event_count","extension_um","largest_avalanche_fraction","right_censored_at_target"]
onset_fields=["provider","material_class","seed","onset_role","pre_event_time_s","pre_event_opening_m","pre_event_native_J_J_per_m2","pre_event_native_KJ_MPa_sqrt_m"]
summary_fields=["provider","material_class","temperature_K","seed","prefix_target_um","long_target_um","launch_status","classification","blocking_gate","parameter_status"]
empty_csv("oneD_v2_baseline_event_transactions.csv",event_fields)
empty_csv("oneD_v2_baseline_physical_avalanches.csv",av_fields)
empty_csv("oneD_v2_baseline_onset_candidates.csv",onset_fields)
empty_csv("oneD_v2_baseline_summary.csv",summary_fields)

cases=[{"provider":p,"material_class":m,"seed":s,"status":"NOT_LAUNCHED"}
       for p in ("PF_CONSISTENT","FEMCZM_CONSISTENT")
       for m,s in (("Peak",8666),("DBTT",1008666))]
decision={
 "schema":"oneD_v2_baseline_decision_v1",
 "mechanics_maps":{"PF":"PF_PRODUCTION_DISCRETE_MAP_QUALIFIED","FEM_native":"QUALIFIED","FEM_structural_G":"QUALIFIED"},
 "typed_provider_interface":"QUALIFIED",
 "forward_launch_gate":"BLOCKED_NATIVE_METRIC_TO_LOCAL_PROCESS_ZONE_CLOSURE_ABSENT",
 "blocking_evidence":{
   "hazard_kernel_requires":["cleavage_stress_Pa","emission_stress_Pa"],
   "PF_map_provides":["reaction","PF_native_J","PF_native_KJ"],
   "FEM_map_provides":["reaction","FEM_native_J","FEM_native_KJ","qualified_G","qualified_KG","G_uncertainty"],
   "missing_source_qualified_contract":"backend-native mechanics plus process-zone state to local cleavage/emission stresses",
   "prohibited_substitution":"native KJ or qualified KG cannot be passed as a local stress"
 },
 "cases":cases,"prefixes_launched":0,"long_cases_launched":0,"campaign_launched":False,
 "canonical_parameters_changed":False,"parameter_status":"INSUFFICIENT_EVIDENCE",
 "refit_justified":False,"new_2d_stochastic_run_required":False,
 "next_required_work":"derive and source-test backend-specific local stress/state closures; deterministic frozen-state diagnostics are sufficient",
 "physical_or_stochastic_law_altered":False,"production_trajectory_altered":False,
 "map_hashes":{"PF":sha(pf_path),"FEM_native":sha(fn_path),"FEM_structural_G":sha(fg_path)}
}
(OUT/"oneD_v2_baseline_decision.json").write_text(json.dumps(decision,indent=2,sort_keys=True)+"\n")

common="""The PF production-discrete, FEM/CZM native, and FEM/CZM structural-G maps pass their deterministic mechanics qualifications. The typed provider interface also passes and keeps native J/KJ separate from qualified G/KG.\n\nThe four 1000-K reduced baselines were **not launched**. The forward gate fails closed because the shared hazard kernel consumes local cleavage and emission stresses, while neither qualified mechanics map nor a source-qualified adapter defines how backend-native global mechanics and the evolving process-zone state produce those local stresses. Substituting KJ, KG, J, or G would change the model and its units. Therefore no 100-µm prefix or 1000-µm continuation result is claimed.\n\nAll canonical parameter rows remain unchanged. Their plausibility cannot be judged from an unexecuted baseline, and no refit is justified. No new 2-D stochastic run is required; source-level closure work and bounded deterministic frozen-state tests are sufficient.\n"""
(OUT/"ONE_D_V2_1000K_FORWARD_BASELINE_RESULTS.md").write_text("# One-dimensional V2 1000-K forward baseline results\n\n"+common)
(OUT/"ONE_D_V2_PF_1000K_BASELINE_COMPARISON.md").write_text("# PF 1000-K reduced-versus-2-D comparison\n\nNo reduced PF baseline was launched, so onset and avalanche agreement are **INSUFFICIENT_EVIDENCE**. The authoritative 2-D reference remains one 208-event target-censored avalanche for Peak and two reload-separated avalanches for DBTT (11 and 197 events).\n\n"+common)
(OUT/"ONE_D_V2_FEMCZM_1000K_BASELINE_COMPARISON.md").write_text("# FEM/CZM 1000-K reduced-versus-2-D comparison\n\nNo reduced FEM/CZM baseline was launched, so comparison with the corrected 2-D Peak and DBTT topologies is **INSUFFICIENT_EVIDENCE**.\n\n"+common)
(OUT/"ONE_D_V2_BASELINE_PARAMETER_DECISION.md").write_text("# V2 baseline parameter decision\n\n**INSUFFICIENT_EVIDENCE.** Parameters are unchanged; no refit is justified. The remaining discrepancy is an unimplemented, source-unqualified mechanics-to-local-state closure, not evidence of parameter error.\n")

for name,title in (("ONE_D_V2_PF_BASELINES_VS_2D.png","PF baseline gate"),("ONE_D_V2_FEMCZM_BASELINES_VS_2D.png","FEM/CZM baseline gate"),("ONE_D_V2_BASELINE_AVALANCHE_TOPOLOGY.png","Reduced baseline topology gate")):
    fig,ax=plt.subplots(figsize=(8,3.4)); ax.axis("off")
    ax.text(.5,.62,"NOT LAUNCHED",ha="center",va="center",fontsize=20,weight="bold")
    ax.text(.5,.37,"Native mechanics → local cleavage/emission stress closure absent",ha="center",va="center",wrap=True)
    ax.set_title(title); fig.tight_layout(); fig.savefig(OUT/name,dpi=180); plt.close(fig)
