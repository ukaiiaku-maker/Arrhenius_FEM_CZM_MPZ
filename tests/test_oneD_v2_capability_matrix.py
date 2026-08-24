from pathlib import Path
import json,pandas as pd,pytest
from reduced_fracture_v2.lifecycle import PFLateGeometryVetoPolicy,pf_fail_closed_veto
OUT=Path(__file__).parents[1]/"analysis_outputs/oneD_v2_capability_v2"
CAP=pd.read_csv(OUT/"oneD_v2_backend_capabilities.csv")
def status(b,c):return CAP[(CAP.backend==b)&(CAP.capability==c)].status.iloc[0]
def test_01_unsupported_rollback_distinct():assert status("PF","LATE_VETO_ROLLBACK_AND_CONTINUE")=="UNSUPPORTED_BY_PRODUCTION" and status("PF","NORMAL_ACCEPTED_EVENT")=="QUALIFIED"
def test_02_fail_closed_explicit():assert PFLateGeometryVetoPolicy.FAIL_CLOSED_TERMINATE.value=="FAIL_CLOSED_TERMINATE"
def test_03_tentative_not_authoritative():assert not json.loads((OUT/"pf_late_geometry_veto_fail_closed_audit.json").read_text())["tentative_geometry_published"]
def test_04_pf_commits_once():assert all(x==[1,1,1] for x in [z["geometry_commit_counts"] for z in json.loads((OUT/"oneD_v2_diagnostic_neutrality.json").read_text())["cases"]])
def test_05_subdivision_honest():assert status("PF","SUBDIVISION_THRESHOLD_PRESERVATION")=="PARTIALLY_QUALIFIED"
def test_06_event_length_once():assert json.loads((OUT/"oneD_v2_rng_threshold_audit.json").read_text())["PF"]["event_length_draws_per_event"]==1
def test_07_renewal_once():assert status("PF","POST_EVENT_PROCESS_ZONE_RENEWAL")=="QUALIFIED"
def test_08_neutral():assert json.loads((OUT/"oneD_v2_diagnostic_neutrality.json").read_text())["pass"]
def test_09_fem_independent():assert status("FEMCZM","NORMAL_ACCEPTED_EVENT")=="QUALIFIED"
def test_10_maps_independent_pending():assert status("PF","MECHANICS_MAP")=="NOT_EVALUATED"
def test_11_campaign_blocked():assert status("PF","FORWARD_REDUCED_MODEL")=="BLOCKED_BY_REQUIRED_CAPABILITIES"
def test_12_pf_not_G():assert json.loads((OUT/"oneD_v2_mechanics_map_manifest.json").read_text())["PF"]["continuum_G_status"]=="UNQUALIFIED"
def test_13_fem_maps_distinct():
 m=json.loads((OUT/"oneD_v2_mechanics_map_manifest.json").read_text())["FEMCZM"];assert "native_J_map_status" in m and "qualified_G_map_status" in m
def test_14_metadata_separate():
 m=json.loads((OUT/"oneD_v2_mechanics_map_manifest.json").read_text());assert m["PF"]["bulk_constraint"]!="FINITE_WIDTH_STIFFNESS_KILLED_WAKE"
def test_15_pz_not_interface():assert json.loads((OUT/"oneD_v2_mechanics_map_manifest.json").read_text())["FEMCZM"]["local_process_zone"]=="FINITE_SIZE_REDUCED_STATE"
def test_16_bounded_interpolation_policy():assert json.loads((OUT/"oneD_v2_mechanics_map_manifest.json").read_text())["interpolation"].startswith("BOUNDED_PIECEWISE_LINEAR")
def test_17_no_extrapolation():assert "NO_EXTRAPOLATION" in json.loads((OUT/"oneD_v2_mechanics_map_manifest.json").read_text())["interpolation"]
def test_18_no_parameter_or_campaign_change():assert not json.loads((OUT/"oneD_v2_mechanics_map_manifest.json").read_text())["campaign_outputs_populated"]
