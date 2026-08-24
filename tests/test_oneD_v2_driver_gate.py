from pathlib import Path
import json,pandas as pd
OUT=Path(__file__).parents[1]/"analysis_outputs/oneD_v2_driver_diagnostics"
def j(n):return json.loads((OUT/n).read_text())
def test_01_neutrality_not_claimed():assert not j("oneD_v2_diagnostic_neutrality.json")["PF"]["diagnostic_pair_run"]
def test_02_threshold_subdivision_unqualified():assert j("oneD_v2_rng_threshold_audit.json")["PF"]["threshold_subdivision"]=="NOT_QUALIFIED"
def test_03_rng_rollback_unqualified():assert j("oneD_v2_rng_threshold_audit.json")["PF"]["rng_rollback"]=="NOT_QUALIFIED"
def test_04_rejected_trial_time_not_inferred():assert "NOT_QUALIFIED" in j("oneD_v2_rng_threshold_audit.json")["PF"].values()
def test_05_exact_commit_not_claimed():assert j("oneD_v2_rollback_audit.json")["PF"]["classification"]=="UNQUALIFIED"
def test_06_event_length_phase_not_inferred():assert "event-length" not in (OUT/"oneD_v2_driver_lifecycle_records.jsonl").read_text()
def test_07_renewal_not_claimed():assert "LIFECYCLE_QUALIFIED" not in (OUT/"ONE_D_V2_PF_DRIVER_LIFECYCLE_QUALIFICATION.md").read_text()
def test_08_adapter_owned_fields_still_audited():assert (OUT.parent/"oneD_v2_backend_qualification/oneD_v2_source_state_sufficiency_matrix.csv").exists()
def test_09_driver_state_not_inferred():assert j("oneD_v2_rng_threshold_audit.json")["FEMCZM"]["status"].startswith("NOT_RUN")
def test_10_metadata_separates_constraint_and_crack():
 x=j("oneD_v2_mechanics_map_manifest.json");assert x["bulk_constraint"]=="PLANE_STRAIN" and "crack_representation" in x["FEMCZM"]
def test_11_process_zone_not_interface_thickness():assert j("oneD_v2_mechanics_map_manifest.json")["FEMCZM"]["local_process_zone"]=="FINITE_SIZE_REDUCED_STATE"
def test_12_pf_map_not_candidate_tuned():assert j("oneD_v2_mechanics_map_manifest.json")["PF"]["map_status"]=="NOT_GENERATED"
def test_13_fem_map_not_candidate_tuned():assert j("oneD_v2_mechanics_map_manifest.json")["FEMCZM"]["map_status"]=="NOT_GENERATED"
def test_14_pf_no_G_claim():assert "continuum-G claim is made" in (OUT/"ONE_D_V2_PF_MECHANICS_MAP_QUALIFICATION_V2.md").read_text()
def test_15_fem_native_and_G_files_separate():assert (OUT/"oneD_v2_fem_native_mechanics_map.csv").exists() and (OUT/"oneD_v2_fem_qualified_G_map.csv").exists()
def test_16_scaling_not_fabricated():assert pd.read_csv(OUT/"oneD_v2_mechanics_map_interpolation_tests.csv").availability_status.iloc[0].startswith("NOT_GENERATED")
def test_17_interpolation_not_fabricated():assert len(pd.read_csv(OUT/"oneD_v2_mechanics_map_interpolation_tests.csv"))==1
def test_18_extrapolation_not_enabled():assert "extrapolat" not in (OUT/"oneD_v2_mechanics_map_manifest.json").read_text().lower()
def test_19_campaign_rows_not_produced():assert not j("oneD_v2_mechanics_map_manifest.json")["campaign_tables_populated"]
def test_20_no_long_2d_run_launched():assert not j("oneD_v2_diagnostic_neutrality.json")["physical_trajectory_changed"]
