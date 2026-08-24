#!/usr/bin/env python3
from pathlib import Path
import csv,json,hashlib

OUT=Path("analysis_outputs/oneD_v2_driver_diagnostics");OUT.mkdir(parents=True,exist_ok=True)
PF_COMMIT="ab6279331d050919f78e2e7f278bc332466e34f8"
FEM_COMMIT="931bed66913afc970117bce900805ecc9b6225f8"
records=[
 {"backend":"PF","material_class":"Peak","temperature_K":1000,"seed":8666,"phase":"ROLLBACK_BEGIN","source_commit":PF_COMMIT,"status":"IRREDUCIBLE_BLOCKER_BEFORE_RUN","detail":"production restore_geometry_veto raises; exact coupled-step replay unavailable"},
 {"backend":"PF","material_class":"DBTT","temperature_K":1000,"seed":1008666,"phase":"ROLLBACK_BEGIN","source_commit":PF_COMMIT,"status":"IRREDUCIBLE_BLOCKER_BEFORE_RUN","detail":"same authoritative engine and rollback contract"},
 {"backend":"FEMCZM","material_class":"Peak","temperature_K":1000,"seed":8666,"phase":"INITIAL_AUTHORITATIVE_STATE","source_commit":FEM_COMMIT,"status":"NOT_LAUNCHED_PF_COMMON_GATE_ALREADY_FAILED","detail":"no long or unnecessary diagnostic launched"},
 {"backend":"FEMCZM","material_class":"DBTT","temperature_K":1000,"seed":1008666,"phase":"INITIAL_AUTHORITATIVE_STATE","source_commit":FEM_COMMIT,"status":"NOT_LAUNCHED_PF_COMMON_GATE_ALREADY_FAILED","detail":"no long or unnecessary diagnostic launched"}]
with (OUT/"oneD_v2_driver_lifecycle_records.jsonl").open("w") as f:
 for r in records:f.write(json.dumps(r,sort_keys=True)+"\n")
with (OUT/"oneD_v2_driver_state_boundaries.csv").open("w",newline="") as f:
 w=csv.DictWriter(f,fieldnames=records[0].keys());w.writeheader();w.writerows(records)
audit={"PF":{"threshold_subdivision":"NOT_QUALIFIED","rng_rollback":"NOT_QUALIFIED","reason":"late-veto full-state restore absent in authoritative driver"},"FEMCZM":{"status":"NOT_RUN_AFTER_PF_GATE_FAILURE"}}
(OUT/"oneD_v2_rng_threshold_audit.json").write_text(json.dumps(audit,indent=2)+"\n")
rollback={"PF":{"source_commit":PF_COMMIT,"classification":"UNQUALIFIED","direct_source_evidence":{"audit_payload":"fail_closed_no_partial_state_rollback","restore_geometry_veto":"raises RuntimeError","driver_callsite":"sharp_front.py late geometry veto"}},"FEMCZM":{"source_commit":FEM_COMMIT,"classification":"NOT_EVALUATED"}}
(OUT/"oneD_v2_rollback_audit.json").write_text(json.dumps(rollback,indent=2)+"\n")
neutral={"PF":{"diagnostic_pair_run":False,"reason":"instrumentation cannot repair missing production rollback"},"FEMCZM":{"diagnostic_pair_run":False,"reason":"stopped at failed common forward gate"},"physical_trajectory_changed":False}
(OUT/"oneD_v2_diagnostic_neutrality.json").write_text(json.dumps(neutral,indent=2)+"\n")
docs={
"ONE_D_V2_BOUNDED_DRIVER_DIAGNOSTICS.md":"The authorized diagnostics were provenance-resolved to PF commit `ab627933` and FEM/CZM commit `931bed6`. Before launch, direct production-source inspection found an irreducible PF rollback failure: late geometry veto invokes `restore_geometry_veto`, which intentionally raises because coupled-step replay is unavailable. A forced-veto diagnostic therefore cannot satisfy full-state restoration. No run was launched or extended, and no trajectory result is claimed.",
"ONE_D_V2_PF_DRIVER_LIFECYCLE_QUALIFICATION.md":"Classification: **UNQUALIFIED**. The production audit payload explicitly records `fail_closed_no_partial_state_rollback`. The mandatory forced late-veto path cannot restore mechanics, MPZ, hazard, RNG, and geometry state. Instrumentation alone cannot change this behavior without becoming a production-physics transaction fix.",
"ONE_D_V2_FEMCZM_DRIVER_LIFECYCLE_QUALIFICATION.md":"Classification: **NOT EVALUATED**. FEM/CZM source commit `931bed6` was resolved, but bounded runs were not launched after the shared forward gate had already failed on PF. Its independent lifecycle remains unqualified; object-level snapshot tests do not substitute for driver qualification."}
docs["ONE_D_V2_PF_MECHANICS_MAP_QUALIFICATION_V2.md"]="Classification: **NOT GENERATED**. Deterministic map work is scientifically independent, but no exact production-operator map was completed in this milestone. No continuum-G claim is made."
docs["ONE_D_V2_FEMCZM_MECHANICS_MAP_QUALIFICATION_V2.md"]="Classification: **NOT GENERATED**. No native-J or qualified energy/compliance/VCCT map is claimed. Plane strain, crack representation, failed-interface state, and local process-zone representation are stored as independent metadata."
for n,t in docs.items():(OUT/n).write_text("# "+n[:-3].replace("_"," ").title()+"\n\n"+t+"\n")
manifest={"bulk_constraint":"PLANE_STRAIN","PF":{"crack_representation":"FINITE_WIDTH_STIFFNESS_KILLED_WAKE","local_process_zone":"FINITE_SIZE_REDUCED_STATE","map_status":"NOT_GENERATED"},"FEMCZM":{"crack_representation":"CODIMENSION_ONE_DISPLACEMENT_DISCONTINUITY","failed_interface_state":"TRACTION_FREE","local_process_zone":"FINITE_SIZE_REDUCED_STATE","map_status":"NOT_GENERATED"},"campaign_tables_populated":False}
(OUT/"oneD_v2_mechanics_map_manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
for n in ("oneD_v2_pf_mechanics_map.csv","oneD_v2_fem_native_mechanics_map.csv","oneD_v2_fem_qualified_G_map.csv","oneD_v2_mechanics_map_interpolation_tests.csv"):
 (OUT/n).write_text("availability_status\nNOT_GENERATED_LIFECYCLE_GATE_FAILED\n")
