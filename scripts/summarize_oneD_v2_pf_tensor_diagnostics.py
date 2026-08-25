#!/usr/bin/env python3
from __future__ import annotations
import hashlib,json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'analysis_outputs/oneD_v2_local_drive';RUN=Path('/private/tmp/pf-v2-tensor-diagnostics')
PHYSICAL=['anisotropic_emission_audit_v10174.json','branch_diagnostics_1000K.csv','crack_path_1000K.csv','crack_path_front0_1000K.csv','fronts_1000K.csv','kinetic_tip_cell_audit_v101.json','selected_material_manifest_v10_2_22.csv','sharp_wake_advance_log.csv','steps_1000K.csv','stochastic_avalanche_geometry_events.json','summary.json','v10_2_22_persistent_site_model.json']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 result={'schema':'oneD_v2_pf_tensor_diagnostic_neutrality_v1','observer_default_off':True,'cases':{}}
 for m,label in [('peak','Peak'),('dbtt','DBTT')]:
  off,on=RUN/f'{m}_off',RUN/f'{m}_on';pairs={n:{'off_sha256':sha(off/n),'on_sha256':sha(on/n),'identical':sha(off/n)==sha(on/n)} for n in PHYSICAL}
  summary=json.loads((on/'summary.json').read_text())[0];tensor=on/'tensor_drive.jsonl'
  first=json.loads(tensor.open().readline())
  result['cases'][label]={'seed':8666 if m=='peak' else 1008666,'temperature_K':1000,'target_extension_um':10,
   'accepted_geometry_events':summary['n_geometry_events'],'final_extension_um':summary['geometry_projected_extension_m']*1e6,
   'tensor_evaluation_count':sum(1 for _ in tensor.open()),'tensor_jsonl_sha256':sha(tensor),'all_physical_artifacts_identical':all(v['identical'] for v in pairs.values()),
   'physical_artifact_hashes':pairs,'first_drive':{k:first[k] for k in ('opening_tensor_Pa','channel_tensors_Pa','tau_signed_Pa','drive_factors','front_direction','front_normal','tip_xy_m','reliable')}}
 result['qualification']='PASS_OFF_ON_PHYSICAL_BYTE_IDENTITY' if all(c['all_physical_artifacts_identical'] for c in result['cases'].values()) else 'FAIL_OBSERVER_NOT_NEUTRAL'
 (OUT/'oneD_v2_pf_tensor_diagnostic_neutrality.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
 fingerprints={p.name:sha(p) for p in sorted(OUT.iterdir()) if p.is_file() and p.name!='oneD_v2_local_drive_fingerprints.json'}
 (OUT/'oneD_v2_local_drive_fingerprints.json').write_text(json.dumps(fingerprints,indent=2,sort_keys=True)+'\n')
if __name__=='__main__':main()
