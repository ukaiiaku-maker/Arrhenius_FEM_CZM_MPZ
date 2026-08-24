#!/usr/bin/env python3
"""Exercise actual backend process-zone classes in isolated source imports."""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parents[1];sys.path.insert(0,str(HERE))
from reduced_fracture_v2.adapters import FEMCZMProcessZoneStateAdapter,PFProcessZoneStateAdapter

def pf(repo):
 sys.path.insert(0,str(repo))
 from arrhenius_fracture.config import ElasticProperties
 from arrhenius_fracture.material_manifest import MaterialManifest,default_manifest_path
 from arrhenius_fracture.sharp_front import FrontConfig,default_cleavage_barrier,default_emission_barrier
 from arrhenius_fracture.unified_front import UnifiedMPZFrontEngine
 from arrhenius_fracture.unified_mpz import MPZConfig
 mat=ElasticProperties();f=FrontConfig();f.r0=1e-6;f.sigma_cap=0;f.da=5e-6;f.L_pz=100e-6
 e=UnifiedMPZFrontEngine(f,default_cleavage_barrier(),default_emission_barrier(mat.b),mat.G,mat.nu,mat.b,MaterialManifest.from_csv(default_manifest_path("DBTT")),MPZConfig(length_m=100e-6,n_bins=40,wake_length_m=100e-6))
 ad=PFProcessZoneStateAdapter();s0=ad.snapshot(e.mpz);rates={"cleavage":e.lambda_cleave(2e9,1000)[0],"emission":e.lambda_emit(2e9,1000)[0]}
 no=e.step(0,1000,0);e.mpz.mobile[:,:2]=1;before=e.N_em;em=e.mpz.evolve(1e-9,1000,3e9,mat.b);emission_changed=e.N_em!=before
 ad.restore(e.mpz,s0);restored=ad.snapshot(e.mpz).fingerprint==s0.fingerprint
 e.B=1.25;cross=e.step(0,1000,0);renewal=(e.n_adv==1 and np.isclose(e.B,.25))
 snap=ad.snapshot(e.mpz);e.mpz.mobile[0,0]+=9;ad.restore(e.mpz,snap);rollback=ad.snapshot(e.mpz).fingerprint==snap.fingerprint
 return {"backend":"PF","production_class":"UnifiedMPZFrontEngine/UnifiedMPZState","frozen_rates":rates,"no_event_update_fired":bool(no["fired"]),"emission_only_state_changed":bool(emission_changed),"cleavage_crossing_fired":bool(cross["fired"]),"post_event_renewal_pass":bool(renewal),"rollback_pass":bool(rollback),"snapshot_restore_pass":bool(restored),"subdivision_status":"SOURCE_EXERCISED_NOT_CLAIMED_EXACT","rng_threshold_status":"B_IS_SERIALIZABLE_BASE_ENGINE_HAS_NO_STOCHASTIC_THRESHOLD_STREAM"}

def fem(repo):
 sys.path.insert(0,str(repo))
 from types import SimpleNamespace
 from arrhenius_fracture.mpz_parameterization_v911 import build_mpz_config,load_selected_row
 from arrhenius_fracture.moving_process_zone_v911 import MovingProcessZoneState
 row=load_selected_row(repo/"mpz_v9_11_parameters/DBTT/spatial_promotion_manifest.csv","DBTT");state=MovingProcessZoneState(build_mpz_config(SimpleNamespace(mpz_length_um=100.,mpz_n_bins=40,r_pz=1e-6),row));ad=FEMCZMProcessZoneStateAdapter();s0=ad.snapshot(state)
 rates=state._pt_model().rates(2e9,np.full(state.n_bins,1e14),1000,2.74e-10);no=state.evolve(0,1000,0,2.74e-10);state.mobile[:,:2]=1;before=state.mobile_count;em=state.evolve(1e-9,1000,3e9,2.74e-10);changed=state.mobile_count!=before
 ad.restore(state,s0);restored=ad.snapshot(state).fingerprint==s0.fingerprint;snap=ad.snapshot(state);state.advance(2*state.dx);renewal=state.advance_total_m>0;ad.restore(state,snap);rollback=ad.snapshot(state).fingerprint==snap.fingerprint
 return {"backend":"FEMCZM","production_class":"MovingProcessZoneState v9.11","frozen_rates_available":bool(rates),"no_event_update_available":bool(no),"emission_only_state_changed":bool(changed),"cleavage_crossing_status":"DRIVER_OWNED_NOT_MPZ_STATE","post_event_renewal_pass":bool(renewal),"rollback_pass":bool(rollback),"snapshot_restore_pass":bool(restored),"subdivision_status":"SOURCE_EXERCISED_NOT_CLAIMED_EXACT","rng_threshold_status":"DRIVER_OWNED_NOT_MPZ_STATE"}

def main():
 p=argparse.ArgumentParser();p.add_argument("backend",choices=("pf","fem"));p.add_argument("repo",type=Path);p.add_argument("out",type=Path);a=p.parse_args();result=pf(a.repo) if a.backend=="pf" else fem(a.repo);a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,indent=2)+"\n")
if __name__=="__main__":main()
