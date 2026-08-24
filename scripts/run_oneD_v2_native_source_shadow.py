#!/usr/bin/env python3
"""Isolated source-object shadow fixture; invoke once per backend package."""
from __future__ import annotations
import argparse,copy,csv,json,sys,tempfile
from pathlib import Path
from types import SimpleNamespace
import numpy as np

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from reduced_fracture_v2 import CanonicalParameters
from reduced_fracture_v2.adapters import PFProcessZoneStateAdapter,FEMCZMProcessZoneStateAdapter
from reduced_fracture_v2.mechanics import PFMechanicsProvider,FEMCZMMechanicsProvider
from reduced_fracture_v2.native_closure import PFNativeStateClosure,FEMCZMNativeStateClosure

IDS={"v913_zeroD_sobol_0242980":"Peak","v913_zeroD_sobol_0202500":"DBTT","v913_zeroD_sobol_0129902":"weak-T","v913_zeroD_sobol_0077080":"ceramic-like"}
OUT=ROOT/"analysis_outputs/oneD_v2_mechanics_maps_and_baselines"

def rows(repo):
 p=Path("/Volumes/Data/Data/Nanopillar_calculation/PF-fracture-fatigue_v10_2_21_persistent_sites_top1/arrhenius_fracture/data/materials/v10_2_27_v913_four_class_paper_registry.csv")
 with p.open(newline="") as f:return [r for r in csv.DictReader(f) if r["candidate_id"] in IDS]

def pf_engine(repo,row,exact=False):
 from arrhenius_fracture.config import ElasticProperties
 from arrhenius_fracture.material_manifest import MaterialManifest
 from arrhenius_fracture.parameter_registry_v9111 import select_option,write_compatibility_manifest
 from arrhenius_fracture.sharp_front import FrontConfig,default_cleavage_barrier,default_emission_barrier
 from arrhenius_fracture.unified_front import UnifiedMPZFrontEngine
 from arrhenius_fracture.unified_mpz import MPZConfig
 options={"v913_zeroD_sobol_0242980":"v913_paper_peak01_0242980_persistent_sites","v913_zeroD_sobol_0202500":"v913_paper_dbtt01_0202500_persistent_sites","v913_zeroD_sobol_0129902":"v913_paper_weakT01_0129902_persistent_sites","v913_zeroD_sobol_0077080":"v913_paper_ceramic01_0077080_persistent_sites"}
 registry=repo/"arrhenius_fracture/data/materials/v10_2_27_v913_four_class_paper_registry.csv"
 selected=select_option(options[row["candidate_id"]],registry,canonical_stage3_only=False)
 name=tempfile.NamedTemporaryFile(suffix=".csv",delete=False).name
 write_compatibility_manifest(selected,name)
 manifest=MaterialManifest.from_csv(name);mat=ElasticProperties();L=float(row["L_pz_um_recommended"])*1e-6;n=int(float(row["n_bins_recommended"]));fc=FrontConfig();fc.r0=1e-6;fc.sigma_cap=0.;fc.da=5e-6;fc.L_pz=L;fc.c_blunt=float(row["c_blunt"]);fc.max_advances_per_step=1
 cls=UnifiedMPZFrontEngine
 if exact:
  from arrhenius_fracture.persistent_site_source_v10221 import PersistentSiteConfig,PersistentSiteStateResolvedTipEngine
  from arrhenius_fracture.persistent_site_physical_width_v10222 import install_physical_front_width
  family=repo/"runs/v10_2_28_kernel_cache/1447653d199f0b43cb475951092d69444c9b785f6fdf518c723792abb3b1f5e5/family.json"
  PersistentSiteStateResolvedTipEngine.configure_state_resolved_physics(family)
  PersistentSiteStateResolvedTipEngine.configure_persistent_sites(PersistentSiteConfig(rho_site0_m2=float(row["rho_source0_m2"]),reference_source_area_m2=float(row["reference_source_area_um2"])*1e-12,reference_front_width_m=float(row["reference_front_width_um"])*1e-6,reference_density_m2=float(row["rho_forest_floor_m2"]),source_zone_length_m=float(row["source_zone_length_um"])*1e-6,maximum_front_width_m=L))
  install_physical_front_width();fc.sigma_cap=30e9;cls=PersistentSiteStateResolvedTipEngine
 return cls(fc,default_cleavage_barrier(),default_emission_barrier(mat.b),mat.G,mat.nu,mat.b,manifest,MPZConfig(length_m=L,n_bins=n,wake_length_m=L))

def fem_engine(repo,row,seed):
 from arrhenius_fracture.config import FractureBarrier
 from arrhenius_fracture.mpz_front_engine_v911 import MovingProcessZone2DFrontEngine
 from arrhenius_fracture.mpz_parameterization_v911 import build_mpz_config
 from arrhenius_fracture.sharp_front import FrontConfig
 cfg=build_mpz_config(SimpleNamespace(mpz_length_um=100.,mpz_n_bins=40,r_pz=1e-6),row)
 cfg.stochastic_seed=seed;cfg.event_statistics="deterministic"
 f=FrontConfig();f.r0=1e-6;f.sigma_cap=0.;f.da=5e-6;f.max_advances_per_step=1
 return MovingProcessZone2DFrontEngine(f,FractureBarrier(),FractureBarrier(),160e9,.28,2.74e-10,cfg)

def compare(backend,repo,row,T):
 p=CanonicalParameters.from_mapping(row);seed=8666 if IDS[row["candidate_id"]]=="Peak" else 1008666
 if backend in {"pf","pf_exact"}:
  eng=pf_engine(repo,row,exact=backend=="pf_exact");adapter=PFProcessZoneStateAdapter();provider=PFMechanicsProvider.from_csv(OUT/"oneD_v2_pf_mechanics_map.csv","pinned");closure=PFNativeStateClosure(provider,eng,adapter,p)
 else:
  eng=fem_engine(repo,row,seed);adapter=FEMCZMProcessZoneStateAdapter();provider=FEMCZMMechanicsProvider.from_csv(OUT/"oneD_v2_fem_native_mechanics_map.csv",OUT/"oneD_v2_fem_qualified_G_map.csv","pinned");closure=FEMCZMNativeStateClosure(provider,eng,adapter,p)
 source=copy.deepcopy(eng);records=[]
 for i,(opening,dt,force_event) in enumerate(((2e-6,0.,False),(8e-6,1e-9,False),(10e-6,0.,True),(11e-6,0.,False),(11e-6,0.,True),(11e-6,0.,True))):
  if force_event: eng.B=source.B=1.1
  pre=closure.state(opening);K=pre.native_KJ_Pa_sqrt_m
  try:
   result=closure.advance_interval(opening,T,dt)
   direct=source.step(K,T,dt) if backend in {"pf","pf_exact"} else source.step_drives(K,K,T,dt)
  except RuntimeError as exc:
   return [{"backend":"PF_EXACT","candidate_id":row["candidate_id"],"material_class":IDS[row["candidate_id"]],"temperature_K":T,"interval_index":i,"fixture_engine":type(source).__name__,"production_wrapper_exact":True,"pre_extension_m":pre.extension_m,"post_extension_m":pre.extension_m,"geometry_commits":0,"post_event_renewals":0,"process_zone_fingerprint_match":False,"action_match":False,"event_count_match":False,"event_length_match":False,"tip_radius_match":False,"tip_stress_match":False,"fired_match":False,"shared_core_source_rate_status":"UNAVAILABLE_RELIABLE_2D_TENSOR_DRIVE_REQUIRED","required_field_status":"UNAVAILABLE","closure_status":"EXACT_PRODUCTION_WRAPPER_REJECTED_SCALAR_NATIVE_KJ: "+str(exc)}]
  a=adapter.snapshot(closure._pz());b=adapter.snapshot(getattr(source,"mpz",getattr(source,"mpz_state",None)))
  records.append({"backend":backend.upper(),"candidate_id":row["candidate_id"],"material_class":IDS[row["candidate_id"]],"temperature_K":T,"interval_index":i,
   "fixture_engine":type(source).__name__,"production_wrapper_exact":bool(backend=="pf_exact" and type(source).__name__=="PersistentSiteStateResolvedTipEngine"),
   "pre_extension_m":result.pre_state.extension_m,"post_extension_m":result.post_state.extension_m,"geometry_commits":result.geometry_commits,"post_event_renewals":result.post_event_renewals,
   "process_zone_fingerprint_match":a.fingerprint==b.fingerprint,"action_match":float(eng.B)==float(source.B),"event_count_match":int(eng.n_adv)==int(source.n_adv),"event_length_match":float(eng.a_adv)==float(source.a_adv),
   "tip_radius_match":float(eng.r_eff())==float(source.r_eff()),"tip_stress_match":float(eng.sigma_tip(K))==float(source.sigma_tip(K)),"fired_match":bool(result.fired)==bool(direct["fired"]),
   "shared_core_source_rate_status":"UNAVAILABLE_EXACT_BARRIER_RATE_COMPARISON_NOT_IMPLEMENTED","required_field_status":"UNAVAILABLE","closure_status":"SOURCE_OBJECT_SHADOW_EXACT_BUT_RATE_PARITY_UNAVAILABLE"})
 return records

def main():
 q=argparse.ArgumentParser();q.add_argument("backend",choices=("pf","pf_exact","fem"));q.add_argument("repo",type=Path);q.add_argument("out",type=Path);a=q.parse_args();sys.path.insert(0,str(a.repo))
 result=[]
 for row in rows(a.repo):
  for T in (300.,1000.,1200.):result.extend(compare(a.backend,a.repo,row,T))
 a.out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
if __name__=="__main__":main()
