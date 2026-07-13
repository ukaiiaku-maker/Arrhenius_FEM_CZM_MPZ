#!/usr/bin/env python3
"""Deterministic anisotropic traction-controlled first-passage campaign."""
from __future__ import annotations
import argparse,csv,json,math,shlex,subprocess
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
import pandas as pd
REQ={"target_class","exp_G00_eV","exp_gT_eV_per_K","exp_sigc0_GPa","exp_sT_MPa_per_K","exp_a","exp_n","exp_floor_frac","cleave_G00_eV","cleave_gT_eV_per_K","cleave_sigc0_GPa","cleave_sT_MPa_per_K","cleave_exp_a","cleave_exp_n","cleave_floor_frac","cleave_S_hs_kB","chi_shield","N_sat"}
def items(s,cast=str):return [cast(x) for x in str(s).replace(',',' ').split() if x]
def fs(x):return 'inf' if math.isinf(float(x)) else f'{float(x):.16g}'
def tag(x):return ('m' if x<0 else 'p')+f'{abs(x):05.1f}'.replace('.','p')
def pyenv(e):
 cp=subprocess.run(['conda','run','-n',e,'python','-c','import sys;print(sys.executable)'],capture_output=True,text=True)
 if cp.returncode:raise SystemExit(cp.stderr)
 return [x for x in cp.stdout.splitlines() if x.strip()][-1]
def main():
 p=argparse.ArgumentParser();p.add_argument('--parameter-table',default='four_class_exp_floor_exact_model_inputs.csv');p.add_argument('--calibration-csv',required=True);p.add_argument('--classes',default='ceramic DBTT');p.add_argument('--target-psi-deg',default='-60 -45 -30 -15 0 15 30 45 60');p.add_argument('--T-K',type=float,default=500);p.add_argument('--outroot',default='runs/mixed_mode_fem_czm_v4_anisotropic_500K');p.add_argument('--conda-env',default='arrhenius-fem-czm');p.add_argument('--max-jobs',type=int,default=1);p.add_argument('--force',action='store_true');p.add_argument('--nx',type=int,default=24);p.add_argument('--ny',type=int,default=48);p.add_argument('--tip-h-fine',type=float,default=3e-6);p.add_argument('--tip-ratio',type=float,default=1.25);p.add_argument('--dU',type=float,default=2e-7);p.add_argument('--dt',type=float,default=8.4);p.add_argument('--steps',type=int,default=3000);p.add_argument('--print-every',type=int,default=100);p.add_argument('--crystal-theta-deg',type=float,default=0);p.add_argument('--crystal-C11',type=float,default=523e9);p.add_argument('--crystal-C12',type=float,default=203e9);p.add_argument('--crystal-C44',type=float,default=160e9);p.add_argument('--cleave-gamma-aniso',type=float,default=.3);p.add_argument('--traction-probe-radius-m',type=float,default=10e-6)
 a=p.parse_args();py=pyenv(a.conda_env);df=pd.read_csv(a.parameter_table);miss=REQ-set(df.columns)
 if miss:raise SystemExit(f'parameter table missing {sorted(miss)}')
 df=df.set_index('target_class',drop=False);cal={round(float(r['target_psi_deg']),6):r for r in csv.DictReader(open(a.calibration_csv))};out=Path(a.outroot);out.mkdir(parents=True,exist_ok=True)
 cp=subprocess.run([py,'-c','from arrhenius_fracture.mixed_mode_first_passage_v4 import MODEL_ID;print(MODEL_ID)'],capture_output=True,text=True)
 if cp.returncode:raise SystemExit(cp.stderr);print(cp.stdout)
 def run(k,target):
  c=cal[round(target,6)]
  if str(c.get('phase_converged','')).lower() not in {'1','true','yes'}:raise RuntimeError(f'calibration not converged {target}')
  r=df.loc[k];alpha=float(c['loading_angle_deg']);sign=float(c['traction_shear_sign']);case=out/k/f'psi_{tag(target)}';case.mkdir(parents=True,exist_ok=True);summ=case/'anisotropic_mixed_mode_first_passage_summary.json'
  if summ.exists() and not a.force:return json.loads(summ.read_text())|{'class':k,'target_psi_deg':target}
  emitG=.75*float(r.exp_G00_eV);emitg=.75*float(r.exp_gT_eV_per_K)
  cmd=[py,'-m','arrhenius_fracture.mixed_mode_first_passage_v4','--mixity-loading-angle-deg',fs(alpha),'--target-traction-phase-deg',fs(target),'--traction-shear-sign',fs(sign),'--traction-probe-radius-m',fs(a.traction_probe_radius_m),'--mode','2d','--nx',str(a.nx),'--ny',str(a.ny),'--tip-h-fine',fs(a.tip_h_fine),'--tip-ratio',fs(a.tip_ratio),'--dU',fs(a.dU),'--dt',fs(a.dt),'--steps',str(a.steps),'--n-stagger','2','--print-every',str(a.print_every),'--stop-after-first-fire','--max-fronts','1','--adaptive-events','--adaptive-event-target','.25','--adaptive-min-frac','1e-8','--adaptive-grow','4','--da-phys','5e-6','--j-decomposition','cluster','--rJ-cluster','20e-6','--rJ-outer','25e-6','--temperatures',fs(a.T_K),'--crack-backend','adaptive_czm','--czm-max-angle-error-deg','35','--crystal-aniso','--crystal-compete','--crystal-theta-deg',fs(a.crystal_theta_deg),'--crystal-C11',fs(a.crystal_C11),'--crystal-C12',fs(a.crystal_C12),'--crystal-C44',fs(a.crystal_C44),'--cleave-gamma-aniso',fs(a.cleave_gamma_aniso),'--crystal-material','w','--emit-barrier-kind','exp_floor','--emit-G00-eV',fs(emitG),'--emit-gT-eV-per-K',fs(emitg),'--emit-sigc0-GPa',fs(r.exp_sigc0_GPa),'--emit-sT-GPa-per-K',fs(float(r.exp_sT_MPa_per_K)/1000),'--emit-exp-a',fs(r.exp_a),'--emit-exp-n',fs(r.exp_n),'--emit-floor-frac',fs(r.exp_floor_frac),'--emit-Tref-K','300','--cleave-barrier-kind','exp_floor','--cleave-exp-T-mode','linear','--cleave-G00-eV',fs(r.cleave_G00_eV),'--cleave-gT-eV-per-K',fs(r.cleave_gT_eV_per_K),'--cleave-sigc0-GPa',fs(r.cleave_sigc0_GPa),'--cleave-sT-GPa-per-K',fs(float(r.cleave_sT_MPa_per_K)/1000),'--cleave-exp-a',fs(r.cleave_exp_a),'--cleave-exp-n',fs(r.cleave_exp_n),'--cleave-floor-frac',fs(r.cleave_floor_frac),'--cleave-S-hs-kB',fs(r.cleave_S_hs_kB),'--cleave-sigma-S-GPa','6','--cleave-S-hs-power','2','--cleave-S-hs-Tref-K','300','--cleave-Tref-K','300','--cleave-shield-chi',fs(r.chi_shield),'--n-sat',fs(r.N_sat),'--multihit-m','3','--multihit-tau','1e-6','--emb-sat-frac','1','--save-snapshots','0','--no-plots','--out',str(case)]
  (case/'command.txt').write_text(shlex.join(cmd)+'\n')
  with (case/'run.log').open('w') as log:rc=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT).returncode
  if rc:raise RuntimeError(f'case failed rc={rc}: {case}')
  z=json.loads(summ.read_text());z.update({'class':k,'target_psi_deg':target,'calibrated_loading_angle_deg':alpha,'calibration_phase_deg':float(c['achieved_traction_phase_deg']),'calibration_phase_error_deg':float(c['traction_phase_error_deg']),'calibration_zener_ratio':float(c['zener_ratio'])});return z
 jobs=[(k,t) for k in items(a.classes) for t in items(a.target_psi_deg,float)];res=[]
 with ThreadPoolExecutor(max_workers=max(1,a.max_jobs)) as ex:
  fut={ex.submit(run,k,t):(k,t) for k,t in jobs}
  for f in as_completed(fut):
   k,t=fut[f]
   try:z=f.result();z['status']='event' if z.get('control_state')=='first_passage' else 'right_censored'
   except Exception as e:z={'class':k,'target_psi_deg':t,'status':'failed','error':repr(e)}
   res.append(z);print({q:z.get(q) for q in ('class','target_psi_deg','status','KJ_reference_first_MPa_sqrt_m','traction_phase_first_deg','candidate_angle_first_deg','mode_classification')})
 pd.DataFrame(res).to_csv(out/'campaign_status_v4.csv',index=False);good=[x for x in res if x['status']!='failed']
 if good:pd.DataFrame(good).to_csv(out/'mixed_mode_v4_anisotropic_all_cases.csv',index=False)
 if any(x['status']=='failed' for x in res):raise SystemExit('one or more v4 cases failed')
if __name__=='__main__':main()
