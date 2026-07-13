#!/usr/bin/env python3
"""Run short mixed-mode FEM/CZM first-passage cases."""
from __future__ import annotations
import argparse,csv,json,math,os,shlex,subprocess,sys
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
import pandas as pd

REQ={"target_class","exp_G00_eV","exp_gT_eV_per_K","exp_sigc0_GPa","exp_sT_MPa_per_K","exp_a","exp_n","exp_floor_frac","cleave_G00_eV","cleave_gT_eV_per_K","cleave_sigc0_GPa","cleave_sT_MPa_per_K","cleave_exp_a","cleave_exp_n","cleave_floor_frac","cleave_S_hs_kB","chi_shield","N_sat"}
def items(s,cast=str): return [cast(x) for x in str(s).replace(',',' ').split() if x]
def fs(x):
    x=float(x); return 'inf' if math.isinf(x) else f'{x:.16g}'
def resolve_table(path):
    if path:
        p=Path(path); 
        if not p.exists(): raise SystemExit(f'parameter table not found: {p}')
        return p
    for n in ['four_class_exp_floor_exact_model_inputs.csv','four_class_exp_floor_final_parameters.csv','four_class_exp_floor_model_inputs.csv']:
        if Path(n).exists(): return Path(n)
    raise SystemExit('No four-class parameter table found. Supply --parameter-table PATH.')
def resolve_python(a):
    if a.python_bin:return str(Path(a.python_bin).expanduser().resolve())
    cp=subprocess.run(['conda','run','-n',a.conda_env,'python','-c','import sys;print(sys.executable)'],text=True,capture_output=True)
    if cp.returncode: raise SystemExit(cp.stderr)
    return [x for x in cp.stdout.splitlines() if x.strip()][-1]
def tag_angle(x): return ('m' if x<0 else 'p')+f'{abs(x):05.1f}'.replace('.','p')

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--parameter-table',default='');p.add_argument('--calibration-csv',default='')
    p.add_argument('--classes',default='ceramic DBTT');p.add_argument('--target-psi-deg',default='-60 -45 -30 -15 0 15 30 45 60')
    p.add_argument('--seeds',default='1101 1102 1103');p.add_argument('--T-K',type=float,default=500)
    p.add_argument('--outroot',default='runs/mixed_mode_fem_czm_v1_first_passage_500K')
    p.add_argument('--conda-env',default='arrhenius-fem-czm');p.add_argument('--python-bin',default='')
    p.add_argument('--max-jobs',type=int,default=1);p.add_argument('--force',action='store_true')
    p.add_argument('--nx',type=int,default=18);p.add_argument('--ny',type=int,default=36)
    p.add_argument('--tip-h-fine',type=float,default=5e-6);p.add_argument('--tip-ratio',type=float,default=1.30)
    p.add_argument('--dU',type=float,default=2e-7);p.add_argument('--dt',type=float,default=8.4);p.add_argument('--steps',type=int,default=3000)
    p.add_argument('--theta-deg',type=float,default=45);p.add_argument('--isotropic',action='store_true')
    p.add_argument('--crystal-material',default='w');p.add_argument('--cleave-gamma-aniso',type=float,default=0.3);p.add_argument('--shear-emission-weight',type=float,default=1.0)
    p.add_argument('--deterministic-threshold',action='store_true')
    p.add_argument('--rJ-cluster',type=float,default=20e-6);p.add_argument('--rJ-outer',type=float,default=25e-6)
    p.add_argument('--print-every',type=int,default=100);p.add_argument('--save-snapshots',type=int,default=0)
    a=p.parse_args();py=resolve_python(a);table=resolve_table(a.parameter_table)
    df=pd.read_csv(table);miss=REQ-set(df.columns)
    if miss:raise SystemExit(f'parameter table missing: {sorted(miss)}')
    df=df.set_index('target_class',drop=False)
    cal={}
    if a.calibration_csv:
        for r in csv.DictReader(open(a.calibration_csv)):
            cal[round(float(r['target_psi_deg']),6)]=float(r['loading_angle_deg'])
    cp=subprocess.run([py,'-m','arrhenius_fracture.mixed_mode_first_passage_v1','--help'],text=True,capture_output=True)
    if cp.returncode:raise SystemExit('mixed-mode module preflight failed:\n'+cp.stderr)
    outroot=Path(a.outroot);outroot.mkdir(parents=True,exist_ok=True)
    jobs=[]
    for klass in items(a.classes):
        if klass not in df.index:raise SystemExit(f'class {klass!r} absent from {table}')
        r=df.loc[klass]
        for psi in items(a.target_psi_deg,float):
            alpha=cal.get(round(psi,6),psi)
            for seed in items(a.seeds,int):
                out=outroot/klass/f'psi_{tag_angle(psi)}'/f'seed_{seed}'
                jobs.append((klass,psi,alpha,seed,r,out))
    def run(job):
        klass,psi,alpha,seed,r,out=job;out.mkdir(parents=True,exist_ok=True)
        if (out/'mixed_mode_first_passage_summary.json').exists() and not a.force:return {'status':'skipped','case':str(out)}
        emitG=.75*float(r.exp_G00_eV);emitg=.75*float(r.exp_gT_eV_per_K)
        cmd=[py,'-m','arrhenius_fracture.mixed_mode_first_passage_v1',
             '--mixity-loading-angle-deg',fs(alpha),'--solver-seed',str(seed),'--shear-emission-weight',fs(a.shear_emission_weight),
             '--mode','2d','--nx',str(a.nx),'--ny',str(a.ny),'--tip-h-fine',fs(a.tip_h_fine),'--tip-ratio',fs(a.tip_ratio),
             '--dU',fs(a.dU),'--dt',fs(a.dt),'--steps',str(a.steps),'--n-stagger','2','--print-every',str(a.print_every),
             '--stop-after-first-fire','--max-fronts','1','--adaptive-events','--adaptive-event-target','0.25','--adaptive-min-frac','1e-8','--adaptive-grow','4',
             '--da-phys','5e-6','--j-decomposition','cluster','--rJ-cluster',fs(a.rJ_cluster),'--rJ-outer',fs(a.rJ_outer),
             '--temperatures',fs(a.T_K),'--crystal-theta-deg',fs(a.theta_deg),'--crack-backend','adaptive_czm','--czm-max-angle-error-deg','35',
             '--emit-barrier-kind','exp_floor','--emit-G00-eV',fs(emitG),'--emit-gT-eV-per-K',fs(emitg),'--emit-sigc0-GPa',fs(r.exp_sigc0_GPa),
             '--emit-sT-GPa-per-K',fs(float(r.exp_sT_MPa_per_K)/1000),'--emit-exp-a',fs(r.exp_a),'--emit-exp-n',fs(r.exp_n),'--emit-floor-frac',fs(r.exp_floor_frac),'--emit-Tref-K','300',
             '--cleave-barrier-kind','exp_floor','--cleave-exp-T-mode','linear','--cleave-G00-eV',fs(r.cleave_G00_eV),'--cleave-gT-eV-per-K',fs(r.cleave_gT_eV_per_K),
             '--cleave-sigc0-GPa',fs(r.cleave_sigc0_GPa),'--cleave-sT-GPa-per-K',fs(float(r.cleave_sT_MPa_per_K)/1000),'--cleave-exp-a',fs(r.cleave_exp_a),'--cleave-exp-n',fs(r.cleave_exp_n),
             '--cleave-floor-frac',fs(r.cleave_floor_frac),'--cleave-S-hs-kB',fs(r.cleave_S_hs_kB),'--cleave-sigma-S-GPa','6','--cleave-S-hs-power','2','--cleave-S-hs-Tref-K','300','--cleave-Tref-K','300',
             '--cleave-shield-chi',fs(r.chi_shield),'--n-sat',fs(r.N_sat),'--multihit-m','3','--multihit-tau','1e-6','--emb-sat-frac','1',
             '--out',str(out)]
        if not a.isotropic:cmd += ['--crystal-aniso','--crystal-compete','--crystal-material',a.crystal_material,'--cleave-gamma-aniso',fs(a.cleave_gamma_aniso)]
        if a.deterministic_threshold:cmd += ['--deterministic-threshold']
        if a.save_snapshots>0:cmd += ['--save-snapshots',str(a.save_snapshots),'--snapshot-cols','4']
        else:cmd += ['--save-snapshots','0','--no-plots']
        (out/'command.txt').write_text(shlex.join(cmd)+'\n')
        (out/'campaign_metadata.json').write_text(json.dumps({'class':klass,'target_psi_deg':psi,'loading_angle_deg':alpha,'solver_seed':seed,'T_K':a.T_K,'parameter_table':str(table)},indent=2))
        with (out/'run.log').open('w') as log: rc=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT).returncode
        return {'status':'complete' if rc==0 else 'failed','returncode':rc,'case':str(out)}
    rows=[]
    with ThreadPoolExecutor(max_workers=max(1,a.max_jobs)) as ex:
        fut={ex.submit(run,j):j for j in jobs}
        for f in as_completed(fut):
            z=f.result();rows.append(z);print(z)
    with (outroot/'campaign_status.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=sorted({k for r in rows for k in r}));w.writeheader();w.writerows(rows)
    if any(r.get('status')=='failed' for r in rows):raise SystemExit('one or more mixed-mode cases failed')
if __name__=='__main__':main()
