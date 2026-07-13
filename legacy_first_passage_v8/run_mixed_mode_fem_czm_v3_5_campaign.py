#!/usr/bin/env python3
"""One-run-per-condition v3.5 campaign runner for the validated v3.3 mechanics."""
from __future__ import annotations
import argparse,csv,json,math,shlex,subprocess
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
import pandas as pd

REQ={"target_class","exp_G00_eV","exp_gT_eV_per_K","exp_sigc0_GPa","exp_sT_MPa_per_K","exp_a","exp_n","exp_floor_frac","cleave_G00_eV","cleave_gT_eV_per_K","cleave_sigc0_GPa","cleave_sT_MPa_per_K","cleave_exp_a","cleave_exp_n","cleave_floor_frac","cleave_S_hs_kB","chi_shield","N_sat"}
def items(s,cast=str):return [cast(x) for x in str(s).replace(',',' ').split() if x]
def fs(x):return 'inf' if math.isinf(float(x)) else f'{float(x):.16g}'
def tag(x):return ('m' if x<0 else 'p')+f'{abs(x):05.1f}'.replace('.','p')
def resolve_python(env):
    cp=subprocess.run(['conda','run','-n',env,'python','-c','import sys;print(sys.executable)'],text=True,capture_output=True)
    if cp.returncode:raise SystemExit(cp.stderr)
    return [x for x in cp.stdout.splitlines() if x.strip()][-1]
def resolve_table(path):
    p=Path(path)
    if not p.exists():raise SystemExit(f'parameter table not found: {p}')
    return p
def load_cal(path):
    out={}
    for r in csv.DictReader(open(path)):
        out[round(float(r['target_psi_deg']),6)]=r
    return out

EXPECTED_MECHANICS_MODEL = "FEM_CZM_mixed_mode_first_passage_v3_3_J_consistent_circular_phase"

def driver_probe_command(py):
    """Import-check the installed mechanics without invoking its required CLI."""
    code = (
        "from arrhenius_fracture.mixed_mode_first_passage_v3_3 import MODEL_ID;"
        "print(MODEL_ID)"
    )
    return [py, "-c", code]

def required_mixed_args(alpha, target, shear_emission_weight):
    """Arguments that must precede the base sharp-front arguments."""
    return [
        "--mixity-loading-angle-deg", fs(alpha),
        "--target-mode-phase-deg", fs(target),
        "--deterministic-threshold",
        "--shear-emission-weight", fs(shear_emission_weight),
    ]

def probe_driver(py):
    cp = subprocess.run(driver_probe_command(py), text=True, capture_output=True)
    if cp.returncode:
        raise SystemExit(cp.stderr or cp.stdout)
    lines = [x.strip() for x in cp.stdout.splitlines() if x.strip()]
    model = lines[-1] if lines else ""
    if model != EXPECTED_MECHANICS_MODEL:
        raise SystemExit(
            f"wrong mixed-mode mechanics model: {model!r}; "
            f"expected {EXPECTED_MECHANICS_MODEL!r}"
        )
    print("mixed-mode mechanics:", model)

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--parameter-table',default='four_class_exp_floor_exact_model_inputs.csv')
    p.add_argument('--calibration-csv',required=True)
    p.add_argument('--classes',default='ceramic DBTT');p.add_argument('--target-psi-deg',default='-60 -45 -30 -15 0 15 30 45 60')
    p.add_argument('--T-K',type=float,default=500);p.add_argument('--outroot',default='runs/mixed_mode_fem_czm_v3_5_campaign_500K')
    p.add_argument('--conda-env',default='arrhenius-fem-czm');p.add_argument('--max-jobs',type=int,default=1);p.add_argument('--force',action='store_true')
    p.add_argument('--nx',type=int,default=24);p.add_argument('--ny',type=int,default=48);p.add_argument('--tip-h-fine',type=float,default=3e-6);p.add_argument('--tip-ratio',type=float,default=1.25)
    p.add_argument('--dU',type=float,default=2e-7);p.add_argument('--dt',type=float,default=8.4);p.add_argument('--steps',type=int,default=3000)
    p.add_argument('--shear-emission-weight',type=float,default=1.0);p.add_argument('--print-every',type=int,default=100);p.add_argument('--save-snapshots',type=int,default=0)
    a=p.parse_args();py=resolve_python(a.conda_env);table=resolve_table(a.parameter_table);cal=load_cal(a.calibration_csv)
    df=pd.read_csv(table);miss=REQ-set(df.columns)
    if miss:raise SystemExit(f'parameter table missing: {sorted(miss)}')
    df=df.set_index('target_class',drop=False);outroot=Path(a.outroot);outroot.mkdir(parents=True,exist_ok=True)
    probe_driver(py)

    def run_case(klass,target):
        key=round(target,6)
        if key not in cal:raise RuntimeError(f'no calibration row for target {target}')
        c=cal[key]
        if str(c.get('phase_converged','')).lower() not in {'1','true','yes'}:
            raise RuntimeError(f'isotropic phase-ratio calibration not converged for target {target}')
        if klass not in df.index:raise RuntimeError(f'class {klass!r} absent from table')
        r=df.loc[klass];alpha=float(c['loading_angle_deg']);case=outroot/klass/f'psi_{tag(target)}';case.mkdir(parents=True,exist_ok=True)
        summary=case/'mixed_mode_first_passage_summary.json'
        if summary.exists() and not a.force:return json.loads(summary.read_text())|{'class':klass,'target_psi_deg':target}
        emitG=.75*float(r.exp_G00_eV);emitg=.75*float(r.exp_gT_eV_per_K)
        cmd=[py,'-m','arrhenius_fracture.mixed_mode_first_passage_v3_3'] + required_mixed_args(alpha,target,a.shear_emission_weight) + [
             '--mode','2d','--nx',str(a.nx),'--ny',str(a.ny),'--tip-h-fine',fs(a.tip_h_fine),'--tip-ratio',fs(a.tip_ratio),'--dU',fs(a.dU),'--dt',fs(a.dt),'--steps',str(a.steps),'--n-stagger','2','--print-every',str(a.print_every),
             '--stop-after-first-fire','--max-fronts','1','--adaptive-events','--adaptive-event-target','0.25','--adaptive-min-frac','1e-8','--adaptive-grow','4','--da-phys','5e-6','--j-decomposition','cluster','--rJ-cluster','20e-6','--rJ-outer','25e-6',
             '--temperatures',fs(a.T_K),'--crack-backend','adaptive_czm','--czm-max-angle-error-deg','35',
             '--emit-barrier-kind','exp_floor','--emit-G00-eV',fs(emitG),'--emit-gT-eV-per-K',fs(emitg),'--emit-sigc0-GPa',fs(r.exp_sigc0_GPa),'--emit-sT-GPa-per-K',fs(float(r.exp_sT_MPa_per_K)/1000),'--emit-exp-a',fs(r.exp_a),'--emit-exp-n',fs(r.exp_n),'--emit-floor-frac',fs(r.exp_floor_frac),'--emit-Tref-K','300',
             '--cleave-barrier-kind','exp_floor','--cleave-exp-T-mode','linear','--cleave-G00-eV',fs(r.cleave_G00_eV),'--cleave-gT-eV-per-K',fs(r.cleave_gT_eV_per_K),'--cleave-sigc0-GPa',fs(r.cleave_sigc0_GPa),'--cleave-sT-GPa-per-K',fs(float(r.cleave_sT_MPa_per_K)/1000),'--cleave-exp-a',fs(r.cleave_exp_a),'--cleave-exp-n',fs(r.cleave_exp_n),'--cleave-floor-frac',fs(r.cleave_floor_frac),'--cleave-S-hs-kB',fs(r.cleave_S_hs_kB),'--cleave-sigma-S-GPa','6','--cleave-S-hs-power','2','--cleave-S-hs-Tref-K','300','--cleave-Tref-K','300','--cleave-shield-chi',fs(r.chi_shield),'--n-sat',fs(r.N_sat),'--multihit-m','3','--multihit-tau','1e-6','--emb-sat-frac','1','--out',str(case)]
        cmd += ['--save-snapshots',str(a.save_snapshots)] if a.save_snapshots>0 else ['--save-snapshots','0','--no-plots']
        (case/'command.txt').write_text(shlex.join(cmd)+'\n')
        with (case/'run.log').open('w') as log:rc=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT).returncode
        if rc:raise RuntimeError(f'case failed rc={rc}: {case}')
        z=json.loads(summary.read_text());z.update(class_=klass,target_psi_deg=target,calibrated_loading_angle_deg=alpha)
        z.update(
            calibration_achieved_psi_deg=float(c.get('achieved_psi_deg','nan')),
            calibration_psi_error_deg=float(c.get('psi_error_deg','nan')),
            calibration_phase_spread_deg=float(c.get('projection_psi_spread_deg','nan')),
            calibration_phase_spread_warning=str(c.get('phase_spread_warning','')).lower() in {'1','true','yes'},
            calibration_confidence=str(c.get('calibration_confidence','')),
            calibration_basis_condition=float(c.get('basis_condition','nan')),
        )
        z['class']=z.pop('class_');return z

    jobs=[(k,p) for k in items(a.classes) for p in items(a.target_psi_deg,float)];results=[]
    with ThreadPoolExecutor(max_workers=max(1,a.max_jobs)) as ex:
        fut={ex.submit(run_case,k,p):(k,p) for k,p in jobs}
        for f in as_completed(fut):
            k,p=fut[f]
            try:
                z=f.result();z['status']='event' if z.get('control_state')=='first_passage' else 'right_censored'
            except Exception as exc:z={'class':k,'target_psi_deg':p,'status':'failed','error':repr(exc)}
            results.append(z);print({q:z.get(q) for q in ['class','target_psi_deg','status','KJ_first_MPa_sqrt_m','Kopen_maxhoop_first_MPa_sqrt_m','mode_classification']})
    pd.DataFrame(results).to_csv(outroot/'campaign_status_v3_5.csv',index=False)
    good=[r for r in results if r['status']!='failed']
    if good:pd.DataFrame(good).to_csv(outroot/'mixed_mode_v3_5_all_cases.csv',index=False)
    if any(r['status']=='failed' for r in results):raise SystemExit('one or more v3 cases failed')

if __name__=='__main__':main()
