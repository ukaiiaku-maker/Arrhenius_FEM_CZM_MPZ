#!/usr/bin/env python3
"""Deterministic, event-state-controlled mixed-mode FEM/CZM campaign."""
from __future__ import annotations
import argparse,csv,json,math,shlex,subprocess
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
import numpy as np
import pandas as pd

REQ={"target_class","exp_G00_eV","exp_gT_eV_per_K","exp_sigc0_GPa","exp_sT_MPa_per_K","exp_a","exp_n","exp_floor_frac","cleave_G00_eV","cleave_gT_eV_per_K","cleave_sigc0_GPa","cleave_sT_MPa_per_K","cleave_exp_a","cleave_exp_n","cleave_floor_frac","cleave_S_hs_kB","chi_shield","N_sat"}
def items(s,cast=str):return [cast(x) for x in str(s).replace(',',' ').split() if x]
def fs(x):
    x=float(x);return 'inf' if math.isinf(x) else f'{x:.16g}'
def tag_angle(x):return ('m' if x<0 else 'p')+f'{abs(x):05.1f}'.replace('.','p')
def tag_alpha(x):return ('m' if x<0 else 'p')+f'{abs(x):07.3f}'.replace('.','p')
def angle_error(a,b):return (float(a)-float(b)+180.0)%360.0-180.0

def resolve_table(path):
    if path:
        p=Path(path)
        if not p.exists():raise SystemExit(f'parameter table not found: {p}')
        return p
    for n in ['four_class_exp_floor_exact_model_inputs.csv','four_class_exp_floor_final_parameters.csv','four_class_exp_floor_model_inputs.csv']:
        if Path(n).exists():return Path(n)
    raise SystemExit('No four-class parameter table found. Supply --parameter-table PATH.')

def resolve_python(a):
    if a.python_bin:return str(Path(a.python_bin).expanduser().resolve())
    cp=subprocess.run(['conda','run','-n',a.conda_env,'python','-c','import sys;print(sys.executable)'],text=True,capture_output=True)
    if cp.returncode:raise SystemExit(cp.stderr)
    return [x for x in cp.stdout.splitlines() if x.strip()][-1]

def load_calibration(path):
    out={}
    if not path:return out
    for r in csv.DictReader(open(path)):
        out[round(float(r['target_psi_deg']),6)]={
            'alpha':float(r['loading_angle_deg']),
            'slope':float(r.get('dpsi_dalpha',1.0) or 1.0),
            'reliable':str(r.get('projection_reliable','True')).lower() in {'1','true','yes'},
        }
    return out

def next_alpha(history,initial_slope,amin,amax,max_step):
    good=[h for h in history if np.isfinite(h['error_deg'])]
    if not good:return 0.0
    latest=good[-1]
    # Prefer a bracket and use regula falsi.
    brackets=[]
    for i,a in enumerate(good):
        for b in good[i+1:]:
            if a['error_deg']*b['error_deg']<=0 and abs(a['alpha_deg']-b['alpha_deg'])>1e-6:
                brackets.append((abs(a['alpha_deg']-b['alpha_deg']),a,b))
    if brackets:
        _,a,b=min(brackets,key=lambda x:x[0])
        den=b['error_deg']-a['error_deg']
        raw=0.5*(a['alpha_deg']+b['alpha_deg']) if abs(den)<1e-10 else a['alpha_deg']-a['error_deg']*(b['alpha_deg']-a['alpha_deg'])/den
    elif len(good)>=2:
        a,b=good[-2],good[-1];den=b['error_deg']-a['error_deg']
        slope=den/(b['alpha_deg']-a['alpha_deg']) if abs(b['alpha_deg']-a['alpha_deg'])>1e-8 else initial_slope
        if not np.isfinite(slope) or abs(slope)<0.08:slope=initial_slope
        raw=b['alpha_deg']-b['error_deg']/slope
    else:
        slope=initial_slope if np.isfinite(initial_slope) and abs(initial_slope)>=0.08 else math.copysign(0.5,initial_slope or 1.0)
        raw=latest['alpha_deg']-latest['error_deg']/slope
    step=float(np.clip(raw-latest['alpha_deg'],-max_step,max_step))
    candidate=float(np.clip(latest['alpha_deg']+step,amin,amax))
    if any(abs(candidate-h['alpha_deg'])<0.02 for h in good):
        candidate=float(np.clip(candidate+math.copysign(0.25,-latest['error_deg']*(initial_slope or 1.0)),amin,amax))
    return candidate

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--parameter-table',default='');p.add_argument('--calibration-csv',required=True)
    p.add_argument('--classes',default='ceramic DBTT');p.add_argument('--target-psi-deg',default='-60 -45 -30 -15 0 15 30 45 60')
    p.add_argument('--T-K',type=float,default=500);p.add_argument('--solver-seed',type=int,default=1)
    p.add_argument('--outroot',default='runs/mixed_mode_fem_czm_v2_event_controlled_500K')
    p.add_argument('--conda-env',default='arrhenius-fem-czm');p.add_argument('--python-bin',default='')
    p.add_argument('--max-jobs',type=int,default=1);p.add_argument('--force',action='store_true')
    p.add_argument('--psi-tol-deg',type=float,default=2.0);p.add_argument('--max-control-iters',type=int,default=5)
    p.add_argument('--alpha-min-deg',type=float,default=-72.0);p.add_argument('--alpha-max-deg',type=float,default=72.0);p.add_argument('--max-alpha-step-deg',type=float,default=15.0)
    p.add_argument('--allow-unreliable-projection',action='store_true')
    p.add_argument('--nx',type=int,default=18);p.add_argument('--ny',type=int,default=36)
    p.add_argument('--tip-h-fine',type=float,default=5e-6);p.add_argument('--tip-ratio',type=float,default=1.30)
    p.add_argument('--dU',type=float,default=2e-7);p.add_argument('--dt',type=float,default=8.4);p.add_argument('--steps',type=int,default=3000)
    p.add_argument('--theta-deg',type=float,default=45);p.add_argument('--isotropic',action='store_true')
    p.add_argument('--crystal-material',default='w');p.add_argument('--cleave-gamma-aniso',type=float,default=0.3);p.add_argument('--shear-emission-weight',type=float,default=1.0)
    p.add_argument('--rJ-cluster',type=float,default=20e-6);p.add_argument('--rJ-outer',type=float,default=25e-6)
    p.add_argument('--mode-projection-angle-deg',type=float,default=105.0);p.add_argument('--mode-projection-damage-cutoff',type=float,default=0.85)
    p.add_argument('--print-every',type=int,default=100);p.add_argument('--save-snapshots',type=int,default=0)
    a=p.parse_args();py=resolve_python(a);table=resolve_table(a.parameter_table);cal=load_calibration(a.calibration_csv)
    df=pd.read_csv(table);miss=REQ-set(df.columns)
    if miss:raise SystemExit(f'parameter table missing: {sorted(miss)}')
    df=df.set_index('target_class',drop=False)
    cp=subprocess.run([py,'-m','arrhenius_fracture.mixed_mode_first_passage_v2','--help'],text=True,capture_output=True)
    if cp.returncode:raise SystemExit('mixed-mode v2 module preflight failed:\n'+cp.stderr)
    outroot=Path(a.outroot);outroot.mkdir(parents=True,exist_ok=True)

    def run_trial(klass,target,alpha,it,r,case_root):
        out=case_root/'trials'/f'iter_{it:02d}_alpha_{tag_alpha(alpha)}';out.mkdir(parents=True,exist_ok=True)
        summary_path=out/'mixed_mode_first_passage_summary.json'
        if summary_path.exists() and not a.force:return json.loads(summary_path.read_text())|{'trial_dir':str(out)}
        emitG=.75*float(r.exp_G00_eV);emitg=.75*float(r.exp_gT_eV_per_K)
        cmd=[py,'-m','arrhenius_fracture.mixed_mode_first_passage_v2',
             '--mixity-loading-angle-deg',fs(alpha),'--solver-seed',str(a.solver_seed),'--deterministic-threshold','--shear-emission-weight',fs(a.shear_emission_weight),
             '--mode-projection-angle-deg',fs(a.mode_projection_angle_deg),'--mode-projection-damage-cutoff',fs(a.mode_projection_damage_cutoff),
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
             '--cleave-shield-chi',fs(r.chi_shield),'--n-sat',fs(r.N_sat),'--multihit-m','3','--multihit-tau','1e-6','--emb-sat-frac','1','--out',str(out)]
        if not a.isotropic:cmd += ['--crystal-aniso','--crystal-compete','--crystal-material',a.crystal_material,'--cleave-gamma-aniso',fs(a.cleave_gamma_aniso)]
        if a.save_snapshots>0:cmd += ['--save-snapshots',str(a.save_snapshots),'--snapshot-cols','4']
        else:cmd += ['--save-snapshots','0','--no-plots']
        (out/'command.txt').write_text(shlex.join(cmd)+'\n')
        (out/'control_trial_metadata.json').write_text(json.dumps({'class':klass,'target_psi_deg':target,'loading_angle_deg':alpha,'iteration':it,'deterministic_threshold_H':1.0,'T_K':a.T_K},indent=2))
        with (out/'run.log').open('w') as log:rc=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT).returncode
        if rc!=0:raise RuntimeError(f'trial failed rc={rc}: {out}')
        payload=json.loads(summary_path.read_text());payload['trial_dir']=str(out);return payload

    def control_case(klass,target):
        if klass not in df.index:raise RuntimeError(f'class {klass!r} absent from {table}')
        r=df.loc[klass];case_root=outroot/klass/f'psi_{tag_angle(target)}';case_root.mkdir(parents=True,exist_ok=True)
        final_path=case_root/'mixed_mode_v2_final_summary.json'
        if final_path.exists() and not a.force:return json.loads(final_path.read_text())
        c=cal.get(round(target,6),{'alpha':-target,'slope':-1.0,'reliable':False})
        alpha=float(np.clip(c['alpha'],a.alpha_min_deg,a.alpha_max_deg));slope=float(c['slope'])
        history=[]
        for it in range(a.max_control_iters):
            payload=run_trial(klass,target,alpha,it,r,case_root)
            achieved=float(payload.get('mode_phase_first_deg',np.nan));reliable=bool(payload.get('projection_reliable',False))
            err=angle_error(achieved,target) if np.isfinite(achieved) else float('nan')
            row={'iteration':it,'alpha_deg':alpha,'target_psi_deg':target,'achieved_psi_deg':achieved,'error_deg':err,
                 'projection_reliable':reliable,'projection_rel_rmse':payload.get('projection_rel_rmse'),'projection_psi_spread_deg':payload.get('projection_psi_spread_deg'),
                 'control_state':payload.get('control_state'),'mode_classification':payload.get('mode_classification'),'trial_dir':payload['trial_dir']}
            history.append(row)
            if np.isfinite(err) and abs(err)<=a.psi_tol_deg and (reliable or a.allow_unreliable_projection):break
            alpha=next_alpha(history,slope,a.alpha_min_deg,a.alpha_max_deg,a.max_alpha_step_deg)
        candidates=[h for h in history if np.isfinite(h['error_deg']) and (h['projection_reliable'] or a.allow_unreliable_projection)]
        if not candidates:candidates=[h for h in history if np.isfinite(h['error_deg'])]
        if not candidates:raise RuntimeError(f'no valid mode projection for {klass}, target {target}')
        best=min(candidates,key=lambda h:abs(h['error_deg']))
        selected=json.loads((Path(best['trial_dir'])/'mixed_mode_first_passage_summary.json').read_text())
        converged=bool(abs(best['error_deg'])<=a.psi_tol_deg and (best['projection_reliable'] or a.allow_unreliable_projection))
        final={**selected,'class':klass,'target_psi_deg':target,'achieved_psi_deg':best['achieved_psi_deg'],'psi_error_deg':best['error_deg'],
               'event_state_control_converged':converged,'control_iterations':len(history),'selected_loading_angle_deg':best['alpha_deg'],
               'selected_trial_dir':best['trial_dir'],'deterministic_threshold_H':1.0,'elastic_initial_alpha_deg':c['alpha'],'elastic_initial_slope_dpsi_dalpha':slope}
        (case_root/'mixed_mode_control_history.csv').write_text(pd.DataFrame(history).to_csv(index=False))
        final_path.write_text(json.dumps(final,indent=2,default=str))
        pd.DataFrame([final]).to_csv(case_root/'mixed_mode_v2_final_summary.csv',index=False)
        return final

    jobs=[(k,psi) for k in items(a.classes) for psi in items(a.target_psi_deg,float)]
    results=[]
    with ThreadPoolExecutor(max_workers=max(1,a.max_jobs)) as ex:
        fut={ex.submit(control_case,k,psi):(k,psi) for k,psi in jobs}
        for f in as_completed(fut):
            k,psi=fut[f]
            try:z=f.result();z['status']='converged' if z['event_state_control_converged'] else 'not_converged'
            except Exception as exc:z={'class':k,'target_psi_deg':psi,'status':'failed','error':repr(exc)}
            results.append(z);print({q:z.get(q) for q in ['class','target_psi_deg','status','achieved_psi_deg','psi_error_deg','control_iterations','mode_classification']})
    pd.DataFrame(results).to_csv(outroot/'campaign_status_v2.csv',index=False)
    finals=[r for r in results if r.get('status')!='failed']
    if finals:pd.DataFrame(finals).to_csv(outroot/'mixed_mode_v2_all_final_cases.csv',index=False)
    if any(r.get('status')=='failed' for r in results):raise SystemExit('one or more v2 controlled cases failed')

if __name__=='__main__':main()
