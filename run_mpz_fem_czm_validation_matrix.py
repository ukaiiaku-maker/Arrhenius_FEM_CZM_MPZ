#!/usr/bin/env python3
"""Small full-physics FEM/CZM validation matrix for v9 moving-PZ parameters.

The solver remains the production anisotropic multifront code.  Single-front is
the default validation protocol only to isolate constitutive persistence; use
``--enable-branching`` to exercise branch birth/coalescence with the same model.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
from pathlib import Path
import subprocess
import sys

import pandas as pd

from mpz_run_utils import moving_pz_cli, fs, check_parameter_status


def build_cmd(py,row,klass,T,case,a):
    cmd=[py,'-m','arrhenius_fracture.sharp_front','--mode','2d',
         '--temperatures',str(int(T)),'--out',str(case),
         '--nx',str(a.nx),'--ny',str(a.ny),'--tip-h-fine',fs(a.tip_h_fine),
         '--tip-ratio',fs(a.tip_ratio),'--dU',fs(a.dU),'--dt',fs(a.dt),
         '--steps',str(a.steps),'--n-stagger',str(a.n_stagger),'--print-every',str(a.print_every),
         '--target-crack-extension-um',fs(a.target_ext_um),
         '--crystal-aniso','--crystal-compete','--crystal-material',a.crystal_material,
         '--crystal-theta-deg',fs(a.theta),'--cleave-gamma-aniso',fs(a.cleave_gamma_aniso),
         '--adaptive-events','--adaptive-event-target',fs(a.adaptive_event_target),
         '--adaptive-min-frac','1e-8','--adaptive-grow','4',
         '--da-phys',fs(a.da_phys),'--j-decomposition','cluster',
         '--rJ-cluster',fs(a.rJ_cluster),'--rJ-outer',fs(a.rJ_outer),
         '--crack-backend','adaptive_czm','--czm-max-angle-error-deg','35',
         '--save-snapshots',str(a.save_snapshots),'--snapshot-cols',str(a.snapshot_cols),
         '--snapshot-by-crack-extension-um',fs(a.snapshot_by_ext_um),
         '--coalesce-cracks'] + moving_pz_cli(row)
    if a.enable_branching:
        cmd += ['--crystal-branch','--max-fronts',str(a.max_fronts),
                '--branch-spacing',fs(a.branch_spacing),
                '--branch-fp-min-ratio',fs(a.branch_fp_min_ratio),
                '--branch-secondary-min-K-ratio',fs(a.branch_secondary_min_K_ratio)]
    else:
        cmd += ['--max-fronts','1']
    return cmd


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--parameters',default='mpz_three_class_initial_guesses.csv')
    ap.add_argument('--classes',default='ceramic weakT DBTT')
    ap.add_argument('--temperatures',default='300 700 1100')
    ap.add_argument('--out',default='runs/mpz_v9_fem_czm_validation')
    ap.add_argument('--python',default=sys.executable)
    ap.add_argument('--max-jobs',type=int,default=1)
    ap.add_argument('--dry-run',action='store_true')
    ap.add_argument('--skip-existing',action='store_true', help='Skip cases with an existing summary.json.')
    ap.add_argument('--require-fitted',action='store_true', help='Refuse parameter rows not marked FITTED_MPZ_V9.')
    ap.add_argument('--nx',type=int,default=60); ap.add_argument('--ny',type=int,default=120)
    ap.add_argument('--tip-h-fine',type=float,default=1e-6); ap.add_argument('--tip-ratio',type=float,default=1.25)
    ap.add_argument('--dU',type=float,default=2e-7); ap.add_argument('--dt',type=float,default=8.4)
    ap.add_argument('--steps',type=int,default=3000); ap.add_argument('--n-stagger',type=int,default=2)
    ap.add_argument('--print-every',type=int,default=10); ap.add_argument('--target-ext-um',type=float,default=100.0)
    ap.add_argument('--theta',type=float,default=45.0); ap.add_argument('--crystal-material',default='w')
    ap.add_argument('--cleave-gamma-aniso',type=float,default=0.3)
    ap.add_argument('--adaptive-event-target',type=float,default=0.25)
    ap.add_argument('--da-phys',type=float,default=5e-6); ap.add_argument('--rJ-cluster',type=float,default=20e-6)
    ap.add_argument('--rJ-outer',type=float,default=25e-6)
    ap.add_argument('--save-snapshots',type=int,default=8); ap.add_argument('--snapshot-cols',type=int,default=4)
    ap.add_argument('--snapshot-by-ext-um',type=float,default=20.0)
    ap.add_argument('--enable-branching',action='store_true')
    ap.add_argument('--max-fronts',type=int,default=16); ap.add_argument('--branch-spacing',type=float,default=10.0)
    ap.add_argument('--branch-fp-min-ratio',type=float,default=0.95)
    ap.add_argument('--branch-secondary-min-K-ratio',type=float,default=0.85)
    a=ap.parse_args()
    root=Path(a.out); root.mkdir(parents=True,exist_ok=True)
    params=pd.read_csv(a.parameters).set_index('target_class',drop=False)
    check_parameter_status(params, a.parameters, a.require_fitted)
    jobs=[]
    for klass in a.classes.replace(',',' ').split():
      row=params.loc[klass]
      for T in [float(x) for x in a.temperatures.replace(',',' ').split()]:
        case=root/klass/f'T{int(T)}K'
        cmd=build_cmd(a.python,row,klass,T,case,a)
        skipped=bool(a.skip_existing and (case/'summary.json').exists())
        jobs.append({'class':klass,'T_K':T,'case':str(case),'cmd':cmd,'skipped_existing':skipped})
        print(' '.join(cmd))
    (root/'run_manifest.json').write_text(json.dumps(jobs,indent=2))
    if a.dry_run: return
    jobs_to_run=[j for j in jobs if not j.get('skipped_existing',False)]
    for j in jobs:
        if j.get('skipped_existing',False):
            print('SKIP existing',j['class'],j['T_K'],j['case'])
    def run(j):
        Path(j['case']).mkdir(parents=True,exist_ok=True)
        log=Path(j['case'])/'run.log'
        with log.open('w') as fp:
            cp=subprocess.run(j['cmd'],stdout=fp,stderr=subprocess.STDOUT,text=True)
        return j,cp.returncode
    with ThreadPoolExecutor(max_workers=max(a.max_jobs,1)) as ex:
      futs=[ex.submit(run,j) for j in jobs_to_run]
      for f in as_completed(futs):
        j,rc=f.result(); print(j['class'],j['T_K'],'returncode',rc)
        if rc: raise SystemExit(f"failed: {j['case']}")

if __name__=='__main__': main()
