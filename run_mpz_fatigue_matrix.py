#!/usr/bin/env python3
"""Run cyclic loading with the same v9 material parameters used in monotonic/dwell tests."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys

import pandas as pd

from mpz_run_utils import moving_pz_cli, fs, check_parameter_status


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--parameters',default='mpz_four_class_initial_guesses.csv')
    ap.add_argument('--classes',default='ceramic peak weakT DBTT')
    ap.add_argument('--temperatures',default='300 700 1100')
    ap.add_argument('--Kmax-values',default='6 10 14')
    ap.add_argument('--R',type=float,default=0.1)
    ap.add_argument('--frequency-Hz',type=float,default=1000.0)
    ap.add_argument('--cycles-max',type=float,default=1e10)
    ap.add_argument('--block-cycles',type=float,default=1e5)
    ap.add_argument('--max-blocks',type=int,default=10000)
    ap.add_argument('--n-advances',type=int,default=20)
    ap.add_argument('--out',default='runs/mpz_v9_fatigue_matrix')
    ap.add_argument('--python',default=sys.executable)
    ap.add_argument('--dry-run',action='store_true')
    ap.add_argument('--skip-existing',action='store_true', help='Skip cases with an existing summary_all.json.')
    ap.add_argument('--require-fitted',action='store_true', help='Refuse parameter rows not marked FITTED_MPZ_V9.')
    a=ap.parse_args()
    params=pd.read_csv(a.parameters).set_index('target_class',drop=False)
    check_parameter_status(params, a.parameters, a.require_fitted)
    root=Path(a.out); root.mkdir(parents=True,exist_ok=True)
    jobs=[]
    for klass in a.classes.replace(',',' ').split():
      row=params.loc[klass]
      for T in [float(x) for x in a.temperatures.replace(',',' ').split()]:
       for K in [float(x) for x in a.Kmax_values.replace(',',' ').split()]:
        case=root/klass/f'T{int(T)}K'/f'K{K:g}'
        cmd=[a.python,'-m','arrhenius_fracture.fatigue_sharp_front',
             '--temperatures',fs(T),'--Kmax-MPa-sqrt-m',fs(K),'--R',fs(a.R),
             '--frequency-Hz',fs(a.frequency_Hz),'--cycles-max',fs(a.cycles_max),
             '--cycle-block-mode','hazard_limited','--block-cycles',fs(a.block_cycles),
             '--max-block-cycles',fs(a.block_cycles),'--max-blocks',str(a.max_blocks),
             '--n-advances',str(a.n_advances),'--continue-after-fire','--no-plots',
             '--out',str(case)] + moving_pz_cli(row)
        skipped=bool(a.skip_existing and (case/'summary_all.json').exists())
        jobs.append({'class':klass,'T_K':T,'Kmax':K,'case':str(case),'cmd':cmd,'skipped_existing':skipped})
        print(' '.join(cmd))
        if skipped:
            print('SKIP existing',case)
        if not a.dry_run and not skipped:
            case.mkdir(parents=True,exist_ok=True)
            cp=subprocess.run(cmd,text=True)
            if cp.returncode: raise SystemExit(cp.returncode)
    (root/'run_manifest.json').write_text(json.dumps(jobs,indent=2))

if __name__=='__main__': main()
