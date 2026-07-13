#!/usr/bin/env python3
"""Constant-K dwell/creep-fracture driver for the moving process-zone engine."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from arrhenius_fracture import sharp_front as sf
from arrhenius_fracture.config import make_emergent_config
from fit_mpz_three_classes import row_to_args
from mpz_run_utils import check_parameter_status


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--parameters',default='mpz_three_class_initial_guesses.csv')
    ap.add_argument('--classes',default='ceramic weakT DBTT')
    ap.add_argument('--temperatures',default='300 700 1100')
    ap.add_argument('--K-MPa-sqrt-m',type=float,default=15.0,dest='K_MPa_sqrt_m')
    ap.add_argument('--hold-s',type=float,default=1e5)
    ap.add_argument('--dt-initial-s',type=float,default=1.0)
    ap.add_argument('--dt-max-s',type=float,default=1e4)
    ap.add_argument('--target-dB',type=float,default=0.1)
    ap.add_argument('--target-emission-hazard',type=float,default=0.1,
                    help='Maximum fraction-like finite-site emission hazard while a non-negligible source inventory remains.')
    ap.add_argument('--source-active-fraction-min',type=float,default=1e-6,
                    help='Stop emission-driven substepping after the finite source inventory falls below this fraction.')
    ap.add_argument('--max-blocks',type=int,default=100000,
                    help='Hard safety bound preventing a stalled dwell integration.')
    ap.add_argument('--n-advances',type=int,default=10)
    ap.add_argument('--out',default='runs/mpz_v9_dwell')
    ap.add_argument('--skip-existing',action='store_true', help='Skip cases with an existing summary.json.')
    ap.add_argument('--require-fitted',action='store_true', help='Refuse parameter rows not marked FITTED_MPZ_V9.')
    a=ap.parse_args()
    out=Path(a.out); out.mkdir(parents=True,exist_ok=True)
    params=pd.read_csv(a.parameters).set_index('target_class',drop=False)
    check_parameter_status(params, a.parameters, a.require_fitted)
    rows=[]; summaries=[]
    for klass in a.classes.replace(',',' ').split():
      r=params.loc[klass]
      for T in [float(x) for x in a.temperatures.replace(',',' ').split()]:
        cdir=out/klass/f'T{int(T)}K'
        if a.skip_existing and (cdir/'summary.json').exists():
          print('SKIP existing',cdir)
          summaries.append(json.loads((cdir/'summary.json').read_text()))
          continue
        args=row_to_args(r,dK=0.1,Kdot=0.005,n_advances=a.n_advances,Kmax=100.0,da_um=5.0)
        eng=sf.build_engine(args,make_emergent_config().material)
        t=0.0; dt=a.dt_initial_s; ib=0
        case=[]
        while t<a.hold_s and eng.n_adv<a.n_advances:
          if ib >= a.max_blocks:
            raise RuntimeError(f'dwell integration exceeded --max-blocks={a.max_blocks} for {klass} at {T:g} K')
          sig=eng.sigma_tip(a.K_MPa_sqrt_m*1e6)
          le=eng.lambda_emit(sig,T)[0]
          db_rate=eng.lambda_cleave(sig,T)[0]
          limits=[dt,a.dt_max_s,a.hold_s-t]
          if db_rate>0:
            limits.append(a.target_dB/db_rate)
          # Finite-site emission is integrated exactly as 1-exp(-H).  Limit the
          # depletion fraction only while a meaningful source inventory remains;
          # unlike the legacy lambda*dt cap, this cannot trap the controller in
          # vanishing timesteps after the sites have been exhausted.
          cap=max(float(np.sum(eng.mpz_state.site_capacity)),1e-300)
          avail_frac=float(np.sum(eng.mpz_state.available_sites))/cap
          if (le>0 and a.target_emission_hazard>0
                  and avail_frac>a.source_active_fraction_min):
            frac=min(max(a.target_emission_hazard,1e-12),1.0-1e-12)
            limits.append(-np.log1p(-frac)/le)
          dt=max(min(limits),1e-12)
          info=eng.step(a.K_MPa_sqrt_m*1e6,T,dt)
          t+=dt; ib+=1
          rec={'class':klass,'T_K':T,'block':ib,'time_s':t,'dt_s':dt,
               'K_MPa_sqrt_m':a.K_MPa_sqrt_m,'n_adv':eng.n_adv,'a_adv_um':eng.a_adv*1e6,
               **{k:v for k,v in info.items() if isinstance(v,(int,float,bool,np.integer,np.floating))}}
          case.append(rec); rows.append(rec)
          if info.get('fired'): dt=max(a.dt_initial_s,dt*0.25)
          else: dt=min(a.dt_max_s,dt*1.5)
        cdir.mkdir(parents=True,exist_ok=True)
        pd.DataFrame(case).to_csv(cdir/'dwell_history.csv',index=False)
        sm={'class':klass,'T_K':T,'K_MPa_sqrt_m':a.K_MPa_sqrt_m,'time_s':t,
            'n_adv':int(eng.n_adv),'a_adv_um':eng.a_adv*1e6,'failed':bool(eng.n_adv>0)}
        summaries.append(sm); (cdir/'summary.json').write_text(json.dumps(sm,indent=2))
        print(sm)
    pd.DataFrame(summaries).to_csv(out/'dwell_summary.csv',index=False)
    (out/'run_config.json').write_text(json.dumps(vars(a),indent=2))

if __name__=='__main__': main()
