#!/usr/bin/env python3
"""Small V1-only cleavage EXP-floor shape sweep.

This is intended as the first-pass tuning/classification tool before expensive
v8 2-D validation.  It sweeps the cleavage free-energy magnitude/shape while
using the same K-controlled V1 fatigue controller and writes one row per
parameter set with a coarse regime label.
"""
from __future__ import annotations
import argparse, itertools, json, shutil, subprocess, sys
from pathlib import Path
import pandas as pd


def run(cmd, log):
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open('w') as f:
        f.write('$ '+' '.join(cmd)+'\n\n'); f.flush()
        p = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)
    if p.returncode:
        raise RuntimeError(f'command failed: {cmd}; see {log}')


def classify(df, cycles_max):
    fired = df['cycles_to_first_fire'].notna()
    if not fired.any():
        return 'inactive_or_threshold_above_range'
    if fired.all():
        # if every K fires and the lowest-K life is far below the horizon, it is
        # more direct-fracture-like than endurance-like for this K window.
        nlo = float(df.sort_values('Kmax_MPa_sqrtm').iloc[0]['cycles_to_first_fire'])
        if nlo < 0.1 * cycles_max:
            return 'direct_fracture_no_limit_in_window'
        return 'all_fire_near_horizon'
    # Some fire and some do not: a threshold/window exists in this K range.
    return 'apparent_threshold_in_window'


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', default='runs/v1_cleave_exp_floor_sweep')
    ap.add_argument('--Kmax-MPa-sqrt-m', nargs='+', type=float, default=[8,9,9.5,9.75,10,10.25,10.5,11])
    ap.add_argument('--T', type=float, default=300.0)
    ap.add_argument('--R', type=float, default=0.1)
    ap.add_argument('--frequency-Hz', type=float, default=1000.0)
    ap.add_argument('--cycles-max', type=float, default=1e11)
    ap.add_argument('--max-blocks', type=int, default=160)
    ap.add_argument('--target-dB', type=float, default=0.02)
    ap.add_argument('--target-dN-store', type=float, default=0.025)
    ap.add_argument('--target-dN-emit', type=float, default=0.25)
    ap.add_argument('--target-dN-mobile', type=float, default=0.25)
    ap.add_argument('--storage-model', default='fixed_fraction', choices=['fixed_fraction','all_retained','escape_limited'])
    ap.add_argument('--fixed-retained-fraction', type=float, default=0.1)
    ap.add_argument('--cleave-G00-eV', nargs='+', type=float, default=[2.0,2.4,2.8])
    ap.add_argument('--cleave-sigc0-GPa', nargs='+', type=float, default=[3.5,5.0,7.0])
    ap.add_argument('--cleave-exp-a', nargs='+', type=float, default=[0.25,0.5,1.0])
    ap.add_argument('--cleave-exp-n', nargs='+', type=float, default=[1.5,2.0,3.0])
    ap.add_argument('--cleave-floor-frac', nargs='+', type=float, default=[0.02,0.05,0.1])
    ap.add_argument('--cleave-exp-T-mode', default='mu_scale', choices=['linear','mu_scale'])
    ap.add_argument('--keep-existing', action='store_true')
    args = ap.parse_args()
    out = Path(args.out)
    if out.exists() and not args.keep_existing:
        shutil.rmtree(out)
    out.mkdir(parents=True, exist_ok=True)
    rows=[]
    combos=list(itertools.product(args.cleave_G00_eV,args.cleave_sigc0_GPa,args.cleave_exp_a,args.cleave_exp_n,args.cleave_floor_frac))
    for i,(G00,sigc,a,n,ff) in enumerate(combos):
        case=out/f'case_{i:04d}_G{G00:g}_sc{sigc:g}_a{a:g}_n{n:g}_ff{ff:g}'.replace('.','p')
        K_rows=[]
        for K in args.Kmax_MPa_sqrt_m:
            kout=case/f'K{K:g}'.replace('.','p')
            cmd=[sys.executable,'-m','arrhenius_fracture.fatigue_sharp_front',
                 '--temperatures',str(args.T),'--Kmax-MPa-sqrt-m',str(K),
                 '--R',str(args.R),'--frequency-Hz',str(args.frequency_Hz),
                 '--cycles-max',str(args.cycles_max),'--max-blocks',str(args.max_blocks),
                 '--block-cycles','1e5','--max-block-cycles','inf','--cycle-block-mode','hazard_limited',
                 '--target-dB',str(args.target_dB),'--target-dN-store',str(args.target_dN_store),
                 '--target-dN-emit',str(args.target_dN_emit),'--target-dN-mobile',str(args.target_dN_mobile),
                 '--storage-model',args.storage_model,'--dN-cap','inf','--sigma-cap-GPa','0','--no-plots',
                 '--cleave-barrier-kind','exp_floor','--cleave-exp-T-mode',args.cleave_exp_T_mode,
                 '--cleave-G00-eV',str(G00),'--cleave-sigc0-GPa',str(sigc),
                 '--cleave-exp-a',str(a),'--cleave-exp-n',str(n),'--cleave-floor-frac',str(ff),
                 '--out',str(kout)]
            if args.storage_model=='fixed_fraction':
                cmd += ['--fixed-retained-fraction', str(args.fixed_retained_fraction)]
            run(cmd, case/f'K{K:g}.log'.replace('.','p'))
            hist=pd.read_csv(kout/f'T{int(args.T)}K'/'fatigue_v1_history.csv')
            fired=hist[hist['n_fire']>0]
            first=float(fired.iloc[0]['cycles_total']) if len(fired) else None
            last=hist.iloc[-1]
            K_rows.append({'Kmax_MPa_sqrtm':K,'cycles_to_first_fire':first,'cycles_total':float(last['cycles_total']),
                           'B_final':float(last['B']),'mu_cleave':float(last['mu_cleave_pred']),
                           'S_cleave_kB':float(last.get('S_cleave_kB', float('nan'))),
                           'dGds_eV_per_GPa':float(last.get('dGcleave_dsigma_eV_per_GPa', float('nan')))})
        kdf=pd.DataFrame(K_rows); kdf.to_csv(case/'K_summary.csv',index=False)
        reg=classify(kdf,args.cycles_max)
        rows.append({'case':i,'regime':reg,'G00_eV':G00,'sigc0_GPa':sigc,'a':a,'n':n,'floor_frac':ff,
                     'n_fire_conditions':int(kdf['cycles_to_first_fire'].notna().sum()),
                     'min_fire_cycles':kdf['cycles_to_first_fire'].dropna().min() if kdf['cycles_to_first_fire'].notna().any() else None,
                     'case_dir':str(case)})
        pd.DataFrame(rows).to_csv(out/'sweep_summary.csv',index=False)
        print(f'[{i+1}/{len(combos)}] {reg}: G00={G00:g} sigc={sigc:g} a={a:g} n={n:g} floor={ff:g}')
    with (out/'sweep_settings.json').open('w') as f: json.dump(vars(args),f,indent=2)
    print(f'Wrote {out/"sweep_summary.csv"}')

if __name__=='__main__': main()
