#!/usr/bin/env python3
"""Re-parameterize four response classes with the v9 moving process zone.

The objective contains both initiation and repeated-growth terms.  No fitted
``chi_shield``, scalar ``N_sat``, emission-per-step cap, or stored-energy
cleavage offset exists in this model.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
from scipy.optimize import differential_evolution

from arrhenius_fracture import sharp_front as sf
from arrhenius_fracture.config import make_emergent_config


def fs(x):
    return f"{float(x):.16g}"


def row_to_args(row, dK, Kdot, n_advances, Kmax):
    cli = [
        '--mode','1d','--front-state-model','moving_pz','--temperatures','300',
        '--Kdot',fs(Kdot),'--Kmax',fs(Kmax),'--dt',fs(dK/Kdot),
        '--continue-after-init','--n-advances',str(n_advances),'--sigma-cap-GPa','0',
        '--r-pz',fs(row.r_pz_m),'--c-blunt',fs(row.c_blunt),
        '--L-pz',fs(row.mpz_length_m),'--mpz-length-m',fs(row.mpz_length_m),
        '--mpz-n-bins',str(int(row.mpz_n_bins)),'--mpz-n-systems',str(int(row.mpz_n_systems)),
        '--mpz-source-sites-per-system',fs(row.mpz_source_sites_per_system),
        '--mpz-source-recovery-rate-s',fs(row.mpz_source_recovery_rate_s),
        '--mpz-source-refresh-length-m',fs(row.mpz_source_refresh_length_m),
        '--mpz-shielding-factors',str(row.mpz_shielding_factors),
        '--mpz-glide-barrier-eV',fs(row.mpz_glide_barrier_eV),
        '--mpz-glide-activation-volume-b3',fs(row.mpz_glide_activation_volume_b3),
        '--mpz-trap-barrier-eV',fs(row.mpz_trap_barrier_eV),
        '--mpz-detrap-barrier-eV',fs(row.mpz_detrap_barrier_eV),
        '--mpz-retained-recovery-barrier-eV',fs(row.mpz_retained_recovery_barrier_eV),
        '--mpz-pair-annihilation-rate-per-count-s',fs(row.mpz_pair_annihilation_rate_per_count_s),
        '--cleave-barrier-kind','exp_floor','--cleave-exp-T-mode','linear',
        '--cleave-G00-eV',fs(row.cleave_G00_eV),
        '--cleave-gT-eV-per-K',fs(row.cleave_gT_eV_per_K),
        '--cleave-sigc0-GPa',fs(row.cleave_sigc0_GPa),
        '--cleave-sT-GPa-per-K',fs(row.cleave_sT_GPa_per_K),
        '--cleave-exp-a',fs(row.cleave_exp_a),'--cleave-exp-n',fs(row.cleave_exp_n),
        '--cleave-floor-frac',fs(row.cleave_floor_frac),
        '--cleave-S-hs-kB',fs(row.cleave_S_hs_kB),
        '--cleave-sigma-S-GPa','6','--cleave-S-hs-power','2',
        '--emit-barrier-kind','exp_floor',
        '--emit-G00-eV',fs(row.emit_G00_eV),
        '--emit-gT-eV-per-K',fs(row.emit_gT_eV_per_K),
        '--emit-sigc0-GPa',fs(row.emit_sigc0_GPa),
        '--emit-sT-GPa-per-K',fs(row.emit_sT_GPa_per_K),
        '--emit-exp-a',fs(row.emit_exp_a),'--emit-exp-n',fs(row.emit_exp_n),
        '--emit-floor-frac',fs(row.emit_floor_frac),
    ]
    return sf._build_parser().parse_args(cli)


def simulate(row, T, dK, Kdot, n_advances, Kmax):
    args = row_to_args(row, dK, Kdot, n_advances, Kmax)
    eng = sf.build_engine(args, make_emergent_config().material)
    dt = dK / Kdot
    events = []
    max_shield_fraction = 0.0
    event_shield_fractions = []
    min_site_fraction = 1.0
    max_retained = 0.0
    for i in range(int(math.ceil(Kmax / dK))):
        Kmpa = (i + 1) * dK
        info = eng.step(Kmpa * 1e6, T, dt)
        Ksh = float(info.get('mpz_K_shield_Pa_sqrt_m', 0.0)) / 1e6
        min_site_fraction = min(min_site_fraction, float(info.get('mpz_available_site_fraction', 1.0)))
        max_retained = max(max_retained, float(info.get('mpz_retained_count', info.get('N_em', 0.0))))
        nf = int(info.get('n_fire', 0))
        if nf:
            frac = Ksh / max(Kmpa, 1e-12)
            event_shield_fractions.extend([frac] * nf)
            max_shield_fraction = max(max_shield_fraction, frac)
            events.extend([Kmpa] * nf)
        if len(events) >= n_advances:
            break
    ev = np.asarray(events[:n_advances], dtype=float)
    return {
        'K_init': float(ev[0]) if len(ev) else np.nan,
        'K_late': float(np.median(ev[-min(5, len(ev)):])) if len(ev) else np.nan,
        'delta_KR': float(np.median(ev[-min(5, len(ev)):]) - ev[0]) if len(ev) else np.nan,
        'n_events': int(len(ev)),
        'max_shield_fraction': max_shield_fraction,
        'min_site_fraction': min_site_fraction,
        'max_retained': max_retained,
        'events': ev.tolist(),
    }


def vector_to_row(base, x):
    r = base.copy()
    (r.cleave_G00_eV, r.cleave_gT_eV_per_K,
     r.emit_G00_eV, r.emit_gT_eV_per_K,
     log_sites, r.mpz_glide_barrier_eV, r.mpz_trap_barrier_eV,
     r.mpz_detrap_barrier_eV, r.mpz_retained_recovery_barrier_eV) = x
    r.mpz_source_sites_per_system = math.exp(log_sites)
    return r


def bounds_for(base):
    cg = float(base.cleave_G00_eV); eg = float(base.emit_G00_eV)
    return [
        (max(0.3, 0.55*cg), 1.6*cg),
        (-0.006, 0.012),
        (max(0.25, 0.45*eg), 1.8*eg),
        (-0.006, 0.014),
        (math.log(5.0), math.log(2000.0)),
        (0.25, 2.5), (0.15, 2.5), (0.35, 3.0), (0.6, 3.5),
    ]


def objective_for_class(base, targets, klass, opt):
    Ts = targets.T_K.to_numpy(float)
    Ktar = targets.K_init_target_MPa_sqrt_m.to_numpy(float)

    def obj(x, details=False):
        row = vector_to_row(base, x)
        sims = [simulate(row, T, opt.dK, opt.Kdot, opt.n_advances, opt.Kmax) for T in Ts]
        Ki = np.array([q['K_init'] for q in sims], float)
        Kl = np.array([q['K_late'] for q in sims], float)
        valid = np.isfinite(Ktar)
        observed = valid & np.isfinite(Ki)
        miss = int(np.count_nonzero(valid & ~np.isfinite(Ki)))
        if np.any(observed):
            first = float(np.sqrt(np.mean(((Ki[observed] - Ktar[observed]) /
                                           max(opt.K_scale, 1e-9))**2)))
        else:
            first = 0.0
        if miss:
            first += 20.0 * miss

        finite = np.isfinite(Ki) & np.isfinite(Kl)
        persistence = 0.0
        if finite.sum() >= 3:
            a = Ki[finite]; b = Kl[finite]
            an = (a-a.mean()) / max(a.std(), 1e-9)
            bn = (b-b.mean()) / max(b.std(), 1e-9)
            persistence += float(np.sqrt(np.mean((an-bn)**2)))
            if klass == 'ceramic':
                persistence += max(float(np.nanmax(Kl-Ki)) - 2.0, 0.0) / 2.0
                persistence += max(max(q['max_shield_fraction'] for q in sims) - 0.15, 0.0) * 5.0
            elif klass == 'weakT':
                persistence += np.nanstd(Kl) / 2.0
            elif klass == 'DBTT':
                n = max(1, len(a)//3)
                ci = float(np.nanmean(a[-n:]) - np.nanmean(a[:n]))
                cl = float(np.nanmean(b[-n:]) - np.nanmean(b[:n]))
                persistence += max(8.0-ci, 0.0)/4.0 + max(0.7*ci-cl, 0.0)/4.0
            elif klass == 'peak':
                if len(a) >= 5:
                    pi = float(np.nanmax(a[1:-1]) - max(a[0], a[-1]))
                    pl = float(np.nanmax(b[1:-1]) - max(b[0], b[-1]))
                    persistence += max(1.0-pi, 0.0) + max(0.5*pi-pl, 0.0)

        state_pen = 0.0
        for q in sims:
            state_pen += max(q['max_shield_fraction'] - 0.95, 0.0) * 20.0
            state_pen += 2.0 if q['min_site_fraction'] < 1e-6 else 0.0
            state_pen += 5.0 if q['n_events'] == 0 else 0.0
        score = opt.w_first*first + opt.w_persistence*persistence + opt.w_state*state_pen
        if details:
            return score, row, sims, {'first': first, 'persistence': persistence, 'state': state_pen}
        return score
    return obj


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--initial', default='mpz_four_class_initial_guesses.csv')
    ap.add_argument('--targets', default='legacy_first_passage_targets_generated.csv')
    ap.add_argument('--out', default='runs/mpz_v9_four_class_fit')
    ap.add_argument('--classes', default='ceramic peak weakT DBTT')
    ap.add_argument('--temperatures', default='300 500 700 900 1100')
    ap.add_argument('--dK', type=float, default=0.10)
    ap.add_argument('--Kdot', type=float, default=0.005)
    ap.add_argument('--Kmax', type=float, default=65.0)
    ap.add_argument('--n-advances', type=int, default=20)
    ap.add_argument('--K-scale', type=float, default=2.0)
    ap.add_argument('--w-first', type=float, default=1.0)
    ap.add_argument('--w-persistence', type=float, default=2.0)
    ap.add_argument('--w-state', type=float, default=1.0)
    ap.add_argument('--popsize', type=int, default=12)
    ap.add_argument('--maxiter', type=int, default=40)
    ap.add_argument('--seed', type=int, default=9107)
    ap.add_argument('--smoke', action='store_true', help='evaluate initial guesses without optimization')
    a = ap.parse_args()
    out=Path(a.out); out.mkdir(parents=True, exist_ok=True)
    init=pd.read_csv(a.initial).set_index('target_class', drop=False)
    targets=pd.read_csv(a.targets)
    wantedT={float(x) for x in a.temperatures.replace(',',' ').split()}
    classes=a.classes.replace(',',' ').split()
    final=[]; predictions=[]
    for klass in classes:
        base=init.loc[klass].copy()
        tg=targets[(targets.target_class.astype(str)==klass) & targets.T_K.astype(float).isin(wantedT)].sort_values('T_K')
        if tg.empty: raise SystemExit(f'no targets for {klass}')
        fun=objective_for_class(base,tg,klass,a)
        x0=np.array([base.cleave_G00_eV,base.cleave_gT_eV_per_K,base.emit_G00_eV,
                     base.emit_gT_eV_per_K,math.log(base.mpz_source_sites_per_system),
                     base.mpz_glide_barrier_eV,base.mpz_trap_barrier_eV,
                     base.mpz_detrap_barrier_eV,base.mpz_retained_recovery_barrier_eV],float)
        if a.smoke:
            xbest=x0
        else:
            res=differential_evolution(fun,bounds_for(base),seed=a.seed,popsize=a.popsize,
                                       maxiter=a.maxiter,polish=True,workers=1,updating='immediate')
            xbest=res.x
        score,row,sims,parts=fun(xbest,details=True)
        row['fit_score']=score; row['fit_first_component']=parts['first']; row['fit_persistence_component']=parts['persistence']; row['fit_state_component']=parts['state']; row['status']='FITTED_MPZ_V9' if not a.smoke else 'INITIAL_SMOKE_ONLY'
        final.append(dict(row))
        for T,q in zip(tg.T_K,sims):
            predictions.append({'target_class':klass,'T_K':T,**q})
        print(klass,'score',score,parts)
    pd.DataFrame(final).to_csv(out/'mpz_four_class_parameters.csv',index=False)
    pd.DataFrame(predictions).to_json(out/'mpz_fit_predictions.json',orient='records',indent=2)
    (out/'run_config.json').write_text(json.dumps(vars(a),indent=2))


if __name__=='__main__':
    main()
