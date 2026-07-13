"""Entropy form & magnitude sweep for the sharp-front fracture model.

Motivation
----------
The activation entropy S*(sigma,T) of the EMISSION channel sets the slope of the
barrier with temperature and therefore controls the ductile-to-brittle
transition.  Its functional FORM and its MAGNITUDE are both uncertain, and they
matter very differently for the two channels:

  * crack-tip / cleavage:  small S*, a few k_B.  An entropy-enthalpy
    compensation (Meyer-Neldel / Schoeck) form gives ~3-4 k_B and is bounded.

  * dislocation emission / nucleation:  LARGE S*.  Nanopillar fits across
    several systems (G0(T) linear, S* = -dG0/dT) give a remarkably universal
    -36 to -46 k_B, e.g. W[100] = -45.7 k_B with G0(481 K)=1.94 eV.

This driver sweeps the emission-entropy magnitude across that whole range and,
separately, a few candidate functional forms, and reports how each reshapes the
toughness-temperature curve Kc(T) and the DBTT.  It runs on the 1D K-ramp
(single advance law, no FEM): fast, and the cleavage clock is dt-independent so
the result is trustworthy at coarse dt.

Experimental emission-entropy anchors (S* = -dG0/dT from nanopillar fits):
    system              S*/k_B   G0(Tref) eV   sigc(Tref) GPa   Tref K
    Al0.7CoCrFeNi-BCC   -40.8    0.854         2.847            159
    Al0.7CoCrFeNi-FCC   -46.4    0.684         1.662            159
    Cu                  -38.5    2.114         0.549            481
    Si                  -45.6    2.471         6.017            527
    Ta[111]             -35.6    2.116         0.713            381
    W[100]              -45.7    1.940         2.298            481

Usage
-----
    python -m arrhenius_fracture.entropy_form_sweep --mode magnitude
    python -m arrhenius_fracture.entropy_form_sweep --mode form
    python -m arrhenius_fracture.entropy_form_sweep --mode both --out runs/Ssweep
"""
import argparse
import csv
import json
import os

import numpy as np

from arrhenius_fracture.sharp_front import (_build_parser, build_engine,
                                            make_emergent_config)

EV = 1.602176634e-19
KB_EV = 8.617333262e-5  # eV/K

# --- W[100] pivot anchor (so all magnitude curves share G0(T_ref)) -----------
W_G0_REF_EV = 1.94
W_T_REF_K = 481.33
W_SIGC_GPA = 2.298

# experimental emission-entropy anchors (kB) for marking on plots
EXP_ANCHORS = {
    'Al0.7CoCrFeNi-BCC': -40.8, 'Al0.7CoCrFeNi-FCC': -46.4, 'Cu': -38.5,
    'Si': -45.6, 'Ta[111]': -35.6, 'W[100]': -45.7,
}


# ----------------------------------------------------------------------------
# Kc(T) for one emission-barrier configuration (reuses the engine machinery)
# ----------------------------------------------------------------------------
def _kc_vs_T(configure, temps, base_args, configure_cb=None):
    """Run the 1D K-ramp at each T after `configure(eb)` sets the emission
    barrier (and optional `configure_cb(cb)` sets cleavage).  Returns list of
    per-T dicts (Kc, sigma_back, r_eff, mode)."""
    mat = make_emergent_config().material
    Kdot = base_args.Kdot * 1e6
    dt = base_args.dt
    t_max = base_args.Kmax * 1e6 / Kdot
    nstep = int(np.ceil(t_max / dt))
    out = []
    for T in temps:
        eng = build_engine(base_args, mat)
        configure(eng.eb)
        eng.eb.__dict__.pop('_argmin_cache', None)   # form changed: drop cache
        if configure_cb is not None:
            configure_cb(eng.cb)
            eng.cb.__dict__.pop('_argmin_cache', None)
        Kc = None
        info = {}
        for i in range(nstep):
            K = Kdot * (i + 1) * dt
            info = eng.step(K, T, dt)
            if info['fired'] and Kc is None:
                Kc = K
                break
        sig_tip = max(info.get('sigma_tip', 0.0), 1.0)
        shield = info.get('sigma_back', 0.0) / sig_tip
        blunt = info.get('r_eff', eng.f.r0) / eng.f.r0
        ductile = (shield > base_args.ductile_shield) or (blunt > base_args.ductile_blunt)
        out.append({
            'T': float(T),
            'Kc_MPa_sqrt_m': None if Kc is None else Kc / 1e6,
            'sigma_back_GPa': info.get('sigma_back', 0.0) / 1e9,
            'r_eff_over_r0': blunt, 'N_em': info.get('N_em', 0.0),
            'mode': 'no-fracture' if Kc is None else ('ductile' if ductile else 'brittle'),
        })
    return out


def _dbtt(rows):
    """Onset T of the contiguous high-T ductile block (None if no block)."""
    rows = sorted(rows, key=lambda r: r['T'])
    dbtt = None
    for r in rows:
        if r['mode'] == 'ductile':
            if dbtt is None:
                dbtt = r['T']
        else:
            dbtt = None  # block broken; reset (want the FINAL contiguous block)
    # re-scan for the *lowest* T whose tail is all ductile
    dbtt = None
    for i in range(len(rows)):
        if all(r['mode'] == 'ductile' for r in rows[i:]):
            dbtt = rows[i]['T']
            break
    return dbtt


# ----------------------------------------------------------------------------
# emission-barrier configurators
# ----------------------------------------------------------------------------
# exact EXP_floor parameters from the BarrierModel export (nanopillar fits).
# system: (G00_eV, gT_eV/K, sigc0_Pa, sT_Pa/K, Tref_K, a, n, Tmin_K, Tmax_K)
EXP_FLOOR_SYSTEMS = {
    'W[100]':            (1.94022, 0.00393367, 2.29797e9, -656405,  481.33, 0.0845685, 1.0,      298, 673),
    'Ta[111]':           (2.11553, 0.003072,   7.13105e8, -738603,  381.33, 0.367026,  0.1,      298, 473),
    'Cu':                (2.11351, 0.00331867, 5.4872e8,  -922577,  481.33, 0.228066,  0.100063, 298, 673),
    'Si':                (2.47105, 0.00392574, 6.01655e9, -1.08048e7, 527.17, 0.189927, 1.0,     298, 773),
    'Al0.7CoCrFeNi-BCC': (0.854276, 0.00351484, 2.84708e9, -2.7131e6, 159.33, 0.654958, 0.1,     40, 295),
    'Al0.7CoCrFeNi-FCC': (0.683717, 0.00399529, 1.66175e9, -2.33278e6, 159.33, 0.130071, 0.357397, 40, 295),
}


def _cfg_exp_floor_emit(system='W[100]', a=None, n=None, gt_sign=1.0,
                        S_hs_kB=0.0, sigma_S_GPa=6.0, S_hs_power=3.0,
                        S_hs_dT=0.0):
    """Emission barrier = experimental EXP_floor model for `system`, using the
    exact fitted (G00,gT,sigc0,sT,Tref,a,n) from the BarrierModel export.
    Pass a/n to override the fitted exponents (a,n sensitivity sweep).

    High-stress entropy crossover (fatigue-paper hypothesis): S_hs_kB shifts S*
    toward less-negative/positive above sigma_S_GPa (preserving the fit at low
    stress and at Tref).  S_hs_kB=0 is the pure data fit; S_hs_kB > ~45 makes
    S* positive at the crack tip, letting emission ramp with T -> a DBTT.

    gt_sign multiplies the fitted gT (=+1 literal fit, S*<0)."""
    G00, gT, sigc0, sT, Tref, a_fit, n_fit, Tmin, Tmax = EXP_FLOOR_SYSTEMS[system]
    a_use = a_fit if a is None else a
    n_use = n_fit if n is None else n

    def configure(eb):
        eb.barrier_kind = 'exp_floor'
        eb.ef_G00_eV, eb.ef_gT_eV_per_K = G00, gt_sign * gT
        eb.ef_sigc0_Pa, eb.ef_sT_Pa_per_K = sigc0, sT
        eb.ef_Tref_K = Tref
        eb.ef_a, eb.ef_n = a_use, n_use
        eb.ef_S_hs_kB = S_hs_kB
        eb.ef_sigma_S_GPa = sigma_S_GPa
        eb.ef_S_hs_power = S_hs_power
        eb.ef_S_hs_dT_per_K = S_hs_dT
        eb.ef_S_hs_Tref_K = Tref
    return configure


def _cfg_cleave_modulus(S_kB_scale=3.6):
    """Crack-opening (cleavage) barrier: WEAK entropy at the modulus-softening
    scale (Meyer-Neldel), per the instruction that crack-opening entropy is
    ~ the temperature dependence of the modulus, not the large emission value.
    T_MN chosen so S*(sigma=0) ~ S_kB_scale at H0~2 eV."""
    def configure(cb):
        cb.barrier_kind = 'classic'
        cb.use_negative_entropy = True
        cb.entropy_stress_form = 'meyer_neldel'
        # S*(0) = H0/(kB*T_MN); pick T_MN to hit the requested few-kB scale
        H0_eV = cb.H0_eV
        cb.S_MN_T_MN_K = H0_eV / (KB_EV * max(S_kB_scale, 1e-6))
        cb.S_MN_sign = -1.0   # constrained crack-tip TS (toughening sign)
    return configure


def run_exp_floor(base_args, temps, systems, an_list, out_dir,
                  cleave_scale=3.6, gt_sign=1.0):
    print('=' * 72)
    print('  EXP_floor EMISSION + weak-modulus CLEAVAGE')
    print(f'  systems: {systems}   (a,n) pairs: {an_list}')
    print(f'  cleavage entropy scale: ~{cleave_scale} kB (modulus softening)')
    sgn = ('+1 (literal fit, S*<0)' if gt_sign > 0
           else '-1 (entropic-nucleation, S*>0)')
    print(f'  emission gT sign: {sgn}')
    print('  NOTE: a,n are PLACEHOLDERS pending the BarrierModel export.')
    print('=' * 72)
    series = {}
    cb_cfg = _cfg_cleave_modulus(cleave_scale)
    for system in systems:
        for (a, n) in an_list:
            rows = _kc_vs_T(_cfg_exp_floor_emit(system, a, n, gt_sign), temps,
                            base_args, configure_cb=cb_cfg)
            a_lab, n_lab, *_ = EXP_FLOOR_SYSTEMS[system][5:7] if a is None else (a, n)
            a_show = EXP_FLOOR_SYSTEMS[system][5] if a is None else a
            n_show = EXP_FLOOR_SYSTEMS[system][6] if n is None else n
            key = f'{system} a={a_show:.3g} n={n_show:.3g}'
            series[key] = rows
            dbtt = _dbtt(rows)
            kcs = ' '.join(f'{(r["Kc_MPa_sqrt_m"] or float("nan")):5.2f}' for r in rows)
            print(f'  {key:26}  DBTT={("none" if dbtt is None else f"{dbtt:.0f}K"):>6}  '
                  f'Kc(T): {kcs}')
    _save_series(out_dir, 'exp_floor', temps, series, key_label='config')
    _plot_series(out_dir, 'exp_floor', temps, series,
                 title='Kc(T): EXP_floor emission + weak-modulus cleavage',
                 legend_fmt=lambda k: str(k))
    return series


def run_crossover(base_args, temps, system, S_hs_list, sigma_S_list, out_dir,
                  cleave_scale=3.6, S_hs_power=3.0):
    """Sweep the high-stress entropy crossover (amplitude S_hs_kB x crossover
    stress sigma_S) and report where a DBTT-like transition emerges.  Emission =
    EXP_floor fit + crossover; cleavage = weak modulus entropy."""
    print('=' * 72)
    print(f'  HIGH-STRESS ENTROPY CROSSOVER SWEEP  (system={system})')
    print('  emission = EXP_floor fit + crossover; cleavage = weak modulus')
    print(f'  S_hs/kB: {S_hs_list}   sigma_S/GPa: {sigma_S_list}  (power={S_hs_power})')
    print('  S_hs=0 is the pure data fit (S*<0 everywhere).')
    print('=' * 72)
    series = {}
    cb_cfg = _cfg_cleave_modulus(cleave_scale)
    for S_hs in S_hs_list:
        for sigS in sigma_S_list:
            cfg = _cfg_exp_floor_emit(system, S_hs_kB=S_hs, sigma_S_GPa=sigS,
                                      S_hs_power=S_hs_power)
            rows = _kc_vs_T(cfg, temps, base_args, configure_cb=cb_cfg)
            key = f'S_hs={S_hs:g}kB sigS={sigS:g}GPa'
            series[key] = rows
            # emission-onset DBTT: lowest T whose tail all has N_em>5
            rs = sorted(rows, key=lambda r: r['T'])
            onset = None
            for i in range(len(rs)):
                if all(r['N_em'] > 5.0 for r in rs[i:]):
                    onset = rs[i]['T']
                    break
            kcs = ' '.join(f'{(r["Kc_MPa_sqrt_m"] or float("nan")):5.2f}' for r in rs)
            nems = ' '.join(f'{r["N_em"]:5.0f}' for r in rs)
            print(f'  {key:24} emit-onset={("none" if onset is None else f"{onset:.0f}K"):>6}')
            print(f'    {"Kc(T)":>10}: {kcs}')
            print(f'    {"N_em(T)":>10}: {nems}')
    _save_series(out_dir, 'crossover', temps, series, key_label='config')
    _plot_series(out_dir, 'crossover', temps, series,
                 title=f'Kc(T): high-stress entropy crossover ({system})',
                 legend_fmt=lambda k: str(k))
    return series


def _cfg_magnitude(S0_kB, pivot=True):
    """Constant-S baseline at magnitude S0_kB (matching the linear-G0 data).
    If pivot, co-set H0 so G*(sigma=0, T_ref) = W_G0_REF_EV (isolates the
    entropy T-slope from an overall level shift)."""
    def configure(eb):
        eb.use_negative_entropy = True
        eb.entropy_stress_form = 'physical'
        eb.S_T_c0_kB = S0_kB
        eb.S_T_c1_kB_per_K = 0.0
        eb.S_T_c2_kB_per_K2 = 0.0
        eb.S_T_min_kB = min(-60.0, S0_kB - 1.0)
        eb.S_T_max_kB = 0.0
        eb.S_sigma_max_kB = 0.0          # isolate the baseline; no stress gate
        if pivot:
            # G*(0,Tref) = H0 - Tref*S  with S = S0_kB*kB  (S0_kB<0)
            #   ->  H0 = G0_ref + Tref*S0_kB*kB_eV
            # For W[100] (S0=-45.7) this gives H0 ~ 0.05 eV: the barrier at Tref
            # is almost entirely entropic.  Beyond ~-45 kB the pivot enthalpy
            # would go NEGATIVE (unphysical) -- floor it at a small positive
            # value and note that the pivot can no longer be held there.
            H0 = W_G0_REF_EV + W_T_REF_K * S0_kB * KB_EV
            eb._pivot_H0_floored = bool(H0 < 0.02)
            eb.H0_eV = max(H0, 0.02)
    return configure


def _cfg_form(name):
    """Named functional forms for the emission entropy (fixed nominal scale)."""
    def configure(eb):
        eb.use_negative_entropy = True
        if name == 'const':
            eb.entropy_stress_form = 'physical'
            eb.S_T_c0_kB, eb.S_T_c1_kB_per_K, eb.S_T_c2_kB_per_K2 = -40.0, 0.0, 0.0
            eb.S_T_min_kB, eb.S_T_max_kB, eb.S_sigma_max_kB = -60.0, 0.0, 0.0
        elif name == 'linear_T':
            eb.entropy_stress_form = 'physical'
            eb.S_T_c0_kB, eb.S_T_c1_kB_per_K, eb.S_T_c2_kB_per_K2 = -46.0, 0.012, 0.0
            eb.S_T_min_kB, eb.S_T_max_kB, eb.S_sigma_max_kB = -60.0, 0.0, 0.0
        elif name == 'gated_schoeck':
            eb.entropy_stress_form = 'physical'
            eb.S_T_c0_kB, eb.S_T_c1_kB_per_K, eb.S_T_c2_kB_per_K2 = -40.0, 0.0, 0.0
            eb.S_T_min_kB, eb.S_T_max_kB = -60.0, 0.0
            eb.S_sigma_max_kB, eb.sigma0_S_GPa, eb.entropy_gate_power = 8.0, 3.0, 1.0
        elif name == 'meyer_neldel':
            eb.entropy_stress_form = 'meyer_neldel'
            eb.S_MN_T_MN_K, eb.S_MN_sign = 6500.0, 1.0
        elif name == 'poly_veverka':
            # MATLAB yield-fit T-quadratic baseline (A-coeffs), clipped; NO
            # exploding stress gain.  The low-stress fit example.
            eb.entropy_stress_form = 'physical'
            eb.S_T_c0_kB = -9.51
            eb.S_T_c1_kB_per_K = -2.55e-1
            eb.S_T_c2_kB_per_K2 = 3.45e-4
            eb.S_T_min_kB, eb.S_T_max_kB, eb.S_sigma_max_kB = -60.0, 0.0, 0.0
        else:
            raise ValueError(f'unknown form {name}')
    return configure


# ----------------------------------------------------------------------------
# drivers
# ----------------------------------------------------------------------------
def run_magnitude(base_args, temps, S0_list, out_dir):
    print('=' * 72)
    print('  EMISSION-ENTROPY MAGNITUDE SWEEP (pivot on W[100]: G0(481K)=1.94 eV)')
    print(f'  S0/kB levels: {S0_list}')
    print(f'  experimental emission anchors: {", ".join(f"{k}={v}" for k,v in EXP_ANCHORS.items())}')
    print('=' * 72)
    series = {}
    for S0 in S0_list:
        rows = _kc_vs_T(_cfg_magnitude(S0, pivot=True), temps, base_args)
        series[S0] = rows
        dbtt = _dbtt(rows)
        kcs = ' '.join(f'{(r["Kc_MPa_sqrt_m"] or float("nan")):5.2f}' for r in rows)
        print(f'  S0={S0:6.1f} kB  DBTT={("none" if dbtt is None else f"{dbtt:.0f}K"):>6}  '
              f'Kc(T): {kcs}')
    _save_series(out_dir, 'magnitude', temps, series, key_label='S0_kB')
    _plot_series(out_dir, 'magnitude', temps, series,
                 title='Kc(T) vs emission-entropy magnitude (W[100] pivot)',
                 legend_fmt=lambda k: f'S0={k:.1f} kB'
                 + ('  (W[100])' if abs(k + 45.7) < 0.6 else ''))
    return series


def run_forms(base_args, temps, forms, out_dir):
    print('=' * 72)
    print('  ENTROPY FORM COMPARISON (emission channel)')
    print(f'  forms: {forms}')
    print('=' * 72)
    series = {}
    for nm in forms:
        rows = _kc_vs_T(_cfg_form(nm), temps, base_args)
        series[nm] = rows
        dbtt = _dbtt(rows)
        kcs = ' '.join(f'{(r["Kc_MPa_sqrt_m"] or float("nan")):5.2f}' for r in rows)
        print(f'  {nm:14}  DBTT={("none" if dbtt is None else f"{dbtt:.0f}K"):>6}  '
              f'Kc(T): {kcs}')
    _save_series(out_dir, 'form', temps, series, key_label='form')
    _plot_series(out_dir, 'form', temps, series,
                 title='Kc(T) vs entropy functional form (emission channel)',
                 legend_fmt=lambda k: str(k))
    return series


# ----------------------------------------------------------------------------
# output helpers
# ----------------------------------------------------------------------------
def _save_series(out_dir, tag, temps, series, key_label):
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, f'entropy_sweep_{tag}.csv')
    with open(path, 'w', newline='') as fp:
        w = csv.writer(fp)
        w.writerow([key_label, 'T_K', 'Kc_MPa_sqrt_m', 'sigma_back_GPa',
                    'r_eff_over_r0', 'N_em', 'mode'])
        for k, rows in series.items():
            for r in rows:
                w.writerow([k, r['T'], r['Kc_MPa_sqrt_m'], r['sigma_back_GPa'],
                            r['r_eff_over_r0'], r['N_em'], r['mode']])
    with open(os.path.join(out_dir, f'entropy_sweep_{tag}.json'), 'w') as fp:
        json.dump({str(k): v for k, v in series.items()}, fp, indent=2)


def _plot_series(out_dir, tag, temps, series, title, legend_fmt):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except Exception:
        return None
    fig, ax = plt.subplots(figsize=(7.5, 5))
    for k, rows in series.items():
        rows = sorted(rows, key=lambda r: r['T'])
        Ts = [r['T'] for r in rows]
        Kc = [r['Kc_MPa_sqrt_m'] for r in rows]
        emph = isinstance(k, float) and abs(k + 45.7) < 0.6
        ax.plot(Ts, Kc, 'o-', lw=2.4 if emph else 1.4,
                color='k' if emph else None, label=legend_fmt(k))
        for r in rows:
            if r['mode'] == 'ductile' and r['Kc_MPa_sqrt_m'] is not None:
                ax.plot(r['T'], r['Kc_MPa_sqrt_m'], 's', ms=8,
                        mfc='none', mec='tab:red', mew=1.5)
    ax.set_xlabel('Temperature [K]')
    ax.set_ylabel(r'$K_{c,\,first}$ [MPa$\sqrt{m}$]')
    ax.set_title(title)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, ncol=2)
    fig.text(0.5, -0.02, 'red squares: ductile classification',
             ha='center', fontsize=8, color='tab:red')
    fig.tight_layout()
    p = os.path.join(out_dir, f'entropy_sweep_{tag}.png')
    fig.savefig(p, dpi=130, bbox_inches='tight')
    plt.close(fig)
    print(f'  saved {p}')
    return p


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--mode', choices=['magnitude', 'form', 'exp_floor',
                                       'crossover', 'both'], default='both')
    ap.add_argument('--temps', type=float, nargs='+',
                    default=[300, 400, 500, 600, 700, 800, 900, 1000])
    ap.add_argument('--S0-list', type=float, nargs='+', dest='S0_list',
                    default=[-5, -10, -20, -30, -40, -45.7, -50])
    ap.add_argument('--forms', nargs='+',
                    default=['const', 'linear_T', 'gated_schoeck',
                             'meyer_neldel', 'poly_veverka'])
    ap.add_argument('--ef-systems', nargs='+', dest='ef_systems',
                    default=['W[100]'])
    ap.add_argument('--ef-an', nargs='+', dest='ef_an',
                    default=['fit'],
                    help="comma-separated a,n pairs, or 'fit' to use the fitted "
                         "per-system exponents from the export")
    ap.add_argument('--cleave-scale', type=float, default=3.6, dest='cleave_scale')
    ap.add_argument('--ef-gt-sign', type=float, default=1.0, dest='ef_gt_sign',
                    help='+1 literal fit (S*<0); -1 entropic-nucleation (S*>0)')
    ap.add_argument('--xover-S-hs', type=float, nargs='+', dest='xover_S_hs',
                    default=[0, 30, 50, 70, 90],
                    help='(crossover) high-stress entropy shifts [kB]')
    ap.add_argument('--xover-sigma-S', type=float, nargs='+', dest='xover_sigma_S',
                    default=[5.0],
                    help='(crossover) crossover stresses [GPa]')
    ap.add_argument('--xover-power', type=float, default=3.0, dest='xover_power')
    ap.add_argument('--xover-system', default='W[100]', dest='xover_system')
    ap.add_argument('--out', default='runs/entropy_sweep')
    ap.add_argument('--Kdot', type=float, default=0.005)
    ap.add_argument('--dt', type=float, default=1.0)
    args = ap.parse_args(argv)

    # base 1D engine args (multihit ESSENTIAL or no DBTT resolves)
    base = _build_parser().parse_args(
        ['--mode', '1d', '--multihit-m', '3', '--multihit-tau', '1e-6',
         '--dt', str(args.dt), '--Kdot', str(args.Kdot)])

    os.makedirs(args.out, exist_ok=True)
    if args.mode in ('magnitude', 'both'):
        run_magnitude(base, args.temps, args.S0_list, args.out)
    if args.mode in ('form', 'both'):
        run_forms(base, args.temps, args.forms, args.out)
    if args.mode == 'exp_floor':
        an = [(None, None) if p == 'fit'
              else tuple(float(z) for z in p.split(',')) for p in args.ef_an]
        run_exp_floor(base, args.temps, args.ef_systems, an, args.out,
                      cleave_scale=args.cleave_scale, gt_sign=args.ef_gt_sign)
    if args.mode == 'crossover':
        run_crossover(base, args.temps, args.xover_system, args.xover_S_hs,
                      args.xover_sigma_S, args.out, cleave_scale=args.cleave_scale,
                      S_hs_power=args.xover_power)
    print(f'\n  outputs in {args.out}')


if __name__ == '__main__':
    main()
