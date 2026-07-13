"""
Regime map for the Arrhenius-hazard sharp-front engine.
================================================================================
The temperature trend of fracture toughness is set by a competition, at every T,
between two emission-driven channels acting on the SAME cleavage barrier. Both
grow with T (emission turns on with T):

    shielding      sigma_eff = sigma_tip - chi_shield * sigma_back   RAISES barrier
                   bounded above by the zero-stress barrier H0
    embrittlement  dG_emb ~ N_em (stored energy)                     LOWERS barrier
                   as shipped UNBOUNDED -> always wins the high-T limit;
                   capped at emb_sat_frac * barrier when emb_sat_frac < 1

Two knobs span every observed behavior, and they answer "what drives each regime":

    emb_sat_frac (x H0)  -> SHELF SURVIVAL.  Residual barrier under full shielding
                            is ~ (1 - emb_sat_frac)*H0.  Large (unbounded
                            embrittlement) -> high-T always collapses (ceramic/peak).
                            Small (bounded) -> a tough upper shelf survives (DBTT).
    chi_shield           -> ONSET / LOW-T toughening.  The only channel that
                            toughens WITH T; cancels the intrinsic cleavage
                            softening and sets the transition temperature.

Regimes:
    ceramic : embrittlement dominates           -> Kc falls with T
    peak    : shielding wins mid-T, unbounded
              embrittlement reclaims high-T      -> Kc rises then collapses
    weak-T  : chi cancels intrinsic softening    -> flat lower branch
    DBTT    : shielding + bounded embrittlement   -> brittle low-T, tough high-T

Kc = applied K at first cleavage advance in a K-ramp; no fire => arrest (tough).
"""
from __future__ import annotations
import os
import numpy as np
from .sharp_front import _build_parser, build_engine
from .config import make_emergent_config

KDOT_SI, KMAX_SI, DT = 0.005e6, 120.0e6, 4.0
TS = list(range(500, 1301, 100))


def base_args():
    a = _build_parser().parse_args(['--mode', '1d', '--temperatures', '700'])
    a.emit_S_T_c0_kB, a.emit_S_T_c1, a.emit_S_sigma_max_kB = -20.0, 0.02, 8.0
    a.multihit_m, a.multihit_tau = 3.0, 1e-6
    return a


def kc_curve(chi, nsat, H0, mat, reck=0.0):
    a = base_args(); a.chi_shield = chi; a.emb_sat_frac = 1.0; a.cleave_H0_eV = H0
    a.N_sat = nsat; a.recover_k = reck
    out = []; nstep = int(np.ceil(KMAX_SI / (KDOT_SI * DT)))
    for T in TS:
        eng = build_engine(a, mat); kc = None
        for i in range(nstep):
            K = KDOT_SI * (i + 1) * DT
            if eng.step(K, T, DT)['fired']:
                kc = K / 1e6; break
        out.append(kc)
    return out


def classify(Kcs):
    def m(lo, hi):
        v = [k for t, k in zip(TS, Kcs) if lo <= t <= hi and k is not None]
        return float(np.mean(v)) if v else None
    n_hi_none = sum(1 for t, k in zip(TS, Kcs) if t >= 1000 and k is None)
    lo, mid, hi = m(500, 700), m(800, 1000), m(1050, 1300)
    if all(k is None for k in Kcs):
        return 'tough'
    if lo is None:
        return 'tough'
    if hi is None or n_hi_none >= 2:
        return 'DBTT'                       # high-T didn't fire within window -> very tough shelf
    if mid is not None and mid > 1.3 * lo and hi < 0.65 * mid:
        return 'peak'
    if hi > 1.3 * lo:
        return 'DBTT'
    if hi < 0.7 * lo:
        return 'ceramic'
    return 'weak-T'


CODE = {'ceramic': 0, 'peak': 1, 'weak-T': 2, 'DBTT': 3, 'tough': 4, 'indet': -1}
SYM = {0: 'C', 1: 'P', 2: '~', 3: 'D', 4: 'T', -1: '?'}


def main():
    out = os.environ.get('REGIME_OUT', '/tmp/afx/regime_map')
    os.makedirs(out, exist_ok=True)
    mat = make_emergent_config().material

    # ---- Panel A: one verified representative per regime ----
    INF = float('inf')
    # (name, color, chi, n_sat, H0); emb cap OFF -> shelf controlled by recovery (n_sat)
    reps = [
        ('ceramic',  '#b22222', 0.00, INF,   2.6),   # no recovery -> embrittlement wins high-T
        ('peak',     '#e07b00', 0.10, INF,   3.6),   # no recovery + mild shield -> rise then collapse
        ('weak-T',   '#caa800', 0.20, 1.5e3, 4.0),   # recovery + light shield -> flat across range
        ('DBTT',     '#2e8b57', 0.60, 2.0e3, 6.0),   # recovery bounds dG_emb -> tough shelf survives
    ]
    print("=" * 88)
    print("  PANEL A: representative Kc(T) per regime")
    print("=" * 88)
    repcurves = []
    for name, c, chi, ns, H0 in reps:
        Kc = kc_curve(chi, ns, H0, mat); repcurves.append((name, c, chi, ns, H0, Kc))
        print(f"  {name:8s} (chi={chi}, n_sat={ns}, H0={H0}) [{classify(Kc)}]  "
              + " ".join('arr' if k is None else f'{k:5.1f}' for k in Kc))

    # ---- Panel B: (chi_shield x emb_sat_frac) phase grid at fixed H0 ----
    H0_B = 4.0
    chis = [0.0, 0.2, 0.4, 0.6, 0.9, 1.3]
    sats = [float("inf"), 8.0e3, 3.0e3, 1.5e3, 8.0e2]   # n_sat: unbounded -> tightly bounded
    print("\n" + "=" * 88)
    print(f"  PANEL B: phase grid (chi x n_sat) at H0={H0_B} eV, emb cap OFF (physical recovery)")
    print("=" * 88)
    grid = np.full((len(sats), len(chis)), -1, int)
    for ic, sat in enumerate(sats):
        row = []
        for ix, chi in enumerate(chis):
            grid[ic, ix] = CODE[classify(kc_curve(chi, sat, H0_B, mat))]
            row.append(SYM[grid[ic, ix]])
        tag = "inf  " if sat==float("inf") else f"{sat:5.0f}"
        print(f"  n_sat={tag} | " + "  ".join(f"chi{chi}:{ss}" for chi, ss in zip(chis, row)))

    # ------------------------------- figure -------------------------------
    import matplotlib
    matplotlib.use('Agg'); import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, BoundaryNorm
    from matplotlib.patches import Patch
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(14.2, 5.6))
    KmaxM = KMAX_SI / 1e6
    for name, c, chi, ns, H0, Kc in repcurves:
        Tf = [t for t, k in zip(TS, Kc) if k is not None]
        Kf = [k for k in Kc if k is not None]
        Ta = [t for t, k in zip(TS, Kc) if k is None]
        axA.plot(Tf, Kf, 'o-', color=c, ms=5, lw=2, label=name)
        if Ta:
            axA.plot(Ta, [KmaxM] * len(Ta), '^', color=c, ms=9, mfc='none')
    axA.axhline(KmaxM, color='0.6', ls=':', lw=1)
    axA.text(510, KmaxM - 2.4, 'arrest (tough / no advance)', color='0.45', fontsize=9)
    axA.set_xlabel('Temperature [K]'); axA.set_ylabel(r'$K_c$ [MPa$\sqrt{\mathrm{m}}$]')
    axA.set_title('A.  three target behaviors (+ peak)\nshelf set by recovery, no artificial cap')
    axA.legend(fontsize=10, loc='center right'); axA.grid(alpha=0.3)

    reg_cmap = ListedColormap(['#b22222', '#e07b00', '#caa800', '#2e8b57', '#15543a'])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5, 4.5], reg_cmap.N)
    axB.imshow(np.clip(grid, 0, 4), origin='lower', aspect='auto', cmap=reg_cmap, norm=norm,
               extent=[-0.5, len(chis) - 0.5, -0.5, len(sats) - 0.5])
    axB.set_xticks(range(len(chis))); axB.set_xticklabels(chis)
    axB.set_yticks(range(len(sats))); axB.set_yticklabels(['inf' if v==float('inf') else f'{v:.0f}' for v in sats])
    axB.set_xlabel(r'$\chi_{shield}$  (onset / low-T toughening $\to$)')
    axB.set_ylabel('n_sat (saturation density)\n(no recovery inf $\\to$ strong recovery)')
    axB.set_title(f'B.  phase diagram, physical recovery ($H_0$={H0_B} eV)')
    for ic in range(len(sats)):
        for ix in range(len(chis)):
            axB.text(ix, ic, SYM[grid[ic, ix]], ha='center', va='center',
                     color='white', fontweight='bold', fontsize=12)
    axB.legend(handles=[Patch(color='#b22222', label='C ceramic (softens with T)'),
                        Patch(color='#e07b00', label='P toughness peak'),
                        Patch(color='#caa800', label='~ weakly T-dependent'),
                        Patch(color='#2e8b57', label='D DBTT (brittle $\\to$ tough)'),
                        Patch(color='#15543a', label='T fully tough')],
               loc='lower right', fontsize=7.5, framealpha=0.95)
    fig.suptitle('Arrhenius-hazard sharp front: shielding (chi) vs density recovery (n_sat) '
                 'select the T-dependence of toughness  [no artificial caps]', fontsize=11, y=1.02)
    fig.tight_layout()
    p = os.path.join(out, 'regime_map.png'); fig.savefig(p, dpi=130, bbox_inches='tight'); plt.close(fig)
    print(f"\n  saved {p}")


if __name__ == '__main__':
    main()
