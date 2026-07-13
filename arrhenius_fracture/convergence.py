"""
Convergence / robustness tests for the sharp-front regime results.

The suspicious numerics are (i) the cleavage first-passage CLOCK integration over
the step (dt), and (ii) the emission-rate throttle dN_cap, which interacts with dt
because dN = min(lam_e*dt, dN_cap). We also check the loading window (Kmax) and the
model-sensitivity to the correlation window tau_c.

A result is trustworthy only where Kc is flat against dt and dN_cap. Run at a
PHYSICAL regime point (recovery on, artificial barrier cap OFF: emb_sat_frac=1).
"""
from __future__ import annotations
import numpy as np
from .sharp_front import _build_parser, build_engine
from .config import make_emergent_config

mat = make_emergent_config().material


def base():
    a = _build_parser().parse_args(['--mode', '1d', '--temperatures', '700'])
    a.emit_S_T_c0_kB, a.emit_S_T_c1, a.emit_S_sigma_max_kB = -20.0, 0.02, 8.0
    a.multihit_m, a.multihit_tau = 3.0, 1e-6
    a.emb_sat_frac = 1.0                 # artificial barrier cap OFF
    a.cleave_H0_eV = 6.0
    a.chi_shield = 0.6
    a.N_sat = 2.0e3                      # physical density ceiling
    return a


def kc(a, T, Kdot=0.005e6, Kmax=40e6, dt=1.0):
    eng = build_engine(a, mat)
    nstep = int(np.ceil(Kmax / (Kdot * dt)))
    chi = getattr(a, 'chi_shield', 0.0)
    sh = 0
    for i in range(nstep):
        K = Kdot * (i + 1) * dt
        sig = eng.sigma_tip(K); seff = max(sig - chi * eng.sigma_back(), 0.0)
        if eng.step(K, T, dt)['fired']:
            return K / 1e6
        if K > 8e6 and seff < 0.02 * max(sig, 1.0):
            sh += 1
            if sh > int(250 / dt):
                return None
        else:
            sh = 0
    return None


def fmt(v):
    return 'arr' if v is None else f'{v:6.2f}'


def main():
    Ttest = [600, 800, 900]    # cold shelf, pre-onset, transition
    print("Regime point: chi=0.6, H0=6.0 eV, n_sat=2e3, emb_sat_frac=1 (cap OFF)\n")

    print("== (1) timestep dt  [clock integration]  default dN_cap=50 ==")
    print("  dt[s]   " + "  ".join(f"T={T}" for T in Ttest))
    for dt in [4.0, 2.0, 1.0, 0.5, 0.25]:
        a = base()
        print(f"  {dt:5.2f}   " + "  ".join(fmt(kc(a, T, dt=dt)) for T in Ttest))

    print("\n== (1b) dt with dN_cap=1e9 [remove emission throttle so dt is isolated] ==")
    print("  dt[s]   " + "  ".join(f"T={T}" for T in Ttest))
    for dt in [4.0, 2.0, 1.0, 0.5, 0.25]:
        a = base(); a.dN_cap = 1e9
        # dN_cap is set on cfg in build_engine via f.dN_cap default; override:
        eng_probe = build_engine(a, mat); eng_probe.f.dN_cap = 1e9
        # rebuild path: patch kc to use a fresh engine each T with dN_cap set
        def kc_big(T, dt=dt):
            eng = build_engine(a, mat); eng.f.dN_cap = 1e9
            nstep = int(np.ceil(40e6 / (0.005e6 * dt))); chi = a.chi_shield; sh = 0
            for i in range(nstep):
                K = 0.005e6 * (i + 1) * dt
                sig = eng.sigma_tip(K); seff = max(sig - chi * eng.sigma_back(), 0.0)
                if eng.step(K, T, dt)['fired']:
                    return K / 1e6
                if K > 8e6 and seff < 0.02 * max(sig, 1.0):
                    sh += 1
                    if sh > int(250 / dt):
                        return None
                else:
                    sh = 0
            return None
        print(f"  {dt:5.2f}   " + "  ".join(fmt(kc_big(T)) for T in Ttest))

    print("\n== (2) emission throttle dN_cap [per-step cap on emitted count] ==")
    print("  dN_cap   " + "  ".join(f"T={T}" for T in Ttest))
    for cap in [25.0, 50.0, 100.0, 200.0, 1e9]:
        def kc_cap(T, cap=cap):
            eng = build_engine(base(), mat); eng.f.dN_cap = cap
            nstep = int(np.ceil(40e6 / (0.005e6 * 1.0))); chi = 0.6; sh = 0
            for i in range(nstep):
                K = 0.005e6 * (i + 1)
                sig = eng.sigma_tip(K); seff = max(sig - chi * eng.sigma_back(), 0.0)
                if eng.step(K, T, 1.0)['fired']:
                    return K / 1e6
                if K > 8e6 and seff < 0.02 * max(sig, 1.0):
                    sh += 1
                    if sh > 250:
                        return None
                else:
                    sh = 0
            return None
        print(f"  {cap:7.0f}  " + "  ".join(fmt(kc_cap(T)) for T in Ttest))

    print("\n== (3) loading window Kmax [Kc must be < Kmax to be resolved] ==")
    print("  Kmax    " + "  ".join(f"T={T}" for T in Ttest))
    for Kmax in [30e6, 40e6, 60e6, 80e6]:
        a = base()
        print(f"  {Kmax/1e6:5.0f}   " + "  ".join(fmt(kc(a, T, Kmax=Kmax)) for T in Ttest))

    print("\n== (4) correlation window tau_c [MODEL sensitivity, not convergence] ==")
    print("  tau_c    " + "  ".join(f"T={T}" for T in Ttest))
    for tau in [1e-7, 1e-6, 1e-5]:
        a = base(); a.multihit_tau = tau
        print(f"  {tau:.0e}  " + "  ".join(fmt(kc(a, T)) for T in Ttest))


if __name__ == '__main__':
    main()
