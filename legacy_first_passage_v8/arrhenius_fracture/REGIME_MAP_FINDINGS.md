# Temperature-dependent fracture regimes from the Arrhenius-hazard sharp front

## What was added (both backward-compatible; defaults reproduce prior results byte-for-byte)

Two knobs in `FrontConfig` / `_build_parser` / `build_engine`, used only in `lambda_cleave`:

1. `--cleave-shield-chi` (`chi_shield`, default 0.0)
   Wires the introduction's `sigma -> sigma - sigma_back` substitution into the
   crack-opening hazard:  `sigma_eff = max(sigma_tip - chi_shield*sigma_back, 0)`,
   and the cleavage barrier is evaluated at `sigma_eff`. This is the only channel
   that toughens WITH temperature (sigma_back grows with the emitted-dislocation
   ledger N_em, which turns on with T). Linear in N_em, so it competes head-to-head
   with the linear embrittlement.

2. `--emb-sat-frac` (`emb_sat_frac`, default 1.0)  [SUPERSEDED for physics by `--n-sat`/`--recover-k`; keep emb_sat_frac=1 in studies]
   Caps the stored-energy embrittlement: `dG_emb -> min(dG_emb, emb_sat_frac*Gstar)`.
   1.0 = uncapped (as shipped). <1 lets a tough upper shelf survive (see below).

`chi_shield = 0` and `emb_sat_frac = 1` reproduce the as-shipped (ceramic-only) model.

## The mechanism (answers "these phenomena arise for different reasons - which?")

Toughness trend is a competition, at every T, between two emission-driven channels
acting on the SAME cleavage barrier, both growing with T:

  shielding      sigma_eff = sigma_tip - chi*sigma_back   RAISES the barrier
                 (bounded above by the zero-stress barrier H0)
  embrittlement  dG_emb ~ N_em                            LOWERS the barrier
                 (as shipped UNBOUNDED -> floors Geff at 0 -> always wins high-T)

Key structural asymmetry, established by direct sweeps:
  * Shielding only raises the stress->barrier mapping, capped by H0.
  * Embrittlement subtracts barrier energy and is floored at Geff=0, i.e. unbounded
    downward. Because emission runs away at high T, UNBOUNDED embrittlement ALWAYS
    wins the high-T limit -- raising cleave_H0 from 3 to 10 eV barely lifts the
    high-T floor (0.4 -> 2.2). So a tough upper shelf is impossible unless dG_emb
    saturates (emb_sat_frac < 1). This is the single most important finding.

## The two governing axes

  emb_sat_frac (x H0)  -> SHELF SURVIVAL. Residual barrier under full shielding is
                          ~ (1 - emb_sat_frac)*H0. Unbounded -> high-T collapse
                          (ceramic / peak). Bounded -> tough shelf (DBTT).
  chi_shield           -> ONSET / low-T toughening. Cancels the intrinsic cleavage
                          entropy softening (cb.S_T_c0_kB = -2) and sets the
                          transition temperature.

## Regimes and the parameter that drives each

  ceramic (reduced toughness with T): embrittlement dominates (emb_sat ~ 1, low chi).
      Kc falls monotonically. Steeper for larger H0 (wider cold-hot spread).
  toughness peak: emb_sat ~ 1 + moderate chi. Shielding wins mid-T, unbounded
      embrittlement reclaims high-T -> Kc rises then collapses.
  weakly T-dependent: chi tuned to cancel the intrinsic softening over the low-T
      branch (flat lower shelf). The narrow boundary between ceramic and DBTT.
  DBTT (brittle -> tough): bounded embrittlement (emb_sat < ~0.5) + shielding.
      Brittle/flat low-T plateau, tough (arrest) above the onset.

## Honest limitations / next steps for a textbook-clean DBTT

  * The cleavage first-passage clock is effectively bimodal at high T (fires fast
    or arrests) because the Arrhenius rate is exponentially sensitive, so the DBTT
    shows as plateau->arrest rather than a smooth low->high upturn. A finite
    high-T crack-velocity law (saturating at v_R, per the introduction's
    v_crack relation) would give a graded upper shelf instead of hard arrest.
  * The emission entropy crossover (eTc0=-20, eTc1=0.02, eSsig=8) is switch-like,
    producing a sharp spike at the onset T. A smoother S_emit(T) would smooth the
    peak / weak-T curves.
  * For clean regime-boundary mapping, tighten the T grid (currently 100 K) and
    fix a single first-passage Kc definition (done here: first advance).

## Reproduce

    cd <parent of arrhenius_fracture>
    REGIME_OUT=./regime_map python3 -m arrhenius_fracture.regime_map

Backward-compat (unchanged from as-shipped):
    python3 -m arrhenius_fracture.sharp_front --mode 1d --temperatures 700 800 900 1000 \
        --Kdot 0.005 --Kmax 40 --dt 1.0 --multihit-m 3 --multihit-tau 1e-6 \
        --emit-S-T-c0-kB=-20 --emit-S-T-c1=0.02 --emit-S-sigma-max-kB=8 --out /tmp/bc

## Update: finite upper shelf + crack-velocity law (resolving the "bimodal" caveat)

Testing showed the earlier "lower-shelf-to-arrest" behavior was a MEASUREMENT
artifact: the loading window (Kmax=40) sat below the upper shelf, and the driver's
early-arrest break censored it. Extending the window reveals a finite, two-shelf
DBTT, e.g. at chi=0.6, H0=6, n_sat=2e3: lower shelf ~24, transition ~860-900 K,
finite upper shelf ~46 MPa*sqrt(m) (with mild upper-shelf softening 46.6->44.9).
The transition spans ~40 K (graded, not a step); its width/position are set by the
emission entropy crossover (eTc0, eTc1).

Two engine changes (both backward-compatible; defaults reproduce prior results):
  * `--v-rayleigh` : finite crack-velocity ceiling (Region III). Bounds unstable
    fast advance; reports v_crack. Does NOT change first-passage Kc.
  * regime_map.py now uses a wide window (Kmax=120) and no censoring break, so the
    finite upper shelf is captured. The scan runs at dt=4 (dt-converged to 0.2%).

Net: the regime mechanism was never bimodal; the model already yields a finite
two-shelf DBTT. Remaining kinetic refinement (for a broader, R-curve-like upper
shelf) is a smoother S_emit(T) form, not a change to the selection mechanism.
