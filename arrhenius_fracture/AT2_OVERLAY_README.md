# AT2 overlay model (prototype)  --  at2_overlay.py

A parallel fracture model that keeps the sharp-engine Arrhenius kinetics but lets
an AT2 damage field carry geometry, so branching, fragmentation, and spontaneous
nucleation come for free. Built from the design we converged on.

## The coupling, in one paragraph
The phi field is the source of truth for GEOMETRY; per-tip FrontEngine instances
are the source of truth for KINETICS. Each step: extract tips from phi; give each
tip its isolated K plus the PAIRWISE elastic (Kachanov) interaction from OTHER
cracks (collinear shields, parallel amplifies) -> K_eff; each tip's first-passage
clock AUTHORIZES an advance da_auth (its own labeled purse); the AT2 field SPENDS
exactly that, advancing each tip by da_auth along the AT2-preferred direction
(magnitude from the ledger, shape from AT2). A branch SPLITS its parent's purse
(each lobe gets half da_auth) so branching cannot manufacture dissipation. Surface
the field wants where no tip purse reaches is referred to a BULK nucleation hazard.

## Conservation gate (type A, length form)
The handshake is asserted as a single scalar per step:
    residual = (realized advance length) - (sum of authorized da_auth)  -> 0.
Energy is Gc*length on both sides; matching lengths matches energy. In the demo
the residual median ~0.06 (discretization: da_auth rounded to whole cells; spikes
when an authorized advance lands on already-cracked cells). The AT2 *functional*
value differs from Gc*length by the AT2 normalization constant c_w -- see TODO.

## Run
    python3 -c "from arrhenius_fracture.at2_overlay import run_demo; run_demo()"
Outputs at2_overlay.png: phi evolution (branching + nucleated secondary cracks),
Gc(x) microstructure, the advance handshake (authorized vs realized), and the
crack-inventory count. Demo: edge notch at 900 K, rising load; inventory grows
1 -> ~12 cracks with the main crack branching and secondaries nucleating on weak
(low-Gc) features.

## What is faithful vs reduced (all flagged in code with [TAGS])
Faithful: per-tip Arrhenius first-passage kinetics (the regime physics), labeled
local purses, branch-splits-parent-purse, pairwise elastic interaction, bulk
nucleation hazard, propagation/nucleation as a financed decision, length-form
conservation gate, spatial engine persistence across unstable tip relabeling.

Reduced / deferred (revisit if they bite):
  [FEM]    reduced analytic K (DrivingForce.k_iso) and a reduced sigma proxy for
           nucleation. Hook: swap in fem.py + j_integral.py for K and the FEM
           principal-stress field. THIS IS THE BIGGEST APPROXIMATION.
  [ITER]   one pass per load step; the Delta2 fixed-point (re-read K after AT2
           moves the geometry, because compliance changed) is stubbed (n_inner=1).
  [c_w]    the AT2 surface FUNCTIONAL value is not yet calibrated to Gc*length by
           the standard AT2 normalization constant; the length-form gate sidesteps
           this, but reconcile before quantitative energy claims.
  [PZMERGE]pairwise elastic K is used up to process-zone contact; AT2 label-merge
           handles coalescence; the smooth handoff is crude.
  [3BODY]  pairwise interaction only; degrades in dense knots (the fragmentation
           regime) -- per our decision, accepted for now.
  [KERNEL] isotropic L_pz purse kernel; anisotropic (forward-biased) deferred.
  [PLAST]  no spatial bulk plastic field fed back into phi; per-tip engines still
           carry their internal shielding/embrittlement ledgers, so tip KINETICS
           are intact, but there is no plastic-zone field shaping phi yet.
  [AMR]    uniform grid; fragmentation needs adaptive refinement following active
           fronts (ell <~ L_pz only where needed). Without AMR, fine-everywhere is
           the cost ceiling.

## Validation ladder (next, in order)
 1. wire FEM K (replace [FEM]) and check a single isolated nucleus reproduces the
    sharp-2D regime TRENDS (DBTT rises, ceramic falls) -- directional, not numeric.
 2. turn on the Delta2 iteration ([ITER]); confirm the length gate tightens.
 3. nucleation density vs stress / rate / microstructure (no branching).
 4. branching morphology + fragment-size statistics; AMR for density.

## Fix log: tip propagation vs runaway nucleation
First demo nucleated cracks ~uniformly and the main tip stalled. Root cause: the
reduced stress field was the UNIFORM far-field sigma_app (no tip concentration),
so the bulk nucleation hazard fired everywhere and the nucleated cracks then
SHIELDED the main tip (pairwise interaction), stalling it. Fixes (all in the
reduced proxy -- the real cure is the FEM stress field, [FEM]):
  * near-tip K-field concentration sigma ~ K_eff/sqrt(2 pi r) ahead of each tip,
    capped at the cohesive ceiling (not singular) -> the tip is the dominant
    stress concentrator and propagates.
  * nucleation excluded within a process-zone MARGIN of all damage -> the K-field
    halo feeds tip propagation, not satellite cracks; nucleation fires only in the
    genuine far field at weak features.
  * selective barrier (dG0=1.7 eV, sigma_star=6 GPa) -> a few secondary cracks at
    the weakest stressed sites, not a swarm.
Result: main crack propagates across the domain; ~5 spontaneous nucleations at
weak features away from the tip; length-gate residual median ~0.06. The DENSITY of
secondary cracks is only qualitative until the FEM stress field replaces the proxy.
