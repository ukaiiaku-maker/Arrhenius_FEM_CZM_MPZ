# AT2 overlay -- remaining updates before the code is "complete"

## STATUS (read this first)
- PROPAGATION is now robust + reproducible: tip direction is geometric/PCA-based,
  skimage-free, with a FIXED per-crack growth axis. The main crack crosses the
  domain for every tested seed and is environment-independent. Run --nucleation-off
  for the clean, low-residual baseline.
- MEANDER is OFF. Seven heuristic meander schemes were tried; all stall on some
  microstructures because a hand-coded path bypasses the energy minimization that
  should select it. TRUE meander requires the FEM-driven AT2 evolution (#1).
- NUCLEATION placement is qualitative and, with nucleation ON, the conservation
  residual rises (overlap of many straight advances). Same root cause: no real
  stress field. Both are fixed by #1.

## RESOLVED since first cut
- tip propagation reproducibility: growth direction now geometric (component
  centroid -> tip), independent of the skimage skeletonize version. (Was: tip
  authorized advance but never realized it on some machines; resid=1.0.)
- straight cracks: advance now WALKS cell-by-cell steering toward lower Gc, so the
  crack meanders along weak paths instead of running dead straight.
- corner nucleation: boundary band excludes seeding near the domain edges.
  (Placement away from edges is improved but not yet correct -- see #1.)


Ordered by priority. Tags match the [FLAGS] in at2_overlay.py.

## Tier 1 -- structural (the model is not quantitative until these are done)
1. [FEM] Replace the reduced analytic driving force with the real FEM solve +
   J-integral / configurational force (fem.py + j_integral.py, already used by
   sharp_front.run_2d). Hooks to replace: DrivingForce.k_iso and
   AT2Overlay._sigma_proxy. THIS IS THE KEY STEP -- a real stress field gives
   genuine tip concentration and far-field stress, and removes THREE current
   crutches at once: the near-tip K-field proxy, the pairwise-shielding patch,
   and the hand-tuned nucleation stress gate. Everything downstream rides on it.
2. Validation gate (do immediately after #1): a single isolated nucleus with
   branching suppressed must reproduce the sharp-2D regime TRENDS (DBTT rises,
   ceramic falls, peak/weak appear) -- directional, not numeric (branching
   dissipates more, so single-crack agreement is trend-only). If this fails the
   multi-crack results are meaningless.
3. [ITER] Implement the Delta2 fixed-point iteration (AT2Config.n_inner>1):
   re-read K after AT2 moves the geometry (compliance changed) within a step, so
   the budget handshake is a true fixed point rather than a single pass.

## Tier 2 -- energy accounting (needed for dissipation claims)
4. [c_w] Calibrate the AT2 surface FUNCTIONAL to Gc*length via the AT2
   normalization constant, so the functional-energy gate and the length-form gate
   agree. Until then only the length-form conservation gate is trustworthy.
5. Global energy balance: fold the nucleation channel (currently logged separately
   as dS_nucl) into a single conservation statement = propagation ledger +
   nucleation hazard ledger. Assert total dissipation closes per Delta2.
6. Rate-matching: align dt and dsigma so the loading rate equals the 1-D/2-D
   Kdot convention (0.005 MPa*sqrt(m)/s), so T-sweeps compare at fixed rate and
   line up with the sharp-model regime map.

## Tier 3 -- physics fidelity
7. [PLAST] Add a spatial bulk plastic field and feed it back into phi and the
   stress field. Currently only the per-tip engine ledgers carry shielding /
   embrittlement (tip kinetics are intact, but there is no plastic-zone field).
8. Energetic branch criterion: replace the heuristic "weak flanks -> split"
   trigger with a real bifurcation test (branch only when splitting the authorized
   advance lowers total energy). Current branching is microstructure-cued only.
9. [KERNEL] Anisotropic (forward-biased) purse kernel for cleavage, plus a real
   cleavage-orientation field in the microstructure (deferred with GBs).
10. a_eff: replace the bounding-box-diagonal crack-length proxy with crack length
    measured along the loading normal from the actual geometry.

## Tier 4 -- robustness / scale (needed for fragmentation statistics)
11. [AMR] Adaptive mesh refinement following active fronts so ell <~ L_pz only
    where needed. Uniform fine grid is the cost ceiling for many cracks.
12. [PZMERGE] Smooth handoff from pairwise elastic K to AT2 label-merge when two
    process zones physically contact (currently crude: elastic K up to contact,
    then label-merge).
13. [3BODY] Three-body interaction terms for dense crack clusters (currently
    pairwise Kachanov only; degrades in the most-fragmented knots).

## Tier 5 -- usability / studies
14. Expose the per-tip engine regime preset (chi_shield / N_sat / cleave_H0) on
    the CLI so ceramic / peak / weak-T / DBTT can be selected per run (currently
    fixed in _default_engine_args).
15. Output: per-crack tracking (length, orientation, growth rate) and
    fragment-size distribution vs strain rate / T as first-class diagnostics.
16. Seed-ensemble runner: average response over microstructure realizations
    (crack density, dissipation) since single seeds are noisy.

## Validation targets (no numeric targets -- trends only)
- single nucleus reproduces sharp-2D Kc(T) trends (gate, #2);
- nucleation density increases with stress and loading rate, sites on weak/stressed
  features;
- apparent toughness / total dissipation increases monotonically with branch
  density;
- fragment size decreases with increasing rate (qualitatively).
