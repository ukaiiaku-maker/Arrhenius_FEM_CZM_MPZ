# v9.11 full 2-D FEM/CZM integration contract

## Correction of scope

This branch modifies the existing `Arrhenius_FEM_CZM_MPZ` anisotropic 2-D
FEM/CZM solver. It does not replace that solver with a reduced front surrogate.

The separate `FEM-CZM-fatigue-fracture` repository is not an implementation
source for this branch. Its reduced constitutive driver is retained only as an
abandoned diagnostic artifact and must not be used for validation or production.

## Non-negotiable retained 2-D features

The following existing features remain active and must pass regression tests:

- elastic-plastic finite-element equilibrium on the evolving cracked geometry;
- domain/interaction-integral driving force and mixed-mode decomposition;
- anisotropic elasticity and crystallographic directional competition;
- `sharp_wake`, `edge_split_czm`, and `adaptive_czm` geometry backends;
- zero-thickness cohesive topology insertion after Arrhenius first passage;
- unilateral crack-face contact where available in the active backend;
- multi-front inventory, branch clocks, cluster/local-J handoff, retirement,
  coalescence bookkeeping, and conservative branch-state splitting;
- Gauss-point plastic strain, density, residual-stress, and process-zone history;
- topology/history transfer during local refinement;
- monotonic adaptive event stepping without hazard discard;
- fatigue cycle quadrature, adaptive cycle blocks, checkpoint/restart, and
  publication diagnostics.

## Constitutive target

The selected ceramic and weakT parameterizations are the accepted v9.10.2
independent-shape candidates. The selected DBTT parameterization is the accepted
v9.10.3 target-aware candidate, using the same v9.10.2 constitutive state.

Each front retains:

- finite source-site inventory and crack-advance source refresh;
- independent cleavage and emission EXP-floor surfaces;
- independent Peierls and Taylor `H0`, activation entropy, `alpha`, and `n`;
- detailed-balance Peierls transport;
- geometric forest-encounter retention;
- detailed-balance correlated Taylor release with uncapped natural hit order;
- mobile and retained spatial fields;
- accumulated-slip blunting;
- direct retained-line `K_shield` evaluated once;
- moving-frame state translation and wake accounting.

No active constitutive path may reintroduce:

- scalar `N_sat` production saturation;
- emission increment caps;
- mobile-density saturation or floor;
- minimum jump-length floor;
- Taylor hit-order or stress-amplification cap;
- empirical scalar backstress from unsigned density;
- stored-energy lowering of the cleavage barrier;
- a second cohesive/traction failure criterion.

## Current full-repository integration gaps

1. `build_engine(..., front_state_model='moving_pz')` imports the package-level
   `moving_process_zone` module rather than the final independent-shape v9.10.2
   state.
2. The 2-D bulk plasticity initialization in `sharp_front.run_2d` still builds
   Peierls/Taylor kinetics from emission-derived energy/entropy ratios and
   legacy defaults for saturation, density floor, jump floor, and Taylor
   amplification.
3. The final three-class exporter and promotion manifests are not yet the single
   authoritative parameter source for both front-local and bulk 2-D kinetics.
4. The 2-D process-zone sampler/coupler must pass the FEM-derived local stress
   shape and scalar forest-density baseline into each front without converting
   unsigned scalar density into signed shielding or subtracting bulk shielding
   twice.
5. Full 2-D parity tests and one-angle/four-temperature validation have not yet
   been run with the final parameter sets.

## Required implementation sequence

1. Add a versioned v9.11 front-state adapter around the final v9.10.2 moving
   process zone; do not overwrite legacy modules.
2. Add a versioned v9.11 bulk Peierls/Taylor adapter using independent shapes and
   the same selected parameter record.
3. Wire both adapters into `sharp_front.run_2d` behind an explicit
   `--front-state-model moving_pz_v911`/production preset while preserving all
   geometry and solver code.
4. Add static checks that legacy caps/floors and empirical scalar shielding are
   inactive.
5. Add 1-D constitutive parity tests at prescribed stress histories.
6. Add 2-D coupling tests for stress-profile normalization, density-floor
   bookkeeping, retained-only direct shielding, and state conservation through
   advance/refinement.
7. Run `sharp_wake` versus `adaptive_czm` first-passage parity with branching
   disabled.
8. Run ceramic/weakT/DBTT Mode-I gates at 300, 700, 900, and 1200 K.
9. Run one mixed-mode angle per class, then long R-curve and fatigue gates.
10. Re-enable and validate branching/coalescence only after the single-front
    constitutive and geometric parity gates pass.

No source-law redesign or additional fitting is permitted before these gates.
