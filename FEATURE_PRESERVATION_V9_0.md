# Feature-preservation map — v8 to v9

The v9 change is intentionally localized to the front-local constitutive state. The spatial solver is not replaced by a reduced or single-front implementation.

| Existing capability | v9 status | Integration point |
|---|---|---|
| Cubic anisotropic elasticity and crystal rotation | Preserved | Existing FEM/material construction in `sharp_front.py` |
| Cleavage-system and slip-system directional competition | Preserved | Existing crystal candidate generation and mixed-mode directional factors |
| Signed vector/domain J-integral | Preserved | Existing `j_integral.py` and directional projection logic |
| Local and clustered J evaluation | Preserved | Existing multifront J selection and cluster handoff |
| Multiple active fronts | Preserved | Existing front inventory; each front now owns its own MPZ state |
| Crack deflection | Preserved | Existing directional continuation |
| First-passage branch birth | Preserved | Existing branch clocks and geometric checks |
| Parent/daughter state transfer | Strengthened | `clone_split()` now conserves finite sites, mobile/retained fields, slip, and cumulative inventories |
| Adaptive sharp-wake backend | Preserved | Existing backend interface unchanged |
| Adaptive cohesive/topological backend | Preserved | Existing CZM insertion/remeshing unchanged |
| Geometry-event veto | Strengthened | Full MPZ state is rolled back, not only scalar `N_em` |
| Crack coalescence | Preserved | Existing coalescence/network logic unchanged |
| Side-branch retirement | Preserved | Existing front retirement criteria unchanged |
| Monotonic loading | Preserved | Same driver, new selectable front state |
| Cyclic mechanics | Preserved | Existing 2-D cyclic field updates remain; front hazard state is unified |
| Fatigue cycle blocking | Preserved | Existing adaptive controller delegates MPZ state evolution to the front |
| Dwell/creep loading | Added | Uses the same material row and MPZ kinetics |
| Mixed-mode v8 production wrapper | Preserved and patched | Dynamic wrapper retains the actual MPZ engine subclass |
| Field snapshots and diagnostics | Preserved and extended | Added MPZ profile JSON histories and aggregate columns |
| Legacy regression | Preserved | `legacy_scalar` remains the default; exact input archives are frozen |
| Stateful local-peridynamics crack formation | Preserved and repaired | Restored the missing shared intact-FEM helper module; candidate-site, embryo, stabilization/healing, bond-softening, and graph-connectivity tests pass |

## Deliberate validation restriction

`run_mpz_fem_czm_validation_matrix.py` defaults to `--max-fronts 1` only for the first constitutive validation gate. This isolates initiation and propagation-resistance persistence from changes in branch count. Passing `--enable-branching` uses the same full solver, same material parameters, branch clocks, coalescence, and multifront bookkeeping.

## Historical modules

Mixed-mode modules v1–v7 are retained unchanged for historical reproducibility. The active production mixed-mode module, v8, is MPZ-compatible. No historical module is deleted.
