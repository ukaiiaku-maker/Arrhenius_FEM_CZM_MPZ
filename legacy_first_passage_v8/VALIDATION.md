# Migration branch validation status

## Automated checks completed

### 1. Python compilation

```bash
python -m compileall -q arrhenius_fracture
```

Status: passed.

### 2. CZM topology and interface tests

```bash
pytest -q tests/test_czm_backend.py
```

Status: `5 passed`.

The tests verify:

- a crack event can split node topology and insert a cohesive interface while preserving bulk element count and Gauss-point indexing;
- a fully failed interface transmits no tensile cohesive force;
- an intact interface transmits traction and satisfies internal force balance;
- `adaptive_czm` inserts a 23.5 degree requested off-axis increment at 23.5 degrees with zero angular error while preserving element indexing;
- four successive changing off-axis increments are inserted at their requested angles without changing bulk element count.

### 3. Legacy backend regression smoke

The legacy `sharp_wake` backend completed a short 2-D run after the architecture changes.

### 4. Monotonic CZM event smoke

A low-barrier controlled test completed an Arrhenius first-passage event, inserted a topological cohesive crack segment, and wrote:

- `cohesive_elements.csv`;
- `czm_advance_log.json`.

### 5. Full inherited-physics integration smoke

```bash
OUT=/tmp/fem_czm_git_smoke KMAX=7.0 bash run_fem_czm_fullphysics_smoke.sh
```

Status: completed successfully.

This test ran the existing V1/V8 comparison path and the V8 2-D model with:

- fatigue cycle integration;
- adaptive cycle blocks;
- cyclic FEM mechanics;
- Arrhenius plasticity;
- spatial process-zone bookkeeping;
- anisotropic/branch-capable front inventory;
- the `edge_split_czm` geometry backend.

The smoke is an architecture/integration test, not a calibrated physics comparison. The script deliberately uses only three blocks and disables the expensive 2-D K calibration.

## Oriented DBTT integration checks

The canonical DBTT-like 300 K test was run with anisotropic elasticity and crystal orientation enabled.

- 30 degree crystal case: hazard-selected direction = 22.6671007 degrees; inserted adaptive-CZM direction = 22.6671007 degrees; first-passage toughness = 23.463963 MPa sqrt(m).
- 45 degree crystal case: hazard-selected direction = 31.3852164 degrees; inserted adaptive-CZM direction = 31.3852164 degrees; first-passage toughness = 23.459 MPa sqrt(m) in the pilot run.

The 30 degree first-passage toughness is unchanged from the prior edge-steered pilot while the inserted direction changed from 45 degrees to the exact 22.667 degree hazard-selected direction. This isolates the correction to post-event geometry rather than pre-event kinetics.

## Validation still required before broad production campaigns

1. Crack-path convergence with tip mesh refinement and `da_phys` refinement.
2. Multi-event monotonic propagation over substantial crack extension with mesh-quality diagnostics.
3. Branch-birth and daughter-topology regression tests.
4. Energy accounting across topology insertion and local r-adaptation events.
5. Comparison of `sharp_wake` and `adaptive_czm` for the four monotonic fracture classes.
6. Comparison of fatigue `da/dN` and threshold trends for the six canonical barrier systems.
7. Add full local patch retriangulation only if long tortuous paths exhaust the r-adaptive quality safeguards.
