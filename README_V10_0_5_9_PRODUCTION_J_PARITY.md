# v10.0.5.9 — production-path J parity audit

v10.0.5.8 established a converged fixed-grip elastic reference for the finite
killed notch and showed that the standalone full-energy domain J agrees with
`-dU_el/da`.

v10.0.5.9 asks whether the actual v10.0.5.x production initialization reproduces
that reference on the same post-equilibrium elastic state.

This release does not change Arrhenius barriers, stochastic thresholds, MPZ
transport, cohesive evolution, crack-growth logic, material parameters, or the
production mesh generator.

## Active production semantics

Each one-step probe traverses the audited v10.0.5.5 stochastic/VHCF source-
transform stack and therefore uses the production:

- geometry and killed-notch stamp;
- graded mesh and physical refinement extent;
- symmetric displacement boundary conditions;
- anisotropic elasticity assembly;
- kinetic-trial adaptive-CZM backend initialization;
- line-of-sight segment inventory;
- maximum-load FEM equilibrium state.

The v10.0.3 progressive adapter deliberately preserves anisotropic FEM and J but
sets the crack-path selector to `deflect=False` for the straight, single-front
checkpoint lifecycle. Consequently, the active v10.0.5.x production J contour is
`straight_progressive_cluster_no_exclusion`; the root-front `2*kill_r` exclusion
branch is not active for this parity run.

Plastic evolution is replaced by an elastic no-op only in this audit entry. The
recorder is inserted immediately after the maximum-load FEM equilibrium solve.

## Primary metric

The anisotropy-safe comparison is

```text
J / sigma_gross^2
```

rather than an isotropic `J -> K` conversion. Two grip openings are run so that
`J/sigma^2` must remain constant under elastic load scaling.

At 180, 240, and 300 micrometers the probe records:

- full stored-energy J with exact production domain semantics;
- tensile-filtered J with exact production domain semantics;
- no-exclusion J at the identical state;
- active element counts, effective notch tip, mesh spacing, refinement radius,
  reaction stress, elastic-energy closure, and production contour metadata.

For the active straight progressive path, production and no-exclusion values
should be identical. The root-front exclusion calculation remains represented in
the source-transform audit only for compatibility with non-progressive paths.

## Validation

```bash
python -m py_compile \
  arrhenius_fracture/production_j_parity_v10059.py \
  arrhenius_fracture/mode_i_first_passage_v10_0_5_9_production_j_probe.py \
  run_v10_0_5_9_production_j_parity.py

pytest -q \
  tests/test_v10056_kj_audit_bracket.py \
  tests/test_v10056_audited_wrappers.py \
  tests/test_v10058_fixed_grip_elastic_audit.py \
  tests/test_v10058_fixed_grip_audit_mesh.py \
  tests/test_v10059_production_j_parity.py \
  tests/test_v10059_v911_probe_contract.py \
  tests/test_v10059_straight_path_recorder.py
```

## Run

```bash
REFERENCE=/Volumes/Data/Data/Nanopillar_calculation/Arrhenius_FEM_CZM_MPZ_v10_0_5_8_fixed_grip_elastic_convergence/runs/v10_0_5_8_fixed_grip_elastic_convergence_v2/fixed_grip_elastic_convergence_v10_0_5_8.json

python run_v10_0_5_9_production_j_parity.py \
  --reference-json "$REFERENCE" \
  --openings-um "1 2" \
  --contour-outer-um "180 240 300" \
  --selected-outer-um 240 \
  --nx 36 \
  --ny 72 \
  --tip-h-um 2.5 \
  --tip-ratio 1.15 \
  --cluster-J-outer-um 240 \
  --local-J-outer-um 100 \
  --mpz-length-um 100 \
  --mpz-n-bins 80 \
  --da-phys-um 5 \
  --crack-backend adaptive_czm \
  --crystal-aniso \
  --crystal-theta-deg 45 \
  --out runs/v10_0_5_9_production_j_parity_v5
```

The probe entry inserts the v9.11-required `--crystal-compete` switch and enforces
branching disabled with `--max-fronts 1`. It also sends the authoritative
`--mpz-length-um` value through the legacy `--L-pz` interface so the two lengths
cannot conflict.

## Outputs

```text
production_j_parity_v10_0_5_9.json
production_j_probe_summary_v10_0_5_9.csv
production_j_contours_v10_0_5_9.csv
opening_1um/production_j_probe_v10_0_5_9.json
opening_1um/production_j_probe.log
opening_2um/production_j_probe_v10_0_5_9.json
opening_2um/production_j_probe.log
```

A launch failure also writes:

```text
production_j_probe_launch_failure_v10_0_5_9.json
```

## Fail-closed interpretations

```text
production_J_parity_passed
production_exclusion_disk_controls_J_mismatch
production_refinement_extent_or_mesh_support_mismatch
production_J_not_quadratic_in_load
production_elastic_energy_closure_failed
production_J_fixed_grip_mismatch_unresolved
```

For the active straight progressive path,
`production_exclusion_disk_controls_J_mismatch` should not occur because the
production exclusion radius is zero. Only `production_J_parity_passed` unblocks
promotion of the v10.0.5.8 geometry-specific reference into a later stochastic
first-passage wrapper.
