# v10.0.5.9 — production-path J parity audit

v10.0.5.8 established a converged fixed-grip elastic reference for the actual
finite killed notch and showed that the standalone full-energy domain J agrees
with `-dU_el/da`.

v10.0.5.9 asks the remaining question: does the **actual production initialization
path** reproduce that reference on the same post-equilibrium elastic state?

This release does not change Arrhenius barriers, stochastic thresholds, MPZ
transport, cohesive evolution, crack-growth logic, material parameters, or the
production mesh generator.

## What the probe reuses

Each one-step probe passes through the audited v10.0.5.5 stochastic/VHCF source
transform stack and therefore uses the production:

- geometry and killed-notch stamp;
- graded mesh and physical refinement extent;
- symmetric displacement boundary conditions;
- anisotropic elasticity assembly;
- crack backend initialization;
- root-front versus straight-path contour policy;
- line-of-sight segment inventory;
- `2*kill_r` exclusion rule when the anisotropic root-front path is active.

Plastic evolution is replaced by an elastic no-op only in this audit entry.
The recorder is inserted immediately after the maximum-load FEM equilibrium solve.

## Primary metric

The anisotropy-safe comparison is

```text
J / sigma_gross^2
```

rather than an isotropic `J -> K` conversion. Two grip openings are run so that
`J/sigma^2` must also remain constant under elastic load scaling.

At 180, 240, and 300 micrometers the probe records:

- full stored-energy J with exact production domain semantics;
- tensile-filtered J with exact production domain semantics;
- full stored-energy J with the exclusion disk removed;
- tensile-filtered J with the exclusion disk removed;
- active element counts, effective notch tip, mesh spacing, refinement radius,
  reaction stress, elastic-energy closure and production contour metadata.

## Install

```bash
BASE=/Volumes/Data/Data/Nanopillar_calculation/Arrhenius_FEM_CZM_MPZ_v10_0_5_8_fixed_grip_elastic_convergence
NEW=/Volumes/Data/Data/Nanopillar_calculation/Arrhenius_FEM_CZM_MPZ_v10_0_5_9_production_j_parity
BRANCH=v10.0.5.9-production-j-parity-audit

git -C "$BASE" fetch origin \
  "+refs/heads/${BRANCH}:refs/remotes/origin/${BRANCH}"

git -C "$BASE" worktree add \
  -b "$BRANCH" \
  "$NEW" \
  "refs/remotes/origin/${BRANCH}"

cd "$NEW"
python -m pip install -e . --no-deps
```

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
  tests/test_v10059_production_j_parity.py
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
  --L-pz-um 20 \
  --mpz-length-um 100 \
  --mpz-n-bins 80 \
  --da-phys-um 5 \
  --crack-backend adaptive_czm \
  --crystal-aniso \
  --crystal-theta-deg 45 \
  --out runs/v10_0_5_9_production_j_parity_v1
```

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

## Fail-closed interpretations

```text
production_J_parity_passed
production_exclusion_disk_controls_J_mismatch
production_refinement_extent_or_mesh_support_mismatch
production_J_not_quadratic_in_load
production_elastic_energy_closure_failed
production_J_fixed_grip_mismatch_unresolved
```

Only `production_J_parity_passed` unblocks promotion of the v10.0.5.8
geometry-specific reference into a later stochastic first-passage wrapper.
All other statuses keep the bracket blocked and identify the next mechanics audit.
