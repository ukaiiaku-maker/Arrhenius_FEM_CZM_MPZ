# v10.0.5.10 — physical refinement-support audit

v10.0.5.9 established that the actual v10.0.5.x production path is linear,
energy-consistent, and quadratic in load, but its 240 micrometer J/sigma^2 is only
0.2767 of the converged v10.0.5.8 fixed-grip reference. The active straight
progressive path uses no exclusion disk, so the mismatch cannot come from the
legacy `2*kill_r` exclusion rule.

The v10.0.5.9 production mesh used a 100 micrometer refined disk because the
legacy production rule is `max(40*h_tip, 50 micrometers)`. The selected J annulus
extended to 240 micrometers, and J decreased as the contour grew. v10.0.5.10 tests
that mesh-support hypothesis directly.

## Scope

This release is audit-only. It does not change:

- elastic constants or anisotropic orientation;
- Arrhenius barriers or attempt frequencies;
- MPZ state, transport, source inventory, or shielding;
- cohesive penalties, topology rules, or crack-event criteria;
- stochastic thresholds or fatigue integration.

The audit installs a mesh constructor through the existing delayed live-binding
stack. It reuses the production radial-ring spacing law, coarsening ratio,
Delaunay triangulation, pre-refined Mode-I corridor, and boundary clipping. The
only mesh change is an explicit 330 micrometer physical refinement radius.

## Acceptance gates

Two grip openings verify elastic load scaling. Contours are evaluated at 100,
140, 180, 240, and 300 micrometers. The publication gate requires:

1. elastic-energy closure within the existing tolerance;
2. J/sigma^2 load scaling within 2 percent;
3. physical refinement radius larger than every accepted contour;
4. peak-to-peak J/sigma^2 variation over 180, 240, and 300 micrometers no larger
   than 10 percent of their median;
5. the selected 240 micrometer J/sigma^2 within 10 percent of the v10.0.5.8
   fixed-grip reference.

Only

```text
production_refinement_support_parity_passed
```

unblocks a later stochastic bracket. Other fail-closed statuses are:

```text
production_refinement_support_inadequate
production_J_contour_instability
production_J_parity_failed_with_adequate_support
production_J_not_quadratic_in_load
production_elastic_energy_closure_failed
```

## New diagnostics

Each probe records:

- the actual installed physical refinement radius and centers;
- radial element count and mean/median/p95/max characteristic size;
- annulus width divided by p95 element size;
- de-duplicated J contours;
- production, tensile-filtered, and no-exclusion J values;
- the complete v10.0.5.9 parity analysis plus the new contour-support gate.

## Worktree

```bash
BASE=/Volumes/Data/Data/Nanopillar_calculation/Arrhenius_FEM_CZM_MPZ_v10_0_5_9_production_j_parity
NEW=/Volumes/Data/Data/Nanopillar_calculation/Arrhenius_FEM_CZM_MPZ_v10_0_5_10_refinement_support
BRANCH=v10.0.5.10-physical-refinement-support-audit

git -C "$BASE" fetch origin \
  "+refs/heads/${BRANCH}:refs/remotes/origin/${BRANCH}"

git -C "$BASE" worktree add \
  -b "$BRANCH" \
  "$NEW" \
  "refs/remotes/origin/${BRANCH}"

cd "$NEW"
python -m pip install -e . --no-deps
```

If the local branch already exists, omit `-b` and use the remote-tracking ref in a
new version-specific folder.

## Validation

```bash
python -m py_compile \
  arrhenius_fracture/physical_refinement_mesh_v100510.py \
  arrhenius_fracture/production_j_refinement_support_v100510.py \
  arrhenius_fracture/mode_i_first_passage_v10_0_5_10_refinement_probe.py \
  run_v10_0_5_10_refinement_support_audit.py

pytest -q \
  tests/test_v10056_kj_audit_bracket.py \
  tests/test_v10056_audited_wrappers.py \
  tests/test_v10058_fixed_grip_elastic_audit.py \
  tests/test_v10058_fixed_grip_audit_mesh.py \
  tests/test_v10059_production_j_parity.py \
  tests/test_v10059_v911_probe_contract.py \
  tests/test_v10059_straight_path_recorder.py \
  tests/test_v100510_refinement_support.py
```

## Run

```bash
REFERENCE=/Volumes/Data/Data/Nanopillar_calculation/Arrhenius_FEM_CZM_MPZ_v10_0_5_8_fixed_grip_elastic_convergence/runs/v10_0_5_8_fixed_grip_elastic_convergence_v2/fixed_grip_elastic_convergence_v10_0_5_8.json

python run_v10_0_5_10_refinement_support_audit.py \
  --reference-json "$REFERENCE" \
  --openings-um "1 2" \
  --contour-outer-um "100 140 180 240 300" \
  --accepted-contour-um "180 240 300" \
  --selected-outer-um 240 \
  --tip-refinement-radius-um 330 \
  --contour-stability-rel-tol 0.10 \
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
  --out runs/v10_0_5_10_refinement_support_v1
```

A final return code of 2 is a valid fail-closed mechanics diagnosis. A missing
probe JSON is a runtime failure and writes a launch-failure manifest with the log
tail.

## Outputs

```text
production_j_refinement_support_v10_0_5_10.json
production_j_refinement_cases_v10_0_5_10.csv
production_j_refinement_contours_v10_0_5_10.csv
production_mesh_radial_support_v10_0_5_10.csv
opening_1um/production_j_refinement_probe_v10_0_5_10.json
opening_1um/production_j_refinement_probe.log
opening_2um/production_j_refinement_probe_v10_0_5_10.json
opening_2um/production_j_refinement_probe.log
```
