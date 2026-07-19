# Arrhenius FEM/CZM moving-process-zone model — v9.0

## Scope

Version 9.0 preserves the full spatial fracture solver and replaces only the front-local scalar crack-tip closure. The production code still contains:

- anisotropic cubic elasticity and crystal rotation;
- crystallographic cleavage and slip-path competition;
- signed directional configurational driving forces;
- multiple independently evolving crack fronts;
- crack deflection and first-passage branch birth;
- parent/daughter state partitioning;
- local and clustered J-integral evaluation;
- adaptive sharp-wake and adaptive cohesive/topological crack backends;
- crack coalescence, branch retirement, and front-network bookkeeping;
- monotonic, cyclic, and mixed-mode loading;
- cyclic mechanical-field updates and spatial diagnostics;
- adaptive physical time/cycle integration;
- snapshots, event histories, and existing postprocessing scripts.

The default remains `legacy_scalar`, so archived v8 commands retain their prior constitutive closure. New work must explicitly select:

```bash
--front-state-model moving_pz
```

## Constitutive replacement

The moving-process-zone (`moving_pz`) engine removes the following quantities as material-fit dimensions:

- scalar fitted shielding coefficient `chi_shield`;
- scalar emission saturation `N_sat`;
- per-numerical-step emission cap `dN_cap`;
- stored-energy subtraction from the cleavage barrier;
- default crack-tip stress ceiling.

It replaces them with a front-owned, moving one-dimensional process-zone state:

1. **Finite source sites.** Each source site emits through the same Arrhenius emission barrier used in all loading protocols. Exact finite-site renewal, `1-exp(-H)`, makes the source update independent of timestep partitioning.
2. **Mobile and retained defect fields.** Each active front carries slip-system-resolved arrays over distance ahead of the current tip.
3. **Activated transport, trapping, detrapping, and recovery.** These are material kinetics, not fatigue-only parameters.
4. **Conservative moving frame.** Crack advance translates the process zone; content left behind is recorded in the wake and virgin material enters ahead.
5. **Direct shielding.** `K_sh` is evaluated from the signed spatial defect distribution with an elastic kernel. No scalar backstress is subtracted from the opening stress.
6. **Slip-based blunting.** Effective tip radius is determined by accumulated local slip rather than retained-event count.
7. **Common loading-protocol parameters.** Monotonic fracture, fatigue, and dwell use the same opening, emission, transport, trapping, and recovery parameters.

## Legacy freeze

The exact source package supplied for the prior first-passage implementation is copied under:

```text
legacy_first_passage_v8/arrhenius_fracture/
```

The legacy implementation is not modified by v9 development. A SHA-256 manifest is included in that folder.

## Important status

The software integration is complete enough for audit, reduced-model fitting, fatigue/dwell tests, and small FEM/CZM validation runs. The four rows in `mpz_four_class_initial_guesses.csv` are **initial guesses only**. They are not a completed calibration and must not be used as publication parameter sets.

The full differential-evolution refit and long FEM/CZM production matrix have intentionally not been run in this package. The intended workflow is to fit and screen in 1-D first, then run only the low/transition/high-temperature FEM/CZM validation matrix.

## Installation

Use the isolated project environment:

```bash
conda activate arrhenius-fem-czm
cd Arrhenius_FEM_CZM_MPZ_v9_0
python -m pip install -e .
```

## Verification

```bash
bash verify_mpz_v9.sh
```

The default verifier compiles the active modules and runs the conservation/integration unit suite. To add the reduced audit, calibration-objective, fatigue, and dwell protocol smokes:

```bash
RUN_PROTOCOL_SMOKES=1 bash verify_mpz_v9.sh
```

To include a minimal anisotropic 2-D CZM solve:

```bash
RUN_2D_SMOKE=1 bash verify_mpz_v9.sh
```

## 1. Audit the frozen model

A small audit:

```bash
python audit_legacy_caps_and_ablations.py \
  --classes "ceramic peak weakT DBTT" \
  --temperatures "300 500 700 900 1100" \
  --dK-values "0.25 0.05 0.02 0.005" \
  --ablations "baseline no_dN_cap no_N_sat no_chi_shield no_stored_energy no_sigma_cap all_removed" \
  --n-advances 20 \
  --out runs/legacy_cap_ablation_audit_v9
```

The audit records the fraction of steps for which the emission cap, scalar saturation, and stress cap are active, together with first and repeated crack events.

## 2. Reduced moving-process-zone fit

First verify the objective and data paths:

```bash
python fit_mpz_four_classes.py \
  --smoke \
  --temperatures "300 500 700 900 1100" \
  --n-advances 20 \
  --out runs/mpz_v9_four_class_fit_smoke
```

Then run the actual refit:

```bash
python fit_mpz_four_classes.py \
  --classes "ceramic peak weakT DBTT" \
  --temperatures "300 400 500 600 700 800 900 1000 1100 1200" \
  --dK 0.05 \
  --Kdot 0.005 \
  --Kmax 65 \
  --n-advances 20 \
  --popsize 16 \
  --maxiter 80 \
  --out runs/mpz_v9_four_class_fit_production
```

The objective combines initiation, repeated-growth class persistence, and state regularization. It does not optimize `chi_shield`, `N_sat`, a per-step cap, or stored-energy cleavage lowering.

Before any spatial validation, inspect:

```text
runs/mpz_v9_four_class_fit_production/mpz_four_class_parameters.csv
runs/mpz_v9_four_class_fit_production/mpz_fit_predictions.json
```

## 3. Unified fatigue and dwell tests

Use the fitted table from the previous step:

```bash
python run_mpz_fatigue_matrix.py \
  --parameters runs/mpz_v9_four_class_fit_production/mpz_four_class_parameters.csv \
  --temperatures "300 700 1100" \
  --Kmax-values "6 8 10 12 14" \
  --R 0.1 \
  --frequency-Hz 1000 \
  --cycles-max 1e10 \
  --require-fitted \
  --out runs/mpz_v9_fatigue_screen
```

```bash
python run_mpz_dwell.py \
  --parameters runs/mpz_v9_four_class_fit_production/mpz_four_class_parameters.csv \
  --temperatures "300 700 1100" \
  --K-MPa-sqrt-m 15 \
  --hold-s 1e5 \
  --require-fitted \
  --out runs/mpz_v9_dwell_screen
```

No separate fatigue or creep material parameter table is introduced.

## 4. Small full-physics FEM/CZM validation

The default validation protocol uses one active front to isolate whether each constitutive response class persists after repeated crack growth. This is a protocol choice, not a simplified solver; the same production executable retains all multifront functionality.

```bash
python run_mpz_fem_czm_validation_matrix.py \
  --parameters runs/mpz_v9_four_class_fit_production/mpz_four_class_parameters.csv \
  --classes "ceramic peak weakT DBTT" \
  --temperatures "300 700 1100" \
  --target-ext-um 100 \
  --max-jobs 1 \
  --require-fitted \
  --out runs/mpz_v9_fem_czm_validation_single_front
```

After the single-front constitutive screen passes, exercise branch birth and coalescence without changing material parameters:

```bash
python run_mpz_fem_czm_validation_matrix.py \
  --parameters runs/mpz_v9_four_class_fit_production/mpz_four_class_parameters.csv \
  --classes "peak DBTT" \
  --temperatures "700 1100" \
  --enable-branching \
  --max-fronts 16 \
  --target-ext-um 100 \
  --max-jobs 1 \
  --require-fitted \
  --out runs/mpz_v9_fem_czm_validation_branching
```

## 5. Production decision gates

Do not launch the prior full 300–1200 K, multi-rate, long-extension campaign until all of the following are true:

- the reduced model reproduces initiation without timestep dependence;
- repeated-growth curves preserve the intended material class;
- shielding remains below the applied driving force and has plausible spatial support;
- source exhaustion is not being used as an implicit threshold;
- results converge with process-zone bins, `dK`, crack increment, time step, and cycle block;
- the small FEM/CZM matrix preserves the class at low, transition, and high temperature;
- fatigue and dwell use the same parameter row and remain numerically stable;
- branch-enabled runs conserve parent-plus-daughter process-zone state.

## Main files added or modified

```text
arrhenius_fracture/moving_process_zone.py
arrhenius_fracture/mpz_front_engine.py
arrhenius_fracture/sharp_front.py
arrhenius_fracture/fatigue_v1.py
arrhenius_fracture/fatigue_sharp_front.py
arrhenius_fracture/mixed_mode_first_passage_v8.py
audit_legacy_caps_and_ablations.py
fit_mpz_four_classes.py
run_mpz_fatigue_matrix.py
run_mpz_dwell.py
run_mpz_fem_czm_validation_matrix.py
mpz_run_utils.py
tests/test_moving_process_zone.py
arrhenius_fracture/sn_intact_fem.py
pytest.ini
```

The supplied source referenced `sn_intact_fem.py` from the stateful local-peridynamics workflow but did not include that module. v9 restores it by factoring the existing intact-FEM algorithms from the retired full-field initiation driver into a shared implementation. The active regression suite completes with 139 tests and 16 subtests.
