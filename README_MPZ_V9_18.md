# MPZ v9.18: persistent signed plastic wake

Branch: `v9.18-persistent-plastic-wake`

## Physical scope

v9.18 does not assign toughness to a wake merely because plastic state exists
behind the crack.  The only new crack-driving contribution is the elastic
Mode-I-equivalent field of the signed retained/mobile dislocation population.
There is no bridging traction, transformation strain, empirical wake toughness,
or extra fracture-energy term.

Accumulated wake slip is conserved and exported for later eigenstrain/FEM work,
but it does not enter the current-tip blunting radius in this branch.

## Changes from v9.17.1

1. Active MPZ state crossed by a committed crack increment is remapped into
   finite-volume bins behind the current tip.
2. Existing wake bins translate farther behind the tip as the crack advances.
3. Wake retained/mobile populations continue Peierls--Taylor exchange, recovery,
   annihilation, and mobile transport at the unloaded wake stress.
4. Total shielding is the sum of active-zone and persistent-wake signed line
   fields. `ARRHENIUS_WAKE_SHIELDING=0` provides an ablation.
5. Event pause/resume uses `ARRHENIUS_NOMINAL_LOADING_DT_S`, not the adaptive
   event-localization fraction.
6. Event audits distinguish precommit active retention, postcommit active state,
   persistent wake state, discarded state, wake shielding, and carryover to the
   next event.
7. Once the committed target is reached, the active front is set to zero allowed
   renewals so no third trial interface is created.

## Important limitation

The inherited `sharp_front.run_2d` stopping condition compares the requested
extension with raw topology extension. v9.18 suppresses post-target renewals and
therefore prevents an extra topology event, but it does not yet add a generic
committed-extension callback to that shared solver loop. A short gate may remain
in an idle state until its configured step limit after the committed target.
The audit remains authoritative for completion.

## New files

- `arrhenius_fracture/moving_process_zone_v918.py`
- `arrhenius_fracture/mode_i_first_passage_v9_18.py`
- `arrhenius_fracture/coupled_event_audit_v918.py`
- `run_mpz_v9_18_mode_i_rcurve.py`
- `run_mpz_v9_18_persistent_plastic_wake.py`
- `run_mpz_v9_18_persistent_plastic_wake_700K.sh`
- `tests/test_persistent_wake_v918.py`
- `tests/test_persistent_wake_audit_v918.py`

## Required gate

```bash
CONDA_ENV=arrhenius-fem-czm \
SEEDS="1" \
CLASSES="ceramic weakT DBTT" \
T_K=700 \
TARGET_EXT_UM=10 \
STEPS=10000 \
NX=36 NY=72 \
TIP_H_FINE=1e-6 TIP_RATIO=1.20 \
DU=2e-7 DT=8.4 \
MPZ_LENGTH_UM=100 MPZ_N_BINS=200 \
WAKE_LENGTH_UM=100 WAKE_N_BINS=0 \
WAKE_SHIELDING=1 WAKE_SHIELD_PROJECTION=1 \
BULK_PLASTICITY_MODE=tip_only \
EVENT_TARGET_DQ=0.05 \
EVENT_MIN_DT_S=1e-12 \
EVENT_MAX_FIXED_HOLD_S=inf \
OUTROOT=runs/mpz_v9_18_persistent_plastic_wake_700K_10um_v1 \
AUDIT_OUT=runs/mpz_v9_18_matched_stress_audit_700K_v1 \
bash run_mpz_v9_18_persistent_plastic_wake_700K.sh
```

## Interpretation gate

A mechanically meaningful persistent wake requires all of the following:

- retained-state conservation on commit;
- nonzero `wake_retained_postcommit`;
- nonzero `wake_K_shield_postcommit_Pa_sqrt_m`;
- nonzero wake retention at the next event nucleation;
- class-dependent propagation response that survives the
  `ARRHENIUS_WAKE_SHIELDING=0` ablation.

If the state is conserved but the wake shielding remains negligible, the correct
conclusion is that wake plasticity is not a meaningful toughness mechanism for
that parameterization. The code must not compensate by inserting an arbitrary
wake-toughness term.
