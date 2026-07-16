# MPZ v9.16: kinetic trial-cohesive opening

This branch replaces the prescribed v9.15 cohesive smoothstep with a kinetic
trial-interface protocol while retaining Arrhenius first passage as the sole
fracture nucleation criterion.

## Constitutive changes

1. A completed cleavage clock inserts an initially intact trial CZM segment.
2. The MPZ moving-frame renewal is consumed but its spatial translation is
   deferred.
3. At fixed remote displacement, cohesive progress obeys

   `dq/dt = (lambda_c / lambda_c,nucleation)^p / tau_event`.

4. Emission, transport, shielding, FEM equilibrium, and cohesive opening evolve
   over the same physical substeps.
5. The MPZ renewal and wake transfer are committed only when `q -> 1`.
6. If the cleavage-rate ratio remains below the arrest threshold, the trial
   interface remains partially intact and loading resumes. The event can restart
   when the rate ratio exceeds the resume threshold.

The current solver still inserts the trial topology before the relaxation loop.
Accordingly, v9.16 distinguishes **topology extension** from **committed physical
extension** in its summaries. The authoritative v9.16 extension is the sum of
committed kinetic events, not the raw topology path.

## New outputs

Each case writes:

- `kinetic_trial_event_relaxation_v916.json`
- `kinetic_trial_event_audit_v916.json`
- `analysis_window_v916.json`
- `v9_16_case_summary.csv/json`

The event history includes instantaneous cleavage-rate ratios, damage
increments, MPZ populations, shielding, cohesive opening/traction, a cohesive
softening-work estimate, and tip-emission work.

## Initial validation

Run the short three-class gate first:

```bash
CONDA_ENV=arrhenius-fem-czm \
SEEDS="1" \
CLASSES="ceramic weakT DBTT" \
T_K=700 \
TARGET_EXT_UM=10 \
BULK_PLASTICITY_MODE=tip_only \
OUTROOT=runs/mpz_v9_16_kinetic_trial_opening_700K_three_class_10um_v1 \
bash run_mpz_v9_16_kinetic_trial_opening_700K.sh
```

Then repeat with explicit bulk plasticity:

```bash
CONDA_ENV=arrhenius-fem-czm \
SEEDS="1" \
CLASSES="ceramic weakT DBTT" \
T_K=700 \
TARGET_EXT_UM=10 \
BULK_PLASTICITY_MODE=bulk_same_pt_km \
OUTROOT=runs/mpz_v9_16_kinetic_trial_opening_700K_three_class_10um_bulk_v1 \
bash run_mpz_v9_16_kinetic_trial_opening_700K.sh
```

Do not start a long-growth campaign until each class has
`all_trial_commit_requirements_passed=true`. The separate field
`physical_coupling_relevance_observed` determines whether event-window
plasticity is large enough to affect the propagation response rather than being
merely present in the software path.
