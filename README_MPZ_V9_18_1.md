# MPZ v9.18.1: active-event renewal rollback

Branch: `v9.18.1-active-event-renewal-rollback`

## Failure reproduced from v9.18

A 700 K ceramic run inserted nine 5 um topology segments (45 um raw topology
extension) before the inherited emission-only continuation path attempted to
defer another cleavage renewal while a cohesive event was already active:

```text
RuntimeError: cannot defer a second renewal while a trial event is active
```

Because the exception occurred before the event payload was written, the archive
does not establish that every inserted segment had reached physical cohesive
commit. The raw topology count is therefore not used as committed extension.

The continuation path sets the instantaneous cleavage prefactor and action to
zero, but the first trapezoidal hazard update can still inherit the previous
loading-step value of `_lambda_c_prev`. That interpolated half-step can cross a
threshold even though cleavage is supposed to be disabled during the
emission-only continuation.

## Correction

`RenewalRollbackPersistentWakeController` treats a threshold crossing during an
already active event as a numerical continuation artifact. It restores the exact
pre-renewal:

- MPZ and persistent-wake state;
- cleavage action;
- stochastic/deterministic threshold-stream state;
- crack-advance counters;
- reload state.

The current cohesive event remains active. A normal one-fire renewal when no
trial event is active follows the unchanged v9.18/v9.17 path.

The payload reports:

```text
active_event_renewal_transactional_rollback_enabled
active_event_renewals_rolled_back
active_event_thresholds_rolled_back
```

## Temperature syntax

`T_K` is scalar. Do not use:

```bash
T_K=300 700 1100
```

For a sweep, use:

```bash
TEMPS="300 700 1100" \
TARGET_EXT_UM=10 \
OUTROOT_BASE=runs/mpz_v9_18_1_persistent_wake_3T_10um_v1 \
bash run_mpz_v9_18_1_persistent_plastic_wake_sweep.sh
```

Run the 10 um gate before increasing to 100 or 500 um.
