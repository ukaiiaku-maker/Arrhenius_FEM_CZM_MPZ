# MPZ v9.18.2: committed-event completion handshake

Branch: `v9.18.2-committed-completion-handshake`

## Scope

This branch changes campaign metadata only. It does not change:

- absolute cleavage-hazard integration;
- source refresh;
- Peierls--Taylor transport/retention;
- persistent wake remapping;
- wake shielding;
- cohesive opening;
- crack-advance distance.

## Why v9.18.1 reported `rc=1`

The inherited v9.16 protocol requests one additional topology quantum beyond the
analysis target. For a 10 um analysis with 5 um events, the inner topology target
is 15 um. This guard prevents the solver from stopping immediately after the
second trial interface is inserted but before it completes physical opening.

v9.18 correctly disables further renewals once 10 um has been physically
committed. The inner legacy driver therefore sees only 10 of its requested 15 um
of raw topology and exits `right_censored`, even though:

- the inner FEM solver return code is zero;
- the v9.18 committed-event audit reaches 10 um;
- no trial event remains active;
- persistent-wake conservation passes.

## Promotion rule

A topology-guard exit is promoted to successful completion only when all are true:

1. the inner FEM solver return code is zero;
2. `v918_target_committed` is true;
3. `v918_no_uncommitted_trial_at_exit` is true;
4. `v918_persistent_wake_commit_gate_passed` is true.

The old topology target and return code are retained as provenance in:

- `legacy_topology_guard_target_extension_um`;
- `legacy_topology_guard_returncode`.

The analysis target remains the physically committed target.

## Repair the completed 700 K v9.18.1 gate

No FEM rerun is required:

```bash
python repair_v9181_committed_completion.py \
  --campaign-root runs/mpz_v9_18_1_persistent_wake_700K_10um_v1/T700K \
  --T-K 700 \
  --target-extension-um 10 \
  --seeds "1" \
  --classes "ceramic weakT DBTT"
```

## Future sweeps

```bash
TEMPS="300 700 1100" \
TARGET_EXT_UM=10 \
OUTROOT_BASE=runs/mpz_v9_18_2_persistent_wake_3T_10um_v1 \
bash run_mpz_v9_18_2_persistent_plastic_wake_sweep.sh
```

The 10 um result contains only two physical events and is a software/physics
transfer gate, not a publication R-curve. Identical straight crack paths and
strongly similar two-point normalized shapes are expected and do not establish
that the material parameterizations are equivalent.

## 700 K physical result from the completed gate

- Ceramic: wake state is conserved but retained population and shielding are
  effectively zero.
- weakT: retained wake persists into event 2 with approximately 0.023--0.024
  MPa sqrt(m) wake shielding.
- DBTT: retained wake persists into event 2 with approximately 0.0015--0.0051
  MPa sqrt(m) wake shielding.

These wake contributions are small relative to initiation toughness. A longer
run and a matched `ARRHENIUS_WAKE_SHIELDING=0` ablation are required before
claiming that the wake materially raises the R-curve.
