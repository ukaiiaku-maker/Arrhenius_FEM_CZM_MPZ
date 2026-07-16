# MPZ v9.18.5.1: startup-safe committed-target horizon

Branch: `v9.18.5.1-safe-target-stop-horizon`

## Failure addressed

The first v9.18.5 ceramic diagnostic passed all unit/regression tests but failed
before creating any physical event or field snapshot. The campaign summary
contained zero committed events, NaN endpoint metrics, and subprocess return
code 1.

The pasted terminal output did not include the per-case traceback from
`matrix_logs/ceramic_tip_only_seed1_700K.log`, so the exact failing initialization
line was not available. The failure boundary nevertheless isolates the new
startup-only feature in v9.18.5: replacement of `args.steps` by a custom
numeric proxy.

The 2-D driver and downstream configuration paths may use strict integer/index
semantics before the first FEM solve. A class that merely implements common
arithmetic operators is not a safe substitute for an actual integer.

## Correction

v9.18.5.1 replaces the proxy with `SafeDynamicStepHorizon`, an `int` subclass.
It therefore behaves as the ordinary requested step count for:

- `range()` and `operator.index()`;
- arithmetic and comparisons;
- NumPy integer/index conversion;
- JSON serialization and formatting;
- assignment to integer configuration fields.

Only the reflected comparison used by

```python
while accepted_step < args.steps:
```

is dynamic. Before the target commit it has ordinary integer behavior. Once the
v9.18.5 controller sets `v9185_stop_requested`, the comparison returns false and
the accepted-step loop exits.

No v9.18.5 geometry, quality, corridor, mechanics, hazard, cohesive-opening,
MPZ, wake, shielding, or analysis behavior is changed.

## Required diagnostic

Use a new worktree and output root, then repeat ceramic at 700 K to 60 um. The
first acceptance criterion is now that the solver reaches its startup banner,
initial K/J evaluation, and first physical event rather than failing with zero
events.

If a subprocess still fails, inspect the actual traceback immediately with:

```bash
LOG=runs/mpz_v9_18_5_1_ceramic_700K_60um_v1/T700K/matrix_logs/ceramic_tip_only_seed1_700K.log
sed -n '1,240p' "$LOG"
tail -n 160 "$LOG"
```

The campaign summary alone does not contain the Python exception.
