# MPZ v9.1 validation report

## Regression

```text
142 passed, 16 subtests passed
```

## Initial-guess smoke

All three initial rows completed 41 crack events at a 25 micrometre increment, spanning 0–1000 micrometres at 300, 700, and 1100 K. This verifies the repeated-growth integration and robust metric outputs.

The inherited guesses are not acceptable final parameterizations:

- ceramic already has the intended decreasing-temperature topology and a small R-curve increment, but its absolute targets still require fitting;
- weakT decreases strongly with temperature and develops only about 1 MPa sqrt(m) of resistance, rather than the required temperature-independent moderate R-curve;
- DBTT decreases with temperature and therefore fails both the first-passage and developed-resistance transition constraints.

This is the intended starting condition: the new optimization must generate the response through the moving process zone rather than inherit the old scalar fit.

## Optimizer-path smoke

A one-generation differential-evolution run was completed for all three classes. This test verifies parameter transforms, negative thermal-slope serialization, checkpointing, class-specific objectives, and output generation. The resulting rows are not retained as calibrated parameters.

## Numerical safeguards

- A completed renewal is resolved by adaptive subdivision; `max_advances_per_step=1` is an event-resolution mechanism, not a cap on physical event count.
- Residual crack-opening clock is retained exactly.
- Finite-site emission uses exact renewal probability and is invariant to timestep partitioning.
- Event-level shielding constraints use pre-renewal state and do not divide a transient finite shielding field by a vanishing applied K.
- Incomplete 1000 micrometre histories are penalized rather than interpreted as brittle response.
