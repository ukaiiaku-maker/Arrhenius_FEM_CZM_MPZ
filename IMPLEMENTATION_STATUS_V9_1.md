# MPZ v9.1 implementation status

## Complete

- Three-class target system: ceramic, weakT/FCC-like, and DBTT.
- Peak class removed from the new fitting workflow.
- First-passage and 0–1000 micrometre repeated-growth simulation.
- Adaptive crack-clock integration without a physical event cap.
- Robust early-slope, late-window median, terminal-slope, and saturation-fit metrics.
- Pre-renewal process-zone diagnostics at every crack event.
- Staged first-passage, R-curve, and joint optimization.
- Checkpoint/restart and completed-class skipping.
- DBTT constraints requiring a microstructural low/high-temperature switch.
- Shared fitted parameter table for monotonic, fatigue, dwell, and FEM/CZM protocols.
- Convergence audit over loading increment, MPZ bin count, and crack increment.
- Full architecture preservation from v9.0.

## Deliberately not complete

- No production three-class parameter fit is claimed in this package.
- The CSV target values are synthetic design anchors and remain editable.
- No long FEM/CZM production campaign has been run with v9.1 parameterizations.
- No separate fatigue or creep parameter table has been introduced.

## Validation completed in the package build

- Initial-guess three-class smoke completed all requested 1000 micrometre histories.
- Differential-evolution optimizer entry was exercised for all three classes.
- Full regression result: 142 tests and 16 subtests passed.
