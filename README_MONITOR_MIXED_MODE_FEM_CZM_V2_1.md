# Mixed-mode FEM/CZM v2.1 read-only monitor

Copy both monitor files into the top-level `Arrhenius_FEM_CZM` directory.
The monitor reads campaign output files only; it never imports or modifies the solver.

## Live use

```bash
OUTROOT=runs/mixed_mode_fem_czm_v2_event_controlled_500K \
CLASSES="ceramic DBTT" \
TARGET_PSI="-60 -45 -30 -15 0 15 30 45 60" \
bash monitor_mixed_mode_fem_czm_v2_1.sh
```

Press Ctrl-C to stop the monitor. The campaign continues running.

## One snapshot

```bash
OUTROOT=runs/mixed_mode_fem_czm_v2_event_controlled_500K \
bash monitor_mixed_mode_fem_czm_v2_1.sh --once
```

## Show recent log lines

```bash
SHOW_ACTIVE_LOG_LINES=5 \
bash monitor_mixed_mode_fem_czm_v2_1.sh
```

## Status meanings

- `PENDING`: no trial directory exists yet.
- `STARTING`: trial directory exists but its log has not started.
- `RUNNING`: the latest trial log was updated recently.
- `STALE?`: no log update within `STALE_MINUTES`; inspect before assuming failure.
- `TRIAL_DONE`: the deterministic solve finished and the controller is selecting the next angle.
- `CONVERGED`: final event-state phase-angle error is within tolerance and the projection is reliable.
- `NOT_CONV`: the controller finalized the best available result without meeting all acceptance criteria.
- `FAILED`: the campaign status file reports a failed case.

The `fit` column is the robust mode-projection reliability flag. `sig` is the most recent printed local tip-stress scale in GPa, not the nominal applied stress.
