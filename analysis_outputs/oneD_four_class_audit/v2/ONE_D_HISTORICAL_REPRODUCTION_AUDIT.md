# One D Historical Reproduction Audit

## Result

The unchanged four-class, five-temperature, 100-µm matrix completed: 20/20 cases at 300, 600, 900, 1100, and 1200 K. It used the archived 23-event loading map (111.53753370116257 µm coverage), seed 3621, nominal `dU=2e-7 m`, nominal `dt=8.4 s`, and translation exponent 0.95.

Weak-T and ceramic legacy developed endpoints reproduce the archived selection values to floating-point serialization. Peak and DBTT were historically selected through different temperature grids/50-µm screening followed by 2-D transfer; no matched archived five-temperature/100-µm score exists. Their deterministic rerun is complete, but a bitwise historical score claim would be false, so those cells are `NO_ARCHIVED_MATCHED_PROTOCOL`.

The V2 analysis is derived from completed ledgers and does not advance RNG. Original ranks/scores are retained only where the selection manifest actually provides them.
