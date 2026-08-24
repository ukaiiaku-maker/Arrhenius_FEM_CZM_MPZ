# One-dimensional V2 PF normal-path qualification V2

Classification: **PF_NORMAL_PATH_QUALIFIED_FAIL_CLOSED_ON_VETO**.

The controlled Peak fixture used the actual `v10.2.30` PF entry, canonical
seed 8666, the production adaptive predictor/retry loop, and stopped after one
accepted 5 µm event. Tightening the adaptive target to 0.01 forced
**86 rejected trials**. All rejected-trial state
fingerprints exactly matched their transaction snapshots. Diagnostics OFF and
ON produced byte-identical authoritative physical artifacts.

Exactly-once results: first-passage accept = **1**,
geometry commit = **1**, post-event renewal =
**1**.

This source does not draw stochastic cleavage or emission thresholds. Cleavage
uses a deterministic unit-action renewal surface; emission is a continuous
rate. Event length is fixed by `f.da`, and direction is selected deterministically
from the cleavage planes. Their RNG/draw/restore requirements are therefore
**NOT_APPLICABLE**, not inferred from later rows.

Late-veto termination remains **SOURCE_CONTRACT_QUALIFIED** and rollback-and-
continue remains **UNSUPPORTED_BY_PRODUCTION**. This normal-path result does not
alter either status.
