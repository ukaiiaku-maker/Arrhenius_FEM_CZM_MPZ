# One-dimensional V2 current parameter diagnostic

## Classification

| Class | Current-row classification | Reason |
|---|---|---|
| Peak | QUALITATIVELY_VALID_BUT_PROVIDER_SENSITIVE | Both providers show an intermediate/high-temperature enhancement and later weakening; PF alone adds minor low-temperature subdivision. |
| DBTT | PARAMETER_RESEARCH_REQUIRED | Both providers reproduce a low shelf and upper transition, but high-temperature state growth is excessive and two cases per provider hit the radius-map bound. The 1000 K native topology does not match either own 2-D target. |
| weak-T | QUALITATIVELY_VALID_BUT_PROVIDER_SENSITIVE | Initial onset remains limited in span, but PF has 2–5 avalanches while FEM/CZM has one. |
| ceramic-like | ROBUST_ACROSS_BOTH_PROVIDERS | Low/declining onset and essentially one long avalanche under both providers. |

| material_class | provider | complete_cases | cases | avalanche_min | avalanche_max | onset_min | onset_max |
|---|---|---|---|---|---|---|---|
| DBTT | FEMCZM | 3 | 5 | 1 | 10 | 22.5043 | 206.427 |
| DBTT | PF | 3 | 5 | 1 | 3 | 22.3521 | 107.17 |
| Peak | FEMCZM | 5 | 5 | 1 | 1 | 22.7055 | 36.9371 |
| Peak | PF | 5 | 5 | 1 | 2 | 22.6378 | 57.9284 |
| ceramic-like | FEMCZM | 5 | 5 | 1 | 1 | 9.8214 | 14.1671 |
| ceramic-like | PF | 5 | 5 | 1 | 1 | 10.1574 | 14.005 |
| weak-T | FEMCZM | 5 | 5 | 1 | 1 | 14.3469 | 15.6601 |
| weak-T | PF | 5 | 5 | 2 | 5 | 18.8517 | 25.2876 |

Absolute PF-native and FEM-qualified values were not forced to agree; classification uses trend, topology, and state evolution.
