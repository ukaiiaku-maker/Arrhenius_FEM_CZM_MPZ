# One-dimensional V2 current-parameter diagnostic

The historical rows are reference states, not certified material constants. Correcting model form supersedes the e2ed84c search conclusion: Peak and DBTT remain useful, while weak-T and ceramic-like require shared-row replacement.

| material_class | current_classification | decision |
| --- | --- | --- |
| Peak | ROBUST_ACROSS_BOTH_PROVIDERS | retain historical row |
| DBTT | QUALITATIVELY_VALID_BUT_PROVIDER_SENSITIVE | retain historical row; low-T PF grouping mismatch |
| weak-T | NEW_PARAMETER_SEARCH_REQUIRED | historical FEM/CZM had six avalanches and a large re-initiation envelope |
| ceramic-like | NEW_PARAMETER_SEARCH_REQUIRED | historical FEM/CZM had three avalanches and excessive onset |

At 900 K the historical weak-T FEM/CZM lane produced six avalanches and the historical ceramic lane produced three; the final shared rows reduce these to two and one respectively for the canonical seed without changing backend lifecycle parameters.
