# Multi-front v2 branch-propagation correction

This patch addresses the observed production-sweep behavior where branch inventories were born but the rendered branches remained short stubs.

## Root cause

The previous multi-front refactor split a single global crack-extension budget across every front that happened to fire in the same accepted load/time step. That is too restrictive. Energy sharing is required at a local bifurcation event, where one parent tip creates a daughter. It should not be applied globally to unrelated active fronts that each completed their own local renewal clock.

With several active fronts, global sharing gave each front a sub-mesh increment, often much smaller than the damage-kill width and far smaller than the J-integral contour. New daughter branches were therefore seeded as very short stubs; on the next step the branch J contour did not see a meaningful resolved crack and the branch did not continue to propagate.

## Correction

The advance budget is now local-family based:

- An independent fired front gets its own local budget `da_phys * n_fire`.
- If that front branches on the same event, the parent and newborn daughter split that parent budget using the branch hazard weights.
- Other independently fired fronts do not divide the same budget.

This keeps the conservation rule where it belongs, at branch birth, but stops unrelated tips from starving one another numerically.

## Diagnostics/rendering

Field snapshots now overlay the explicit sharp-front polylines. The damage field alone can make short branches look absent because the crack is represented by element/node killing; the polyline overlay shows what the crack inventory actually did.

`branch_diagnostics_*K.csv` now reports `advance1_m` and `advance2_m` from the actual last advance of the primary and leading daughter, instead of zeros.

## Still to check

If branches remain short after this correction, the next issue is likely scale separation between `da_phys` and the J-contour radius. A newborn branch with length much less than the J-integral inner contour is not yet an independently resolved sharp crack. The physically consistent closure would be a process-zone handoff: unresolved daughter branches use the parent-tip kink hazard until they grow beyond a handoff length of order the process-zone/J-contour inner radius. That is not included in this patch.
