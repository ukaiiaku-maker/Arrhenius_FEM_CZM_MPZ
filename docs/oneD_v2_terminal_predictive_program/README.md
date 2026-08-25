# oneD V2 terminal predictive program: results record

This directory records the completed one-dimensional V2 predictive-program result generated on the local branch:

- local branch: `codex/oneD-v2-terminal-predictive-program`
- reported local HEAD: `5045fee`
- producer-code commit recorded by the output manifest: `cddc51605aee93d8ecbaa0dda76c12085e08f9de`
- diagnostic baseline: `e2ed84c`

## Important repository note

This GitHub record contains the final reports and compact machine-readable outputs. It does **not** reconstruct or replace the unpublished local implementation commit history. The local Codex worktree must still push the full source branch if the implementation code is to be reviewed or merged.

## Final result

The PF-consistent and FEM/CZM-consistent reduced models are fit for screening over the declared domain. At 1000 K, onset errors were +4.2% (PF Peak), +13.1% (PF DBTT), -7.3% (FEM/CZM Peak), and -9.1% (FEM/CZM DBTT), with matching physical-avalanche topology for all four native baselines.

A 972-case sensitivity study produced 84.8% overall cross-provider onset-sign agreement and agreement for the dominant parameter directions. The final shared screening rows are:

- Peak: `v913_zeroD_sobol_0242980`
- DBTT: `v913_zeroD_sobol_0202500`
- weak-T: `oneD_v2_focused_weak_T_0016`
- ceramic-like: `oneD_v2_focused_ceramic_like_0018`

Six new bounded two-dimensional PF transfer cases validated the new weak-T and ceramic-like rows at 300, 1000, and 1200 K. No new two-dimensional FEM/CZM run was launched, no canonical production trajectory was changed, and no production physical formula was altered.

## Intended use

These rows are **V2 reduced-model screening parameterizations**, not uniquely calibrated material constants. Use them for trend exploration, parameter screening, response-class studies, and selection of higher-fidelity validation cases.

## Files

- `docs/oneD_v2_terminal_predictive_program/`: final scientific reports
- `data/oneD_v2_terminal_predictive_program/`: final registry, baseline results, domain, PF transfer results, decision, and provenance
- `docs/oneD_v2_terminal_predictive_program/ONE_D_V2_TO_FATIGUE_INTEGRATION_HANDOFF.md`: integration contract for cyclic fatigue code
