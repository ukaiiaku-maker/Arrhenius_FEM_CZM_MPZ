# MPZ v9.6 changelog

## Peierls–Taylor correction

Version 9.6 removes the exploratory caps and algebraic saturation functions that generated the artificial rise–plateau–rise flow-stress curve in v9.4/v9.5.

Production Peierls–Taylor kinetics now use:

- the active emission EXP-floor surface as the common parent;
- scaled Peierls and Taylor EXP-floor barriers;
- exact forward-minus-reverse detailed balance;
- the natural forest spacing `delta = 1/(2 sqrt(rho_f))` without `phi_max`;
- an unbounded correlated hit order derived from a physical correlation length;
- the inverse mean gamma waiting time `lambda_T/m` for constant-condition Taylor completion;
- an exact sequential Peierls–Taylor rate;
- explicit mobile density when available, with an unsaturated linear fallback only for the current one-field bulk FEM;
- the natural forest spacing as the event travel length without a minimum jump length.

The production path no longer uses `pt_taylor_m_cap`, `pt_taylor_phi_max`, `pt_mobile_saturation_density_m2`, `pt_mobile_density_floor_m2`, `pt_jump_length_min_m`, the Taylor renewal-time conversion, or a constitutive rate cap. Legacy fields remain parseable for restart compatibility but are ignored by the v9.6 model.

## Search correction

The previous top-five-per-region and strict-common-closure down-selection is no longer the entry gate for DBTT development. The new broad analytical map:

- evaluates the full refined intrinsic atlas;
- adds the exact prior four-class EXP-floor references;
- samples Peierls–Taylor closures independently before asking whether a common family exists;
- uses the prior DBTT `N_sat` and shielding coefficient only as benchmark coordinates for the retained-state scale;
- never applies those benchmark values as production caps;
- scores intrinsic and developed DBTT trends together.

The canonical reference rows are attached to the nearest fully evaluated atlas first-passage curve, and the proxy candidate ID and normalized distance are written explicitly. The proxy is used only because the exact historical reference file stores barrier/state parameters rather than the complete first-passage trajectory.

## Workflow gate

Do not resume the v9.5 continuation or developed-state optimization until:

1. the uncapped PT audit has been reviewed across temperature, strain rate, forest density, and mobile density;
2. the broad DBTT-capacity map has identified viable intrinsic/PT families;
3. those families have been checked against the prior canonical DBTT benchmark;
4. mobile and forest densities are split explicitly in the bulk FEM state update.
