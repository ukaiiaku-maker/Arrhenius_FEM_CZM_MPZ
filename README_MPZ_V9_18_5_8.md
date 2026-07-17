# MPZ v9.18.5.8: state-coupled material differentiation

## Why this branch exists

The v9.18.5.7 weakT/DBTT endpoint campaign remained nearly self-similar after normalizing the R-curves by their initiation magnitude. The cause is structural: the campaign was run with `BULK_PLASTICITY_MODE=tip_only`. In that mode the continuum FEM is elastic. The MPZ/wake state changes the Arrhenius hazard through direct scalar K shielding, but it does not change the equilibrium compliance or the domain-integral J response.

The audited v9.18.5.7 archive showed:

- all four `KJ/Uapp` geometry-factor histories were identical to floating-point precision;
- weakT at 300 and 1100 K differed by less than 0.04% after K normalization;
- direct MPZ/wake shielding was below 0.4% of the applied K in the nonzero cases;
- many events consumed the full source inventory refreshed by the fixed 5 um crack increment, making emission per event nearly extension-controlled rather than hazard-controlled.

Therefore magnitude shifts were being misidentified as material-response differentiation.

## Changes

1. `run_mpz_v9_18_5_8_state_coupled_material_sweep.sh` defaults to and requires `bulk_same_pt_km` for material-comparison campaigns. `tip_only` requires the explicit diagnostic override `ALLOW_ELASTIC_COLLAPSE=1`.
2. The existing independently calibrated class Peierls/Taylor surfaces and Kocks--Mecking storage/recovery law are used. No barrier rescaling, artificial backstress multiplier, reload gate, shielding multiplier, or emission cap is added.
3. `audit_v91858_state_coupled_differentiation.py` requires:
   - explicit mobile/retained bulk state;
   - nonzero bulk-state update calls;
   - nonzero accepted accumulated plastic strain;
   - target-complete and quality-certified cases;
   - a non-collapsed normalized R-curve or geometry-factor history for every same-temperature material pair.
4. The audit reports source-inventory-limited event fractions rather than silently treating constant refreshed-site consumption as temperature-dependent kinetics.

## Promotion rule

A run is not promoted merely because the classes have different initiation K values. It must pass both:

- v9.18.5.7 subsegment-aware geometry certification;
- v9.18.5.8 normalized material-shape differentiation.

Default collapse thresholds are 2% maximum normalized-K separation and 1% maximum relative geometry-factor separation. A pair is rejected only when both measures remain below their thresholds.

## Scope

This branch activates an existing state-coupled constitutive path; it does not claim that the bulk scalar mobile/retained transfer is the final slip-system-resolved formulation. The first run should remain a short weakT/DBTT 700 K gate before any temperature sweep.
