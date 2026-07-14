# MPZ v9.2 changelog

## Analytical intrinsic first-passage atlas

- Replaced the numerical first-passage search as the primary intrinsic-screening stage with direct analytical integration of the EXP-floor cleavage and emission hazards.
- Added a six-dimensional Sobol atlas over the intrinsic cleavage/emission barrier parameters for each ceramic, weak-T, and DBTT shape family.
- Added explicit region labels for ceramic intrinsic behavior, weak-temperature intrinsic behavior, DBTT precursors, emission-saturated cases, and unclassified mixed cases.
- Added fine-increment re-evaluation of shortlisted candidates.
- Added complete material-row exports for the next steady/transient moving-process-zone stage.
- Added analytical evaluation of the v9.1 initial and completed Stage-1 rows as regression anchors.
- Preserved the attached Panel-A workflow's raw-data/metrics/replot pattern while removing its legacy scalar `chi_shield` and `N_sat` continuation path from the active calibration.

## Interpretation

A `DBTT_precursor` indicates a low-temperature cleavage-dominated tip and a strong high-temperature increase in emission exposure. It is not promoted to a DBTT material until the moving-process-zone model demonstrates transport, trapping, retention, shielding/blunting, and convergence toward the developed state over a physically relevant crack extension.
