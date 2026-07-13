# Mixed-mode FEM/CZM v3.1: J-consistent amplitude with phase-ratio calibration

## Why v3 failed

v3 correctly moved the authoritative mixed-mode amplitude to the domain J integral, but its elastic calibration still required the full Williams amplitude fit to satisfy `projection_rel_rmse <= 0.40` and annular K-magnitude stability. On the active blunt-notch mesh, those amplitude diagnostics fail even for nominal Mode I, so every target was rejected before a fracture run began.

That was the wrong gate. In v3/v3.1, the Williams projection is used only to determine the elastic phase ratio `KII/KI`; the fracture amplitude comes from `KJ`. A poor absolute amplitude residual should therefore be reported, not used to reject phase calibration.

## v3.1 correction

1. Solve two elastic basis cases: unit opening and unit sliding.
2. Construct the measured 2x2 response matrix mapping boundary displacement components to projected `[KI,KII]`.
3. Solve this matrix for each requested phase angle, retaining geometry-induced cross-coupling.
4. Verify the requested phase with an actual elastic solve and refine the boundary angle if needed.
5. Gate calibration only on finite phase ratio, point count, fit count, conditioning, and annulus-to-annulus phase spread.
6. Retain amplitude residual and K-magnitude spread as diagnostics only.
7. Continue to use the domain J integral as the authoritative amplitude during fracture calculations.

The v3 files remain untouched. All v3.1 files have unique names and output folders.
