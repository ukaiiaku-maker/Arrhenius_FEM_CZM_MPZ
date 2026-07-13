# V5 validation record

Completed in the available source snapshot:

* Python compilation of all v5 modules and scripts.
* Fifteen regression tests.
* Actual anisotropic elastic calibration at crystal theta = 45 degrees.
* Calibration recovered exact Mode-I normalization:
  * cleavage factor = 1.0
  * emission factor = 1.0
* Nine-angle calibration passed from -60 to +60 degrees.
* Directional factors were symmetric and finite; no factor cap activated.
* One-step coupled 2D anisotropic sharp-front smoke test completed through FEM,
  J evaluation, directional partition, calibrated tip-stress kinetics, and
  summary writing.

The available source snapshot predates the active adaptive-CZM backend, so the
full adaptive-CZM production campaign could not be executed in this environment.
That backend-facing path is exercised on the user's active project by the
recommended preflight.
