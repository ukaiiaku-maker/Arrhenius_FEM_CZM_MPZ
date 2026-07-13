# Validation performed

- Python compilation: passed for all new modules/scripts.
- Shell syntax: passed for `run_mixed_mode_fem_czm_v1_campaign.sh`.
- Unit tests: 4 passed.
- Elastic calibration smoke test on a coarse mesh:
  - targets: `-30, 0, +30 deg`;
  - achieved errors: approximately `+0.19, +0.01, -0.14 deg`.
- Full sharp-front wrapper smoke test:
  - completed two mixed-boundary FEM/J steps;
  - wrote mixed-mode projection records and summary;
  - deterministic threshold audit written correctly.

The full adaptive-CZM campaign was not run in this build environment because the uploaded source snapshot does not contain the newer adaptive-CZM backend used by the user's active `Arrhenius_FEM_CZM` project. The verification script checks that backend-facing CLI options exist before a campaign is launched.
