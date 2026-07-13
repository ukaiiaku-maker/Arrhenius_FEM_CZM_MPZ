# Validation report — MPZ v9.0

Validation was performed against the active v9 tree, not the frozen legacy copy.

## Automated regression

```text
PYTHONPATH=. pytest -q
139 passed, 16 subtests passed
```

The suite includes the moving-process-zone conservation/invariance tests, mixed-mode wrapper tests, adaptive CZM and multifront infrastructure tests, and the stateful local-peridynamics core tests.

## Package verifier

```text
python -m compileall -q arrhenius_fracture *.py
bash verify_mpz_v9.sh
8 passed
package version: 0.9.0
MPZ v9 verification passed
```

## Entry-point checks

The following new drivers successfully parsed their complete command-line interfaces with `--help`:

- `audit_legacy_caps_and_ablations.py`
- `fit_mpz_four_classes.py`
- `run_mpz_fatigue_matrix.py`
- `run_mpz_dwell.py`
- `run_mpz_fem_czm_validation_matrix.py`

Separate smoke invocations completed for the 1-D MPZ front, unified fatigue, dwell, anisotropic 2-D adaptive CZM, branch-enabled anisotropic 2-D entry, and mixed-mode v8 right-censored entry.

## Scope boundary

These checks establish implementation integrity and feature compatibility. They do **not** constitute production calibration of the four material classes. `mpz_four_class_initial_guesses.csv` remains explicitly uncalibrated, and publication-level FEM/CZM, fatigue, and dwell campaigns must wait for the repeated-growth fit and convergence gates described in `README_MPZ_V9_0.md`.
