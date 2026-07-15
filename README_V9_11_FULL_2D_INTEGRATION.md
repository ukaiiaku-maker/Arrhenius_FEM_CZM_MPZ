# FEM/CZM–MPZ v9.11 full 2-D integration

This branch integrates the selected v9.10.2 ceramic/weakT and v9.10.3 DBTT parameterizations into the existing `Arrhenius_FEM_CZM_MPZ` solver. It does not replace the FEM/CZM mechanics with a reduced surrogate.

## Preserved solver infrastructure

The implementation calls the existing 2-D elastic–plastic FEM, domain/interaction integral, anisotropic mixed-mode directionality, adaptive CZM crack backend, adaptive event stepping, process-zone history, remeshing, fatigue controller, and restart machinery. Existing v8/v9.10 modules are not overwritten.

The initial v9.11 validation runner intentionally uses one active front with branching disabled. This is a validation gate, not removal of the underlying branch/coalescence code. Multi-front profile ownership must be validated separately before branching is enabled in v9.11 production.

## Constitutive changes

- Exact selected parameter manifests are the only source of ceramic, weakT, and DBTT parameters.
- Cleavage, emission, Peierls, and Taylor retain independent EXP-floor shapes.
- The finite source-site inventory and crack-advance refresh used during spatial promotion are retained.
- Peierls controls mobile transport; geometric forest encounters retain mobile lines; correlated Taylor completion releases them.
- The natural Taylor hit order is uncapped.
- Mobile saturation, mobile density floor, jump-length floor, constitutive rate cap, scalar `N_sat`, emission increment cap, stored-energy cleavage lowering, and unsigned-density backstress are inactive.
- The FEM-derived finite-radius stress supplies only a dimensionless process-zone shape. The calibrated sharp-tip magnitude remains `K/sqrt(2*pi*r_eff)`.
- Bulk scalar density supplies the Taylor forest-density baseline and is never interpreted as signed shielding.
- Bulk plastic redistribution is already included in the FEM/domain-integral drive. Only unresolved retained-line `K_shield` is subtracted explicitly, once.

## Install

```bash
cd /Volumes/Data/Data/Nanopillar_calculation/Arrhenius_FEM_CZM_MPZ

git fetch origin \
  refs/heads/v9.11-full-2d-three-class-integration:refs/remotes/origin/v9.11-full-2d-three-class-integration

git switch -c v9.11-full-2d-three-class-integration \
  --track origin/v9.11-full-2d-three-class-integration

conda activate arrhenius-fem-czm
python -m pip install -e . --no-deps
```

## Preflight

```bash
python verify_mpz_v9_11_install.py .
python verify_mpz_v9_11_physics.py --parameter-root mpz_v9_11_parameters
python -m pytest -q \
  tests/test_mpz_v9_11_2d_coupling.py \
  tests/test_mpz_v9_10_2_independent_shapes.py \
  tests/test_bulk_pt_plasticity.py
```

## First full 2-D smoke

The smoke runner requires a verified v8 production-backend calibration CSV containing the `psi=0` row.

```bash
CALIBRATION_CSV=runs/<verified_v8_calibration>/mixed_mode_loading_calibration_v8.csv \
CLASS=DBTT \
T_K=700 \
OUTROOT=runs/mpz_v9_11_full2d_DBTT_700K_smoke_v1 \
bash run_mpz_v9_11_full2d_smoke.sh
```

Set `RUN_SOLVER=0` to run only the installation, physics, and unit-test preflight.

## Four-temperature gate

```bash
CALIBRATION_CSV=runs/<verified_v8_calibration>/mixed_mode_loading_calibration_v8.csv \
CLASSES="ceramic weakT DBTT" \
TEMPS="300 700 900 1200" \
ROOT=runs/mpz_v9_11_modeI_temperature_gate_v1 \
bash run_mpz_v9_11_modeI_temperature_gate.sh
```

No parameter fitting is permitted during these gates. Evaluate first-passage toughness, source availability, mobile/retained populations, direct unresolved `K_shield`, 2-D profile coverage, bulk plastic work, and J-contour stability before long R-curve, fatigue, or branching runs.
