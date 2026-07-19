# v10.0.5.8 — fixed-grip elastic FEM convergence artifact

This point release adds the independent elastic reference that was intentionally
absent from the v10.0.5.6 remote-stress / `K_J` audit.

It does **not** change Arrhenius barriers, process-zone kinetics, cohesive damage,
crack-renewal clocks, stochastic thresholds, plasticity, or production crack
growth. The artifact is a fail-closed mechanics benchmark.

## Why this is required

The v10.0.5.6 audit compared the domain-integral result with the standard
finite-width sharp single-edge-crack relation

```text
K = sigma_gross sqrt(pi a) Y_edge(a/W).
```

The modeled starter defect is instead a finite killed band, and the body is loaded
by symmetric prescribed displacement. A mismatch with the sharp-crack formula can
therefore come from the actual notch geometry and fixed-grip compliance, from the
J implementation, or from both. The sharp-edge formula alone cannot distinguish
those possibilities.

v10.0.5.8 computes the geometry-specific reference directly from elastic FEM:

```text
G_fixed_grip = -d U_el / da at fixed total grip opening.
```

The derivative is evaluated from three separate equilibrated elastic states with
nearby killed-notch lengths on the same mesh. No hazard or crack-growth update is
executed.

## Independent checks

The benchmark requires:

1. stored elastic energy from `sum(psi_stored A)` to agree with `0.5 u^T K u`;
2. the fixed-grip energy derivative to converge with crack-length increment;
3. the energy derivative to converge with tip-mesh refinement;
4. a stable domain-J contour plateau using the **full stored elastic energy**;
5. the plateau median to agree with `G_fixed_grip`.

The current tensile-filtered energy J is also written as an ablation. It is not the
reference used to pass the gate.

For anisotropic elasticity, the primary comparison is the energy ratio
`J_full/G_fixed_grip`. Any `K` and geometry-factor values are isotropic-equivalent
common-unit displays and are not an anisotropic Stroh conversion.

## Installation

```bash
BRANCH=v10.0.5.8-fixed-grip-elastic-fem-convergence

git fetch origin "$BRANCH"
git switch --track "origin/$BRANCH"
python -m pip install -e . --no-deps
```

## Validation

```bash
python -m py_compile \
  arrhenius_fracture/fixed_grip_elastic_audit_v10058.py \
  run_v10_0_5_8_fixed_grip_elastic_convergence.py

pytest -q \
  tests/test_v10056_kj_audit_bracket.py \
  tests/test_v10056_audited_wrappers.py \
  tests/test_v10058_fixed_grip_elastic_audit.py
```

## Production audit

```bash
python run_v10_0_5_8_fixed_grip_elastic_convergence.py \
  --tip-h-um "10 5 2.5" \
  --crack-increment-um "40 20 10" \
  --contour-outer-um "100 140 180 240 300" \
  --width-mm 2 \
  --height-mm 4 \
  --crack-mm 0.5 \
  --notch-half-thickness-um 80 \
  --grip-opening-um 2 \
  --anisotropic \
  --crystal-theta-deg 45 \
  --out runs/v10_0_5_8_fixed_grip_elastic_convergence_v1
```

The crack increments are intentionally no smaller than the coarsest 10 µm tip
spacing. If a requested perturbation does not move the actual killed-node tip, the
audit fails closed instead of fabricating a zero-width derivative.

The process returns zero only when every convergence and J/energy agreement gate
passes. A nonzero return is an audit result, not a solver crash.

## Outputs

```text
fixed_grip_energy_release_v10_0_5_8.csv
fixed_grip_J_contours_v10_0_5_8.csv
fixed_grip_elastic_convergence_v10_0_5_8.json
fixed_grip_elastic_convergence_v10_0_5_8.png
```

The JSON records the derived fixed-grip geometry factor, but no expected numerical
value is hard-coded. In particular, the previously quoted `Y_fg = 1.2003` is not
accepted as an input; it must be reproduced by the present repository, geometry,
mesh sequence, elastic constants, degradation law, and displacement constraints.

## Promotion rule

This release does not automatically unblock the stochastic first-passage bracket.

- If `J_full/G_fixed_grip` converges near one, the geometry-specific fixed-grip
  artifact can replace the sharp-edge LEFM value in a subsequent audited bracket.
- If only the full-energy J agrees while the tensile-filtered production J does
  not, the production J energy input must be corrected first.
- If neither J agrees, the domain-integral implementation or starter-notch
  representation requires a separate mechanics correction.
