# v10.0.5.14.3 — PF persistent-site and signed-kernel-family parity

This point release ports the audited PF v10.2.22 persistent-site signed MPZ and
the actual v10.2.14 crack-extension-indexed shielding atlas onto the validated
v10.0.5.13.5 full 2-D FEM/CZM mechanics.

## Preserved mechanics

- plane-strain FEM and cubic anisotropic elasticity;
- domain-integral J and equivalent KJ;
- validated 330/240/100 µm mesh/J policy;
- long-corridor Euclidean node deduplication;
- adaptive-CZM topology transactions and quality vetoes;
- stochastic cleavage/first-passage architecture;
- elastic continuum bulk under `bulk_plasticity_mode=tip_only`.

## Persistent-site physics

- two signed reduced BCC slip channels;
- separate positive/negative mobile, retained, and accumulated-slip fields;
- persistent areal source density `rho_source0_m2`;
- no source depletion, temporal source recovery, or crack-advance refresh;
- geometry-controlled multiplicity `M_s = rho_source0 c_arc r_tip w_eff`;
- physical front-width floor `max(w_min_physical,b)`, never MPZ `dx`;
- implicit emission/backstress complementarity;
- unsigned mobile+retained content for Taylor backstress;
- signed retained content for shielding, with mobile shielding fraction zero;
- fractional moving-frame translation and natural resharpening;
- trial/commit and restart preservation of the complete signed state.

## PF signed-kernel family

Production uses:

```text
--signed-kernel-family v10_2_14_active_only_campaign_family.json
```

The loader accepts schema
`v10.2.14_active_only_real_signed_2d_shielding_atlas`, resamples each measured
2×40 active kernel onto the runtime 2×80 MPZ grid, and interpolates by committed
cumulative crack-path extension using the atlas-specified inverse-distance
method. The family-level activation-to-line conversion is retained. Wake
shielding and constitutive shielding caps remain disabled, and extrapolation
outside the measured extension envelope is rejected.

## Transport corrections

v10.0.5.14.1 retained explicit CFL microstepping and failed at the admissible
Peierls rates. v10.0.5.14.2 replaced this with backward Euler, but step-doubling
continued to resolve the numerical stiff tail even after nearly all mobile
content had escaped. The candidate-0118 700 K smoke therefore exhausted 2,000
linear solves at an interval of approximately 5.4e-4 s.

v10.0.5.14.3 preserves the same finite-volume equations but advances each
frozen transport generator with a dense scaling-and-squaring matrix exponential.
The frozen linear Peierls/encounter/Taylor/escape system is therefore solved
without CFL or backward-Euler stiff-tail error. Step doubling remains only for
the nonlinear change of state-dependent coefficients between intervals. Error
normalization uses the line content at the beginning of the accepted macrostep,
so an insignificant residual tail does not force unbounded refinement.

The transport equations and constitutive rates are unchanged. Runtime records:

- `transport_integrator`;
- `transport_attempted_exponentials`;
- `transport_rejected_intervals`;
- `transport_nonlinear_error_max`;
- `max_frozen_courant`;
- `line_content_conservation_error`.

## Install or update

```bash
ROOT=/Volumes/Data/Data/Nanopillar_calculation/Arrhenius_FEM_CZM_MPZ_v10_0_5_14_persistent_sites
cd "$ROOT"
conda activate arrhenius-fem-czm
git pull --ff-only origin v10.0.5.14-persistent-sites-v10222-parity
python -m pip install -e "$ROOT" --no-deps
```

Expected package version: `10.0.5.14.3`.

## Verification

```bash
python -m pytest -q \
  tests/test_persistent_site_transport_v1005143.py \
  tests/test_persistent_site_transport_v1005142.py \
  tests/test_signed_kernel_family_v1005141.py \
  tests/test_persistent_site_v100514.py \
  tests/test_v100513_barrier_only.py \
  tests/test_v1005131_preserved_state.py \
  tests/test_v1005132_startup_resolution_warning.py \
  tests/test_v1005133_tip_only_ramp.py \
  tests/test_v1005134_tip_only_policy_propagation.py \
  tests/test_v1005135_long_corridor.py \
  tests/test_v1005123_phase_c_repairs.py \
  tests/test_mpz_v9_10_unified_transport.py \
  tests/test_v100510_refinement_support.py \
  tests/test_v100511_same_mesh_energy.py
```

The v10.0.5.14.3-specific tests include exact semigroup behavior for a frozen
stiff generator, a Courant number far above the explicit limit, the actual
candidate-0118 Peierls/Taylor parameters at 700 K over 840 s, nonnegativity,
line-content conservation, and zero-content clock advancement.

## Candidate-0118 smoke

The versioned runner creates the output directory before opening its log:

```bash
ROOT=/Volumes/Data/Data/Nanopillar_calculation/Arrhenius_FEM_CZM_MPZ_v10_0_5_14_persistent_sites
cd "$ROOT"
conda activate arrhenius-fem-czm
bash run_v10_0_5_14_3_0118_family_smoke.sh
```

Optional environment overrides include `PFROOT`, `FAMILY_JSON`, `OUTROOT`,
`T_K`, `TARGET_EXT_UM`, `STEPS`, `DU`, `DT`, and `PRINT_EVERY`.

## Mandatory runtime invariants

- `available_site_fraction = 1`;
- `persistent_source_inventory_active = false`;
- `source_depletion_active = false`;
- `source_refresh_active = false`;
- `source_sites_refreshed = 0`;
- `front_width_grid_independent = true`;
- `ahead_of_tip_dx_used_as_front_width_floor = false`;
- `two_channel_drive_reliable = true`;
- `transport_integrator = adaptive_frozen_generator_exponential_v10_0_5_14_3`;
- `transport_cfl_limited = false`.

The remaining promotion gate is the full candidate-0118, 700 K, 20 µm 2-D
smoke using the production PF atlas. The PR remains draft until that runtime
output is reviewed.
