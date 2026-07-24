# v10.0.5.15 — PF update-map and continuous moving-tip parity

This release restores constitutive time-update parity with the active PF
v10.2.22 benchmark at commit
`198ece3aeb1d193a8c1c4857676fba720c088d27` while retaining the validated full
2-D FEM/J/adaptive-CZM mechanics.

## Local signed-state update

For every constitutive interval, FEM/CZM now uses the PF sequence:

1. implicit persistent-site emission against Taylor backstress;
2. analytic mobile/retained exchange in every MPZ bin;
3. zero explicit mobile and retained recovery;
4. one population-weighted scalar Peierls velocity for the active MPZ;
5. fractional conservative mobile advection with absorbing escape;
6. the corresponding zero-stress wake exchange and advection.

The forest density used by transport is the unsigned retained population, as in
the PF scalar transport path. The v10.0.5.14.2–14.4 coupled backward-Euler and
step-doubling transport solver is not installed in this release.

## Crack-opening and microstructure coupling

The outer FEM equilibrium and two-channel tensor drive remain fixed during one
accepted mechanics interval. Inside that interval the front engine uses the PF
kinetic moving-tip Strang sequence:

1. evolve emission, transport, retention, escape, shielding, and blunting for a
   half interval;
2. recompute the cleavage rate from the evolved tip state;
3. advance the virtual crack tip and translate the MPZ by
   `da_checkpoint*dB`;
4. evolve the microstructure for the second half interval.

The operation repeats on internal substeps until the accepted time is consumed
or the cleavage progress reaches one checkpoint. At that point adaptive CZM is
asked to insert one physical cohesive crack increment. The MPZ is not translated
again at the checkpoint because it has already moved continuously during crack
opening.

The outer adaptive predictor performs the same coupled calculation on a trial
copy. A rapidly developing plastic zone can therefore reduce the predicted
cleavage increment before the mechanical step or cohesive event is accepted.

## Fail-closed cohesive topology

If adaptive CZM rejects a completed continuous-tip checkpoint, the constitutive
state is restored to the renewal origin and the calculation stops. Continuing
would leave the MPZ ahead of the physical cohesive crack and is forbidden.

## Install

```bash
ROOT=/Volumes/Data/Data/Nanopillar_calculation/Arrhenius_FEM_CZM_MPZ_v10_0_5_15_pf_parity

cd /Volumes/Data/Data/Nanopillar_calculation

git clone \
  --branch v10.0.5.15-pf-update-map-moving-tip-parity \
  --single-branch \
  https://github.com/ukaiiaku-maker/Arrhenius_FEM_CZM_MPZ.git \
  "$ROOT"

cd "$ROOT"
conda activate arrhenius-fem-czm
python -m pip install -e "$ROOT" --no-deps
```

Expected package version: `10.0.5.15`.

## Tests

```bash
python -m pytest -q \
  tests/test_pf_update_map_v100515.py \
  tests/test_moving_tip_coupling_v100515.py \
  tests/test_signed_kernel_family_v1005141.py \
  tests/test_persistent_site_v100514.py \
  tests/test_adaptive_czm_tip_support_v1005144.py \
  tests/test_persistent_site_diagnostics_v1005144.py
```

## Production-atlas campaign

```bash
MAX_JOBS=2 \
bash run_v10_0_5_15_0118_300_1200K_200um_family_campaign.sh
```

The v10.0.5.14.4 and v10.0.5.15 outputs must not be combined. All temperatures
must be rerun because both the local transport map and crack-opening/
microstructure coupling changed.
