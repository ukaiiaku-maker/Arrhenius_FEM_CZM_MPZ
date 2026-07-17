# v10.0.1 PF-equivalent kinetic CZM foundation hardening

Branch: `v10.0.1-pf-equivalent-kinetic-czm-foundation-hardening`

Base: `v10.0-pf-equivalent-kinetic-czm`

Reference constitutive implementation: `ukaiiaku-maker/PF-fracture-fatigue`, branch `v10.1.7.1-final-production-temperature-sweep`.

## Correction in this point release

The v10.0 engine installed `CampaignKineticMPZState` after construction but did not override the inherited v9.11 `reset()` method. An explicit reset could therefore replace the campaign state with the v9.11 state while leaving the engine label unchanged.

v10.0.1 adds a reset-safe engine class and versioned abrupt/progressive entry points. Reset now preserves the established v9.11 threshold/control initialization and then reconstructs the same virgin PF campaign state from the immutable material manifest. No temperature-dependent initial state is introduced.

## Scope

This is a foundation hardening release. The following remain blocked:

1. progressive runs longer than one 5 um segment;
2. rejected trial-damage increments are not yet retried automatically at the recommended smaller physical interval;
3. unused physical time after a checkpoint is reported but is not yet carried through a new equilibrium solve;
4. PF/CZM Stage-A parity is not certified until matching PF reference traces are supplied;
5. penalty convergence and the reduced three-class matrix have not run locally.

Do not interpret the existence of the progressive entry point as production certification.

## Checkout

```bash
cd /Volumes/Data/Data/Nanopillar_calculation/Arrhenius_FEM_CZM_MPZ_v10_0_pf_equivalent_kinetic_czm

git fetch origin \
  refs/heads/v10.0.1-pf-equivalent-kinetic-czm-foundation-hardening:refs/remotes/origin/v10.0.1-pf-equivalent-kinetic-czm-foundation-hardening

git worktree add \
  ../Arrhenius_FEM_CZM_MPZ_v10_0_1_pf_equivalent_foundation \
  -b v10.0.1-pf-equivalent-kinetic-czm-foundation-hardening \
  origin/v10.0.1-pf-equivalent-kinetic-czm-foundation-hardening
```

## First gate

```bash
cd /Volumes/Data/Data/Nanopillar_calculation/Arrhenius_FEM_CZM_MPZ_v10_0_1_pf_equivalent_foundation
conda activate arrhenius-fem-czm
CONDA_ENV=arrhenius-fem-czm bash run_v10_0_1_foundation_tests.sh
```

Do not run the progressive smoke until this command passes and its isolated traces have been inspected.
