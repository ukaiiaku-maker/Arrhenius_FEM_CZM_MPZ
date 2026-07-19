# v10.0.2 PF-equivalent kinetic CZM event lifecycle

Branch: `v10.0.2-pf-equivalent-event-lifecycle`

Base: `v10.0.1-pf-equivalent-kinetic-czm-foundation-hardening`

## Purpose

This point release closes the two physical-time lifecycle gaps that blocked the first progressive FEM/CZM smoke:

1. a trial cohesive increment rejected because `delta damage` is too large is restored and retried at a reduced physical time interval;
2. physical time left after a committed 5 um checkpoint is retained at the same applied load, a new trial segment is inserted, mechanics/J is recomputed, and only then is the remaining interval consumed.

The release also includes `dot_ep` in the complete trial rollback history. The prior progressive transform remapped this field during topology insertion but did not restore it during a full trial rollback.

## Preserved constraints

- PF v10.1.7.1 material rows remain unchanged.
- The reset-safe `CampaignKineticMPZState` remains active.
- Opening, cleavage, and emission stress channels remain separated.
- The active MPZ is translated continuously before commit.
- Commit performs zero additional MPZ translation or source refresh.
- One trial topology event is admitted per re-equilibrated mechanics state.
- No cohesive critical-opening, traction, or energy failure criterion is added.
- Wake shielding, smeared variational fracture, temporal source recycling, empirical source caps, and temperature-dependent fitted controls remain disabled.

## Validation status

Implemented but not locally certified by the repository authoring environment:

- lifecycle unit tests for rejection/retry, commit-time carry, target stop, and tiny physical intervals;
- guarded compilation against the actual production `sharp_front.run_2d`;
- versioned v10.0.2 progressive entry point;
- package provenance bumped to `10.0.2`.

Longer progressive runs remain blocked until the lifecycle foundation tests and one-segment FEM smoke pass locally. Penalty convergence is the next gate after the one-segment smoke.

## Checkout

```bash
cd /Volumes/Data/Data/Nanopillar_calculation/Arrhenius_FEM_CZM_MPZ_v10_0_1_pf_equivalent_foundation

BRANCH=v10.0.2-pf-equivalent-event-lifecycle
WORKTREE=../Arrhenius_FEM_CZM_MPZ_v10_0_2_event_lifecycle

git fetch origin \
  refs/heads/${BRANCH}:refs/remotes/origin/${BRANCH}

git worktree add \
  -b ${BRANCH} \
  ${WORKTREE} \
  origin/${BRANCH}
```

## First gate

```bash
cd /Volumes/Data/Data/Nanopillar_calculation/Arrhenius_FEM_CZM_MPZ_v10_0_2_event_lifecycle
conda activate arrhenius-fem-czm
CONDA_ENV=arrhenius-fem-czm bash run_v10_0_2_lifecycle_tests.sh
```

Do not launch a progressive FEM run until this gate passes.
