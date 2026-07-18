# v10.0.5.5 stochastic VHCF hybrid integration

This point release extends the validated v10.0.5.4 cycle-integrated
first-passage implementation. The physical cycle horizon is user-selected and
acts as an experimental censoring cutoff. `1e14` cycles is supported, but it is
not a mandatory target and is not the default full-run horizon.

## Active stochastic mechanisms

- Cleavage renewal uses the existing exponential integrated-hazard threshold.
- Finite crack-tip source emission is realized by reproducible bounded
  Bernoulli/binomial sampling of the existing Arrhenius source probabilities.
- The one-cycle predictor is mean-field and transactional; it does not consume
  random numbers.
- Peierls/Taylor transport, trapping, release and recovery remain deterministic
  conditional on the realized source-emission history in this release.

No barrier surface, source capacity, source refresh length, shielding law,
blunting law, cohesive law or FEM constitutive parameter is changed.

## Hybrid cycle blocks

The scheduler first applies the existing cleavage and MPZ state-change limits.
It then limits accepted cycle blocks by the expected number of state-changing
plastic events:

- rare-event blocks target `RARE_EVENT_TARGET=0.25` expected events;
- high-activity bounded tau-leaps target `TAU_LEAP_TARGET=3` events;
- tau-leaping is selected when the unconstrained block would contain more than
  `TAU_SWITCH_EXPECTED_EVENTS=10` expected events.

Finite source emission is sampled with a bounded binomial distribution, rather
than an unbounded Poisson distribution, because each source site can emit only
once until physical source refresh occurs after crack advance.

## FEM caching

The maximum-load FEM cache is implemented for tip-only, fixed-amplitude,
non-cyclic-mechanics runs. The cache key includes:

- remote displacement;
- mesh size;
- bulk damage sum;
- crack-tip position;
- every cohesive element's topology, damage and clock.

Caching is **off by default** until a cache-off/cache-on equivalence smoke passes.
Enable it with `VHCF_FEM_CACHE=1` only for that validation.

## Local worktree

```bash
BASE="/Volumes/Data/Data/Nanopillar_calculation/Arrhenius_FEM_CZM_MPZ_v10_0_5_4_vhcf"
NEW="/Volumes/Data/Data/Nanopillar_calculation/Arrhenius_FEM_CZM_MPZ_v10_0_5_5_stochastic_vhcf"
BRANCH="v10.0.5.5-stochastic-vhcf-hybrid"

git -C "$BASE" fetch origin \
  "+refs/heads/${BRANCH}:refs/remotes/origin/${BRANCH}"

git -C "$BASE" worktree add \
  -b "$BRANCH" \
  "$NEW" \
  "refs/remotes/origin/${BRANCH}"

cd "$NEW"
```

If the local branch already exists, omit `-b "$BRANCH"` and use the existing
branch in the worktree command.

## Compilation and tests

```bash
python -m py_compile \
  arrhenius_fracture/stochastic_campaign_v10055.py \
  arrhenius_fracture/mode_i_first_passage_v10_0_5_5_stochastic_vhcf.py \
  arrhenius_fracture/mode_i_first_passage_v10_0_5_5_stochastic_vhcf_audited.py \
  run_v10_0_5_5_stochastic_vhcf_delta_sigma.py \
  run_v10_0_5_5_stochastic_vhcf_delta_sigma_compat.py

bash -n run_v10_0_5_5_stochastic_vhcf_delta_sigma.sh

pytest -q \
  tests/test_v10053_audited_source_preflight.py \
  tests/test_kinetic_fatigue_v10053.py \
  tests/test_v10054_vhcf_first_passage.py \
  tests/test_v10055_stochastic_vhcf.py
```

## First stochastic smoke

```bash
MODE=smoke \
MATERIAL=DBTT \
TEMPERATURES="700" \
DELTA_SIGMA_MPA="350" \
R=0.1 \
FREQUENCY_HZ=1000 \
CYCLES_MAX=1e5 \
MAX_BLOCK_CYCLES=inf \
MAX_BLOCKS=1000 \
TARGET_EXTENSION_UM=5 \
STOCHASTIC_SEED=1 \
RARE_EVENT_TARGET=0.25 \
TAU_LEAP_TARGET=3 \
TAU_SWITCH_EXPECTED_EVENTS=10 \
VHCF_FEM_CACHE=0 \
RESOLVE_CYCLIC_MECHANICS=0 \
OUTROOT=runs/v10_0_5_5_DBTT_700K_350MPa_stochastic_smoke_seed1_v1 \
bash run_v10_0_5_5_stochastic_vhcf_delta_sigma.sh
```

Required audit fields:

```text
status = complete or right_censored only at the requested physical limit
legacy_scalar_predictor_calls = 0
predictor_mean_field_calls > 0
event_statistics = stochastic
stochastic_emission_active = 1
constitutive_surfaces_changed = false
```

A low-stress smoke is not required to produce a cleavage event. It must show a
finite, reproducible stochastic source history and correct physical-horizon
accounting.

## Reproducibility and seed divergence

Run seed 1 twice in distinct output folders. The complete stochastic audit,
source history and cycle blocks must match. Then run seeds 2-5. At least one
source-event trajectory should differ while all conservation and completion
gates remain valid.

Using the same seed at different stress amplitudes intentionally supplies common
random numbers for paired stress comparisons. Independent replicate curves are
obtained by repeating the complete stress sweep with different seeds.

## Cache equivalence gate

Use the same seed and identical physical controls for two runs:

```bash
VHCF_FEM_CACHE=0 OUTROOT=runs/v10055_cache_off_seed1 ...
VHCF_FEM_CACHE=1 OUTROOT=runs/v10055_cache_on_seed1  ...
```

The following must agree to numerical tolerance:

- cycle blocks and limiter sequence;
- stochastic source events per block;
- `KJmax` and `Delta KJ`;
- cleavage action and stochastic threshold state;
- mobile, retained, escaped and available-source populations;
- crack events and final termination.

Only after this gate should `VHCF_FEM_CACHE=1` be used for production sweeps.

## Flexible cycle horizons

Suggested horizons are operational, not physical constants:

- smoke: `1e5`;
- pilot: `1e8`;
- standard full run: `1e12`;
- maximum experimental censoring study: explicitly set `CYCLES_MAX=1e14`.

A case terminates earlier whenever stochastic cleavage first passage and the
requested crack-extension criterion are reached.
