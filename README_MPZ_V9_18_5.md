# MPZ v9.18.5: committed-target stop, strict mesh quality, and Mode-I corridor

Branch: `v9.18.5-target-stop-quality-corridor`

## Why v9.18.5 is needed

The v9.18.4 ceramic 700 K gate physically reached 60 um and passed all
committed-event and persistent-wake audits. It also exposed three long-growth
problems:

1. the FEM loop continued for more than 13,000 steps after the committed target;
2. the final topology updates accepted triangle qualities of order 1e-4--1e-3;
3. the tip-local mesh scale grew above the 5 um physical crack quantum after the
   crack left the notch-centered refinement patch.

The v9.18.4 angular endpoint regularization was therefore useful as a software
gate but is not retained for production long growth.

## Exact endpoint and rigid-body constraint

v9.18.5 returns to the exact edge-aware insertion introduced in v9.18.3. The
singular matrix seen there is treated as a boundary-condition problem rather
than a reason to perturb the crack path.

Once a fully failed crack separates the upper and lower bulk bodies, each body
must have a minimal horizontal rigid-body restraint. The original model fixed x
only on bottom-boundary nodes. v9.18.5 identifies connected bulk components and,
for any component without an x restraint, holds the current x displacement of
one grip-boundary node. This removes only incremental rigid translation and does
not add cohesive stiffness or change the crack-opening law.

## Immediate physical-target stop

The accepted-step horizon is now dynamic. After the final cohesive event commits
its MPZ translation and the physical target is reached, the next evaluation of
the accepted-step loop terminates. The original nominal `--steps` value remains
available for reporting and all normal output writers execute.

## Production mesh gates

Every accepted topology transaction must satisfy:

- finite positive element areas;
- minimum affected-element triangle quality;
- minimum child-to-parent area ratio;
- no orphan bulk nodes;
- bulk support for every cohesive endpoint;
- maximum tip-local `h/da` ratio.

Defaults:

```text
minimum triangle quality       0.035
minimum child area ratio       0.08
maximum tip h / physical da    0.75
```

A failed transaction is rolled back and the Arrhenius renewal remains available.
The v9.16 retry loop cannot silently relax these production floors.

## Pre-refined Mode-I corridor

The initial graded mesh is built from overlapping refinement centers along the
full requested Mode-I growth length plus a guard distance. Defaults:

```text
corridor center spacing   25 um
corridor guard             10 um
```

For a 100 um target this keeps the tip inside the fine corridor instead of
letting it enter the coarse far-field mesh after approximately 50 um.

## Physics unchanged

v9.18.5 does not change:

- cleavage or emission barriers;
- absolute-hazard integration;
- stochastic/deterministic event statistics;
- cohesive-opening progress;
- source refresh;
- Peierls--Taylor state evolution;
- persistent-wake transport;
- wake shielding;
- the 5 um physical crack quantum.

## Required diagnostic

First rerun only ceramic to 60 um in a fresh output root:

```bash
TEMPS="700" \
CLASSES="ceramic" \
TARGET_EXT_UM=60 \
STEPS=15000 \
OUTROOT_BASE=runs/mpz_v9_18_5_ceramic_700K_60um_v1 \
bash run_mpz_v9_18_5_persistent_plastic_wake_sweep.sh
```

Acceptance requires:

- 60 um committed extension;
- `v9185_stop_requested=true`;
- final accepted step close to the target-commit step rather than 15,000;
- no singular matrix or NaN mechanics;
- no quality-gate relaxation;
- minimum accepted triangle quality >= 0.035;
- minimum accepted child area ratio >= 0.08;
- tip `h/da <= 0.75`;
- complete v9.16/v9.17/v9.18 audits.

The case directory writes `target_stop_quality_corridor_v9185.json` with the
corridor centers, component anchors, quality vetoes, and target-stop state.

Only after this gate passes should the three-class 100 um shielding-on campaign
be run.
