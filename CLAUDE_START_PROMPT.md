Read `CLAUDE.md` and `FEM_CZM_HANDOFF.md` completely before changing files.

Work continuously in the isolated local workspace through audit, implementation,
focused tests, local simulations, diagnosis, and evidence-based correction. Do
not stop after producing a plan.

Begin with Phase 0. Prove the repository root, branch, HEAD, Git status,
environment, package import path, exact PF reference artifacts, current
production entry point, and archived slow-ramp baseline. Then reproduce the
reported low-J loading rate.

The immediate implementation objective is a transactional PF-reference
driving-force controller based on `J_target(t)` or `KJ_target(t)`. It must
preserve physical time and localize stochastic events inside each target
interval without consuming or changing production hazard, MPZ, cohesive,
geometry, bulk-plastic, or RNG state during rejected trial solves.

Use the smallest physics-neutral changes, add a regression for each discovered
bug, run a real production-entry smoke after focused tests, make narrow local
commits, and update `CLAUDE_PROGRESS.md` at every milestone.

Continue automatically through prefracture parity and then a 20–50 micrometer
short-growth qualification unless a stop condition in `CLAUDE.md` is reached.
Do not launch 400 or 1000 micrometer production growth without first reporting
the short-growth acceptance evidence and exact launch contract.
