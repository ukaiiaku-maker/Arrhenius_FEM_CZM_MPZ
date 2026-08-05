# Full-Physics Baseline Audit

Produced in response to an explicit project-direction correction: return to
the last known-good full-production FEM/CZM implementation and integrate
only targeted R-curve/advection fixes, rather than continuing the
`legacy_scalar` controlled-K_J diagnostic path. See `FULL_PHYSICS_BASELINE_AUDIT.json`
for the machine-readable per-mechanism table this file summarizes.

## 1. Authoritative baseline commit

```text
293491157063484bb10df6adc479f847a6a08ba6
"Make optional plot flag compatible with macOS Bash 3.2"
2026-08-03 15:58:31 -0700
tag: claude-v10.0.5.18.4.0-source-baseline
== HEAD of origin/v10.0.5.18.4.0-theta0-pf-parity
```

This is not a filename-based guess. It is:

- the exact commit CLAUDE.md pins as "verified source commit at workspace
  creation";
- the exact commit CLAUDE_PROGRESS.md's "HEAD before this session's work"
  section identifies as `ad06e4c~1`;
- confirmed via `git rev-parse` to be simultaneously the tip of the
  `claude-v10.0.5.18.4.0-source-baseline` tag and of
  `origin/v10.0.5.18.4.0-theta0-pf-parity`.

A sibling branch, `origin/v10.0.5.18.4.0-theta0-pf-parity-validation`
(`3791cab9`), was checked and rejected as an alternative candidate: it is
**not** a descendant of the source-baseline commit (they diverge from a
common ancestor `6352a4c4`), and it only adds a CI-isolation GitHub Actions
workflow on top of that ancestor -- it is **missing** the kernel-compatibility
layer, the production entry-point wiring, the preflight launcher, and the
parity test file that the source-baseline commit has. The source-baseline
commit is the more physics-complete of the two.

## 2. Current-branch divergence from the baseline

`git diff --stat` between the baseline commit and the current branch HEAD
(`6100d50`, 16 commits ahead) touches 27 files: **25 are entirely new files**
(documentation, standalone `pf_theta0_*.py` modules, new scripts, new
tests) that do not modify any existing physics file, plus exactly **one**
existing physics file (`arrhenius_fracture/sharp_front.py`) and one existing
shell launcher.

A full hunk-by-hunk read of the `sharp_front.py` diff confirms every change
is a new, additively-gated code path (`if pf_target_active:` / `if
pf_target_trajectory is not None:`), with the `else:` branch in every case
being the pre-existing statement, unindented and byte-identical -- this was
a deliberate design constraint from the commit that fixed an earlier
AST-patcher anchor regression (`c998041`), precisely to avoid touching any
text that other modules' `inspect.getsource(run_2d)`-based patchers depend
on. No anisotropic elasticity, remeshing, MPZ, persistent-site, shielding,
bulk Peierls-Taylor, or cohesive code was touched by this session's prior
scalar-controller work.

**Conclusion: the current branch has not diverged materially from the
known-good full-physics baseline.**

## 3. Recovery-strategy decision (Phase 2)

Per the given decision rule ("If the current branch still contains the full
implementation intact, continue on it and bypass/remove only the scalar
workflow"): **continue on the current branch.** No new recovery branch, no
cherry-picking, and no reset are needed -- there is nothing to recover from,
because nothing was lost. The scalar-controller code in `sharp_front.py`
remains present but is dead code (default-off) for every production entry
point; per the correction's own allowance, it may remain as an isolated
root-finder regression and will not be extended, used as a scientific
result, or built upon.

A safety reference exists regardless: the previous local commits made
during the rejected scalar-diagnostic continuation
(`2bb128c`/`4fcf08a`/`a147774`) are recorded in `origin`'s ref for this
branch and in reflog, should any of that (non-physics) documentation or
script content ever be wanted again -- they are not being deleted, only not
built upon.

## 4. Production composition chain (verified by direct source read + one archived empirical run)

```text
mode_i_first_passage_v10_0_5_18_4_0_theta0_pf_full_field_production.py:main
  -> installs active_only_kernel_family_compat, then
mode_i_first_passage_v10_0_5_18_4_0_theta0_pf_full_field_parity.py:main
  -> monkeypatches persistent-site policy, bulk-PT model (exact v10.4.1
     port), and continuum-bulk field metadata; hard-gates
     --signed-kernel-family against REFERENCE_KERNEL_SHA256; then calls
mode_i_first_passage_v10_0_5_18_3_9_atomic_path_corridor.py:main
  -> the anisotropic remeshing / atomic path-corridor lineage
     (FEM_CZM_HANDOFF.md: demonstrated ~400 um growth without the
     previous mesh failure); ultimately reaches
arrhenius_fracture/sharp_front.py:run_2d
  -> the shared step loop every wrapper AST/text-patches or otherwise
     composes with.
```

The persistent-site engine (`PersistentSiteMovingProcessZoneFrontEngineV100514`
and its `PersistentSitePFMovingTipFrontEngineV100515` subclass) is wired in
by `mode_i_first_passage_v10_0_5_14_persistent_site.py`, which monkeypatches
`mixed_mode_first_passage_v9_11.MovingProcessZone2DFrontEngine` to the
persistent-site class; that name is what `mixed_mode_first_passage_v8`'s own
`sf.build_engine = _engine_factory(...)` mechanism actually instantiates --
**not** `sharp_front.build_engine`'s own direct `state_model == 'moving_pz'`
branch, which is bypassed once this monkeypatch is installed.

The exact intermediate hops from `mode_i_first_passage_v9_18_5_6.py` down
through the v9.18.x/v9.11 lineage to the point where that `build_engine`
swap is actually installed were **not exhaustively re-traced hop-by-hop**
this session -- doing so by static reading alone, across roughly ten nested
legacy wrapper layers, was judged lower-value than empirical confirmation.
Instead, an **archived, already-completed run** from the read-only
`Arrhenius_FEM_CZM_MPZ_v10_0_5_17_paper7_self_contained` reference repo
(`runs/v10_0_5_18_4_0_theta0_full_field_preflight/`) was inspected directly:
its `persistent_site_production_manifest_v10_0_5_18_3_9.json` records
`front_engine.continuous_fractional_MPZ_translation=True` and
`front_engine.kinetic_tip_cell=True` -- literal fields from
`persistent_site_moving_tip_v100515.py`'s own `audit_payload()` -- and its
`theta0_pf_full_field_parity_contract_v10_0_5_18_4_0.json` records
`bulk_parity.model = "v10.4.1_bulk_peierls_taylor_detailed_balance_exact_port"`
and `exact_selected_registry_row=true`. This is direct, empirical proof the
full chain executes correctly end-to-end, not an inference from filenames.

## 5. Kernel-provenance status (unchanged, not weakened)

The hard SHA-256 gate in `mode_i_first_passage_v10_0_5_18_4_0_theta0_pf_parity.py`
(`REFERENCE_KERNEL_SHA256 = a85b57ad9e...`) remains exactly as coded. No
locally-available kernel file currently hashes to that value (see
`CLAUDE_PROGRESS.md`'s "EXTERNAL BLOCKER" section for the incident record),
so the fully historical-parity-asserting entry point cannot currently be
run end-to-end with a byte-verified claim of exact provenance.

This does **not** block a native-loading run (Phase 6): that gate's entire
purpose is asserting numerical identity with one specific historical PF
case (it also hard-forces `nx`, `ny`, `dU`, `dt`, etc. to match that case's
own controls). A native-loading run is FEM's own independent prediction
under its own loading schedule and does not need to satisfy this specific
historical-parity claim. It can invoke the same full physics through a
lower-level entry point (e.g. `mode_i_first_passage_v10_0_5_18_3_9_atomic_path_corridor.py`
directly, or `mode_i_first_passage_v10_0_5_14_1_persistent_site_family.py`)
supplying a documented, currently-available kernel snapshot with its hash
recorded and the provenance caveat stated explicitly -- never claimed as
historically verified, and never routed through the hash-gated wrapper's
exact-match assertion.

## 6. Mechanism-by-mechanism status

All fourteen required mechanisms were located in the current working tree
and confirmed to be ancestors of the baseline commit. Full detail (file,
symbol, composition module, launcher flags, focused test, archived-run
evidence, current-branch status) is in `FULL_PHYSICS_BASELINE_AUDIT.json`.
Summary:

| Mechanism | Status |
|---|---|
| Anisotropic elasticity and crystal orientation | unchanged |
| Anisotropic crack direction | unchanged |
| Anisotropic remeshing / mesh-quality corrections | unchanged (must not be replaced) |
| Atomic cohesive-zone trial/commit/veto | unchanged |
| Signed active shielding | unchanged (kernel-gated run currently blocked, see above) |
| Persistent-site moving MPZ engine | unchanged |
| Crack-tip Peierls transport (continuous MPZ translation) | unchanged -- **primary Phase 3 focus** |
| Bulk full-field Peierls-Taylor evolution | unchanged |
| Emission-derived Peierls-Taylor barriers | unchanged |
| Stochastic cleavage and emission | unchanged |
| Event-length sampling | unchanged |
| Energy-gated crack extension | **gap**: no separate FEM-side energy-budget gate module analogous to PF's; the FEM analog is hazard-threshold crossing + mesh-quality veto. A genuine architectural difference, not a bug -- not to be papered over by inventing a new gate. |
| Exact rollback | unchanged |
| Moving-tip state transfer | unchanged -- **secondary Phase 3 focus** (see distinction below) |
| Restart support | **gap**: no restart-from-checkpoint mechanism located for this production chain (only snapshot-writing, not snapshot-resuming). Not built speculatively; flagged for when actually needed. |

## 7. Where Phase 3 should look first

Two distinct files govern crack-tip microstructure transport, and they must
not be conflated:

1. `persistent_site_transport_v1005144.py` (`installed_split_transport_v1005144`)
   governs **within-step** population exchange/advection along the
   tip-relative x-coordinate inside a fixed MPZ frame (glide/encounter/
   Taylor-release, an exact physical generator matrix, backward-Euler with
   step doubling). It already has its own conservation invariants and a
   three-file test suite (`test_persistent_site_transport_v1005142/3/4.py`).
   This is very likely **not** the site of the R-curve bug.
2. `persistent_site_moving_tip_v100515.py`'s `_integrate_coupled` is what
   actually **re-anchors the MPZ frame's origin** when the crack tip
   advances: it accumulates `dB` from the integrated cleavage hazard rate,
   translates by `da = self.f.da * dB` via `self.mpz_state.advance(da)`,
   and only requests a cohesive-geometry checkpoint once accumulated `B`
   reaches 1.0 (`self.a_adv += self.f.da`). This is a **continuous,
   physically-integrated** translation, not a fixed nominal event length --
   which is reassuring, but Phase 3 must still verify this fractional
   `da` accumulated *before* a commit is consistent with the *actual
   committed* geometric advance `a_new - a_old` once the cohesive backend's
   atomic remesh actually executes (a discrepancy between the two is
   exactly the kind of gap the governing instructions describe).

## 8. Explicit non-actions this phase

No physics code was written or modified. No recovery branch was created
(none was needed). No PF campaign files were touched. No kernel-hash guard
was weakened.
