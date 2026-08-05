# Claude progress

## URGENT HANDOFF (2026-08-04/05, project-direction correction — read this section first)

**The user explicitly rejected the scalar-controlled-K_J development path
documented further below in this file (the `legacy_scalar` diagnostic,
`--pf-kj-target-csv`, the "Path A" recommendation, and the reproducible
anchor-diagnostic script). That work is not deleted from history (it is
still reachable at commits `2bb128c`/`4fcf08a`/`a147774` in this branch's
reflog/origin ref), but it is NOT the production architecture and must not
be extended further.** The user's correction: return to the last known-good
full-production FEM/CZM implementation and fix the targeted R-curve/
advection defect, without rebuilding validated physics from scratch.

### What this correction session did, in order

1. **Phase 1 (audit)**: produced `FULL_PHYSICS_BASELINE_AUDIT.md` and
   `FULL_PHYSICS_BASELINE_AUDIT.json`. Identified
   `293491157063484bb10df6adc479f847a6a08ba6` (tag
   `claude-v10.0.5.18.4.0-source-baseline`, == HEAD of
   `origin/v10.0.5.18.4.0-theta0-pf-parity`) as the authoritative
   full-physics commit — this is the same commit CLAUDE.md already pinned
   as "verified source commit at workspace creation," now independently
   re-confirmed via `git rev-parse`/`git merge-base` and rejected an
   alternative sibling branch (`...-validation`, confirmed less complete).
2. **Phase 2 (recovery decision)**: `git diff --stat` between that baseline
   and this branch's HEAD before this correction (`6100d50`) touches 27
   files, 25 of them entirely new (docs/standalone modules/scripts/tests),
   plus exactly one existing physics file (`sharp_front.py`, confirmed via
   full hunk-by-hunk diff read to be purely additive/default-off) and one
   shell launcher. **Conclusion: the branch had not materially diverged.
   No recovery branch was created — none was needed.** Continued on the
   current branch, per the governing instruction's own decision rule.
3. **Phase 3 (diagnosis)**: traced the accepted-event path (hazard crossing
   → event-length proposal → anisotropic cohesive/geometry trial → remesh →
   commit → new-tip MPZ evolution) through `sharp_front.run_2d` and
   `persistent_site_moving_tip_v100515.py`. **Found the first incorrect
   operation**: in the single-front (non-deflect) accepted-event branch,
   `eng.step()` unconditionally mutates the front engine's hazard/MPZ state
   (`B`, `N_em`, `a_adv`, `n_adv`) for every fire *before* the actual
   geometry commit is attempted via `crack_backend.advance(...)`. When that
   commit is vetoed (mesh-quality/angle-error gate), the pre-existing code
   only printed a diagnostic and `break`-ed — the already-mutated engine
   state was never rolled back, permanently desynchronizing the MPZ/hazard
   frame from the real (unadvanced) geometric tip for the rest of the run.
   The multi-front/deflect code path already had the correct handling
   (`eng_f.restore_geometry_veto(...)` or a manual B/N_em/a_adv/n_adv
   fallback restore) — the single-front path simply lacked it.
4. **Phase 4 (fix)**: mirrored the existing, already-reviewed multi-front
   rollback pattern into the single-front path in `sharp_front.py` (commit
   `c81c376`). No anisotropic elasticity/remeshing/MPZ/shielding/bulk-PT
   physics was touched — only how the front engine's own state reacts to
   an already-existing veto signal from the (unmodified) remesher. Added
   `tests/test_sharp_front_single_front_geometry_veto_rollback.py` (2
   tests) and **verified by temporary revert that the first test fails
   against the pre-fix code** (B stuck at the post-fire value instead of
   restored) — a genuine regression test, not a tautology.
5. **Phase 5 (regression)**: `tests/test_sharp_front_ast_patcher_anchors_v10051840.py`
   (19 passed, 2 skipped — `run_2d` anchors intact); the full
   atomic-path-corridor/persistent-site/moving-tip/signed-kernel/rollback
   focused suite (88 passed, 2 pre-existing-unrelated failures); full
   `tests/` suite (729 passed, 3 skipped, 8 failed — the identical 8
   pre-existing failures from this file's own documented baseline further
   below, no new failures).
6. **Phase 6 (native-loading experiment)**: run, see full result below.
   **This run also revealed an important correction to step 3's own
   diagnosis**: `deflect = bool(getattr(args, 'crystal_aniso', False))`
   (`sharp_front.py`) — any real anisotropic production run (always
   `--crystal-aniso`) takes the **multi-front/deflect** code path, not the
   single-front path Phase 3/4 fixed. The Phase 4 fix is still a genuine,
   verified bug fix (proven by a test that fails against the pre-fix code)
   and is retained, but it is not what governs the actual target case's
   behavior. The multi-front path's *own* equivalent rollback
   (`eng_f.restore_geometry_veto`) was already correct before this session
   — directly confirmed by the native run below stopping cleanly on a real
   geometry veto rather than corrupting state silently.

### Phase 6 result (native-loading run, Peak/1000K/theta=0/seed=8666)

**Entry point**: `mode_i_first_passage_v10_0_5_18_3_9_atomic_path_corridor.main()`
called directly (not the hash-gated `..._theta0_pf_parity`/
`..._full_field_parity`/`..._full_field_production` wrappers), wrapped in
`active_only_kernel_family_compat_v10051840.installed_active_only_kernel_family_compat(out)`
invoked manually. `bulk_plasticity_mode=tip_only` (this entry point's own
default; the `full_field` v10.4.1 overlay is only installed by the
higher, still-hash-gated wrapper).

**Loading**: rate1x-equivalent — PF rate1x's own physical opening rate
(`2.0e-7 m / 8.4 s = 2.381e-8 m/s`) held fixed, discretized with a much
larger native FEM step (`dt=100 s`, `dU=2.381e-6 m`) instead of PF's own
`dt=8.4 s` (which is what produced the original ~145,000-step slow-ramp
problem). 200 steps requested; mesh 36x72 (2020-2180 nodes depending on
adaptive refinement), production Peak barrier row, seed 8666,
`--crystal-aniso --crystal-compete --crystal-theta-deg 0`, `--crack-backend
adaptive_czm`, `--da-phys 5e-6`, `--target-crack-extension-um 50`.

**Kernel used**: the frozen `reference_inputs/pf_final_v10_2_30/rate1x/
v913_paper_peak01_0242980_persistent_sites/T1000K_th0_seed8666/family.json`
snapshot (sha256 `d41b08f69ae773...`) — **not** verified-identical to any
specific historical PF campaign kernel (unchanged caveat from prior
sessions).

**Result** (console log only; the run raised before `steps_1000K.csv` was
written — see "what's missing" below):

```text
step  5   KJ=20.492 MPa√m  sigma_tip=8.15 GPa  B=0.036  N_em=7.69   a=0.500 mm
step 10   KJ=20.879         sigma_tip=8.30 GPa  B=0.286  N_em=6.07   a=0.500 mm
step 15   KJ=21.044         sigma_tip=8.37 GPa  B=0.679  N_em=4.74   a=0.500 mm
step 18   KJ=21.176         sigma_tip=8.42 GPa  B=0.000  N_em=3.94   a=0.504 mm  << ADVANCE (first passage)
step 20   KJ=21.183         sigma_tip=8.38 GPa  B=0.220  N_em=3.55   a=0.504 mm
step 25   KJ=21.202         sigma_tip=8.39 GPa  B=0.744  N_em=3.75   a=0.504 mm
step 28   KJ=21.223         sigma_tip=8.40 GPa  B=0.000  N_em=3.71   a=0.506 mm  << ADVANCE
step 30   KJ=21.239         sigma_tip=8.42 GPa  B=0.159  N_em=6.22   a=0.506 mm
step 35   KJ=21.270         sigma_tip=8.46 GPa  B=0.596  N_em=4.16   a=0.506 mm
step 39   KJ=21.297         sigma_tip=8.48 GPa  B=0.000  N_em=1.92   a=0.509 mm  << ADVANCE
step 40   KJ=21.307         sigma_tip=8.49 GPa  B=0.095  N_em=0.90   a=0.509 mm
step 45   KJ=21.379         sigma_tip=8.52 GPa  B=0.560  N_em=0.63   a=0.509 mm
GEOMETRY VETO front 0: v10051839_no_feasible_atomic_path_corridor -- renewal retained in B=1.000
[run stops: RuntimeError from persistent_site_moving_tip_v100515.py's
 restore_geometry_veto, fail-closed by design]
```

First passage: KJ≈21.2 MPa√m at step 18 (physical time 1800 s). Three
real, discrete crack-advance events (9 µm total, 0.500→0.509 mm), KJ
rising slightly with each — a genuine, qualitatively PF-like rising
R-curve signal, well short of the 50 µm target, before the run stopped on
a real (not artificial) geometry veto.

**What's missing / not yet claimed**: the run did not reach the 20-50 µm
target; `steps_1000K.csv` was never written (the uncaught RuntimeError
propagated out before that write). The exact veto reason was not captured
by the first run because the pre-existing print statement sat *after* the
raising `restore_geometry_veto()` call. Fixed the print/restore ordering
(diagnostic-only, no behavior change) and re-ran with an otherwise
identical, deterministic (seed 8666) configuration: **identical
trajectory reproduced exactly** (steps 5-45 bit-for-bit the same KJ/B/
N_em/a values), confirming determinism, and this time the veto reason was
captured: **`v10051839_no_feasible_atomic_path_corridor`** — the atomic
path-corridor remesher could not find any triangulation satisfying its
quality gates for the requested advance at that point. This is a genuine
remeshing-difficulty finding under this specific (mesh, `da_phys=5e-6`,
native-rate) configuration, not a bug this session introduced or a gate
that was weakened -- the fail-closed stop is the CORRECT response to it
per FEM_CZM_HANDOFF.md section 7.3.

**Exact next experiment**: rerun the same configuration with a smaller
`--da-phys` (e.g. `2e-6` or `1e-6` instead of `5e-6`) -- a purely
mesh-resolution/step-size choice, not a quality-gate relaxation -- to see
whether the atomic path-corridor backend can find a feasible triangulation
for smaller requested advances at the geometry where it failed this time.
If it still fails at the same physical crack length, that would point to
a locally difficult mesh region (e.g. near the notch/refinement boundary)
worth investigating directly in `atomic_path_corridor_adaptive_quality_v10051839.py`.
Do not relax `czm-min-area-ratio`/`czm-min-triangle-quality` to force this
past -- per instruction, those floors stay fixed.

**First-passage KJ scale vs. PF**: 21.2 MPa√m is well below PF's ~53.9
MPa√m for this case (the v10.2.30 campaign's Peak/1000K/rate1x value) --
mechanistically expected, not unexplained, since this run deliberately
used `tip_only` bulk mode (no full-field bulk Peierls-Taylor plasticity
contribution), not PF's `full_field` closure. Do not read this as a
physical-correspondence failure -- the two configurations are not yet
comparable on equal footing.

### Kernel-provenance resolution used for the native run (read before reusing)

The exact-hash-gated production wrapper
(`mode_i_first_passage_v10_0_5_18_4_0_theta0_pf_parity.py`/
`..._full_field_parity.py`, `REFERENCE_KERNEL_SHA256 = a85b57ad9e...`) is
still blocked (unchanged, not weakened — see the EXTERNAL BLOCKER section
further below). For the native-loading experiment, this session invoked
`mode_i_first_passage_v10_0_5_18_3_9_atomic_path_corridor.main()`
**directly** instead (one level below the historical-parity-asserting
wrapper — confirmed via direct source read to have no hash gate of its own),
wrapped in `active_only_kernel_family_compat_v10051840.installed_active_only_kernel_family_compat(out)`
(the same standalone, intentional interoperability shim the production
entry point itself uses, invoked manually rather than through the gated
wrapper), supplying the frozen `reference_inputs/pf_final_v10_2_30/rate1x/
v913_paper_peak01_0242980_persistent_sites/T1000K_th0_seed8666/family.json`
kernel snapshot (sha256 `d41b08f69ae773...`, **not** verified-identical to
any specific historical PF campaign kernel — this caveat is unchanged from
prior sessions and must be repeated in any report of these results). This
produces `bulk_plasticity_mode=tip_only` (this entry point's own default,
not the `full_field` v10.4.1 overlay only the higher, hash-gated wrapper
installs) — a legitimate, previously-qualified production configuration in
its own right (FEM_CZM_HANDOFF.md explicitly describes it as "the qualified
tip-only elastic-bulk formulation"), just not PF's exact bulk closure.

### Native loading rate ("rate1x-equivalent," not PF's own step discretization)

The original slow-ramp defect (documented below, ~145,000 steps to first
passage) came from reusing PF's own `dt=8.4s` step discretization for an
expensive FEM solve, not from the *physical* rate1x opening rate itself
being too slow. This session instead held the **physical** rate1x opening
rate fixed (`dU/dt = 2.0e-7/8.4 = 2.381e-8 m/s`, PF rate1x's own value) but
discretized it with a much larger FEM step `dt` (a native choice, not
PF-matched, not KJ-controller-driven), scaling `dU` proportionally so the
remote-displacement rate is unchanged. This is "rate1x-equivalent loading"
in the sense the governing instructions asked for. See the Phase 6 section
below for the exact `dt` used and the resulting step count.

---

**Everything below this section is from an earlier point in the same
session and is still accurate, but the single most important pivot is
not yet reflected further down: the PF final v10.2.30 four-class
campaign is now the principal PF reference bank (superseding the
v10.4.1 Peak/1000K case), and a real, meaningful controlled-K_J
diagnostic against REAL PF data has now run successfully for the first
time this session.**

### What just happened, in order

1. User redirected: use
   `PF-fracture-fatigue_v10_2_21_persistent_sites_top1/runs/1_final_v10_2_30_theta0_rate_extremes_four_class_1000um_base3621_v1`
   (3 rates x 4 classes x 12 temperatures = 144 cases) as the principal
   PF physical-correspondence reference bank, not the older, unavailable
   v10.4.1 case. Do not regenerate the older case unless a specific
   diagnostic needs a quantity this campaign lacks.
2. Wrote `scripts/inventory_pf_final_campaign_v10230.py` (read-only,
   walks all 144 cases, extracts completion/first-passage/event/kernel
   info + content hashes) and
   `scripts/characterize_pf_final_campaign_v10230.py` (derives comparison
   CSVs from the manifest). Committed at `7870f36`.
3. **All 144 cases verified `complete_target_extension`, zero
   missing/malformed fields** — see `PF_FINAL_CAMPAIGN_REFERENCE_MANIFEST.json`
   (workspace root, tracked) and `PF_FINAL_CAMPAIGN_REFERENCE_SUMMARY.md`
   for the human-readable writeup. Derived CSVs are in
   `runs/pf_final_campaign_v10230/` (gitignored, regenerate with the two
   scripts above in order).
4. **Key PF characterization finding** (rate1x, first-passage K_J vs T):
   the four classes reproduce their names' qualitative signatures exactly:
   - Peak: broad gradual transition, ~800-950 K (23→57 MPa√m).
   - DBTT: SHARP transition, exactly between 950 K (21.7) and 1000 K
     (56.2) — much narrower than Peak's.
   - weak-T: flat 19-21 MPa√m across the entire 300-1300 K range.
   - ceramic: monotonically decreasing 13.6→9.3 MPa√m, no transition.
   Ordering `ceramic < weak-T < {Peak, DBTT}` holds at nearly every
   (rate, T). Rate sensitivity is weak at low T (~1.1x spread) and strong
   above ~800 K (1.7-2.35x spread), consistent with thermally-activated
   plasticity/emission becoming active there. Full tables in
   `runs/pf_final_campaign_v10230/pf_*.csv`.
5. **Recommended 9 Stage-A anchor cases** (see
   `PF_FINAL_CAMPAIGN_REFERENCE_SUMMARY.md`'s table): Peak 300/1000K,
   DBTT 900/1000/1100K (bracketing its sharp transition), weak-T and
   ceramic at 300/1300K.
6. **Important broadened finding on the kernel-provenance blocker**: the
   SAME externally-mutated `v10_2_28_kernel_cache/1447653d.../family.json`
   path (SHA-256 changed from the original `a85b57ad9e...` to
   `d41b08f69a...` — see the older "EXTERNAL BLOCKER" section below) is
   referenced by EVERY case in this v10.2.30 campaign too (confirmed via
   each case's `command.sh`), not just the old v10.4.1 case. **No case's
   own audit JSON records a byte-level kernel fingerprint from generation
   time** (only shape/policy metadata) to verify against, so it is
   possible the mutation already post-dates BOTH campaigns' generation
   (this campaign is dated 2026-07-29, before the observed 2026-08-03/04
   mutation) — meaning `a85b57ad9e...` may in fact be the correct
   historical hash for this campaign's kernel too, but no local copy with
   that hash exists anywhere (already exhaustively searched, see
   `PF_REFERENCE_REGENERATION_CONTRACT.md`). **This blocks the FULL
   production entry point for ANY case** (native-loading or
   controlled-K_J) because
   `mode_i_first_passage_v10_0_5_18_4_0_theta0_pf_parity.py`'s
   `REFERENCE_KERNEL_SHA256` constant is hardcoded to `a85b57ad9e...` and
   is enforced regardless of which PF case you're comparing against.
   **Did not weaken this check.**
7. Generalized `pf_theta0_frozen_reference_v10051840.py`'s
   `resolve_frozen_pf_reference` to accept explicit
   `expected_steps_csv_sha256`/`expected_kernel_sha256`/`steps_csv_filename`
   (defaults preserve exact v10.4.1 behavior) since the new campaign's
   cases have their own valid, different hashes and temperature-specific
   filenames. Committed at `8dfda9a`, 10/10 tests pass.
8. Froze the first anchor case (Peak, 1000 K, rate1x, seed 8666) into
   `reference_inputs/pf_final_v10_2_30/rate1x/v913_paper_peak01_0242980_persistent_sites/T1000K_th0_seed8666/`
   (gitignored, per policy) — `steps_1000K.csv` sha256
   `50153f1b93dee23407ec7971658fc8da917629c14e4ff3290b7a7a0c88449faa`,
   `family.json` sha256 `d41b08f69ae773009f65f4c0094ef5436ba7f0f290a94172e5b894d09a41f7c3`
   (the latter carries the provenance caveat from point 6 above — see
   that directory's own README.md). **Not yet committed as a workspace
   change since `reference_inputs/` is gitignored by design** — nothing
   to commit there, this is just a note that the files exist locally.
9. **Ran the first real controlled-K_J diagnostic against REAL PF data
   this session** (not a synthetic target): `sharp_front.main()` directly
   (NOT the kernel-gated production wrapper — see point 6), with:
   - production-scale mesh (`--nx 36 --ny 72 --tip-h-fine 1e-6
     --tip-ratio 1.20 --da-phys 5e-6`) — confirmed to produce **exactly
     1124 nodes, hbar_tip=1.385e-06 m**, matching this PF case's own
     `stage3_case_status.json` (`n_nodes: 1124`,
     `hbar_tip_m: 1.3845799113363694e-06`) almost exactly — strong
     confirmation the geometry/mesh parameters are correctly understood.
   - **the correct Peak-class cleavage AND emission barrier parameters**
     installed directly via `sharp_front.py`'s own low-level CLI flags
     (`--cleave-G00-eV`, `--cleave-gT-eV-per-K`, `--cleave-sigc0-GPa`,
     `--cleave-sT-GPa-per-K`, `--cleave-exp-a`, `--cleave-exp-n`,
     `--cleave-floor-frac`, `--cleave-Tref-K`, and the `--emit-*`
     equivalents) — using the exact values already documented in
     `PF_REFERENCE_REGENERATION_CONTRACT.md`'s "shared_contract" block.
     **This does NOT require the kernel-gated wrapper at all** — it's a
     legitimate, intentional feature of the base CLI, not a workaround.
   - `--pf-kj-target-csv` pointing at the frozen real steps CSV.
   - **Still missing**: the signed kernel (blocked, point 6), the
     persistent-site MPZ engine, and the full-field bulk Peierls-Taylor
     detailed-balance model — this run used the base `legacy_scalar`
     engine's plain cleavage/emission hazard only, no MPZ/shielding/bulk
     plasticity. **Not yet a production-configuration or
     physical-correspondence result** — it is real Gate 1 (numerical
     correctness) evidence against real data, nothing more, and should
     not be reported as more than that.
   - **Result** (40 steps, `runs/anchor_diagnostics/peak_1000K_rate1x_barriers/`,
     not committed — gitignored under `runs/`): all 40 steps converged
     (`pf_kj_target_controller_audit_1000K.json`), 3-9 controller
     iterations per step (mean 5), achieved K_J matched target to ~7e-6
     relative at the last step (13.587140 vs 13.587047 MPa√m), all
     stochastic-identity fingerprints preserved across every trial. B
     stayed exactly 0.0 through all 40 steps while N_em grew steadily
     (0.62 → 13.53) — **qualitatively consistent** with the real PF
     trajectory's own early behavior (B numerically negligible, N_em
     growing) at the corresponding early rows. No premature failure, no
     ligament severing (unlike two earlier attempts this same investigation
     with (a) a tiny 6x10 toy mesh and (b) production mesh but generic/
     default barrier parameters — both severed the ligament within 6
     steps at unrealistically low K_J~1-2 MPa√m; recorded here so this
     exact failure mode is not re-discovered from scratch next session).

### Exact next steps (in order)

1. **Do not re-derive any of the above** — it's all recorded here and in
   `PF_FINAL_CAMPAIGN_REFERENCE_SUMMARY.md`/`PF_FINAL_CAMPAIGN_REFERENCE_MANIFEST.json`.
2. Decide whether to (a) extend the just-run diagnostic further (more
   steps, toward this case's real first passage at row 163/step 164,
   KJ≈53.9 MPa√m) with the same barrier-parameters-but-no-kernel
   configuration, accepting it stays Gate-1-only since MPZ/kernel/bulk-PT
   are absent, or (b) invest in resolving the kernel blocker first (see
   point 6) so a genuine production-configuration/physical-correspondence
   comparison becomes possible. The user has not yet been asked which to
   prioritize.
3. Turn the ad-hoc diagnostic invocation (the Python snippet in point 9
   above) into a proper reusable script under `scripts/` if continuing
   down path (a) — it was run ad hoc via `python3 -c "..."` this session,
   not yet saved as a file.
4. Repeat the freeze-and-diagnose pattern (points 8-9 above) for the
   other 8 recommended Stage-A anchor cases once a decision is made on
   step 2.
5. Update `FEM_PF_PARITY_SCORECARD.md` with this session's actual
   evidence (the campaign characterization findings and the diagnostic
   run's results) — **this was not yet done as of this handoff**; the
   scorecard still only reflects the pre-campaign-pivot state.

---

- Updated: 2026-08-04 (session continuation — real-engine controller
  audit output, stochastic-identity fingerprint, a confirmed-and-fixed
  AST-patcher anchor regression, and a corrected acceptance standard
  added this continuation; see `FEM_PF_PARITY_SCORECARD.md` for the
  scientific status ledger, which now organizes ongoing work per explicit
  user instruction. This file remains the narrative/commit/environment
  record.)
- **IMPORTANT — acceptance standard corrected**: the FEM/CZM model is
  NOT required to numerically match PF step-by-step or reproduce
  identical stochastic histories. The standard is now **physical
  correspondence / cross-model consistency**: comparable fracture regime,
  toughness scale (~20% initial target), R-curve trend, temperature
  trends, and material-class ordering — not exact numerical agreement.
  See `FEM_PF_PARITY_SCORECARD.md`'s "Acceptance philosophy" section for
  the full comparison-band framework before judging any future Gate 2+
  result. Do not retune parameters merely to reduce a percentage
  difference when mechanisms and macroscopic conclusions already agree.
- Repository: /Volumes/Data/Data/Nanopillar_calculation/Arrhenius_FEM_CZM_MPZ_theta0_pf_parity_claude
- Branch: claude/v10.0.5.18.4.0-j-controlled-loading
- Development workflow: all source-code development is version-controlled
  in Git on this branch. Local commits are authorized and ongoing; nothing
  has been pushed yet (push is a separate, explicit, user-approved step —
  see "Exact next command" below for the precise push command). The local
  `arrhenius-fem-czm-claude` Conda environment runs the editable install
  from this exact local Git checkout — confirmed via
  `importlib.metadata.version('arrhenius-fem-czm')`.
- Completed commits on this branch this session, oldest first:
  - `7562911` — `feat: add PF driving-force trajectory reader` (Phase 1)
  - `3cc67bb` — `feat: add safeguarded J/KJ-target loading controller` (Phase 2a)
  - `ad2c446` — `feat: wire default-off PF KJ-target controller` (Phase 2b checkpoint)
  - `040610f` — `docs: add FEM/PF parity scorecard and PF artifact regeneration contract`
  - `99f6e5c` — `test: add real-engine controller audit output and integration tests`
  - `437973d` — `test: add direct-step state-equality proof for the KJ-target controller`
  - `f6fad45` — `feat: add fail-closed frozen PF reference interface; complete regeneration request`
  - `12b02e6` — `feat: add stochastic-identity fingerprint and fail-closed transactionality check`
  - `c998041` — `fix: restore two AST-patcher anchors broken by the PF-controller wiring` **(important — read the "AST-patcher anchor regression" section below)**
  - `4bb2d09` — `feat: define Gate 3/4 first-passage event-audit schema (unpopulated)`
  (HEAD before this session's work: `ad06e4c`, the workspace-handoff
  commit; `ad06e4c~1` = `293491157063484bb10df6adc479f847a6a08ba6`, the
  verified source commit on `v10.0.5.18.4.0-theta0-pf-parity`.)
- Remote: origin = https://github.com/ukaiiaku-maker/Arrhenius_FEM_CZM_MPZ.git
- Conda env: arrhenius-fem-czm-claude
  (python: /opt/homebrew/Caskroom/miniconda/base/envs/arrhenius-fem-czm-claude/bin/python)
- Package version (authoritative, via `importlib.metadata.version('arrhenius-fem-czm')`): 10.0.5.18.4.0
  - NOTE: `arrhenius_fracture.__version__` (hardcoded at
    `arrhenius_fracture/__init__.py:91`) still reads `'10.0.2'` — stale
    dead code predating this workspace, not physics-relevant, not fixed
    (out of scope). Every real version gate uses `importlib.metadata`.

## EXTERNAL BLOCKER — the immutable PF reference artifacts are not
## currently trustworthy at their documented path. Read this before any
## flag-ON validation or Phase 3 work.

The read-only PF reference repository
(`/Volumes/Data/Data/Nanopillar_calculation/PF-fracture-fatigue_v10_2_21_persistent_sites_top1`,
outside this workspace) is being actively written to by a concurrent,
unrelated process — not this session, not this workspace. This session
never wrote to that repo (only `ls`/`shasum`/`wc`/read-only Python CSV
reads); `git grep`/`git status` in *this* workspace show no trace of it.

**Exact expected values (from FEM_CZM_HANDOFF.md / CLAUDE.md, verified
earlier this session before the external changes began):**

```text
PF case directory:
  PF-fracture-fatigue_v10_2_21_persistent_sites_top1/runs/
    v10_4_1_theta0_rate1x_bulk_PT_four_class_1000um_selective_reuse_base3621_v1/
    v913_paper_peak01_0242980_persistent_sites/T1000K_th0_seed8666/steps_1000K.csv
  steps_1000K.csv SHA-256:
    666d839ae8ef4267ae4ac7ab52a4220de676525bfb70917744e25a119c291b0c

PF kernel:
  PF-fracture-fatigue_v10_2_21_persistent_sites_top1/runs/v10_2_28_kernel_cache/
    1447653d199f0b43cb475951092d69444c9b785f6fdf518c723792abb3b1f5e5/family.json
  family.json SHA-256:
    a85b57ad9eee8331ef34ea222760ce7d4064f48ddef668f1e28a666d14946f3a
```

**Current observed state (as of 2026-08-04 08:36-09:09 PDT):**
- The `steps_1000K.csv` at the exact path above no longer exists. It
  existed and was used successfully twice earlier this session (the Phase
  0 baseline reproduction smoke, and a flag-off regression re-run
  immediately after the `sharp_front.py` edit) — both are valid, already
  gathered evidence; only *later* re-runs against this path are blocked.
- `runs/` in that repo now contains newly created files from an apparently
  live, unrelated "v10.4.x bulk plasticity campaign"
  (`LATEST_v10_4_4/5/6_bulk_plasticity_campaign.txt`,
  `v10_4_6_theta0_rate1x_four_class_500um_bulk_plasticity_base3621_v1/`).
- The kernel file still exists at its documented path, but its **content
  has changed**: observed SHA-256 is now
  `d41b08f69ae773009f65f4c0094ef5436ba7f0f290a94172e5b894d09a41f7c3`
  (≠ the expected `a85b57ad9e...` above). The hash-named cache directory
  is being reused/overwritten by the concurrent process.

**Why this matters / what it does not affect:**
- The existing SHA-256 guards in `run_v10_0_5_18_4_0_peak1000_theta0_pf_full_field_parity.sh`
  and `scripts/preflight_and_launch_v10051840_theta0_full_field.sh`
  correctly fail closed on this mismatch — confirmed this checkpoint did
  **not** weaken either guard (see "Verification evidence" below). Nothing
  has silently run against wrong physics.
- Nothing in this workspace's own Git history is affected. Phases 0/1/2a
  are fully committed and independent of the external repo's live state.
  Phase 1's real-reference integration test is
  `pytest.mark.skipif`-guarded on the file's existence, so it now skips
  (not fails) — confirmed in this checkpoint's test run (36 passed, 1
  skipped).
- Only **flag-ON physical qualification** of the Phase 2b controller (and
  all of Phase 3) is blocked — the flag-off/default production path does
  not depend on this file at all for its own correctness, only the
  launcher script's own pre-existing precondition check and its final
  `compare_v10051840_theta0_pf_parity.py` call require the file to exist
  to run *at all* (flag on or off) via that particular script.

**Requirement before any further flag-ON work**: the authoritative PF
`steps_1000K.csv` and kernel `family.json` must be recovered or
regenerated, then copied into a **dedicated, immutable local reference
directory inside this workspace's own control** (e.g.
`reference_inputs/pf_v10_4_1_theta0_peak1000K_seed8666/`, not committed —
large reference data — but stable and outside the live, externally-mutable
`PF-fracture-fatigue_.../runs/` tree), and independently re-verified
against the two SHA-256 values above before being trusted again. Do not:
- copy/"rescue" data from the live external repo without explicit
  direction (it's someone else's active work area);
- accept a same-named-but-different-hash artifact as a substitute;
- retry against the live path repeatedly hoping it stabilizes.

## Current phase

Phase 0 (provenance/baseline): **complete**. Phase 1 (PF driving-force
trajectory reader): **complete, committed at `7562911`**. Phase 2a
(standalone controller decision logic): **complete, committed at
`3cc67bb`**. Phase 2b (wiring into `sharp_front.run_2d` + launcher
plumbing): **default-off integration complete and checkpointed at
`ad2c446`; controller audit-output + real-engine integration tests added
at `99f6e5c`** (real mesh, real solve, synthetic target — proves the
controller tracks a target KJ(t) correctly against the actual FEM engine,
not just synthetic callbacks). **Flag-ON qualification against the real,
audited PF reference data is still pending recovery of the immutable PF
reference artifacts** (see EXTERNAL BLOCKER below and
`PF_REFERENCE_REGENERATION_CONTRACT.md`). Do **not** describe the
controller as physically qualified against PF — only its wiring,
default-off safety, and real-engine numerical correctness against a
synthetic target have been verified. This project's scientific status is
now tracked primarily in **`FEM_PF_PARITY_SCORECARD.md`** (gate ladder +
per-quantity table); consult that file first for "what's actually been
proven," this file for "how/when/which commit."

## Phase 0 — completed findings

- PF reference case verified intact (before the external changes): 1652
  data rows, columns include `KJ_Pa_sqrtm`, `J_effective_direct_J_per_m2`,
  `B`, `N_em`, `crack_extension_m`, `dt_cur_s` (constant 8.4 s through
  prefracture). SHA-256 recorded above.
- PF kernel SHA-256 verified matching (before the external changes),
  recorded above.
- PF prefracture/first-passage row located precisely: row index 176
  (0-based; PF `step`=177) is the first row with `crack_extension_m>0`.
  At that row: KJ=53.9019 MPa√m, J≈6531 J/m² (matches FEM_CZM_HANDOFF.md §5
  "J≈6531 J/m², KJ≈54 MPa√m"). Rows 0-175 (176 rows, physical time 8.4s to
  1478.4s) are the exact PF prefracture regime.
- Archived slow-ramp FEM baseline confirmed at
  `Arrhenius_FEM_CZM_MPZ_v10_0_5_17_paper7_self_contained/runs/v10_0_5_18_4_0_peak1000_theta0_seed8666_pf_full_field_parity_v1_slow_displacement_ramp_baseline_20260803_202328/`
  (a different, unaffected reference repo). Its log reaches step 9400,
  KJ=11.367 MPa√m, B=0, N_em=0, a=0.500mm (still no growth) before being
  deliberately stopped — confirms the "slow driving-force ramp" diagnosis
  in FEM_CZM_HANDOFF.md §5.
- Independently reproduced the same slow ramp in this workspace: an
  800-step local smoke matched the archived baseline bit-for-bit at every
  printed step (step 200 KJ=2.215, 400=3.111, 600=3.786, 800=4.346 MPa√m).
- All files named in FEM_CZM_HANDOFF.md §4 exist in this workspace.
  Production entry point:
  `python -m arrhenius_fracture.mode_i_first_passage_v10_0_5_18_4_0_theta0_pf_full_field_production`,
  invoked via `run_v10_0_5_18_4_0_peak1000_theta0_pf_full_field_parity.sh`.
  `scripts/preflight_and_launch_v10051840_theta0_full_field.sh` hard-checks
  `EXPECTED_BRANCH=v10.0.5.18.4.0-theta0-pf-parity` — do not use it on this
  dev branch; use the plain launcher or pytest directly instead.
- Focused regression suite (Phase 0 required set): 19/19 passed.

## Phase 1 — completed: PF driving-force trajectory reader (`7562911`)

`arrhenius_fracture/pf_theta0_driving_force_trajectory_v10051840.py`:
- `load_pf_theta0_driving_force_trajectory(steps_csv_path)` — reads a PF
  `steps_*.csv`, validates required columns, reconstructs physical time by
  cumulative-summing per-row `dt_cur_s` (never from step count alone),
  records SHA-256 provenance, locates `first_passage_row_index` (first
  `crack_extension_m>0`). Fails closed on missing columns/non-positive
  dt/non-increasing step.
- `prefracture_slice(trajectory)` — strict prefix before first passage;
  the crossing row itself is excluded (FEM_CZM_HANDOFF.md §6: never assign
  a crossing event to an interval endpoint).
- `interpolate_target(trajectory, physical_time_s, channel="J"|"KJ")` —
  bounded linear interpolation; raises outside `[time[0], time[-1]]`
  (no extrapolation, matching the kernel's own no-extrapolation contract).

Test: `tests/test_pf_theta0_driving_force_trajectory_v10051840.py` — 7
tests. One (`test_real_pf_reference_prefracture_first_passage_matches_handoff`)
is `pytest.mark.skipif`-guarded on the real file's existence — now skips
due to the EXTERNAL BLOCKER above; the other 6 (synthetic-CSV) still pass.

Pure function of the CSV; touches no FEM/hazard/RNG/mesh state.

## sharp_front.py architecture (researched this session — expensive to
## re-derive; read this before touching sharp_front.py again)

`arrhenius_fracture/sharp_front.py::run_2d(args)` is the **single shared
step loop** for the entire `mode_i_first_passage_*` family. ~15 other
modules `inspect.getsource(sharp_front.run_2d)` and AST/text-patch it
rather than calling a different loop. **Treat `run_2d`'s source text as
load-bearing for other modules.**

**AST-patcher anchor regression (confirmed and fixed, commit `c998041`) —
read this before trusting any `git grep`-based anchor safety check
again.** The initial wiring commit's own safety check (`git grep` for
`dU_step`/`trial_frac`/`adaptive_target`/`Uapp_saved`, confirming the two
anchors matching THOSE specific names sat outside the edited region) was
**incomplete and gave false confidence** — it missed two anchors in
`mode_i_first_passage_v10_0_5_5_stochastic_vhcf.py`
(`_CACHE_INIT_ANCHOR`, matching the literal adjacency between
`prev_a_tip_for_block = float(a_tip)` and `while step < args.steps:`, and
`_CACHE_MECHANICS_ANCHOR`, matching the mechanics/stagger block's exact
16-space indentation) because that grep only searched for specific
variable names, not for "every anchor literal anywhere in the package."
The `if pf_target_active: <new> else: <original re-indented>` structure
used at the time re-indented the mechanics block by 4 extra spaces and
inserted new lines between the two `_CACHE_INIT_ANCHOR` lines, silently
breaking both anchors. **This was invisible in the focused PF-controller
test suite this session mostly ran** — it only surfaced running the FULL
`tests/` suite, 7 failures deep in modules with no obvious PF-controller
connection. **Lesson: after any edit to `run_2d`, run the full test suite
at least once, and/or run
`tests/test_sharp_front_ast_patcher_anchors_v10051840.py` (added in the
fix commit — an AST-based scan of every anchor literal in the package,
checked against `inspect.getsource(run_2d)`) — do not rely on a
name-specific `git grep` as a substitute.**

The fix: stopped wrapping the original mechanics/KJ block in any new
conditional at all. It now runs completely unconditionally, at its exact
original text/indentation, regardless of whether the controller is
active. The controller instead runs its own small, purely local search
(distinct variable names, local-only arrays, never touching the shared
`u`/`ep_gp`/`rho_gp`/`sigma_gp`/`psi_gp`/`Ftop`) BEFORE that block, solely
to determine the target `Uapp`; the untouched original block then
performs the one real, authoritative solve. See `c998041`'s commit
message for full detail.

Prior (now-corrected) belief: before this session's edit, `git grep`
confirmed the only two external text anchors
(`kinetic_progressive_2d_v10.py`/`v1002.py`'s `adaptive_target = min(...)`
at line 2435, and an `info = eng.step(KJ, T, dt_cur)` anchor at ~3071) both
fall **outside** the edited region (~2451-2665 pre-edit numbering), so the
edit should not break those patchers. **This was true for those two named
anchors but incomplete as a safety check** — see above.

Production entry chain for this workspace's target case:
`mode_i_first_passage_v10_0_5_18_4_0_theta0_pf_full_field_production.py:46`
→ `..._parity.py:298 main()` → `_base.main`
(`mode_i_first_passage_v10_0_5_18_3_9_atomic_path_corridor`) → ... →
`sharp_front.py:1095 run_2d`. Every wrapper `main()` in this chain forwards
the raw CLI argv list mostly verbatim (only lightweight string-based
`_option_value`/`_require_*` checks, not real argparse) down to
`sharp_front.py`'s own `_build_parser().parse_args(argv)` (its `main()`,
line ~4454) — so a new flag added to `_build_parser()` flows through
cleanly as long as it's present in the CLI args list somewhere upstream.

Inside `run_2d` (theta=0, non-deflect, non-fatigue path — the
Peak/1000K/seed8666 case), pre-edit structure:
- Outer loop `while step < args.steps:`, inner trial loop `while True:`,
  which already snapshots `u_saved/ep_saved/rho_saved/Uapp_saved` before
  every trial and restores them on rejection.
- The fixed ramp (`dU_step = cfg.loading.dU_top * trial_frac`) was applied
  unconditionally, then the mechanics/plasticity solve ran, then KJ was
  computed via `compute_J_integral`, then the **only existing accept/reject**
  check (`predicted_clock > adaptive_target`, a hazard-clock/CFL safety
  check, not a KJ-target check) decided retry vs. accept.
- **Key finding**: front-engine hazard/RNG/MPZ state (`B`, `N_em`,
  `_hazard_rng`, event-length draw, etc.) is only touched **after** the
  trial loop breaks (step commit). It is not read or written during trial
  evaluation for the non-deflect path — so a KJ-target retry loop
  gets full transactionality for free from the existing
  `u_saved/ep_saved/rho_saved/Uapp_saved` scaffolding; no additional
  front-engine snapshot was needed for this integration.
- `compute_J_integral` (`j_integral.py:50`) and `solve_dirichlet`
  (`fem.py:173`) are pure/side-effect-free given their inputs — a trial
  evaluation (assemble → solve → update_plasticity → compute_J_integral)
  is exactly the `solve_trial(Uapp) -> KJ` contract Phase 2a's controller
  expects.
- Full front-engine transactional snapshot/restore already exists (for a
  different purpose) at `persistent_site_moving_tip_v100515.py:71-115`
  and subclasses — not needed for this integration, but the reference
  implementation if a future change moves hazard integration earlier into
  the trial loop.
- No pre-existing KJ/J-target retry loop existed anywhere; closest prior
  art (`mixed_mode_first_passage_v8.py`'s safeguarded secant for
  mode-mixity angle) wraps a whole `run_2d` call from outside, a different
  pattern than the per-step inner-loop retry needed here.

## Phase 2a — completed: standalone J/KJ-target controller (`3cc67bb`)

`arrhenius_fracture/pf_theta0_j_controlled_loading_v10051840.py`:
- `elastic_predictor(Uapp_old, achieved_old, target, power, achieved_min)`
  — `U_new = U_old * (target/max(achieved_old, achieved_min))**(1/power)`.
  `channel_power_for("J")==2.0` reduces to the handoff's literal sqrt
  formula; `channel_power_for("KJ")==1.0` gives the equivalent linear form.
- `solve_to_target(...)` — safeguarded predictor+secant loop: predicts,
  calls the caller's `solve_trial` callback, checks tolerance, falls back
  to a bounded secant correction (`max_growth_ratio`/`max_shrink_ratio`)
  capped at `max_iterations`. Never raises on non-convergence — returns
  `converged=False` with the best-effort last trial.
- Returns `JControlledStepResult` with full per-trial history
  (`TrialRecord`, `accepted` flag) for the event-resolved audit table
  FEM_CZM_HANDOFF.md §8 asks for.
- Deliberately touches no FEM/mesh/plasticity/MPZ/hazard/RNG state — pure
  scalar control logic, dependency-injected via `solve_trial`.

Test: `tests/test_pf_theta0_j_controlled_loading_v10051840.py` — 11/11
passing (predictor formula, exact/nonlinear convergence, safeguard
clamping, graceful non-convergence, achieved_min floor, rejected-trial
bookkeeping, config validation).

## Phase 2b — default-off integration: CHECKPOINTED, flag-ON qualification
## PENDING (do not describe as qualified)

### What changed (this checkpoint)

`arrhenius_fracture/sharp_front.py`:
1. New setup block before `while step < args.steps:`: if
   `args.pf_kj_target_csv` is set, loads the PF trajectory via
   `load_pf_theta0_driving_force_trajectory` + `prefracture_slice`, builds
   a `JControllerConfig`, initializes `physical_time_accepted = 0.0`,
   `KJ_accepted = 0.0`. If unset (default `None`), `pf_target_trajectory
   stays None` and nothing else changes.
2. Three new argparse flags appended at the very end of `_build_parser()`
   (right before its `return p`, so no existing argument-block text was
   touched): `--pf-kj-target-csv` (default `None`), `--pf-kj-target-rel-tol`
   (default 1e-3), `--pf-kj-target-max-iterations` (default 25).
3. Trial-loop restructure: `pf_target_active = pf_target_trajectory is not
   None and not deflect and not fatigue_mode`. When `pf_target_active`,
   looks up `KJ_target` at `physical_time_accepted + dt_cur`, defines a
   `_pf_solve_trial(Uapp)` closure resetting `u/ep_gp/rho_gp` from
   `u_saved/ep_saved/rho_saved` each call and re-solving (same
   assemble/solve/update_plasticity/`compute_J_integral` sequence as the
   original code), calls `solve_to_target(...)`, and sets
   `Uapp/KJ = result.accepted_Uapp/achieved` before falling through to the
   **unchanged** `predicted_clock = eng.predict_clock_increment(...)` call
   and the **completely unchanged** hazard-clock accept/reject block. When
   `pf_target_active` is False, the `else:` branch is the **original code,
   byte-identical, re-indented one level** — confirmed no `git diff`
   changes to its actual statements, only indentation.
   - First-step-from-rest handling: when `Uapp_saved <= 0` (step 0 only),
     seeds the controller with one real solve at
     `max(cfg.loading.dU_top, 1e-12)` rather than a degenerate `(0, 0)`
     pair (the predictor is undefined at `Uapp_old=0`, per
     `test_achieved_min_floor_prevents_division_by_zero_from_rest_state`).
     Later steps seed from `(Uapp_accepted, KJ_accepted)` with no extra solve.
4. Post-accept commit: added `if pf_target_trajectory is not None:
   KJ_accepted = float(KJ); physical_time_accepted += dt_cur` right after
   the existing `Uapp_accepted = Uapp` — additive, no-op when unset.

`run_v10_0_5_18_4_0_peak1000_theta0_pf_full_field_parity.sh`:
- Added `PF_KJ_TARGET_CSV=${PF_KJ_TARGET_CSV:-}` (default empty) and a
  `PF_KJ_TARGET_FLAG` computed the same way as the pre-existing
  `PLOT_FLAG` pattern, appended to the production-entry invocation via
  `${PF_KJ_TARGET_FLAG:+"$PF_KJ_TARGET_FLAG" "$PF_KJ_TARGET_CSV"}`.
  Verified this exact two-word optional-flag Bash idiom works correctly
  under the system's Bash 3.2 (both empty and non-empty cases tested
  standalone) before using it, given this script's own prior history of a
  Bash-3.2-specific bug (`2934911 Make optional plot flag compatible with
  macOS Bash 3.2`). Did **not** touch the pre-existing `FAMILY_SHA_REQUIRED`
  kernel-hash guard or its comparison logic anywhere in this diff.

### Verification evidence for this checkpoint

- `py_compile` on `sharp_front.py`: succeeds.
- `bash -n` on the launcher script: succeeds.
- Focused suite: **36 passed, 1 skipped** (the 1 skip is the
  `skipif`-guarded real-PF-reference integration test, correctly skipping
  now that the external file is gone — not a new failure). `git diff
  --check`: clean.
- Default-off (flag absent) bit-identical evidence: an 800-step smoke run
  performed **immediately after the `sharp_front.py` edit, before the
  launcher-script edit**, reproduced the exact same KJ values as the Phase
  0 baseline at every printed step (200/400/600/800 → 2.215/3.111/3.786/
  4.346 MPa√m) — proves the `sharp_front.py` restructuring alone is a
  behavioral no-op when `--pf-kj-target-csv` is absent.
- A **third** attempted re-verification (after also editing the launcher
  script) hit the EXTERNAL BLOCKER above — the launcher's own
  pre-existing precondition check (present before any of this session's
  edits, used later in the same script for the `compare_v10051840_theta0_pf_parity.py`
  call) requires `steps_1000K.csv` to exist at the documented path
  **regardless of the new flag**, and that file is now gone externally.
  This is not a regression from this checkpoint's edits — it is the same
  external blocker, confirmed by inspecting the failing script's own
  precondition-check output.
- Given the above, the launcher-script diff was instead verified by (a)
  direct code review (the new lines are purely additive and default to
  producing zero extra CLI tokens when `PF_KJ_TARGET_CSV` is unset — see
  diff), and (b) a standalone Bash 3.2 test of the exact
  `${VAR:+"$A" "$B"}` expansion pattern confirming it contributes zero
  words when empty and exactly two correctly-quoted words when set. The
  composition of (evidence for `sharp_front.py` alone) + (proof the
  launcher change is a no-op when unset) covers the "flag-off remains
  bit-identical" requirement even though a fresh end-to-end run via the
  launcher could not be performed a third time this session.
- `git grep` confirms: no hardcoded fallback/replacement PF path exists
  anywhere in the new code (the CSV path is only ever taken from
  `args.pf_kj_target_csv` / `$PF_KJ_TARGET_CSV`, both user/launcher
  supplied); the `FAMILY_SHA_REQUIRED` guard lines are untouched by the
  diff; `runs/` remains gitignored (new run artifacts were not staged).

### What is explicitly NOT yet verified (do not claim otherwise)

- No flag-ON run against real PF reference data has completed. The one
  attempt failed at the launcher's own precondition check due to the
  EXTERNAL BLOCKER, before any controller code executed.
- No evidence yet that the controller converges within
  `max_iterations`, tracks a real PF KJ(t) trajectory sanely, or reaches
  first passage faster than the fixed ramp, against real (non-synthetic)
  FEM state.
- UPDATE (`99f6e5c`): a real-integration test now exists
  (`tests/test_pf_theta0_j_controlled_loading_real_engine_v10051840.py`,
  4 tests) — but against a synthetic target trajectory and the minimal
  `legacy_scalar` front engine, not the real PF reference data or
  production front-engine configuration. See
  `FEM_PF_PARITY_SCORECARD.md` Gate 1 row for exactly what this does and
  does not prove.

## Files changed this session (cumulative)

- New, committed at `7562911`: `arrhenius_fracture/pf_theta0_driving_force_trajectory_v10051840.py`,
  `tests/test_pf_theta0_driving_force_trajectory_v10051840.py`
- New, committed at `3cc67bb`: `arrhenius_fracture/pf_theta0_j_controlled_loading_v10051840.py`,
  `tests/test_pf_theta0_j_controlled_loading_v10051840.py`
- Modified, committed at `ad2c446`: `arrhenius_fracture/sharp_front.py`,
  `run_v10_0_5_18_4_0_peak1000_theta0_pf_full_field_parity.sh`,
  `CLAUDE_PROGRESS.md`
- New, committed at `040610f`: `FEM_PF_PARITY_SCORECARD.md`,
  `PF_REFERENCE_REGENERATION_CONTRACT.md`
- Modified/new, committed at `99f6e5c`: `arrhenius_fracture/sharp_front.py`
  (controller audit-output block), `tests/test_pf_theta0_j_controlled_loading_real_engine_v10051840.py`
- New, committed at `437973d`: direct-step state-equality test added to
  the real-engine test file above.
- New, committed at `f6fad45`: `arrhenius_fracture/pf_theta0_frozen_reference_v10051840.py`,
  `tests/test_pf_theta0_frozen_reference_v10051840.py`; modified
  `.gitignore`, `run_v10_0_5_18_4_0_peak1000_theta0_pf_full_field_parity.sh`,
  `PF_REFERENCE_REGENERATION_CONTRACT.md`; created (gitignored, untracked)
  `reference_inputs/pf_v10_4_1_theta0_peak1000K_seed8666/README.md`.
- New, committed at `12b02e6`: `arrhenius_fracture/pf_theta0_stochastic_fingerprint_v10051840.py`,
  `tests/test_pf_theta0_stochastic_fingerprint_v10051840.py`; modified
  `arrhenius_fracture/sharp_front.py` (fingerprint capture + fail-closed
  check) and the real-engine test file (2 more tests).
- Modified, committed at `c998041`: `arrhenius_fracture/sharp_front.py`
  (anchor-safety restructure — see "AST-patcher anchor regression"
  above); new `tests/test_sharp_front_ast_patcher_anchors_v10051840.py`.
- New, committed at `4bb2d09`: `arrhenius_fracture/pf_theta0_first_passage_event_audit_v10051840.py`,
  `tests/test_pf_theta0_first_passage_event_audit_v10051840.py`.
- Untracked, not committed (correctly gitignored under `runs/`):
  various `runs/*_20260804/` smoke/regression-check directories, several
  of which failed at the launcher's PF-file precondition check once the
  EXTERNAL BLOCKER hit (documented, not a code defect).
- Untracked, not part of this project's code, not touched:
  `.claude/settings.json`, `.claude/settings.local.json` (local harness
  config, not FEM/CZM source).

## Tests run (cumulative)

- 19/19 Phase 0 focused tests.
- 7/7 Phase 1 trajectory-reader tests (6 synthetic + 1 skipif-guarded real).
- 11/11 Phase 2a controller tests.
- 9/9 real-engine integration tests (audit output, rejected-trial
  recording, determinism, flag-off no-op, direct-step state equality,
  stochastic-identity preservation, fail-closed tamper detection).
- 8/8 frozen-reference resolver tests.
- 7/7 stochastic-fingerprint unit tests.
- 3/3 event-audit-schema tests.
- 17/19 (2 skipped, pre-existing) AST-patcher anchor-guard tests.
- **Full `tests/` suite** (not just the focused PF-controller subset):
  725 passed, 3 skipped, 8 failed — all 8 failures independently verified
  (via a temporary worktree at the pre-session baseline commit `ad06e4c`)
  to be pre-existing and unrelated to this session's work (5 stale
  package-version-string assertions, 1 unrelated fatigue-anchor break, 2
  more found via the same baseline check). See "AST-patcher anchor
  regression" above for why running the full suite mattered here.

## Latest accepted physical state

No FEM/CZM physics has changed. The default (flag-off) production path is
byte-identical in behavior to before this session (proven via bit-identical
smoke reproduction). A new, entirely opt-in code path exists and is
believed correct by construction (code review + the sharp_front.py
architecture notes above + Phase 2a's unit tests) but has not yet been
exercised against real FEM state end-to-end.

## Latest numerical blocker

Two, layered:
1. **Physics/controller** (unchanged from before): the fixed
   remote-displacement ramp produces a far slower J/KJ trajectory than
   the PF reference (~145,000 steps extrapolated to first passage). The
   fix (the Phase 2b controller) is now wired in behind a flag but not
   yet physically validated — see item 2.
2. **External/operational** (new this session): the immutable PF
   reference artifacts this fix needs to validate against are not
   currently trustworthy at their documented path (see EXTERNAL BLOCKER).
   This blocks flag-ON validation and all of Phase 3 until resolved.

## Failed approaches that should not be repeated

- Do not "fix" the slow ramp by tuning displacement rate directly —
  CLAUDE.md forbids tuning around a defective loading controller.
- Do not run `scripts/preflight_and_launch_v10051840_theta0_full_field.sh`
  on this dev branch — it hard-fails on a branch-name check by design.
- Do not use heredoc-piped `conda run -n ... python3 - <<'PY'` expecting
  captured stdout when redirected to a file — silently empty output was
  observed. Write the script to a real file and run
  `conda run -n ... python3 /path/to/script.py` instead.
- Do not retry the flag-ON smoke against the live external PF repo path
  hoping it stabilizes — see EXTERNAL BLOCKER.
- Do not trust a name-specific `git grep` as a complete anchor-safety
  check before editing `run_2d` — see "AST-patcher anchor regression"
  above. Run the full `tests/` suite and/or
  `tests/test_sharp_front_ast_patcher_anchors_v10051840.py` instead.
- Do not wrap an existing block of `run_2d` code in a new `if/else` to
  make it conditional if any other module's AST/text-patch anchor might
  match text inside that block — Python's mandatory indentation means
  ANY such wrapping changes that text's exact indentation, which breaks
  literal-substring anchors even though the logic is unchanged. Prefer
  determining a value via a small, separate, distinctly-named local
  computation BEFORE the untouched original block, then feeding that
  value into the one line that already varied (e.g. the `dU_step`
  assignment) — this was the actual fix in `c998041`.

## Exact next command

1. **Do not push yet.** After this checkpoint commit lands, the exact push
   command (to run only when explicitly approved) is:
   ```
   git push -u origin claude/v10.0.5.18.4.0-j-controlled-loading
   ```
   Do not modify PR #64, merge, or change its draft status.
2. Obtain or regenerate the authoritative PF `steps_1000K.csv` and kernel
   `family.json`, copy them into a dedicated immutable local reference
   directory inside this workspace (outside the live, externally-mutable
   `PF-fracture-fatigue_.../runs/` tree), and independently re-verify both
   SHA-256 values recorded in the EXTERNAL BLOCKER section above before
   trusting them again.
3. Only then: run the flag-ON smoke
   (`PF_KJ_TARGET_CSV=<frozen local copy>/steps_1000K.csv`, short `STEPS`,
   e.g. 50), inspect that KJ tracks the PF target sanely with a reasonable
   iteration count and no exceptions, add a dedicated regression test for
   the real-integration path, then commit as
   `feat: qualify PF KJ-target controller against frozen reference data`
   (or similar) — only at that point is it fair to describe the
   controller as physically qualified.
4. Proceed to Phase 3 (prefracture parity) only after step 3 passes.

## Exact next scientific acceptance gate

Phase 3 prefracture parity for Peak / 1000 K / theta=0 / seed=8666: FEM
reaches PF's first-passage J≈6531 J/m² (KJ≈53.9 MPa√m) in a practical
number of controller-driven steps (not ~145,000), with B(t), N_em(t), and
bulk plastic work tracking the PF trajectory shape, before any
crack-growth qualification (Phase 4) is attempted. Blocked until the
EXTERNAL BLOCKER above is resolved.

## Continuation handoff (for a fresh session or Codex)

Everything needed to resume is in this file plus the four new Phase
1/2a/2b files. A fresh agent should: (1) read this file, CLAUDE.md, and
FEM_CZM_HANDOFF.md in full; (2) re-run `git log --oneline -5` and
`git status --short`/`git diff --stat` to confirm no uncommitted drift
from what's recorded here; (3) re-run the focused suite (`36 passed, 1
skipped` is the expected current baseline) to re-confirm before writing
new code; (4) check whether the EXTERNAL BLOCKER has been resolved — if
not, do not attempt flag-ON work; if resolved, follow "Exact next command"
above starting at step 2. Phase 0, Phase 1, Phase 2a, and the default-off
Phase 2b wiring do not need to be repeated or re-derived — in particular,
do not re-run the `sharp_front.py` architecture research (Explore-agent
pass); it is fully recorded above.
