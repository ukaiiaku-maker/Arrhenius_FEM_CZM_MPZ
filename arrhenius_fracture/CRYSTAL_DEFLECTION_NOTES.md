# Crystallographic deflection / branching — build status & notes

Goal (unchanged): emergent crack GEOMETRY — deflection now, branching later —
that varies with temperature, rate, and material/orientation, built ON the
existing validated sharp-front 2D engine (`run_2d` in `sharp_front.py`), using
the group's own hazard / activation-work formulation (paper Eqs 26–30:
`W = σ:A`, `A = φV*`, Schmid + non-Schmid, hazards add `Λ = Σ_s Λ_s`).

Decision log: deflection first, branching later (energy-shared, see below).
Nucleation deferred. Method is the SHARP FRONT (not AT2/AT1 phase-field — those
were explored and set aside; see "What was tried" at the bottom).

---

## DONE and VERIFIED this round

### `crystal.py` (new) — anisotropy foundation
- `cubic_plane_strain_D(C11,C12,C44,θ)` — cubic elasticity by full 4th-order
  tensor rotation. VERIFIED: exact isotropic limit (matches `fem.plane_strain_D`),
  rotation-invariant to 1e-16, W → Zener A = 1.000 (W is elastically isotropic).
- `bcc_cleavage_traces(θ)` / `bcc_slip_traces(θ)` — {100} cleavage and {110} slip
  in-plane traces (candidate advance directions, each with t̂ advance dir + n̂
  crack-normal), carrying crystal orientation θ.
- W single-crystal constants `W_C11,W_C12,W_C44` (523/203/160 GPa).

### `sharp_front.py` — anisotropic elasticity wired into `run_2d`
- New flags: `--crystal-aniso --crystal-theta-deg θ --crystal-C11/C12/C44`.
  Default OFF → isotropic → existing validated runs UNCHANGED.
- VERIFIED: model BCC (A=2) near-tip σ_yy swings 190→305→190 MPa over θ=0→45→90°
  (correct cubic 90° periodicity); W (A=1) is orientation-independent (194 MPa).
  i.e. the near-tip driving field now carries the crystal orientation.

### J-integral is already tip-general (key finding for the next step)
- `j_integral.compute_J_integral(..., crack_tip: (2,))` takes an ARBITRARY tip
  location; `find_crack_tip(mesh,d,a0)` returns 2D tip position AND growth
  direction from the damage field. So deflection needs NO J-contour rewrite —
  the contour re-centers on whatever tip we pass.

---

## HOW THE FRONT WORKS NOW (so the next dev knows what to change)
- Crack = sharp front; state is `a_tip` (x-position on the y=0 ligament).
- Elements with `x ≤ a_tip` (and |y| ≤ half-thickness) are stiffness-killed (d=1).
- `FrontEngine` runs ONE first-passage cleavage clock from the analytic
  `σ_tip = K/√(2π r_eff)`; when `B ≥ 1` the front fires and `a_tip += da_phys`
  along +x. K comes from the J-integral at the current tip.
- So TODAY the advance DIRECTION is hard-fixed to +x. That single scalar is
  what the directional hazard generalizes.

---

## NEXT STEPS (deflection), in order

1. **Per-plane activation work at the tip.** At the sharp tip read the near-tip
   stress σ (already anisotropic). For each candidate plane s in
   `bcc_cleavage_traces(θ)` (+ later slip), form the Schmid resolved drive
   `τ_s = σ:M_s` with `M_s = sym(m_s⊗n_s)` → activation work `W_s = σ:A_s`
   (Schmid only first; non-Schmid `a₁σ_nn+a₂σ_mm` deferred — paper Eq 28).
   Effective per-plane tip stress drives the SAME `FrontEngine` clock per plane.

2. **Directional advance (deflection).** Generalize `a_tip` (scalar) → a tip
   `(x,y)` + heading. The plane with the largest `W_s` (highest hazard / first
   `B_s≥1`) wins; advance `da_phys` along its t̂_s; kill the element in that
   direction (the killed region becomes the swept path, not `x≤a_tip`). Pass the
   moving tip to `compute_J_integral` (already supported).
   - Use `find_crack_tip` for the 2D tip; keep a path polyline for the kill mask.

3. **Render** the deflected path at two θ and two T to SEE it bend with
   orientation and temperature — the deliverable that was missing all along.

## LATER (branching) — design fixed, not yet built
- When TWO well-separated planes are co-critical (comparable `Λ_s` within one
  advance window) → emit two fronts.
- CONSTRAINT (Hank): the two branches SHARE the single-front energy budget — the
  total advance work at the bifurcation must not exceed the single-front energy
  release (no surface energy created from nothing). Each branch gets a fraction
  of `da` summing to ≤ the single-front authorization. This is the (B)
  conservation ledger applied at the branch. Branching is then naturally rare
  (it must be "paid for"), and density should scale with rate/speed and
  microstructure contrast — an emergent trend.

## LATER (anisotropic plasticity) — required for consistency
- Plastic flow must use the SAME slip geometry ({110}/{112}) as the cleavage
  advance, or cleavage and plasticity anisotropy are inconsistent. Slip-system
  resolved shear τ_s = σ:M_s drives the plastic hazard per system; sum hazards.

---

## REFERENCE TARGET (from the sharp-tip code, for validation of trends not values)
- weak-T regime Kc(T): ~13.1, 11.4, 15.8, 14.7, 13.7 MPa√m (flat ~11–16).
- DBTT regime Kc(T): 12.9, 14.0, 37.9, 37.5, 35.9 (lower shelf ~13 → upper ~38 at 900 K).
- One engine spans ceramic→DBTT (1D K-ramp + regime_map confirm; convergence in
  dt / emission-throttle / Kmax). Priority is emergent TRENDS, not exact Kc.

## KEY PITFALLS ALREADY DIAGNOSED (don't repeat)
- AT2/AT1 phase-field smears unless: adaptive ~1µm TIP mesh (`tip_h_fine≈1.2e-6`,
  NOT uniform coarsening) AND `ell` sized to the tip (`ell_absolute_m≈3e-6`, NOT
  `ell_factor×global_hbar≈330µm`).
- AT1 strength is tied to ell: `σ_c=√(3EGc/(8ell))`. At resolvable ell (~3µm),
  σ_c≈280 MPa ≈ the crack-driving far-field (~290 MPa at K=13) → bulk damages.
  Localizing needs sub-µm ell along the WHOLE path (tip-only refinement is not
  enough). This is why the sharp front (no ell) is the chosen route, and why
  PF-CZM (decouples σ_c from ell) is the noted fallback if phase-field is revisited.
- Over-load was the "no T-dependence" cause: loading to KJ~35 when Kc~13 cleaves
  identically at every T. Keep the load near Kc(T) so the first-passage clock
  (which carries the DBTT) can discriminate.

## What was tried and set aside (so it isn't re-litigated)
- Reduced structured-grid AT2 overlay (`at2_overlay.py`): heuristic tip advance
  (skeleton/centroid/PCA/steering) — environment-fragile, straight or stalling;
  DEAD END for geometry. Kept only as an architecture sketch.
- FEM AT2 (`main.py run_simulation`) + source-mode: works but the source-mode is
  the artificial-gate machinery to avoid; over-load gave straight, T-independent
  cracks.
- AT1 (`at1.py`): implemented + verified (analytic threshold), but strength–ell
  coupling above. Kept as a clean baseline; PF-CZM noted for later.
- CURRENT chosen path: sharp front + crystallographic directional hazard (this
  file).

---

## UPDATE — STEP 1-2 DONE: crystallographic deflection works

`crystal.py` additions: `near_tip_stress_tensor`, `pick_cleavage_plane` (max
resolved opening stress sigma_nn,s = n_s.sigma.n_s == W=sigma:A_s for the
cleavage normal projector).

`sharp_front.py` `run_2d` (all behind `--crystal-aniso`):
- 2D tip state `tip_xy` + heading `fwd`; per step picks the {100} cleavage plane
  with max opening stress, takes K from the J-integral around the ACTUAL 2D tip
  along that plane, advances `da_phys` along it on fire, and kills elements along
  the swept path segment (path-following kill, not the y=0 band).
- Engine/kinetics UNCHANGED -> Kc preserved.

VERIFIED (model BCC A=2, C44=320, T=700K, dU=8e-6, dt=1.0, v_rayleigh=1.2e-5):
- theta=30 deg: crack DEFLECTS, climbing ~30 deg along the {100} plane to
  y~0.7mm as it crosses; Kc_first=13.4.  (out -> deflect_theta30.png)
- theta=0 deg : crack stays straight; Kc_first=11.9.  (deflect_theta0.png)
- => path is orientation-controlled, AND toughness is orientation-dependent
  (11.9 vs 13.4) — an emergent bonus from the anisotropic near-tip field.

KEY RUN KNOBS for a GRADUAL (visible) advance instead of one-step runaway:
brittle cleavage is unstable (runs away once K>Kc), so cap the advance velocity:
`--v-rayleigh 1.2e-5` (= ~3*da/dt increments/step) with `--dt 1.0` and `--dU`
large enough to reach Kc. Without the cap the crack crosses in one unstable step
(still correct Kc, but no visible path). For PHYSICAL dynamic fracture use the
real v_rayleigh (~2600 m/s for W) with a co-designed small dt (loading is
per-step `Uapp=step*dU`, clock is per-time `B+=lambda*dt` — they must be balanced;
see the `co-design nu0/tau_c/Kdot` note). The 1.2e-5 value here is a
visualization cap, not physical.

NEXT: (3) branching — when two {100} planes are co-critical (near theta=45, the
two opening stresses are equal: the picker already shows the tie there), emit two
fronts SHARING the single-front energy budget (Hank's constraint). (4) anisotropic
plasticity on the same {110}/{112} slip geometry. Non-Schmid a1,a2 (Eq 28) and
strict-W (A=1, cleavage-plane-only) deflection are incremental add-ons.

---

## UPDATE — gate concern CLOSED: deflection is real with the cap OFF

The `--v-rayleigh` cap is OFF by default (inf); it is NOT in the advance physics.
Confirmed by running v_rayleigh=inf:
- theta=30, T=900, DBTT regime, dU=4e-6: crack DEFLECTS at ~30 deg and crosses
  (deflect_capoff_theta30.png). PATH GEOMETRY IDENTICAL to the capped run -> the
  cap only spread the same path over more frames, it never changed the path.
=> deflection is genuine physics, not a throttle artifact. The cap can be ignored.

FINDING (matters for branching): with the cap off the advance is near-BISTABLE
under displacement-per-step loading (`Uapp=step*dU`):
  dU<=2e-6 -> advance 1 increment then ARREST (deflection drops K, stable);
  dU>=4e-6 -> cross fast (~4 load steps, ~45k increments).
There is no wide stable-tearing window with this loading scheme. A physical
stable-advance window needs the loading to track the R-curve (load-control or a
co-designed dynamic dt with the real v_rayleigh). This is a LOADING/time-stepping
issue, separate from the (verified, cap-free) deflection physics. It will matter
for branching: two fronts can only COMPETE if advance is gradual, so step 3 needs
either the stable-load window or the dynamic time-stepping.

---

## UPDATE — STEP 3 DEMONSTRATED: branching at co-critical {100} planes

`crystal.cleavage_branch_candidates(...)`: returns 1 plane (deflect) or 2
(branch) — a second {100} plane branches when its opening stress >= branch_ratio
* winner AND it is directionally separated (Eq 30 competing hazards). VERIFIED:
single front for theta=0..40, BRANCH (planes @45 & 135) at theta=45.

`sharp_front.run_2d` (flags `--crystal-branch [--branch-ratio 0.92]`): on a
co-critical step a SECOND front is spawned with its OWN engine/clock; both fronts
advance on the ENERGY-SHARED budget (each gets n_fire//2 — total surface created
<= single-front release, per Hank's constraint). Each front: own near-tip stress,
own plane pick, own J-integral around its 2D tip, own kill-along-path.

VERIFIED (theta=45, T=700, model BCC A=2): crack SPLITS at the notch into two
divergent fronts (+45 and -45 along the two {100} planes); Kc_first=13.4.
=> branch_theta45.png.

HONEST caveats (refinement, not blockers):
- The two branches are ASYMMETRIC (one dominates) even at exactly theta=45.
  Expected causes: once one front advances it shields/relieves the other (real
  competition), plus the energy-split + close-front J-integral interaction is
  approximate. Physically branches do compete, but near-symmetric theta=45 should
  be closer to symmetric — worth refining (e.g. simultaneous shared-budget
  allocation by relative hazard, not fixed half/half).
- Binary branch only (cap at 2 fronts); no re-branch / re-merge yet.
- Same loading caveats as deflection (near-bistable advance; v_rayleigh visual
  cap used for frame-spreading — path geometry is cap-independent).

STATUS: deflection (steps 1-2) + branching (step 3) both demonstrated in the
sharp-front engine via the activation-work hazard. Remaining: (a) symmetric/
multi-way branch refinement + shared-budget allocation by hazard; (b) anisotropic
plasticity on {110}/{112}; (c) stable-advance loading (load-control/dynamic dt)
to replace the visualization cap and give a clean stable-tearing window;
(d) non-Schmid a1,a2 and strict-W (A=1) cleavage-plane-only cases.

---

## UPDATE — K-COLLAPSE STALL FIXED + branch-safe J-integral

BUG (Hank's "tip evolves into a large square, not progressing"): after a
deflected advance, KJ collapsed to 0 and stayed there for 40+ steps (B frozen),
while kill_r painted a square. ROOT CAUSE: kill_r = max(notch_half_thickness,
1.5 hbar_tip) = 80 um, while the J-annulus is r_inner..r_outer = 20..80 um -> the
kill blob SWALLOWED the entire contour -> every domain element was d>0.95 ->
skipped -> J=0. The "square" was the 80 um kill disk; the stall was K=0 freezing B.

THREE FIXES (together):
1. THIN kill trail: kill_r = max(hbar_tip, 0.5um) (~1 um crack width, << r_inner),
   so the annulus is no longer swallowed. The trail is now a physical crack width
   (reads as a thin dashed line at plot resolution, not a fat band).
2. exclude_radius (= 2*kill_r) passed to compute_J_integral: a hard disk around
   the tip forced out of the domain, so the contour always stays OUTSIDE the
   freshly-killed material.
3. Branch-safe DOMAIN SEGMENTATION in compute_J_integral (new section,
   future-proofing for branching per Hank): new params
   `crack_segments=[(p0,p1),...]` (all crack polyline segments) and a
   LINE-OF-SIGHT mask -- an element contributes to a tip's J only if the straight
   path tip->element does not cross a crack segment that is not incident to the
   tip. This (a) keeps the contour out of the wake, (b) auto-segments the domain
   between neighboring cracks so the J-contour never runs THROUGH another crack,
   and (c) routes around a second tip when close. Segments are pre-filtered to
   within 3*r_outer of the tip for speed. Both fronts pass _all_segments().
   `find_crack_tip`/`compute_J_integral` already took an arbitrary 2D tip; these
   params are optional (default None -> old behavior), so isotropic runs unchanged.

VERIFIED:
- theta=30, T=950, cap off (Hank's failing case): KJ now RECOVERS after the
  advance (38.3, not 0) -- no square, no frozen B. (Severs in one step here =
  correct brittle fast fracture at that load.)
- theta=30, T=700, gradual (v_rayleigh cap for frame-spreading): clean ~30 deg
  climb, KJ healthy throughout (17.3 -> 11.2 -> 508, never 0), Kc=13.4.
  (deflect_theta30_fixed.png)
- theta=45 branch still splits with the segment-aware J (no error), Kc=13.4.
  (branch_theta45_fixed.png)

The square/stall was purely in the 2D bookkeeping (kill mask + contour), NOT the
engine or the deflection physics (Kc~37 at 950 K matches the DBTT upper shelf and
mesh-converges 36.0->37.6->37.4->38.0). Engine + Kc untouched.

STILL OPEN (unchanged): near-bistable advance (arrest vs fast crossing) under
displacement-per-step loading -> stable-tearing window needs load-control or
dynamic dt; branch symmetry/shared-budget-by-hazard; anisotropic plasticity on
{110}/{112}; non-Schmid a1,a2; strict-W (A=1) case.

---

## HOW TO RUN (canonical calls — update these whenever the interface changes)

NO special flag enables the deflection fixes; they are internal to the
`--crystal-aniso` path. There is NO `kill_r=80um` square in the current code; if
you see a square, you are running an OLD package — re-extract the zip.

### Visible crystallographic deflection (look at the geometry)
The velocity ceiling `--v-rayleigh` is a frame-spreader: with it OFF (inf) the
brittle crack severs in ONE step (correct physics, but a 1-element-wide invisible
diagonal). Add it (sized to dt) to spread the SAME path over steps.
n_cap = v_rayleigh*dt/da_phys; for dt=84, da_phys=5e-6, v_rayleigh=1.2e-7 -> ~2/step.

    python3 -m arrhenius_fracture.sharp_front --mode 2d --nx 50 --ny 100 \
      --tip-h-fine 0.6e-6 --tip-ratio 1.25 --n-stagger 2 --save-snapshots 6 \
      --crystal-aniso --crystal-theta-deg 30 --crystal-C44 320e9 \
      --emit-S-T-c0-kB=-20 --emit-S-T-c1=0.02 --emit-S-sigma-max-kB=8 \
      --multihit-m 3 --multihit-tau 1e-6 --cleave-H0-eV 3.0 --cleave-shield-chi 0.6 \
      --emb-sat-frac 1 --n-sat 2000 \
      --v-rayleigh 1.2e-7 \
      --temperatures 950 --dU 2e-6 --dt 84 --steps 260 \
      --out run_deflect_theta30

### Branching (theta=45, two co-critical {100} planes)
Same as above plus `--crystal-branch` and `--crystal-theta-deg 45`.

### Cap OFF (physical, but severs in one step -> path invisible)
Drop `--v-rayleigh 1.2e-7` (defaults to inf). Useful to confirm Kc and that there
is no square; not for viewing the path.

## This round's fixes (all internal, no new flags)
1. THIN + ADAPTIVE kill trail: kill_r = max(hbar_tip, 0.5um) in the fine zone, and
   per-element max(kill_r, 0.7*sqrt(area_e)) so the path is laid down even where the
   mesh coarsens (the old 80um kill swallowed the J-annulus -> the square).
2. exclude_radius + branch-safe line-of-sight J-mask (see earlier entry): K stays
   valid along the deflected path (verified KJ 41->19->78, never 0).
3. Heading-preserving domain clip `_clip_to_domain`: a runaway step is limited
   ALONG the cleavage plane, not clipped per-coordinate to a corner.

VERIFIED: theta=0/15/30/45 at T=950 -> NO square (all sever cap-off, or climb
visibly with the cap). theta=30 + cap: clean ~30deg climb, KJ healthy, Kc=37.3.

## KNOWN LIMITATION (next real task): mesh refinement must follow the PATH
The mesh is refined at the NOTCH TIP only. A crack that deflects/branches leaves
the fine zone into coarse mesh, where the J-contour (needs local h < ~r_inner ~
2*ell) is under-resolved even though the adaptive kill keeps the path visible.
Proper deflection/branching into the bulk needs refinement ALONG the 2D crack
path (dynamic remeshing following the tip, or a pre-refined wedge covering the
deflection-angle fan ahead of the notch). Until then, keep deflection studies
near the refined zone and/or use gentle coarsening (tip-ratio ~1.1-1.15).

---

## UPDATE — TIP-FOLLOWING REMESH added; branching status clarified

### Tip-following adaptive remesh (the requested capability) -- DONE
The mesh now re-refines at the CURRENT tip wherever the crack goes, so the
J-contour and kill stay resolved through deflection/branching (previously the
fine disk sat at the notch tip only; a deflected crack left it into coarse mesh).

How it works:
- `make_tri_mesh(..., tip_center=(xc,yc))` now refines the radial disk at an
  arbitrary point (was hardcoded (a0,0)). `_radial_ring_nodes` already supported
  this; it just wasn't moved.
- `_remesh_following_tip(geom,mesh_cfg,seed,tip_xy,old_mesh,rho,ep,u,paths,kill_r,
  half_h)`: rebuilds the mesh centered on the tip, transfers per-element history
  (rho, ep) and nodal u by NEAREST NEIGHBOUR (cKDTree), and RE-LAYS the crack
  EXACTLY from the saved polyline(s) (no damage interpolation -> no smearing).
- run_2d trigger: after an advance, if the primary (or branch) tip has moved more
  than R_trigger = 0.4*R_fine from the current refine center, remesh and refresh
  all mesh-derived locals (mesh,bnd,d,rho,ep,u,x,y,cxe,cye,cx_e,cy_e,elem_rad,
  kill_r,adj,refine_center). Default ON for --crystal-aniso; `--no-tip-remesh`
  disables (debug).
- Snapshot renderer now stores per-snapshot nodes/elems and builds the
  triangulation per panel (the mesh changes during the run).

VERIFIED: theta=30, T=950, gradual (cap): crack deflects ~30deg and stays RESOLVED
all the way across (continuous trail to y~0.9mm), KJ healthy throughout
(60->99->177), Kc=37.3. (deflect_theta30_remesh.png)

### Branching -- IS implemented end-to-end (functional first cut)
Chain present and run-verified: detection (cleavage_branch_candidates, fires at
theta=45) -> spawn front2 w/ own engine+path -> per-step front2 near-tip stress,
plane pick, branch-safe J, own clock, ENERGY-SHARED advance (n_fire//2), kill,
path. VERIFIED: theta=45 splits into +/-45 fronts (branch_theta45_remesh.png),
Kc=38.4.
HONEST limits (not yet done): binary only (caps at 2, neither re-branches);
fixed half/half split (-> asymmetric branches; should allocate by relative
hazard); and the remesh centers on the PRIMARY tip only, so the second branch is
less resolved -- DUAL-CENTER refinement (refine around BOTH tips) is the branch
follow-up.

## HOW TO RUN (updated; remeshing is automatic, no flag)
Same visible-deflection call as before -- now stays resolved across the whole
domain because the mesh follows the tip:

    python3 -m arrhenius_fracture.sharp_front --mode 2d --nx 50 --ny 100 \
      --tip-h-fine 0.6e-6 --tip-ratio 1.25 --n-stagger 2 --save-snapshots 6 \
      --crystal-aniso --crystal-theta-deg 30 --crystal-C44 320e9 \
      --emit-S-T-c0-kB=-20 --emit-S-T-c1=0.02 --emit-S-sigma-max-kB=8 \
      --multihit-m 3 --multihit-tau 1e-6 --cleave-H0-eV 3.0 --cleave-shield-chi 0.6 \
      --emb-sat-frac 1 --n-sat 2000 --v-rayleigh 1.2e-7 \
      --temperatures 950 --dU 2e-6 --dt 84 --steps 260 \
      --out run_deflect_theta30

Branching: add `--crystal-branch` and set `--crystal-theta-deg 45`.
Disable remeshing (debug): add `--no-tip-remesh`.

---

## UPDATE — sigma cap diagnosis + {100} VARIANT COMPETITION (plane-gate-global)

### sigma_tip 30 GPa cap: NOT the DBTT mechanism
Verified at T=300, theta=30: Kc_first = 15.54 with cap = 30 / 80 / OFF -- IDENTICAL.
The first advance always fires with sigma_tip BELOW the cap (6-13 GPa), so the cap
never sets the reported toughness. The DBTT (lower shelf ~13 -> upper shelf ~40)
comes from EMISSION-SHIELDING saturation (n_sat, k_shield reducing K_eff), exactly
the "single tip-stress controller + saturation" route. The cap only governs the
deep-overload PROPAGATION stress: uncapped, sigma_tip reaches 90 GPa (~2x W cohesive
strength -> unphysical). RECOMMENDATION: set the cap to a physical cohesive ceiling
(~50-60 GPa for W) rather than 30 or off; since Kc is cap-independent the value is a
propagation-realism choice, not a toughness knob. `--sigma-cap-GPa`.

### Why the crack never switched/branched onto another {100}: the gate
`pick_cleavage_plane` gated the no-reversal constraint against the LOCAL heading
(min_forward vs `fwd`). A crack climbing at +45 deg sees the bend-back -45 deg {100}
variant project cos(90)=0 on its heading -> BLOCKED -> it locks onto the first plane.
Also at theta=30 one variant out-opens the other 3:1 (0.75 vs 0.25 sigma_yy) so there
is no competition there regardless. The two {100} traces in 2D are only competitive
near theta=45 (the +/-45 variants tie by symmetry).

### Fix: `--plane-gate-global` (gate against the macroscopic crack axis)
New flag gates plane choice against the FIXED mode-I crack axis (+x = initial
heading, `fwd0`) instead of the local heading. A -45 kink after a +45 climb is then
"forward" (projects +0.71 on +x) and ALLOWED, while true backward (wake re-entry,
negative +x) stays blocked. Applied to both the primary branch-candidate pick and
the front2 pick. Default OFF (preserves the locked-plane behavior).

VERIFIED at T=300 (cap 60), branch OFF:
- theta=45 -> HORIZONTAL ZIG-ZAG between +/-45 {100} variants: a macroscopically
  mode-I (horizontal) crack built from alternating {100} facets (classic {100}/{100}
  cleavage zig-zag). Was: monotonic climb to y~0.9mm. (zigzag_theta45.png)
- theta=43 -> ASYMMETRIC STAIR-STEP: broken symmetry, net ~-45 deg descent built from
  alternating-variant facets. (zigzag_theta43.png)
- theta=30 -> UNCHANGED monotonic +30 climb (one variant dominates 3:1, no competition)
Path geometry is crystallographic (T-independent); toughness still varies via shielding.

HONEST: near theta=45 the path is now genuinely competitive and HISTORY-SENSITIVE
(which variant is picked first biases the net drift up vs down). theta=45 is symmetric
(horizontal); small offsets bias the drift. The sign of the drift feedback for a given
small offset should be sanity-checked, but qualitatively variant competition is now present.

## HOW TO RUN (variant competition)
Add `--plane-gate-global` and use theta near 45:

    python3 -m arrhenius_fracture.sharp_front --mode 2d --nx 50 --ny 100 \
      --tip-h-fine 0.6e-6 --tip-ratio 1.25 --n-stagger 2 --save-snapshots 5 \
      --crystal-aniso --crystal-theta-deg 45 --crystal-C44 320e9 --plane-gate-global \
      --emit-S-T-c0-kB=-20 --emit-S-T-c1=0.02 --emit-S-sigma-max-kB=8 \
      --multihit-m 3 --multihit-tau 1e-6 --cleave-H0-eV 3.0 --cleave-shield-chi 0.6 \
      --emb-sat-frac 1 --n-sat 2000 --v-rayleigh 1.2e-7 --sigma-cap-GPa 60 \
      --temperatures 300 --dU 2e-6 --dt 84 --steps 200 --out run_zigzag45

Add `--crystal-branch` to let it ALSO spawn a branch where two variants are co-critical.

---
## DETERMINISM FIX (cross-temperature mesh carryover) + cap clarification

**Bug:** tip-following remeshing reassigns the mesh-derived globals (mesh, bnd, x,
y, cx_e, cy_e, adj) that are built ONCE before the temperature loop.  Only `mesh`
was being rebuilt per temperature, so every temperature AFTER the first inherited
the previous crack's far-field-refined mesh AND its stale node coordinates -- the
notch stamp `d[(x<=a0)&(|y|<=half_h)]` then landed on the wrong nodes, corrupting
Kc for all T>first.  Symptom: same temperature run twice gave fire vs no-fire
(amplified by the bistable emission/cleavage race).
**Fix:** at the top of the T loop, rebuild mesh+bnd (notch-centered, seed 42) AND
refresh x, y, cx_e, cy_e, adj.  Verified: `--temperatures 900 900` now bit-identical
(31.221, 123 advances) with remeshing ON.
**Corrected DBTT (theta=30):** 300=15.5, 500=12.5, 700=13.0, 900=31.2, 1100=36.5
(was 15.5/12.4/12.7/40.1/39.9).  Lower shelf ~12.5-13, upper shelf ~31-37 (rising,
NOT the flat 40 -- the flat 40 was the carryover plus saturation).  ANY prior
multi-temperature curve from this engine is suspect and should be re-run.

**sigma cap (`--sigma-cap-GPa`, default 30):** pure clamp on sigma_tip; only other
use is the shelf_audit diagnostic.  cap=30 and cap=200 give IDENTICAL Kc_first
(=31.22 at 900K) -> the cap does not set toughness.  Uncapped, sigma_tip rides
~45 GPa during propagation (at/above W theoretical strength ~E/10~41 GPa), so the
30 GPa cap keeps propagation physical.  The flat upper shelf is emission/shielding
SATURATION (n_sat, emb_sat_frac), i.e. the single-tip-stress DBTT mechanism -- not
a toughness cap.

## PLANE COMPETITION (why theta=30 doesn't switch/branch)
{100} pick is max opening sigma_nn = n.sigma.n.  For a <001> section the two {100}
variants have normals at theta and theta+90, so under mode-I the opening ratio is
sigma_nn(010):sigma_nn(100) = cot^2(theta):  theta=0 -> straight; theta=30 -> 3:1
(one plane dominant, single climb); theta=45 -> 1:1 (degenerate -> BRANCH).  Also
the two forward traces are always 90 deg apart, so once climbing the no-reversal
gate kills the other variant every step -> a single crack CANNOT zig-zag between
two {100}; competition is expressed as BRANCHING near theta=45, not deflection.
To get richer competition (3+ variants, plane switching) needs a different ZONE
AXIS (e.g. <011>) that exposes more cleavage traces in-plane, or allowing {110}
cleavage alongside {100}.

### Determinism: fixed vs residual (platform FP)
- FIXED: cross-temperature mesh-global carryover (gross bug: same T gave fire/no-fire,
  40 vs 31). Verified on one machine: mesh checksum + RNG probe + Kc identical at 900K
  regardless of how many temperatures precede it; `900 900` bit-identical.
- RESIDUAL (NOT a bug): Kc_first carries a ~3% platform band (e.g. 900K = 31.22 on the
  Linux sandbox vs 30.41 on macOS Accelerate). Cause = platform/BLAS floating-point in
  the FEM solve, AMPLIFIED by the near-bistable first-advance race (advance caught on a
  steep coarse-dU ramp). Each platform is internally reproducible. Cure = stable-advance
  loading (pending #2: load-control / dynamic-dt); capturing the advance off the
  bistable knife-edge removes the round-off sensitivity and the platforms converge.
- Corrected DBTT carries this band: lower shelf ~13 (500-700K), upper shelf ~30-36
  RISING (900->1100K). Robust on both machines to within ~3%.
