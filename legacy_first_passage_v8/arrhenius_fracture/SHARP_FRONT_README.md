# Sharp-front dual-hazard fracture (`sharp_front.py`)

A rebuild of the crack-advance layer for the tungsten DBTT model. It replaces
the AT2 ↔ hazard coupling (source mode, Griffith license, consume halo,
frontier gates, fired-memory relief) with a **single crack-advance law** on a
sharp ligament front. The FEM, plasticity, J-integral, and energy-audit
modules are reused unchanged.

## Why the rebuild

The previous `--preset emergent` model had three structural problems
(diagnosed from the `dbtt_v5` run):

1. **Kinetics were gated out of relevance.** With the constant-entropy shelf,
   the saturated cleavage rate was ~0.7/s; with `dt=1 s` and `B_relax=5 s` the
   clock sat permanently past `B_target` from early loading. What actually
   controlled advance at 300 K was the Griffith license `K_J ≥ K_G` — a
   hand-set input. The "emergent Kc" was `K_G`, not the Arrhenius physics.
2. **Both plasticity couplings were stabilizing.** Plastic-zone growth shielded
   *and* relaxed the tip stress for both channels, so emission and cleavage
   both got harder — the system could only store elastically forever (500–700 K)
   or avalanche via the consume hack (300 K). The intended *divergence*
   (PZ growth suppresses emission but embrittles cleavage) was not in the code.
3. **ρ was decoupled from dissipation.** `Wp/Wext = 0%` while `ρ_max` climbed
   to 5×10¹⁵ — a free regularizer, and the energy audit was off by ~85×
   (`K_J` domain vs global).

`sharp_front` removes the two-criterion conflict entirely: the phase field is
only a stiffness-kill indicator for broken material; it never evolves
variationally and there is no second fracture criterion to reconcile.

## The model

One sharp front on the `y=0` ligament. Per-front state: emitted-dislocation
ledger `N_em`, cleavage first-passage action `B`, blunted tip radius `r_eff`.

```
sigma_tip = K / sqrt(2*pi*r_eff)            analytic; never de-smeared FEM
r_eff     = r0 + c_blunt*b*N_em             emission blunts          (shields, transient)
sigma_back= beta*G*b/(2*pi*(1-nu)*L_pz)*N_em pile-up back-stress, LINEAR in N
                                            -> suppresses EMISSION only   (–)
dG_emb    = e_stored(N_em) * v_emb          stored PZ energy LOWERS
                                            the CLEAVAGE barrier          (+)

lambda_e = nu0_e * exp(-G*_e(sigma_tip - sigma_back, T)/kT)
lambda_c = nu0_c * exp(-(G*_f(sigma_tip, T) - dG_emb)/kT)
           multi-hit renewal: lambda_eff = gammainc(m, lambda_c*tau_c)/tau_c
B += lambda_eff*dt ;  front advances floor(B) increments/step (renewal count)
```

The advance is a **renewal count**, not one fire per step: in time `dt` at rate
`lambda_eff` the first-passage process completes ~`lambda_eff*dt` events, so the
front advances `floor(B)` increments and keeps the remainder. Once K exceeds Kc
this is >> 1 and the crack **runs away** in a few load steps (unstable fast
fracture), as brittle cleavage must. A single-fire-per-step cap would instead
throttle the crack to `da_phys/dt` -- a numerical speed limit -- so a 300 K
crack at 3x Kc would crawl for hundreds of steps rather than severing.

The two couplings have **opposite signs** (emission self-suppresses via back
stress; cleavage is embrittled by stored energy), so the channels diverge and
cannot deadlock. On advance the crack enters fresh material: the tip
re-sharpens, a wake fraction of the ledger is retained, and `B` resets
(renewal).

The DBTT is the **race** between the two clocks under a K-ramp, plus the
blunting feedback (emission → blunting → σ_tip suppression → emission runs
further ahead). It is not a gate or a hand-set crossover.

## The stress-gated entropy (the key physics)

A negative *affine* activation entropy, `S = -S0*(1 + sigma/sigma0_S)`, is
`-S0` even at zero stress, so it collapses the **zero-stress** barrier with T —
cleavage then fires thermally at any stress at high T and there is no ductile
regime. The **gated** form

```
S(sigma) = -S0 * x^n / (1 + x^n),   x = |sigma| / sigma0_S
```

is `0` at `sigma=0` (cold zero-stress barrier preserved) and saturates at `-S0`
under load (threshold rises with T only where the tip is stressed). This is the
structure required for a DBTT — and, per the fatigue draft, the same physics
implicated in **low-T fatigue limits**: the constrained transition-state
entropy exists only under load.

Roles in the calibrated defaults:
- **Cleavage** = cold-easy channel: high `H0` (2.2 eV), *weak* gated entropy →
  nearly T-independent threshold. Wins the race below the DBTT (brittle).
- **Emission** = cold-hard / hot-soft: high `H0` (1.8 eV) + *strong* gated
  entropy (`S0=10 kB`) → frozen below the DBTT, softens steeply above it.
  Pre-empts cleavage and blunts the tip above the DBTT (ductile).

Sweep the entropy form with `--entropy-form {gated,affine}`. `affine` gives no
ductile regime (the pathology); `gated` gives the DBTT. `--cleave-S0-gate-GPa`
and `--emit-S0-gate-GPa` expose `sigma0_S` for each channel.

## Result (1D, `--Kdot 0.005 --multihit-m 3 --multihit-tau 1e-6`)

```
T     Kc [MPa√m]   N_em    sigma_back   r_eff     mode
300     10.4        26      0.25 GPa    1007 nm   brittle
400      8.8        27      0.27 GPa    1008 nm   brittle
500      7.6        29      0.28 GPa    1008 nm   brittle
600      6.6        32      0.31 GPa    1009 nm   brittle
700      5.7        37      0.36 GPa    1010 nm   brittle
800      4.6       125      1.21 GPa    1034 nm   ductile   <- DBTT
900      1.1       552      5.35 GPa    1151 nm   ductile
```

Brittle Kc falls smoothly with T (cleavage-controlled); at the DBTT emission
accumulates explosively, the back stress climbs to GPa, the tip blunts, and the
mode flips. The crossover emerges from the dynamics, not a switch.

## Usage

```bash
# 1D validation / calibration sweep (seconds):
python -m arrhenius_fracture.sharp_front --mode 1d \
    --temperatures 300 400 500 600 700 800 900 \
    --Kdot 0.005 --multihit-m 3 --multihit-tau 1e-6 --out runs/sf1d

# contrast the entropy form (the fatigue-limit lever):
python -m arrhenius_fracture.sharp_front --mode 1d --entropy-form affine ... # no DBTT
python -m arrhenius_fracture.sharp_front --mode 1d --entropy-form gated  ... # DBTT

# 2D FEM-coupled run (FEM supplies K_J; engine is the only advance law):
python -m arrhenius_fracture.sharp_front --mode 2d \
    --temperatures 300 700 --steps 120 --nx 60 --ny 120 \
    --multihit-m 3 --multihit-tau 1e-6 --out runs/sf2d
```

## Physically-grounded activation entropy (`entropy_stress_form='physical'`)

The activation entropy is no longer a tuned shape. It is a composite of three
documented pieces, used for BOTH channels but with physically asymmetric
magnitudes:

```
S*(sigma,T) = S_T(T) + S_sigma(sigma)
S_T(T)      = clip(c0 + c1*T + c2*T^2, S_T_min, S_T_max)      [kB]
S_sigma(s)  = -S_sigma_max * x^n/(1+x^n),  x = sigma/sigma0_S  [kB]
```

- `S_T(T)`: the experimentally/atomistically-derived BASELINE (Allera et al.,
  Nat. Commun. 2025; Veverka & Dillon fatigue draft). For W dislocation glide
  this is a saturating polynomial in T -- the piece the ORIGINAL emission model
  had (T-dependent, no stress dependence). The dislocation core carries ~kB per
  b (Schoeck 1980, eq. 34), so the EMISSION channel gets a large baseline.
- `S_sigma(sigma)`: the Schoeck thermoelastic / interaction entropy (Schoeck
  1980, eqs. 36/40; draft eq. 37). Zero at zero stress (a homogeneous strain
  gives no harmonic frequency change -- Schoeck sec. 4) and more negative under
  load. This is the NEW stress dependence, added to (not replacing) S_T.
- CLEAVAGE gets S_T ~ 0 and a tiny S_sigma: bond rupture has no dislocation-core
  mode reorganization, only the weak -(1/mu)dmu/dT modulus term. This asymmetry
  (large emission entropy, negligible cleavage entropy) is the physically-
  grounded replacement for the old symmetric gated/gated setup.

What the sweep recovers: with PHYSICAL entropy magnitudes the DBTT is controlled
jointly by the emission baseline magnitude (S_T_c0) and the emission ENTHALPY H0
(a measurable dislocation-nucleation activation energy), not by a tuned entropy
shape. The default baseline runs S_T(T) from ~-14 kB (cold) to ~-2 kB (900 K)
-- more negative at low T, per the W glide data. Less-negative baselines lower
the DBTT (e.g. emit-S-T-c0 -18 -> ~900 K, -12 -> ~700 K); a very negative
baseline (|S_T| large across the range) suppresses emission entirely and gives
no ductile regime -- a genuine constraint, not a free knob.

IMPORTANT (fixed): the S_T saturation window defaults are now [-40, 0] kB. An
earlier build clamped the ceiling at -5 kB with a positive slope, which clipped
the whole baseline to -5 across 300-900 K and made the emit-S-T-c0 sweep axis
INERT (identical results for c0 = -2/-4/-6). The axis is now live.

Single-run flags: `--emit-S-T-c0-kB`, `--emit-S-T-c1`, `--emit-S-sigma-max-kB`,
`--cleave-S-sigma-max-kB`, plus `--emit-H0-eV`/`--cleave-H0-eV`.
Sweep axes: `--sweep-emit-S-T-c0-kB`, `--sweep-emit-S-T-c1`,
`--sweep-emit-S-sigma-max-kB` (combine with the H0 and form axes).

```bash
python -m arrhenius_fracture.sharp_front --mode sweep \
    --temperatures 300 400 500 600 700 800 900 \
    --Kdot 0.005 --multihit-m 3 --multihit-tau 1e-6 \
    --sweep-emit-S-T-c0-kB -25 -18 -12 \
    --sweep-emit-S-sigma-max-kB 0 4 8 \
    --out runs/phys_sweep
```

CAVEAT: the polynomial coefficients here are illustrative placeholders. To make
this predictive, fit S_T(T) to the actual W activation-entropy data (Allera
2025 / the draft's polynomial) and set S_sigma_max from the Schoeck
-(1/mu)dmu/dT estimate for W. The machinery is in place; the numbers are not yet
the measured ones.

## Entropy sweep (`--mode sweep`)

The stress dependence of the activation entropy is the main open physics lever
(and the link to low-T fatigue limits). Because 1D runs in seconds, the full
parameter space is cheap to sweep. `--mode sweep` runs the complete 1D
temperature sweep for every entropy condition in a grid, saving each condition
to its own subdirectory plus a combined table of derived scalars (DBTT, brittle
and ductile Kc shelves).

```bash
# strong levers: form (gated vs affine) x emission gate stress
python -m arrhenius_fracture.sharp_front --mode sweep \
    --temperatures 300 400 500 600 700 800 900 \
    --Kdot 0.005 --multihit-m 3 --multihit-tau 1e-6 \
    --sweep-form gated affine \
    --sweep-emit-S0-gate-GPa 2.0 3.0 4.0 \
    --out runs/entropy_sweep
```

Sweep axes (each takes a list; the Cartesian product is run):
- `--sweep-form gated affine` — entropy stress form. THE dominant lever: affine
  collapses the zero-stress barrier and gives NO ductile regime; gated gives the
  DBTT.
- `--sweep-emit-S0-gate-GPa ...` — emission gate stress sigma0_S; shifts where
  emission unfreezes in the stress ramp (DBTT location).
- `--sweep-gate-power ...` — Hill exponent n; gate sharpness (DBTT abruptness).
- `--sweep-emit-S0-kB ...` — emission entropy magnitude (weak lever on DBTT
  location at fixed gate stress, in the tested regime).
- `--sweep-cleave-S0-kB ...` — cleavage entropy magnitude.

Axes left unset fall back to the calibrated default for that parameter. The
same overrides exist as single-value flags for one-off runs: `--entropy-form`,
`--gate-power`, `--emit-S0-kB`, `--emit-S0-gate-GPa`, `--cleave-S0-kB`,
`--cleave-S0-gate-GPa`.

Outputs:
- `runs/entropy_sweep/cond_NNN_<tag>/` — full 1D results per condition
  (`kc_vs_T.json`, `kc_vs_T.png`, per-T `trace_*.csv`).
- `runs/entropy_sweep/sweep_summary.csv` and `.json` — combined table:
  (form, gate_power, emit_S0, emit_S0_gate, cleave_S0) ->
  (DBTT_K, Kc_brittle_shelf, Kc_ductile_shelf, n_brittle/ductile/nofracture).

NOTE: resolve the DBTT with a dense temperature grid (the example uses 7
points); a coarse grid pins the apparent DBTT to the nearest sampled T. The
DBTT scalar is the ONSET OF THE CONTIGUOUS HIGH-T DUCTILE BLOCK -- walking down
from the highest T, the lowest T such that it and all higher-T points are
ductile -- so its resolution is the temperature spacing. (This replaced an
earlier "lowest ductile T anywhere" definition that a single low-T point
grazing the ductility threshold could corrupt.)

### Brittle/ductile classification (continuous criterion)

A point is **ductile** when the emitted pile-up significantly shields the crack
tip, by either of two continuous, scale-free measures:

```
sigma_back / sigma_tip > ductile_shield   (default 0.3)   OR
r_eff / r0             > ductile_blunt    (default 1.02)
```

i.e. the back-stress is shielding the tip, or the tip has substantially
blunted. Both are the "is emission winning the blunting race" question. This
replaced a raw `N_em > 50` count cutoff, which sat in the MIDDLE of the brittle
cluster (N_em ~ 40-55 at all T in the affine case) and flipped labels on a
rounding-level difference (a 300 K point at N_em=50.1 was spuriously called
ductile). The genuine ductile onset is a wide gap away (N_em > ~100,
sigma_back > 1 GPa, r_eff/r0 > 1.03), so the continuous criterion is stable.
Tunable with `--ductile-shield` and `--ductile-blunt`; the legacy `--ductile-N`
still exists but is no longer the discriminant.

## Mesh independence (2D)

Refining the mesh changes only the *resolution of the picture*, never the
physics. All physical length scales are absolute and decoupled from the element
size `hbar`:

- **`r0`, `L_pz`, barrier lengths** — absolute by construction (the tip field
  is analytic, `sigma_tip = K/sqrt(2*pi*r0)`), so the 1D physics and the DBTT
  are exactly mesh-independent.
- **Crack advance per cleavage event** is a PHYSICAL increment `da_phys`
  (`--da-phys`, default ~5*r_pz), NOT one element. Sub-element advances
  accumulate and kill elements only as the physical tip crosses them, so the
  renewal/first-passage dynamics do not depend on element size. (Previously one
  event advanced exactly one element, so `hbar` silently set the cleavage step.)
- **J-integral contour radius** is an ABSOLUTE length `r_J` (`--rJ`, default
  ~10*L_pz), floored at 3 elements. K_J is then comparable across refinement
  rather than measured over a `4*hbar` radius that shrank with the mesh. (Quick
  check: K_J at fixed step changed only ~3% between a 30x60 and 60x120 mesh,
  i.e. genuine discretization convergence.)
- **Crack band half-width** is the physical `notch_half_thickness`, not
  `max(half_h, 0.6*hbar)`, so a coarse mesh no longer draws an artificially
  thick crack.
- **Wake density deposit** spreads over the physical pile-up area `L_pz`,
  floored to one element so a coarse mesh still registers (not drops) the
  deposit.

The driver prints both `hbar` (resolution) and the physical lengths, and emits
explicit NOTEs when `hbar > da_phys` (advance under-resolved) or `hbar > L_pz`
(process zone under-resolved) — so an unconverged run is a reported condition,
not a silent smear.

**Caveat:** with the default 1 um process zone on a ~mm domain, even a 60x120
mesh has `hbar ~ 40x L_pz`, so the pile-up/back-stress terms are still
mesh-floored. To converge those you need `hbar < L_pz` (very fine or locally
refined mesh), or raise `L_pz`/`r0` to physical values you can afford to
resolve. A mesh-convergence check (same T, three densities, compare `Kc_first`
and the snapshot fields) is still worth running before trusting any 2D number
quantitatively — the code now makes that check meaningful, but does not replace
it.

Outputs: per-T `trace_*.csv` / `steps_*.csv`, `kc_vs_T.json` / `summary.json`,
and a `kc_vs_T.png` (1D).

### 2D field snapshots (restored)

The 2D driver renders a per-temperature `field_snapshots_<T>K.png` panel, the
sharp-front analogue of the original emergent-model snapshots. Four rows over
columns sampled across the run (step 1, evenly spaced steps, every advance, and
the final step):

- **damage d** — the sharp ligament front (stiffness-kill indicator); shows
  where the crack actually is, not a smeared phase field.
- **log10 rho** — dislocation density. Uniform in the brittle regime (no
  emission); a growing tip cloud in the ductile regime, deposited where
  emission did work.
- **sigma1 FEM (MPa)** — raw principal stress; the tip concentration and its
  plastic relaxation. Shared color scale across columns.
- **eq. plastic strain** — crack-tip plastic lobes; flat in the brittle
  regime, butterfly lobes developing in the ductile regime. Shared scale,
  floored at 0 (a near-zero field renders as zero, not an autoscaled artifact).

The title of each panel annotates the front scalars (`KJ`, `N_em`, `a_tip`)
since those are front quantities, not fields. Control with
`--save-snapshots N` (columns to capture) and `--snapshot-cols N` (max panel
width). The brittle-vs-ductile contrast is directly visible: compare a 300 K
panel (uniform rho, flat plastic strain, stationary sharp tip) with a 900 K
panel (growing rho cloud, plastic lobes, climbing N_em).

## Diagnostic plots (2D)

Each 2D temperature writes a diagnostic suite (sharp-front analogues of the
legacy main.py figures), plus a cross-temperature `toughness_vs_temperature.png`:

- `load_displacement_<T>K.png` — reaction force vs applied opening.
- `toughness_<T>K.png` — K_J (domain integral) with the crack length a(U) on a
  twin axis (the resistance curve).
- `dislocations_<T>K.png` — rho distribution (mean/p95/p99/max) vs opening.
- `hazard_clocks_<T>K.png` — the renewal-clock competition: cleavage rate
  lambda_c, emission rate lambda_e, and advances-per-step. THIS is the
  sharp-front story (cleavage wins cold, emission wins hot) and has no legacy
  analogue.
- `tip_state_<T>K.png` — sigma_tip (drive), sigma_back (pile-up), r_eff/r0
  (blunting), N_em. The sharp-front analogue of the old tip-memory plot.
- `energetics_<T>K.png` — U_el, W_p, W_emit growth. See caveat below.

**Deliberately omitted** (no sharp-front analogue — these were artifacts of the
AT2-coupling layer that was removed):
- local-Gc toughening state (`toughening_state`) — there is no degradable Gc
  field; cleavage is a hazard, not a Gc.
- M_tip amplification / z_shield (`tip_memory` upper curves) — replaced by the
  physical blunting r_eff/r0 in `tip_state`.
- phase-field energy E_pf — no variational damage energy exists.
- yield / plastic / thermo-admissible fractions — AT2 projection diagnostics.

**Energy balance (fixed):** the `energetics` panel IS a closed balance.
A FEM bug previously left the elastic energy density `psi_e_gp`
UNDEGRADED while the stress and stiffness were scaled by g_d=(1-d)^2,
so stiffness-killed (cracked) elements reported spurious full energy
(sum psi*A was ~3.6x the exact 1/2 u^T K u with a notch present).
psi_e_gp is now degraded by g_d in both assemble_mechanics and
stress_state. Verified: on a clean elastic block all four measures
(1/2 u^T K u, masked/unmasked psi integral, 1/2 F*delta) agree to 6
digits, and with a notch U_el tracks W_ext exactly in the brittle
pre-fracture regime. The panel now also plots the residual
W_ext - U_el - W_p - W_emit (the fracture surface work + discretization
error). The J-integral consumes the same corrected psi; Kc_first is
unchanged (10.98 vs 1D 10.4 at 300 K).

## Wake-ledger conservation (renewal advance)

When the front advances `n_fire` increments in a step, the emitted-dislocation
ledger is split between what the re-sharpened tip RETAINS and what is SHED into
the cracked wake:

```
N_retained = N_pre * wake_retain**n_fire
N_shed     = N_pre * (1 - wake_retain**n_fire)   # deposited as wake density
```

The 2D wake density deposit uses `N_shed` (the density actually left behind),
NOT a fraction of the post-shed remainder. The two agree when `n_fire == 1`,
but during runaway (`n_fire > 1`) `wake_retain**n_fire -> 0`, so the post-shed
remainder is ~0 and depositing from it would drop essentially the entire
plastic wake. The engine now exposes the full ledger for audit:
`N_em_pre_renewal`, `N_em_retained`, `N_em_shed_to_wake`, plus the pre-renewal
driving state `sigma_back_pre_renewal`, `r_eff_pre_renewal`,
`dG_emb_pre_renewal_eV` (the tip plasticity that actually drove the advance,
captured before the renewal reset). The 1D Kc summary is unaffected.

## Rate-shelf audit

Every run prints a shelf audit: the saturated clock rate `lambda_c_max` and the
minimum time-to-fire. If the clock cannot complete within the loading window it
says so **explicitly** — a dead clock is a reported condition, not a silent
mystery (the failure mode of the old runs). Co-design `nu0`, `tau_c`, and the
loading rate against the printed shelf.

## Open knob (yours to decide)

The barrier parameters here produce a DBTT near 700–800 K by construction. They
are **not** fit to data. To pin the model, fit the two barriers to W `Kc(T)`
(DFT or experiment) or to a target DBTT + brittle/ductile Kc. The stress
dependence of the entropy (`entropy_gate_power`, `sigma0_S` per channel, and
affine-vs-gated) is the main lever and is worth a dedicated sweep — it is also
the link to the low-T fatigue-limit behavior.

## What was NOT changed

`mesh.py`, `fem.py`, `plasticity.py`, `j_integral.py`, `materials.py`,
`config.py` (except the additive `entropy_stress_form` / `entropy_gate_power`
fields on `FractureBarrier`). The original `main.py` emergent driver is left
intact for comparison.
