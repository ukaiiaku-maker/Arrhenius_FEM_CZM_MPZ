# CODEX HANDOFF

## FEM/CZM–PF Sharp-Interface Parity for Arrhenius Hazard Fracture at (\theta=0^\circ)

## 1. Project objective

The primary objective is to make the two-dimensional FEM/cohesive-zone fracture model reproduce the same **qualitative physical behavior** as the existing sharp-interface Arrhenius hazard model in the `PF-fracture-fatigue` repository.

For the present development phase, restrict all calculations to:

[
\theta=0^\circ
]

The crack should advance along a straight nominal interface. There is no immediate requirement to return to the earlier (\theta=30^\circ) or other inclined-crack configurations.

The FEM/CZM implementation should reproduce the main trends of the PF model for four calibrated material parameterizations:

1. `v913_paper_peak01_0242980_persistent_sites`
2. `v913_paper_dbtt01_0202500_persistent_sites`
3. `v913_paper_weakT01_0129902_persistent_sites`
4. `v913_paper_ceramic01_0077080_persistent_sites`

The final FEM/CZM model does not need to match every stochastic event or every point of the PF trajectory exactly. It must, however, reproduce the same principal physical trends:

* comparable first-passage fracture resistance;
* qualitatively similar temperature dependence;
* qualitatively similar R-curve shape and magnitude;
* comparable ordering among the four parameterizations;
* stable and interpretable evolution of the cleavage hazard state (B);
* stable and interpretable evolution of the emission population (N_{\mathrm{em}});
* correct translation of the crack-tip microstructure as the crack advances;
* consistent crack-growth statistics and event lengths;
* no unphysical resets, loss of state, or remeshing-induced changes in constitutive history.

The FEM/CZM model should use the same Arrhenius hazard physics as the PF sharp-interface model. The cohesive-zone representation should be treated as a numerical representation of crack geometry and separation, not as an independent empirical fracture criterion.

---

## 2. Repositories and recommended new workspace

### Existing FEM/CZM repository

Current working repository:

```text
/Volumes/Data/Data/Nanopillar_calculation/Arrhenius_FEM_CZM_MPZ_v10_0_5_17_paper7_self_contained
```

GitHub repository:

```text
ukaiiaku-maker/Arrhenius_FEM_CZM_MPZ
```

Current development branch:

```text
v10.0.5.18.4.0-theta0-pf-parity
```

Associated draft pull request:

```text
PR #64
```

PR #64 must remain draft and unmerged until the full parity campaign has been completed and audited.

### Existing PF reference repository

```text
/Volumes/Data/Data/Nanopillar_calculation/PF-fracture-fatigue_v10_2_21_persistent_sites_top1
```

### Recommended new folder

Create a fresh clone rather than modifying the previous FEM/CZM working directory:

```bash
cd /Volumes/Data/Data/Nanopillar_calculation

git clone \
  https://github.com/ukaiiaku-maker/Arrhenius_FEM_CZM_MPZ.git \
  Arrhenius_FEM_CZM_MPZ_theta0_pf_parity_codex

cd Arrhenius_FEM_CZM_MPZ_theta0_pf_parity_codex

git fetch origin
git switch v10.0.5.18.4.0-theta0-pf-parity
```

Before beginning development, verify the current remote head of PR #64. Do not assume an older commit hash from a previous terminal transcript is still authoritative.

Open the new folder in VS Code:

```bash
code .
```

Use the existing environment initially:

```bash
conda activate arrhenius-fem-czm
python -m pip install -e '.[dev]'
```

Do not copy old `runs/` directories into the new source tree. Reference them through absolute paths where needed.

---

## 3. Authoritative PF reference case

The first required comparison is the Peak parameterization at 1000 K:

```text
Parameterization:
v913_paper_peak01_0242980_persistent_sites

Temperature:
1000 K

Orientation:
theta = 0 degrees

Hazard seed:
8666

Target projected crack extension:
1000 micrometers
```

Reference case:

```text
/Volumes/Data/Data/Nanopillar_calculation/PF-fracture-fatigue_v10_2_21_persistent_sites_top1/runs/v10_4_1_theta0_rate1x_bulk_PT_four_class_1000um_selective_reuse_base3621_v1/v913_paper_peak01_0242980_persistent_sites/T1000K_th0_seed8666
```

The reference step table is:

```text
steps_1000K.csv
```

The PF reference uses:

```text
dU                  = 2.0e-7 m
dt                  = 8.4 s
opening rate        = 2.3809523809523807e-8 m/s
n_stagger           = 2
tip_h_fine          = 1.0e-6 m
tip_ratio           = 1.20
da_phys             = 5.0e-6 m
adaptive target     = 0.15
maximum fronts      = 1
bulk mode           = full_field
```

The exact signed shielding-kernel family used by this PF case is located at:

```text
/Volumes/Data/Data/Nanopillar_calculation/PF-fracture-fatigue_v10_2_21_persistent_sites_top1/runs/v10_2_28_kernel_cache/1447653d199f0b43cb475951092d69444c9b785f6fdf518c723792abb3b1f5e5/family.json
```

Expected SHA-256:

```text
a85b57ad9eee8331ef34ea222760ce7d4064f48ddef668f1e28a666d14946f3a
```

The source PF artifact is an active-only kernel family. Its wake channel is disabled. The FEM compatibility layer may adapt the empty disabled-wake representation for the legacy loader, but it must not modify the active kernel, enable wake shielding, or alter the kernel normalization.

---

## 4. Current FEM/CZM architecture

The current FEM/CZM stack contains the following relevant capabilities:

* persistent crack-tip source population;
* moving crack-tip microstructure or MPZ;
* stochastic cleavage first passage;
* stochastic emission;
* signed shielding kernel;
* full-field Peierls–Taylor bulk plasticity;
* forward-minus-reverse detailed balance;
* exact zero net plastic rate at zero stress;
* atomic cohesive crack extension;
* exact preservation of the hazard-selected event endpoint;
* constrained path-corridor remeshing;
* rollback when a full event cannot be committed;
* minimum triangle-quality and local-area-ratio gates;
* quantitative PF–FEM comparison utilities.

The v10.0.5.18.3.9 remeshing work demonstrated that the cohesive crack can propagate through approximately 400 µm without the previous mesh failure. The remesher is therefore no longer the primary known blocker.

The current physics-parity work is in the v10.0.5.18.4.0 branch.

Relevant files include, or should include equivalents of:

```text
arrhenius_fracture/
  bulk_pt_detailed_balance_v1041_exact.py
  mode_i_first_passage_v10_0_5_18_4_0_theta0_pf_parity.py
  mode_i_first_passage_v10_0_5_18_4_0_theta0_pf_full_field_parity.py
  mode_i_first_passage_v10_0_5_18_4_0_theta0_pf_full_field_production.py
  theta0_pf_full_field_metadata_v10051840.py

run_v10_0_5_18_4_0_peak1000_theta0_pf_full_field_parity.sh

scripts/
  compare_v10051840_theta0_pf_parity.py
  preflight_and_launch_v10051840_theta0_full_field.sh
```

The exact current file list should be checked from the branch rather than inferred from this handoff.

---

## 5. Current status and identified problem

The full-field (\theta=0^\circ) preflight has successfully demonstrated:

* correct installation of package version `10.0.5.18.4.0`;
* passage of focused and inherited tests;
* successful construction of a real FEM/CZM model;
* successful execution of a one-step full-field calculation;
* use of the Peak `0242980` parameter row;
* use of seed 8666;
* use of the atomic path-corridor cohesive backend;
* active full-field bulk plasticity;
* correct detailed-balance bulk formulation;
* use of the PF active-only kernel through an audited compatibility layer.

A longer calculation was then launched using the PF remote-opening rate.

That calculation was stopped and archived because it traversed the prefracture loading regime too slowly.

Archived baseline:

```text
/Volumes/Data/Data/Nanopillar_calculation/Arrhenius_FEM_CZM_MPZ_v10_0_5_17_paper7_self_contained/runs/v10_0_5_18_4_0_peak1000_theta0_seed8666_pf_full_field_parity_v1_slow_displacement_ramp_baseline_20260803_202328
```

At approximately step 3600, the calculation had reached only:

```text
KJ approximately 8.48 MPa sqrt(m)
B  = 0
N_em approximately 0
crack extension = 0
```

The PF first-passage state is at a much larger driving force, approximately:

```text
J approximately 6531 J/m2
KJ approximately 54 MPa sqrt(m)
```

A direct extrapolation suggested that the FEM calculation could require roughly 145,000 load steps and several days merely to reach first passage.

The important conclusion is:

> Matching the PF remote displacement rate does not guarantee matching the PF driving-force rate.

The PF and FEM/CZM geometries have different compliance. Therefore, the same (dU/dt) gives different (J(t)) and (K_J(t)).

This must be addressed before evaluating the R-curve.

---

## 6. Immediate development priority: driving-force-controlled loading

The next implementation should replace, or supplement, the fixed remote-displacement ramp with a PF-reference driving-force controller.

### Desired behavior

Read the prefracture PF trajectory from the archived `steps_1000K.csv` and construct a target history such as:

[
J_{\mathrm{target}}(t)
]

or, if more robust,

[
K_{J,\mathrm{target}}(t).
]

At each nominal PF time interval:

1. evaluate the current FEM/CZM (J);
2. predict the boundary displacement required to reach the next PF target (J);
3. solve the FEM system;
4. update the displacement using a safeguarded secant, proportional, or square-root correction;
5. converge to the target (J) within a specified tolerance;
6. integrate all Arrhenius hazard and plasticity state over the actual physical time interval;
7. localize any stochastic event inside the interval;
8. preserve the exact event threshold, event length, and event endpoint.

The controller must not simply jump the hazard state from one target point to the next.

### Suggested predictor

In an approximately elastic prefracture regime:

[
J \propto U^2.
]

A useful first predictor is therefore:

[
U_{\mathrm{new}}
================

U_{\mathrm{old}}
\sqrt{
\frac{J_{\mathrm{target}}}
{\max(J_{\mathrm{old}},J_{\min})}
}.
]

This should be followed by a safeguarded correction because full-field plasticity and crack growth will invalidate exact square-root scaling.

### Event localization

If the integrated cleavage or emission hazard crosses a stochastic threshold between two target-(J) states, the code must localize the event within the interval.

Do not:

* assign the event to the end of the interval;
* redraw a new threshold;
* split one physical event into several artificial physical events;
* reset (B), (N_{\mathrm{em}}), or the MPZ because the load controller retries a solve.

A trial target-(J) solve must be transactional. Rejected load-controller iterations must restore:

* bulk plastic strain;
* bulk forest density;
* tip MPZ populations;
* (B);
* current cleavage threshold;
* current emission thresholds;
* event-length draw;
* crack geometry;
* cohesive state;
* RNG state.

---

## 7. Non-negotiable physical requirements

### 7.1 Fracture remains hazard controlled

There must be no independent athermal fracture threshold inserted to force crack advance.

The crack event is selected by the Arrhenius cleavage hazard.

The cohesive-zone model must not introduce a competing empirical initiation criterion that overrides the stochastic hazard.

### 7.2 Exact event geometry

For every stochastic cleavage renewal:

* preserve the selected total event length;
* preserve the selected direction;
* preserve the exact endpoint;
* allow numerical subsegments only as a representation of one physical event;
* do not replace one event with equal-length 2, 4, or 8-event partitions.

### 7.3 Atomic commit and rollback

If the cohesive geometry cannot be created while satisfying all mesh gates, the complete physical event must fail atomically and restore the prior state.

### 7.4 Mesh-quality gates

Do not weaken the existing acceptance thresholds:

```text
minimum triangle quality >= 0.035
minimum local child-area ratio >= 0.08
```

The exact local-area metric should remain the current path-support-cell metric.

### 7.5 Crack-tip microstructure

The moving MPZ must:

* translate with realized physical crack-path length;
* remain synchronized with the current physical crack tip;
* preserve active population, retained population, accumulated slip, and backstress correctly;
* discard only state that physically exits through the wake boundary;
* not reset during remeshing;
* not be advanced during a rejected trial;
* not be translated twice for one event;
* not be translated only after geometry commit if the constitutive split requires fractional translation during hazard evolution.

### 7.6 Tip and bulk populations remain distinct

The PF v10.4.1 model uses:

* a moving crack-tip MPZ;
* a full-field Peierls–Taylor bulk model;
* no direct deposition of the tip population into the continuum bulk density;
* no artificial bulk-density transport from crack advance.

Preserve that separation.

### 7.7 Full-field bulk law

The FEM implementation should reproduce the PF v10.4.1 bulk closure:

* selected-row Peierls barrier;
* selected-row Taylor barrier;
* correlated Taylor hit-order law;
* 1 ns Taylor renewal window;
* PF mobile-density law;
* PF jump-length law;
* forward-minus-reverse detailed balance;
* exactly zero net plastic rate at zero stress;
* `bulk_mult_frac = 1`.

---

## 8. Required crack-tip state audit

The current concern is that (B), (N_{\mathrm{em}}), and the moving MPZ were less stable in the FEM result than in the PF result.

Add an event-resolved audit for every physical cleavage event.

Record:

```text
event index
solver step
physical time
target J
achieved J
target KJ
achieved KJ
B before increment
B after first constitutive half-step
B after translation
B before event commit
B after event renewal
current cleavage threshold
event length draw
fractional translation distance
committed crack-path distance
projected crack extension
tip coordinates before and after
N_em before and after
mobile MPZ population before and after
retained MPZ population before and after
accumulated-slip measure before and after
wake-discarded mobile population
wake-discarded retained population
wake-discarded accumulated-slip history
bulk plastic work increment
cohesive work increment
elastic-energy change
external-work increment
energy-balance residual
```

The audit should make it possible to determine whether an apparent population collapse is caused by:

* legitimate MPZ translation;
* a renewal reset;
* remeshing;
* rollback failure;
* double translation;
* diagnostic mixing of physical line content and accumulated-slip history.

---

## 9. Required mechanics and (J)-integral audit

A prior (\theta=30^\circ) FEM result produced a decreasing R-curve while the PF result produced a rising R-curve. The discrepancy appeared more strongly in (K_J/F) than in external force itself.

The parity model must therefore output:

```text
signed J
effective positive J
KJ
top reaction force
KJ / abs(force)
all available J contours
selected contour or cluster
current J-contour center
current physical crack-tip coordinates
distance between contour center and crack tip
external work
recoverable elastic energy
bulk plastic dissipation
tip plastic dissipation
cohesive-interface work
fracture-surface work
energy-balance residual
```

Determine whether the cohesive crack representation causes differences through:

* contour intersection with cohesive faces;
* omitted cohesive work;
* use of a stale tip location;
* incorrect contour recentering after remeshing;
* wake compliance;
* path-dependent state transfer;
* local-to-cluster (J) handoff.

Do not tune the Arrhenius parameters to compensate for a defective (J)-integral.

---

## 10. Comparison strategy

### Stage A: Peak, 1000 K, (\theta=0^\circ)

This is the primary debugging case.

Required result:

* reach first passage in a practical number of FEM solves;
* complete at least 400 µm before assessing long-run stability;
* ultimately complete 1000 µm;
* produce a rising or otherwise qualitatively PF-like R-curve;
* maintain stable MPZ and hazard state evolution.

### Stage B: Peak temperature sweep

After the 1000 K case is physically credible, run the Peak model over the PF reference temperatures:

```text
300
600
800
900
950
1000
1050
1100
1150
1200
1250
1300 K
```

Compare:

* first-passage (J);
* first-passage (K_J);
* maximum or plateau R-curve resistance;
* crack-growth trajectory;
* number of events;
* total emitted population;
* temperature of the toughness maximum or transition region.

### Stage C: Four parameterizations

Repeat for:

```text
Peak
DBTT
weakT
ceramic
```

The FEM model should reproduce the qualitative ordering and temperature dependence of the PF model.

The goal is not an arbitrary fit of each FEM result independently. All four cases must use the same numerical and physical interpretation, with differences coming from the selected parameter rows.

---

## 11. Quantitative comparison outputs

Use, extend, or replace:

```text
scripts/compare_v10051840_theta0_pf_parity.py
```

Align PF and FEM trajectories by projected crack extension, not by solver step.

Generate:

```text
aligned_R_curve.csv
parity_summary.json
KJ_R_curve_comparison.png
KJ_over_force_comparison.png
B_comparison.png
N_em_comparison.png
event_state_comparison.csv
energy_balance_comparison.csv
temperature_summary.csv
four_class_summary.csv
```

Suggested metrics:

### First passage

[
\epsilon_{J,\mathrm{FP}}
========================

## \frac{J_{\mathrm{FP}}^{\mathrm{FEM}}

J_{\mathrm{FP}}^{\mathrm{PF}}}
{J_{\mathrm{FP}}^{\mathrm{PF}}}
]

### R-curve error

* RMSE in (K_J) over common extension;
* median absolute relative error;
* slope sign agreement;
* slope magnitude ratio;
* resistance at 50, 100, 200, 400, and 1000 µm.

### Hazard-state stability

* maximum negative reset in (B) outside valid renewal;
* median (|\Delta B|);
* (N_{\mathrm{em}}) dynamic range;
* large isolated population-jump count;
* line-content conservation residual;
* MPZ translation-to-crack-path error.

### Mechanics

* (K_J/F) evolution;
* contour-to-contour spread in (J);
* energy-balance residual;
* cohesive-work fraction;
* plastic-work fraction.

---

## 12. Acceptance criteria

A case should not be described as parity-qualified merely because the solver completes.

### Numerical acceptance

* no uncaught exceptions;
* no geometry vetoes accepted as physical events;
* every committed stochastic event succeeds atomically;
* endpoint error within the existing v3.9 tolerance;
* event-length error within the existing v3.9 tolerance;
* triangle quality remains at least 0.035;
* local area ratio remains at least 0.08;
* no unresolved protected-node remapping;
* no constitutive state change during failed remesh trials;
* requested projected extension completed.

### Crack-tip state acceptance

* MPZ translation equals realized physical path advance within numerical tolerance;
* no unexplained MPZ reset;
* no double translation;
* no remeshing-dependent population loss;
* (B) changes only through integrated hazard and valid renewal;
* (N_{\mathrm{em}}) changes only through the defined emission, escape, retention, recovery, or wake-discard mechanisms.

### Physical parity acceptance

For each parameterization:

* same qualitative temperature trend as PF;
* same qualitative R-curve type;
* same sign of the developed R-curve slope over the principal growth range;
* comparable first-passage resistance;
* comparable resistance scale over common crack extension;
* correct relative ordering of the four model classes;
* no parameter-specific numerical workaround.

Exact stochastic event-by-event identity is not required unless the two implementations are intentionally driven with an identical random stream and identical hazard history.

---

## 13. Development rules for Codex

1. Read the existing handoff, branch history, tests, and production manifests before editing.
2. Do not weaken existing guards to make a run complete.
3. Do not change physical parameters without identifying the exact selected-row field and explaining why the PF implementation differs.
4. Do not introduce an athermal crack-initiation criterion.
5. Do not silently replace the PF kernel.
6. Do not describe metadata normalization as a physics fix.
7. Do not treat numerical subsegments as separate physical events.
8. Do not delete or overwrite prior run directories.
9. Create a new branch for each substantial redesign.
10. Add a regression test for every bug found.
11. Run a real production-entry smoke, not only unit tests.
12. Keep PR #64 draft until long-run parity evidence exists.
13. Record the Git commit, package version, parameter-row hash, kernel hash, seed, and exact launch contract in every run.
14. Prefer a small number of decisive controlled calculations over a broad campaign with unresolved physics.

---

## 14. Recommended immediate task sequence

### Task 1: Establish a reproducible baseline in the new folder

* install the current branch;
* run all focused tests;
* verify the exact PF reference files;
* parse the archived slow-ramp FEM result;
* reproduce the previously observed low-(J) ramp rate.

### Task 2: Implement PF (J(t))-controlled loading

* add a target-(J) trajectory reader;
* add transactional trial solves;
* add safeguarded displacement correction;
* add event localization within the target-(J) interval;
* verify no state advances during rejected controller trials.

### Task 3: Run prefracture parity only

Before attempting 1000 µm growth:

* reproduce the PF prefracture (J(t));
* compare (B(t));
* compare (N_{\mathrm{em}}(t));
* compare bulk plastic work;
* compare first-passage time and resistance.

### Task 4: Run a short crack-growth qualification

Use a 20 or 50 µm target:

* verify exact event geometry;
* verify MPZ translation;
* verify post-event (B) renewal;
* verify (J)-contour recentering;
* verify energy balance.

### Task 5: Run 400 µm and then 1000 µm

Only proceed when the short qualification passes.

### Task 6: Expand to the four-class temperature campaign

Do not tune the four classes independently. Preserve the registry rows and common numerical treatment.

---

## 15. Definition of success

The project is successful when the FEM/CZM model can be run at (\theta=0^\circ) for all four calibrated parameterizations and produces:

* a physically credible temperature-dependent first-passage response;
* a qualitatively PF-like R-curve;
* stable Arrhenius hazard and emission-state evolution;
* correct moving-tip microstructure transport;
* reliable 1000 µm crack propagation;
* no remeshing-induced constitutive artifacts;
* no empirical fracture threshold added to force agreement;
* a transparent accounting of the remaining difference between a PF sharp wake and a cohesive crack.

The cohesive-zone implementation is allowed to produce modest quantitative differences from the PF sharp-interface model. It is not acceptable for it to reverse the R-curve trend, produce unexplained state collapses, or require parameter changes that compensate for an incorrect loading, (J)-integral, MPZ, or remeshing implementation.
