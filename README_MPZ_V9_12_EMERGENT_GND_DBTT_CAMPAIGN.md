# MPZ v9.12 emergent-GND DBTT parameterization campaign

## Scope

This branch is a **parallel research formulation** based on
`v9.11.1-dbtt-mechanism-options`. It does not overwrite or silently activate the
v9.11.1 moving-process-zone model.

The new campaign asks whether a sharp temperature dependence of the
**post-initiation microstructural resistance** can emerge from one coupled
physical state:

1. crack-tip emission creates signed mobile dislocations;
2. the Arrhenius Peierls mechanism controls transport;
3. ordinary forest encounters store mobile content as retained content;
4. the Arrhenius Taylor completion mechanism releases retained content;
5. opposite-sign populations annihilate through a capture-radius law;
6. the signed retained GND field generates both local resolved internal stress
   and crack-tip shielding through fixed elastic kernels.

There is no fitted `N_sat`, shielding coefficient, independent Taylor-like
backstress, constitutive `K_shield` cap, or target on `K0(T)`.

## Campaign quantity

The screening quantity is

```text
Delta_K_micro(T) = median over 10-30 um of
                   [K_required(stateful) - K_required(neutral)]
```

The stateful and neutral evaluations use the same crack geometry, temperature,
cleavage surface, and target cleavage hazard. Only the evolved signed
microstructural feedback is removed in the neutral calculation.

`K0(T)` is not an argument of the objective and carries no target or penalty.

Default acceptance:

```text
Delta_K_micro amplitude        >= 8 MPa sqrt(m)
largest-jump localization      >= 0.50
10-90 percent transition width <= 200 K
```

A candidate also fails unless it develops nonzero signed GND content, GND
stress, and signed shielding.

## State variables and equations

For slip system `alpha`, sign `q`, and physical MPZ cell `i`, the model carries

```text
rho_m[alpha,q,i]  mobile dislocation density [m^-2]
rho_r[alpha,q,i]  retained dislocation density [m^-2]
rho_src[alpha,i]  available physical source density [m^-2]
```

The signed GND density is

```text
kappa[alpha,i] = rho_r[alpha,+,i] - rho_r[alpha,-,i].
```

The forest density is

```text
rho_f[alpha] = rho_floor + sum_beta a[alpha,beta] rho_total[beta].
```

The mean free path, Peierls velocity, and storage rate are

```text
l_mfp = c_mfp / sqrt(rho_f)
jump  = c_jump / [2 sqrt(rho_f)]
v_P   = jump * lambda_P
k_enc = abs(v_P) / l_mfp.
```

Taylor release uses the candidate Taylor EXP-floor barrier and the natural,
uncapped obstacle order

```text
m = 1 + 2 L_corr sqrt(rho_f)
lambda_T_completion = lambda_T_single / m.
```

Retained content evolves through encounter storage, Taylor release, recovery,
and opposite-sign annihilation. There is no logistic factor such as
`1-rho_r/rho_cap`. Any saturation-like state is an output of the kinetic
balance and residence time.

## One signed field, two mechanical projections

The same `kappa` field generates the internal resolved stress and shielding:

```text
tau_GND[alpha,i] = sum_beta,j K_tau[alpha,beta,i,j] kappa[beta,j]
K_shield          = sum_alpha,j K_K[alpha,j] kappa[alpha,j].
```

The 1-D screen includes regularized analytical edge-dislocation kernels. The
2-D adapter accepts mechanically measured kernels sampled at the same physical
MPZ cell centers. Kernels, slip interaction coefficients, source-zone width,
mean-free-path coefficient, and annihilation capture radius are common physics,
not candidate parameters.

## Physical source convention

The new row requires

```text
rho_source0_m2
```

rather than `source_sites_per_system`. Total source inventory is density times
physical MPZ cell area and is therefore invariant to the number of bins. Source
availability depletes through the emission hazard and refreshes only after
accepted crack advance exposes new material.

## Files

```text
arrhenius_fracture/emergent_gnd_types_v912.py
    Barrier, candidate, common-physics and protocol types.

arrhenius_fracture/emergent_gnd_state_v912.py
    Signed mobile/retained evolution, physical source density, Peierls
    transport, forest storage, Taylor release, annihilation, GND stress and
    shielding.

arrhenius_fracture/emergent_gnd_campaign_v912.py
    Matched stateful/neutral resistance, temperature objective and registry I/O.

arrhenius_fracture/emergent_gnd_dbtt_v912.py
    Stable public API.

arrhenius_fracture/emergent_gnd_2d_adapter_v912.py
    Trial/commit interface for a single active 2-D crack front.

scripts/generate_mpz_v9_12_sobol_candidates.py
scripts/build_mpz_v9_12_protocol_from_2d.py
scripts/run_mpz_v9_12_emergent_gnd_screen.py
mpz_v9_12_emergent_gnd_common_physics.json
mpz_v9_12_emergent_gnd_search_bounds.json
mpz_v9_12_protocol_example.csv
```

## Install in a new folder

```bash
cd /Volumes/Data/Data/Nanopillar_calculation

DEST=Arrhenius_FEM_CZM_MPZ_v9_12_emergent_gnd_dbtt

test ! -e "$DEST" || {
  echo "ERROR: destination already exists: $DEST"
  exit 1
}

git clone \
  --branch v9.12-emergent-gnd-dbtt-campaign \
  --single-branch \
  https://github.com/ukaiiaku-maker/Arrhenius_FEM_CZM_MPZ.git \
  "$DEST"

cd "$DEST"
conda activate arrhenius-fem-czm
python -m pip install -e . --no-deps
python -m pytest -q tests/test_emergent_gnd_dbtt_v912.py
```

## Build a matched protocol from a neutral 2-D trajectory

Run the 2-D model with microstructural feedback disabled but with the same
geometry, loading schedule, crack increments and signed-J convention. Convert
its trajectory:

```bash
python scripts/build_mpz_v9_12_protocol_from_2d.py \
  --input-csv runs/<neutral_case>/event_history.csv \
  --extension-column extension_um \
  --K-column K_J_MPa_sqrt_m \
  --time-column time_s \
  --out inputs/v9_12_neutral_2d_protocol.csv
```

The supplied `mpz_v9_12_protocol_example.csv` is only a smoke/example protocol.
It is not a calibrated elastic backbone.

## Generate a Sobol family

The v9.11.1 rows may be used as barrier seeds. Their legacy source counts are
not used by the v9.12 constitutive law.

```bash
python scripts/generate_mpz_v9_12_sobol_candidates.py \
  --base-registry mpz_v9_11_dbtt_option_rows.csv \
  --base-candidate-id DBTT_restart04_candidate03 \
  --bounds-json mpz_v9_12_emergent_gnd_search_bounds.json \
  --n 4096 \
  --seed 912 \
  --out candidates/v9_12_emergent_gnd_sobol_4096.csv
```

## Stage 1: homogeneous 0-D rejection gate

```bash
python -u scripts/run_mpz_v9_12_emergent_gnd_screen.py \
  --stage 0d \
  --candidate-registry candidates/v9_12_emergent_gnd_sobol_4096.csv \
  --protocol-csv inputs/v9_12_neutral_2d_protocol.csv \
  --physics-json mpz_v9_12_emergent_gnd_common_physics.json \
  --temperatures 300 400 500 600 700 800 900 1000 1100 1200 \
  --window-um 10 30 \
  --out runs/v9_12_emergent_gnd_0d_sobol_4096_v1
```

`0d` collapses the MPZ to one physical cell while preserving the same source,
storage, release, recovery and objective definitions. It is a topology screen,
not a spatial validation.

## Stage 2: full 1-D moving-PZ screen

```bash
python -u scripts/run_mpz_v9_12_emergent_gnd_screen.py \
  --stage 1d \
  --candidate-registry candidates/v9_12_emergent_gnd_0d_shortlist.csv \
  --protocol-csv inputs/v9_12_neutral_2d_protocol.csv \
  --physics-json mpz_v9_12_emergent_gnd_common_physics.json \
  --temperatures 300 400 500 600 700 800 900 1000 1100 1200 \
  --window-um 10 30 \
  --out runs/v9_12_emergent_gnd_1d_shortlist_v1
```

The runner prints `CAMPAIGN_START`, `CANDIDATE_START`, `CASE_RESULT`,
`CANDIDATE_COMPLETE`, and `CAMPAIGN_COMPLETE` directly to the terminal.

## Stage 3: 2-D coupling contract

For each mechanical trial:

1. the FEM/CZM driver supplies `K_applied`, `T`, and `dt`;
2. `EmergentGND2DAdapter.evolve_trial` advances a copy of the state;
3. the solver uses returned `K_shield` and `tau_GND` in the existing fixed point;
4. a rejected mechanical trial discards the copy;
5. an accepted trial calls `commit_trial`;
6. accepted crack advance is passed as `accepted_crack_advance_m`, translating
   the state and exposing pristine source-bearing material.

The first 2-D gate should retain one active front, branching off, 50 um / 80
bins, and the measured signed active-zone kernel.

## Required diagnostics

At each extension checkpoint retain

```text
Delta_K_micro
K_shield
tau_GND at the source/tip cells
signed GND density by system and cell
mobile and retained densities by sign
forest density
source available fraction
Pi_store = k_enc * residence_time
Pi_release = lambda_T_completion * residence_time
```

A candidate fails if a sharp resistance increment is not accompanied by a
coincident GND/internal-stress/shielding transition.

## Current limitations

- The 1-D analytical kernels are isotropic regularized edge-dislocation kernels.
  Production 2-D validation should replace them with mechanically measured
  tungsten/crystal-specific kernels.
- The initial common slip-interaction matrix is simple and fixed. It should be
  replaced by a geometry-derived BCC matrix before publication use, not fitted
  per candidate.
- Distributed full-field crystal plasticity is not introduced.
- Activation in the production FEM/CZM driver remains a separate validation
  gate; this branch provides the parallel campaign and trial/commit interface.
