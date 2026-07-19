# Audit: caps / floors / gates in the 1D sharp-front engine

Goal: confirm no artificial clamp *imposes* the regime selection. Each item is
classified PHYSICAL (keep), NUMERICAL-SAFE (keep, guards only), or
ARTIFICIAL-IMPOSING (must be justified or replaced).

| # | location | clamp | class | effect on regime physics |
|---|----------|-------|-------|--------------------------|
| 1 | sigma_tip | `K_eff = max(K - k_shield*N_em*..., 0)` | PHYSICAL | k_shield=0 by default; floor just prevents K<0. |
| 2 | sigma_tip | `s = min(s, sigma_cap=30 GPa)` | PHYSICAL | cohesive/theoretical-strength ceiling. Sets the cold brittle Kc scale; same for all regimes. Convergence checked (Kc flat vs Kmax). |
| 3 | sigma_back | (none) linear in N_em | — | unbounded by itself; bounded only through N_em. |
| 4 | e_stored/dG_emb | (none) linear in N_em | — | unbounded by itself; bounded only through N_em. |
| 5 | lambda_emit | `s_eff = max(sig_tip - sigma_back, 0)` | PHYSICAL | back-stress self-limits emission; cannot drive opening stress negative. Does NOT bound N_em at high T (zero-stress thermal channel stays open). |
| 6 | lambda_emit/cleave | `exp(clip(x,-700,0))` | NUMERICAL-SAFE | underflow guard + rate<=nu0 (barrierless cap). |
| 7 | lambda_cleave | `sig_eff = max(sig_tip - chi*sigma_back, 0)` | PHYSICAL | shielding reduces opening stress to at most zero (barrier to at most H0). This is the bounded shielding channel. |
| 8 | lambda_cleave | `Geff = max(Gstar - dGe, 0)` | PHYSICAL | barrier can't be negative (barrierless = max rate). |
| 9 | lambda_cleave | `dGe = min(dG_emb, emb_sat_frac*Gstar)` | **ARTIFICIAL-IMPOSING** | caps embrittlement ASYMMETRICALLY (bounds dG_emb but not sigma_back). With default emb_sat_frac=1 it still floors Geff at 0. This was the band-aid for the missing recovery; **superseded by N_sat/recover_k** (item 14). |
| 10| lambda_cleave | `gammainc(m, min(lam_raw*tau,1e12))/tau` | NUMERICAL-SAFE + MODEL | overflow guard is safe; m_hits & tau_c are genuine model params (multihit first-passage). tau_c sensitivity ~5%/decade (reported, not convergence). |
| 11| step | `dN <= dN_cap (=50)` | NUMERICAL (was entangled) | per-step emission throttle. **With recovery on it is irrelevant** (Kc identical 25->1e9). Without recovery it set the unbounded-growth rate, contaminating dt. |
| 12| step | `n_fire = floor(min(B,1e7))` | NUMERICAL-SAFE | renewal count (multi-advance/step, no throttle); cap only guards overflow on runaway. |
| 13| step | `if not isfinite(B): B=0` | NUMERICAL-SAFE | NaN guard from singular solves. |
| 14| step | **N_em recovery: `prod*=max(1-N_em/N_sat,0)`, `-recover_k*N_em*dt`** | PHYSICAL (NEW) | the fix. Defaults (inf, 0) reproduce as-shipped. Finite values bound N_em so dG_emb AND sigma_back saturate together -> shelf survival becomes a physical condition (saturated dG_emb vs H0), no per-barrier cap. |
| 15| driver | early-arrest break (`sig_eff<0.02*sig` 250 steps) | DIAGNOSTIC | declares arrest only when shielding has driven sig_eff to ~0; verified against full runs (does not cut a curve that would fire). |

## Verdict
The only clamp that was imposing the key physics was **#9 (emb_sat_frac)**, an
asymmetric barrier cap standing in for the absent dislocation recovery. Replacing
it with **#14 (physical N_sat / recover_k)** reproduces the regimes (verified: the
DBTT shelf survives with emb_sat_frac=1 and a finite N_sat), so the
ceramic/peak/weak-T/DBTT span is real physics. Recommended practice: run all
regime studies with `--emb-sat-frac 1` (cap off) and bound the ledger with
`--n-sat` or `--recover-k`.

## Still to revisit if/when moving to 2D
The deleted legacy full-field variational driver had its own clamps (`rho_cap`, `mobile_rho_floor`,
`peierls_floor`, mesh `da_phys` floor). They do not touch the 1D regime results
but must be audited before 2D regime claims.
