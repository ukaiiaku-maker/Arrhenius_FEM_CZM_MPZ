# OneD V2 Taylor/Peierls Search Contract

## Scope and lineage

This campaign is a monotonic-fracture search rooted at focused Peak/DBTT result commit `3d64c918431874f1ee1c215e4441b86672a52cc6`. It does not alter the completed four-class registry or the completed Peak/DBTT option registry. The authoritative control rows are loaded at full precision from the V2 material registry rather than transcribed:

| Class | Control candidate |
|---|---|
| Peak | `v913_zeroD_sobol_0242980` |
| DBTT | `v913_zeroD_sobol_0202500` |

Only these ten active coordinates may vary:

- Peierls: `peierls_H0_eV`, `peierls_activation_entropy_kB`, `peierls_exp_a`, `peierls_exp_n`.
- Taylor: `taylor_H0_eV`, `taylor_activation_entropy_kB`, `taylor_exp_a`, `taylor_exp_n`, `taylor_corr_rho_c_m2`, `taylor_corr_scale`.

All other active candidate coordinates are copied bit-for-bit from the class control. In particular, `peierls_nu0_s = 1e12 s^-1`, `taylor_nu0_s = 1e11 s^-1`, `rho_source0_m2`, and `c_blunt` are fixed. The source-owned Peierls and Taylor stress projections are both fixed at `1/sqrt(3)`. Common physics, mechanics-provider data, event lifecycle policy, source geometry, and inactive historical registry coordinates are not candidate variables.

No memory, shielding multiplier, wake-decay parameter, event counter, lifecycle parameter, or provider-specific material row is introduced. The PF and FEM/CZM reduced providers evaluate the same material row. No FEM/CZM calculation is authorized by this campaign.

## Immutable barrier fingerprints

The fingerprints hash IEEE-754 hexadecimal coordinate encodings, so a one-bit change fails the contract.

| Class | Cleavage barrier SHA-256 | Emission barrier SHA-256 |
|---|---|---|
| Peak | `50ead03bdc1dc5240eb2fea00eda4bcd24535887379f356a3105f3405852ed76` | `8197ce770333504c47b05ecaeb7497cf66ae6399cd8524360e3146bf2a9cc1f6` |
| DBTT | `ca51dcc632ec9d96dcda209413dbf3b7a914ecb488433a281d6359ddf663c156` | `75bc5a555c5412dfb9c8337c0c0686a3fdc3e9f5a4774cfe4671664f61a40186` |

The corresponding full active-vector control hashes are Peak `d004ea60768cfb90435d8d10f9badd180e614c0fc23950b24554b5a0b4bbe3a2` and DBTT `54104a30a31aa19517810249cb3ea42a346b92670ce3bd4f007110ba3babb1c2`.

## Direct source-equation audit

All ten whitelisted fields are active. Peierls and Taylor each construct an EXP-floor surface from their four searched shape/barrier coordinates and the fixed emission surface. Their fixed attempt frequencies then multiply the Arrhenius exponent. Searching activation entropy while keeping attempt frequency fixed avoids searching two partly degenerate prefactors.

The exact zero-dimensional persistent-state exchange is:

`Peierls rate -> Peierls velocity -> encounter rate`,

with the encounter rate creating retained content. The source's quantity named Taylor completion is the competing completion/release channel. Consequently:

- `tau_transport = current source/process-zone transport distance / abs(Peierls velocity)`;
- `tau_retention = 1 / encounter rate`;
- `chi_ret = tau_transport / tau_retention`;
- `tau_taylor_completion = 1 / Taylor completion rate`;
- `chi_taylor_completion = tau_transport / tau_taylor_completion`.

Both ratios are reported. Taylor completion is not mislabeled as the forward retention rate. The diagnostic calls the same source routines and current state geometry used by the evolution update; it adds no constitutive term and does not feed back into the trajectory.

At fixed geometry the encounter rate is proportional to Peierls velocity, so `chi_ret` can algebraically cancel the absolute velocity. The campaign therefore retains the dimensional Peierls rate/velocity, transport time, encounter time, Taylor completion time, both ratios, and the retained-equilibrium fraction. Selection will not rely on `chi_ret` alone.

## Resistance and decomposition policy

Only the first pre-event state and later pre-event states separated by a resolved reload are resistance candidates. Native drive values inside a physical avalanche remain trajectory diagnostics and are never labeled R-curve points.

Every newly produced onset record preserves remote/common-reference K, native K, local source stress, resolved emission drives, fixed-radius local-equivalent K, current radius, backstress, mobile and retained state, and the source-owned kinetic diagnostics. The zero-dimensional predictive state does not contain a separately identifiable signed `K_shield`; it is recorded as `NOT_REPRESENTED_IN_ZEROD_PREDICTIVE_STATE`, not imputed as zero. Direct PF field diagnostics are required before any shielding-dominated claim.

## Fail-closed controls

The producer validates each candidate before evaluation. It rejects any change outside the whitelist or either barrier fingerprint. Candidate identity is the SHA-256 of the canonical full active material vector, independent of search entry path. Tests cover the whitelist, bit-level barrier identity, fixed attempt frequencies/projections, source-rate reproduction, reload-separated onset extraction, and separation of apparent/native/local metrics.

Search bounds, transformed-coordinate definitions, population counts, source hashes, and committed producer revision will be appended to machine-readable provenance when the deterministic population is generated.
