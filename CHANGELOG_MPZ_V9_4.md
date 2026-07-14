# MPZ v9.4 changelog

## Detailed-balance correction

The v9.3 Peierls–Taylor screening run exposed a zero-stress plastic-flow ratchet. The v9.3 constitutive code evaluated forward Peierls and Taylor Arrhenius rates without subtracting the reverse rate. Because an EXP-floor barrier is finite at zero stress, the model could carry the prescribed plastic strain rate at essentially zero applied stress.

Version 9.4 restores the signed kinetic construction used in the prior DDD implementation:

```text
R_P,net = R_P,forward - R_P,reverse
R_T,net = R_T,forward - R_T,reverse
R_PT = (1/R_P,net + 1/R_T,net)^(-1)
```

The reverse reference uses the corresponding zero-stress EXP-floor barrier. Net Peierls motion, Taylor completion, the sequential Peierls–Taylor rate, and plastic strain rate are therefore exactly zero at zero effective stress.

## Search safeguards

- Added forward, reverse, and net-rate diagnostics for Peierls and Taylor branches.
- Added a strict zero-stress detailed-balance acceptance test.
- Added rejection when an unclamped scaled Peierls or Taylor zero-stress free energy becomes non-positive over the screened temperature range.
- Narrowed the plastic temperature-slope multiplier prior from `0.25–128` to `0.25–8`.
- Reduced the broad maximum-stress acceptance ceiling from `80 GPa` to `40 GPa`.
- Retained the correlated multi-hit Taylor closure, high-density monotonicity audit, and prohibition on a constitutive total-dislocation-density cap.

## Interpretation of the v9.3 v2 run

The completed v9.3 v2 archive is retained as an audit dataset. It demonstrates that monotonicity alone is insufficient: a constitutive branch can be monotonic in density while still being unphysical because it supports nonzero net plastic flow at zero stress. Its shortlist is not a production parameter set and must not be passed to developed-state or transient MPZ calculations.
